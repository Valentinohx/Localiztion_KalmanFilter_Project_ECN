/*
            kalmanfilter.cpp
            created at 2018-05-15
            Xiao SHI && Huaxin LIU
            objective: use the robot's encoder and the sensor as input, develope an hybrid localazation system to
                       localize the robot in the environment that with tiles floor, the main algorithm used is EKF
                       
            Subscribe : 1. /capteur_carrelage, in which the sensors output has been transformed into binary
                        2. /mobile_base/sensors/core , in which the robot's encoder is read from its internal sensor
                        
            Publish   : 1. /filted_state, which is the output of the filter, include its x, y and theta
                        2. /P_variance, which is the filter's variance that described the noises of the system
                        
            Before using: 1. modify the path of matplotlibcpp.h, which can download on https://github.com/lava/matplotlib-cpp
                          2. modify the path of saving graph
                          3. modify the length of tiles
                          4. modify the sensor position in robot frame
                          5. modify the geometry parameters of robot
        */

#include "ros/ros.h"
#include <std_msgs/Bool.h>
#include <geometry_msgs/Pose2D.h>
#include <std_msgs/Float32MultiArray.h>
#include <kobuki_msgs/SensorState.h>
#include <vector>
#include <math.h>
#include <Eigen/Dense>				//matrix library
#include <iostream>
#include <stdlib.h>
#include <std_msgs/Int8MultiArray.h>

//the absolute path of matplotlibcpp.h, which can download on
#include </home/valentinohx/ros/src/kalmanfileterwithrosbag/src/matplotlibcpp.h>
//#include </home/xiao/ros_lab_ws/src/kalman_new/src/kalmanfileterwithrosbag/src/matplotlibcpp.h>

using namespace std;
using namespace Eigen;
namespace plt = matplotlibcpp;


const char *absolutepath = "/home/valentinohx/ros/src/kalmanfileterwithrosbag/";
const char *folder = "plot_results_E/";
const char *bagName = "circle_1";
const char *comment = "406";
double tile_length = 0.30;					//tile size (bat. D : 0.1017, bat. E : 0.30)

/* **********************************for filter tuning********************************* */
//
//P
double initial_uncertainty_xy  =  0.007;    //0.007                 //Uncertainty about the initial position x y coordinates
double initial_uncertainty_theta  =  5;     //5                //Uncertainty about the initial position
//Qalpha
double state_uncertainty_xy = 0.003;        //0.003                    //state noise, in m 
double state_uncertainty_theta = 0.2;       //0.2                     //in degree

//Qbeta
double thegma_wheels = 0.00006;   // in m    //0.00006                  //input noise

//Qgama
double measurement_uncertainty = 0.005;      //0.005               //measurment uncertainty, in m
/* **********************************for filter tuning******************************** */
double initial_theta = 0; //in degrees 

const double PI = 3.14159265359;

ros::Publisher publisher_state;
ros::Publisher publisher_variance;
ros::Publisher publisher_capteur;

//Main variables
vector<VectorXd> measurement_equation;			//estimated measurements of sensors, sensor measurement
VectorXd X;                         //state, with dynamic size, type is double, column vector
VectorXd U;                         //moving under the wheels
VectorXd X1;                         //state, with dynamic size, type is double, column vector
VectorXd U1;                         //moving under the wheels
VectorXi Encoders;                  //encodage des positions of wheels, interger with dynamic size
VectorXd deltaq;
//Matrices de Kalman
MatrixXd A;                         //matrice A of Kalman, with double type, and dynamic matrix
MatrixXd B;                         //matrice B of Kalman
MatrixXd C_x;                       //matrice C of Kalman of X
MatrixXd C_y;                       //matrice C of Kalman of Y

MatrixXd P;                     	//propagation error
MatrixXd K;
MatrixXd jointToCartesian ;


//Matrices of variances
MatrixXd Qalpha;
MatrixXd Qbeta;
double  Qgamma;

//Threshold of Mahalanobis
// use chi2inv to calculate it
double mahaThreshold = 2.7;
//Geometric and temporal variables
double T = 0.02;							//update period of odometry
//double T = 0.1;							//update period of odometry
double t0;									//system time of the beginning of the program
double encoder_coef = 2578.33;				//conversion factor encoder -> angle (increment / wheel circumference)
double rr = 0.035;							//the radius of the wheel
double trackGauge = 0.227;

//sensor number and position definition
double mys = 0.075;							//distance centre-capteur selon y, in robot fram
double mxs = 0;                             //distance centre-capteur selon x, in robot fram

vector< vector<double> > m_sensor_position = {{mxs, mys}, {mxs,-mys}};  //order should be left to right strictly
vector< vector<double> > m_position_active_sensor;

vector<double> filted_x;
vector<double> filted_y;
vector<double> odom_x;
vector<double> odom_y;

vector< int> dMaha_dots_index;

vector< double> dMaha_green1_dots;
vector< double> dMaha_green2_dots;

vector< double> dMaha_red1_dots;
vector< double> dMaha_red2_dots;

vector<double> dMaha_neighboor_red_dots1;
vector<double> dMaha_neighboor_red_dots2;
vector<double> dMaha_neighboor_red_dots3;
vector<double> dMaha_neighboor_red_dots4;

vector<double> travel_distance;
double last_travel_distance = 0;

vector<double> thegema_x;
vector<double> thegema_y;
vector<double> thegema_theta;

int sensor_number = m_sensor_position.size();

inline VectorXd EvolutionModelCore( VectorXd &X, const VectorXd &U )
{
    X[0] = X[0] + U[0] * cos(X[2]);
    X[1] = X[1] + U[0] * sin(X[2]);
    X[2] = X[2] + U[1];
    return X;
}

//---------------------------------------- EVOLUTION MODELE----------------------------------------//
void evolutionModel()
{
    //store the last state of the encoder to get the change of angle, and avoid to miss every increment
    static int left_encoder_prec;    //static file scope: entire file, and it remains in memory until the program ends
    static int right_encoder_prec;
    //the first value of the encoder is abnormal, so we do not take it consideration
    //This modification imposes a pause at the beginning of the robot's tajectory in the command
    if(ros::Time::now().toSec()-t0> 0.4)
    {
        //Encoder increment variation variables
        int delta_left = Encoders[0]-left_encoder_prec;
        int delta_right = Encoders[1]-right_encoder_prec;
        
        //We take into account the limit of the coders (2^16) at the end of which they loop to 0
        if (delta_left > 1<<15)   //priority : << 7, > 9, so here << first, then >  2^15
            delta_left += -1<<16;  //+=: 16, -1<<16  -2^16 = -65536
        else if (delta_left < -1<<15)
            delta_left += 1<<16;
        
        if (delta_right > 1<<15)
            delta_right += -1<<16;
        else if (delta_right < -1<<15)
            delta_right += 1<<16;
        
        //increment radians of the wheel
        VectorXd deltaq(2);   //you must specify the size before you use it
        deltaq[0] = delta_right * 2 * PI / encoder_coef;
        deltaq[1] = delta_left  * 2 * PI / encoder_coef;
        U = jointToCartesian * deltaq;       //from localazation.pdf equation (5.19)
        X = EvolutionModelCore( X, U );  //prediction phase
        //pure odometry
        U1 = jointToCartesian * deltaq;
        X1 = EvolutionModelCore(X1,U1);
    }
    //update the encoder
    left_encoder_prec  = Encoders[0];
    right_encoder_prec = Encoders[1];
}

//-------------------------------------- CALLBACKS FUNCTION ------------------------------------------//

// Measurement from the sensor
void mesureCallback(std_msgs::Int8MultiArray IR_state)
{
    VectorXd g(2);
    for(int i = 0; i < sensor_number; i++ )
    {
        // if there is a new line detected, then update the measurements
        if ( IR_state.data[i])
        {
            //we do not know which kind of line is detected, so calculate for both, then use Mahalanobis distance to tell
            g[0] = X[0] + m_sensor_position[i][0]*cos(X[2]) - m_sensor_position[i][1]*sin(X[2]); // corresponds to the current line in world fram, estimated X line
            g[1] = X[1] + m_sensor_position[i][1]*cos(X[2]) + m_sensor_position[i][0]*sin(X[2]); // estimated Y line cordinates in the world fram
            m_position_active_sensor.push_back({m_sensor_position[i][0], m_sensor_position[i][1]});  //just keep the coordinates of the actived sensor
            measurement_equation.push_back(g);
        }
    }
}
// Update the raw data of the encoders
void vitesseCallback(kobuki_msgs::SensorState encoder_state)
{
    Encoders[0] = encoder_state.left_encoder;
    Encoders[1] = encoder_state.right_encoder;
}

//-------------------------------------- MAIN FONCTION ------------------------------------------//
int main(int argc, char **argv)
{
    //Initialise variables
    X = VectorXd::Zero(3);
    X1 = VectorXd::Zero(3);
    U = VectorXd(2);
    U1 = VectorXd(2);
    
    Encoders = VectorXi(2);
    A = MatrixXd(3,3);
    B = MatrixXd(3,2);
    C_x = MatrixXd(1,3);
    C_y = MatrixXd(1,3);
    P = MatrixXd(3,3);
    measurement_equation = vector<VectorXd>();
    
    jointToCartesian = MatrixXd(2,2);
    jointToCartesian <<         rr/2,             rr/2,
            rr/trackGauge,   -rr/trackGauge;
    //Pose of depart, modify it based on real trajectory
    X << 0, 0, initial_theta*PI/180;
    X1 = X;
    
    //Uncertainty about the initial position
    P <<    pow(initial_uncertainty_xy,2), 0, 0,
            0, pow(initial_uncertainty_xy,2), 0,
            0, 0, pow(initial_uncertainty_theta*PI/180, 2);
    
    //Matrices of variance
    Qalpha = MatrixXd(3,3);  //state noises
    Qalpha <<  pow(state_uncertainty_xy,2)*T, 0, 0,
               0, pow(state_uncertainty_xy,2)*T, 0,
               0, 0, pow(state_uncertainty_theta*PI/180,2)*T;
    
    Qbeta = MatrixXd(2,2);				//input noises
    Qbeta << pow(thegma_wheels,2), 0,
              0, pow(thegma_wheels,2);
    
    Qgamma = pow(measurement_uncertainty,2)/12;			//measurment noises
    
    //Initialisation ROS
    ros::init(argc, argv, "conversion_capteur");
    ROS_INFO("Node conversion_capteur_node connected to roscore");
    ros::NodeHandle nh;
    
    //Subscribing
    ROS_INFO("Subscribing to topics\n");
    ros::Subscriber subscriber_capteur = nh.subscribe<std_msgs::Int8MultiArray>("/capteur_carrelage", 1, mesureCallback);
    ros::Subscriber subscriber_roue = nh.subscribe<kobuki_msgs::SensorState>("/mobile_base/sensors/core" , 1, vitesseCallback);
    
    //Publishing
    publisher_state = nh.advertise<geometry_msgs::Pose2D>("/filted_state", 1);
    publisher_variance = nh.advertise<geometry_msgs::Pose2D>("/P_variance", 1);
    //publisher_capteur = nh.advertise<geometry_msgs::Pose2D>("/position_capteur", 3);
    
    //frequency
    ros::Rate loop_rate(1/T);
    
    //time and index initial
    t0 = ros::Time::now().toSec();
    int index = 0;
    while (ros::ok())
    {
        ros::spinOnce();
        //odometry
        evolutionModel();
        //use to store travel distance
        last_travel_distance = U(0) + last_travel_distance;
        travel_distance.push_back( last_travel_distance );
        
        odom_x.push_back(X1[0]);  //store the pure odometry
        odom_y.push_back(X1[1]);
        
        //calculate the Jacobian matrix of the system
        A << 1, 0, -U[0]*sin(X[2]),
                0, 1, U[0]*cos(X[2]),
                0, 0, 1;
        B << cos(X[2]), 0,
                sin(X[2]), 0,
                0,     1;
        //---------------------------------------------- FILTRE OF KALMAN ------------------------------------------//
        //propagation error
        P = A*P*A.transpose()+B*Qbeta*B.transpose()+Qalpha;
        /*When a measurement is taken by one of the sensors, two distances of Mahalanobis are calculated:
                  one considering that the line detected is horizontal, the other vertical
                  if only one is below the threshold, it is retained, one uses the corresponding matrix C */
        
        for (int i = 0; i < measurement_equation.size(); i++) //measurement_equation is used to store the output of measurment equation
        {
            /*transfer the sensor measurments to exactly in the line, either a horizontal line or a vertical line, the O_IR_P_exact_tileed should be the
                      times of the length of the tiles, this is true based on that the robot always starts from the limits of a tile*/
            VectorXd O_IR_P_exact_tile(2);
            O_IR_P_exact_tile[0] = round(measurement_equation[i][0]/tile_length)*tile_length;		//Sensor position rounded to a tile according to x
            O_IR_P_exact_tile[1] = round(measurement_equation[i][1]/tile_length)*tile_length;		//position du capteur arrondie à un carreau selon y
            
            //get the neighboors of the rounded line, this four more lines are used to calculate the  Mahalanobis distances
            //right or above line
            VectorXd O_IR_P_exact_tile_higher(2);
            O_IR_P_exact_tile_higher[0] = ( round(measurement_equation[i][0]/tile_length) + 1 )*tile_length;   //x higher neighboor
            O_IR_P_exact_tile_higher[1] = ( round(measurement_equation[i][1]/tile_length) + 1 )*tile_length;   //y right neighboor
            
            //left or lower line
            VectorXd O_IR_P_exact_tile_lower(2);
            O_IR_P_exact_tile_lower[0] = ( round(measurement_equation[i][0]/tile_length ) - 1)*tile_length;    // x lower neighboor
            O_IR_P_exact_tile_lower[1] = ( round(measurement_equation[i][1]/tile_length ) - 1)*tile_length;    // y left neighboor
            
            //C matrix of the Kalman filter coresponds to horizontal and verticle line
            C_x << 1, 0, -m_position_active_sensor[i][0]*sin(X[2]) - m_position_active_sensor[i][1]*cos(X[2]);
            C_y << 0, 1, -m_position_active_sensor[i][1]*sin(X[2]) + m_position_active_sensor[i][0]*cos(X[2]);
            
            //Mahalanobis distance calculations for x and y line
            double delta_mesure_X = O_IR_P_exact_tile[0] - measurement_equation[i][0];      //Yk - Yk_hat
            double dMaha_X = abs(delta_mesure_X)/sqrt((C_x*P*C_x.transpose())(0,0)+Qgamma);
            
            double delta_mesure_Y = O_IR_P_exact_tile[1] - measurement_equation[i][1];      //Yk - Yk_hat
            double dMaha_Y = abs(delta_mesure_Y)/sqrt((C_y*P*C_y.transpose())(0,0)+Qgamma);
            
            /* ********************************************************************** */
            // calculate the Mahalanobis for its neighboor x and y line
            double delta_mesure_right_higher_X = O_IR_P_exact_tile_higher[0] - measurement_equation[i][0];      //Yk - Yk_hat
            double dMaha_right_higher_X = abs(delta_mesure_right_higher_X)/sqrt((C_x*P*C_x.transpose())(0,0)+Qgamma);
            
            double delta_mesure_right_higher_Y = O_IR_P_exact_tile_higher[1] - measurement_equation[i][1];      //Yk - Yk_hat
            double dMaha_right_higher_Y = abs(delta_mesure_right_higher_Y)/sqrt((C_y*P*C_y.transpose())(0,0)+Qgamma);
            
            /* ****************** */
            double delta_mesure_left_lower_X = O_IR_P_exact_tile_lower[0] - measurement_equation[i][0];      //Yk - Yk_hat
            double dMaha_left_lower_X = abs(delta_mesure_left_lower_X)/sqrt((C_x*P*C_x.transpose())(0,0)+Qgamma);
            
            double delta_mesure_left_lower_Y = O_IR_P_exact_tile_lower[1] - measurement_equation[i][1];      //Yk - Yk_hat
            double dMaha_left_lower_Y = abs(delta_mesure_left_lower_Y)/sqrt((C_y*P*C_y.transpose())(0,0)+Qgamma);
            /* ******************************************************************************* */
            //vote the dots to plot
             if( abs(delta_mesure_X) < 0.01 || abs(delta_mesure_Y) < 0.01 )      //when the error smaller than 0.01m vote it to be correct one
            {
                if( abs(delta_mesure_X) < 0.01 && abs(delta_mesure_Y) < 0.01 )
                {
                    dMaha_green1_dots.push_back(dMaha_X);
                    dMaha_green2_dots.push_back(dMaha_Y);
                    
                    dMaha_red1_dots.push_back(-10);     //use -10 to take up the index(useless value)
                    dMaha_red2_dots.push_back(-10);
                    
                    dMaha_neighboor_red_dots1.push_back(dMaha_right_higher_X);
                    dMaha_neighboor_red_dots2.push_back(dMaha_right_higher_Y);
                    dMaha_neighboor_red_dots3.push_back(dMaha_left_lower_X);
                    dMaha_neighboor_red_dots4.push_back(dMaha_left_lower_Y);
                    
                    dMaha_dots_index.push_back(index);
                    index = index + 1;      //calculate index to plot
                }
                
                else if(abs(delta_mesure_X) < abs(delta_mesure_Y) )  // if the x line rounded error less than y line rounded error, then the most likely line actually detected is x line
                    //  we should vote the dM as green dots, then calculate its neighboor dM
                {
                    dMaha_green1_dots.push_back(dMaha_X);
                    dMaha_green2_dots.push_back(-10);
                    
                    dMaha_red1_dots.push_back(-10); // y
                    dMaha_red2_dots.push_back(-10);
                    
                    dMaha_neighboor_red_dots1.push_back(dMaha_right_higher_X);
                    dMaha_neighboor_red_dots2.push_back(-10);
                    dMaha_neighboor_red_dots3.push_back(dMaha_left_lower_X);
                    dMaha_neighboor_red_dots4.push_back(-10);
                    
                    dMaha_dots_index.push_back(index);
                    index = index + 1;
                }
                
                else // if the x line rounded error less than y line rounded error, then the most likely line actually detected is x line
                    //  we should vote the dM as green dots, then calculate its neighboor dM
                {
                    dMaha_green1_dots.push_back(dMaha_Y);
                    dMaha_green2_dots.push_back(-10);
                    
                    dMaha_red1_dots.push_back(-10);
                    dMaha_red2_dots.push_back(-10);
                    
                    dMaha_neighboor_red_dots1.push_back(-10);
                    dMaha_neighboor_red_dots2.push_back(dMaha_right_higher_Y);
                    dMaha_neighboor_red_dots3.push_back(-10);
                    dMaha_neighboor_red_dots4.push_back(dMaha_left_lower_Y);
                    
                    dMaha_dots_index.push_back(index);
                    index = index +1 ;
                }
            }
            else
            {
            
                dMaha_red1_dots.push_back(-10);
                dMaha_red2_dots.push_back(-10);
                
                dMaha_green1_dots.push_back(-10);
                dMaha_green2_dots.push_back(-10);
                
                dMaha_neighboor_red_dots1.push_back(dMaha_right_higher_X);
                dMaha_neighboor_red_dots2.push_back(dMaha_right_higher_Y);
                dMaha_neighboor_red_dots3.push_back(dMaha_left_lower_X);
                dMaha_neighboor_red_dots4.push_back(dMaha_left_lower_Y);
                
                dMaha_dots_index.push_back(index);
                index = index +1;
            }
            
            //Kalman wrt. horizantal line
            if (dMaha_X < mahaThreshold && dMaha_Y > mahaThreshold)
            {
                K = P*C_x.transpose()/((C_x*P*C_x.transpose())(0,0)+Qgamma);
                X = X + K*delta_mesure_X;
                P = (MatrixXd::Identity(3,3)-K*C_x)*P;
            }
            
            //Kalman wrt. vertical line
            else if (dMaha_Y < mahaThreshold && dMaha_X > mahaThreshold)
            {
                K = P*C_y.transpose()/((C_y*P*C_y.transpose())(0,0)+Qgamma);
                X = X + K*delta_mesure_Y;
                P = (MatrixXd::Identity(3,3)-K*C_y)*P; 
            }
            //for detect the corner
            else if(dMaha_X < mahaThreshold && dMaha_Y < mahaThreshold)
            {
                K = P*C_x.transpose()/((C_x*P*C_x.transpose())(0,0)+Qgamma);
                X = X + K*delta_mesure_X;
                P = (MatrixXd::Identity(3,3)-K*C_x)*P;
                
                K = P*C_y.transpose()/((C_y*P*C_y.transpose())(0,0)+Qgamma);
                X = X + K*delta_mesure_Y;
                P = (MatrixXd::Identity(3,3)-K*C_y)*P;
                
                /*dMaha_green1_dots.push_back(dMaha_X);
                dMaha_green2_dots.push_back(dMaha_Y);
                
                dMaha_red1_dots.push_back(-10);     //use -10 to take up the index(useless value)
                dMaha_red2_dots.push_back(-10);
                
                dMaha_neighboor_red_dots1.push_back(dMaha_right_higher_X);
                dMaha_neighboor_red_dots2.push_back(dMaha_right_higher_Y);
                dMaha_neighboor_red_dots3.push_back(dMaha_left_lower_X);
                dMaha_neighboor_red_dots4.push_back(dMaha_left_lower_Y);
                
                dMaha_dots_index.push_back(index);
                index = index + 1;      //calculate index to plot*/
            }
            
            /*else
            {
                dMaha_red1_dots.push_back(-10);
                dMaha_red2_dots.push_back(-10);
                
                dMaha_green1_dots.push_back(-10);
                dMaha_green2_dots.push_back(-10);
                
                dMaha_neighboor_red_dots1.push_back(dMaha_right_higher_X);
                dMaha_neighboor_red_dots2.push_back(dMaha_right_higher_Y);
                dMaha_neighboor_red_dots3.push_back(dMaha_left_lower_X);
                dMaha_neighboor_red_dots4.push_back(dMaha_left_lower_Y);
                
                dMaha_dots_index.push_back(index);
                index = index +1;
            }*/
        }
        
        measurement_equation.clear();		//empty the vector
        m_position_active_sensor.clear();
        
        //Publication of the state and the matrix of uncertainties (diagonal coefficients)
        //        geometry_msgs::Pose2D state;
        //        state.x = X[0];
        //        state.y = X[1];
        //        state.theta = X[2];
        //        publisher_state.publish(state);
        // plot the filted trajectory
        filted_x.push_back(X[0]);
        filted_y.push_back(X[1]);
        
        //        geometry_msgs::Pose2D variance;
        //        variance.x = P(0,0);
        //        variance.y = P(1,1);
        //        variance.theta = P(2,2);
        //        publisher_variance.publish(variance);
        //plot the varience of P
        thegema_x.push_back( sqrt(P(0,0)) );
        thegema_y.push_back( sqrt(P(1,1)) );
        thegema_theta.push_back( sqrt(P(2,2))*180/PI );
        
        loop_rate.sleep();
    }
    
    //use matplotcpp.h, which can check from https://github.com/lava/matplotlib-cpp
    plt::figure();
    plt::subplot(3,1,1);
    plt::plot(travel_distance, thegema_x);
    plt::ylabel("thegma_x");
    
    plt::subplot(3,1,2);
    plt::plot(travel_distance, thegema_y);
    plt::ylabel("thegma_y");
    
    plt::subplot(3,1,3);
    plt::plot(travel_distance, thegema_theta);
    plt::ylim(0.0,1.5);
    plt::ylabel("thegma_theta");
    plt::xlabel("travel_distance");
    //the absolute path to save the graph
    
    const char *extention = ".svg";
    const char *underScope = "_";
    
    char* path = new char[strlen(absolutepath)+strlen(folder) + strlen(bagName) + strlen(underScope) +1];
    *path = '\0';
    strcat(path,absolutepath);
    strcat(path,folder);
    strcat(path,bagName);
    strcat(path,underScope);
    
    const char *typePlot = "thegma";
    char* filename1 = new char[strlen(path)+strlen(typePlot) + strlen(comment) + strlen(extention) +1];
    *filename1 = '\0';
    strcat(filename1,path);
    strcat(filename1,typePlot);
    strcat(filename1,comment);
    strcat(filename1,extention);
    plt::save(filename1);
    
    //plot maha distance
    plt::figure();
    vector<double> threshold_line(dMaha_dots_index.size(), mahaThreshold);
    plt::plot(dMaha_dots_index, dMaha_green1_dots,"g.", dMaha_dots_index, dMaha_green2_dots,"g.",
              dMaha_dots_index, dMaha_red1_dots,"r.", dMaha_dots_index, dMaha_red2_dots,"r.", dMaha_dots_index, threshold_line,"b-",
              dMaha_dots_index, dMaha_neighboor_red_dots1, "r.", dMaha_dots_index, dMaha_neighboor_red_dots2, "r.",
              dMaha_dots_index, dMaha_neighboor_red_dots3, "r.", dMaha_dots_index, dMaha_neighboor_red_dots4, "r.");
    
    plt::ylim(0, 160);
    const char *typePlot2 = "Dmaha";
    char* filename2 = new char[strlen(path)+strlen(typePlot2) + strlen(comment) + strlen(extention) +1];
    *filename2 = '\0';
    strcat(filename2,path);
    strcat(filename2,typePlot2);
    strcat(filename2,comment);
    strcat(filename2,extention);
    plt::save(filename2);
    
    // plot trajectory
    plt::figure();
    plt::plot(filted_x, filted_y, "r-", odom_x, odom_y, "k-");
    
    const char *typePlot3 = "trajectory";
    char* filename3 = new char[strlen(path)+strlen(typePlot3) + strlen(comment) + strlen(extention) +1];
    *filename3 = '\0';
    strcat(filename3,path);
    strcat(filename3,typePlot3);
    strcat(filename3,comment);
    strcat(filename3,extention);
    plt::save(filename3);
    return 0;
}
