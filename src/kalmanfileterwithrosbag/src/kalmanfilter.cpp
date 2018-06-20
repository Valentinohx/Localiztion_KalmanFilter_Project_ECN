/*
    kalmanfilter.cpp
    created at 2018-05-15
    Xiao SHI && Huaxin LIU
    objective: use the robot odometry  and the sensor as input, develope an hybrid localazation system to
               localize the robot in the environment that with tiles floor, the main algorithm used is EKF
               
    Subscribe : 1. /capteur_carrelage, in which the sensors output has been transformed into binary
                2. /mobile_base/sensors/core , in which the robot's odometry is read from its internal sensor,
                    more detailed, its wheels' encoder
    Publish   : 1. /odom_carrelage, in which publised the robots' odometry information, include its x,y,theta
                2. /position_capteur, which is the output of the filter, include its x, y and theta
                3. /variance_carrelage, which is the filter's variance that described the noises of the system
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

//#include </home/valentinohx/ros/src/kalmanfileterwithrosbag/src/matplotlibcpp.h>
#include </home/xiao/ros_lab_ws/src/kalman_new/src/kalmanfileterwithrosbag/src/matplotlibcpp.h>


using namespace std;
using namespace Eigen;
namespace plt = matplotlibcpp;

const double PI = 3.14159265359;

ros::Publisher publisher_state;
ros::Publisher publisher_variance;
ros::Publisher publisher_capteur;


//Variables principales
vector<VectorXd> sensor_mesurements;			//mesures estimees des capteurs, sensor measurement
VectorXd X;                         //état  , with dynamic size, type is double, column vector
VectorXd U;                         //déplacement sous les roues
VectorXd X1;                         //état  , with dynamic size, type is double, column vector
VectorXd U1;                         //déplacement sous les roues
VectorXi Encoders;                  //encodage des positions des roues, interger with dynamic size
VectorXd deltaq;
//Matrices de Kalman
MatrixXd A;                         //matrice A de Kalman, with double type, and dynamic matrix
MatrixXd B;                         //matrice B de Kalman
MatrixXd C_x;                       //matrice C de Kalman du X
MatrixXd C_y;                       //matrice C de Kalman du Y
MatrixXd P;                     	//propagation error
MatrixXd K;
MatrixXd jointToCartesian ;


//Matrices de variances
MatrixXd Qalpha;
MatrixXd Qbeta;
double  Qgamma;

//Seuil de Mahalanobis
// use chi2inv to calculate it
double mahaThreshold = 2.7;
//Variables géométriques et temporelles
//double T = 0.02;							//période de mise à jour de l'odométrie
double T = 0.1;							//période de mise à jour de l'odométrie
double t0;									//temps système du début du programme
double encoder_coef = 2578.33;				//facteur de conversion donnée codeur -> angle (en incrément/tour de roue)
double lr = 0.227*1.103;					//longueur entre les roues, track gauch
double rr = 0.035;							//rayon des roues, the radius of the wheel
double trackGauge = 0.227*1.103;
double t_carreau = 0.30;					//taille des carreaux (bat. D : 0.1017, bat. E : 0.30)

//sensor number and position definition
double mys = 0.075;							//distance centre-capteur selon y, in robot fram
double mxs = 0;                             //distance centre-capteur selon x, in robot fram

vector< vector<double> > m_sensor_position = {{mxs, mys}, {mxs,-mys}};  //order should be left to right strictly
vector< vector<double> > m_sensor_position_filted;

vector<double> filted_x;
vector<double> filted_y;
vector<double> odom_x;
vector<double> odom_y;

vector< double> dMaha_green_dots;

vector< double> dMaha_red1_dots;
vector< double> dMaha_red2_dots;

 vector<double> dMaha_neighboor_red_dots1;
  vector<double> dMaha_neighboor_red_dots2;
   vector<double> dMaha_neighboor_red_dots3;
    vector<double> dMaha_neighboor_red_dots4;

vector< int> dMaha_dots_index;

int update_counter = 0;
int sensor_detect_counter = 0;

int sensor_number = m_sensor_position.size();

inline VectorXd EvolutionModelCore( VectorXd &X, const VectorXd &U )
{
    X[0] = X[0] + U[0] * cos(X[2]);
    X[1] = X[1] + U[0] * sin(X[2]);
    X[2] = X[2] + U[1];
    return X;
}

//---------------------------------------- MODELE D'EVOLUTION ----------------------------------------//
void evolutionModel()
{
    //store the last state of the encoder to get the change of angle, and avoid to miss every increment
    static int left_encoder_prec;    //static file scope: entire file, and it remains in memory until the program ends
    static int right_encoder_prec;
    //the first value of the encoder is abnormal, so we do not take it consideration
    //This modification imposes a pause at the beginning of the robot's tajectory in the command
    if(ros::Time::now().toSec()-t0>1.5)
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
        //        ROS_INFO("deltaq is: /n ");
        //        cout<<deltaq<<endl;
        
        U = jointToCartesian * deltaq;       //from localazation.pdf equation (5.19)
        U1 = jointToCartesian * deltaq;
        X = EvolutionModelCore( X, U );  //prediction phase
        X1 = EvolutionModelCore(X1,U1);  //pure odometry
        
    }
    //update the encoder
    left_encoder_prec  = Encoders[0];
    right_encoder_prec = Encoders[1];
}

//-------------------------------------- CALLBACKS MESURE ------------------------------------------//

// Mesure issue du capteur
void mesureCallback(std_msgs::Int8MultiArray etat_capteur)
{
    VectorXd C(2);
    for(int i = 0; i < sensor_number; i++ )
    {
		
        //Test s'il s'agit d'un joint de carrelage
        if ( etat_capteur.data[i])// if there is a new line detected, then update the measurements
        {
            //position estimée du capteur, we do not know which kind of line is detected, so calculate for both, then use Mahalanobis distance to tell
            C[0] = X[0] + m_sensor_position[i][0]*cos(X[2]) - m_sensor_position[i][1]*sin(X[2]); // corresponds to the current line in world fram, estimated X line  
            C[1] = X[1] + m_sensor_position[i][1]*cos(X[2]) + m_sensor_position[i][0]*sin(X[2]); // estimated Y line cordinates in the world fram
             
            sensor_detect_counter ++;
            sensor_mesurements.push_back(C);
            m_sensor_position_filted.push_back({m_sensor_position[i][0], m_sensor_position[i][1]});  //just keep the coordinates of the actived sensor
        }
    }
    cout<<"sensor_detect_counter:  " << sensor_detect_counter<<endl;
}
// Mise à jour les données brutes des encodeurs
void vitesseCallback(kobuki_msgs::SensorState etat_capteur)
{
    Encoders[0] = etat_capteur.left_encoder;
    Encoders[1] = etat_capteur.right_encoder;
}
//-------------------------------------- FONCTION MAIN ------------------------------------------//
int main(int argc, char **argv)
{
    //Initialisation des variables
    X = VectorXd::Zero(3);
    U = VectorXd(2);
    X1 = VectorXd::Zero(3);
    U1 = VectorXd(2);
    Encoders = VectorXi(2);
    A = MatrixXd(3,3);
    B = MatrixXd(3,2);
    C_x = MatrixXd(1,3);
    C_y = MatrixXd(1,3);
    P = MatrixXd(3,3);
    sensor_mesurements = vector<VectorXd>();
    
    jointToCartesian = MatrixXd(2,2);
    jointToCartesian <<         rr/2,             rr/2,
                           rr/trackGauge,   -rr/trackGauge;
    //Pose of depart, modify it based on real trajectory
    X << 0, 0, 45*PI/180;
    X1 << 0, 0, 45*PI/180;
    
    //Incertitude sur la position initiale
    P << pow(0.007,2), 0, 0,
            0, pow(0.007,2), 0,
            0, 0, pow(5*PI/180,2);
    
    //Matrices de variance
    Qalpha = MatrixXd(3,3);				//state noises 
    Qalpha <<  pow(0.003,2)*T, 0, 0,
               0, pow(0.003,2)*T, 0,
               0, 0, pow(0.2*PI/180,2)*T;
            
    Qbeta = MatrixXd(2,2);				//input noises 
    Qbeta << pow(0.00006,2), 0,
             0, pow(0.00006,2);
            
            
    Qgamma = pow(0.005,2)/12;			//measurment noises
    
    //Initialisation ROS
    ros::init(argc, argv, "conversion_capteur");
    ROS_INFO("Node conversion_capteur_node connected to roscore");
    ros::NodeHandle nh;
    
    //Subscribing
    ROS_INFO("Subscribing to topics\n");
    ros::Subscriber subscriber_capteur = nh.subscribe<std_msgs::Int8MultiArray>("/capteur_carrelage", 1, mesureCallback);
    ros::Subscriber subscriber_roue = nh.subscribe<kobuki_msgs::SensorState>("/mobile_base/sensors/core" , 1, vitesseCallback);
    
    //Publishing
    publisher_state = nh.advertise<geometry_msgs::Pose2D>("/odom_carrelage", 1);
    publisher_variance = nh.advertise<geometry_msgs::Pose2D>("/variance_carrelage", 1);
    publisher_capteur = nh.advertise<geometry_msgs::Pose2D>("/position_capteur", 3);
    
    //frequency
    ros::Rate loop_rate(1/T);
    
    //Temps initial
    t0 = ros::Time::now().toSec();
    int index = 0;
    while (ros::ok())
    {
        ros::spinOnce();
        //odometry
        evolutionModel();  
        odom_x.push_back(X1[0]);  //store the pure odometry 
        odom_y.push_back(X1[1]);
        
        //calculate the Jacobian matrix of the system
        A << 1, 0, -U[0]*sin(X[2]),
                0, 1, U[0]*cos(X[2]),
                0, 0, 1;
        B << cos(X[2]), 0,
                sin(X[2]), 0,
                0,     1;
        //---------------------------------------------- FILTRE DE KALMAN ------------------------------------------//
        //propagation error
        P = A*P*A.transpose()+B*Qbeta*B.transpose()+Qalpha;
        /*When a measurement is taken by one of the sensors, two distances of Mahalanobis are calculated:
          one considering that the line detected is horizontal, the other vertical
          if only one is below the threshold, it is retained, one uses the corresponding matrix C */
               
        for (int i = 0; i < sensor_mesurements.size(); i++) //sensor_mesurements is used to store the output of measurment equation
        {
            /*transfer the sensor measurments to exactly in the line, either a horizontal line or a vertical line, the c_rounded should be the 
              times of the length of the tiles, this is true based on that the robot always starts from the limits of a tile*/
            VectorXd C_round(2);
            C_round[0] = round(sensor_mesurements[i][0]/t_carreau)*t_carreau;		//Sensor position rounded to a tile according to x
            C_round[1] = round(sensor_mesurements[i][1]/t_carreau)*t_carreau;		//position du capteur arrondie à un carreau selon y
            
            //get the neighboors of the rounded line, this four more lines are used to calculate the  Mahalanobis distances
            //right or above line
            VectorXd C_round_right_higher(2);
            C_round_right_higher[0] = ( round(sensor_mesurements[i][0]/t_carreau) + 1 )*t_carreau;
            C_round_right_higher[1] = ( round(sensor_mesurements[i][1]/t_carreau) + 1 )*t_carreau;
            
            //left or lower line
            VectorXd C_round_left_lower(2);
            C_round_left_lower[0] = ( round(sensor_mesurements[i][0]/t_carreau ) - 1)*t_carreau;
            C_round_left_lower[1] = ( round(sensor_mesurements[i][1]/t_carreau ) - 1)*t_carreau;
           
           
           // cout<<"c round" << C_round[0] << "  "<<C_round[1] <<endl;
            //C matrix of the Kalman filter coresponds to horizontal and verticle line
            C_x << 1, 0, -m_sensor_position_filted[i][0]*sin(X[2]) - m_sensor_position_filted[i][1]*cos(X[2]);// maybe un peu de problem
            C_y << 0, 1, -m_sensor_position_filted[i][1]*sin(X[2]) + m_sensor_position_filted[i][0]*cos(X[2]);
            
            //Calculs des distances de Mahalanobis for x and y line
            double delta_mesure_X = C_round_right_higher[0] - sensor_mesurements[i][0];      //Yk - Yk_hat
            double dMaha_X = abs(delta_mesure_X)/sqrt((C_x*P*C_x.transpose())(0,0)+Qgamma);
           
            double delta_mesure_Y = C_round[1] - sensor_mesurements[i][1];  
            double dMaha_Y = abs(delta_mesure_Y)/sqrt((C_y*P*C_y.transpose())(0,0)+Qgamma);
           
           
           
           /* ********************************************************************** */
            // calculate the Mahalanobis for its neighboor x and y line 
            double delta_mesure_right_higher_X = C_round_right_higher[0] - sensor_mesurements[i][0];      //Yk - Yk_hat
            double dMaha_right_higher_X = abs(delta_mesure_right_higher_X)/sqrt((C_x*P*C_x.transpose())(0,0)+Qgamma);
            
            double delta_mesure_right_higher_Y = C_round_right_higher[1] - sensor_mesurements[i][0];      //Yk - Yk_hat
            double dMaha_right_higher_Y = abs(delta_mesure_right_higher_Y)/sqrt((C_x*P*C_x.transpose())(0,0)+Qgamma);
           
           /* ****************** */
            double delta_mesure_left_lower_X = C_round_left_lower[0] - sensor_mesurements[i][0];      //Yk - Yk_hat
            double dMaha_left_lower_X = abs(delta_mesure_right_higher_X)/sqrt((C_x*P*C_x.transpose())(0,0)+Qgamma);
            
            double delta_mesure_left_lower_Y = C_round_left_lower[1] - sensor_mesurements[i][0];      //Yk - Yk_hat
            double dMaha_left_lower_Y = abs(delta_mesure_left_lower_Y)/sqrt((C_x*P*C_x.transpose())(0,0)+Qgamma);
            /* ******************************************************************************* */
             
             //store the four irelavent Mahalanobis for the neighboor line
            
             
             dMaha_neighboor_red_dots1.push_back(dMaha_right_higher_X);
             dMaha_neighboor_red_dots2.push_back(dMaha_right_higher_Y);
             dMaha_neighboor_red_dots3.push_back(dMaha_left_lower_X);
             dMaha_neighboor_red_dots4.push_back(dMaha_left_lower_Y);
            
            
           // cout << delta_mesure_X <<"  " << delta_mesure_Y << endl;
            
            if( abs(delta_mesure_X) < 0.01 || abs(delta_mesure_Y) < 0.01 )
            {
                if(abs(delta_mesure_X) < abs(delta_mesure_Y) )  // if the x line rounded error less than y line rounded error, then the most likely line actually detected is x line
                    //  we should vote the dM as green dots, then calculate its neighboor dM
                {
                    dMaha_green_dots.push_back(dMaha_X);
                    dMaha_red1_dots.push_back(dMaha_Y);
					dMaha_red2_dots.push_back(-10);
					 
                    dMaha_dots_index.push_back(index);
                    index = index + 1;
                }
                /*else if( abs(delta_mesure_X) == abs(delta_mesure_Y))
                {
					
				}*/
                
                else // if the x line rounded error less than y line rounded error, then the most likely line actually detected is x line
                    //  we should vote the dM as green dots, then calculate its neighboor dM
                {
                    dMaha_green_dots.push_back(dMaha_Y);
                    dMaha_red1_dots.push_back(dMaha_X);
                    dMaha_red2_dots.push_back(-10);
                    
                    dMaha_dots_index.push_back(index);
                    index = index +1 ;
                }
            }
            else
            {
				
                dMaha_red1_dots.push_back(dMaha_X);
                dMaha_red2_dots.push_back(dMaha_Y);
                dMaha_green_dots.push_back(-10);
                
                dMaha_dots_index.push_back(index);
                index = index +1;
            }
            
               
            //Publication des positions estimées du capteur (pour l'analyse de résultats)
            //Le champ theta renseigne sur la distance de Mahalanobis calculée pour les deux types de ligne
            geometry_msgs::Pose2D capteur;
            capteur.x = sensor_mesurements[i][0];
            capteur.y = sensor_mesurements[i][1];
            capteur.theta = dMaha_X;  // what happens here?
            publisher_capteur.publish(capteur);
            
            capteur.x = sensor_mesurements[i][0];
            capteur.y = sensor_mesurements[i][1];
            capteur.theta = 1000+dMaha_Y;
            publisher_capteur.publish(capteur);
            
            
            //Kalman du X
            if (dMaha_X < mahaThreshold && dMaha_Y > mahaThreshold)
            {
                ROS_INFO("X update");
                K = P*C_x.transpose()/((C_x*P*C_x.transpose())(0,0)+Qgamma);
                X = X + K*delta_mesure_X;
                P = (MatrixXd::Identity(3,3)-K*C_x)*P;
                update_counter ++;
            }
            
            //Kalman du Y
            else if (dMaha_Y < mahaThreshold && dMaha_X > mahaThreshold)
            {
                ROS_INFO("Y update");
                K = P*C_y.transpose()/((C_y*P*C_y.transpose())(0,0)+Qgamma);
                X = X + K*delta_mesure_Y;
                P = (MatrixXd::Identity(3,3)-K*C_y)*P;
                update_counter ++;
            }
            cout<<"update_counter:  "<<update_counter<<endl;
            
        }
        sensor_mesurements.clear();		//penser à vider le vecteur de mesure
        m_sensor_position_filted.clear();
        
        //Publication de l'état et de la matrice d'incertitudes (coefficients diagonaux)
        geometry_msgs::Pose2D state;
        state.x = X[0];
        state.y = X[1];
        state.theta = X[2];
        publisher_state.publish(state);
        
        filted_x.push_back(X[0]);
        filted_y.push_back(X[1]);
        
        geometry_msgs::Pose2D variance;
        variance.x = P(0,0);
        variance.y = P(1,1);
        variance.theta = P(2,2);
        publisher_variance.publish(variance);
        
        /*plt::plot(filted_x,filted_y,"r-",odom_x,odom_y,"k-");
        plt::xlim(-1.5, 1.5);
        plt::ylim(-2.0, 2.5);
        plt::pause(0.01);*/
          vector<double> w(dMaha_dots_index.size(), mahaThreshold);
          
          /*plt::plot(dMaha_dots_index, dMaha_green_dots,"g*", dMaha_dots_index, dMaha_red1_dots,"r.",  dMaha_dots_index, dMaha_red2_dots,"r.", dMaha_dots_index, w,"r-",
          dMaha_dots_index, dMaha_neighboor_red_dots1, "r.", dMaha_dots_index, dMaha_neighboor_red_dots2, "r.",
          dMaha_dots_index, dMaha_neighboor_red_dots3, "r.", dMaha_dots_index, dMaha_neighboor_red_dots4, "r.");
          plt::ylim(0, 160);
          plt::pause(0.01);*/
          
          /*plt::plot(dMaha_dots_index, dMaha_green_dots,"g*", dMaha_dots_index, dMaha_red1_dots,"r.",  dMaha_dots_index, dMaha_red2_dots,"r.", dMaha_dots_index, w,"r-",
          dMaha_dots_index, dMaha_neighboor_red_dots1, "r.",  dMaha_dots_index, dMaha_neighboor_red_dots4, "r.");
          plt::pause(0.01);*/
          
          plt::plot(filted_x,filted_y,"r-",odom_x,odom_y,"k-");
          plt::pause(0.01);
          
        
    //    const char* filename = "/home/plot_test.png";
    //    plt::save(filename);*/
             
      //  ROS_INFO("1");
        loop_rate.sleep();
    }
       
      /* plt::plot(filted_x,filted_y,"r-",odom_x,odom_y,"k-");
       //plt::pause(0.01);
       for(int i = 0; i < 100000000000000 ; )*/
      //  string filename = "/home/xiao/plot_test.png";
    /*const char* filename = "/home/xiao/plot_test.png";
    plt::save(filename);
    ROS_INFO("store done");
    //    plt::show();
    //    while(true);*/
    return 0;
}
