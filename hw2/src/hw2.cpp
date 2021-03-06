// include ros library
#include "ros/ros.h"

// include msg library
#include <turtlesim/Pose.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/Point.h>
#include <time.h>
#include <cstdio>
#include <unistd.h>
#include <termios.h>
#include <fcntl.h>
#include <stdio.h>
// include math 
#include <math.h>
#define PI = 3.1415926;
using namespace std;

turtlesim::Pose pose;
geometry_msgs::Twist vel_msg;
geometry_msgs::Point goal_point;
//Define a data structure to 3D
struct XYZ{
  float x;
  float y;
  float z;
  float abs;
  float theta;
};
//Declare a variable.Its name is pos_err with XYZ data type
struct XYZ pos_err;

// declare call back function(call back the pose of robot)
void pos_cb(const turtlesim::Pose::ConstPtr& msg)
{
  pose = *msg;
} 


void transfer(float &x, float &y, float theta )
{
  float tempx = x;
  float tempy = y;

  x = cos(theta) * tempx - sin(theta) * tempy ;
  y = sin(theta) * tempx + cos(theta) * tempy ;
}

void control(double theta_trans)
{
   if (pos_err.abs > 2)
    {
      pos_err.abs = 2;
    }
    
    vel_msg.linear.x = pos_err.abs;

    vel_msg.angular.z = theta_trans;
    cout << "vel.z : " <<vel_msg.angular.z<<endl;
    
    
    
}


int main(int argc, char **argv)
{
  ros::init(argc, argv, "tutorial_1");
  ros::NodeHandle n;

  // declare publisher & subscriber
  ros::Subscriber pos_sub = n.subscribe<turtlesim::Pose>("turtle1/pose", 10, pos_cb);
  ros::Publisher turtlesim_pub = n.advertise<geometry_msgs::Twist>("/turtle1/cmd_vel", 100);
  //input your desired position
  ROS_INFO("Please input (x,y). x>0,y>0");
  cout<<"desired_X:";
  cin>>goal_point.x;
  cout<<"desired_Y:";
  cin>>goal_point.y;
  // setting frequency as 10 Hz
  ros::Rate loop_rate(10);

  int count = 0;
  while (ros::ok()){

    printf("\ncount : %d\n",count);
    printf("goal x : %f \t y : %f\n",goal_point.x,goal_point.y);
    printf("pose x : %f \t y : %f\n",pose.x,pose.y);

    // Calculate position error(feedback term)
    pos_err.x = goal_point.x - pose.x;
    pos_err.y = goal_point.y - pose.y;
    pos_err.abs = sqrt(pow(pos_err.x, 2) + pow(pos_err.y, 2));
    transfer(pos_err.x, pos_err.y, -pose.theta);
    float theta_trans = atan2(pos_err.y,pos_err.x);
    cout << "theta_trans: " << theta_trans <<endl; 

    control(theta_trans);
    turtlesim_pub.publish(vel_msg);

    
    

    

    /*Your error-driven controller design
    /*vel_msg.linear.x = (Your control input design);
     *vel_msg.angular.z = (Your control input design);*/

   

    count ++;
    ros::spinOnce();
    loop_rate.sleep();

  }
  return 0;
}



