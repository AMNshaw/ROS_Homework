/*
  Author : Ben You
  Target : Mechanical 2-order system simulation
  Date : 03/25 2020
*/
#include "ros/ros.h"
#include "std_msgs/Float64.h"
#include <stdlib.h>
#include <iostream>
using namespace std;
/*The dynamic system can be represented as "M*y_doubledot + b*y_dot + k*y = u ";
* Control canonical form:
* [x1_dot] = [0       1][x1] + [0]u;
* [x2_dot] = [-k/M -b/M][x2] + [1/M]u;
* x1 = y;
* x2 = y_dot;
* y = [1 0][x1];
*          [x2];
*/

/*You may not use this struct! Depend on your demand*/
struct Dynamic_State {
	float last = 0;
	float current = 0;
};
int main(int argc, char **argv)
{
	string mes;
	ros::init(argc,argv,"dynamic");
	ros::NodeHandle nh;
	ros::Publisher state_pub = nh.advertise<std_msgs::Float64>("/Position",10);
	/*Parameters initialization*/
	float M = 2;//mass
	float k = 1;//spring constant
	float b = 0.707;//damping coefficient
	float u = 1;//control input is unit step function(r = u,without controller design)constant force,
	Dynamic_State x2;//velocity
	Dynamic_State x1;//position
	Dynamic_State x1_dot;
	Dynamic_State x2_dot;
	ros::Rate loop_rate(10);
	
	double dt;//sampling time
	std_msgs::Float64 position;
	/*Compute the position*/
	position.data = 0;
	state_pub.publish(position);
	cin >> mes;
	double past = ros::Time::now().toSec();
	while(ros::ok()) {
		double now = ros::Time::now().toSec();
		/*Discrete-Time Linear State-Space*/
		dt = now - past;
		
		/*Following codes are implemented with numerical integration
		 *x(k+1) = x(k) + x_dot * dt;
		 *
		 *Please implement your codes here.
		 */
		x2_dot.last = -k/M * x1.last - b/M * x2.last + u/M;
		x1_dot.last = x2.last;

		x2.current = x2.last + x2_dot.last * dt;
		x1.current = x1.last + x1_dot.last * dt;


		x2.last = x2.current;
		x1.last = x1.current;
		/*Current moment is last moment in the future*/
		past = now;
		/*You can check the position data by PlotJuggler*/
		position.data = x1.current; 
		cout << position.data << endl;
		state_pub.publish(position);
		ros::spinOnce();
		loop_rate.sleep();
	}
	return 0;
}

