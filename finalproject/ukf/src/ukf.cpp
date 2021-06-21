#include "ukf.h"
#include <tf/transform_datatypes.h>

ukf::ukf(int state_size , int measurement_size){

  x_size = state_size;
  y_size = measurement_size;
  alpha = 1e-3;
  kappa = 0;
  beta = 2;
  lambda = -13;

  L = x_size;
  x_sigmavector_size = 2*x_size+1;////////////////////////////////////////////////

  x.setZero(x_size);
  y.setZero(y_size);

  x_hat.setZero(x_size);
  y_hat.setZero(y_size);

  x_a.setZero(x_size+x_size+y_size);

  x_sigmavector.setZero(x_size,x_sigmavector_size);
  y_sigmavector.setZero(x_sigmavector_size,y_size);

  H.setZero(y_size,x_size);  // measurement matrix

  y = H*x;

  w_c.setZero(x_sigmavector_size);
  w_m.setZero(x_sigmavector_size);

  w_m(0) = lambda/(lambda+L);///////////////////////////////////////////////////
  w_c(0) = lambda/(lambda+L)+(1-alpha*alpha+beta);//////////////////////////////////////////////////////
  
  for(int i=1 ; i < x_sigmavector_size ; i++){
    w_c(i) = 1/(2*(L+lambda));////////////////////////////////////////////////////////////
    w_m(i) = 1/(2*(L+lambda));//////////////////////////////////////////////////////////
  }
  // default Q R P matrix
  Q = 5e-7*Eigen::MatrixXd::Identity(x_size, x_size);
  R = 5e-4*Eigen::MatrixXd::Identity(y_size, y_size);
  P = 1e-3*Eigen::MatrixXd::Identity(x_size, x_size);

  P_.setZero(x_size,x_size);
  P_yy.setZero(y_size,y_size);
  P_xy.setZero(x_size,y_size);

  last_quat << 0,0,0,1;

  /*  for(int i=0; i<2; i++)
  {
    std::cout << "w_m(" << i <<") " << w_m(i) << std::endl;
    std::cout << "w_c(" << i <<") " << w_c(i) << std::endl;
  }*/
}

//time update
void ukf::predict(){
  //find sigma point
  P=(lambda+L)*P;
  Eigen::MatrixXd M = (P).llt().matrixL();

  x_sigmavector.setZero();
  x_sigmavector.col(0) = x;
  for(int i=0;i<x_size;i++){
    Eigen::VectorXd sigma = (M.row(i)).transpose();
    x_sigmavector.col(i+1) = x_sigmavector.col(0) + (M.row(i)).transpose();//////////////////////////////////////////////////////////////
    x_sigmavector.col(i+x_size+1) = x_sigmavector.col(0) - (M.row(i)).transpose();/////////////////////////////////////////////////////////////
  }


  std::cout << "x(9) before " << std::endl << x_sigmavector.col(9) << std::endl;
  //std::cout << "x(19) before " << std::endl << x_sigmavector.col(19) << std::endl;
  //std::cout << "x(38) before " << std::endl << x_sigmavector.col(38) << std::endl;
  // process model
  x_sigmavector = dynamics(x_sigmavector);/////////////////////////////////////////////

  /*for(int i=0; i<x_sigmavector_size; i++)
  {
    std::cout << "x(i) after " << i << std::endl << x_sigmavector.col(i) << std::endl;
  }*/

  //x_hat (mean)
  x_hat.setZero(x_size);   //initialize x_hat
  for(int i=0;i<x_sigmavector_size;i++){
    
    x_hat = x_hat + w_m(i)*x_sigmavector.col(i);/////////////////////////////////////////////////////////
    //std::cout << "x_hat_plus " << i << std::endl << x_hat << std::endl;
  }

  //covariance
  P_.setZero(x_size,x_size);

  for(int i=0 ; i<x_sigmavector_size ;i++){
    Eigen::MatrixXd x_err;
    x_err = x_sigmavector.col(i)-x_hat;
    P_ = P_ + w_c(i)*(x_err*x_err.transpose());////////////////////////
  }
  //add process noise covariance
  P_+= Q;
  // measurement model
  y_sigmavector = H*x_sigmavector;//////////////////////////////////////////////////////
  //y_hat (mean)
  y_hat.setZero(y_size);
  for(int i=0;i< x_sigmavector_size;i++){
    y_hat += w_m(i)*y_sigmavector.col(i);//////////////////////////////////////////////////////
  }

  std::cout << "x(9) after " << std::endl << x_sigmavector.col(9) << std::endl;
  //std::cout << "x(4) after " << std::endl << x_sigmavector.col(4) << std::endl;
  //std::cout << "x(23) after " << std::endl << x_sigmavector.col(23) << std::endl;
  //std::cout << "x_sigmavector.col(19) " << std::endl << x_sigmavector.col(19) << std::endl;
  //std::cout << "x_sigmavector.col(38) " << std::endl << x_sigmavector.col(38) << std::endl;
  std::cout << "x_hat1" << std::endl << x_hat << std::endl;
  //std::cout << "y_hat1" << std::endl << y_hat << std::endl;

}

//measurement update
void ukf::correct(Eigen::VectorXd measure){

  Eigen::Vector4d delta_q_m_k;
  Eigen::Vector3d delta_q_m_k3;
  double delta_q_m_k4;
  Eigen::Matrix4d phi_q_m_k;
  Eigen::Vector4d q_m_k_inv;
  Eigen::Vector3d e_m_k;
  double quat_m_value;
  Eigen::Vector4d delta_q_k;
  Eigen::Vector3d delta_q_kv;
  double delta_q_ks;
  double a;
  double f;
  double p_value;
  Eigen::Vector3d p;
  Eigen::Matrix4d phi_q_k;
  Eigen::Vector4d q_k1;

  a = 3;
  f = 2 * a +1;
  quat_m_value = quat_m(0) * quat_m(0) + quat_m(1) * quat_m(1) + quat_m(2) * quat_m(2) + quat_m(3) * quat_m(3);

  q_m_k_inv << -quat_m(0), -quat_m(1), -quat_m(2), quat_m(3);

  phi_q_m_k << last_quat(3), -last_quat(2), last_quat(1), last_quat(0),
               last_quat(2), last_quat(3), -last_quat(0), last_quat(1),
              -last_quat(1), last_quat(0), last_quat(3), last_quat(2),
              -last_quat(0), -last_quat(1), -last_quat(2), last_quat(3);

  delta_q_m_k = phi_q_m_k * q_m_k_inv;
  delta_q_m_k3 << delta_q_m_k(0), delta_q_m_k(1), delta_q_m_k(2);
  delta_q_m_k4 = delta_q_m_k(3);
  e_m_k = f*delta_q_m_k3/(a+delta_q_m_k4);
  measure[6] = e_m_k(0);
  measure[7] = e_m_k(1);
  measure[8] = e_m_k(2);

  //----------------
  y=measure;
  P_yy.setZero(y_size,y_size);
  P_xy.setZero(x_size,y_size);

  for(int i=0;i<x_sigmavector_size;i++){
    Eigen::MatrixXd y_err;
    y_err = y_sigmavector.col(i) - y_hat;//////////////////////////////////////////////////////////////
    P_yy += w_c(i)*(y_err*y_err.transpose());/////////////////////////////////////////////////////////
  }

  //add measurement noise covarinace
  P_yy +=R;
 //std::cout << "P_yy1" << std::endl << P_yy << std::endl;
  for(int i=0;i<x_sigmavector_size;i++){
    Eigen::VectorXd y_err , x_err;

    y_err = y_sigmavector.col(i) - y_hat;//////////////////////////////////////////////
    x_err = x_sigmavector.col(i) - x_hat;/////////////////////////////////////////
    P_xy += w_c(i)*(x_err*y_err.transpose());///////////////////////////////////////////
    
    //std::cout << "x_err" << std::endl << x_err << std::endl;
    //std::cout << "y_err" << std::endl << y_err << std::endl;
    //std::cout << "P_xy(i)" << std::endl << P_xy << std::endl;
  }
  
  Kalman_gain = P_xy * P_yy.inverse(); ///////////////////////////////
  //std::cout << "P_xy1" << std::endl << P_xy << std::endl;
  //std::cout << "Kal_1" << std::endl << Kalman_gain << std::endl;
  // correct states and covariance
  x = x_hat + Kalman_gain*(y-y_hat);////////////////////////////////

  
  P = P_- Kalman_gain*P_yy*Kalman_gain.transpose();//////////////////////////////////

  std::cout << "x correct" << std::endl << x << std::endl;
  //std::cout << "P correct" << std::endl << P << std::endl;
  //----------------
  p_value = x[6]*x[6]+x[7]*x[7]+x[8]*x[8];
  p << x[6], x[7], x[8];
  delta_q_ks = (-a*(p_value) + f*sqrt(f*f+(1-a*a)*(p_value)))/(f*f+p_value);
  delta_q_kv = (a + delta_q_ks)*p/f;
  delta_q_k << delta_q_kv, delta_q_ks;

  phi_q_k << qk11(3), -qk11(2), qk11(1), qk11(0), //x_hat[16], x_hat[17], x_hat[18], x_hat[19]
             qk11(2), qk11(3), -qk11(0), qk11(1), //qk11(0), qk11(1), qk11(2), qk11(3)
            -qk11(1), qk11(0), qk11(3), qk11(2),
            -qk11(0), -qk11(1), -qk11(2), qk11(3);
  q_k1 = phi_q_k*delta_q_k;

  last_quat = q_k1;

  tf::Quaternion quat_transform(q_k1(0),q_k1(1),q_k1(2),q_k1(3));
  double roll,pitch,yaw;
  tf::Matrix3x3(quat_transform).getRPY(roll,pitch,yaw);
  euler_angle << roll*(180/3.1415926),pitch*(180/3.1415926),yaw*(180/3.1415926);
  quaternion << q_k1(0),q_k1(1),q_k1(2),q_k1(3);
}

Eigen::MatrixXd ukf::dynamics(Eigen::MatrixXd sigma_state){
  Eigen::MatrixXd a = sigma_state;
  return a;
}

Eigen::MatrixXd ukf::rotate(double roll, double yaw, double pitch){

  Eigen::MatrixXd frame;
  frame.setZero(3,3);
  double c_r = cos(roll) , s_r = sin(roll);
  double c_p = cos(pitch) , s_p = sin(pitch);
  double c_y = cos(yaw) , s_y = sin(yaw);
  frame << c_p*c_y ,c_y*s_p*s_r-s_y*c_r , c_y*s_p*c_r+s_y*s_r,
           c_p*s_y ,s_y*s_p+c_r*c_y     , s_y*s_p*c_r-c_y*s_r,
           -1*s_p  ,s_r*c_p             , c_r*c_p            ;
  return frame;
}

void ukf::set_measurement_matrix(Eigen::MatrixXd matrix){
  H = matrix;
}
void ukf::set_process_noise(Eigen::MatrixXd matrix){
  this->Q = matrix;
}

void ukf::set_measurement_noise(Eigen::MatrixXd matrix){
  this->R = matrix;
}

ukf::~ukf(){
}

