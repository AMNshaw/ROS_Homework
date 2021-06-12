#include "ukf.h"


ukf::ukf(int state_size , int measurement_size){
  x_size = state_size;
  y_size = measurement_size;
  alpha = 1e-3;
  kappa = 0;
  beta = 2;
  lambda = 0.0;
  C << 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0; 
  
  L=x_size;
  x_sigmavector_size = 2*x_size;////////////////////////////////////////////////

  lambda= alpha * alpha * (L + kappa) -L;

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

  w_c(0) = lambda/(lambda+L);//////////////////////////////////////////////////////
  w_m(0) = lambda/((L+lambda)+(1-pow(alpha,2)+beta));///////////////////////////////////////////////////

  for(int i=1 ; i<x_sigmavector_size ; i++){
    w_c(i) = 1/(2*(L+lambda));////////////////////////////////////////////////////////////
    w_m(i) = 1/(2*(L+lambda));//////////////////////////////////////////////////////////
  }

  // default Q R P matrix
  Q =5e-7*Eigen::MatrixXd::Identity(x_size, x_size);
  R =5e-4*Eigen::MatrixXd::Identity(y_size,y_size);
  P=1e-3*Eigen::MatrixXd::Identity(x_size, x_size);

  P_.setZero(x_size,x_size);
  P_yy.setZero(y_size,y_size);
  P_xy.setZero(x_size,y_size);
}

//time update
void ukf::predict(){

  //find sigma point
  P=(lambda+L)*P;
  Eigen::MatrixXd M= (P).llt().matrixL();
  Eigen::MatrixXd buffer;
  x_sigmavector.col(0) = x;

  for(int i=0;i<x_size;i++){
    Eigen::VectorXd sigma =(M.row(i)).transpose();
    x_sigmavector.col(i+1) = x_sigmavector.col(0) + sigma;//////////////////////////////////////////////////////////////
    x_sigmavector.col(i+x_size+1) = x_sigmavector.col(0) - sigma;/////////////////////////////////////////////////////////////
  }

  // process model
  x_sigmavector = dynamics(x_sigmavector);

  //x_hat (mean)
  x_hat.setZero(x_size);   //initialize x_hat

  for(int i=0;i<x_sigmavector_size;i++){
    x_hat += w_m(i)*x_sigmavector(i);////////////////////////////////////////////
  }
  //covariance
  P_.setZero(x_size,x_size);

  for(int i=0 ; i<x_sigmavector_size ;i++){
    P_ += w_c(i)*(x_sigmavector(i)-x_hat)*(x_sigmavector(i)-x_hat).transpose();
  }

  //add process noise covariance
  P_+= Q;

  // measurement model
  y_sigmavector = C.transpose()*x_sigmavector;//////////////////////////////////////////////////////


  //y_hat (mean)
  y_hat.setZero(y_size);

  for(int i=0;i< x_sigmavector_size;i++){
    y_hat += w_m(i)*y_sigmavector(i);//////////////////////////////////////////////////////
  }
}

//measurement update
void ukf::correct(Eigen::VectorXd measure){

    y=measure;

    P_yy.setZero(y_size,y_size);
    P_xy.setZero(x_size,y_size);

    for(int i=0;i<x_sigmavector_size;i++){
      Eigen::MatrixXd y_err;
      Eigen::MatrixXd y_err_t;
      y_err = y_sigmavector(i) - y_hat;//////////////////////////////////////////////////////////////
      y_err_t = y_err.transpose();

      P_yy += w_c(i)*(y_err*y_err_t);/////////////////////////////////////////////////////////
    }

    //add measurement noise covarinace
    P_yy +=R;

    for(int i=0;i<x_sigmavector_size;i++){
      Eigen::VectorXd y_err , x_err;
      y_err = y_sigmavector(i) - y_hat;//////////////////////////////////////////////
      x_err = x_sigmavector(i) - x_hat;/////////////////////////////////////////
      P_xy += w_c(i)*(x_err*y_err.transpose());///////////////////////////////////////////
    }

    Kalman_gain = ?;

    // correct states and covariance
    x = x_hat+Kalman_gain*(y-y_hat);////////////////////////////////
    P = P_+Kalman_gain*P_yy*Kalman_gain.transpose();//////////////////////////////////
}

Eigen::MatrixXd ukf::dynamics(Eigen::MatrixXd sigma_state){
  Eigen::MatrixXd a = sigma_state;
  return a;
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
