#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  //Predict the state
  x_ = F_* x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_*P_*Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  //update the state by using Kalman Filter equations
  MatrixXd y = z - H_*x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd P_Ht = P_*Ht;
  MatrixXd S = H_*P_Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_Ht*Si;
  x_ += K*y;
  P_ = (MatrixXd::Identity(4,4) - K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  //update the state by using Extended Kalman Filter equations
  VectorXd z_pred=VectorXd(3);
  double px=x_(0);
  double py=x_(1);
  double vx=x_(2);
  double vy=x_(3);

  double roh_pred = sqrt(px*px + py*py);
  double psi_pred = atan(py/px);
  float dot_psi_pred = (px*vx+py*vy)/roh_pred;


  if (roh_pred < 0.01) {
    roh_pred=0.01;
    psi_pred = atan(py/roh_pred)*px/(fabs(px)+0.0001);
    dot_psi_pred = (px*vx+py*vy)/roh_pred;
  }
  z_pred << roh_pred, psi_pred, dot_psi_pred;
  //cout<<"z pred"<<z_pred<<endl;


  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  //cout<<"x="<<x_<<endl;
  //cout<<"y="<<y<<"\n"<<endl;
  //cout<<"K="<<K<<endl;
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}
