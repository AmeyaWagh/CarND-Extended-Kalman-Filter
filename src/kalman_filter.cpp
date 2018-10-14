#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void
KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                   MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void
KalmanFilter::Predict()
{

   x_ = F_*x_;
   MatrixXd Ft = F_.transpose();
   P_ = F_ * P_ * Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z)
{

   // Calculate error
    Eigen::VectorXd y = z - ( H_ * x_ );

    // Update Kalman filter
    MatrixXd Ht = H_.transpose();
    MatrixXd S  = H_*P_*Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K  = P_ * Ht * Si;

    // estimate
    x_ = x_ + ( K * y );
    P_ -= K * H_ * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{

    double rho     = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
    double theta   = atan2(x_(1), x_(0));
    double rho_dot = (x_(0)*x_(2) + x_(1)*x_(3)) / rho;

    VectorXd h = VectorXd(3);
    h << rho, theta, rho_dot;

    VectorXd y = z - h;

    // normalize phi
    y(1) = atan2(sin(y(1)), cos(y(1)));

    // Update Kalman filter
    MatrixXd Ht = H_.transpose();
    MatrixXd S  = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K  = P_ * Ht * Si;

    //estimate
    x_ = x_ + ( K * y );
    P_ -= K * H_ * P_;

}
