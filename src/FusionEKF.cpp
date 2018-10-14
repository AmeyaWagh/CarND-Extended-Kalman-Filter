#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;
  noise_ax  =   9.0;
  noise_ay  =   9.0;

  // initializing matrices
  R_laser_  = MatrixXd(2, 2);
  R_radar_  = MatrixXd(3, 3);
  H_laser_  = MatrixXd(2, 4);
  Hj_       = MatrixXd(3, 4);
  x_        = VectorXd(4);
  Q_        = MatrixXd(4,4);
  P_        = MatrixXd(4,4);
  F_        = MatrixXd(4,4);

  //measurement covariance matrix
  //LASER
  R_laser_ << 0.0225,   0,
              0,        0.0225;

  // RADAR
  R_radar_ << 0.09,     0,      0,
              0,        0.0009, 0,
              0,        0,      0.09;


  // initializing states
  x_ << 0,0,0,0;

  // laser measurement matrix
  H_laser_ << 1,0,0,0,
              0,1,0,0;

  // transition matrix
  F_ << 1,  0,  1,  0,
        0,  1,  0,  1,
        0,  0,  1,  0,
        0,  0,  0,  1;

  P_ << 1,  0,    0,      0,
        0,  1,    0,      0,
        0,  0,    1000,   0,
        0,  0,    0,      1000;

  ekf_.Init(x_,P_,F_,H_laser_,R_laser_,Q_);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_)
  {

    cout << "EKF: " << endl;
    // first measurement    
//    ekf_.x_ = VectorXd(4);
//    ekf_.x_ << 1, 1, 1, 1;
    ekf_.P_ = P_;

    // initialize previous time
    if(     measurement_pack.sensor_type_ == MeasurementPackage::RADAR ||
            measurement_pack.sensor_type_ == MeasurementPackage::LASER)
      {
          previous_timestamp_ = measurement_pack.timestamp_;
          // done initializing, no need to predict or update
          is_initialized_ = true;
      }

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
      {
          //Convert radar from polar to cartesian coordinates and initialize state.
          polarToCartesian(measurement_pack);
      }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
      {

          // double x = measurement_pack.raw_measurements_[0];
          // double y = measurement_pack.raw_measurements_[1];
          // double v_x = 0;
          // double v_y = 0;
          // ekf_.x_ << x,y,v_x,v_y;
          ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0; 
      }

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  updateTime(measurement_pack);

  updateStateTransitionMat();

  updateProcessCovarianceMat();

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/


  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar updates
      ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
  {
    // Laser updates
      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

void
FusionEKF::updateTime(const MeasurementPackage &measurement_pack)
{
    dt = measurement_pack.timestamp_ - previous_timestamp_;
    dt /= 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
}

void
FusionEKF::updateStateTransitionMat()
{
    ekf_.F_ = Eigen::MatrixXd(4,4);
    ekf_.F_ <<    1,0,dt,0,
                  0,1,0, dt,
                  0,0,1, 0,
                  0,0,0, 1;
}

void
FusionEKF::updateProcessCovarianceMat()
{

    float dt_2 = dt * dt;       //  dt^2
    float dt_3 = dt_2 * dt;     //  dt^3
    float dt_4 = dt_3 * dt;     //  dt^4

    float dt_4by4 = dt_4 / 4;    //  dt^4/4
    float dt_3by2 = dt_3 / 2;    //  dt^3/2

    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt_4by4 * noise_ax,   0,                  dt_3by2 * noise_ax,  0,
               0,                   dt_4by4 * noise_ay,  0,                  dt_3by2 * noise_ay,
               dt_3by2 * noise_ax,   0,                  dt_2 * noise_ax,    0,
               0,                   dt_3by2 * noise_ay,  0,                  dt_2 * noise_ay;
}

void
FusionEKF::polarToCartesian(const MeasurementPackage &measurement_pack)
{
    double rho       = measurement_pack.raw_measurements_[0];
    double phi       = measurement_pack.raw_measurements_[1];
    double rho_dot   = measurement_pack.raw_measurements_[2];

    double x    = rho*cos(phi);
    double y    = rho*sin(phi);
    double v_x  = rho_dot*cos(phi);
    double v_y  = rho_dot*sin(phi);
    ekf_.x_ << x,y,v_x,v_y;
}
