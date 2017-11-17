#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Tools tools;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices

  //measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
	  0, 0.0225;


  //measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
	  0, 0.0009, 0,
	  0, 0, 0.09;



  //laser measurement matrix
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
	  0, 1, 0, 0;

  //initial transition matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
	  0, 1, 0, 1,
	  0, 0, 1, 0,
	  0, 0, 0, 1;

  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
	  0, 1, 0, 0,
	  0, 0, 1000, 0,
	  0, 0, 0, 1000;

  noise_ax = 9;
  noise_ay = 9;

} 

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    // first measurement
    cout << "EKF: " << endl;

	//state vector
	ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1; //TODO: play around with second two values to help RMSE

	float px;
	float py;
	float vx = 0;
	float vy = 0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
		float ro = measurement_pack.raw_measurements_[0];
		float theta = measurement_pack.raw_measurements_[1];
		px = ro * cos(theta);
		py = ro * sin(theta);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
		px = measurement_pack.raw_measurements_[0]; 
		py = measurement_pack.raw_measurements_[1];
    }
	// Special case initialisation problem (track 2)
	float eps = 0.000001;
	if (fabs(ekf_.x_(0)) < eps && fabs(ekf_.x_(1)) < eps) {
		ekf_.x_(0) = eps;
		ekf_.x_(1) = eps;
	}

	ekf_.x_ << px, py, vx, vy;
    // done initializing, no need to predict or update
	previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/



  // find difference in time between this measurement and the last in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  float dt2 = dt * dt;
  float dt3 = dt2 * dt;
  float dt4 = dt3 * dt;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << (dt4 / 4.0) * noise_ax,	0,					(dt3 / 2)*noise_ax,	0,
				0,						(dt4 / 4)*noise_ay, 0,					(dt3 / 2)*noise_ay,
				(dt3 / 2)*noise_ax,		0,					dt2*noise_ax,		0,
				0,						(dt3 / 2)*noise_ay, 0,					dt2*noise_ay;



  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/



  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  ekf_.R_ = R_radar_;
	  ekf_.H_ = tools.CalculateJacobian(ekf_.x_);

	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	 
  } else {
    // Laser updates
	  ekf_.R_ = R_laser_;
	  ekf_.H_ = H_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
