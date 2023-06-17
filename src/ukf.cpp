#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // Filter settings
  time_us_ = 0;
  is_initialized_ = false;
  use_laser_ = true;
  use_radar_ = true;

  // Measuremnt uncertainity values
  std_a_ = 30;
  std_yawdd_ = 30;

  // Process uncertainity values
  std_laspx_ = 0.15;
  std_laspy_ = 0.15;
  std_radr_ = 0.3;
  std_radphi_ = 0.03;
  std_radrd_ = 0.3;
  
  // Filter variables
  n_x_ = 5; // Dimensions of the state
  n_aug_ = 7; // Dimensions of the augmented state
  lambda_ = 3 - n_aug_; // Spread for sigma points
  n_sigma_points_ = (2 * n_aug_) + 1; // Number of sigma points needed. (2 per dimension + the current mean)
  
  x_ = VectorXd::Zeros(5); // State 
  P_ = MatrixXd::Identity(5, 5); // Process covariance
  Xsig_pred_ = MatrixXd::Zeros(n_aug_, n_sigma_points_);

  weights_ = MatrixXd::Zeros(n_sigma_points_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  double temp_weight = 0.5 / ( lambda_ + n_aug_);
  for(size_t i=1; i < n_sigma_points_; i++)
  {
    weights_(i) = temp_weight;
  }
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
    const bool recieved_radar_data = (meas_package.sensor_type_ == MeasurementPackage::RADAR);
    const bool recieved_lidar_data = (meas_package.sensor_type_ == MeasurementPackage::LASER);
    
    if(is_initialized_)
    {
      
      // UKF is already initialized, do the predict-update cycle

      double dt = meas_package.timestamp_ - time_us_;
      time_us_ = meas_package.timestamp_;

      Predict(dt);

      if(use_laser_ && recieved_lidar_data)
      {
        UpdateLidar(meas_package);
      }
      else if(use_radar_ && recieved_radar_data)
      {
        UpdateRadar(meas_package);
      }

    }
    else
    {
      // UKF is not initialized yet
      
      if(recieved_lidar_data)
      {
        InitializeFromLidar(meas_package);
      }
      else if(recieved_radar_data)
      {
        InitializeFromRadar(meas_package);
      }

      // Initialize the time
      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
    }

}


void UKF::InitializeFromLidar(MeasurementPackage meas_package)
{
  x_(0) = meas_package.raw_measurements_(0); //px
  x_(1) = meas_package.raw_measurements_(1); //py

  P_(0, 0) = pow(std_laspx_, 2);
  P_(1, 1) = pow(std_laspy_, 2);
}


void UKF::InitializeFromRadar(MeasurementPackage meas_package)
{
  double range, turn_angle, range_rate;
  double px, py, radial_velocity, turn_rate;
  range = meas_package.raw_measurements_(0);
  turn_angle = measurement_package.raw_measurements_(1);
  range_rate = meas_package.raw_measurements_(2);

  px = range * cos(turn_angle);
  py = range * sin(turn_angle);

  x_(0) = px;
  x_(1) = py;
  x_(2) = radial_velocity;
  x_(3) = turn_angle;
  x_(4) = 0; // No direct conversion available for turn rate

  P_(2, 2) = pow(std_radrd_, 2);
  P_(3, 3) = pow(std_radphi_, 2);
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  // Augmented state mean and covariance.
  VectorXd::Zeros X_aug(n_aug_);
  X_aug.head(n_x_) = x_;
  MatrixXd::Zeros P_aug(n_aug_, n_aug_);
  P_.topLeftCorner(n_x_, n_x_) = P_;
  P_(n_x_, n_x_) = pow(std_a_, 2);
  P_(n_x_+1, n_x_+1) = pow(std_yawdd_, 2);



  // Compute the square root of augmented covariance and (lambda + n_aug) to be used in sigma point computation.
  MatrixXd Paug_sqrt = P_aug.llt().matrixL();
  double spread_factor = sqrt(lambda_ + n_aug_);
  MatrixXd sigma_factors = spread_factor * Paug_sqrt.colwise();

  // Compute sigma points
  // Xsig_pred_ is already initialized with zeros during construction.
  MatrixXd::Zeros Xsig_aug(n_aug_, n_sigma_points_);
  Xsig_aug(0) = X_aug;
  Xsig_aug.block(0, 1, n_aug_, n_aug_) = X_aug.replicate(1, n_aug) + sigma_factors;
  Xsig_aug.block(0, n_aug_ + 1, n_aug_, n_aug_) = X_aug.replicate(1, n_aug) - sigma_factors;

  
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}