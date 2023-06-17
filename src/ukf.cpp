#include "ukf.h"


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
  n_obs_radar_ = 3;
  n_obs_lidar_ = 2;
  n_x_ = 5; // Dimensions of the state
  n_aug_ = 7; // Dimensions of the augmented state
  lambda_ = 3 - n_aug_; // Spread for sigma points
  n_sigma_points_ = (2 * n_aug_) + 1; // Number of sigma points needed. (2 per dimension + the current mean)
  
  x_ = VectorXd::Zero(5); // State 
  P_ = MatrixXd::Identity(5, 5); // Process covariance
  Xsig_pred_ = MatrixXd::Zero(n_x_, n_sigma_points_);

  weights_ = MatrixXd::Zero(n_sigma_points_);
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

      Prediction(dt);

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
  turn_angle = meas_package.raw_measurements_(1);
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

void UKF::Prediction(double delta_t) 
{
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  // Augmented state mean and covariance.
  VectorXd X_aug = VectorXd::Zero(n_aug_);
  X_aug.head(n_x_) = x_;
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_.topLeftCorner(n_x_, n_x_) = P_;
  P_(n_x_, n_x_) = pow(std_a_, 2);
  P_(n_x_+1, n_x_+1) = pow(std_yawdd_, 2);

  // Compute the square root of augmented covariance and (lambda_ + n_aug_) to be used in sigma point computation.
  MatrixXd Paug_sqrt = P_aug.llt().matrixL();
  double spread_factor = sqrt(lambda_ + n_aug_);
  MatrixXd sigma_factors = spread_factor * Paug_sqrt;

  // Compute sigma points
  // Xsig_pred_ is already initialized with Zero during construction.
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, n_sigma_points_);
  Xsig_aug.col(0) = X_aug;
  Xsig_aug.block(0, 1, n_aug_, n_aug_) = X_aug.replicate(1, n_aug_) + sigma_factors;
  Xsig_aug.block(0, n_aug_ + 1, n_aug_, n_aug_) = X_aug.replicate(1, n_aug_) - sigma_factors;

  // Predict sigma points to k + 1
  PredictSigmaPoints(X_aug, Xsig_aug, delta_t);

  // Update the state mean x_ and covariance P_ using weighted sum
  for(size_t i = 0; i < n_sigma_points_; i++)
  {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }
  
  MatrixXd error = Xsig_pred_ - x_.replicate(1, n_sigma_points_);
  for(size_t i = 0; i < n_sigma_points_; i++)
  {
    P_ += weights_(i) * (error.col(i) * error.col(i).transpose());
  }

}

void UKF::PredictSigmaPoints(VectorXd& x_k, MatrixXd& sigma_x_k, double dt)
{
  double px, py, radial_velocity, yaw, yaw_rate;
  double px_pred, py_pred, radial_velocity_pred, yaw_pred, yaw_rate_pred;
  for(size_t i=0; i < n_sigma_points_; i++)
  {
    px = sigma_x_k(0);
    py = sigma_x_k(1);
    radial_velocity = sigma_x_k(2);
    yaw = sigma_x_k(3);
    yaw_rate = sigma_x_k(4);

    VectorXd noise = VectorXd::Zero(n_x_);
    noise(0) = 0.5 * pow(dt, 2) * cos(yaw) * sigma_x_k(5);
    noise(1) = 0.5 * pow(dt, 2) * sin(yaw) * sigma_x_k(5);
    noise(2) = dt * sigma_x_k(5);
    noise(3) = 0.5 * pow(dt, 2) * sigma_x_k(6);
    noise(4) = dt * sigma_x_k(6);

    VectorXd sigma_pred = VectorXd::Zero(n_x_);
    if(radial_velocity > 0.01)
    {
      double scale = radial_velocity / yaw;
      sigma_pred(0) = scale * (sin(yaw + (yaw_rate *dt)) - sin(yaw));
      sigma_pred(1) = scale * (-cos(yaw + (yaw_rate * dt)) + cos(yaw));
      sigma_pred(3) = yaw_rate * dt;  
    }
    else
    {
      sigma_pred(0) = radial_velocity * cos(yaw) * dt;
      sigma_pred(1) = radial_velocity * sin(yaw) * dt;
      sigma_pred(3) = yaw_rate * dt;
    }
    Xsig_pred_.col(i) = x_k + sigma_pred + noise;
  }  
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

   MatrixXd I = MatrixXd::Identity(n_x_, n_x_);

   MatrixXd H = MatrixXd::Zero(n_obs_lidar_, n_x_);
   H(0, 0) = 1;
   H(1, 1) = 1;

   MatrixXd R = MatrixXd::Zero(n_obs_lidar_, n_obs_lidar_);
   R(0, 0) = pow(std_laspx_, 2);
   R(1, 1) = pow(std_laspy_, 2);

   // Since the measurement transfrom function for the lidar is simply a selection matrix
   // we can directly use the values from the state x_
   VectorXd y = meas_package.raw_measurements_ - x_.head(2);
   MatrixXd S = H * P_* H.transpose()  + R;
   MatrixXd K = P_ * H.transpose() * S.inverse();


   x_ += K * y;
   P_ = (I - (K * H)) * P_;
   
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

   // Project sigma points into measurement space
   MatrixXd Zsig = MatrixXd::Zero(n_obs_radar_, n_sigma_points_);
   double px, py, radial_velocity, yaw;
   for(size_t i = 0; i < n_sigma_points_; i++)
   {
    px = Xsig_pred_(i, 0);
    py = Xsig_pred_(i, 1);
    radial_velocity = Xsig_pred_(i, 2);
    yaw = Xsig_pred_(i, 3);

    double radial_dist = sqrt(pow(px, 2) + pow(py, 2));
    Zsig(i, 0) = radial_dist;
    Zsig(i, 1) = atan2(py/px);
    Zsig(i, 2) = ( radial_velocity * (px * cos(yaw) + py * sin(yaw))) / radial_dist; 
   }

   VectorXd z = VectorXd::Zero(n_obs_radar_);
   MatrixXd S = MatrixXd::Zero(n_obs_radar_, n_obs_radar_);
   for(size_t i = 0; i < n_sigma_points_; i++)
   {
    z += weights_(i) * Zsig.col(i);
   }

   MatrixXd z_diff = MatrixXd::Zero(n_obs_radar_, n_sigma_points_);
   z_diff = Zsig - z.replicate(1, n_sigma_points_);
   for(size_t i = 0; i < n_sigma_points_; i++)
   {
    S += weights_(i) * (z_diff.col(i) * z_diff.col(i).transpose());
   }
   S(0, 0) += pow(std_radr_, 2);
   S(1, 1) += pow(std_radphi_, 2);
   S(2, 2) += pow(std_radrd_, 2);


  MatrixXd T = MatrixXd::Zero(n_x_, n_obs_radar_);

  for(size_t i = 0; i < n_sigma_points_; i++)
  {
    T += weights_(i) * ( (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z).transpose());
  }

  MatrixXd K = T * S.inverse();

  x_ += K * (meas_package.raw_measurements_ - z);

  P_ -= K * S * K.transpose();

}