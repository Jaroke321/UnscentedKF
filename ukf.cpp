#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.5/2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/4;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  //true after first iteration
  is_initialized_ = false;

  // dimension of state vector
  n_x_ = 5;

  //size of augmented state matrix
  n_aug_ = 7;

  //lambda value
  lambda_ = max(3 - n_x_, 0);

  //previous time step
  time_us_ =  0;

  //initialize weights
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2*n_aug_+1; i++){
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_){
    // initialize x_ and P_
    // calculate time_us_
    time_us_ = meas_package.timestamp_;
    x_ = VectorXd(n_x_);

    //set P_ matrix
    P_ = MatrixXd(n_x_, n_x_);
    P_ << 0.5, 0, 0, 0, 0,
          0, 0.5, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    //alternate between radar and lidar data
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
        // set x_ for radar measurements
        float r = meas_package.raw_measurements_[0];
        float theta = meas_package.raw_measurements_[1];

        float px = r * cos(theta);
        float py = r * sin(theta);

        x_ << px, py, 5, 1, 1;

    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
        // set x_ for laser measurement
        x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 5, 1, 1;
    }

    is_initialized_ = true;

  }

  //calculate time_us_
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  //call prediction step
  Prediction(dt);

  //call UpdateLidar if laser measurement
  if(meas_package.sensor_type_ == MeasurementPackage::LASER){
    //call Update
    UpdateLidar(meas_package);
  }
  //call UpdateRadar if radar measurement
  else if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    //call update
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // calculate augmented sigma points
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // create augmented sigma points vector
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);

  // create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;
  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  // calculate augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++){
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  // use augmented sigma points to get predicted sigma points (Xsig_pred_)

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  for(int i = 0; i < 2*n_aug_+1; i++){

    //get all points for readability
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if(fabs(yawd) < 0.001){
        px_p = p_x + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
        py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise to predictions
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma points into correct column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  // get predicted x_
  x_.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++){
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // get predicted P_
  P_.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++){

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }


}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  //dimensions for Zsig
  int n_z = 2;
  //create matrix to convert sigma points into measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  for(int i = 0; i< 2*n_aug_+1; i++){

    //extract values for readability
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);

    Zsig(0, i) = p_x;
    Zsig(1, i) = p_y;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++){

    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++){

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //calculate S
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add noise
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;

  S = S + R;

  // get current measurements
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

  // create cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++){

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //x difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //kalman gain
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  // update
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  // dimension of Zsig
  int n_z = 3;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  //transform current sigma points into measurement space
  for(int i = 0; i < 2*n_aug_+1; i++){

    //extract values for better readability
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    //fill Zsig with new measurement model
    Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);
    Zsig(1, i) = atan2(p_y, p_x);
    Zsig(2, i) = (p_x*v1 + p_y * v2) / sqrt(p_x*p_x + p_y*p_y);
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++){

      z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

  //innovation covariance martix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++){

      //residual
      VectorXd z_diff = Zsig.col(i) - z_pred;

      //angle normalization
      while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
      while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

      S = S + weights_(i) * z_diff * z_diff.transpose();
    }

  //add noise
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;

  S = S + R;

  //create vector z for incoming measurements
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];

  //create cross correlation matrix
  MatrixXd Tc  = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++){

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    //state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //normalize angles
  while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
  while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

  //update x_ and P_
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();


}
