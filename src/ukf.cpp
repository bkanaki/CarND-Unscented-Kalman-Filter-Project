#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>
#include <fstream>

// #define PI 3.14159265

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

const double pi = 3.14159265;
// function to make sure that the phi (second element) in the
// vector is limited between [-pi,pi]
static inline void limitAngle(VectorXd &y, int idx) {
  float val = y(idx);
  // cout << "val: " << val << endl;
  if (val > pi) {
    y(idx) = fmod(val, 2*pi);
  } else if (val < -pi) {
    y(idx) = fmod(val, -2*pi);
  }
}

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // set n_x first and use this instead of hard coding matrix dims
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // This link helps in making a guess about the noise parameters.  
  // https://discussions.udacity.com/t/how-to-make-an-educated-guess-of-process-noise-parameters/351920/2
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;
  
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
  // initialize x with all zeros and P with identity matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);
  x_.fill(0.0);
  n_aug_ = n_x_ + 2;  // 2 from the two noise parameters
  lambda_ = 3 - n_aug_;
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  n_z_r_ = 3;
  n_z_l_ = 2;
  R_r_ = MatrixXd(n_z_r_, n_z_r_);
  R_l_ = MatrixXd(n_z_l_, n_z_l_);
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
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho    = meas_package.raw_measurements_(0);
      float phi    = meas_package.raw_measurements_(1);

      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
      
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }
    // first measurement
    x_ << 1, 1, 1, 1, 0.1;

    // init covariance matrix
    P_ << 0.15,    0, 0, 0, 0,
             0, 0.15, 0, 0, 0,
             0,    0, 1, 0, 0,
             0,    0, 0, 1, 0,
             0,    0, 0, 0, 1;
    // Fill up the R matrix
    R_r_ << 
    std_radr_*std_radr_,                        0,                     0,
                       0, std_radphi_*std_radphi_,                     0,
                       0,                       0, std_radrd_*std_radrd_;
    R_l_ << std_laspx_*std_laspx_,                     0,
                                0, std_laspy_*std_laspy_;
    // store the timestamp
    time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    // cout << "Intialization is good!" << endl;
    return;
  }
  // find the delta time between this measurement time and previous time
  double dt = (meas_package.timestamp_ - time_us_) / 1.0e6;
  time_us_ = meas_package.timestamp_;

  // predict based on measurements
  // cout << "Doing prediction..." << endl;
  Prediction(dt);

  // update based on measurements
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    // cout << "Updating lidar..." << endl;
    UpdateLidar(meas_package);
    // cout << NIS_lidar_ << endl;
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // cout << "Updating radar..." << endl;
    UpdateRadar(meas_package);
    // cout << NIS_radar_ << endl;
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

  // Step 1: Generate sigma points
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.setZero();
  x_aug.head(5) = x_;

  //create augmented covariance matrix
  MatrixXd Q = MatrixXd(n_aug_-n_x_, n_aug_-n_x_);
  Q << std_a_*std_a_, 0,
       0, std_yawdd_*std_yawdd_;
  P_aug.setZero();
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(n_aug_ - n_x_, n_aug_ - n_x_) = Q;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
      Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  // Step 2: Predict sigma points
  //create matrix with predicted sigma points as columns
  // MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  double dt = delta_t;
  double dt2 = dt*dt;

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //predict sigma points
    double p_x     = Xsig_aug(0,i);
    double p_y     = Xsig_aug(1,i);
    double v       = Xsig_aug(2,i);
    double psi     = Xsig_aug(3,i);
    double psidot  = Xsig_aug(4,i);
    double n_a     = Xsig_aug(5,i);
    double n_psidd = Xsig_aug(6,i);
    
    double p_x_pred;
    double p_y_pred;
    
    //avoid division by zero
    if (abs(psidot) < 0.0001) {
      double t = psi + psidot * dt;
      p_x_pred = p_x + v * (sin(t) - sin(psi)) / psidot + 0.5 * dt2 * n_a * cos(psi);
      p_y_pred = p_y + v * (cos(psi) - cos(t)) / psidot + 0.5 * dt2 * n_a * sin(psi);
    } else {
      p_x_pred = p_x + v * dt * cos(psi) + 0.5 * dt2 * n_a * cos(psi);
      p_y_pred = p_y + v * dt * sin(psi) + 0.5 * dt2 * n_a * sin(psi);
    }
    
    double v_pred = v + dt * n_a;
    double psi_pred = psi + psidot*dt + 0.5 * dt2 * n_psidd;
    double psidot_pred = psidot + dt * n_psidd;
    
    //write predicted sigma points into right column

    Xsig_pred_(0,i) = p_x_pred;
    Xsig_pred_(1,i) = p_y_pred;
    Xsig_pred_(2,i) = v_pred;
    Xsig_pred_(3,i) = psi_pred;
    Xsig_pred_(4,i) = psidot_pred;
  }

  // Step 3: Predict mean and covariance
  //set weights
  const double w0 = lambda_ / (lambda_ + n_aug_);
  const double wi = 0.5 / (lambda_ + n_aug_);
  
  //predict state mean
  x_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    if (i != 0) {
      x_ += (wi * Xsig_pred_.col(i));
    } else {
      x_ += (w0 * Xsig_pred_.col(i));
    }
  }

  // cout << "in prediction x_: " << x_ << endl;
  
  //predict state covariance matrix
  P_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // residual
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    limitAngle(x_diff, 3);
    
    if (i != 0) {
      P_ += (wi * (x_diff * x_diff.transpose()));
    } else {
      P_ += (w0 * (x_diff * x_diff.transpose()));
    }
  }
  // cout << "in prediction P_: " << P_ << endl;
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
  // Step 1: Predict Lidar Measurement
  //set weights
  const double w0 = lambda_ / (lambda_ + n_aug_);
  const double wi = 0.5 / (lambda_ + n_aug_);

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_l_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_l_);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    if (i != 0) {
      z_pred = z_pred + wi * Zsig.col(i);
    } else {
      z_pred = z_pred + w0 * Zsig.col(i);
    }
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z_l_, n_z_l_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    if (i != 0) {
      S = S + wi * z_diff * z_diff.transpose();
    } else {
      S = S + w0 * z_diff * z_diff.transpose();
    }
  }

  // add measurement noise covariance matrix
  S += R_l_;

  // Step 2: Do the final updates

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_l_);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    if (i != 0) {
      Tc += wi* x_diff * z_diff.transpose();
    } else {
      Tc += w0* x_diff * z_diff.transpose();
    }
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  // extract measurement as VectorXd
  VectorXd z = meas_package.raw_measurements_;

  //residual
  VectorXd z_diff = z - z_pred;

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ += -1 * K * S * K.transpose();

  // calculate the NIS
  NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;
  // cout << "NIS Lidar: " << NIS_lidar_ << " 95\% Threshold: 0.21 - 5.9" << endl;
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
  // Step 1: Predict Radar Measurement
  //set weights
  const double w0 = lambda_ / (lambda_ + n_aug_);
  const double wi = 0.5 / (lambda_ + n_aug_);

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_r_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  // cout << "transforming sigma points" << endl;
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  // cout << "Calculating the mean update" << endl;
  VectorXd z_pred = VectorXd(n_z_r_);
  z_pred.fill(0.0);
  for (int i=0; i < 2 * n_aug_+1; i++) {
    if (i != 0) {
      z_pred = z_pred + wi * Zsig.col(i);
    } else {
      z_pred = z_pred + w0 * Zsig.col(i);
    }
  }

  //innovation covariance matrix S
  // cout << "calculating covariance update" << endl;
  MatrixXd S = MatrixXd(n_z_r_, n_z_r_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    limitAngle(z_diff, 1);

    if (i != 0) {
      S = S + wi * z_diff * z_diff.transpose();
    } else {
      S = S + w0 * z_diff * z_diff.transpose();
    }
  }

  // add measurement noise covariance matrix
  S += R_r_;

  // Step 2: Do the final updates
  // cout << "Step2: Doing final uopdates" << endl;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_r_);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  // cout << "Calculate cross correlation matrix" << endl;
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    limitAngle(z_diff, 1);

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    limitAngle(x_diff, 3);

    if (i != 0) {
      Tc += wi* x_diff * z_diff .transpose();
    } else {
      Tc += w0* x_diff * z_diff .transpose();
    }
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  // extract measurement as VectorXd
  VectorXd z = meas_package.raw_measurements_;

  //residual
  VectorXd z_diff = z - z_pred;
  limitAngle(z_diff, 1);

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ += -1 * K * S * K.transpose();

  // calculate the NIS
  // cout << "Calculating NIS radar" << endl;
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
  // cout << "NIS Radar: " << NIS_radar_ << " 95\% Threshold: 0.58 - 7.8" << endl;
}
