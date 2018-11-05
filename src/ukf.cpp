#include "ukf.h"
#include "tools.h"
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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
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

  is_initialized_ = false;  // Set to true after we receive our first measurement.
                            // See ProcessMeasurement().

  n_x_    = 5;  // Size of state vector
                // [px, py, v, psi (yaw angle), psi_dot (yaw_rate)]
  n_aug_  = 7;  // Size of augmented state vector
                // [state, nu_a (linear acceleration noise), nu_psi_ddot (angular acceleration noise)]
                // Number of sigma points
  n_sig_  = 2 * n_aug_ + 1;
                // Spreading parameter
  lambda_ = 3 - n_aug_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;


  Xsig_pred_ = MatrixXd(n_x_, n_sig_);  // # of rows == size of state vector
                                        // # of cols == number of sigma points

  // Calculate weights.
  weights_ = VectorXd(n_sig_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < n_sig_; ++i)
    weights_(i) = 0.5 / (lambda_ + n_aug_);
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

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */

    // first measurement
    cout << "UKF: " << endl;

    // Save the current timestamp so we can calculate delta next time.
    time_us_ = meas_package.timestamp_;
    current_step_ = 0;

    // Initial state (location, velocity).
    if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */

      // Set the state with the initial location and zero velocity and yaw angle.
      double  p_x = meas_package.raw_measurements_[0];
      double  p_y = meas_package.raw_measurements_[1];
      double  v   = 0;
      double  yaw = 0;
      double  yaw_rate = 0;

      x_ << p_x, p_y, v, yaw, yaw_rate;
      cout << "Initial measurement: LASER" << endl;
    }
    else if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      // Measurement vector: (rho, phi, rho_dot)
      double  rho     = meas_package.raw_measurements_[0];
      double  phi     = meas_package.raw_measurements_[1];
      double  rho_dot = meas_package.raw_measurements_[2];

      // State vector: (p_x, p_y, v, yaw, yaw_rate) -> same as for laser measurement
      NormalizeAngle(phi);
      double  cos_phi = cos(phi);
      double  sin_phi = sin(phi);

      double  p_x = rho * cos_phi;      // Translate polar coordinates to Cartesian.
      double  p_y = rho * sin_phi;
      double  v   = 0;
      double  yaw = 0;
      double  yaw_rate = 0;

      x_ << p_x, p_y, v, yaw, yaw_rate;
      cout << "Initial measurement: RADAR" << endl;
    }

    // Initial state covariance matrix.
    P_ = MatrixXd::Identity(n_x_, n_x_);

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  cout << "step: " << ++current_step_ << endl;

  // Calculate elapsed time in seconds.
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0F;
  time_us_ = meas_package.timestamp_;

  // Predict
  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
  else if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }

  // print the output
  //cout << "x_ = " << endl << x_ << endl;
  //cout << "P_ = " << endl << P_ << endl;
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

  /*****************************************************************************
   *  Generate sigma points
   ****************************************************************************/

  //MatrixXd  Xsig = MatrixXd(n_x_, n_sig_);

  //// Calculate square root of P
  //// (lower triangular matrix L of the matrix P such that P = L*L^)
  //MatrixXd A = P_.llt().matrixL();

  //// First sigma point is just the state vector.
  //Xsig.col(0) = x_;

  //// Set remaining sigma points based on the sigma point formula.
  //MatrixXd T = sqrt(lambda_ + n_x_) * A;
  //for (int i = 0; i < n_x_; i++)
  //{
  //  Xsig.col(1 + i)        = x_ + T.col(i);
  //  Xsig.col(n_x_ + 1 + i) = x_ - T.col(i);
  //}


  /*****************************************************************************
   *  Generate augmented sigma points
   ****************************************************************************/

  // Create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug[5] = 0; // longitudinal acceleration noise mean
  x_aug[6] = 0; // yaw acceleration noise mean

  // Create augmented state covariance.
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  // Create square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();

  // Create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
  Xsig_aug.col(0) = x_aug;
  MatrixXd T_aug = sqrt(lambda_ + n_aug_) * A_aug;
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(1 + i)          = x_aug + T_aug.col(i);
    Xsig_aug.col(n_aug_ + 1 + i) = x_aug - T_aug.col(i);
  }


  /*****************************************************************************
   *  Predict sigma points
   ****************************************************************************/

  // For each column/sigma point ...
  for (int i = 0; i < Xsig_pred_.cols(); ++i)
  {
    VectorXd  x_aug = Xsig_aug.col(i);

    // Create variables only for readability.
    double  px          = x_aug[0]; // position - x
    double  py          = x_aug[1]; // position - y
    double  v           = x_aug[2]; // speed
    double  psi         = x_aug[3]; // yaw
    double  psi_dot     = x_aug[4]; // yaw acceleration
    double  nu_a        = x_aug[5]; // longitudinal acceleration noise
    double  nu_psi_ddot = x_aug[6]; // yaw acceleration noise

    // Cache values for expensive operations.
    double  psi_dot_delta_t = psi_dot * delta_t;
    double  sin_psi         = sin(psi);
    double  cos_psi         = cos(psi);
    double  delta_t2        = delta_t * delta_t;

    // Calculate the state transition vector (avoid division by zero).
    VectorXd  state_trans = VectorXd(n_x_);
    state_trans <<
      ((fabs(psi_dot) < 0.0001) ? v * cos_psi * delta_t : (v / psi_dot) * ( sin(psi + psi_dot_delta_t) - sin_psi)),
      ((fabs(psi_dot) < 0.0001) ? v * sin_psi * delta_t : (v / psi_dot) * (-cos(psi + psi_dot_delta_t) + cos_psi)),
      0.0,
      psi_dot_delta_t,
      0.0;

    // Calculate the noise vector.
    VectorXd  noise = VectorXd(n_x_);
    noise <<
      0.5 * delta_t2 * cos_psi * nu_a,
      0.5 * delta_t2 * sin_psi * nu_a,
      delta_t * nu_a,
      0.5 * delta_t2 * nu_psi_ddot,
      delta_t * nu_psi_ddot;

    // Update the state vector.
    VectorXd  x = x_aug.head(n_x_); // Slice off the augmented (noise) values.
    VectorXd  x_upd = x + state_trans + noise;

    // Update the state prediction.
    Xsig_pred_.col(i) = x_upd;
  }


  /*****************************************************************************
   *  Predict mean and covarance
   ****************************************************************************/

  // Predict state mean
  x_.fill(0.0);
  for (int i = 0; i < n_sig_; ++i)
    x_ += weights_(i) * Xsig_pred_.col(i);

  // Predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; ++i)
  {
    // State difference
    VectorXd  x_diff = Xsig_pred_.col(i) - x_;

    // Angle normalization (4th element - index 3 - is yaw angle)
    NormalizeAngle(x_diff(3));

    P_ += weights_(i) * x_diff * x_diff.transpose();
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

  VectorXd  z = meas_package.raw_measurements_;
  int   n_z = 2;  // Laser sensor has two measurement values - px, py


  /*****************************************************************************
   *  Transform sigma points into measurement space
   ****************************************************************************/

  MatrixXd  Zsig = MatrixXd(n_z, n_sig_);

  // map process state (px, py) into measurement state (px, py)
  for (int i = 0; i < n_sig_; ++i)
  {
    double  px = Xsig_pred_(0, i);
    double  py = Xsig_pred_(1, i);

    Zsig(0, i) = px;
    Zsig(1, i) = py;
  }


  /*****************************************************************************
   *  Calculate measurement noise matrix
   ****************************************************************************/

  // Measurement noise matrix.
  MatrixXd R = MatrixXd(n_z, n_z);
  R.fill(0.0);
  R(0, 0) = std_laspx_ * std_laspx_;
  R(1, 1) = std_laspy_ * std_laspy_;


  /*****************************************************************************
   *  Update state
   ****************************************************************************/

  int meas_angle_index = -1;  // No measurement angle for laser measurements.
  nis_laser_ = UpdateHelper(z, n_z, meas_angle_index, Zsig, R);
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

  VectorXd  z = meas_package.raw_measurements_;
  int n_z = 3;  // Radar sensor has three measurement values - rho, phi, rho_dot


  /*****************************************************************************
   *  Transform sigma points into measurement space
   ****************************************************************************/

  MatrixXd  Zsig = MatrixXd(n_z, n_sig_);

  // map process state (px, py, v, psi, psi_dot) into measurement state (rho, phi, rho_dot)
  for (int i = 0; i < n_sig_; ++i)
  {
    double  px       = Xsig_pred_(0, i);
    double  py       = Xsig_pred_(1, i);
    double  v        = Xsig_pred_(2, i);
    double  yaw      = Xsig_pred_(3, i);
    double  yaw_rate = Xsig_pred_(4, i);

    double  range           = sqrt(px * px + py * py);
    double  bearing         = atan2(py, px);
    if (range < 0.0001) // Avoid division by zero.
      range = 0.0001;
    double  radial_velocity = (px * cos(yaw) * v + py * sin(yaw) * v) / range;

    // No need to normalize the bearing because atan2 will return a value in the range [-pi, pi]
    Zsig(0, i) = range;
    Zsig(1, i) = bearing;
    Zsig(2, i) = radial_velocity;
  }


  /*****************************************************************************
   *  Calculate measurement noise matrix
   ****************************************************************************/

  // Measurement noise matrix.
  MatrixXd R = MatrixXd(n_z, n_z);
  R.fill(0.0);
  R(0, 0) = std_radr_   * std_radr_;
  R(1, 1) = std_radphi_ * std_radphi_;
  R(2, 2) = std_radrd_  * std_radrd_;


  /*****************************************************************************
   *  Update state
   ****************************************************************************/

  int meas_angle_index = 1; // 2nd element - index 1 - is bearing
  nis_radar_ = UpdateHelper(z, n_z, meas_angle_index, Zsig, R);
}


/**
 * Updates the state and the state covariance matrix using a laser or radar measurement.
 * @param {VectorXd} z
 * @param {int} n_z
 * @param {int} meas_angle_index
 * @param {MatrixXd} Zsig
 * @param {MatrixXd} R
 * @return {double} nis
 */
double UKF::UpdateHelper(VectorXd z, int n_z, int meas_angle_index, MatrixXd Zsig, MatrixXd R) {

  /*****************************************************************************
   *  Predict measurement mean and covarance
   ****************************************************************************/

  // Calculate mean predicted measurement
  VectorXd  z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < n_sig_; ++i)
    z_pred += weights_(i) * Zsig.col(i);

  // Calculate measurement covariance matrix S
  MatrixXd  S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sig_; ++i)
  {
    // Measurement difference
    VectorXd  z_diff = Zsig.col(i) - z_pred;

    // Normalize measurement angle if indicated.
    if (meas_angle_index >= 0)
      NormalizeAngle(z_diff(meas_angle_index));

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // Add noise component to covariance matrix.
  S += R;


  /*****************************************************************************
   *  Update state
   ****************************************************************************/

  // Calculate cross correlation matrix, Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; ++i)
  {
    // State difference
    VectorXd  x_diff = Xsig_pred_.col(i) - x_;
    // Angle normalization (4th element - index 3 - is yaw angle)
    NormalizeAngle(x_diff(3));

    // Measurement difference (residual)
    VectorXd  z_diff = Zsig.col(i) - z_pred;
    // Normalize measurement angle if indicated.
    if (meas_angle_index >= 0)
      NormalizeAngle(z_diff(meas_angle_index));

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Calculate Kalman gain K
  MatrixXd  K = Tc * S.inverse();

  // Measurement difference (residual)
  VectorXd  z_diff = z - z_pred;  // actual measurement minus measurement prediction
  // Normalize measurement angle if indicated.
  if (meas_angle_index >= 0)
    NormalizeAngle(z_diff(meas_angle_index));

  // Update state mean vector and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();


  /*****************************************************************************
   *  Calculate and return NIS
   ****************************************************************************/

  return Tools::CalculateNIS(z, z_pred, S);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {double &} angle
 */
void UKF::NormalizeAngle(double &angle) {

    while (angle >  M_PI) angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
}
