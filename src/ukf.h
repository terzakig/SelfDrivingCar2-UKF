#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  const int dim_x = 5; // the state dimension
  const int dim_x_joint = 7; // the dimension of the joint of x and linear acceleration-angular accelration
  const double inf_variance = .1;
  
  // Sigma-point weights - related parameters...
  //const double beta = 6; 
  const double beta = 0;
  //const double kappa = 3;
  const double kappa = 0;
  
  
  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd Sigma_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  MatrixXd Xsig_joint_;
  
  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;


  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  // previous sample timespamp
  double previous_timestamp_;
  
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  // initilize the filter
  void Init(const MeasurementPackage&);
  
  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(const MeasurementPackage&);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  
  static double UpdateRadar(const MeasurementPackage& pack,
				 const MatrixXd& Xsig,  // The predicted/prior distribution as points along 
						        // the covariance matrix principal directions.
				 const VectorXd& x,     // The prior mean
				 const MatrixXd& Sigma, // The prior covariance
				 double std_radr,      // radar range stdev
				 double std_radphi,     // radar angle-to-object stdev
				 double std_radrd,     // radar range-rate stdev
				 VectorXd& x_post,      // The posterior mean
				 MatrixXd& Sigma_post,  // The posterior covariance
				                       // NOTE: We multiply the entries of Xsig with sqrt(1/(1 + kappa))
				 double beta,
				 double kappa,
				 double gamma = 0      // A regularization parameter (employed by iterative updates) 
				); 
  
				 
// Do the Lidar update				
  static double UpdateLidar(const MeasurementPackage& pack,
				 const MatrixXd& Xsig,  // The predicted/prior distribution as points along 
						        // the covariance matrix principal directions.
				 VectorXd& x,     // The state mean
				 MatrixXd& Sigma, // The state covariance
				 double std_laspx,      // Lidar px stdev
				 double std_laspy,     // Lidar py stdev
				 double beta,
				 double kappa
				 ); 
				 
  

  // The update 
  void Update(const MeasurementPackage&);
  
private:
  
  static void GenerateSigmaPoints(const VectorXd& x,      // state mean
				         const MatrixXd& Sigma, // state covariance
				         MatrixXd& Xsig_out  // reference to sigma point matrix
 				    );
  
  inline static void GenerateJointSigmaPoints(const VectorXd& x,     // prior mean
					      const MatrixXd& Sigma, // prior covariance
					      double std_a,          // linear accelration standard deivation
					      double std_omega_dot,      // angular accelration standard deviation
					      MatrixXd& Xsig_joint  // The joint distribution sigma points
					      );
  


};

#endif /* UKF_H */
