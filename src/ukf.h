#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "tools.h"

#include<math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  const int dim_x = 5; // the state dimension
  const int dim_x_joint = 7; // the dimension of the joint of x and linear acceleration-angular accelration
  const double inf_variance = 1;
  
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
  inline static double UpdateRadar(const MeasurementPackage& pack,
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
				) 
  {
    
    assert(pack.sensor_type_ == MeasurementPackage::RADAR);
    // skip update if degenrate measurement
    if (fabs(pack.raw_measurements_[0]) < 0.00001) {
     
      x_post = x;
      Sigma_post = Sigma;
      
      return 0; 
    }
    
    //set state dimension
    const int dim_x = 5;
    const int dim_joint = dim_x + 2; // this should be the dimension of the joint!!!!
    
    //Joint sigma points
    const int n_joint_sigma_points = 2 * dim_joint + 1;
    
    assert(Xsig.rows() == dim_x && Xsig.cols() ==  n_joint_sigma_points);
    assert(x.size() == dim_x);
    
    //measurement dimension
    const int dim_z = pack.raw_measurements_.size();
    
    // shortcut for the measurement vector
    VectorXd z = pack.raw_measurements_;
    // make suree angle is in the [-π, π] interval
    if (z[1] > 2 * M_PI) z[1] = -(2 * M_PI - z[1]);

    //The lambda par5ameter (need to look into this rfeally...)
    double lambda = 3 - dim_joint;

    double alpha_sq = (lambda + dim_joint) / (dim_joint + kappa);
    
    // Create storage for the predicted measurement sigma-points
    MatrixXd Zsig = MatrixXd(dim_z, n_joint_sigma_points);

    //mean predicted measurement
    VectorXd z_pred = VectorXd::Zero(dim_z);
  
    //Run the sigma points through the measurement model 
    double w_m[n_joint_sigma_points]; // mean weights ...
    double w_c[n_joint_sigma_points]; // covariance weights ...
    
    for (int i = 0; i < n_joint_sigma_points; i++) {
      // get the weight out of the way first...
      if (i == 0) {
	
	w_m[i] = lambda / (lambda + dim_joint);
	w_c[i] = w_m[i] + (1 - alpha_sq) + beta;
      
      } else 
	w_m[i] = w_c[i] = 0.5 / (lambda + dim_joint);
      
      // Npw get the state in named variables
      double px = (1 / sqrt(1 + gamma) ) * Xsig(0, i),
             py = (1 / sqrt(1 + gamma) ) * Xsig(1, i),
             v = (1 / sqrt(1 + gamma) ) * Xsig(2, i),
             phi = (1 / sqrt(1 + gamma) ) * Xsig(3, i);
      double theta = (1 / sqrt(1 + gamma) ) * atan2(py, px);
      
      double vx = v * cos(phi),
             vy = v * sin(phi);
             
   // 2. Now storing ne
      double rho = sqrt(px * px + py * py);
      double rho_dot = rho > 0.0001 ? (vx * px + vy * py) / rho : 0;
      
      Zsig(0, i) = rho;
      Zsig(1, i) = theta;
      Zsig(2, i) = rho_dot;
      
      z_pred += w_m[i] * Zsig.col(i);
      //std::cout<<"Zsig currently : "<<std::endl<<Zsig<<std::endl;
  }

  
  // Now, S is the so-called "predicted measurement covariance".
  // NOTE: What S actually is, is the covariance matrix of the measurfement variable
  //       in the joint distrbution of the measurement likelihood and the prior. 
  MatrixXd S = MatrixXd::Zero(dim_z, dim_z);

  for (int i = 0; i < n_joint_sigma_points; i++) 
    S += w_c[i] * (Zsig.col(i) - z_pred) * (Zsig.col(i) - z_pred).transpose();  
  
  
  // And lastly, add the likelihood noise variances straight to the diagonal
  S(0, 0) += std_radr * std_radr;
  S(1, 1) += std_radphi * std_radphi;
  S(2, 2) += std_radrd * std_radrd;
   
  // Print the darn S...
  //std::cout<<"S : "<<std::endl<<S<<std::endl;
  
   // 2. Now storing new state mean and covariance
  // Cross-correlation matrix Tc 
  MatrixXd Sigma_xz = MatrixXd::Zero(dim_x, dim_z);


  //calculate cross correlation matrix
  for (int i = 0; i < n_joint_sigma_points; i++) 
      Sigma_xz += w_c[i] * (Xsig.col(i) - x) *(Zsig.col(i) - z_pred).transpose(); 
  
  //calculate Kalman gain K;
  MatrixXd K = Sigma_xz * S.inverse();
  
  //update state mean and covariance matrix
  x_post = x + K * (z - z_pred);
  Sigma_post = Sigma - K*S*K.transpose();

  return (z - z_pred).dot( S.inverse() * (z - z_pred) );

  }
				 

				 
  // Do the Lidar update				
  inline static double UpdateLidar(const MeasurementPackage& pack,
				 const MatrixXd& Xsig,  // The predicted/prior distribution as points along 
						        // the covariance matrix principal directions.
				 VectorXd& x,     // The state mean
				 MatrixXd& Sigma, // The state covariance
				 double std_laspx,      // Lidar px stdev
				 double std_laspy,     // Lidar py stdev
				 double beta,
				 double kappa
				 ) 
  {
    
    assert(pack.sensor_type_ == MeasurementPackage::LASER);
    
    // skip if bad measurement
    if (pack.raw_measurements_[0] * pack.raw_measurements_[0] + 
        pack.raw_measurements_[1] * pack.raw_measurements_[1] < 0.00001) return 0;
    
    //set state dimension
    const int dim_x = 5;
 
    
    const int dim_joint = dim_x + 2;
    
    
    //number of joint sigma points
    const int n_joint_sigma_points = 2 * dim_joint + 1;

    assert(x.size() == dim_x);
    assert(Xsig.rows() == dim_x && Xsig.cols() == n_joint_sigma_points);
    
    //measurement dimension
    const int dim_z = pack.raw_measurements_.size();

    // shortcut name for the measurement vector
    VectorXd z = pack.raw_measurements_;
    
    //define spreading parameter
    double lambda = 3 - dim_joint;
    
    // The α^2 parameter
    double alpha_sq = (lambda + dim_joint) / (dim_joint + kappa);
    
    // Create storage for the predicted measurement sigma-points
    MatrixXd Zsig = MatrixXd(dim_z, n_joint_sigma_points);

    //mean predicted measurement
    VectorXd z_pred = VectorXd::Zero(dim_z);
  
    //Now run the sigma points through the measurement model equations
    double w_m[n_joint_sigma_points]; // mean weights
    double w_c[n_joint_sigma_points]; // covariance weights
    
    for (int i = 0; i < n_joint_sigma_points; i++) {
      // compute the weight first thing...
      if (i == 0) {
	
	w_m[i] = lambda / (lambda + dim_joint);
	w_c[i] = w_m[i] + (1 - alpha_sq + beta);
      
      } else 
	w_m[i] = w_c[i] = 0.5 / (lambda + dim_joint);
      
      // Now obtaining the predicted values for position
      double px =  Xsig(0, i),
             py = Xsig(1, i);
             
      Zsig(0, i) = px;
      Zsig(1, i) = py;
      
      z_pred += w_m[i] * Zsig.col(i);
      //std::cout<<"Zsig currently : "<<std::endl<<Zsig<<std::endl;
  }

  // Now do the "predicted measurement" covariance, S (essentially the Z-block in the z, x joint)
  //NOTE: S is the covariance of the measurement variable (say, y to distinguish from the instantiation, z)
  //      in the joint distribution of the measurement likelihood with the prior.
  MatrixXd S = MatrixXd::Zero(dim_z, dim_z);
  for (int i = 0; i < n_joint_sigma_points; i++) 
    S += w_c[i] * (Zsig.col(i) - z_pred) * (Zsig.col(i) - z_pred).transpose();  
  
    
  // And lastly, add the likelihood noise variances straight to the diagonal
  S(0, 0) += std_laspx * std_laspx;
  S(1, 1) += std_laspy * std_laspy;
  
    
  // print the darn thing...
  std::cout<<"S : "<<std::endl<<S<<std::endl;
  
   // 2. Now obtaining new state mean and covariance
  // Cross-correlation matrix Sigma_xz 
  MatrixXd Sigma_xz = MatrixXd::Zero(dim_x, dim_z);
  //calculate cross correlation matrix
  for (int i = 0; i < n_joint_sigma_points; i++) 
      Sigma_xz += w_c[i] * (Xsig.col(i) - x) *(Zsig.col(i) - z_pred).transpose(); 
  
  //calculate Kalman gain K;
  MatrixXd K = Sigma_xz * S.inverse();
  
  //std::cout<<"K : "<<std::endl<<S<<std::endl;
  //std::cout<<"K*S*K' : "<<std::endl<<K*S*K.transpose()<<std::endl;
  //std::cout<<"Sigma : "<<std::endl<<Sigma<<std::endl;
  //update state mean and covariance matrix
  x = x + K * (z - z_pred);
  Sigma = Sigma - K*S*K.transpose();
  
  // return the NIS
  return (z - z_pred).dot( S.inverse() * (z - z_pred) );
  
  }

  // The update 
  void Update(const MeasurementPackage&);
  
private:
  
  // Generate the Sigma points given a prior and covariance matrix
  inline static void GenerateSigmaPoints(const VectorXd& x,      // state mean
				         const MatrixXd& Sigma, // state covariance
				         MatrixXd& Xsig_out  // reference to sigma point matrix
 				    ) 
  {

	    //set state dimension
	    const int dim_x = x.size();
  
	    // number if sigma points
	    const int n_sigma_points = 2 * dim_x + 1;
  
	    //define spreading parameter
	    double lambda = 3 - dim_x;

	    //Cholesky decomposition of Sigma
	    MatrixXd L = Sigma.llt().matrixL();
	    /*Eigen::JacobiSVD<MatrixXd> svd(Sigma, Eigen::ComputeFullU | Eigen::ComputeFullV);
	    MatrixXd L = svd.matrixU();
	    VectorXd ss = svd.singularValues();
	    for (int i = 0; i<dim_x; i++) {
	      L.col(i) *= sqrt(ss[i]);
	      //std::cout<<"singular value : "<<ss[i]<<std::endl;
	      
	    }*/
	    
	    
	    if (Xsig_out.rows() != dim_x || Xsig_out.cols() != n_sigma_points) Xsig_out = MatrixXd(dim_x, n_sigma_points);
	    
	    
	    // fill-in the columns of Xsig
	    Xsig_out.col(0) = x;
	    for (int i = 1; i < dim_x; i++) {
    
	      Xsig_out.col(1 + i) = x + sqrt(lambda + dim_x) * L.col(i);
 
	      Xsig_out.col(1 + dim_x + i) = x - sqrt(lambda + dim_x) * L.col(i);
	    }
  
  }
  
  // Incorporate the "noisy" accelerations in a joint distribution described by new sigma points
  inline static void GenerateJointSigmaPoints(const VectorXd& x,     // prior mean
					      const MatrixXd& Sigma, // prior covariance
					      double std_a,          // linear accelration standard deivation
					      double std_omega_dot,      // angular accelration standard deviation
					      MatrixXd& Xsig_joint  // The joint distribution sigma points
					      )  
  {
    
    const int dim_x = 5;
    const int dim_joint = dim_x + 2;
    const int n_joint_sigma_points = 2 * dim_joint + 1;
    double lambda = 3 - dim_joint;
    
    // 1. Create the joint covariance matrix
    MatrixXd Sigma_joint = MatrixXd::Zero(dim_joint, dim_joint);
    Sigma_joint.block<dim_x, dim_x>(0, 0) = Sigma;
    Sigma_joint(dim_x, dim_x) = std_a * std_a;  // linear acceleration along the vehicle axis ("longitudinal") variance
    Sigma_joint(dim_x + 1, dim_x + 1) = std_omega_dot * std_omega_dot; // angular acceleration variance
    // Joint covariance done...
    
    // 2. Decomposing the joint covariance matrix
    // NOTE: This is redundant, because the ellipsoid gains simply two more axes uncorrelated with the existing ones
    // and everything can be worked-out without the  additional Cholesky, but hey...
    MatrixXd L = Sigma_joint.llt().matrixL();
    // ************** SVD - L ********************
    /*Eigen::JacobiSVD<MatrixXd> svd(Sigma_joint, Eigen::ComputeFullU | Eigen::ComputeFullV);
    MatrixXd L = svd.matrixU();
    VectorXd ss = svd.singularValues();
    for (int i = 0; i<dim_x; i++) {
      L.col(i) *= sqrt(ss[i]);
      //std::cout<<"singular value : "<<ss[i]<<std::endl;
	      
    }*/
    // ************** SVD - L ********************
    if (Xsig_joint.rows() != dim_joint || Xsig_joint.cols() != n_joint_sigma_points) Xsig_joint = MatrixXd(dim_joint, n_joint_sigma_points);
    
    // 3. Now filling the sigma-point matrix of the joint distribution
    VectorXd x_joint(dim_joint);
    x_joint.segment<dim_x>(0) = x;
    x_joint.col(0).segment<2>(dim_x) << 0, 0;
    Xsig_joint.col(0) = x_joint;
    
    // And the rest of the columns now...
    for (int i = 0; i < dim_joint; i++) {
      
      Xsig_joint.col(1 + i) = x_joint + sqrt(lambda + dim_joint) * L.col(i);
      Xsig_joint.col(1 + dim_joint + i) = x_joint - sqrt(lambda + dim_joint) * L.col(i);
      
    }
    // Done ...
  }


};

#endif /* UKF_H */
