#ifndef MEASUREMENT_PACKAGE_H_
#define MEASUREMENT_PACKAGE_H_

#include "Eigen/Dense"

using namespace Eigen;

class MeasurementPackage {
public:
  long long timestamp_;

  enum SensorType{
    LASER = 0,
    RADAR
  } sensor_type_;

  Eigen::VectorXd raw_measurements_;
  
  // Get the Jacobian of the Radar measurement model for
  // a given state vector x = [px py vx vy]
  inline static void RadarJacobian(const VectorXd &x, MatrixXd &J) {
    
    // Now assign the Jacobian
    double px = x[0], py = x[1], vx = x[2], vy = x[3];
    double rho = sqrt(px*px + py*py);
    
    // Case #1: radius is very small for some reason... In this case, we don't trust the sensor and continue on with the current estimate
    //          To this, we zero the Jacobian and the measurement will not affect the state and the Levenberg-Marquardt loop will exit...
    if (rho < 0.0000001)  J.setZero();
    else { 
      double rho_inv = 1.0 / rho;
      double rho_inv2 = rho_inv * rho_inv;
      double rho_inv3 = rho_inv2 * rho_inv;
    
      J(0, 0) = px * rho_inv;                       J(0, 1) = py * rho_inv;                     J(0, 2) = 0;            J(0, 3) = 0;
      J(1, 0) = -py * rho_inv2;                     J(1, 1) = px * rho_inv2;                    J(1, 2) = 0;            J(1, 3) = 0;
      J(2, 0) = py * (vx*py - vy*px) * rho_inv3;    J(2, 1) = px*(vy*px - vx*py) * rho_inv3;    J(2, 2) = px * rho_inv;   J(2, 3) = py * rho_inv; 
    }
    
  }
  
  
  inline static VectorXd EvaluateRadarModel(const VectorXd& x) {
    
    double px = x[0], py = x[1], vx = x[2], vy = x[3];
    double d = sqrt(px*px + py*py);
    if (d < 0.0000001)  return VectorXd::Zero(3);
    
    VectorXd ret(3);
    ret[0] = d;
    ret[1] = atan2(py, px);
    ret[2] = (px*vx + py*vy) / d; 
    
    return ret;
  }
  
};

#endif /* MEASUREMENT_PACKAGE_H_ */
