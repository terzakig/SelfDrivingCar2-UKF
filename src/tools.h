#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

namespace Tools { 

  
  /**
  * A helper method to calculate RMSE
  */
  inline Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimates, const std::vector<Eigen::VectorXd> &ground_truth) {
   
    assert(estimates.size() == ground_truth.size() && estimates.size() > 0);
    
    int n = 0; // I assume ground truth should not contain invalid entries, but just in case...
    VectorXd vRMS = VectorXd::Zero(estimates[0].size());
    for (int i = 0; i < estimates.size(); i++)
      if (estimates[i].size() > 0 && ground_truth[i].size() > 0) {
	
	VectorXd diff = estimates[i] - ground_truth[i];
	for (int j = 0; j<estimates[0].size() ; j++)
	  vRMS[j] += diff[j] * diff[j];
	
	n++;
      }
      for (int j =0; j < estimates[0].size(); j++)
	vRMS[j] = vRMS[j] / n; 
      return vRMS;
  }

  
  
}

#endif /* TOOLS_H_ */
