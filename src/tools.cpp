#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd  rmse(4);  // x, y, vx, vy
  rmse << 0, 0, 0, 0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  if (0 == estimations.size())
    cout << "Empty estimation vector." << endl;
  //  * the estimation vector size should equal ground truth vector size
  else if (estimations.size() != ground_truth.size())
    cout << "Estimation vector size and ground truth vector size are different." << endl;
  else {
    //accumulate squared residuals
    vector<VectorXd>::const_iterator est = estimations.begin();
    vector<VectorXd>::const_iterator truth = ground_truth.begin();
    do {
      VectorXd    residual = *est - *truth;
      residual = residual.array() * residual.array();
      rmse += residual;
    } while (++est != estimations.end() && ++truth != ground_truth.end());

    //calculate the mean
    rmse /= (double) estimations.size();

    //calculate the square root
    rmse = rmse.array().sqrt();
  }

  //return the result
  return rmse;
}