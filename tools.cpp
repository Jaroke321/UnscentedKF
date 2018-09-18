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
  VectorXd rmse = VectorXd(4);
  rmse << 0, 0, 0, 0;

  //check for valid sizes of vectors
  if(estimations.size() != ground_truth.size()
     || estimations.size() == 0)
     {
        std::cout << "invalid size of either ground truth or estimations vector";
        return rmse;
     }
  //accumulate squared residuals
  for(unsigned int i = 0; i < estimations.size(); ++i){

    VectorXd residual = estimations[i] - ground_truth[i];

    //coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }
  //calculate the mean
  rmse = rmse.array() / estimations.size();

  //calculate the square root
  rmse = rmse.array().sqrt();

  //return result
  return rmse;
}
