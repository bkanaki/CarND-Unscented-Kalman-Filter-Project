#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
	rmse << 0,0,0,0;

  // check if the estimations is non empty
  if (estimations.size() < 1) {
	    cout << "Estimations are not present" << endl;
	    return rmse;
	}

  // check if the vector length of estimates and ground truth data match
  if (estimations.size() != ground_truth.size()) {
    cout << "Size of estimations and ground truth do not match" << endl;
    return rmse;
  }

	// accumulate squared residuals
	VectorXd residual(4);
	VectorXd square(4);
	for(int i=0; i < estimations.size(); ++i){
    residual = (estimations[i] - ground_truth[i]);
    square = residual.array() * residual.array();
    rmse = rmse + square;
	}

	// calculate the mean
	rmse = rmse/estimations.size();

	// calculate the squared root
	rmse = rmse.array().sqrt();

	// return the result
	return rmse;
}