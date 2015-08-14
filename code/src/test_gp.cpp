#include <iostream>
#include <vector>

#include <Eigen/Dense> 

#include "gp.h"

using Eigen::MatrixXd;

int main(int argc, char* argv[]){

	// Note: for the rbf for this one, \sigma_f = sqrt(1.40)
	// l = 
	std::vector<double> xs{-1.50,-1.00,-0.75,-0.40,-0.25,0.00};	
	std::vector<double> ys{-1.70,-1.15,-0.40,0.20,0.55,0.9};

	MatrixXd m(2,2);
	m(0,0) = 3;
	m(1,0) = 2.5;
	m(0,1) = -1;
	m(1,1) = m(1,0) + m(0,1);
	std::cout << m << std::endl;	 
	
	return 0;
}
