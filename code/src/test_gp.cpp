#include <iostream>
#include <vector>

#include <Eigen/Dense> 

#include "gp.h"


int main(int argc, char* argv[]){

	// Note: for the rbf for this one, \sigma_f = sqrt(1.40)
	// l = 
// 	std::vector<double> xs{-1.50,-1.00,-0.75,-0.40,-0.25,0.00};	
// 	std::vector<double> ys{-1.70,-1.15,-0.40,0.20,0.55,0.9};

	std::vector<double> xs{-1.00,-0.75,-0.40,-0.25,0.00,0.2,0.5,0.8,1.0};	
	std::vector<double> ys{0.305,-0.007,0.783,0.989,0.869,1.361,0.786,0.175,0.184};
	
	LinearFunction fn(0.5,0.5,0.5);
// 	RadialBasisFunction fn(0.5,0.5);
// 	RadialBasisFunction fn(1.27,1);
	GaussianProcess gp(&fn);

	gp.observe(xs,ys);
//  	std::cout << "Before: " << std::endl << gp.kernel->params.transpose() << std::endl;
	gp.optimize();
/*
	VectorXd d(2);
	d(0) = 1.27;
	d(1) = 1.;

	std::cout << gp.objective(d) << " " << d.transpose() << std::endl;
	exit(-1);
*/
//  	std::cout << "After: " << std::endl << gp.kernel->params.transpose() << std::endl;

//	std::cout << gp.K << std::endl;
//	std::pair<double,double> eval = gp.predict(0.2);
//	std::cout << eval.first << "," << eval.second << std::endl;
//	exit(-1);

	std::vector<double> mean;
	std::vector<double> var;

	const unsigned int num_steps = 100;
	for (unsigned int i = 0; i < num_steps; ++i){
// 		double x = -2.0 + i*(0.3+2.0)/static_cast<double>(num_steps);
		double x = -1.0 + i*(2.0)/static_cast<double>(num_steps);
		std::pair<double,double> eval = gp.predict(x);
//		mean.push_back(eval.first);
//		var.push_back(eval.second);
// 		std::cout << x << ", " << eval.first << ", " << eval.second << std::endl;
	}
	return 0;
}

