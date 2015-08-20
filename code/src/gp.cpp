#include "gp.h"

#include <cassert>
#include <cmath>
#include <iostream>

const double PI = atan(1.0)*4.;


GaussianProcess::GaussianProcess(KernelFunction *fn) : K(0,0),kernel(fn), last_x_size(0), stale(true)
{
}

void GaussianProcess::observe(double x, double y)
{
	stale = true;
	last_x_size = xs.size();
	const unsigned int x_size = xs.size();
	const unsigned int y_size = ys.size();

	xs.resize(x_size+1);
	ys.resize(y_size+1);

	xs(x_size) = x;
	ys(y_size) = y;
}

void GaussianProcess::observe(std::vector<double> xs, std::vector<double> ys)
{
	assert(xs.size() == ys.size());
	stale = true;
	last_x_size = xs.size();
	const unsigned int x_size = this->xs.size();
	const unsigned int y_size = this->ys.size();
	const unsigned int x_len = xs.size();
	const unsigned int y_len = ys.size();

	this->xs.resize(x_size+x_len);
	this->ys.resize(y_size+y_len);

	for (unsigned int i = 0; i < x_len; ++i){
		this->xs(x_size+i) = xs[i];
		this->ys(x_size+i) = ys[i];
	}

// 	std::cout << this->xs.size() << std::endl;
}

void GaussianProcess::update_K(bool all)
{
	//const unsigned int start_idx = all? last_x_size : 0;
	const unsigned int start_idx = 0;
	const unsigned int x_len = xs.size() - start_idx;
// 	const unsigned int num_rows = K.rows();
// 	const unsigned int num_cols = K.cols();
	/*
	if (!all){
		assert(num_rows < last_x_size);
	}
	*/
	//K.resize(num_rows+x_len,num_cols+x_len);
	K.resize(x_len,+x_len);
	for (unsigned int i = 0; i < xs.size(); ++i) {
		for (unsigned int j = start_idx; j < xs.size(); ++j) {
			const double diff = kernel->eval(xs(j),xs(i));
			K(i,j) = diff;
			K(j,i) = diff;
		}
	}
	stale = false;
}

double GaussianProcess::objective(VectorXd const &params){

	// update K with new kernel params
	kernel->params = params;
	update_K(); // commented out true
		
	// This is the objective function from "Gaussian Proces for Regression:
	// a Quick Introduction"
	const double n = static_cast<double>(xs.size());
	const double a = -0.5*ys.transpose()*K.inverse()*ys;
	const double b = -0.5*std::log(K.determinant());
	const double c = -0.5*n*std::log(2.*PI);

	return -(a+b+c);
}

void GaussianProcess::optimize()
{
	const double dp = 0.1;
	const double dp_squared = dp*dp;
	const double step_size = 0.1;
	double err = 1000000000000000.0;

// 	const double beta = 0.8;

	unsigned int iter = 0;
// 	VectorXd params(params_orig.size());
	VectorXd params = kernel->params;
	VectorXd params_orig = kernel->params;
	double base = objective(params_orig);

	VectorXd gradient(params.size());
	MatrixXd hessian(params.size(),params.size());
	double alpha = 10000;
	while (err > 0.0001 && iter < 5000){
// 	while (iter < 100){
		params_orig = params;	
		for (unsigned int i = 0; i < params.size(); ++i){
			params(i) = params_orig(i) + dp;
			const double v1 = objective(params);
			params(i) = params_orig(i) - dp;
			const double v2 = objective(params);
			std::cout << v1 <<","<< v2<<std::endl;

			gradient(i) = (v1 - v2)/(2.*dp);
// 			hessian(i,i) = alpha + (v1 - 2.*base - v2)/dp_squared;
// 			alpha = 0.95*alpha;
			params(i) = params_orig(i);
/*
			for (unsigned j = i+1; j < params.size(); ++j){
				params(i) = params_orig(i) + dp;
				params(j) = params_orig(j) + dp;
				const double h1 = objective(params);

				params(i) = params_orig(i) + dp;
				params(j) = params_orig(j) - dp;
				const double h2 = objective(params);

				params(i) = params_orig(i) - dp;
				params(j) = params_orig(j) + dp;
				const double h3 = objective(params);
				
				params(i) = params_orig(i) - dp;
				params(j) = params_orig(j) - dp;
				const double h4 = objective(params);
	
				const double h_val = (h1 - h2 - h3 + h4)/(4.*dp_squared);
				hessian(i,j) = h_val;
				hessian(j,i) = h_val;
			}
*/
		}

// 		const VectorXd d_param = hessian.inverse()*gradient;
 		std::cout << "Before: " << params.transpose() << std::endl;
		const VectorXd d_param = step_size*gradient;
 		std::cout << d_param.transpose() << std::endl;
		params -= d_param; 
 		std::cout <<  "After: " << params.transpose() << std::endl;
		const double val = objective(params);
		std::cout << val <<std::endl;
		kernel->params = params;
		err = fabs(val - base);///fabs(base);
// 		std::cout << gradient.transpose() << std::endl;
// 		std::cout << err << ", " << val << std::endl;
		base = val;
		/*
		if (iter % 10){
			std::cout << val << " " << kernel->params.transpose() << std::endl;
		}
		*/
		iter += 1;
	} // end while (err > 0.0001)


	
}

// Returns mean and variance
std::pair<double,double> GaussianProcess::predict(double x){
// 	std::pair<double,double> retval(0.,0.);
// 	std::cout << stale << std::endl;
	if (stale){
// 		std::cout << "Updating K" << std::endl;
		update_K();
// 		std::cout << K << std::endl;
// 		exit(-1);
	}
	VectorXd k_star(xs.size());
	for (unsigned int i = 0; i < xs.size(); ++i){
		k_star(i) = kernel->eval(x,xs(i));	
	}
	const double k_starstar = kernel->eval(x,x);
	const double mean = k_star.transpose()*K.inverse()*ys;
	const double var = k_starstar - k_star.transpose()*K.inverse()*k_star;

	std::pair<double,double> retval(mean,var);
	return retval;
}

double GaussianProcess::probability(double f_x,double x){
	std::pair<double,double> p = predict(x);
	p.first = p.first;



	return 0.;
}

	
