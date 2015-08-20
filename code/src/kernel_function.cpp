#include "kernel_function.h"

#include <cmath>
#include <iostream>

const double PI = atan(1.0)*4.;

// Linear Function
LinearFunction::LinearFunction(double slope, double offset, double intercept){
	this->params.resize(3);
	this->params(0) = slope;
	this->params(1) = offset;
	this->params(2) = intercept;
}

double LinearFunction::eval(double x, double xp){
	const double sv = this->params(0);
	const double l  = this->params(1);
	const double sb = this->params(2);

	std::cout << x << "," << xp << "," << sv << "," << l << "," << sb << std::endl;
	return sb*sb + sv*sv*(x-l)*(xp-l);
}

// Radial Basis Function
RadialBasisFunction::RadialBasisFunction(double scale, double length) {
	this->params.resize(2);
	this->params(0) = scale;
	this->params(1) = length;
} 

double RadialBasisFunction::eval(double x, double xp){
	const double dx = x - xp;
	const double scale = this->params(0);
	const double length = this->params(1);	
	const double retval = scale*scale*std::exp(-dx*dx/(2*length*length));
	return retval + 0.3*0.3*(x == xp);
}

// Fluctuating Basis Function
FluctuatingBasisFunction::FluctuatingBasisFunction(double scale1, double scale2, double length1, double noise_scale)
{
	this->params.resize(4);
	this->params(0) = scale1;
	this->params(1) = scale2;
	this->params(2) = length1;
	this->params(3) = noise_scale;

}
double FluctuatingBasisFunction::eval(double x, double xp)
{
	const double dx = x-xp;
	const double s1 = this->params(0);
	const double s2 = this->params(1);
	const double l1 = this->params(2);
	const double l2 = 6.*l1;
	const double s3 = this->params(3);

	const double a = s1*s1*std::exp(-dx*dx/(2.*l1*l1));
	const double b = s2*s2*std::exp(-dx*dx/(2.*l2*l2));
	const double c = s3*s3*static_cast<float>(x == xp);
	const double retval = a+b+c;
	return retval;
}
// Sinusoidal Basis Function
SinusoidalBasisFunction::SinusoidalBasisFunction(double scale, double length, double freq, double noise_scale){
	this->params.resize(4);
	this->params(0) = scale;
	this->params(1) = length;
	this->params(2) = freq;
	this->params(3) = noise_scale;
}

double SinusoidalBasisFunction::eval(double x, double xp){
	const double dx = x - xp;

	const double s1 = this->params(0);
	const double l1 = this->params(1);
	const double v = this->params(2);
	const double s2 = this->params(3);

	const double sin_x = std::sin(v*PI*dx);

	const double a = s1*s1*std::exp(-dx*dx/(2.*l1*l1));
	const double b = std::exp(-2.*sin_x*sin_x);
	const double c =  s2*s2*static_cast<float>(x == xp);
		
	const double retval = a+b+c;
	return retval;
}
