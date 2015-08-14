#include "kernel_function.h"

#include <cmath>

RadialBasisFunction::RadialBasisFunction(double scale, double length) {
	this->params.push_back(scale);
	this->params.push_back(length);
} 

double RadialBasisFunction::eval(double x, double xp){
	const double dx = x - xp;
	const double scale = this->params[0];
	const double length = this->params[1];	
	const double retval = scale*scale*std::exp(-dx/(2*length*length));
	return retval;
}
