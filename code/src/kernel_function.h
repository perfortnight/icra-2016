#ifndef kernel_function_h
#define kernel_function_h

#include "class_modifiers.h"

#include <Eigen/Dense>
#include <vector>

using Eigen::VectorXd;

class KernelFunction : virtual public uncopyable {
public:
	KernelFunction() {};
//	double operator()(double x, double xp){ return eval(x,xp);}
	virtual double eval(double x, double xp) = 0;
	VectorXd params;
private:

};

class LinearFunction : virtual public KernelFunction {
public:	
	LinearFunction(double slope, double offset, double intercept);
	virtual double eval(double x, double xp);
private:
};

class RadialBasisFunction : virtual public KernelFunction {
public:	
	RadialBasisFunction(double scale, double length);
	virtual double eval(double x, double xp);
private:
};

class FluctuatingBasisFunction : virtual public KernelFunction {
	public:
FluctuatingBasisFunction(double scale1, double scale2, double length1, double noise_scale);
	virtual double eval(double x, double xp);
};

class SinusoidalBasisFunction : virtual public KernelFunction {
public:
	SinusoidalBasisFunction(double scale, double length, double freq, double noise_scale);
	virtual double eval(double x, double xp);
};

#endif //kernel_function_h
