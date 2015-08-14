#ifndef kernel_function_h
#define kernel_function_h

#include "class_modifiers.h"
#include <vector>

class KernelFunction : virtual public uncopyable {
public:
	KernelFunction() {};
	double operator()(double x, double xp){ return eval(x,xp);}
	std::vector<double> params;
private:
	virtual double eval(double x, double xp) = 0;

};

class RadialBasisFunction : virtual public KernelFunction {
public:	
	RadialBasisFunction(double scale, double length);
private:
	virtual double eval(double x, double xp);
};

#endif //kernel_function_h
