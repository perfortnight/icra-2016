#ifndef gp_h
#define gp_h

#include <vector>
#include <utility>
#include "kernel_function.h"

#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class GaussianProcess : virtual public uncopyable {
	public:
		GaussianProcess(KernelFunction *fn);

		void observe(double x, double y);
		void observe(std::vector<double> xs, std::vector<double> ys);

		void optimize();
		// Returns mean and variance
		std::pair<double,double> predict(double x);
		double probability(double f_x, double x);

		// Make these public because we may want to access them
		// and we don't really need a getter/setter.
		VectorXd xs;
		VectorXd ys;
		MatrixXd K;

		KernelFunction* kernel;		
		double objective(VectorXd const &params);
	private:
		void update_K(bool all=false);
		unsigned int last_x_size;
		bool stale;

};


#endif //gp_h
