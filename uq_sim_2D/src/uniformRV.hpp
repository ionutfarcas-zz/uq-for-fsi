#ifndef UNIFORMRV_HPP_
#define UNIFORMRV_HPP_

#include "rand_num_gen.hpp"

class UniformRandomVariable : public RandomNumberGenerator
{
private:
	std::random_device dev;
	std::mt19937* rng;
	std::uniform_real_distribution<double>* var_uniform;

public:
	UniformRandomVariable()
	{
		rng = new std::mt19937(dev());
		var_uniform = new std::uniform_real_distribution<double>(0.0, 1.0);
	}

	virtual double get_rn() const
	{
		return (*var_uniform)(*rng);
	}

	virtual std::vector<double> get_samples(const double& x, const double& y, const int& nsamples) const
	{
		std::vector<double> samples_temp;
		double temp = 0.0;

		for (int i = 0  ; i < nsamples ; ++i) 
		{
			temp = fabs(x + y*this->get_rn());
			samples_temp.push_back(temp);
		}

		return samples_temp;
	}

	~UniformRandomVariable() 
	{
		delete rng;
		rng = nullptr;

		delete var_uniform;
		var_uniform = nullptr;
	}
};

#endif /* UNIFORMRV_HPP_ */