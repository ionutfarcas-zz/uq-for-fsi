#ifndef NORMALRV_HPP_
#define NORMALRV_HPP_

#include "rand_num_gen.hpp"

class NormalRandomVariable : public RandomNumberGenerator
{
private:
	std::random_device dev;
	std::mt19937* rng;
	std::normal_distribution<double>* var_normal;

public:
	NormalRandomVariable()
	{
		rng = new std::mt19937(dev());
		var_normal = new std::normal_distribution<double>(0.0, 1.0);
	}

	virtual double get_rn() const
	{
		return (*var_normal)(*rng);
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

	~NormalRandomVariable()
	{
		delete rng;
		rng = nullptr;

		delete var_normal;
		var_normal = nullptr;
	}
};

#endif /* NORMALRV_HPP_ */