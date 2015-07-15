#ifndef UNIFORMRV_HPP_
#define UNIFORMRV_HPP_

#include "rand_num_gen.hpp"

class UniformRandomVariable : public RandomNumberGenerator
{
private:
	boost::mt19937 rng;
	boost::uniform_real<> uniform_distr;
	boost::variate_generator<boost::mt19937&, boost::uniform_real<> >* var_uniform;

public:
	UniformRandomVariable()
	{
		rng = boost::mt19937();
		uniform_distr = boost::uniform_real<>(-1.0, 1.0);
		var_uniform = new boost::variate_generator<boost::mt19937&, boost::uniform_real<> >(rng, uniform_distr);	
	}

	virtual double get_rn() const
	{
		return (*var_uniform)();
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
		if(var_uniform != 0)
			delete var_uniform;
	}
};

#endif /* UNIFORMRV_HPP_ */