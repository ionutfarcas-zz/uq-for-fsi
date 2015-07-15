#ifndef NORMALRV_HPP_
#define NORMALRV_HPP_

#include "rand_num_gen.hpp"

class NormalRandomVariable : public RandomNumberGenerator
{
private:
	boost::mt19937 rng;
	boost::normal_distribution<> normal_distr;
	boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >* var_normal;

public:
	NormalRandomVariable()
	{
		rng = boost::mt19937();
		normal_distr = boost::normal_distribution<>(0.0, 1.0);
		var_normal = new boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >(rng, normal_distr);
	}

	virtual double get_rn() const
	{
		return (*var_normal)();
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
		if(var_normal != 0)
			delete var_normal;
	}
};

#endif /* NORMALRV_HPP_ */