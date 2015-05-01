#ifndef RANDNUMGEN_HPP_
#define RANDNUMGEN_HPP_

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

class RandomNumberGenerator
{
public:
    virtual double get_rn() const = 0;

    virtual std::vector<double> get_samples(const double& x, const double& y, const int& nsamples) const = 0;

    virtual ~RandomNumberGenerator() {}

};

#endif /* RANDNUMGEN_HPP_ */