#ifndef RANDNUMGEN_HPP_
#define RANDNUMGEN_HPP_

#include <iostream>
#include <random>
#include <vector>
#include <cmath>

class RandomNumberGenerator
{
public:
    virtual double get_rn() const = 0;

    virtual std::vector<double> get_samples(const double& x, const double& y, const int& nsamples) const = 0;

    virtual ~RandomNumberGenerator() {}

};

#endif /* RANDNUMGEN_HPP_ */