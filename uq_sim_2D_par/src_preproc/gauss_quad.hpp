#ifndef GAUSSQUAD_HPP_
#define GAUSSQUAD_HPP_

#include <iostream>
#include <vector>
#include <cassert>

class GaussQuadrature
{
public:
	virtual double orthogonal_poly(const int& degree, const double& var) const = 0;

	virtual void quad_nodes_weights(const int& degree, std::vector<double>& weights, std::vector<double>& nodes) const = 0;

	virtual double norm_factor(const int& n) const = 0;

	virtual ~GaussQuadrature() {}
};

#endif /* GAUSSQUAD_HPP_ */