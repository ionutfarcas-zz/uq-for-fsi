#ifndef GAUSSLEGENDREQUAD_HPP_
#define GAUSSLEGENDREQUAD_HPP_

#include "gauss_quad.hpp"

class GaussLegendreQuadrature: public GaussQuadrature
{
public:
	GaussLegendreQuadrature() {}

	virtual double orthogonal_poly(const int& degree, const double& var) const
	{
		assert(degree >=0 );

		double poly_val;

		if (degree == 0)
			poly_val = 1.0;
		else if (degree == 1)
			poly_val =  var;
		else
			poly_val =  ((2*degree - 1)*var*orthogonal_poly(degree-1, var) - (degree - 1)*orthogonal_poly(degree - 2, var))/degree;

		return poly_val;
	}

	virtual void quad_nodes_weights(const int& degree, std::vector<double>& weights, std::vector<double>& nodes) const
	{
		assert(degree >=0 );

		switch(degree)
		{
			case 1:
			nodes.push_back(0.0000);
			weights.push_back(2.0000);
			break;

			case 2:
			nodes.push_back(-0.5773);
			nodes.push_back(0.5773);
			weights.push_back(1.0000);
			weights.push_back(1.0000);
			break;

			case 3:
			nodes.push_back(-0.77459);
			nodes.push_back(0.0000);
			nodes.push_back(0.77459);
			weights.push_back(0.5555);
			weights.push_back(0.8888);
			weights.push_back(0.5555);
			break;

			case 4:
			nodes.push_back(-0.8611);
			nodes.push_back(-0.3399);
			nodes.push_back(0.3399);
			nodes.push_back(0.8611);
			weights.push_back(0.3478);
			weights.push_back(0.6521);
			weights.push_back(0.6521);
			weights.push_back(0.3478);
			break;

			case 5:
			nodes.push_back(-0.9061);
			nodes.push_back(-0.5384);
			nodes.push_back(0.0000);
			nodes.push_back(0.5384);
			nodes.push_back(0.9061);
			weights.push_back(0.2369);
			weights.push_back(0.4786);
			weights.push_back(0.5688);
			weights.push_back(0.4786);
			weights.push_back(0.2369);
			break;

			case 6:
			nodes.push_back(-0.9324);
			nodes.push_back(-0.6612);
			nodes.push_back(-0.2386);
			nodes.push_back(0.2386);
			nodes.push_back(0.6612);
			nodes.push_back(0.9324);
			weights.push_back(0.1713);
			weights.push_back(0.3607);
			weights.push_back(0.4679);
			weights.push_back(0.4679);
			weights.push_back(0.3607);
			weights.push_back(0.1713);
			break;

			case 7:
			nodes.push_back(-0.9491);
			nodes.push_back(-0.7415);
			nodes.push_back(-0.4058);
			nodes.push_back(0.0000);
			nodes.push_back(0.4058);
			nodes.push_back(0.7415);
			nodes.push_back(0.9491);
			weights.push_back(0.1294);
			weights.push_back(0.2797);
			weights.push_back(0.3818);
			weights.push_back(0.4179);
			weights.push_back(0.3818);
			weights.push_back(0.2797);
			weights.push_back(0.1294);

			break;

			case 8:
			nodes.push_back(-0.9602);
			nodes.push_back(-0.7966);
			nodes.push_back(-0.5255);
			nodes.push_back(-0.1834);
			nodes.push_back(0.1834);
			nodes.push_back(0.5255);
			nodes.push_back(0.7966);
			nodes.push_back(0.9602);
			weights.push_back(0.1012);
			weights.push_back(0.2223);
			weights.push_back(0.3137);
			weights.push_back(0.3626);
			weights.push_back(0.3626);
			weights.push_back(0.3137);
			weights.push_back(0.2223);
			weights.push_back(0.1012);
			break;

			case 9:
			nodes.push_back(-0.9681);
			nodes.push_back(-0.8360);
			nodes.push_back(-0.6133);
			nodes.push_back(-0.3242);
			nodes.push_back(0.0000);
			nodes.push_back(0.3242);
			nodes.push_back(0.6133);
			nodes.push_back(0.8360);
			nodes.push_back(0.9681);
			weights.push_back(0.0812);
			weights.push_back(0.1806);
			weights.push_back(0.2606);
			weights.push_back(0.3123);	
			weights.push_back(0.3302);
			weights.push_back(0.3123);
			weights.push_back(0.2606);
			weights.push_back(0.1806);
			weights.push_back(0.0812);
			break;

			case 10:
			nodes.push_back(-0.9739);
			nodes.push_back(-0.8650);
			nodes.push_back(-0.6794);
			nodes.push_back(-0.4333);
			nodes.push_back(-0.1488);
			nodes.push_back(0.1488);
			nodes.push_back(0.4333);
			nodes.push_back(0.6794);
			nodes.push_back(0.8650);
			nodes.push_back(0.9739);
			weights.push_back(0.0666);
			weights.push_back(0.1494);
			weights.push_back(0.2190);
			weights.push_back(0.2692);
			weights.push_back(0.2955);
			weights.push_back(0.2955);
			weights.push_back(0.2692);
			weights.push_back(0.2190);
			weights.push_back(0.1494);
			weights.push_back(0.0666);
			break;

			default:
			std::cout << "Please input a degree less or equal to 10" << std::endl;
			break;
		}
	}
	
	virtual double norm_factor(const int& n) const
	{
		assert(n >= 0);

		return 2.0/(2.0*n + 1.0);
	}

	~GaussLegendreQuadrature() {}
};

#endif /* GAUSSLEGENDREQUAD_HPP_ */