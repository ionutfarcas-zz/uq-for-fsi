#ifndef GAUSSHERMITEQUAD_HPP_
#define GAUSSHERMITEQUAD_HPP_

#include "gauss_quad.hpp"

class GaussHermiteQuadrature : public GaussQuadrature
{
public:
	GaussHermiteQuadrature() {};

	virtual double orthogonal_poly(const int& degree, const double& var) const
	{
		assert(degree >=0 );

		double poly_eval;

		if (degree == 0)
			poly_eval = 1.0;
		else if (degree == 1)
			poly_eval = var;
		else
			poly_eval = var*orthogonal_poly(degree - 1, var) - (degree - 1)*orthogonal_poly(degree - 2, var);

		return poly_eval;
	}

	virtual void quad_nodes_weights(const int& degree, std::vector<double>& weights, std::vector<double>& nodes) const
	{
		assert(degree >=0 );

		switch(degree)
		{
			case 1:
			nodes.push_back(0.0);
			weights.push_back(1.7725);
			break;

			case 2:
			nodes.push_back(-0.7071);
			nodes.push_back(0.7071);
			weights.push_back(0.8862);
			weights.push_back(0.8862);
			break;

			case 3:
			nodes.push_back(-1.2247);
			nodes.push_back(0.0000);
			nodes.push_back(1.2247);
			weights.push_back(0.2954);
			weights.push_back(1.1816);
			weights.push_back(0.2954);
			break;

			case 4:
			nodes.push_back(-1.6507);
			nodes.push_back(-0.5246);
			nodes.push_back(0.5246);
			nodes.push_back(1.6507);
			weights.push_back(0.0813);
			weights.push_back(0.8049);
			weights.push_back(0.8049);
			weights.push_back(0.0813);
			break;

			case 5:
			nodes.push_back(-2.0202);
			nodes.push_back(-0.9586);
			nodes.push_back(0.0);
			nodes.push_back(0.9586);
			nodes.push_back(2.0202);
			weights.push_back(0.0200);
			weights.push_back(0.3936);
			weights.push_back(0.9453);
			weights.push_back(0.3936);
			weights.push_back(0.0200);
			break;

			case 6:
			nodes.push_back(-2.3506);
			nodes.push_back(-1.3358);
			nodes.push_back(-0.4361);
			nodes.push_back(0.4361);
			nodes.push_back(1.3358);
			nodes.push_back(2.3506);
			weights.push_back(0.0045);
			weights.push_back(0.1571);
			weights.push_back(0.7246);
			weights.push_back(0.7246);
			weights.push_back(0.1571);
			weights.push_back(0.0045);
			break;

			case 7:
			nodes.push_back(-2.6520);
			nodes.push_back(-1.6736);
			nodes.push_back(-0.8163);
			nodes.push_back(0.0000);
			nodes.push_back(0.8163);
			nodes.push_back(1.6736);
			nodes.push_back(2.6520);
			weights.push_back(0.0010);
			weights.push_back(0.0545);
			weights.push_back(0.4256);
			weights.push_back(0.8103);
			weights.push_back(0.4256);
			weights.push_back(0.0545);
			weights.push_back(0.0010);
			break;

			case 8:
			nodes.push_back(-2.9306);
			nodes.push_back(-1.9817);
			nodes.push_back(-1.1572);
			nodes.push_back(-0.3812);
			nodes.push_back(0.3812);
			nodes.push_back(1.1572);
			nodes.push_back(1.9817);
			nodes.push_back(2.9306);
			weights.push_back(0.0002);
			weights.push_back(0.0171);
			weights.push_back(0.2078);
			weights.push_back(0.6611);
			weights.push_back(0.6611);
			weights.push_back(0.2078);
			weights.push_back(0.0171);
			weights.push_back(0.0002);
			break;

			case 9:
			nodes.push_back(-3.1910);
			nodes.push_back(-2.2666);
			nodes.push_back(-1.4686);
			nodes.push_back(-0.7236);
			nodes.push_back(0.0);
			nodes.push_back(0.7236);
			nodes.push_back(1.4686);
			nodes.push_back(2.2666);
			nodes.push_back(3.1910);
			weights.push_back(0.0000);
			weights.push_back(0.0049);
			weights.push_back(0.0885);
			weights.push_back(0.4327);
			weights.push_back(0.7202);
			weights.push_back(0.4327);
			weights.push_back(0.0885);
			weights.push_back(0.0049);
			weights.push_back(0.0000);
			break;

			case 10:
			nodes.push_back(-3.4362);
			nodes.push_back(-2.5327);
			nodes.push_back(-1.7567);
			nodes.push_back(-1.0366);
			nodes.push_back(-0.3429);
			nodes.push_back(0.3429);
			nodes.push_back(1.0366);
			nodes.push_back(1.7567);
			nodes.push_back(2.5327);
			nodes.push_back(3.4362);
			weights.push_back(0.0000);
			weights.push_back(0.0013);
			weights.push_back(0.0339);
			weights.push_back(0.2401);
			weights.push_back(0.6109);
			weights.push_back(0.6109);
			weights.push_back(0.2401);
			weights.push_back(0.0339);
			weights.push_back(0.0013);
			weights.push_back(0.0000);
			break;

			default:
			std::cout << "Please input a degree less or equal to 10!" << std::endl;
			break;
		}
	}
	
	virtual double norm_factor(const int& n) const
	{
		assert(n >=0 );

		double factorial;

		if(n == 0)
			factorial = 1.0;
		else
			factorial = n*norm_factor(n-1);

		return factorial;	
	}

	~GaussHermiteQuadrature() {}
};

#endif /* GAUSSHERMITEQUAD_HPP_ */