#include <iostream>
#include <cmath>
#include <memory>
#include "sgpp_base.hpp"

using namespace SGPP::base;

SGPP::float_t f(int dim, SGPP::float_t* x) 
{
  SGPP::float_t res = 1.0;

  for (int i = 0; i < dim; i++) 
  {
    res *= (tan(x[i])*tan(x[i]))*(1.0 + tan(x[i])*tan(x[i]))*exp(-tan(x[i])*tan(x[i])/2.0);
  }

  return res;
}

SGPP::float_t sg_quad(const int& dim, const int& level, const double& a, const double& b)
{
  SGPP::float_t res = 0.0;
  size_t grid_storage_size = 0;
  double no_grid_points = 0.0;
  double volume = 1.0;

  std::unique_ptr<SGPP::float_t[]> p(new SGPP::float_t[dim]);
  GridIndex* gp;

  std::shared_ptr<Grid> grid(Grid::createLinearGrid(dim));
  GridStorage* gridStorage = grid->getStorage(); 
  GridGenerator* gridGen = grid->createGridGenerator();
  gridGen->regular(level);

  grid_storage_size = gridStorage->size();
  no_grid_points = static_cast<double>(gridStorage->size());

  DataVector alpha(no_grid_points);
  alpha.setAll(0.0);

  for(size_t i = 0; i < grid_storage_size; ++i) 
  {
    gp = gridStorage->get(i);

    for(int j = 0 ; j < dim ; ++j)
    {
      p[j] =  a + (b-a)*gp->getCoord(j);
    }

    alpha[i] = f(dim, p.get());
  }

  for(int i = 0 ; i < dim ; ++i)
  {
    volume *= (b - a);
  }

  SGPP::op_factory::createOperationHierarchisation(*grid)->doHierarchisation(alpha);
  std::shared_ptr<OperationQuadrature> opQ(SGPP::op_factory::createOperationQuadrature(*grid));

  res = volume*opQ->doQuadrature(alpha);

  return res;
}

int main() 
{
  std::cout << sg_quad(2, 10, -M_PI/2, M_PI/2) << std::endl;
}