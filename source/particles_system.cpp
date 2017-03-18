#include "particles_system.hpp"

vect F12(vect p1, vect v1, vect p2, vect v2, double sigma, double R)
{
  const double alpha=12.0;
  const double beta=6.0;
   
  double eps=sigma/(pow(2.0,alpha)-pow(2.0,beta));

  double r=p1.dist(p2);
  double F=r>R?0.0:eps*(pow(R/r, alpha)-pow(R/r, beta));
  return (p1-p2)*(F/r);
}