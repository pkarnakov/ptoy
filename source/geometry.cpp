#include "geometry.hpp"

std::ostream& operator<<(std::ostream& out, vect v)
{
  out<<"("<<v.x<<","<<v.y<<")";
  return out;
}