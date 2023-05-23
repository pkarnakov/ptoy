#include <iostream>

#include "geometry.h"

std::ostream& operator<<(std::ostream& out, Vect v) {
  out << "(" << v.x << "," << v.y << ")";
  return out;
}
