#include <array>
#include <iostream>
#include <vector>

#include "span.h"

struct P {
  float x, y;
};

int main() {
  std::array<int, 3> a{0, 1, 2};
  std::array<int, 3> b{3, 4, 5};
  span<int> sa{a};
  span<int> sb{b};

  std::vector<P> p{
      P{0.1, 0.2},
      P{0.1, 0.2},
      P{0.1, 0.2},
  };

  span<P> sp{p};

  sp[1].x = 0;

  span<std::array<float, 2>> sr{
      reinterpret_cast<std::array<float, 2>*>(&(*p.begin())), p.size()};
  sr[1][0] += 1;
  for (auto e : sp) {
    std::cout << e.x << ' ' << e.y << std::endl;
  }
}
