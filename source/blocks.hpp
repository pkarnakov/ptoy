/*
  BLOCKS
*/

#include "geometry.hpp"
#include <cmath>
#include <vector>

using std::min;
using std::max;

template<class T>
class array2D
{
  std::vector<T> D;
  mindex N;
public:
  array2D() {;}
  array2D(mindex size) : N(size)
  {
    D.resize(N.i*N.j);
  }
  void resize(mindex size)
  {
    N=size;
    D.resize(N.i*N.j);
  }
  int getn(mindex m) const
  {
    return N.i*m.j+m.i;
  }
  int getn(int i, int j) const
  {
    return N.i*j+i;
  }
  T& operator[](mindex m)
  {
    return D[getn(m)];
  }
  const T& operator[](mindex m) const
  {
    return D[getn(m)];
  }
  T& operator()(int i, int j)
  {
    return D[getn(i,j)];
  }
  const T& operator()(int i, int j) const
  {
    return D[getn(i,j)];
  }
  T& operator[](int n)
  {
    return D[n];
  }
  const T& operator[](int n) const
  {
    return D[n];
  }
  int size() const
  {
    return N.i*N.j;
  }
  mindex msize() const
  {
    return N;
  }
  mindex getm(int n) const
  {
    return mindex(n%N.i, n/N.j);
  }
  bool valid(mindex m) const
  {
    return m.i>=0 && m.j>=0 && m.i<N.i && m.j<N.j;
  }
};

class blocks
{
  rect_vect domain;
  vect block_size;
  int block_capacity;
  mindex N;
  void init()
  {
    for(int i=-1; i<=1; ++i)
    for(int j=-1; j<=1; ++j)
    {
      NEAR.push_back(mindex(i,j));
    }

    N.i=int(domain.size().x/block_size.x)+1;
    N.j=int(domain.size().y/block_size.y)+1;
    B.resize(N);
    for(int n=0; n<B.size(); ++n)
    {
      B[n].reserve(block_capacity);
    }
  }
public:
  array2D<std::vector<int>> B;
  std::vector<mindex> NEAR;
  int close_packing(vect size, double r) const
  {
    // approximate value of volume occupied by close-packed circles in a rectangle with given size
    double occupied_volume=size.x*size.y*PI/(2.0*sqrt(3.0));
    // approximate number of close-packed circles
    double circles_number=occupied_volume/(PI*r*r);
    // truncate and add a reserve (+1)
    return int(circles_number)+1;
  }
  blocks(rect_vect _domain, vect _block_size, int _block_capacity) : domain(_domain), block_size(_block_size), block_capacity(_block_capacity)
  {
    init();
  }
  blocks(rect_vect _domain, vect _block_size, double particle_radius) : domain(_domain), block_size(_block_size), block_capacity(close_packing(_block_size, particle_radius))
  {
    init();
  }
  mindex get_block(vect p) const
  {
    vect rel=p-domain.A;
    vect size=domain.size();
    return mindex(int(rel.x/size.x), int(rel.y/size.y));
  }
  mindex constraints(mindex m) const
  {
    return mindex(max(0,min(N.i-1,m.i)), max(0,min(N.j-1,m.j)));
  }
  void arrange(std::vector<vect> X)
  {
    for(int n=0; n<B.size(); ++n)
    {
      B[n].resize(0);
    }

    for(size_t k=0; k<X.size(); ++k)
    {
      B[constraints(get_block(X[k]))].push_back(k);
    }
  }

};
