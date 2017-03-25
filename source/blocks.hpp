/*
  BLOCKS
*/

#include "geometry.hpp"
#include <cmath>
#include <vector>
#include "particles_system.hpp"
#include <cassert>

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
  size_t size() const
  {
    return static_cast<size_t>(N.i*N.j);
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
  size_t block_capacity;
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
    for(size_t n=0; n<B.size(); ++n)
    {
      B[n].reserve(block_capacity);
    }
  }
  void reinit() {
    mindex Nnew;
    Nnew.i=int(domain.size().x/block_size.x)+1;
    Nnew.j=int(domain.size().y/block_size.y)+1;
    if (Nnew.i > N.i || Nnew.j > N.j) {
      const auto Nold = N;
      N += N;
      const auto Bold = B;
      B = array2D<std::vector<particle>>(N);
      for (int i = 0; i < Nold.i; ++i) {
        for (int j = 0; j < Nold.j; ++j) {
          const mindex m(i, j);
          B[m] = std::move(Bold[m]);
        }
      }
    }
  }
public:
  mindex N;
  const size_t kBlockNone = static_cast<size_t>(-1);
  array2D<std::vector<particle>> B;
  std::vector<mindex> NEAR;
  size_t num_particles_;
  size_t close_packing(vect size, Scal r) const
  {
    // approximate value of volume occupied by close-packed circles in a rectangle with given size
    Scal occupied_volume=size.x*size.y*PI/(2.0*sqrt(3.0));
    // approximate number of close-packed circles
    Scal circles_number=occupied_volume/(PI*r*r);
    // truncate and add a reserve (+1)
    return int(circles_number)+1;
  }
  blocks(rect_vect _domain, vect _block_size, 
      size_t _block_capacity, bool) : domain(_domain), block_size(_block_size), block_capacity(_block_capacity)
  {
    init();
  }
  blocks(rect_vect _domain, vect _block_size, Scal particle_radius) : domain(_domain), block_size(_block_size), block_capacity(close_packing(_block_size, particle_radius))
  {
    init();
  }
  mindex get_block(vect p) const
  {
    vect rel=p-domain.A;
    return mindex(int(rel.x/block_size.x), int(rel.y/block_size.y));
  }
  mindex constraints(mindex m) const
  {
    return mindex(max(0,min(N.i-1,m.i)), max(0,min(N.j-1,m.j)));
  }
  size_t getn_check(vect p) const {
    const mindex m = get_block(p);
    if (m.i < 0 || m.i >= N.i || m.j < 0 || m.j >= N.j) {
      return kBlockNone;
    }
    return B.getn(m);
  }
  void add_particles(const std::vector<particle>& P) {
    for (auto& part : P) {
      auto n = constraints(get_block(part.p));
      //std::cout << "p=" << part.p.x << " " << part.p.y 
      //  << " b=" << B.getn(n) << std::endl;
      B[n].push_back(part);
      assert(B[n].size() <= block_capacity);
    }
  }
  void SetDomain(rect_vect new_domain) { 
    domain = new_domain;
    reinit();
  }
  size_t GetNumParticles() const { return num_particles_; }
  void arrange()
  {
    reinit();
    size_t lnum_particles_ = 0;
    for (size_t n = 0; n < B.size(); ++n) {
      size_t w = 0;
      size_t we = B[n].size();
      while (w < we) {
        size_t new_block = getn_check(B[n][w].p); 
        if (n != new_block) {
          if (new_block != kBlockNone) {
            B[new_block].push_back(B[n][w]);
            if (new_block < n) {
              ++lnum_particles_;
            }
          }
          std::swap(B[n][w], B[n].back());
          B[n].pop_back();
          --we;
        } else {
          ++w;
          ++lnum_particles_;
        }
      }
    }

    num_particles_ = lnum_particles_;
  }
  void print_status() const {
    std::cout << "Blocks" << std::endl;
    for (size_t n = 0; n < B.size(); ++n) {
      if (B[n].size()) { 
        std::cout << n << " " << B[n].size() << std::endl;
      }
    }
  }

};
