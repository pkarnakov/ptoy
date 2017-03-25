/*
  BLOCKS
*/

#include "geometry.hpp"
#include <cmath>
#include <vector>
#include "particles_system.hpp"
#include <cassert>

class blocks
{
 public:
  static const size_t kNumNeighbors = 9;
  static const size_t kBlockNone = static_cast<size_t>(-1);
  using DataVect = std::vector<std::vector<vect>>;
  using DataInt = std::vector<std::vector<int>>;
  struct BlockData {
    DataVect position, position_tmp, velocity, velocity_tmp, force;
    DataInt id;
    void resize(size_t size) {
      position.resize(size);
      position_tmp.resize(size);
      velocity.resize(size);
      velocity_tmp.resize(size);
      force.resize(size);
      id.resize(size);
    }
    void RemoveParticle(
        size_t src, // source block 
        size_t idx  // particle index within the source block
        ) {
      std::swap(position[src][idx], position[src].back());
      std::swap(position_tmp[src][idx], position_tmp[src].back());
      std::swap(velocity[src][idx], velocity[src].back());
      std::swap(velocity_tmp[src][idx], velocity_tmp[src].back());
      std::swap(force[src][idx], force[src].back());
      std::swap(id[src][idx], id[src].back());

      position[src].pop_back();
      position_tmp[src].pop_back();
      velocity[src].pop_back();
      velocity_tmp[src].pop_back();
      force[src].pop_back();
      id[src].pop_back();
    }
    void MoveParticle(
        size_t src, // source block 
        size_t idx, // particle index within the source block
        size_t dest // destination block
        ) {
      position[dest].push_back(position[src][idx]);
      position_tmp[dest].push_back(position_tmp[src][idx]);
      velocity[dest].push_back(velocity[src][idx]);
      velocity_tmp[dest].push_back(velocity_tmp[src][idx]);
      force[dest].push_back(force[src][idx]);
      id[dest].push_back(id[src][idx]);

      RemoveParticle(src, idx);
    }
    void AddParticle(
        size_t dest, // destination block
        vect particle_position,
        vect particle_velocity,
        int particle_id
        ) {
      position[dest].push_back(particle_position);
      position_tmp[dest].push_back(vect::kNan);
      velocity[dest].push_back(particle_velocity);
      velocity_tmp[dest].push_back(vect::kNan);
      force[dest].push_back(vect::kNan);
      id[dest].push_back(particle_id);
    }
  };
  blocks(rect_vect domain, vect block_size) 
      : domain_(domain), block_size_(block_size), num_particles_(0) {
    // Determine the number of blocks
    dims_.i = int(domain_.size().x / block_size.x) + 1;
    dims_.j = int(domain_.size().y / block_size.y) + 1;
    num_blocks_ = static_cast<size_t>(dims_.i * dims_.j);

    assert(dims_.i > 0 && dims_.j > 0 && num_blocks_ > 0);

    data_.resize(num_blocks_);
    std::cout << "foo" << GetNumBlocks() << std::endl;

    // Calc offsets to neighbors
    size_t n = 0;
    for (int i = -1; i <= 1; ++i) {
      for (int j = -1; j <= 1; ++j) {
        neighbor_offsets_[n] = i * dims_.j + j; 
        ++n;
      }
    }
  }
  size_t FindBlock(vect position) const {
    mindex m(static_cast<int>((position.x - domain_.A.x) / block_size_.x),
             static_cast<int>((position.y - domain_.A.y) / block_size_.y));
    if (m.i >= 0 && m.i < dims_.i && m.j >= 0 && m.j < dims_.j) {
      return m.i * dims_.j + m.j;
    }
    return kBlockNone;
  }
  void AddParticles(
      const std::vector<vect>& position, 
      const std::vector<vect>& velocity,
      const std::vector<int>& id
      ) {
    const size_t size = position.size();
    assert(velocity.size() == size);
    assert(id.size() == size);

    for (size_t p = 0; p < size; ++p) {
      size_t i = FindBlock(position[p]);
      if (i != kBlockNone) {
        data_.AddParticle(i, position[p], velocity[p], id[p]);
      }
    }
  }
  BlockData& GetData() {
    return data_;
  }
  size_t GetNumBlocks() const {
    return num_blocks_;
  }
  size_t GetNumParticles() const { 
    return num_particles_; 
  }
  std::array<int, kNumNeighbors> GetNeighborOffsets() const {
    return neighbor_offsets_;
  }
  void SortParticles() {
    size_t lnum_particles_ = 0;
    for (size_t i = 0; i < num_blocks_; ++i) {
      size_t p = 0;
      size_t pe = data_.position[i].size();
      while (p < pe) {
        const size_t j = FindBlock(data_.position[i][p]);
        if (i != j) {
          if (j != kBlockNone) {
            data_.MoveParticle(i, p, j);
            if (j < i) {
              ++lnum_particles_;
            }
          } else {
            data_.RemoveParticle(i, p);
          }
          --pe;
        } else {
          ++p;
          ++lnum_particles_;
        }
      }
    }

    num_particles_ = lnum_particles_;
  }

 private:
  rect_vect domain_;
  vect block_size_;
  BlockData data_;
  mindex dims_;
  size_t num_blocks_;
  std::array<int, kNumNeighbors> neighbor_offsets_;
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

};
