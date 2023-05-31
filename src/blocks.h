#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <vector>
#include "aligned_allocator.h"
#include "geometry.h"

template <class T>
using Array = std::array<T, 16>;

using ArrayVect = Array<Vect>;
using ArrayInt = Array<int>;

template <class T>
using Data = std::vector<Array<T>, AlignedAllocator<Array<T>, 64>>;

class Blocks {
 public:
  static const size_t kNumNeighbors = 9;
  static const size_t kBlockNone = -1;
  struct BlockData {
   public:
    Data<Vect> position, position_tmp, velocity, velocity_tmp, force;
    Data<int> id;
    std::vector<int> block_size;
    BlockData() = delete;
    BlockData(Blocks* owner) : owner_(owner) {}
    void Clear() {
      position.clear();
      position_tmp.clear();
      velocity.clear();
      velocity_tmp.clear();
      force.clear();
      id.clear();
      block_size.clear();

      for (auto& pair : owner_->block_by_id_) {
        pair.first = kBlockNone;
      }
    }
    void Resize(size_t size) {
      position.resize(size);
      position_tmp.resize(size);
      velocity.resize(size);
      velocity_tmp.resize(size);
      force.resize(size);
      id.resize(size);
      block_size.resize(size, 0);

      // TODO: consider updating block_by_id_ here
      for (auto& pair : owner_->block_by_id_) {
        pair.first = kBlockNone;
      }
    }
    // src: source block
    // idx: particle index within the source block
    void RemoveParticle(size_t src, size_t idx) {
      const int last = block_size[src] - 1;
      std::swap(position[src][idx], position[src][last]);
      std::swap(position_tmp[src][idx], position_tmp[src][last]);
      std::swap(velocity[src][idx], velocity[src][last]);
      std::swap(velocity_tmp[src][idx], velocity_tmp[src][last]);
      std::swap(force[src][idx], force[src][last]);
      std::swap(id[src][idx], id[src][last]);

      owner_->block_by_id_[id[src][idx]] = {src, idx};
      owner_->block_by_id_[id[src][last]].first = kBlockNone;

      --block_size[src];
    }
    // src: source block
    // idx: particle index within the source block
    // dest: destination block
    void MoveParticle(size_t src, size_t idx, size_t dest) {
      const int end = block_size[dest];
      position[dest][end] = position[src][idx];
      position_tmp[dest][end] = position_tmp[src][idx];
      velocity[dest][end] = velocity[src][idx];
      velocity_tmp[dest][end] = velocity_tmp[src][idx];
      force[dest][end] = force[src][idx];
      id[dest][end] = id[src][idx];

      RemoveParticle(src, idx);
      owner_->block_by_id_[id[dest][end]] = {dest, end};
      ++block_size[dest];
    }
    void AddParticle(
        size_t dest, Vect particle_position, Vect particle_velocity,
        int particle_id) {
      const int end = block_size[dest];
      position[dest][end] = particle_position;
      position_tmp[dest][end] = GetNan<Vect>();
      velocity[dest][end] = particle_velocity;
      velocity_tmp[dest][end] = GetNan<Vect>();
      force[dest][end] = GetNan<Vect>();
      id[dest][end] = particle_id;
      const size_t pid = static_cast<size_t>(id[dest][end]);
      if (owner_->block_by_id_.size() <= pid) {
        owner_->block_by_id_.resize(pid + 1);
      }
      owner_->block_by_id_[pid] = {dest, block_size[dest]};
      ++block_size[dest];
    }

   private:
    Blocks* owner_;
  };
  void SetDomain(RectVect proposal) {
    RectVect domain = domain_;
    while (proposal.A.x < domain.A.x) {
      domain.A.x -= block_size_.x;
    }
    while (proposal.A.y < domain.A.y) {
      domain.A.y -= block_size_.y;
    }
    while (proposal.B.x > domain.B.x) {
      domain.B.x += block_size_.x;
    }
    while (proposal.B.y > domain.B.y) {
      domain.B.y += block_size_.y;
    }

    if (domain != domain_) {
      const BlockData old_data = GetData();
      InitEmptyBlocks(domain, block_size_);
      AddParticles(old_data);
      std::cout << "Create blocks for new domain: " << domain_.A << " "
                << domain_.B << std::endl;
    }
  }
  Scal GetCircumRadius() const {
    return block_size_.length() * 0.5;
  }
  Vect GetCenter(size_t block) const {
    const MIdx m(block / dims_.j, block % dims_.j);
    return domain_.A +
           Vect((0.5 + m.i) * block_size_.x, (0.5 + m.j) * block_size_.y);
  }
  Blocks(RectVect domain, Vect block_size)
      : data_(this), num_particles_(0), num_per_cell_(0) {
    InitEmptyBlocks(domain, block_size);
  }
  size_t FindBlock(Vect position) const {
    MIdx m(
        static_cast<int>((position.x - domain_.A.x) / block_size_.x),
        static_cast<int>((position.y - domain_.A.y) / block_size_.y));
    if (m.i >= 0 && m.i < dims_.i && m.j >= 0 && m.j < dims_.j) {
      return m.i * dims_.j + m.j;
    }
    return kBlockNone;
  }
  void AddParticles(const BlockData& other) {
    const size_t size = other.position.size();
    assert(other.velocity.size() == size);
    assert(other.id.size() == size);
    for (size_t k = 0; k < size; ++k) {
      for (size_t p = 0; p < other.position[k].size(); ++p) {
        size_t i = FindBlock(other.position[k][p]);
        if (i != kBlockNone) {
          data_.AddParticle(
              i, other.position[k][p], other.velocity[k][p], other.id[k][p]);
        }
      }
    }
    SortParticles();
  }
  template <class VectorVect, class VectorInt>
  void AddParticles(
      const VectorVect& position, const VectorVect& velocity,
      const VectorInt& id) {
    const size_t size = position.size();
    assert(velocity.size() == size);
    assert(id.size() == size);

    for (size_t p = 0; p < size; ++p) {
      size_t i = FindBlock(position[p]);
      if (i != kBlockNone) {
        data_.AddParticle(i, position[p], velocity[p], id[p]);
      }
    }
    SortParticles();
  }
  const BlockData& GetData() const {
    return data_;
  }
  BlockData& GetData() {
    return data_;
  }
  const std::vector<std::pair<size_t, size_t>> GetBlockById() const {
    return block_by_id_;
  }
  size_t GetNumBlocks() const {
    return num_blocks_;
  }
  size_t GetNumParticles() const {
    return num_particles_;
  }
  size_t GetNumPerCell() const {
    return num_per_cell_;
  }
  std::array<int, kNumNeighbors> GetNeighborOffsets() const {
    return neighbor_offsets_;
  }
  void SortParticles() {
    size_t lnum_particles_ = 0;
    size_t max_per_cell = 0;
    for (size_t i = 0; i < num_blocks_; ++i) {
      size_t p = 0;
      size_t pe = data_.block_size[i];
      while (p < pe) {
        const size_t j = FindBlock(data_.position[i][p]);
        if (i != j) {
          if (j != kBlockNone) {
            data_.MoveParticle(i, p, j);
            // if the new block is already processed, increase the counter
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
      max_per_cell = std::max(max_per_cell, pe);
    }

    num_particles_ = lnum_particles_;
    num_per_cell_ = max_per_cell;
  }

 private:
  RectVect domain_;
  Vect block_size_;
  BlockData data_;
  MIdx dims_;
  size_t num_blocks_;
  std::array<int, kNumNeighbors> neighbor_offsets_;
  std::vector<std::pair<size_t, size_t>> block_by_id_;
  size_t num_particles_;
  size_t num_per_cell_;
  void InitEmptyBlocks(RectVect domain, Vect block_size) {
    domain_ = domain;
    block_size_ = block_size;
    dims_.i = int(domain.size().x / block_size.x) + 1;
    dims_.j = int(domain.size().y / block_size.y) + 1;
    num_blocks_ = static_cast<size_t>(dims_.i * dims_.j);

    assert(dims_.i > 0 && dims_.j > 0 && num_blocks_ > 0);

    data_.Clear();
    data_.Resize(num_blocks_);

    // Calc offsets to neighbors
    size_t n = 0;
    for (int i = -1; i <= 1; ++i) {
      for (int j = -1; j <= 1; ++j) {
        neighbor_offsets_[n] = i * dims_.j + j;
        ++n;
      }
    }
  }
};
