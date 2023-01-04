#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <vector>
#include "aligned_allocator.hpp"
#include "geometry.hpp"

using ArrayVect = std::vector<Vect, AlignedAllocator<Vect, 64>>;
using ArrayInt = std::vector<int, AlignedAllocator<int, 64>>;

class blocks {
 public:
  static const size_t kNumNeighbors = 9;
  static const size_t kBlockNone = static_cast<size_t>(-1);
  using DataVect = std::vector<ArrayVect>;
  using DataInt = std::vector<ArrayInt>;
  struct BlockData {
   private:
    blocks* parent;

   public:
    DataVect position, position_tmp, velocity, velocity_tmp, force;
    DataInt id;
    BlockData() = delete;
    BlockData(blocks* parent_) : parent(parent_) {}
    void clear() {
      position.clear();
      position_tmp.clear();
      velocity.clear();
      velocity_tmp.clear();
      force.clear();
      id.clear();

      for (auto& pair : parent->block_by_id_) {
        pair.first = kBlockNone;
      }
    }
    void resize(size_t size) {
      position.resize(size);
      position_tmp.resize(size);
      velocity.resize(size);
      velocity_tmp.resize(size);
      force.resize(size);
      id.resize(size);

      for (size_t i = 0; i < position.size(); ++i) {
        const size_t kBlockPadding = 16;
        position[i].reserve(kBlockPadding);
        position_tmp[i].reserve(kBlockPadding);
        velocity[i].reserve(kBlockPadding);
        velocity_tmp[i].reserve(kBlockPadding);
        force[i].reserve(kBlockPadding);
        id[i].reserve(kBlockPadding);
      }

      // TODO: consider updating block_by_id_ here
      for (auto& pair : parent->block_by_id_) {
        pair.first = kBlockNone;
      }
    }
    void RemoveParticle(
        size_t src, // source block
        size_t idx // particle index within the source block
    ) {
      std::swap(position[src][idx], position[src].back());
      std::swap(position_tmp[src][idx], position_tmp[src].back());
      std::swap(velocity[src][idx], velocity[src].back());
      std::swap(velocity_tmp[src][idx], velocity_tmp[src].back());
      std::swap(force[src][idx], force[src].back());
      std::swap(id[src][idx], id[src].back());

      parent->block_by_id_[id[src][idx]] = {src, idx};
      parent->block_by_id_[id[src].back()].first = kBlockNone;

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
      assert(src != dest);
      position[dest].push_back(position[src][idx]);
      position_tmp[dest].push_back(position_tmp[src][idx]);
      velocity[dest].push_back(velocity[src][idx]);
      velocity_tmp[dest].push_back(velocity_tmp[src][idx]);
      force[dest].push_back(force[src][idx]);
      id[dest].push_back(id[src][idx]);

      RemoveParticle(src, idx);

      parent->block_by_id_[id[dest].back()] = {dest, position[dest].size() - 1};
    }
    void AddParticle(
        size_t dest, // destination block
        Vect particle_position, Vect particle_velocity, int particle_id) {
      position[dest].push_back(particle_position);
      position_tmp[dest].push_back(Vect::kNan);
      velocity[dest].push_back(particle_velocity);
      velocity_tmp[dest].push_back(Vect::kNan);
      force[dest].push_back(Vect::kNan);
      id[dest].push_back(particle_id);

      assert(id[dest].back() >= 0);
      const size_t pid = static_cast<size_t>(id[dest].back());
      if (parent->block_by_id_.size() <= pid) {
        parent->block_by_id_.resize(pid + 1);
      }
      parent->block_by_id_[pid] = {dest, position[dest].size() - 1};
    }
  };
  void SetDomain(rect_vect proposal) {
    rect_vect domain = domain_;
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
  // TODO: draw blocks
  Vect GetCenter(size_t block) const {
    const mindex m(block / dims_.j, block % dims_.j);
    return domain_.A +
           Vect((0.5 + m.i) * block_size_.x, (0.5 + m.j) * block_size_.y);
  }
  blocks(rect_vect domain, Vect block_size)
      : data_(this), num_particles_(0), num_per_cell_(0) {
    InitEmptyBlocks(domain, block_size);
  }
  size_t FindBlock(Vect position) const {
    mindex m(
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
  void AddParticles(
      const ArrayVect& position, const ArrayVect& velocity,
      const std::vector<int>& id) {
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
      size_t pe = data_.position[i].size();
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
  rect_vect domain_;
  Vect block_size_;
  BlockData data_;
  mindex dims_;
  size_t num_blocks_;
  std::array<int, kNumNeighbors> neighbor_offsets_;
  std::vector<std::pair<size_t, size_t>> block_by_id_;
  size_t num_particles_;
  size_t num_per_cell_;
  void InitEmptyBlocks(rect_vect domain, Vect block_size) {
    domain_ = domain;
    block_size_ = block_size;
    dims_.i = int(domain.size().x / block_size.x) + 1;
    dims_.j = int(domain.size().y / block_size.y) + 1;
    num_blocks_ = static_cast<size_t>(dims_.i * dims_.j);

    assert(dims_.i > 0 && dims_.j > 0 && num_blocks_ > 0);

    data_.clear();
    data_.resize(num_blocks_);

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
