#include <chrono>
#include <cstdint>
#include <limits>
#include <thread>

#include "macros.h"
#include "particles.h"

const Scal kRadius = 0.02;
const Scal kSigma = 1;
const Scal kSigmaWall = 1;
const Scal kSigmaBond = 1e5;
const Scal kSigmaPick = 1e3;
const Scal kSigmaPortalEdge = 1;
const Scal kRadiusPortalEdge = 3 * kRadius;
const Scal kMass = kRadius * kRadius * 100;
const Scal kPointForce = 0.1;
const Scal kPointForceAttractive = 0.1;
const Scal kDissipation = 0.001;
// Damping of the normal relative velocity inside a contact. Unlike
// kDissipation it only acts between overlapping particles and only on their
// approach, so it leaves bulk motion alone. See docs/MODEL.md.
const Scal kDashpot = 100;
const Scal kBlockSize = 4. * kRadius;
const Scal kGravity = 10;
const Scal kPortalThickness = 0.02;
const Scal kVelocityLimit = 10;

const int kParticleIdNone = -1;
const Scal kTimeStep = 0.0005;

Particles::Particles()
    : domain(RectVect(Vect(-1, -1), Vect(1, 1)))
    , Blocks(domain, Vect(kBlockSize, kBlockSize))
    , blocks_buffer_(Blocks) {
  force_enabled = false;
  force_center = Vect(0, 0);
  remove_last_portal_ = false;

  t = 0.0;
  dt = kTimeStep;
  gravity_ = Vect(0, -1) * kGravity;

  SetDomain(domain);
  resize_queue_ = domain;
  ResetEnvObjFrame(domain);

  SetParticleGrid(kGridMedium);
}

// Bonds, frozen particles and portals refer to particles by id, so they cannot
// outlive the particles they were made of and are cleared here. The two
// portals flanking the grid are recreated, as in the initial state. The domain
// is left alone since it follows the window size.
void Particles::SetParticleGrid(size_t size) {
  grid_size_ = size;

  bonds_.clear();
  no_rendering_.clear();
  no_rendering_buffer_.clear();
  frozen_.clear();
  portals_.clear();
  particle_to_move_.clear();
  bonds_prev_particle_id_ = kParticleIdNone;
  bonds_enabled_ = false;
  freeze_enabled_ = false;
  pick_enabled_ = false;
  portal_enabled_ = false;
  portal_mouse_moving_ = false;
  portal_stage_ = 0;
  t = 0.;

  auto& data = Blocks.GetData();
  data.clear();
  data.resize(Blocks.GetNumBlocks());

  // Hexagonal packing resting on the bottom of the domain.
  const Scal width = size * 2. * kRadius;
  const Scal height = size * std::sqrt(3.) * kRadius;
  const RectVect box(Vect(-0.5 * width, -1.), Vect(0.5 * width, -1. + height));
  ArrayVect position;
  ArrayVect velocity;
  std::vector<int> id;
  for (size_t j = 0; j < size; ++j) {
    for (size_t i = 0; i < size; ++i) {
      position.push_back(Vect(
          box.A.x + kRadius * (2. * i + 1. + (j % 2)),
          box.A.y + kRadius * (std::sqrt(3.) * j + 1.)));
      velocity.push_back(Vect(0.));
      id.push_back(id.size());
    }
  }
  Blocks.AddParticles(position, velocity, id);
  SetParticleBuffer();

  const Scal dx = kPortalThickness;
  PortalStart(Vect(box.A.x - dx, box.A.y));
  PortalStop(Vect(box.A.x - dx, box.B.y + 0.2));
  PortalStart(Vect(box.B.x + dx, box.A.y));
  PortalStop(Vect(box.B.x + dx, box.B.y + 0.2));

  std::cout << "Particle grid " << size << "x" << size << " = "
            << size * size << " particles" << std::endl;
}

void Particles::AddParticleBlock(Vect relative, size_t size) {
  const Vect center = domain.A + domain.size() * relative;
  // Hexagonal packing, as in SetParticleGrid().
  const Scal width = size * 2. * kRadius;
  const Scal height = size * std::sqrt(3.) * kRadius;
  const Vect low = center - Vect(width, height) * 0.5;

  // Ids are never reused, so the first free one is past the last known.
  int next_id = static_cast<int>(Blocks.GetBlockById().size());
  ArrayVect position;
  ArrayVect velocity;
  std::vector<int> id;
  for (size_t j = 0; j < size; ++j) {
    for (size_t i = 0; i < size; ++i) {
      position.push_back(Vect(
          low.x + kRadius * (2. * i + 1. + (j % 2)),
          low.y + kRadius * (std::sqrt(3.) * j + 1.)));
      velocity.push_back(Vect(0.));
      id.push_back(next_id);
      ++next_id;
    }
  }
  Blocks.AddParticles(position, velocity, id);
  SetParticleBuffer();

  std::cout << "Added " << size << "x" << size << " particles, "
            << Blocks.GetNumParticles() << " in total" << std::endl;
}

Particles::~Particles() {}
void Particles::SetParticleBuffer() {
  blocks_buffer_ = Blocks;
  std::vector<particle> res;
  for (size_t iblock = 0; iblock < blocks_buffer_.GetNumBlocks(); ++iblock) {
    const auto& data = blocks_buffer_.GetData();
    for (size_t p = 0; p < data.position[iblock].size(); ++p) {
      res.emplace_back(data.position[iblock][p], data.velocity[iblock][p]);
    }
  }
  particle_buffer_ = res;
}
void Particles::AddEnvObj(env_object* env) {
  ENVOBJ.push_back(std::unique_ptr<env_object>(env));
}
void Particles::step(Scal time_target, bool quit) {
#pragma omp parallel
  {
    while (t < time_target && !quit) {
#pragma omp for schedule(dynamic, 8)
      for (size_t iblock = 0; iblock < Blocks.GetNumBlocks(); ++iblock) {
        calc_forces(iblock);
      }

#pragma omp single
      RHS_bonds();

      // Parallelized inside, so called by all threads.
      ApplyPortalsForces();

#pragma omp single
      ApplyFrozen();

#pragma omp for schedule(dynamic, 8)
      for (size_t iblock = 0; iblock < Blocks.GetNumBlocks(); ++iblock) {
        auto& data = Blocks.GetData();
        for (size_t p = 0; p < data.position[iblock].size(); ++p) {
          data.velocity_tmp[iblock][p] = data.velocity[iblock][p];
          data.position_tmp[iblock][p] = data.position[iblock][p];
          data.velocity[iblock][p] +=
              data.force[iblock][p] * (dt * 0.5 / kMass);
          if (data.velocity[iblock][p].length() > kVelocityLimit) {
            data.velocity[iblock][p] *=
                kVelocityLimit / data.velocity[iblock][p].length();
          }
          data.position[iblock][p] += data.velocity[iblock][p] * dt * 0.5;
        }
      }

#pragma omp single
      t += 0.5 * dt;

#pragma omp for schedule(dynamic, 8)
      for (size_t iblock = 0; iblock < Blocks.GetNumBlocks(); ++iblock) {
        calc_forces(iblock);
      }

#pragma omp single
      RHS_bonds();

      // Parallelized inside, so called by all threads.
      ApplyPortalsForces();

#pragma omp single
      ApplyFrozen();

#pragma omp for schedule(dynamic, 8)
      for (size_t iblock = 0; iblock < Blocks.GetNumBlocks(); ++iblock) {
        auto& data = Blocks.GetData();
        for (size_t p = 0; p < data.position[iblock].size(); ++p) {
          data.velocity[iblock][p] = data.velocity_tmp[iblock][p] +
                                     data.force[iblock][p] * (dt / kMass);
          if (data.velocity[iblock][p].length() > kVelocityLimit) {
            data.velocity[iblock][p] *=
                kVelocityLimit / data.velocity[iblock][p].length();
          }
          data.position[iblock][p] =
              data.position_tmp[iblock][p] +
              (data.velocity_tmp[iblock][p] + data.velocity[iblock][p]) *
                  (dt * 0.5);
        }
      }

#pragma omp single
      {
        // Resize the frame if needed (with a limited speed)
        auto limit = [](Scal& current, const Scal target) {
          const Scal limit = 0.01;
          current =
              std::min(current + limit, std::max(current - limit, target));
        };

        RectVect new_domain = domain;
        limit(new_domain.A.x, resize_queue_.A.x);
        limit(new_domain.A.y, resize_queue_.A.y);
        limit(new_domain.B.x, resize_queue_.B.x);
        limit(new_domain.B.y, resize_queue_.B.y);

        if (int(t / dt) % int(0.01 / dt) == 0) {
          if (new_domain != domain) {
            SetDomain(new_domain);
            ResetEnvObjFrame(new_domain);
          }
        }

        DetectPortals();
        ApplyPortals();
      }

      // Re-bin the particles that left their block. The scan is parallel,
      // the moves that follow are not, see blocks::ScanBlock().
#pragma omp for schedule(static)
      for (size_t iblock = 0; iblock < Blocks.GetNumBlocks(); ++iblock) {
        Blocks.ScanBlock(iblock);
      }

#pragma omp single
      {
        Blocks.ApplyMoves();
        CheckBonds();
        CheckFrozen();

        if (remove_last_portal_) {
          if (portals_.size()) {
            portals_.pop_back();
          }
          remove_last_portal_ = false;
        }

        // Pass the data to renderer if ready
        if (renderer_ready_for_next_) {
          SetParticleBuffer();
          RHS_bonds();
          no_rendering_buffer_ = no_rendering_;
          renderer_ready_for_next_ = false;
          // std::this_thread::sleep_for(std::chrono::milliseconds(25));
        }

        // Advance in time (another half)
        t += 0.5 * dt;
      }
    }
  }
}

void Particles::CheckBonds() {
  const auto& bbi = Blocks.GetBlockById();
  for (auto it = bonds_.begin(); it != bonds_.end();) {
    if (bbi[it->first].first == blocks::kBlockNone ||
        bbi[it->second].first == blocks::kBlockNone) {
      it = bonds_.erase(it);
    } else {
      ++it;
    }
  }
}

// Particles that left the domain are removed by blocks::SortParticles(), which
// only marks them in block_by_id_. Forget them here, otherwise ApplyFrozen()
// and the renderer would look up a block that does not exist.
void Particles::CheckFrozen() {
  const auto& bbi = Blocks.GetBlockById();
  for (auto it = frozen_.begin(); it != frozen_.end();) {
    if (bbi[*it].first == blocks::kBlockNone) {
      it = frozen_.erase(it);
    } else {
      ++it;
    }
  }
}

void Particles::SetForce(Vect center, bool enabled) {
  force_center = center;
  force_enabled = enabled;
}
void Particles::SetForce(Vect center) {
  force_center = center;
}
void Particles::SetForce(bool enabled) {
  force_enabled = enabled;
}

// One step changes gravity by half of its default magnitude, which reaches
// zero in two steps and reverses in four.
const Scal kGravityStep = kGravity * 0.5;
const Scal kGravityMax = kGravity * 2;

void Particles::ChangeGravity(int steps) {
  gravity_.y = std::min(
      kGravityMax, std::max(-kGravityMax, gravity_.y + steps * kGravityStep));
  // Changing gravity while it is off would do nothing visible.
  gravity_enable_ = true;
  const Scal down = -gravity_.y;
  const char* direction = down > 0 ? " down" : down < 0 ? " up" : "";
  std::cout << "Gravity: " << std::abs(down) << direction << std::endl;
}

void Particles::PickStart(Vect point) {
  size_t min_block = blocks::kBlockNone;
  size_t min_particle = 0;
  Scal min_dist;

  auto& data = blocks_buffer_.GetData();
  for (size_t iblock = 0; iblock < blocks_buffer_.GetNumBlocks(); ++iblock) {
    for (size_t p = 0; p < data.position[iblock].size(); ++p) {
      if (min_block == blocks::kBlockNone ||
          data.position[iblock][p].dist(point) < min_dist) {
        min_block = iblock;
        min_particle = p;
        min_dist = data.position[iblock][p].dist(point);
      }
    }
  }

  pick_particle_id_ = data.id[min_block][min_particle];
  pick_pointer_ = point;
  pick_enabled_ = true;
}

void Particles::PickMove(Vect point) {
  if (!pick_enabled_) {
    return;
  }
  pick_pointer_ = point;
}

void Particles::PickStop(Vect) {
  if (!pick_enabled_) {
    return;
  }
  pick_enabled_ = false;
}

void Particles::MoveToPortal(
    Vect& position, Vect& velocity, const Portal& src, const Portal& dest) {
  const Vect src_a = src.begin;
  const Vect src_b = src.end;
  const Vect src_r = src_b - src_a;
  const Vect src_n = Vect(-src_r.y, src_r.x).GetNormalized();

  const Vect dest_a = dest.begin;
  const Vect dest_b = dest.end;
  const Vect dest_r = dest_b - dest_a;
  const Vect dest_n = Vect(-dest_r.y, dest_r.x).GetNormalized();

  const Scal lambda_pos = src_r.dot(position - src_a) / src_r.dot(src_r);
  const Scal offset_pos = src_n.dot(position - src_a);
  const Scal lambda_vel = src_r.dot(velocity) / src_r.dot(src_r);
  const Scal offset_vel = src_n.dot(velocity);

  const Scal sign = (offset_pos > 0. ? 1. : -1.);

  position = dest_a + dest_r * lambda_pos +
             dest_n * (offset_pos - sign * 2. * kPortalThickness);
  velocity = dest_r * lambda_vel + dest_n * offset_vel;
}

void Particles::DetectPortals() {
  for (auto& pair : portals_) {
    for (int d = 0; d <= 1; ++d) {
      const auto& portal = pair[d];

      const Vect a = portal.begin;
      const Vect b = portal.end;
      const Vect r = b - a;
      const Vect n = Vect(-r.y, r.x).GetNormalized();

      auto& data = Blocks.GetData();
      for (size_t iblock : portal.blocks) {
        for (size_t p = 0; p < data.position[iblock].size(); ++p) {
          auto id = static_cast<size_t>(data.id[iblock][p]);
          if (particle_to_move_.size() <= id) {
            particle_to_move_.resize(id + 1);
          }
          const Vect curr = data.position[iblock][p];
          const Scal lambda_curr = (curr - a).dot(r) / r.dot(r);
          const Scal offset_curr = (curr - a).dot(n);
          if (lambda_curr > 0. && lambda_curr < 1. &&
              std::abs(offset_curr) <= kPortalThickness) {
            particle_to_move_[id] = 1;
          }
        }
      }
    }
  }
}

void Particles::UpdatePortalCache(const Portal& portal, PortalCache& cache) {
  const Vect a = portal.begin;
  const Vect r = portal.end - a;
  const Vect n = Vect(-r.y, r.x).GetNormalized();
  const auto& data = Blocks.GetData();

  cache.block.clear();
  cache.index.clear();
  cache.position.clear();
  cache.lambda.clear();
  cache.offset.clear();
  cache.groups.clear();
  for (size_t iblock : portal.blocks) {
    if (data.position[iblock].empty()) {
      continue;
    }
    const Scal inf = std::numeric_limits<Scal>::max();
    PortalCache::Group group;
    group.begin = cache.position.size();
    group.low = Vect(inf);
    group.high = Vect(-inf);
    group.lambda_min = inf;
    group.lambda_max = -inf;
    for (size_t p = 0; p < data.position[iblock].size(); ++p) {
      const Vect x = data.position[iblock][p];
      const Scal lambda = r.dot(x - a) / r.dot(r);
      cache.block.push_back(iblock);
      cache.index.push_back(p);
      cache.position.push_back(x);
      cache.lambda.push_back(lambda);
      cache.offset.push_back((x - a).dot(n));
      group.low = Vect(std::min(group.low.x, x.x), std::min(group.low.y, x.y));
      group.high =
          Vect(std::max(group.high.x, x.x), std::max(group.high.y, x.y));
      group.lambda_min = std::min(group.lambda_min, lambda);
      group.lambda_max = std::max(group.lambda_max, lambda);
    }
    group.end = cache.position.size();
    cache.groups.push_back(group);
  }
}

void Particles::ApplyPortalsForces() {
  // Assume that the particles have just been moved
  // with their velocity
  // so that (position - dt * velocity) is the previous position

  // Distance beyond which F12() returns zero.
  const Scal kCutoff = 2. * kRadius;

  // Called from the parallel region of step(), so the directives below are
  // orphaned and bind to it. Every particle only gets a force of its own, so
  // the loop over them is shared between the threads.
  auto& data = Blocks.GetData();
  for (auto& pair : portals_) {
    // The coordinates relative to the portal are the same for every particle
    // of the pair, so they are computed once instead of once per pair of
    // particles below.
#pragma omp single
    {
      UpdatePortalCache(pair[0], portal_cache_[0]);
      UpdatePortalCache(pair[1], portal_cache_[1]);
    }

    for (int d = 0; d <= 1; ++d) {
      auto& portal = pair[d];
      const PortalCache& cache = portal_cache_[d];
      const PortalCache& other_cache = portal_cache_[1 - d];
      const Vect a = portal.begin;
      const Vect b = portal.end;

      const Vect r = b - a;
      const Vect n = Vect(-r.y, r.x).GetNormalized();
      const Scal r_length = r.length();
#pragma omp for schedule(static)
      for (size_t p = 0; p < cache.position.size(); ++p) {
        const Vect curr = cache.position[p];
        const Scal lambda_curr = cache.lambda[p];
        const Scal offset_curr = cache.offset[p];
        if (!(lambda_curr > 0. && lambda_curr < 1.)) {
          continue;
        }
        Vect& force = data.force[cache.block[p]][cache.index[p]];

        // Check particle forces
        for (const auto& group : cache.groups) {
          if (curr.x < group.low.x - kCutoff ||
              curr.x > group.high.x + kCutoff ||
              curr.y < group.low.y - kCutoff ||
              curr.y > group.high.y + kCutoff) {
            continue;
          }
          for (size_t q = group.begin; q < group.end; ++q) {
            const Vect neighbor = cache.position[q];
            const Vect dp = neighbor - curr;
            if (dp.dot(dp) > kCutoff * kCutoff) {
              continue;
            }
            const Scal lambda_neighbor = cache.lambda[q];
            const Scal offset_neighbor = cache.offset[q];
            if (lambda_neighbor > 0. && lambda_neighbor < 1. &&
                offset_neighbor * offset_curr < 0.) {
              force -= F12(curr, neighbor);
            }
          }
        }

        // The projection of a particle of the other portal lies at its
        // position along that portal, so a larger gap along the portal alone
        // already puts the projection beyond the cutoff.
        const Scal window = kCutoff / r_length;
        for (const auto& group : other_cache.groups) {
          if (group.lambda_max < lambda_curr - window ||
              group.lambda_min > lambda_curr + window) {
            continue;
          }
          for (size_t q = group.begin; q < group.end; ++q) {
            const Scal other_lambda_neighbor = other_cache.lambda[q];
            if (std::abs(other_lambda_neighbor - lambda_curr) > window) {
              continue;
            }
            const Scal other_offset_neighbor = other_cache.offset[q];
            if (other_lambda_neighbor > 0. && other_lambda_neighbor < 1. &&
                other_offset_neighbor * offset_curr < 0.) {
              const Scal sign = (offset_curr > 0. ? 1. : -1.);
              const Vect proj =
                  a + r * other_lambda_neighbor +
                  n * (other_offset_neighbor + 2. * sign * kPortalThickness);
              force += F12(curr, proj);
            }
          }
        }
      }
    }
  }
}
void Particles::ApplyPortals() {
  // Assume that the particles have just been moved
  // with their velocity
  // so that (position - dt * velocity) is the previous position

  auto& data = Blocks.GetData();
  for (auto& pair : portals_) {
    for (int d = 0; d <= 1; ++d) {
      auto& portal = pair[d];
      auto& other = pair[1 - d];

      for (size_t iblock : portal.blocks) {
        for (size_t p = 0; p < data.position[iblock].size(); ++p) {
          const auto id = static_cast<size_t>(data.id[iblock][p]);
          if (particle_to_move_.size() <= id) {
            particle_to_move_.resize(id + 1);
          }
          if (particle_to_move_[id]) {
            MoveToPortal(
                data.position[iblock][p], data.velocity[iblock][p], portal,
                other);
            particle_to_move_[id] = 0;
          }
        }
      }
    }
  }
}

void Particles::PortalStart(Vect point) {
  portal_enabled_ = true;
  portal_begin_ = point;
  portal_current_ = point;
  portal_mouse_moving_ = true;
}

void Particles::PortalMove(Vect point) {
  if (!portal_enabled_) {
    return;
  }
  portal_current_ = point;
}

void Particles::PortalStop(Vect point) {
  if (!portal_enabled_) {
    return;
  }
  portal_enabled_ = false;
  portal_mouse_moving_ = false;
  if (portal_stage_ == 0) {
    portal_prev_.first = portal_begin_;
    portal_prev_.second = point;
    portal_stage_ = 1;
  } else {
    std::array<Portal, 2> pair;
    pair[0].begin = portal_prev_.first;
    pair[0].end = portal_prev_.second;
    pair[1].begin = portal_begin_;
    pair[1].end = point;
    pair[1].end =
        pair[1].begin + (pair[1].end - pair[1].begin).GetNormalized() *
                            pair[0].begin.dist(pair[0].end);
    portals_.push_back(pair);
    portal_stage_ = 0;

    UpdatePortalBlocks(portals_.back()[0]);
    UpdatePortalBlocks(portals_.back()[1]);
  }
}

void Particles::BondsStart(Vect point) {
  bonds_enabled_ = true;
  int id = kParticleIdNone;

  for (size_t iblock = 0; iblock < blocks_buffer_.GetNumBlocks(); ++iblock) {
    auto& data = blocks_buffer_.GetData();
    for (size_t p = 0; p < data.position[iblock].size(); ++p) {
      if (data.position[iblock][p].dist(point) < kRadius) {
        id = data.id[iblock][p];
      }
    }
  }

  bonds_prev_particle_id_ = id;
}

void Particles::BondsMove(Vect point) {
  if (!bonds_enabled_) {
    return;
  }
  int id = kParticleIdNone;

  for (size_t iblock = 0; iblock < blocks_buffer_.GetNumBlocks(); ++iblock) {
    auto& data = blocks_buffer_.GetData();
    for (size_t p = 0; p < data.position[iblock].size(); ++p) {
      if (data.position[iblock][p].dist(point) < kRadius &&
          data.id[iblock][p] != bonds_prev_particle_id_) {
        id = data.id[iblock][p];
      }
    }
  }

  if (bonds_prev_particle_id_ != kParticleIdNone && id != kParticleIdNone) {
    assert(bonds_prev_particle_id_ != id);
    std::pair<int, int> bond(bonds_prev_particle_id_, id);
    std::pair<int, int> bond_r(id, bonds_prev_particle_id_);

    if (bonds_.count(bond)) {
      bonds_.erase(bond);
      bonds_.erase(bond_r);
    } else {
      bonds_.insert(bond);
      bonds_.insert(bond_r);
    }
  }

  if (id != kParticleIdNone) {
    bonds_prev_particle_id_ = id;
  }
}

void Particles::BondsStop(Vect) {
  if (!bonds_enabled_) {
    return;
  }
  bonds_enabled_ = false;
}

void Particles::FreezeStart(Vect point) {
  freeze_last_id_ = -1;
  freeze_enabled_ = true;
  FreezeMove(point);
}

void Particles::FreezeMove(Vect mousepos) {
  if (!freeze_enabled_) {
    return;
  }

  const auto kFreezeRadius = kRadius * 3;
  size_t min_iblock = -1;
  size_t min_ip = -1;
  Scal min_dist = kFreezeRadius;

  auto& data = blocks_buffer_.GetData();
  // Find the nearest particle within the radius.
  for (size_t iblock = 0; iblock < blocks_buffer_.GetNumBlocks(); ++iblock) {
    for (size_t ip = 0; ip < data.position[iblock].size(); ++ip) {
      const Scal dist = data.position[iblock][ip].dist(mousepos);
      if (dist < min_dist) {
        min_dist = dist;
        min_iblock = iblock;
        min_ip = ip;
      }
    }
  }
  // Freeze or unfreeze the particle if diffferent from previous.
  if (min_dist < kFreezeRadius) {
    const int id = data.id[min_iblock][min_ip];
    if (id != freeze_last_id_) {
      const auto it = frozen_.find(id);
      if (it == frozen_.end()) {
        frozen_.insert(it, id);
        std::cout << "Freeze particle id=" << id << std::endl;
      } else {
        frozen_.erase(id);
        std::cout << "Unfreeze particle id=" << id << std::endl;
      }
      freeze_last_id_ = id;
    }
  }
}

void Particles::FreezeStop(Vect) {
  if (!freeze_enabled_) {
    return;
  }
  freeze_enabled_ = false;
}

Vect F12(Vect p1, Vect p2, Scal R, Scal sigma) {
  const Vect dp = p1 - p2;
  const Scal r2 = dp.dot(dp);
  const Scal r2inv = 1. / r2;
  const Scal d2 = r2inv * (R * R);
  const Scal d6 = d2 * d2 * d2;
  const Scal d12 = d6 * d6;
  return dp * std::max<Scal>(0., sigma * (d12 - d6) * r2inv);
}

Vect F12wall(Vect p1, Vect p2) {
  const Scal sigma = kSigmaWall;
  const Scal R = kRadius;
  const Vect dp = p1 - p2;
  const Scal r2 = dp.dot(dp);
  const Scal r2inv = 1. / r2;
  const Scal d2 = r2inv * (R * R);
  const Scal d6 = d2 * d2 * d2;
  const Scal d12 = d6 * d6;
  return dp * std::max<Scal>(0., sigma * (d12 - d6) * r2inv);
}

Vect F12(Vect p1, Vect p2) {
  const Scal threshold = std::pow(kRadius, 2) * 1e-3;
  const Scal sigma = kSigma;
  const Scal R = 2. * kRadius;
  const Vect dp = p1 - p2;
  const Scal r2 = std::max(threshold, dp.dot(dp));
  const Scal r2inv = 1. / r2;
  const Scal d2 = r2inv * (R * R);
  const Scal d6 = d2 * d2 * d2;
  const Scal d12 = d6 * d6;
  return dp * std::max<Scal>(0., sigma * (d12 - d6) * r2inv);
}

// Pair force with the contact dashpot. The damping acts along the line of
// centers, so it does not resist shear, and it is scaled by the overlap so
// that it vanishes together with the contact instead of jumping to zero.
Vect F12(Vect p1, Vect p2, Vect v1, Vect v2) {
  const Scal threshold = std::pow(kRadius, 2) * 1e-3;
  const Scal sigma = kSigma;
  const Scal R = 2. * kRadius;
  const Vect dp = p1 - p2;
  const Scal r2 = std::max(threshold, dp.dot(dp));
  const Scal r2inv = 1. / r2;
  const Scal d2 = r2inv * (R * R);
  const Scal d6 = d2 * d2 * d2;
  const Scal d12 = d6 * d6;
  const Scal spring = std::max<Scal>(0., sigma * (d12 - d6) * r2inv);
  const Scal overlap = std::max<Scal>(0., 1. - r2 / (R * R));
  const Scal damping = -kDashpot * overlap * (v1 - v2).dot(dp) * r2inv;
  return dp * (spring + damping);
}

template <bool ApplyThreshold = true>
void CalcForceSerial(
    ArrayVect& force, ArrayVect& position, ArrayVect& position_other,
    ArrayVect& velocity, ArrayVect& velocity_other) {
  for (size_t q = 0; q < position_other.size(); ++q) {
    for (size_t p = 0; p < position.size(); ++p) {
      if (&position[p] != &position_other[q])
        force[p] +=
            F12(position[p], position_other[q], velocity[p], velocity_other[q]);
    }
  }
}

template <bool ApplyThreshold = true>
void CalcForceSerialPadded(
    ArrayVect& force, ArrayVect& position, ArrayVect& position_other,
    ArrayVect& velocity, ArrayVect& velocity_other) {
  for (size_t q = 0; q < position_other.size(); ++q) {
    for (size_t p = 0; p < position.size(); p += 8) {
      for (size_t k = 0; k < 8; ++k) {
        force[p + k] +=
            F12(position[p + k], position_other[q], velocity[p + k],
                velocity_other[q]);
      }
    }
  }
}

#if USEFLAG(AVX)
#include <x86intrin.h>
#define CALC_FORCE CalcForceAvx
// Accesses `force` and `position` in groups of 8 particles, up to the next
// multiple of 8 past their size. This stays inside the allocation since
// `blocks` keeps the capacity of every block padded to a multiple of
// blocks::kVectorWidth with the padding elements initialized.
template <bool ApplyThreshold = true>
void CalcForceAvx(
    ArrayVect& force, ArrayVect& position, ArrayVect& position_other,
    ArrayVect& velocity, ArrayVect& velocity_other) {
  static_assert(blocks::kVectorWidth == 8, "kernel processes 8 particles");
  assert(position.capacity() % blocks::kVectorWidth == 0);
  assert(force.capacity() >= position.capacity());
  assert(velocity.capacity() >= position.capacity());
  assert(reinterpret_cast<uintptr_t>(position.data()) % 32 == 0);
  assert(reinterpret_cast<uintptr_t>(force.data()) % 32 == 0);
  assert(reinterpret_cast<uintptr_t>(velocity.data()) % 32 == 0);

  // sigma = kSigma;
  const __m256 sigma = _mm256_broadcast_ss(&kSigma);
  // R2 = (2. * kRadius) ^ 2;
  const float tmp = std::pow(2. * kRadius, 2);
  const __m256 R2 = _mm256_broadcast_ss(&tmp);
  // threshold = kRadius ^ 2 * 1e-3
  const float tmp_th = std::pow(kRadius, 2) * 1e-3;
  const __m256 threshold = _mm256_broadcast_ss(&tmp_th);
  const float tmp_zero = 0.;
  const __m256 zero = _mm256_broadcast_ss(&tmp_zero);
  // one = 1, used to turn r2 / R2 into the overlap.
  const float tmp_one = 1.;
  const __m256 one = _mm256_broadcast_ss(&tmp_one);
  // R2inv = 1 / (2. * kRadius) ^ 2
  const float tmp_r2inv = 1. / tmp;
  const __m256 R2inv = _mm256_broadcast_ss(&tmp_r2inv);
  // dashpot = -kDashpot, folded into the sign of the damping term
  const float tmp_dashpot = -kDashpot;
  const __m256 dashpot = _mm256_broadcast_ss(&tmp_dashpot);

  // Padding elements are accessed through the data pointers, past the size.
  const float* position_data = (const float*)position.data();
  const float* velocity_data = (const float*)velocity.data();
  float* force_data = (float*)force.data();

  for (size_t q = 0; q < position_other.size(); ++q) {
    const __m256 qx = _mm256_broadcast_ss((float*)&position_other[q].x);
    const __m256 qy = _mm256_broadcast_ss((float*)&position_other[q].y);
    // qxy = (q.x, q.y)
    const __m256 qxy = _mm256_blend_ps(qx, qy, 0xAA);
    const __m256 qvx = _mm256_broadcast_ss((float*)&velocity_other[q].x);
    const __m256 qvy = _mm256_broadcast_ss((float*)&velocity_other[q].y);
    // qvxy = (q.v.x, q.v.y)
    const __m256 qvxy = _mm256_blend_ps(qvx, qvy, 0xAA);
    for (size_t p = 0; p < position.size(); p += 8) {
      // pxy =(p.x, p.y)
      const __m256 pxy_l = _mm256_load_ps(position_data + p * 2);
      const __m256 pxy_h = _mm256_load_ps(position_data + p * 2 + 8);
      // rxy = pxy - qxy
      const __m256 rxy_l = _mm256_sub_ps(pxy_l, qxy);
      const __m256 rxy_h = _mm256_sub_ps(pxy_h, qxy);
      // rxy2 = rxy * rxy
      const __m256 rxy2_l = _mm256_mul_ps(rxy_l, rxy_l);
      const __m256 rxy2_h = _mm256_mul_ps(rxy_h, rxy_h);
      // r2 = (rx * rx + ry * ry)
      // r2 = ([7] [6] [3] [2] [5] [4] [1] [0])
      __m256 r2 = _mm256_hadd_ps(rxy2_l, rxy2_h);
      if (ApplyThreshold) {
        // r2 = max(r2, threshold)
        r2 = _mm256_max_ps(r2, threshold);
      }
      // c2 = 1. / r2
      const __m256 c2 = _mm256_rcp_ps(r2);
      // d2 = c2 * R2 = R2 / r2
      const __m256 d2 = _mm256_mul_ps(c2, R2);
      // d6 = d2 * d2 * d2
      const __m256 d6 = _mm256_mul_ps(d2, _mm256_mul_ps(d2, d2));
      // d12 = d6 * d6
      const __m256 d12 = _mm256_mul_ps(d6, d6);
      // k = (d12 - d6) * sigma * c2
      __m256 k =
          _mm256_mul_ps(sigma, _mm256_mul_ps(c2, _mm256_sub_ps(d12, d6)));

      k = _mm256_max_ps(k, zero);

      // Contact dashpot: k += -kDashpot * overlap * dot(vrel, rxy) / r2,
      // which acts along rxy and so damps only the normal relative motion.
      // vxy = (p.v.x, p.v.y)
      const __m256 vxy_l = _mm256_load_ps(velocity_data + p * 2);
      const __m256 vxy_h = _mm256_load_ps(velocity_data + p * 2 + 8);
      // vrel = vxy - qvxy
      const __m256 vrel_l = _mm256_sub_ps(vxy_l, qvxy);
      const __m256 vrel_h = _mm256_sub_ps(vxy_h, qvxy);
      // rv = dot(vrel, rxy), in the same lane order as r2
      const __m256 rv = _mm256_hadd_ps(
          _mm256_mul_ps(vrel_l, rxy_l), _mm256_mul_ps(vrel_h, rxy_h));
      // overlap = max(0, 1 - r2 / R2), zero once the contact breaks
      const __m256 overlap =
          _mm256_max_ps(_mm256_sub_ps(one, _mm256_mul_ps(r2, R2inv)), zero);
      // kd = -kDashpot * overlap * rv * c2
      const __m256 kd =
          _mm256_mul_ps(_mm256_mul_ps(dashpot, overlap), _mm256_mul_ps(rv, c2));
      k = _mm256_add_ps(k, kd);

      // lo = k([3] [3] [2] [2] [1] [1] [0] [0])
      const __m256 kxy_l = _mm256_unpacklo_ps(k, k);
      // hi = k([7] [7] [6] [6] [5] [5] [4] [4])
      const __m256 kxy_h = _mm256_unpackhi_ps(k, k);

      // load force to fxy
      __m256 fxy_l = _mm256_load_ps(force_data + p * 2);
      __m256 fxy_h = _mm256_load_ps(force_data + p * 2 + 8);
      // fxy += rxy * kxy
      fxy_l = _mm256_add_ps(fxy_l, _mm256_mul_ps(kxy_l, rxy_l));
      fxy_h = _mm256_add_ps(fxy_h, _mm256_mul_ps(kxy_h, rxy_h));
      // store force
      _mm256_store_ps(force_data + p * 2, fxy_l);
      _mm256_store_ps(force_data + p * 2 + 8, fxy_h);
    }
  }
}
#else
#define CALC_FORCE CalcForceSerial
#endif

void Particles::UpdatePortalBlocks(Portal& portal) {
  portal.blocks.clear();
  for (size_t iblock = 0; iblock < Blocks.GetNumBlocks(); ++iblock) {
    if (portal.IsClose(Blocks.GetCenter(iblock), Blocks.GetCircumRadius())) {
      portal.blocks.push_back(iblock);
    }
  }
}

void Particles::UpdateEnvObj() {
  block_envobj_.clear();
  block_envobj_.resize(Blocks.GetNumBlocks());
  for (size_t iblock = 0; iblock < Blocks.GetNumBlocks(); ++iblock) {
    for (size_t k = 0; k < ENVOBJ.size(); ++k) {
      if (ENVOBJ[k]->IsClose(
              Blocks.GetCenter(iblock), Blocks.GetCircumRadius())) {
        block_envobj_[iblock].push_back(k);
      }
    }
  }

  for (auto& pair : portals_) {
    for (int d = 0; d <= 1; ++d) {
      UpdatePortalBlocks(pair[d]);
    }
  }
}

Vect F12_bond(Vect p1, Vect p2) {
  const Scal sigma = kSigmaBond;
  const Scal R = 2. * kRadius;
  const Vect dp = p1 - p2;
  const Scal r = dp.length();
  const Scal kLimit = 5;
  const Scal k = std::max(1. - r / R, 1. - kLimit);

  return dp * (sigma * k);
}

Vect F12_portal_edge(Vect p1, Vect p2) {
  const Scal sigma = kSigmaPortalEdge;
  const Scal R = kRadiusPortalEdge;
  const Vect dp = p1 - p2;
  const Scal r2 = dp.dot(dp);
  const Scal r2inv = 1. / r2;
  const Scal d2 = r2inv * (R * R);
  const Scal d6 = d2 * d2 * d2;
  const Scal d12 = d6 * d6;
  const Scal k = std::max<Scal>(0., sigma * (d12 - d6) * r2inv);
  return dp * k;
}

Vect F12_pick(Vect p1, Vect p2) {
  const Scal sigma = kSigmaPick;
  const Scal R = 0.1 * kRadius;
  const Vect dp = p1 - p2;
  const Scal r = dp.length();
  const Scal k = (R - r) / R;

  return dp * (sigma * k);
}

void Particles::RHS_bonds() {
  const auto& bbi = Blocks.GetBlockById();
  auto& data = Blocks.GetData();
  no_rendering_.clear();
  for (auto bond : bonds_) {
    const auto& p = bbi[bond.first];
    const auto& q = bbi[bond.second];
    assert(p.first != blocks::kBlockNone);
    assert(q.first != blocks::kBlockNone);

    const Vect p_pos = data.position[p.first][p.second];
    const Vect q_pos = data.position[q.first][q.second];

    Vect res;

    // TODO: revise
    const Scal kCutoff = 8. * kRadius;
    // If particles are close, just apply the bond force
    if (p_pos.dist(q_pos) < kCutoff) {
      res = F12_bond(p_pos, q_pos);
    } else {
      bool found = false;
      // Otherwise, apply the force through the portals (if any)
      for (auto& pair : portals_) {
        for (size_t d = 0; d <= 1; ++d) {
          const auto& portal = pair[d];
          const auto& other = pair[1 - d];

          const Vect a = portal.begin;
          const Vect b = portal.end;
          const Vect other_a = other.begin;
          const Vect other_b = other.end;

          const Vect r = b - a;
          const Vect n = Vect(-r.y, r.x).GetNormalized();
          const Vect other_r = other_b - other_a;
          const Vect other_n = Vect(-other_r.y, other_r.x).GetNormalized();

          const Scal p_lambda = r.dot(p_pos - a) / r.dot(r);
          const Scal p_offset = (p_pos - a).dot(n);
          const Scal q_lambda_other =
              other_r.dot(q_pos - other_a) / other_r.dot(other_r);
          const Scal q_offset_other = (q_pos - other_a).dot(other_n);

          if (p_lambda > 0. && p_lambda < 1. && q_lambda_other > 0. &&
              q_lambda_other < 1. && p_offset * q_offset_other < 0. &&
              std::abs(p_offset - q_offset_other) <
                  2. * kPortalThickness + kCutoff) {
            const Scal sign = (p_offset > 0. ? 1. : -1.);
            const Vect proj =
                a + r * q_lambda_other +
                n * (q_offset_other + 2. * sign * kPortalThickness);
            res = F12_bond(p_pos, proj);
            found = true;
            no_rendering_.emplace(bond.first, bond.second);
            no_rendering_.emplace(bond.second, bond.first);
          }
        }
      }
      if (!found) {
        res = F12_bond(p_pos, q_pos);
      }
    }
    data.force[p.first][p.second] += res;
  }
}

void Particles::ApplyFrozen() {
  const auto& bbi = Blocks.GetBlockById();
  auto& data = Blocks.GetData();
  for (int id : frozen_) {
    const auto iblock = bbi[id].first;
    const auto p = bbi[id].second;
    data.velocity[iblock][p] *= 0.;
    data.force[iblock][p] *= 0.;
  }
}

void Particles::calc_forces(size_t iblock) {
  auto& data = Blocks.GetData();
  for (size_t p = 0; p < data.position[iblock].size(); ++p) {
    auto& f = data.force[iblock][p];
    auto& x = data.position[iblock][p];
    auto& v = data.velocity[iblock][p];

    f = Vect(0);

    // Gravity.
    if (gravity_enable_) {
      f = gravity_ * kMass;
    }

    // Attraction or repulsion force.
    if (force_enabled) {
      const Vect r = x - force_center;
      if (r.length() > kRadius) {
        if (force_attractive_) {
          f += r * (-kPointForceAttractive / std::pow(r.length(), 3));
        } else {
          f += r * (kPointForce / std::pow(r.length(), 4));
        }
      }
    }

    // Portal edges.
    for (auto& pair : portals_) {
      for (auto& portal : pair) {
        f += F12_portal_edge(x, portal.begin);
        f += F12_portal_edge(x, portal.end);
      }
    }

    // Force from Pick tool.
    if (pick_enabled_ && data.id[iblock][p] == pick_particle_id_) {
      f += F12_pick(x, pick_pointer_);
    }

    // dissipation
    // f -= v * (0.1 * p1.m);
    f -= v * (kDissipation * kMass);
  }

  // pairwise interactions
  for (int offset : Blocks.GetNeighborOffsets()) {
    const size_t j = iblock + offset;
    // TODO: revise outside block condition
    if (j >= Blocks.GetNumBlocks()) {
      continue;
    }

    if (iblock != j) { // no check for self-force needed
      CALC_FORCE<false>(
          data.force[iblock], data.position[iblock], data.position[j],
          data.velocity[iblock], data.velocity[j]);
    }

    if (iblock == j) { // apply threshold to distance to avoid self-force
      CALC_FORCE<true>(
          data.force[iblock], data.position[iblock], data.position[j],
          data.velocity[iblock], data.velocity[j]);
    }
  }

  // environment objects
  for (size_t k : block_envobj_[iblock]) {
    for (size_t p = 0; p < data.position[iblock].size(); ++p) {
      auto& obj = ENVOBJ[k];
      data.force[iblock][p] +=
          obj->F(data.position[iblock][p], data.velocity[iblock][p]);
    }
  }
}
