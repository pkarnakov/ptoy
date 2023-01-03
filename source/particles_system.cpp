#include "particles_system.hpp"
#include <chrono>
#include <thread>

particles_system::particles_system()
    : domain(rect_vect(Vect(-1., -1.), Vect(1., 1.)))
    , Blocks(domain, Vect(kBlockSize, kBlockSize))
    , blocks_buffer_(Blocks) {
  force_enabled = false;
  force_center = Vect(0., 0.);
  remove_last_portal_ = false;

  // place particles in the domain
  const Scal r = kRadius;
  /*
  const int row = 1. / r;
  const int col = 0.5 / r;
  const int N = row * col;

  std::vector<particle> P;
  for(int i=0; i<N; ++i) {
    const Vect p((i%row*2.0+1.0)*r-1.0, (i/row*2.0+1.0)*r-1.0);
    const Vect v(0., 0.);
    // TODO: adjust sigma so that with r->0 it converges to a solid body
    const Scal sigma = kSigma;
    switch (i % 1) {
      case 0:
        P.push_back(particle(p, v, 0.01, r, sigma, 0x1, rgb(1.,0.,0.)));
        break;
      case 1:
        P.push_back(particle(p, v, 0.02, r, sigma, 0x2, rgb(0.,1.,0.)));
        break;
      case 2:
        P.push_back(particle(p, v, 0.05, r, sigma, 0x4, rgb(0.,0.,1.)));
        break;
    }
  }
  */

  const size_t rows = 40;
  const size_t columns = 40;
  const Scal width = columns * 2. * kRadius;
  const Scal height = rows * std::sqrt(3.) * kRadius;
  std::vector<particle> P;
  rect_vect box(Vect(-0.5 * width, -1.), Vect(0.5 * width, -1. + height));
  for (size_t j = 0; j < rows; ++j) {
    for (size_t i = 0; i < columns; ++i) {
      const Scal x = box.A.x + kRadius * (2. * i + 1. + (j % 2));
      const Scal y = box.A.y + kRadius * (std::sqrt(3.) * j + 1.);
      P.emplace_back(
          Vect(x, y), Vect(0., 0.), kMass, r, kSigma, 0x1, rgb(1., 0., 0.));
    }
  }

  ArrayVect position;
  ArrayVect velocity;
  std::vector<int> id;

  for (auto part : P) {
    position.push_back(part.p);
    velocity.push_back(part.v);
    id.push_back(id.size());
  }

  Blocks.AddParticles(position, velocity, id);

  t = 0.0;
  dt = kTimeStep;
  g = Vect(0.0, -1.0) * kGravity;

  SetDomain(domain);
  SetParticleBuffer();
  resize_queue_ = domain;
  ResetEnvObjFrame(domain);

  const Scal dx = kPortalThickness;
  if (1) {
    PortalStart(Vect(box.A.x - dx, box.A.y));
    PortalStop(Vect(box.A.x - dx, box.B.y + 0.2));
    PortalStart(Vect(box.B.x + dx, box.A.y));
    PortalStop(Vect(box.B.x + dx, box.B.y + 0.2));
  }
}
particles_system::~particles_system() {}
void particles_system::SetParticleBuffer() {
  blocks_buffer_ = Blocks;
  std::vector<particle> res;
  for (size_t iblock = 0; iblock < blocks_buffer_.GetNumBlocks(); ++iblock) {
    const auto& data = blocks_buffer_.GetData();
    for (size_t p = 0; p < data.position[iblock].size(); ++p) {
      res.push_back(particle(
          data.position[iblock][p], data.velocity[iblock][p], 0.01, kRadius,
          kSigma, 0x1, rgb(1., 0., 0.)));
    }
  }
  particle_buffer_ = res;
}
void particles_system::AddEnvObj(env_object* env) {
  ENVOBJ.push_back(std::unique_ptr<env_object>(env));
}
void particles_system::status(std::ostream& out) {
  out << "status N/A";
  // out<<"Particles system"<<std::endl<<"Particles number =
  // "<<P.size()<<std::endl;
}
void particles_system::step(Scal time_target, const std::atomic<bool>& quit) {
#pragma omp parallel
  {
    while (t < time_target && !quit.load()) {
#pragma omp for schedule(dynamic, 8)
      for (size_t iblock = 0; iblock < Blocks.GetNumBlocks(); ++iblock) {
        calc_forces(iblock);
      }

#pragma omp single
      {
        RHS_bonds();
        ApplyPortalsForces();
        ApplyFrozen();
      }

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
      {
        RHS_bonds();
        ApplyPortalsForces();
        ApplyFrozen();
      }

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
              data.position_tmp[iblock][p] + data.velocity[iblock][p] * dt;
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

        rect_vect new_domain = domain;
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
        Blocks.SortParticles();
        CheckBonds();

        if (remove_last_portal_) {
          if (portals_.size()) {
            portals_.pop_back();
          }
          remove_last_portal_ = false;
        }

        // Pass the data to renderer if ready
        if (renderer_ready_for_next_) {
          std::lock_guard<std::mutex> lg(m_buffer_);
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

void particles_system::CheckBonds() {
  const auto bbi = Blocks.GetBlockById();
  for (auto it = bonds_.begin(); it != bonds_.end();) {
    if (bbi[it->first].first == blocks::kBlockNone ||
        bbi[it->second].first == blocks::kBlockNone) {
      it = bonds_.erase(it);
    } else {
      ++it;
    }
  }
}

void particles_system::SetForce(Vect center, bool enabled) {
  force_center = center;
  force_enabled = enabled;
}
void particles_system::SetForce(Vect center) {
  force_center = center;
}
void particles_system::SetForce(bool enabled) {
  force_enabled = enabled;
}

void particles_system::PickStart(Vect point) {
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

void particles_system::PickMove(Vect point) {
  if (!pick_enabled_) {
    return;
  }
  pick_pointer_ = point;
}

void particles_system::PickStop(Vect) {
  if (!pick_enabled_) {
    return;
  }
  pick_enabled_ = false;
}

void particles_system::MoveToPortal(
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

void particles_system::DetectPortals() {
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

void particles_system::ApplyPortalsForces() {
  // Assume that the particles have just been moved
  // with their velocity
  // so that (position - dt * velocity) is the previous position

  for (auto& pair : portals_) {
    for (int d = 0; d <= 1; ++d) {
      auto& portal = pair[d];
      auto& other = pair[1 - d];
      const Vect a = portal.begin;
      const Vect b = portal.end;
      const Vect other_a = other.begin;
      const Vect other_b = other.end;

      const Vect r = b - a;
      const Vect n = Vect(-r.y, r.x).GetNormalized();
      const Vect other_r = other_b - other_a;
      const Vect other_n = Vect(-other_r.y, other_r.x).GetNormalized();
      auto& data = Blocks.GetData();
      for (size_t iblock : portal.blocks) {
        for (size_t p = 0; p < data.position[iblock].size(); ++p) {
          const Vect curr = data.position[iblock][p];
          const Scal lambda_curr = r.dot(curr - a) / r.dot(r);
          const Scal offset_curr = (curr - a).dot(n);

          // Check particle forces
          if (lambda_curr > 0. && lambda_curr < 1.) {
            for (size_t j : portal.blocks) {
              for (size_t q = 0; q < data.position[j].size(); ++q) {
                const Vect neighbor = data.position[j][q];
                const Scal lambda_neighbor = r.dot(neighbor - a) / r.dot(r);
                const Scal offset_neighbor = (neighbor - a).dot(n);
                if (lambda_neighbor > 0. && lambda_neighbor < 1. &&
                    offset_neighbor * offset_curr < 0.) {
                  data.force[iblock][p] -= F12(curr, neighbor);
                }
              }
            }
          }

          if (lambda_curr > 0. && lambda_curr < 1.) {
            for (size_t j : other.blocks) {
              for (size_t q = 0; q < data.position[j].size(); ++q) {
                const Vect neighbor = data.position[j][q];
                const Scal other_lambda_neighbor =
                    other_r.dot(neighbor - other_a) / other_r.dot(other_r);
                const Scal other_offset_neighbor =
                    (neighbor - other_a).dot(other_n);
                if (other_lambda_neighbor > 0. && other_lambda_neighbor < 1. &&
                    other_offset_neighbor * offset_curr < 0.) {
                  const Scal sign = (offset_curr > 0. ? 1. : -1.);
                  const Vect proj = a + r * other_lambda_neighbor +
                                    n * (other_offset_neighbor +
                                         2. * sign * kPortalThickness);
                  data.force[iblock][p] += F12(curr, proj);
                }
              }
            }
          }
        }
      }
    }
  }
}
void particles_system::ApplyPortals() {
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

void particles_system::PortalStart(Vect point) {
  portal_enabled_ = true;
  portal_begin_ = point;
  portal_current_ = point;
  portal_mouse_moving_ = true;
}

void particles_system::PortalMove(Vect point) {
  if (!portal_enabled_) {
    return;
  }
  portal_current_ = point;
}

void particles_system::PortalStop(Vect point) {
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

void particles_system::BondsStart(Vect point) {
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

void particles_system::BondsMove(Vect point) {
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

void particles_system::BondsStop(Vect) {
  if (!bonds_enabled_) {
    return;
  }
  bonds_enabled_ = false;
}

void particles_system::FreezeStart(Vect point) {
  freeze_last_id_ = -1;
  freeze_enabled_ = true;
  FreezeMove(point);
}

void particles_system::FreezeMove(Vect point) {
  if (!freeze_enabled_) {
    return;
  }

  for (size_t iblock = 0; iblock < blocks_buffer_.GetNumBlocks(); ++iblock) {
    auto& data = blocks_buffer_.GetData();
    for (size_t p = 0; p < data.position[iblock].size(); ++p) {
      if (data.position[iblock][p].dist(point) < kRadius) {
        const int id = data.id[iblock][p];
        if (id != freeze_last_id_) {
          if (frozen_.count(id)) {
            frozen_.erase(id);
            std::cout << "Unfreeze particle id=" << id << std::endl;
          } else {
            frozen_.insert(id);
            std::cout << "Freeze particle id=" << id << std::endl;
          }
        }
        freeze_last_id_ = id;
      }
    }
  }
}

void particles_system::FreezeStop(Vect) {
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

template <bool ApplyThreshold = true>
void CalcForceSerial(
    ArrayVect& force, ArrayVect& position, ArrayVect& position_other) {
  for (size_t q = 0; q < position_other.size(); ++q) {
    for (size_t p = 0; p < position.size(); ++p) {
      if (&position[p] != &position_other[q])
        force[p] += F12(position[p], position_other[q]);
    }
  }
}

template <bool ApplyThreshold = true>
void CalcForceSerialPadded(
    ArrayVect& force, ArrayVect& position, ArrayVect& position_other) {
  for (size_t q = 0; q < position_other.size(); ++q) {
    for (size_t p = 0; p < position.size(); p += 8) {
      force[p] += F12(position[p], position_other[q]);
      force[p + 1] += F12(position[p + 1], position_other[q]);
      force[p + 2] += F12(position[p + 2], position_other[q]);
      force[p + 3] += F12(position[p + 3], position_other[q]);
      force[p + 4] += F12(position[p + 4], position_other[q]);
      force[p + 5] += F12(position[p + 5], position_other[q]);
      force[p + 6] += F12(position[p + 6], position_other[q]);
      force[p + 7] += F12(position[p + 7], position_other[q]);
    }
  }
}

#ifdef USE_AVX
#include <x86intrin.h>
#define CALC_FORCE CalcForceAvx
template <bool ApplyThreshold = true>
void CalcForceAvx(
    ArrayVect& force, ArrayVect& position, ArrayVect& position_other) {
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

  for (size_t q = 0; q < position_other.size(); ++q) {
    const __m256 qx = _mm256_broadcast_ss((float*)&position_other[q].x);
    const __m256 qy = _mm256_broadcast_ss((float*)&position_other[q].y);
    // qxy = (q.x, q.y)
    const __m256 qxy = _mm256_blend_ps(qx, qy, 0xAA);
    for (size_t p = 0; p < position.size(); p += 8) {
      // pxy =(p.x, p.y)
      const __m256 pxy_l = _mm256_load_ps((float*)&position[p]);
      const __m256 pxy_h = _mm256_load_ps((float*)&position[p + 4]);
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

      // lo = k([3] [3] [2] [2] [1] [1] [0] [0])
      const __m256 kxy_l = _mm256_unpacklo_ps(k, k);
      // hi = k([7] [7] [6] [6] [5] [5] [4] [4])
      const __m256 kxy_h = _mm256_unpackhi_ps(k, k);

      // load force to fxy
      __m256 fxy_l = _mm256_load_ps((float*)&force[p]);
      __m256 fxy_h = _mm256_load_ps((float*)&force[p + 4]);
      // fxy += rxy * kxy
      fxy_l = _mm256_add_ps(fxy_l, _mm256_mul_ps(kxy_l, rxy_l));
      fxy_h = _mm256_add_ps(fxy_h, _mm256_mul_ps(kxy_h, rxy_h));
      // store force
      _mm256_store_ps((float*)&force[p], fxy_l);
      _mm256_store_ps((float*)&force[p + 4], fxy_h);
    }
  }
}
#else
#define CALC_FORCE CalcForceSerial
#endif

#if 0
void print(const __m256& mm) {
  float a[8];
  _mm256_store_ps(a, mm);
  for (size_t i = 0; i < 8; ++i) {
    std::cout << a[i] << " ";
  }
  std::cout << std::endl;
}
void TestUni() {
  std::cout.precision(8);
  std::cout << std::scientific;
  const Vect p0(0., 0.);
  const Vect q0(1., 0.);
  const Vect zero(0., 0.);

  const size_t N = 8;
  ArrayVect position(N, p0);
  ArrayVect position_other(N, q0);

  for (size_t i = 0; i < N; ++i) {
    const Scal k = 0.1;
    position[i] += Vect(0., k * i);
    position_other[i] += Vect(0., k * i);
  }

  Vect offset(1., 0.);
  for (size_t i = 0; i < N; ++i) {
    position[i] += offset;
    position_other[i] += offset;
  }

  ArrayVect force2(N, zero);
  CalcForceSerial(force2, position, position_other);

  ArrayVect force(N, zero);
  CalcForceAvx(force, position, position_other);

  for (size_t i = 0; i < N; ++i) {
    std::cout << force2[i] << " " << force[i] << " "
              << (force[i] - force2[i]) * (1. / force2[i].length()) << " "
              << std::endl;
  }
}
#endif

void particles_system::UpdatePortalBlocks(Portal& portal) {
  portal.blocks.clear();
  for (size_t iblock = 0; iblock < Blocks.GetNumBlocks(); ++iblock) {
    if (portal.IsClose(Blocks.GetCenter(iblock), Blocks.GetCircumRadius())) {
      portal.blocks.push_back(iblock);
    }
  }
}

void particles_system::UpdateEnvObj() {
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

void particles_system::RHS_bonds() {
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

void particles_system::ApplyFrozen() {
  const auto& bbi = Blocks.GetBlockById();
  auto& data = Blocks.GetData();
  for (int id : frozen_) {
    const auto iblock = bbi[id].first;
    const auto p = bbi[id].second;
    data.velocity[iblock][p] *= 0.;
    data.force[iblock][p] *= 0.;
  }
}

void particles_system::calc_forces(size_t iblock) {
  auto& data = Blocks.GetData();
  for (size_t p = 0; p < data.position[iblock].size(); ++p) {
    auto& f = data.force[iblock][p];
    auto& x = data.position[iblock][p];
    auto& v = data.velocity[iblock][p];

    f = Vect(0);

    // Gravity.
    if (gravity_enable_) {
      f = g * kMass;
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
          data.force[iblock], data.position[iblock], data.position[j]);
    }

    if (iblock == j) { // apply threshold to distance to avoid self-force
      CALC_FORCE<true>(
          data.force[iblock], data.position[iblock], data.position[j]);
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
