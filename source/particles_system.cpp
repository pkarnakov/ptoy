#include "particles_system.hpp"
#include <thread>
#include <chrono>

particles_system::particles_system() : 
    domain(rect_vect(vect(-1.,-1.),vect(1.,1.))),
    Blocks(domain, vect(kBlockSize, kBlockSize)),
    blocks_buffer_(Blocks)
{
  force_enabled = false;
  force_center = vect(0., 0.);

  // place particles in the domain
  const Scal r = kRadius;
  /*
  const int row = 1. / r;
  const int col = 0.5 / r;
  const int N = row * col;

  std::vector<particle> P;
  for(int i=0; i<N; ++i) {
    const vect p((i%row*2.0+1.0)*r-1.0, (i/row*2.0+1.0)*r-1.0);
    const vect v(0., 0.);
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

  const Scal coeff = 1;
  std::vector<particle> P;
  rect_vect box(vect(-0.2, -1.), vect(0.2, -1. + 30. * kRadius));
  for (Scal x = box.A.x + r; x + r * 0.9 < box.B.x; x += 2. * r) {
    for (Scal y = box.A.y + r; y + r * 0.9 < box.B.y; y += 2. * r) {
      P.push_back(particle(
          vect(x*coeff, y), vect(0., 0.), 
          kMass, r, kSigma, 0x1, rgb(1., 0., 0.)));
    }
  }

  const Scal dx = kPortalThickness;
  if (1) {
    PortalStart(vect(box.A.x * coeff-dx, box.A.y));
    PortalStop(vect(box.A.x * coeff -dx, box.B.y + 0.2));
    PortalStart(vect(box.B.x * coeff +dx, box.A.y));
    PortalStop(vect(box.B.x * coeff +dx, box.B.y + 0.2));
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
  g = vect(0.0, -1.0) * kGravity;

  SetDomain(domain);
  SetParticleBuffer();
  resize_queue_ = domain;
  ResetEnvObjFrame(domain);
}
particles_system::~particles_system() {}
void particles_system::SetParticleBuffer() {
  blocks_buffer_ = Blocks;
  std::vector<particle>res;
  for (size_t i = 0; i < blocks_buffer_.GetNumBlocks(); ++i) {
    const auto& data = blocks_buffer_.GetData();
    for (size_t p = 0; p < data.position[i].size(); ++p) {
      res.push_back(particle(
          data.position[i][p], data.velocity[i][p],
          0.01, kRadius, kSigma, 0x1, rgb(1., 0., 0.)));
    }
  }
  particle_buffer_ = res;
}
void particles_system::AddEnvObj(env_object* env) {
  ENVOBJ.push_back(std::unique_ptr<env_object>(env));
}
void particles_system::status(std::ostream& out)
{
  out << "status N/A";
  //out<<"Particles system"<<std::endl<<"Particles number = "<<P.size()<<std::endl;
}
void particles_system::step(Scal time_target, const std::atomic<bool>& quit)
{
  #pragma omp parallel
  {

  while (t < time_target && !quit.load()) {
    #pragma omp for
    for (size_t i = 0; i < Blocks.GetNumBlocks(); ++i) {
      RHS(i);
    }

    #pragma omp single 
    {
      RHS_bonds();
      DetectPortals(dt);
      ApplyPortalsForces();
      DetectPortals(dt);
      ApplyFrozen();
    }

    #pragma omp for 
    for (size_t i = 0; i < Blocks.GetNumBlocks(); ++i) {
      auto& data = Blocks.GetData();
      for (size_t p = 0; p < data.position[i].size(); ++p) {
        data.velocity_tmp[i][p] = data.velocity[i][p]; 
        data.position_tmp[i][p] = data.position[i][p]; 
        data.velocity[i][p] += data.force[i][p] * (dt * 0.5 / kMass);
        data.position[i][p] += data.velocity[i][p] * dt * 0.5;
      }
    }

    #pragma omp single
    t += 0.5 * dt;
    
    #pragma omp for 
    for (size_t i = 0; i < Blocks.GetNumBlocks(); ++i) {
      RHS(i);
    }

    #pragma omp single 
    {
      RHS_bonds();
      DetectPortals(dt * 0.5);
      ApplyPortalsForces();
      DetectPortals(dt * 0.5);
      ApplyFrozen();
    }

    #pragma omp for
    for (size_t i = 0; i < Blocks.GetNumBlocks(); ++i) {
      auto& data = Blocks.GetData();
      for (size_t p = 0; p < data.position[i].size(); ++p) {
        data.velocity[i][p] = 
            data.velocity_tmp[i][p] + data.force[i][p] * (dt / kMass);
        const Scal limit = 10.;
        if (data.velocity[i][p].length() > limit) {
          data.velocity[i][p] *= limit / data.velocity[i][p].length();
        }
        data.position[i][p] = 
            data.position_tmp[i][p] + data.velocity[i][p] * dt;
      }
    }

    #pragma omp single
    {
      DetectPortals(dt);
      // Resize the frame if needed (with a limited speed)
      auto limit = [](Scal& current, const Scal target) {
        const Scal limit = 0.02;
        current = std::min(current + limit,
                           std::max(current - limit, target));
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

      Blocks.SortParticles();
      ApplyPortals();
      Blocks.SortParticles();

      // Pass the data to renderer if ready
      if (renderer_ready_for_next_) {
        std::lock_guard<std::mutex> lg(m_buffer_);
        SetParticleBuffer();
        renderer_ready_for_next_ = false;
        //std::this_thread::sleep_for(std::chrono::milliseconds(25));
      }

      // Advance in time (another half)
      t += 0.5 * dt;
    }
  }
  }
}

void particles_system::SetForce(vect center, bool enabled) {
  force_center = center;
  force_enabled = enabled;
}
void particles_system::SetForce(vect center) {
  force_center = center;
}
void particles_system::SetForce(bool enabled) {
  force_enabled = enabled;
}

void particles_system::PickStart(vect point) {
  size_t min_block = blocks::kBlockNone;
  size_t min_particle;
  Scal min_dist;

  auto& data = blocks_buffer_.GetData();
  for (size_t i = 0; i < blocks_buffer_.GetNumBlocks(); ++i) {
    for (size_t p = 0; p < data.position[i].size(); ++p) {
      if (min_block == blocks::kBlockNone  
          || data.position[i][p].dist(point) < min_dist) {
        min_block = i;
        min_particle = p;
        min_dist = data.position[i][p].dist(point);
      }
    }
  }

  pick_particle_id_ = data.id[min_block][min_particle];
  pick_pointer_ = point;
  pick_enabled_ = true;
}

void particles_system::PickMove(vect point) {
  if (!pick_enabled_) {
    return;
  }
  pick_pointer_ = point;
}

void particles_system::PickStop(vect) {
  if (!pick_enabled_) {
    return;
  }
  pick_enabled_ = false;
}

void particles_system::MoveToPortal(
    vect& position, vect& velocity,
    const Portal& src, const Portal& dest) {
  const vect src_a = src.begin;
  const vect src_b = src.end; 
  const vect src_r = src_b - src_a; 
  const vect src_n = vect(-src_r.y, src_r.x).GetNormalized();

  const vect dest_a = dest.begin;
  const vect dest_b = dest.end; 
  const vect dest_r = dest_b - dest_a; 
  const vect dest_n = vect(-dest_r.y, dest_r.x).GetNormalized();

  auto Move = [src_a, src_r, src_n, dest_a, dest_r, dest_n](
      vect& v) {
    const Scal lambda = src_r.dot(v - src_a) / src_r.dot(src_r);
    const Scal offset = src_n.dot(v - src_a);
    
    Scal nof;

    if (offset < 0.) {
      nof = offset + 2. * kPortalThickness;
    } else {
      nof = offset - 2. * kPortalThickness;
    }

    v = dest_a + dest_r * lambda + dest_n * nof;
  };

  Move(position);
  //Move(velocity); 
}

void particles_system::DetectPortals(const Scal local_dt) {
  for (auto& pair : portals_) {
    for (int d = 0; d <= 1; ++d) {
      const auto& portal = pair[d];

      const vect a = portal.begin;
      const vect b = portal.end;
      const vect r = b - a; 
      const vect n = vect(-r.y, r.x).GetNormalized();

      auto& data = Blocks.GetData();
      for (size_t i : portal.blocks) {
        for (size_t p = 0; p < data.position[i].size(); ++p) {
          auto id = static_cast<size_t>(data.id[i][p]);
          if (particle_to_move_.size() <= id) {
            particle_to_move_.resize(id + 1);
          }
          const vect curr = data.position[i][p];
          const vect prev = curr - data.velocity[i][p] * local_dt;
          const Scal lambda_curr = (curr - a).dot(r) / r.dot(r);
          const Scal offset_curr = (curr - a).dot(n);
          const Scal offset_prev = (prev - a).dot(n);
          if (lambda_curr > 0. && lambda_curr < 1.
              && std::abs(offset_curr) <= kPortalThickness 
              && std::abs(offset_prev) > kPortalThickness) {
            //particle_to_move_[id] = 1 - particle_to_move_[id];
            particle_to_move_[id] = 1;
            //std::cout << "tomove(" << id << "): " 
            //    << particle_to_move_[id] << std::endl;
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
      auto& other = pair[1-d];
      vect a = portal.begin;
      vect b = portal.end;
      vect other_a = other.begin;
      vect other_b = other.end;

      const vect r = b - a; 
      const vect n = vect(-r.y, r.x).GetNormalized();
      const vect other_r = other_b - other_a;
      const vect other_n = vect(-other_r.y, other_r.x).GetNormalized();
      auto& data = Blocks.GetData();
      for (size_t i : portal.blocks) {
        for (size_t p = 0; p < data.position[i].size(); ++p) {
          const vect curr = data.position[i][p];
          const auto id_curr = data.id[i][p];
          const Scal lambda_curr = r.dot(curr - a) / r.dot(r);
          const Scal offset_curr = (curr - a).dot(n);
          const Scal moffset_curr = 
              offset_curr * (particle_to_move_[id_curr] - 0.5);

          // Check particle forces
          if (lambda_curr > 0. && lambda_curr < 1.) {
            for (size_t j : portal.blocks) {
              for (size_t q = 0; q < data.position[j].size(); ++q) {
                const vect neighbor = data.position[j][q];
                const auto id_neighbor = data.id[j][q];
                const Scal lambda_neighbor = r.dot(neighbor - a) / r.dot(r);
                const Scal offset_neighbor = (neighbor - a).dot(n);
                const Scal moffset_neighbor =
                    offset_neighbor * (particle_to_move_[id_neighbor] - 0.5);
                if (lambda_neighbor > 0. && lambda_neighbor < 1.
                    && moffset_neighbor * moffset_curr < 0.) {
                  data.force[i][p] -= F12(curr, neighbor);
                }
              }
            }
          }

          if (lambda_curr > 0. && lambda_curr < 1.) {
            for (size_t j : other.blocks) {
              for (size_t q = 0; q < data.position[j].size(); ++q) {
                const vect neighbor = data.position[j][q];
                const auto id_neighbor = data.id[j][q];
                const Scal other_lambda_neighbor = 
                    other_r.dot(neighbor - other_a) / other_r.dot(other_r);
                const vect other_q_neighbor = 
                    other_a + other_r * other_lambda_neighbor;
                const Scal other_offset_neighbor = 
                    (neighbor - other_a).dot(other_n);
                if (other_lambda_neighbor > 0. && other_lambda_neighbor < 1.
                    && other_offset_neighbor * offset_curr < 0.) {
                  Scal nof;
                  if (offset_curr > 0.) {
                    nof = other_offset_neighbor + 2. * kPortalThickness;
                  } else {
                    nof = other_offset_neighbor - 2. * kPortalThickness;
                  }
                  //data.force[i][p] += 
                  //    F12(curr, 
                  //        a + r * other_lambda_neighbor + 
                  //        n * other_offset_neighbor);
                  const vect proj = a + r * other_lambda_neighbor + n * nof;
                  data.force[i][p] += F12(curr, proj);
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
      auto& other = pair[1-d];

      for (size_t i : portal.blocks) {
        for (size_t p = 0; p < data.position[i].size(); ++p) {
          const auto id = data.id[i][p];
          if (particle_to_move_.size() <= id) {
            particle_to_move_.resize(id + 1);
          }
          if (particle_to_move_[id]) {
            MoveToPortal(data.position[i][p], data.velocity[i][p],
                portal, other);
            particle_to_move_[id] = 0;
            //std::cout << "moved: " << id << std::endl;
          } 
        }
      }
    }
  }
}

void particles_system::PortalStart(vect point) {
  portal_enabled_ = true;
  portal_begin_ = point;
  portal_current_ = point;
  portal_mouse_moving_ = true;
}

void particles_system::PortalMove(vect point) {
  if (!portal_enabled_) {
    return;
  }
  portal_current_ = point;
}

void particles_system::PortalStop(vect point) {
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
    pair[1].end = pair[1].begin +
        (pair[1].end - pair[1].begin).GetNormalized() *
        pair[0].begin.dist(pair[0].end);
    portals_.push_back(pair);
    portal_stage_ = 0;

    UpdatePortalBlocks(portals_.back()[0]);
    UpdatePortalBlocks(portals_.back()[1]);
  }
}

void particles_system::BondsStart(vect point) {
  bonds_enabled_ = true;
  int id = kParticleIdNone;

  for (size_t i = 0; i < blocks_buffer_.GetNumBlocks(); ++i) {
    auto& data = blocks_buffer_.GetData();
    for (size_t p = 0; p < data.position[i].size(); ++p) {
      if (data.position[i][p].dist(point) < kRadius) {
        id = data.id[i][p];
        std::cout << "Found particle id=" << id << std::endl;
      }
    }
  }

  bonds_prev_particle_id_ = id;
}

void particles_system::BondsMove(vect point) {
  if (!bonds_enabled_) {
    return;
  }
  int id = kParticleIdNone;

  for (size_t i = 0; i < blocks_buffer_.GetNumBlocks(); ++i) {
    auto& data = blocks_buffer_.GetData();
    for (size_t p = 0; p < data.position[i].size(); ++p) {
      if (data.position[i][p].dist(point) < kRadius &&
          data.id[i][p] != bonds_prev_particle_id_) {
        id = data.id[i][p];
        std::cout << "Found next particle id=" << id << std::endl;
      }
    }
  }

  if (bonds_prev_particle_id_ != kParticleIdNone && 
      id != kParticleIdNone) {
    assert(bonds_prev_particle_id_ != id);
    bonds_.emplace_back(bonds_prev_particle_id_, id);
  }
  if (id != kParticleIdNone) {
    bonds_prev_particle_id_ = id;
  }
}

void particles_system::BondsStop(vect) {
  if (!bonds_enabled_) {
    return;
  }
  bonds_enabled_ = false;
  for (auto bond : bonds_) {
    std::cout << "(" << bond.first << ", " << bond.second << ") ";
  }
  std::cout << std::endl;
}

void particles_system::FreezeStart(vect point) {
  freeze_last_id_ = -1;
  freeze_enabled_ = true;
  FreezeMove(point);
}

void particles_system::FreezeMove(vect point) {
  if (!freeze_enabled_) {
    return;
  }

  for (size_t i = 0; i < blocks_buffer_.GetNumBlocks(); ++i) {
    auto& data = blocks_buffer_.GetData();
    for (size_t p = 0; p < data.position[i].size(); ++p) {
      if (data.position[i][p].dist(point) < kRadius) {
        const int id = data.id[i][p];
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

void particles_system::FreezeStop(vect) {
  if (!freeze_enabled_) {
    return;
  }
  freeze_enabled_ = false;
}

vect F12(vect p1, vect p2, Scal R, Scal sigma) {
  const vect dp = p1 - p2;
  const Scal r2 = dp.dot(dp);
  const Scal r2inv = 1. / r2;
  const Scal d2 = r2inv * (R * R);
  const Scal d6 = d2 * d2 * d2;
  const Scal d12 = d6 * d6;
  return dp * std::max<Scal>(0., sigma * (d12 - d6) * r2inv);
}

vect F12wall(vect p1, vect p2) {
  const Scal sigma = kSigmaWall;
  const Scal R = kRadius;
  const vect dp = p1 - p2;
  const Scal r2 = dp.dot(dp);
  const Scal r2inv = 1. / r2;
  const Scal d2 = r2inv * (R * R);
  const Scal d6 = d2 * d2 * d2;
  const Scal d12 = d6 * d6;
  return dp * std::max<Scal>(0., sigma * (d12 - d6) * r2inv);
}

vect F12(vect p1, vect p2) {
  const Scal threshold = std::pow(kRadius, 2) * 1e-3;
  const Scal sigma = kSigma;
  const Scal R = 2. * kRadius;
  const vect dp = p1 - p2;
  const Scal r2 = std::max(threshold, dp.dot(dp));
  const Scal r2inv = 1. / r2;
  const Scal d2 = r2inv * (R * R);
  const Scal d6 = d2 * d2 * d2;
  const Scal d12 = d6 * d6;
  return dp * std::max<Scal>(0., sigma * (d12 - d6) * r2inv);
}

//#define CALC_FORCE CalcForceAvx
#define CALC_FORCE CalcForceSerial
//#define CALC_FORCE CalcForceSerialPadded


template <bool ApplyThreshold=true>
void CalcForceSerial(ArrayVect& force,
               ArrayVect& position,
               ArrayVect& position_other) {

  for (size_t q = 0; q < position_other.size(); ++q) {
    for (size_t p = 0; p < position.size(); ++p) {
if (&position[p] != &position_other[q])
      force[p] += F12(position[p], position_other[q]);
    }
  }
}

template <bool ApplyThreshold=true>
void CalcForceSerialPadded(ArrayVect& force,
               ArrayVect& position,
               ArrayVect& position_other) {

  for (size_t q = 0; q < position_other.size(); ++q) {
    for (size_t p = 0; p < position.size(); p+=8) {
      force[p] += F12(position[p], position_other[q]);
      force[p+1] += F12(position[p+1], position_other[q]);
      force[p+2] += F12(position[p+2], position_other[q]);
      force[p+3] += F12(position[p+3], position_other[q]);
      force[p+4] += F12(position[p+4], position_other[q]);
      force[p+5] += F12(position[p+5], position_other[q]);
      force[p+6] += F12(position[p+6], position_other[q]);
      force[p+7] += F12(position[p+7], position_other[q]);
    }
  }
}

void print(const __m256& mm) {
  float a[8];
  _mm256_store_ps(a, mm);
  for (size_t i = 0; i < 8; ++i) {
    std::cout << a[i] << " ";
  }
  std::cout << std::endl;
}

#ifdef asdf
template <bool ApplyThreshold=true>
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
    for (size_t p = 0; p < position.size(); p+=8) {
      // pxy =(p.x, p.y)
      const __m256 pxy_l = _mm256_load_ps((float*)&position[p]);
      const __m256 pxy_h = _mm256_load_ps((float*)&position[p+4]);
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
      __m256 k = _mm256_mul_ps(
          sigma, _mm256_mul_ps(c2, _mm256_sub_ps(d12, d6)));

      k = _mm256_max_ps(k, zero);

      // lo = k([3] [3] [2] [2] [1] [1] [0] [0])
      const __m256 kxy_l = _mm256_unpacklo_ps(k, k);
      // hi = k([7] [7] [6] [6] [5] [5] [4] [4]) 
      const __m256 kxy_h = _mm256_unpackhi_ps(k, k);

      // load force to fxy
      __m256 fxy_l = _mm256_load_ps((float*)&force[p]);
      __m256 fxy_h = _mm256_load_ps((float*)&force[p+4]);
      // fxy += rxy * kxy
      fxy_l =  _mm256_add_ps(fxy_l, _mm256_mul_ps(kxy_l, rxy_l));
      fxy_h =  _mm256_add_ps(fxy_h, _mm256_mul_ps(kxy_h, rxy_h));
      // store force
      _mm256_store_ps((float*)&force[p], fxy_l);
      _mm256_store_ps((float*)&force[p+4], fxy_h);
    }
  }
}
#endif

#ifdef FALSE
void TestUni() {
  std::cout.precision(8);
  std::cout << std::scientific;
  const vect p0(0., 0.);
  const vect q0(1., 0.);
  const vect zero(0., 0.);

  const size_t N = 8;
  ArrayVect position(N, p0);
  ArrayVect position_other(N, q0);

  for (size_t i = 0; i < N; ++i) {
    const Scal k = 0.1;
    position[i] += vect(0., k * i);
    position_other[i] += vect(0., k * i);
  }

  vect offset(1., 0.);
  for (size_t i = 0; i < N; ++i) {
    position[i] += offset;
    position_other[i] += offset;
  }


  ArrayVect force2(N, zero);
  CalcForceSerial(force2, position, position_other);

  ArrayVect force(N, zero);
  CalcForceAvx(force, position, position_other);

  for (size_t i = 0; i < N; ++i) {
    std::cout 
        << force2[i] << " " 
        << force[i] << " "
        << (force[i] - force2[i]) * (1. / force2[i].length()) << " "
        << std::endl;
  }
}
#endif

void particles_system::UpdatePortalBlocks(Portal& portal) {
  portal.blocks.clear();
  for (size_t i = 0; i < Blocks.GetNumBlocks(); ++i) {
    if (portal.IsClose(Blocks.GetCenter(i),
                       Blocks.GetCircumRadius())) {
      portal.blocks.push_back(i);
      std::cout << "add block " << i << std::endl;
    }
  }
}


void particles_system::UpdateEnvObj() { 
  block_envobj_.clear();
  block_envobj_.resize(Blocks.GetNumBlocks());
  for (size_t i = 0; i < Blocks.GetNumBlocks(); ++i) {
    for (size_t k = 0; k < ENVOBJ.size(); ++k) {
      if (ENVOBJ[k]->IsClose(Blocks.GetCenter(i),
                             Blocks.GetCircumRadius())) {
        block_envobj_[i].push_back(k);
      }
    }
  }

  for (auto& pair : portals_) {
    for (int d = 0; d <= 1; ++d) {
      UpdatePortalBlocks(pair[d]);
    }
  }
}

vect F12_bond(vect p1, vect p2) {
  const Scal sigma = kSigmaBond;
  const Scal R = 2. * kRadius;
  const vect dp = p1 - p2;
  const Scal r = dp.length();
  const Scal k = (R - r) / R;

  return dp * (sigma * k);
}

void particles_system::RHS_bonds() {
  const auto& bbi = Blocks.GetBlockById();
  auto& data = Blocks.GetData();
  for (auto bond : bonds_) {
    const auto& a = bbi[bond.first];
    const auto& b = bbi[bond.second];
    const auto f = F12_bond(data.position[a.first][a.second],
                        data.position[b.first][b.second]);
    data.force[a.first][a.second] += f;
    data.force[b.first][b.second] -= f;
  }
}

void particles_system::ApplyFrozen() {
  const auto& bbi = Blocks.GetBlockById();
  auto& data = Blocks.GetData();
  for (int id : frozen_) {
    const auto i = bbi[id].first;
    const auto p = bbi[id].second;
    data.velocity[i][p] *= 0.;
    data.force[i][p] *= 0.;
  }
}

void particles_system::RHS(size_t i)
{
  auto& data = Blocks.GetData();
  for (size_t p = 0; p < data.position[i].size(); ++p) {
    auto& f = data.force[i][p];
    auto& x = data.position[i][p];
    auto& v = data.velocity[i][p];

    // zero
    f = vect(0., 0.);

    // gravity
    if (gravity_enable_) {
      f = g*kMass;
    }

    // point force
    if (force_enabled) {
      const vect r = x - force_center; 
      if (r.length() > kRadius) {
        if (force_attractive_) {
          f += r * (-kPointForceAttractive / std::pow(r.length(), 3));
        } else {
          f += r * (kPointForce / std::pow(r.length(), 4));
        }
      }
    }
    
    // pick force
    if (pick_enabled_ && data.id[i][p] == pick_particle_id_) {
      f += F12_bond(x, pick_pointer_);
    }

    // dissipation
    //f -= v * (0.1 * p1.m);
    f -= v * (kDissipation * kMass);
  }

  // pairwise interactions
  for (int offset : Blocks.GetNeighborOffsets()) {
    const size_t j = i + offset;
    // TODO: revise outside block condition
    if (j >= Blocks.GetNumBlocks()) {
      continue;
    }

    if (i != j) { // no check for self-force needed
      CALC_FORCE<false>(data.force[i], data.position[i], data.position[j]);
    }

    if (i == j) { // apply threshold to distance to avoid self-force 
      CALC_FORCE<true>(data.force[i], data.position[i], data.position[j]);
    }
  }

  // environment objects
  for (size_t k : block_envobj_[i]) {
    for (size_t p = 0; p < data.position[i].size(); ++p) {
      auto& obj=ENVOBJ[k];
      data.force[i][p] += obj->F(
          data.position[i][p], data.velocity[i][p]);
    }
  }
}

