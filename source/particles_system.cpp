#include "particles_system.hpp"

particles_system::particles_system() : 
    domain(rect_vect(vect(-1.,-1.),vect(1.,1.))),
    Blocks(domain, vect(4*kRadius, 4*kRadius))
{
  force_enabled = false;
  force_center = vect(0., 0.);

  // place particles in the domain
  const Scal r = kRadius;
  const int row=1. / r;
  const int N=row * row;

  std::vector<particle> P;
  for(int i=0; i<N; ++i)
  {
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
  const Scal gravity = 10.;
  g = vect(0.0, -1.0) * gravity;

  SetDomain(domain);
  resize_queue_ = domain;
  ResetEnvObjFrame(domain);
}
particles_system::~particles_system() {}
const std::vector<particle>& particles_system::GetParticles() {
  return particle_buffer_;
}
void particles_system::SetParticleBuffer() {
  std::vector<particle> res;
  for (size_t i = 0; i < Blocks.GetNumBlocks(); ++i) {
    const auto& data = Blocks.GetData();
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

  while (t < time_target) {
    #pragma omp for
    for (size_t i = 0; i < Blocks.GetNumBlocks(); ++i) {
      RHS(i);
    }

    #pragma omp for 
    for (size_t i = 0; i < Blocks.GetNumBlocks(); ++i) {
      auto& data = Blocks.GetData();
      for (size_t p = 0; p < data.position[i].size(); ++p) {
        data.position_tmp[i][p] = data.position[i][p]; 
        data.velocity_tmp[i][p] = data.velocity[i][p]; 
        data.position[i][p] += data.velocity[i][p] * dt * 0.5;
        data.velocity[i][p] += data.force[i][p] * dt * 0.5;
      }
    }

    #pragma omp single
    t += 0.5 * dt;

    #pragma omp for 
    for (size_t i = 0; i < Blocks.GetNumBlocks(); ++i) {
      RHS(i);
    }

    #pragma omp for
    for (size_t i = 0; i < Blocks.GetNumBlocks(); ++i) {
      auto& data = Blocks.GetData();
      for (size_t p = 0; p < data.position[i].size(); ++p) {
        data.velocity[i][p] = 
            data.velocity_tmp[i][p] + data.force[i][p] * dt;
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
      auto limited = [](Scal& current, const Scal target) {
        const Scal limit = 0.02;
        current = std::min(current + limit,
                           std::max(current - limit, target));
      };

      rect_vect new_domain = domain;
      limited(new_domain.A.x, resize_queue_.A.x);
      limited(new_domain.A.y, resize_queue_.A.y);
      limited(new_domain.B.x, resize_queue_.B.x);
      limited(new_domain.B.y, resize_queue_.B.y);

      if (int(t / dt) % int(0.01 / dt) == 0) {
        if (new_domain != domain) {
          SetDomain(new_domain);
          ResetEnvObjFrame(new_domain);
        }
      }
      SetParticleBuffer();
      t += 0.5 * dt;
      Blocks.SortParticles();
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

vect F12(vect p1, vect p2)
{
  const Scal sigma = kSigma;
  const Scal R = 2. * kRadius;
  const vect r = p1 - p2;
  const Scal ar2 = r.dot(r);
  const Scal ad2 = 1. / ar2;
  const Scal d2 = ad2 * (R * R);
  const Scal d6 = d2 * d2 * d2;
  const Scal d12 = d6 * d6;
  return r * (sigma * (d12 - d6) * ad2);
}

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
      const __m256 k = _mm256_mul_ps(
          sigma, _mm256_mul_ps(c2, _mm256_sub_ps(d12, d6)));

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

#define CALC_FORCE CalcForceAvx
//#define CALC_FORCE CalcForceSerial
//#define CALC_FORCE CalcForceSerialPadded


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
    //f=g*part.m;
    f = g*kMass;

    // point force
    if (force_enabled) {
      const vect r = x - force_center; 
      f += r * (kPointForce / std::pow(r.length(), 3));
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

  for (size_t p = 0; p < data.position[i].size(); ++p) {
    auto& f = data.force[i][p];
    f *= 1. / kMass;
  }

  // environment objects
  for (size_t k : block_envobj_[i]) {
    for (size_t p = 0; p < data.position[i].size(); ++p) {
      auto& obj=ENVOBJ[k];
      data.force[i][p] += obj->F(
          data.position[i][p], data.velocity[i][p], kRadius,kSigma);
    }
  }
}

vect F12(vect p1, vect /*v1*/, vect p2, vect /*v2*/, 
         Scal /*sigma*/, Scal /*R*/) {
  return F12(p1, p2);
}
