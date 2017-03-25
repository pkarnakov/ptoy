#include "particles_system.hpp"

particles_system::particles_system() : 
    domain(rect_vect(vect(-1.,-1.),vect(1.,1.))),
    Blocks(domain, vect(4*kRadius, 4*kRadius))
{
  std::lock_guard<std::mutex> lg(m_step);

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

  std::vector<vect> position;
  std::vector<vect> velocity;
  std::vector<int> id;

  for (auto part : P) {
    position.push_back(part.p);
    velocity.push_back(part.v);
    id.push_back(id.size());
  }

  Blocks.AddParticles(position, velocity, id);

  t=0.0;
  dt=0.0002;
  const Scal gravity = 10.;
  g=vect(0.0, -1.0) * gravity;
}
particles_system::~particles_system() {}
std::vector<particle> particles_system::GetParticles() {
  //std::lock_guard<std::mutex> lg(m_step);
  std::vector<particle> res;
  for (size_t i = 0; i < Blocks.GetNumBlocks(); ++i) {
    const auto& data = Blocks.GetData();
    for (size_t p = 0; p < data.position[i].size(); ++p) {
      res.push_back(particle(
          data.position[i][p], data.velocity[i][p],
          0.01, kRadius, kSigma, 0x1, rgb(1., 0., 0.)));
    }
  }
  return res;
}
void particles_system::AddEnvObj(env_object* env) {
  ENVOBJ.push_back(std::unique_ptr<env_object>(env));
}
void particles_system::status(std::ostream& out)
{
  out << "status N/A";
  //out<<"Particles system"<<std::endl<<"Particles number = "<<P.size()<<std::endl;
}
void particles_system::step(Scal time_target)
{
  std::lock_guard<std::mutex> lg(m_step);

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
    //f = g*kMass;

    // point force
    if (force_enabled) {
      const Scal intensity = 0.1;
      const vect r = x - force_center; 
      f += r * (intensity / std::pow(r.length(), 3));
    }

    // dissipation
    //f -= v * (0.1 * p1.m);
    f -= v * (0.1 * kMass);
  }

  // pairwise interactions
  for (int offset : Blocks.GetNeighborOffsets()) {
    const size_t j = i + offset;
    // TODO: revise outside block condition
    if (j >= Blocks.GetNumBlocks()) {
      continue;
    }

    if (i != j) // no check for self-force needed
    for (size_t p = 0; p < data.position[i].size(); ++p) {
      for (size_t q = 0; q < data.position[j].size(); ++q) {
        //if(&p1!=&p2 && (p1.layers_mask & p2.layers_mask))
        data.force[i][p] += F12(
            data.position[i][p], data.velocity[i][p],
            data.position[j][q], data.velocity[j][q],
            kSigma, kRadius * 2.);
      }
    }

    if (i == j) // need to exclude self-force
    for (size_t p = 0; p < data.position[i].size(); ++p) {
      for (size_t q = 0; q < data.position[j].size(); ++q) {
        //if(&p1!=&p2 && (p1.layers_mask & p2.layers_mask))
        if (p != q) {
          data.force[i][p] += F12(
              data.position[i][p], data.velocity[i][p],
              data.position[j][q], data.velocity[j][q],
              kSigma, kRadius * 2.);
        }
      }
    }
  }

  for (size_t p = 0; p < data.position[i].size(); ++p) {
    auto& f = data.force[i][p];
    f *= 1. / kMass;
  }

  // environment objects
  for(size_t k=0; k<ENVOBJ.size(); ++k) {
    for (size_t p = 0; p < data.position[i].size(); ++p) {
      auto& obj=ENVOBJ[k];
      //f+=obj->F(p,v,part.r,part.sigma);
      data.force[i][p] += obj->F(
          data.position[i][p], data.velocity[i][p],
          kRadius,kSigma);
    }
  }
}

vect F12(vect p1, vect /*v1*/, vect p2, vect /*v2*/, Scal sigma, Scal R)
{
  //return vect(0.,0.);
  const vect r = p1 - p2;
  const Scal ar2 = r.dot(r);
  const Scal ad2 = 1. / ar2;
  const Scal r2 = ar2 * (1. / (R * R));
  const Scal d2 = ad2 * (R * R);
  const Scal d6 = d2 * d2 * d2;
  const Scal d12 = d6 * d6;
  const Scal cutoff = 2.;
  const Scal cutoff2 = cutoff * cutoff;
  Scal F = 0.;
  if (r2 < cutoff2) {
    F = d2 > cutoff2 ? 0.0 : sigma * (d12 - d6);
    if (r2 > 1.) {
      F *= (cutoff2 - r2) / (cutoff2 - 1.); 
    }
  }
  return r * (F * ad2);
}
