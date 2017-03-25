#include "particles_system.hpp"

particles_system::particles_system() : 
    domain(rect_vect(vect(-1.,-1.),vect(1.,1.))),
    Blocks(domain, vect(4*kRadius, 4*kRadius), 200, false)
{
  force_enabled = false;
  force_center = vect(0., 0.);

  // place particles in the domain
  const double r = kRadius;
  const int row=1. / r;
  const int N=row * row;

  std::vector<particle> P;
  for(int i=0; i<N; ++i)
  {
    const vect p((i%row*2.0+1.0)*r-1.0, (i/row*2.0+1.0)*r-1.0);
    const vect v(0., 0.);
    // TODO: adjust sigma so that with r->0 it converges to a solid body
    const double sigma = kSigma;
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

  Blocks.add_particles(P);
  Blocks.print_status();

  t=0.0;
  dt=0.0002;
  const double gravity = 10.;
  g=vect(0.0, -1.0) * gravity;
}
particles_system::~particles_system()
{}
std::vector<particle> particles_system::GetParticles() const {
  //std::lock_guard<std::mutex> lg(m_step);
  std::vector<particle> res;
  for(int bj=0; bj<Blocks.B.msize().j; ++bj)
  for(int bi=0; bi<Blocks.B.msize().i; ++bi)
  {
    for (auto& part : Blocks.B[mindex(bi, bj)]) {
      res.push_back(part);
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
void particles_system::step(double time_target)
{
  std::lock_guard<std::mutex> lg(m_step);

  #pragma omp parallel
  {

  while (t < time_target) {
    #pragma omp for collapse(2) schedule(dynamic, Blocks.B.msize().i)
    for(int bj=0; bj<Blocks.B.msize().j; ++bj)
    for(int bi=0; bi<Blocks.B.msize().i; ++bi)
    {
      RHS(mindex(bi,bj));
    }

    #pragma omp for collapse(2) schedule(dynamic, Blocks.B.msize().i)
    for(int bj=0; bj<Blocks.B.msize().j; ++bj)
    for(int bi=0; bi<Blocks.B.msize().i; ++bi)
    {
      mindex m(bi,bj);
      std::vector<particle>& b1=Blocks.B[m];
      for (auto& part : b1) {
        part.p0 = part.p;
        part.v0 = part.v;
        part.p += part.v * dt * 0.5;
        part.v += part.f * dt * 0.5;
      }
    }

    #pragma omp single
    t += 0.5 * dt;

    #pragma omp for collapse(2) schedule(dynamic, Blocks.B.msize().i)
    for(int bj=0; bj<Blocks.B.msize().j; ++bj)
    for(int bi=0; bi<Blocks.B.msize().i; ++bi)
    {
      RHS(mindex(bi,bj));
    }

    #pragma omp for collapse(2) schedule(dynamic, Blocks.B.msize().i)
    for(int bj=0; bj<Blocks.B.msize().j; ++bj)
    for(int bi=0; bi<Blocks.B.msize().i; ++bi)
    {
      mindex m(bi,bj);
      std::vector<particle>& b1=Blocks.B[m];
      for (auto& part : b1) {
        part.v = part.v0 + part.f * dt;
        const double limit = 10.;
        if (part.v.length() > limit) {
          part.v *= limit / part.v.length();
        }
        part.p = part.p0 + part.v * dt;
      }
    }

    #pragma omp single
    {
      t += 0.5 * dt;
      Blocks.arrange();
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

void particles_system::RHS(mindex block_m)
{
  const mindex m = block_m;

  std::vector<particle>& b1=Blocks.B[m];
  for(size_t w1=0; w1<b1.size(); ++w1)
  {
    particle& p1 = b1[w1]; // first particle
    auto& part = p1;
    auto& f = part.f;
    auto& p = part.p;
    auto& v = part.v;

    // zero
    f = vect(0., 0.);

    // gravity
    //f=g*part.m;
    f=g*kMass;

    // point force
    if (force_enabled) {
        const double intensity = 0.1;
        const vect r = p - force_center; 
        f += r * (intensity / std::pow(r.length(), 3));
    }

    // dissipation
    //f -= v * (0.1 * p1.m);
    f -= v * (0.1 * kMass);
  }

  // pairwise interactions
  for(size_t k=0; k<Blocks.NEAR.size(); ++k)
  {
    mindex mnear = m + Blocks.NEAR[k];
    if(Blocks.B.valid(mnear))
    {
      std::vector<particle>& b2=Blocks.B[mnear];

      for(size_t w1=0; w1<b1.size(); ++w1)
      {
        particle& p1 = b1[w1]; // first particle

        for(size_t w2=0; w2<b2.size(); ++w2)
        {
          const particle& p2 = b2[w2]; // second particle
          //if(&p1!=&p2 && (p1.layers_mask & p2.layers_mask))
          if (&p1!=&p2) {
            p1.f += F12(p1.p, p1.v, p2.p, p2.v, kSigma, kRadius * 2.);
          }
        }
      }
    }
  }

  for(size_t w1=0; w1<b1.size(); ++w1)
  {
    // mass
    //f *= 1. / p1.m;
    b1[w1].f *= 1. / kMass;
  }

  if (m.i == 0 || m.i == Blocks.N.i - 2 ||
      m.j == 0 || m.j == Blocks.N.j - 2)
  for(size_t w1=0; w1<b1.size(); ++w1)
  {
    particle& p1 = b1[w1]; 
    // environment objects
    {
     // std::lock_guard<std::mutex> lg(m_ENVOBJ);
      for(size_t k=0; k<ENVOBJ.size(); ++k)
      {
        auto& obj=ENVOBJ[k];
        //f+=obj->F(p,v,part.r,part.sigma);
        p1.f+=obj->F(p1.p,p1.v,kRadius,kSigma);
      }
    }
  }
}

vect F12(vect p1, vect /*v1*/, vect p2, vect /*v2*/, double sigma, double R)
{
  const vect r = p1 - p2;
  const double ar2 = r.dot(r);
  const double ad2 = 1. / ar2;
  const double r2 = ar2 * (1. / (R * R));
  const double d2 = ad2 * (R * R);
  const double d6 = d2 * d2 * d2;
  const double d12 = d6 * d6;
  const double cutoff = 2.;
  const double cutoff2 = cutoff * cutoff;
  double F = 0.;
  if (r2 < cutoff2) {
    F = d2 > cutoff2 ? 0.0 : sigma * (d12 - d6);
    if (r2 > 1.) {
      F *= (cutoff2 - r2) / (cutoff2 - 1.); 
    }
  }
  return r * (F * ad2);
}
