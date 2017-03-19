#include "particles_system.hpp"


particles_system::particles_system() : Blocks(rect_vect(vect(-1.,-1.),vect(1.,1.)), vect(0.1, 0.1), 200, false)
{
  force_enabled = false;
  force_center = vect(0., 0.);

  // place particles in the domain
  double r=0.02;
  int N=1000;

  std::vector<particle> P;
  for(int i=0; i<N; ++i)
  {
    int row=37;
    const vect p((i%row*2.0+1.0)*r-1.0, (i/row*2.0+1.0)*r-1.0);
    const vect v(0., 0.);
    const double sigma = 100000.;
    switch (i % 3) {
      case 0:
        P.push_back(particle(p, v, 0.01, r, sigma, 0x3, rgb(1.,0.,0.)));
        break;
      case 1:
        P.push_back(particle(p, v, 0.02, r, sigma, 0x3, rgb(0.,1.,0.)));
        break;
      case 2:
        P.push_back(particle(p, v, 0.05, r, sigma, 0x3, rgb(0.,0.,1.)));
        break;
    }
  }

  Blocks.add_particles(P);
  Blocks.print_status();

  t=0.0;
  dt=0.0003;
  g=vect(0.0, -10.0);
}
particles_system::~particles_system()
{}
std::vector<particle> particles_system::GetParticles() const {
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
void particles_system::step()
{
  RHS();

  #pragma omp parallel for collapse(2) schedule(dynamic, Blocks.B.msize().i)
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
  const auto t0 = t;
  t += 0.5 * dt;

  RHS();

  #pragma omp parallel for collapse(2) schedule(dynamic, Blocks.B.msize().i)
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
  t = t0 + dt;

  Blocks.arrange();
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

void particles_system::RHS()
{
  #pragma omp parallel for collapse(2) schedule(dynamic, Blocks.B.msize().i)
  for(int bj=0; bj<Blocks.B.msize().j; ++bj)
  for(int bi=0; bi<Blocks.B.msize().i; ++bi)
  {
    mindex m(bi,bj);
    std::vector<particle>& b1=Blocks.B[m];
    for(size_t w1=0; w1<b1.size(); ++w1)
    {
      particle& p1 = b1[w1]; // first particle
      auto& part = p1;
      auto& f = part.f;
      auto& p = part.p;
      auto& v = part.v;
      // gravity
      f=g*part.m;
      // environment objects
      for(size_t k=0; k<ENVOBJ.size(); ++k)
      {
        auto& obj=ENVOBJ[k];
        f+=obj->F(p,v,part.r,part.sigma);
      }
      // point force
      if (force_enabled) {
          const double intensity = 0.1;
          const vect r = p - force_center; 
          f += r * (intensity / std::pow(r.length(), 3));
      }

      // dissipation
      f -= v * (0.1 * p1.m);

      // pairwise interactions
      for(size_t k=0; k<Blocks.NEAR.size(); ++k)
      {
        mindex mnear = m + Blocks.NEAR[k];
        if(Blocks.B.valid(mnear))
        {
          std::vector<particle>& b2=Blocks.B[mnear];
          for(size_t w2=0; w2<b2.size(); ++w2)
          {
            particle& p2=b2[w2]; // second particle
            if(&p1!=&p2 && (p1.layers_mask & p2.layers_mask))
            {
              f += F12(p1.p, p1.v, p2.p, p2.v, 0.5*(p1.sigma+p2.sigma), p1.r+p2.r);
            }
          }
        }
      }
      
      // mass
      f *= 1. / p1.m;
    }
  }
}

vect F12(vect p1, vect /*v1*/, vect p2, vect /*v2*/, double sigma, double R)
{
  const double alpha=12.0;
  const double beta=0.0;
   
  double eps=sigma/(pow(2.0,alpha)-pow(2.0,beta));

  double r=p1.dist(p2);
  double F=r>R?0.0:eps*(pow(R/r, alpha)-pow(R/r, beta));
  return (p1-p2)*(F/r);
}
