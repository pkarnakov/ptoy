#include "particles_system.hpp"


particles_system::particles_system() : Blocks(rect_vect(vect(-1.,-1.),vect(1.,1.)), vect(0.2, 0.2), 100)
{
  INT=std::unique_ptr<integrator>(new integrator_Euler);
   
  force_enabled = false;
  force_center = vect(0., 0.);

  // place particles in the domain
  double r=0.02;
  int N=100;
  for(int i=0; i<N; ++i)
  {
    int row=11;
    P.push_back(particle(vect((i%row*2.0+1.0)*r-1.0, (i/row*2.0+1.0)*r-1.0), vect(0.1, 0.0), 0.01, r, 10000.0, 1+i%2*0, rgb(1.0, i%2, 0.0)));
  }

  transfer_data(P, X,V);

  t=0.0;
  dt=0.0001;
  g=vect(0.0, -10.0);

  INT->set_data(X, V);
  INT->set_t(t);
  INT->set_dt(dt);

  using std::placeholders::_1;
  using std::placeholders::_2;
  using std::placeholders::_3;
  using std::placeholders::_4;
  INT->set_RHS(std::bind(&particles_system::RHS, this, _1, _2, _3, _4));
}
particles_system::~particles_system()
{}
std::vector<particle> particles_system::GetParticles() const {
  return P;
}
void particles_system::AddEnvObj(env_object* env) {
  ENVOBJ.push_back(std::unique_ptr<env_object>(env));
}
void particles_system::status(std::ostream& out)
{
  out<<"Particles system"<<std::endl<<"Particles number = "<<P.size()<<std::endl;
}
void particles_system::step()
{
  std::vector<vect> X(P.size()), V(P.size()), F(P.size());
  for(std::size_t i=0; i<P.size(); ++i)
  {
    auto& p = P[i].p;
    auto& v = P[i].v;
    p += v * dt;
  }

  for(std::size_t i=0; i<P.size(); ++i)
  {
    X[i] = P[i].p;
    V[i] = P[i].v;
  }

  RHS(X, V, t, F);

  for(std::size_t i=0; i<P.size(); ++i)
  {
    auto& v = P[i].v;
    auto& f = F[i];
    v += (f + v * (-0.1)) * dt;
    const double limit = 100.;
    if (v.length() > limit) {
      v *= limit / v.length();
    }
  }

  Blocks.arrange(X);
  t += dt;
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

void particles_system::transfer_data(const vector<particle>& P, vector<vect>& X, vector<vect>& V)
{
  size_t N=P.size();
  X.resize(N);
  V.resize(N);
  for(size_t i=0; i<N; ++i)
  {
    X[i]=P[i].p;
    V[i]=P[i].v;
  }
}

void particles_system::RHS(const vector<vect>& X, const vector<vect>& V, double /*t*/, vector<vect>& F) const
{
  F.resize(X.size());

  for(size_t i=0; i<F.size(); ++i)
  {
    // gravity
    F[i]=g*P[i].m;
    // environment objects
    for(size_t k=0; k<ENVOBJ.size(); ++k)
    {
      auto& obj=ENVOBJ[k];
      F[i]+=obj->F(X[i],V[i],P[i].r,P[i].sigma);
    }
    // point force
    if (force_enabled) {
        const double intensity = 0.1;
        const vect r = X[i] - force_center; 
        F[i] += r * (intensity / std::pow(r.length(), 3));
    }

    // dissipation
    F[i] -= V[i] * 0.1;
  }

  // pairwise interactions
  for(int bj=0; bj<Blocks.B.msize().j; ++bj)
  for(int bi=0; bi<Blocks.B.msize().i; ++bi)
  {
    mindex m(bi,bj);
    const std::vector<int>& b1=Blocks.B[m];
    for(size_t w1=0; w1<b1.size(); ++w1)
    {
      int p1=b1[w1]; // first particle

      for(size_t k=0; k<Blocks.NEAR.size(); ++k)
      {
        mindex mnear=m+Blocks.NEAR[k];
        if(Blocks.B.valid(mnear))
        {
          const std::vector<int>& b2=Blocks.B[mnear];
          for(size_t w2=0; w2<b2.size(); ++w2)
          {
            int p2=b2[w2]; // second particle
            if(p1!=p2 && (P[p1].layers_mask & P[p2].layers_mask))
            {
              F[p1]+=F12(X[p1], V[p1], X[p2], V[p2], 0.5*(P[p1].sigma+P[p2].sigma), P[p1].r+P[p2].r);
            }
          }
        }
      }
    }
  }

  for(size_t i=0; i<F.size(); ++i)
  {
    F[i]*=1.0/P[i].m;
  }
}

vect F12(vect p1, vect /*v1*/, vect p2, vect /*v2*/, double sigma, double R)
{
  const double alpha=12.0;
  const double beta=6.0;
   
  double eps=sigma/(pow(2.0,alpha)-pow(2.0,beta));

  double r=p1.dist(p2);
  double F=r>R?0.0:eps*(pow(R/r, alpha)-pow(R/r, beta));
  return (p1-p2)*(F/r);
}
