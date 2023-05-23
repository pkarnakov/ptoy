# New Design

## MVC

* model: particle system
* view: opengl window
* controller: user interface

## Interaction

* model to view: provide element positions (particles, bonds)
* view to model: nothing
* view to controller: pass mouse and keyboard events
* controller to view: provide geometry of control elements (buttons)
*  controller to model: pass user commands
* model to controller: provide particle properties

```
using Idx = size_t;
using IdxInt = int;
using Scal = float;
using Vect = geom::Vect<Scal, 2>;
using Id = size_t;
using ParticleId = Id;

struct Particle {
  Vect position;
  Vect velocity;
  Scal mass;
  // Family family;
};

struct Bond {
  ParticleId first, second;
};

struct VisibleState {
  std::vector<Particle> particles;
  std::vector<Bond> bonds;
};

class Model {
 public:
  void GetVisibleState();
};

class Interface {
 public:
  Interface(Model* model);
  virtual void ForceRedraw() = 0;
};
```

## Particles

Evolution of particles.

Getters:
- positions and velocity of particles
- environment objects
- constaints and bonds
- custom fields
- domain size

Modifiers:
- add particles at given positions
- resize domain
- add environment object

Internal state:
- cell lists
