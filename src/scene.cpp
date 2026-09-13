#include "scene.h"

#include "particles.h"

void UpdateScene(const Particles& partsys, SceneData& data, Scene& scene) {
  { // Particles.
    const auto& particles = partsys.GetParticles();
    data.position.resize(particles.size());
    data.velocity.resize(particles.size());
    for (size_t i = 0; i < particles.size(); ++i) {
      data.position[i] = particles[i].p;
      data.velocity[i] = particles[i].v;
    }
    scene.particles.p = data.position;
    scene.particles.v = data.velocity;
  }

  { // Portals.
    using Portal = Scene::Portal;
    using Pair = std::array<Portal, 2>;
    const auto& portals = partsys.GetPortals();
    data.portals.resize(portals.size());
    for (size_t i = 0; i < portals.size(); ++i) {
      data.portals[i] = Pair{
          Portal{portals[i][0].begin, portals[i][0].end},
          Portal{portals[i][1].begin, portals[i][1].end},
      };
    }
    // Append the incomplete pair currently drawn, leaving the portal that is
    // not drawn yet empty.
    const auto drawing = partsys.GetPortalDrawing();
    const Portal empty{Vect(0), Vect(0)};
    if (drawing.stage == 0) {
      if (drawing.mouse_moving) {
        data.portals.emplace_back(
            Pair{Portal{drawing.begin, drawing.current}, empty});
      }
    } else {
      data.portals.emplace_back(
          Pair{Portal{drawing.prev.first, drawing.prev.second}, empty});
      if (drawing.mouse_moving) {
        data.portals.back()[1] = Portal{drawing.begin, drawing.current};
      }
    }
    scene.portals = data.portals;
  }

  { // Bonds.
    data.bonds.clear();
    const auto& norend = partsys.GetNoRendering();
    const auto& pos = partsys.GetBlockData().position;
    const auto& bbi = partsys.GetBlockById();
    for (auto bond : partsys.GetBonds()) {
      const auto& a = bbi[bond.first];
      const auto& b = bbi[bond.second];
      if (!norend.count(bond)) {
        data.bonds.push_back({pos[a.first][a.second], pos[b.first][b.second]});
      }
    }
    scene.bonds = data.bonds;
  }

  { // Frozen particles.
    data.frozen.clear();
    const auto& pos = partsys.GetBlockData().position;
    const auto& bbi = partsys.GetBlockById();
    for (auto id : partsys.GetFrozen()) {
      const auto& a = bbi[id];
      data.frozen.push_back(pos[a.first][a.second]);
    }
    scene.frozen = data.frozen;
  }

  scene.domain = partsys.GetDomain();
}
