#pragma once

#include <string>

#include "geometry.h"
#include "particles.h"

class Control {
 public:
  using Keysym = char;

  enum class EventType {
    key_down,
    mouse_motion,
    mouse_down,
    mouse_up,
  };

  enum class MouseMode {
    None,
    Attraction,
    Repulsion,
    Bond,
    Pick,
    Freeze,
    Portal,
  };

  struct Event {
    EventType type;
    Keysym keysym; // ASCII value.
    Vect mousepos; // Mouse position in domain units.
  };

  Control(Particles*);
  void Send(Event);
  void SendKeyDown(Keysym);
  void SendMouseMotion(Vect mousepos);
  void SendMouseDown(Vect mousepos);
  void SendMouseUp(Vect mousepos);
  void SetMouseMode(MouseMode s);

  static const char* EventTypeToStr(EventType type);
  static const char* MouseModeToStr(MouseMode s);

 public:
  bool debug = false;
  MouseMode mouse_mode = MouseMode::Repulsion;

 private:
  void Handle(Event);

  Particles* partsys_;
};
