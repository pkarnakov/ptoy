#include <iostream>

#include "control.h"

const char* Control::EventTypeToStr(EventType type) {
  switch (type) {
    case EventType::key_down:
      return "key_down";
    case EventType::mouse_motion:
      return "mouse_motion";
    case EventType::mouse_down:
      return "mouse_down";
    case EventType::mouse_up:
      return "mouse_up";
    default:
      return "";
  }
}

const char* Control::MouseModeToStr(MouseMode s) {
  switch (s) {
    case MouseMode::None:
      return "none";
    case MouseMode::Attraction:
      return "attraction";
    case MouseMode::Repulsion:
      return "repulsion";
    case MouseMode::Bond:
      return "bonds";
    case MouseMode::Pick:
      return "pick";
    case MouseMode::Freeze:
      return "freeze";
    case MouseMode::Portal:
      return "portal";
    default:
      return nullptr;
  }
}

Control::Control(Particles* partsys_) : partsys_(partsys_) {}

void Control::Send(Event e) {
  Handle(e);
  if (!debug) {
    return;
  }
  const char keysym = std::isprint(e.keysym) ? e.keysym : ' ';
  std::cout << EventTypeToStr(e.type);
  switch (e.type) {
    case EventType::key_down:
      std::cout << " keysym='" << keysym << "'"
                << " dec=" << int(e.keysym);
      break;
    case EventType::mouse_motion:
    case EventType::mouse_down:
    case EventType::mouse_up:
      std::cout << " mousepos=" << e.mousepos;
      break;
    default:
      break;
  }
  std::cout << std::endl;
}

void Control::SetMouseMode(MouseMode s) {
  mouse_mode = s;
  std::cout << std::string("Mouse switched to mode: ") + MouseModeToStr(s)
            << std::endl;
  switch (s) {
    case MouseMode::Repulsion:
      partsys_->SetForceAttractive(false);
      break;
    case MouseMode::Attraction:
      partsys_->SetForceAttractive(true);
      break;
    default:
      break;
  }
}

void Control::Handle(Event e) {
  if (e.type == EventType::key_down) {
    switch (e.keysym) {
      case 'g':
        partsys_->SetGravity(!partsys_->GetGravity());
        std::cout << (partsys_->GetGravity() ? "Gravity on" : "Gravity off")
                  << std::endl;
        break;
      case 'n':
        SetMouseMode(MouseMode::None);
        break;
      case 'r':
        SetMouseMode(MouseMode::Repulsion);
        break;
      case 'f':
        SetMouseMode(MouseMode::Freeze);
        break;
      case 'p':
        SetMouseMode(MouseMode::Pick);
        break;
      case 'a':
        SetMouseMode(MouseMode::Attraction);
        break;
      case 'b':
        SetMouseMode(MouseMode::Bond);
        break;
      case 'o':
        SetMouseMode(MouseMode::Portal);
        break;
      case 'i':
        std::cout << "Remove last portal" << std::endl;
        partsys_->RemoveLastPortal();
        break;
      default:
        break;
    }
  } else if (e.type == EventType::mouse_motion) {
    switch (mouse_mode) {
      case MouseMode::Attraction:
      case MouseMode::Repulsion:
        partsys_->SetForce(e.mousepos);
        break;
      case MouseMode::Bond:
        partsys_->BondsMove(e.mousepos);
        break;
      case MouseMode::Pick:
        partsys_->PickMove(e.mousepos);
        break;
      case MouseMode::Freeze:
        partsys_->FreezeMove(e.mousepos);
        break;
      case MouseMode::Portal:
        partsys_->PortalMove(e.mousepos);
        break;
      case MouseMode::None:
        break;
      default:
        break;
    }
  } else if (e.type == EventType::mouse_down) {
    switch (mouse_mode) {
      case MouseMode::Attraction:
      case MouseMode::Repulsion:
        partsys_->SetForce(e.mousepos, true);
        break;
      case MouseMode::Bond:
        partsys_->BondsStart(e.mousepos);
        break;
      case MouseMode::Pick:
        partsys_->PickStart(e.mousepos);
        break;
      case MouseMode::Freeze:
        partsys_->FreezeStart(e.mousepos);
        break;
      case MouseMode::Portal:
        partsys_->PortalStart(e.mousepos);
        break;
      case MouseMode::None:
        break;
      default:
        break;
    }
  } else if (e.type == EventType::mouse_up) {
    switch (mouse_mode) {
      case MouseMode::Attraction:
      case MouseMode::Repulsion:
        partsys_->SetForce(false);
        break;
      case MouseMode::Bond:
        partsys_->BondsStop(e.mousepos);
        break;
      case MouseMode::Pick:
        partsys_->PickStop(e.mousepos);
        break;
      case MouseMode::Freeze:
        partsys_->FreezeStop(e.mousepos);
        break;
      case MouseMode::Portal:
        partsys_->PortalStop(e.mousepos);
        break;
      case MouseMode::None:
        break;
      default:
        break;
    }
  }
}

void Control::SendKeyDown(Keysym keysym) {
  Event e;
  e.type = EventType::key_down;
  e.keysym = keysym;
  Send(e);
}

void Control::SendMouseMotion(Vect mousepos) {
  Event e;
  e.type = EventType::mouse_motion;
  e.mousepos = mousepos;
  Send(e);
}

void Control::SendMouseDown(Vect mousepos) {
  Event e;
  e.type = EventType::mouse_down;
  e.mousepos = mousepos;
  Send(e);
}

void Control::SendMouseUp(Vect mousepos) {
  Event e;
  e.type = EventType::mouse_up;
  e.mousepos = mousepos;
  Send(e);
}
