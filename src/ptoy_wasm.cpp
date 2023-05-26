#include <stdint.h>
#include <stdlib.h>
#include <emscripten.h>
#include <emscripten/html5.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>

static void main_loop() {
  EM_ASM_({ Draw(); });
}

extern "C" {
int SetConfig(const char* str) {
  std::stringstream buf(str);
  try {
    throw std::runtime_error("SetConfig()");
  } catch (const std::runtime_error& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  return 0;
}
const char* GetConfig() {
  return "GetConfig()";
}
}

int main() {
  emscripten_set_main_loop(main_loop, 1, 0);
}
