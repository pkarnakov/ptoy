// Created by Petr Karnakov on 25.02.2020

#include <sstream>
#include <stdexcept>
#include <string>

#define FILELINE (std::string() + __FILE__ + ":" + std::to_string(__LINE__))

// force assert
#define fassert_1(x)                                                      \
  do {                                                                    \
    if (!(x)) {                                                           \
      throw std::runtime_error(FILELINE + ": assertion failed '" #x "'"); \
    }                                                                     \
  } while (0);

#define fassert_2(x, msg)                                      \
  do {                                                         \
    if (!(x)) {                                                \
      throw std::runtime_error(                                \
          FILELINE + ": assertion failed '" #x "'\n" + (msg)); \
    }                                                          \
  } while (0);

#define GET_COUNT(_1, _2, _3, COUNT, ...) COUNT
#define VA_SIZE(...) GET_COUNT(__VA_ARGS__, 3, 2, 1, 0)
#define PTOY_CAT(x, y) x##y
#define PTOY_XCAT(x, y) PTOY_CAT(x, y)

#define fassert(...) PTOY_XCAT(fassert##_, VA_SIZE(__VA_ARGS__))(__VA_ARGS__)

#define NAMEVALUE(x)                 \
  ([&]() -> std::string {            \
    std::stringstream namevalue_s;   \
    namevalue_s << #x << '=' << (x); \
    return namevalue_s.str();        \
  }())

#define fassert_equal(x, y)                                                   \
  do {                                                                        \
    if (!((x) == (y))) {                                                      \
      std::stringstream fasrteq_s;                                            \
      fasrteq_s << FILELINE << ": assertion failed, expected equal ";         \
      fasrteq_s << #x << "='" << (x) << "' and " << #y << "='" << (y) << "'"; \
      throw std::runtime_error(fasrteq_s.str());                              \
    }                                                                         \
  } while (0);
