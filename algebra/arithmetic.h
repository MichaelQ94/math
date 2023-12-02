#ifndef HOI_MATH_ALGEBRA_ARITHMETIC_H
#define HOI_MATH_ALGEBRA_ARITHMETIC_H

namespace hoi::algebra {

// Defines static functions allowing retrieval of additive and multiplicative identities
// for an arithmetic type T. Be default, these are assumed to be constructible as
// `T(0)` and `T(1)`, respectively. Custom behavior can be substituted by specializing
// this struct.
template<typename T>
struct Arithmetic {
  // For an arithmetic type T, returns the additive identity satisfying:
  // `zero() + t == t == t + zero()`
  static const T& zero() {
    static const T ZERO(0);
    return ZERO;
  }

  // For an arithmetic type T, returns the multiplicative identity satisfying:
  // `one() * t == t == t * one()`
  static const T& one() {
    static const T ONE(1);
    return ONE;
  }
};

}  // namespace hoi::algebra

#endif
