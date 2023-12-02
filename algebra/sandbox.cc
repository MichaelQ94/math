#include <iostream>
#include <sstream>
#include <string>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "algebra/polynomial.h"

ABSL_FLAG(double, x, 0, "Value at which to evaluate the polynomial");

using ::hoi::algebra::Polynomial;
using ::hoi::algebra::div_mod;

template<typename S>
std::string ToString(const Polynomial<S>& p) {
  std::stringstream ss;
  for (const auto& [deg, coeff] : p.terms_desc()) {
    ss << "(" << deg << ", " << coeff << ") ";
  }
  return ss.str();
}

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);

  Polynomial<double> a;
  a[3] = 1;

  Polynomial<double> b;
  b[2] = 3;

  Polynomial<double> prod = a * b;
  const auto [quot, rem] = div_mod(a, b);

  double x = absl::GetFlag(FLAGS_x);
  double p_x = prod(x);

  std::cout << "prod = " << ToString(prod) << std::endl;
  std::cout << "quot = " << ToString(quot) << std::endl;
  std::cout << "rem = " << ToString(rem) << std::endl;
  std::cout << "qb + r = " << ToString((quot * b) + rem) << std::endl;
  std::cout << "p(" << x << ") = " << p_x << std::endl;
  std::cout << "b(a) = " << ToString(b(a)) << std::endl;

  return 0;
}
