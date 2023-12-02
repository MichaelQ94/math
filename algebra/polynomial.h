#ifndef HOI_MATH_ALGEBRA_POLYNOMIAL_H
#define HOI_MATH_ALGEBRA_POLYNOMIAL_H

#include <concepts>
#include <iterator>
#include <map>
#include <type_traits>
#include <utility>

#include "algebra/arithmetic.h"
#include "absl/log/check.h"

namespace hoi::algebra {

// Represents a polynomial of the form `p(x) = p[0] + p[1]*x + p[2]*x^2 + ...`
// A default-constructed instance represents the zero polynomial `p(x) = 0`.
// Coefficients can be accessed and mutated using `operator[]`. Any unset
// coefficient is implicitly 0, and any non-zero coefficient which is replaced
// with a value of 0 is removed from the internal representation. A Polynomial
// can be evaluated as a function using `operator()`, which substitutes the
// given argument into the encapsulated polynomial expression.
template<typename S>
class Polynomial {
 public:
  template<typename BackingIt> class Terms;
  using TermsAsc = Terms<typename std::map<int, S>::const_iterator>;
  using TermsDesc = Terms<typename std::map<int, S>::const_reverse_iterator>;
  template<typename BackingIt> class ConstTermIterator;
  class CoeffRef;

  // A default-constructed Polynomial represents `p(x) = 0`.
  Polynomial() = default;
  Polynomial(const Polynomial&) = default;
  Polynomial(Polynomial&&) = default;

  // Allows implicit creation of constant polynomials `p(x) = c`.
  Polynomial(const S& c);
  Polynomial(S&& c);

  // Create a polynomial representing `p(x) = a*(x^n)`.
  Polynomial(int n, const S& a);
  Polynomial(int n, S&& a);

  bool operator==(const Polynomial&) const = default;
  bool operator!=(const Polynomial&) const = default;

  Polynomial& operator=(const Polynomial&) = default;
  Polynomial& operator=(Polynomial&&) = default;

  // Accesses the coefficient of the degree `deg` term (`deg` must be positive). The
  // returned CoeffRef is implicitly convertible to `const S&` if value inspection is
  // desired, and also has an overloaded `operator=(const S&)` which can be used to
  // set the value of the referenced coefficient.
  CoeffRef operator[](int deg);

  // Returns the coefficient of the degree `deg` term (`deg` must be positive).
  const S& operator[](int deg) const;

  // Returns the number of nonzero terms in the polynomial.
  int num_terms() const;

  // Returns an object which can be ued to iterate over the terms in this polynomial,
  // in ascending order by degree.
  TermsAsc terms_asc() const;

  // Returns an object which can be ued to iterate over the terms in this polynomial,
  // in descending order by degree.
  TermsDesc terms_desc() const;

  // Substitutes the given value of `x` into the encapsulated polynomial expression
  // and computes the result. `X` must be constructible from `S`.
  template<typename X> requires std::convertible_to<S, X>
  X operator()(const X& x) const;

  // Arithmetic operators.
  template<typename T>
  friend Polynomial<T> operator+(const Polynomial<T>& lhs, const Polynomial<T>& rhs);

  template<typename T>
  friend Polynomial<T> operator-(const Polynomial<T>& lhs, const Polynomial<T>& rhs);

  template<typename T>
  friend Polynomial<T> operator*(const Polynomial<T>& lhs, const Polynomial<T>& rhs);

  template<typename T>
  friend Polynomial<T> operator/(const Polynomial<T>& lhs, const Polynomial<T>& rhs);

  template<typename T>
  friend Polynomial<T> operator%(const Polynomial<T>& lhs, const Polynomial<T>& rhs);

  // Proxy object allowing polynomial coeffecients to be accessed and assigned,
  // which is intended to be used as though it is of type `S&`.
  //
  // Implementation detail: reassigning coefficient values sometimes also requires
  // modifying the internal representation of the host Polynomial. See the
  // implementation of operator= for details.
  class CoeffRef {
    using const_ref = const S&;

   public:
    // Assign/reassign the referenced coefficient.
    CoeffRef& operator=(const S& s);
    CoeffRef& operator=(S&& s);

    // Allow implicit conversion to `const S&`, so that coefficients can be
    // inspected using:
    // ```
    // Polynomial<int> p;
    // S p_0 = p[0];
    // ```
    operator const_ref() const;

   private:
    friend class Polynomial;

    CoeffRef(const int deg, Polynomial& host);

    // Used to implement assignment operators.
    template<typename T> requires std::same_as<std::remove_cvref_t<T>, std::remove_cvref_t<S>>
    CoeffRef& assign_coeff(T&& t);  // T&& is a forwarding reference

    const int deg_;
    Polynomial& host_;
  };

  // Provides read access to the non-zero monomial terms in a Polynomial.
  struct ConstMonomialRef {
    const int deg;
    const S& coeff;
  };

  class ArrowHelper {
   public:
    ConstMonomialRef* operator->() { return &ref_; }

   private:
    template<typename BackingIt> friend class ConstTermIterator;
  
    explicit ArrowHelper(ConstMonomialRef&& ref) : ref_(std::move(ref)) {}

    ConstMonomialRef ref_;
  };

  template<typename BackingIt>
  class ConstTermIterator {
   public:
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::iterator_traits<BackingIt>::difference_type;
    using value_type = ConstMonomialRef;
    using pointer = value_type*;
    using reference = value_type&;

    ConstTermIterator() = default;
    ConstTermIterator(const ConstTermIterator&) = default;
    ConstTermIterator(ConstTermIterator&&) = default;
    ConstTermIterator& operator=(const ConstTermIterator&) = default;
    ConstTermIterator& operator=(ConstTermIterator&&) = default;

    ConstMonomialRef operator*() const {
      return ConstMonomialRef{.deg = backing_it_->first, .coeff = backing_it_->second};
    }

    ArrowHelper operator->() const {
      return ArrowHelper(**this);
    }

    bool operator==(const ConstTermIterator& it) const = default;
    bool operator!=(const ConstTermIterator& it) const = default;

    ConstTermIterator& operator++() {
      ++backing_it_;
      return *this;
    }

    ConstTermIterator operator++(int) {
      ConstTermIterator cur = *this;
      ++(*this);
      return cur;
    }

   private:
    friend class Polynomial;

    explicit ConstTermIterator(const BackingIt& backing_it) : backing_it_(backing_it) {}

    BackingIt backing_it_;
  };

  template<typename BackingIt>
  class Terms {
   public:
    ConstTermIterator<BackingIt> begin() const {
      return ConstTermIterator<BackingIt>(backing_begin_);
    }

    ConstTermIterator<BackingIt> end() const {
      return ConstTermIterator<BackingIt>(backing_end_);
    }

   private:
    friend class Polynomial;

    Terms(BackingIt&& backing_begin, BackingIt&& backing_end)
        : backing_begin_(std::move(backing_begin))
        , backing_end_(std::move(backing_end)) {}

    const BackingIt backing_begin_;
    const BackingIt backing_end_;
  };

 private:
  std::map<int, S> deg_to_coeff_;
};


// Returns the degree of the highest-order term in the polynomial with a non-zero
// coefficient. In the special case of the zero polynomial, this will return -1.
// For any two polynomials `p` and `q`, `deg(p + q) == max(deg(p), deg(q))`.
// For any two *non-zero* polynomials `p` and `q`, `deg(p * q) == deg(p) + deg(q)`.
template<typename S>
int deg(const Polynomial<S>& p) {
  return p.num_terms() > 0 ? p.terms_desc().begin()->deg : -1;
}

// Computes the unique polynomials `q` and `r` satisfying `(a == (q * b) + r) &&
// (deg(r) < deg(b))`, and returns the pair `(q, r)`. `b` must be nonzero.
template<typename S>
std::pair<Polynomial<S>, Polynomial<S>> div_mod(const Polynomial<S>& a, const Polynomial<S>& b) {
  const int deg_b = deg(b);
  CHECK_GE(deg_b, 0) << "Cannot divide by zero polynomial.";

  std::pair<Polynomial<S>, Polynomial<S>> result;
  auto& [quot, rem] = result;
  rem = a;

  for (int deg_rem = deg(rem); deg_rem >= deg_b; deg_rem = deg(rem)) {
    const int deg_diff = deg_rem - deg_b;
    const S leading_coeff_ratio = rem[deg_rem] / b[deg_b];

    quot[deg_diff] = leading_coeff_ratio;
    rem = rem - (Polynomial(deg_diff, leading_coeff_ratio) * b);
    rem[deg_rem] = Arithmetic<S>::zero();  // Ensure exact cancellation of the leading-order term.
  }
  return result;
}

template<typename S>
Polynomial<S>::Polynomial(const S& c) : Polynomial(0, c) {}

template<typename S>
Polynomial<S>::Polynomial(S&& c) : Polynomial(0, std::move(c)) {}

template<typename S>
Polynomial<S>::Polynomial(int n, const S& a) {
  if (a != Arithmetic<S>::zero()) {
    (*this)[n] = a;
  }
}

template<typename S>
Polynomial<S>::Polynomial(int n, S&& a) {
  if (a != Arithmetic<S>::zero()) {
    (*this)[n] = std::move(a);
  }
}

template<typename S>
Polynomial<S>::CoeffRef Polynomial<S>::operator[](int deg) {
  CHECK_GE(deg, 0) << "Polynomial terms must have non-negative degree.";
  return CoeffRef(deg, *this);
}

template<typename S>
const S& Polynomial<S>::operator[](int deg) const {
  CHECK_GE(deg, 0) << "Polynomial terms must have non-negative degree.";
  auto it = deg_to_coeff_.find(deg);
  return it != deg_to_coeff_.end() ? it->second : Arithmetic<S>::zero();
}

template<typename S>
int Polynomial<S>::num_terms() const {
  return deg_to_coeff_.size();
}

template<typename S>
Polynomial<S>::TermsAsc Polynomial<S>::terms_asc() const {
  return TermsAsc(deg_to_coeff_.begin(), deg_to_coeff_.end());
}

template<typename S>
Polynomial<S>::TermsDesc Polynomial<S>::terms_desc() const {
  return TermsDesc(deg_to_coeff_.rbegin(), deg_to_coeff_.rend());
}


template<typename S>
template<typename X> requires std::convertible_to<S, X>
X Polynomial<S>::operator()(const X& x) const {
  int n = 0;
  X x_nth_power = Arithmetic<X>::one();
  X result = Arithmetic<X>::zero();

  for (const auto& [deg, coeff] : this->terms_asc()) {
    for (; n < deg; ++n) {
      x_nth_power = x_nth_power * x;
    }
    result = result + (static_cast<const X&>(coeff) * x_nth_power);
  }
  return result;
}

template<typename S>
Polynomial<S> operator+(const Polynomial<S>& lhs, const Polynomial<S>& rhs) {
  Polynomial<S> sum;
  for (const auto& [deg, coeff] : lhs.terms_asc()) {
    sum[deg] = coeff;
  }
  for (const auto& [deg, coeff] : rhs.terms_asc()) {
    sum[deg] = lhs[deg] + coeff;
  }
  return sum;
}

template<typename S>
Polynomial<S> operator-(const Polynomial<S>& lhs, const Polynomial<S>& rhs) {
  Polynomial<S> diff;
  for (const auto& [deg, coeff] : lhs.terms_asc()) {
    diff[deg] = coeff;
  }
  for (const auto& [deg, coeff] : rhs.terms_asc()) {
    diff[deg] = lhs[deg] - coeff;
  }
  return diff;
}

template<typename S>
Polynomial<S> operator*(const Polynomial<S>& lhs, const Polynomial<S>& rhs) {
  Polynomial<S> prod;
  for (const auto& [deg_l, coeff_l] : lhs.terms_asc()) {
    for (const auto& [deg_r, coeff_r] : rhs.terms_asc()) {
      prod[deg_l + deg_r] = prod[deg_l + deg_r] + (coeff_l * coeff_r);
    }
  }
  return prod;
}

template<typename S>
Polynomial<S> operator/(const Polynomial<S>& lhs, const Polynomial<S>& rhs) {
  return div_mod(lhs, rhs).first;
}

template<typename S>
Polynomial<S> operator%(const Polynomial<S>& lhs, const Polynomial<S>& rhs) {
  return div_mod(lhs, rhs).second;
}

template<typename S>
Polynomial<S>::CoeffRef::CoeffRef(const int deg, Polynomial<S>& host)
    : deg_(deg), host_(host) {}

template<typename S>
Polynomial<S>::CoeffRef& Polynomial<S>::CoeffRef::operator=(const S& s) {
  return assign_coeff(s);
}

template<typename S>
Polynomial<S>::CoeffRef& Polynomial<S>::CoeffRef::operator=(S&& s) {
  return assign_coeff(s);
}

template<typename S>
Polynomial<S>::CoeffRef::operator Polynomial<S>::CoeffRef::const_ref() const {
  return static_cast<const Polynomial<S>&>(host_)[deg_];
}

template<typename S>
template<typename T> requires std::same_as<std::remove_cvref_t<T>, std::remove_cvref_t<S>>
Polynomial<S>::CoeffRef& Polynomial<S>::CoeffRef::assign_coeff(T&& t) {
  // By ensuring that the host Polynomial's deg_to_coeff_ map contains only
  // entries describing non-zero coefficients, we can compute the degree of
  // the polynomial in O(1) time by simply extracting the highest key.
  // Exposing direct references of type `S&` to map values would not allow us
  // to do this.
  if (t == Arithmetic<S>::zero()) {
    host_.deg_to_coeff_.erase(deg_);
  } else {
    host_.deg_to_coeff_.insert_or_assign(deg_, std::forward<T>(t));
  }
  return *this;
}

}  // namespace hoi::algebra

#endif
