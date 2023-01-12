#pragma once

#include "./GreatestCommonDivisor.hpp"
#include "./Math.hpp"
#include <cstddef>
#include <cstdint>
#include <llvm/ADT/SmallVector.h>
#include <llvm/Support/raw_ostream.h>
#include <optional>

template <class T, int Bits>
concept is_int_v = std::signed_integral<T> && sizeof(T) == (Bits / 8);

template <is_int_v<32> T> constexpr auto widen(T x) -> int64_t { return x; }
template <is_int_v<64> T> constexpr auto widen(T x) -> __int128_t { return x; }
template <is_int_v<32> T> constexpr auto splitInt(T x) -> int64_t { return x; }

template <typename T> struct Rational {
  [[no_unique_address]] T numerator{0};
  [[no_unique_address]] T denominator{1};

  constexpr Rational<T>() = default;
  constexpr Rational<T>(T coef) : numerator(coef){};
  constexpr Rational<T>(T n, T d)
    : numerator(d > 0 ? n : -n), denominator(n ? (d > 0 ? d : -d) : 1) {}
  constexpr static auto create(T n, T d) -> Rational<T> {
    if (n) {
      T sign = 2 * (d > 0) - 1;
      T g = gcd(n, d);
      n *= sign;
      d *= sign;
      if (g != 1) {
        n /= g;
        d /= g;
      }
      return Rational<T>{n, d};
    } else {
      return Rational<T>{0, 1};
    }
  }
  constexpr static auto createPositiveDenominator(T n, T d) -> Rational<T> {
    if (n) {
      T g = gcd(n, d);
      if (g != 1) {
        n /= g;
        d /= g;
      }
      return Rational<T>{n, d};
    } else {
      return Rational<T>{0, 1};
    }
  }

  [[nodiscard]] constexpr auto safeAdd(Rational<T> y) const
    -> std::optional<Rational<T>> {
    auto [xd, yd] = divgcd(denominator, y.denominator);
    T a, b, n, d;
    bool o1 = __builtin_mul_overflow(numerator, yd, &a);
    bool o2 = __builtin_mul_overflow(y.numerator, xd, &b);
    bool o3 = __builtin_mul_overflow(denominator, yd, &d);
    bool o4 = __builtin_add_overflow(a, b, &n);
    if ((o1 | o2) | (o3 | o4)) {
      return {};
    } else if (n) {
      auto [nn, nd] = divgcd(n, d);
      return Rational<T>{nn, nd};
    } else {
      return Rational<T>{0, 1};
    }
  }
  constexpr auto operator+(Rational<T> y) const -> Rational<T> {
    return *safeAdd(y); // NOLINT(bugprone-unchecked-optional-access)
  }
  constexpr auto operator+=(Rational<T> y) -> Rational<T> & {
    std::optional<Rational<T>> a = *this + y;
    assert(a.has_value());
    *this = *a;
    return *this;
  }
  [[nodiscard]] constexpr auto safeSub(Rational<T> y) const
    -> std::optional<Rational<T>> {
    auto [xd, yd] = divgcd(denominator, y.denominator);
    T a, b, n, d;
    bool o1 = __builtin_mul_overflow(numerator, yd, &a);
    bool o2 = __builtin_mul_overflow(y.numerator, xd, &b);
    bool o3 = __builtin_mul_overflow(denominator, yd, &d);
    bool o4 = __builtin_sub_overflow(a, b, &n);
    if ((o1 | o2) | (o3 | o4)) {
      return {};
    } else if (n) {
      auto [nn, nd] = divgcd(n, d);
      return Rational<T>{nn, nd};
    } else {
      return Rational<T>{0, 1};
    }
  }
  constexpr auto operator-(Rational<T> y) const -> Rational<T> {
    return *safeSub(y); // NOLINT(bugprone-unchecked-optional-access)
  }
  constexpr auto operator-=(Rational<T> y) -> Rational<T> & {
    std::optional<Rational<T>> a = *this - y;
    assert(a.has_value());
    *this = *a;
    return *this;
  }
  [[nodiscard]] constexpr auto safeMul(T y) const
    -> std::optional<Rational<T>> {
    auto [xd, yn] = divgcd(denominator, y);
    T n;
    if (__builtin_mul_overflow(numerator, yn, &n)) return {};
    else return Rational<T>{n, xd};
  }
  [[nodiscard]] constexpr auto safeMul(Rational<T> y) const
    -> std::optional<Rational<T>> {
    if ((numerator != 0) & (y.numerator != 0)) {
      auto [xn, yd] = divgcd(numerator, y.denominator);
      auto [xd, yn] = divgcd(denominator, y.numerator);
      T n, d;
      bool o1 = __builtin_mul_overflow(xn, yn, &n);
      bool o2 = __builtin_mul_overflow(xd, yd, &d);
      if (o1 | o2) return {};
      else return Rational<T>{n, d};
    } else {
      return Rational<T>{0, 1};
    }
  }
  constexpr auto operator*(T y) const -> Rational<T> {
    return *safeMul(y); // NOLINT(bugprone-unchecked-optional-access)
  }
  constexpr auto operator*(Rational<T> y) const -> Rational<T> {
    return *safeMul(y); // NOLINT(bugprone-unchecked-optional-access)
  }
  constexpr auto operator*=(Rational<T> y) -> Rational<T> & {
    if ((numerator != 0) & (y.numerator != 0)) {
      auto [xn, yd] = divgcd(numerator, y.denominator);
      auto [xd, yn] = divgcd(denominator, y.numerator);
      numerator = xn * yn;
      denominator = xd * yd;
    } else {
      numerator = 0;
      denominator = 1;
    }
    return *this;
  }
  [[nodiscard]] constexpr auto inv() const -> Rational<T> {
    if (numerator < 0) {
      // make sure we don't have overflow
      assert(denominator != std::numeric_limits<int64_t>::min());
      return Rational<T>{-denominator, -numerator};
    } else {
      return Rational<T>{denominator, numerator};
    }
    // return Rational<T>{denominator, numerator};
    // bool positive = numerator > 0;
    // return Rational<T>{positive ? denominator : -denominator,
    //                 positive ? numerator : -numerator};
  }
  [[nodiscard]] constexpr auto safeDiv(Rational<T> y) const
    -> std::optional<Rational<T>> {
    return (*this) * y.inv();
  }
  constexpr auto operator/(Rational<T> y) const -> Rational<T> {
    return *safeDiv(y); // NOLINT(bugprone-unchecked-optional-access)
  }
  // *this -= a*b
  constexpr auto fnmadd(Rational<T> a, Rational<T> b) -> bool {
    if (std::optional<Rational<T>> ab = a.safeMul(b)) {
      if (std::optional<Rational<T>> c = safeSub(*ab)) {
        *this = *c;
        return false;
      }
    }
    return true;
  }
  constexpr auto div(Rational<T> a) -> bool {
    if (std::optional<Rational<T>> d = safeDiv(a)) {
      *this = *d;
      return false;
    }
    return true;
  }
  // Rational<T> operator/=(Rational<T> y) { return (*this) *= y.inv(); }
  constexpr operator double() {
    return double(numerator) / double(denominator);
  }

  constexpr auto operator==(Rational<T> y) const -> bool {
    return (numerator == y.numerator) & (denominator == y.denominator);
  }
  constexpr auto operator!=(Rational<T> y) const -> bool {
    return (numerator != y.numerator) | (denominator != y.denominator);
  }
  [[nodiscard]] constexpr auto isEqual(int64_t y) const -> bool {
    if (denominator == 1) return (numerator == y);
    else if (denominator == -1) return (numerator == -y);
    else return false;
  }
  constexpr auto operator==(int y) const -> bool { return isEqual(y); }
  constexpr auto operator==(int64_t y) const -> bool { return isEqual(y); }
  constexpr auto operator!=(int y) const -> bool { return !isEqual(y); }
  constexpr auto operator!=(int64_t y) const -> bool { return !isEqual(y); }
  constexpr auto operator<(Rational<T> y) const -> bool {
    return (widen(numerator) * widen(y.denominator)) <
           (widen(y.numerator) * widen(denominator));
  }
  constexpr auto operator<=(Rational<T> y) const -> bool {
    return (widen(numerator) * widen(y.denominator)) <=
           (widen(y.numerator) * widen(denominator));
  }
  constexpr auto operator>(Rational<T> y) const -> bool {
    return (widen(numerator) * widen(y.denominator)) >
           (widen(y.numerator) * widen(denominator));
  }
  constexpr auto operator>=(Rational<T> y) const -> bool {
    return (widen(numerator) * widen(y.denominator)) >=
           (widen(y.numerator) * widen(denominator));
  }
  constexpr auto operator>=(int y) const -> bool {
    return *this >= Rational<T>(y);
  }
  [[nodiscard]] constexpr auto isInteger() const -> bool {
    return denominator == 1;
  }
  constexpr void negate() { numerator = -numerator; }
  constexpr operator bool() const { return numerator != 0; }

  friend inline auto operator<<(llvm::raw_ostream &os, const Rational<T> &x)
    -> llvm::raw_ostream & {
    os << x.numerator;
    if (x.denominator != 1) os << " // " << x.denominator;
    return os;
  }
  void dump() const { llvm::errs() << *this << "\n"; }
};
template <typename T>
constexpr auto gcd(Rational<T> x, Rational<T> y) -> std::optional<Rational<T>> {
  return Rational<T>{gcd(x.numerator, y.numerator),
                     lcm(x.denominator, y.denominator)};
}
inline auto denomLCM(PtrVector<Rational<int64_t>> x) -> int64_t {
  int64_t l = 1;
  for (auto r : x) l = lcm(l, r.denominator);
  return l;
}

static_assert(AbstractVector<PtrVector<Rational<int64_t>>>);
static_assert(AbstractVector<LinearAlgebra::ElementwiseVectorBinaryOp<
                LinearAlgebra::Sub, PtrVector<Rational<int64_t>>,
                PtrVector<Rational<int64_t>>>>);
