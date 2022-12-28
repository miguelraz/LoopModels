#pragma once
#include <concepts>
#include <cstddef>
#include <limits>
#include <type_traits>

template <typename T, T N> struct Static {
  using type = T;
  static constexpr bool isStatic = true;
  constexpr Static() = default;
  template <typename O, O M>
  constexpr auto operator+(const Static<O, M> &) const {
    return Static<decltype(N + M), N + M>{};
  }
  template <typename O, O M>
  constexpr auto operator-(const Static<O, M> &) const {
    return Static<decltype(N - M), N - M>{};
  }
  template <typename O, O M>
  constexpr auto operator*(const Static<O, M> &) const {
    return Static<decltype(N * M), N * M>{};
  }
  template <typename O, O M>
  constexpr auto operator/(const Static<O, M> &) const {
    static_assert(M != 0, "Division by zero");
    return Static<decltype(N / M), N / M>{};
  }
  template <typename O, O M>
  constexpr auto operator%(const Static<O, M> &) const {
    static_assert(M != 0, "Remainder in division by zero");
    return Static<decltype(N % M), N % M>{};
  }
  template <typename O, O M>
  constexpr auto operator&(const Static<O, M> &) const {
    return Static<decltype(N & M), N & M>{};
  }
  template <typename O, O M>
  constexpr auto operator|(const Static<O, M> &) const {
    return Static<decltype(N | M), N | M>{};
  }
  template <typename O, O M>
  constexpr auto operator^(const Static<O, M> &) const {
    return Static<decltype(N ^ M), N ^ M>{};
  }
  template <typename O, O M>
  constexpr auto operator<<(const Static<O, M> &) const {
    return Static<decltype(N << M), (N << M)>{};
  }
  template <typename O, O M>
  constexpr auto operator>>(const Static<O, M> &) const {
    return Static<decltype(N >> M), (N >> M)>{};
  }
  template <typename O, O M>
  constexpr auto operator==(const Static<O, M> &) const {
    return Static<decltype(N == M), N == M>{};
  }
  template <typename O, O M>
  constexpr auto operator!=(const Static<O, M> &) const {
    return Static<decltype(N != M), N != M>{};
  }
  template <typename O, O M>
  constexpr auto operator<(const Static<O, M> &) const {
    return Static<decltype(N < M), (N < M)>{};
  }
  template <typename O, O M>
  constexpr auto operator>(const Static<O, M> &) const {
    return Static<decltype(N > M), (N > M)>{};
  }
  template <typename O, O M>
  constexpr auto operator<=(const Static<O, M> &) const {
    return Static<decltype(N <= M), (N <= M)>{};
  }
  template <typename O, O M>
  constexpr auto operator>=(const Static<O, M> &) const {
    return Static<decltype(N >= M), (N >= M)>{};
  }
  constexpr operator T() const { return N; }
};
template <typename T> using static_type_t = typename T::type;

template <typename T>
concept SignedStaticInt = T::isStatic && std::is_signed_v<static_type_t<T>>;
template <typename T>
concept UnsignedStaticInt = T::isStatic && std::is_unsigned_v<static_type_t<T>>;

static_assert(std::is_convertible_v<Static<size_t, 16>, size_t>);
static_assert(!std::is_convertible_v<size_t, Static<size_t, 16>>);

static_assert(SignedStaticInt<Static<int, 0>>);
static_assert(!UnsignedStaticInt<Static<int, 0>>);
static_assert(!SignedStaticInt<Static<unsigned int, 0>>);
static_assert(UnsignedStaticInt<Static<unsigned int, 0>>);

static_assert(std::same_as<int, decltype(Static<int, 0>{} + 0)>);
static_assert(
  std::same_as<ptrdiff_t, decltype(Static<int, 0>{} + ptrdiff_t(0))>);

template <typename T>
concept IntConvertible =
  std::is_integral_v<T> || std::is_convertible_v<T, size_t>;
template <typename T>
concept SignedIntConvertible = std::is_signed_v<T> ||
                               (std::is_convertible_v<T, ptrdiff_t> &&
                                SignedStaticInt<T>);
;
template <typename T>
concept UnsignedIntConvertible = std::is_unsigned_v<T> ||
                                 (std::is_convertible_v<T, size_t> &&
                                  UnsignedStaticInt<T>);
;

static_assert(IntConvertible<size_t>);
static_assert(IntConvertible<ptrdiff_t>);
static_assert(IntConvertible<Static<size_t, 0>>);
static_assert(IntConvertible<Static<ptrdiff_t, 0>>);

static_assert(!SignedIntConvertible<size_t>);
static_assert(SignedIntConvertible<ptrdiff_t>);
static_assert(!SignedIntConvertible<Static<size_t, 0>>);
static_assert(SignedIntConvertible<Static<ptrdiff_t, 0>>);

static_assert(UnsignedIntConvertible<size_t>);
static_assert(!UnsignedIntConvertible<ptrdiff_t>);
static_assert(UnsignedIntConvertible<Static<size_t, 0>>);
static_assert(!UnsignedIntConvertible<Static<ptrdiff_t, 0>>);

static_assert(
  std::same_as<Static<ptrdiff_t, 12>,
               decltype(Static<ptrdiff_t, 5>{} + Static<ptrdiff_t, 7>{})>);

template <typename T> struct SignedType {};
template <SignedIntConvertible T> struct SignedType<T> {
  using type = T;
};
template <UnsignedStaticInt T> struct SignedType<T> {
  static_assert(
    T::value <=
    std::numeric_limits<std::make_signed_t<typename T::type>>::max());
  using type = Static<std::make_signed_t<typename T::type>, T::value>;
};
template <std::unsigned_integral T> struct SignedType<T> {
  using type = std::make_signed_t<T>;
};
template <typename T> struct UnsignedType {};
template <UnsignedIntConvertible T> struct UnsignedType<T> {
  using type = T;
};
template <SignedStaticInt T> struct UnsignedType<T> {
  static_assert(T::value >= 0);
  using type = Static<std::make_unsigned_t<typename T::type>, T::value>;
};
template <std::signed_integral T> struct UnsignedType<T> {
  using type = std::make_unsigned_t<T>;
};
template <typename T> using signed_type_t = typename SignedType<T>::type;
template <typename T> using unsigned_type_t = typename UnsignedType<T>::type;

template <typename T, T U> constexpr auto signedOf(Static<T, U>) {
  using I = std::make_signed_t<T>;
  return Static<I, I(U)>{};
}
template <typename T, T U> constexpr auto unsignedOf(Static<T, U>) {
  using I = std::make_unsigned_t<T>;
  return Static<I, I(U)>{};
}
template <std::integral T> constexpr auto signedOf(T i) {
  using I = std::make_signed_t<T>;
  return I(i);
}
template <std::integral T> constexpr auto unsignedOf(T i) {
  using I = std::make_unsigned_t<T>;
  return I(i);
}
