#pragma once
#include <llvm/Support/raw_ostream.h>

namespace LinearAlgebra {
enum class AxisType {
  Row,
  Column,
  RowStride,
};
inline auto operator<<(llvm::raw_ostream &os, AxisType x)
  -> llvm::raw_ostream & {
  switch (x) {
  case AxisType::Row:
    return os << "Row";
  case AxisType::Column:
    return os << "Column";
  case AxisType::RowStride:
    return os << "RowStride";
  }
  llvm_unreachable("Unknown AxisType");
  return os;
}

// strong typing
template <AxisType T> struct AxisInt {
  using V = size_t;
  [[no_unique_address]] V value{0};
  // [[no_unique_address]] unsigned int value{0};
  constexpr AxisInt() = default;
  constexpr AxisInt(V v) : value(v) {}
  explicit constexpr operator size_t() const { return value; }
  explicit constexpr operator ptrdiff_t() const { return value; }
  explicit constexpr operator unsigned() const { return value; }
  explicit constexpr operator bool() const { return value; }

  constexpr auto operator+(V i) const -> AxisInt<T> { return value + i; }
  constexpr auto operator-(V i) const -> AxisInt<T> { return value - i; }
  constexpr auto operator*(V i) const -> AxisInt<T> { return value * i; }
  constexpr auto operator/(V i) const -> AxisInt<T> { return value / i; }
  constexpr auto operator%(V i) const -> AxisInt<T> { return value % i; }
  constexpr auto operator==(V i) const -> bool { return value == i; }
  constexpr auto operator!=(V i) const -> bool { return value != i; }
  constexpr auto operator<(V i) const -> bool { return value < i; }
  constexpr auto operator<=(V i) const -> bool { return value <= i; }
  constexpr auto operator>(V i) const -> bool { return value > i; }
  constexpr auto operator>=(V i) const -> bool { return value >= i; }
  constexpr auto operator++() -> AxisInt<T> & {
    ++value;
    return *this;
  }
  constexpr auto operator++(int) -> AxisInt<T> { return value++; }
  constexpr auto operator--() -> AxisInt<T> & {
    --value;
    return *this;
  }
  constexpr auto operator--(int) -> AxisInt<T> { return value--; }
  constexpr auto operator+=(AxisInt<T> i) -> AxisInt<T> & {
    value += V(i);
    return *this;
  }
  constexpr auto operator+=(V i) -> AxisInt<T> & {
    value += i;
    return *this;
  }
  constexpr auto operator-=(AxisInt<T> i) -> AxisInt<T> & {
    value -= V(i);
    return *this;
  }
  constexpr auto operator-=(V i) -> AxisInt<T> & {
    value -= i;
    return *this;
  }
  constexpr auto operator*=(AxisInt<T> i) -> AxisInt<T> & {
    value *= V(i);
    return *this;
  }
  constexpr auto operator*=(V i) -> AxisInt<T> & {
    value *= i;
    return *this;
  }
  constexpr auto operator/=(AxisInt<T> i) -> AxisInt<T> & {
    value /= V(i);
    return *this;
  }
  constexpr auto operator/=(V i) -> AxisInt<T> & {
    value /= i;
    return *this;
  }
  constexpr auto operator%=(AxisInt<T> i) -> AxisInt<T> & {
    value %= V(i);
    return *this;
  }
  constexpr auto operator%=(V i) -> AxisInt<T> & {
    value %= i;
    return *this;
  }
  constexpr auto operator*() const -> V { return value; }

  friend inline auto operator<<(llvm::raw_ostream &os, AxisInt<T> x)
    -> llvm::raw_ostream & {
    return os << T << "{" << *x << "}";
  }
};
template <typename T, AxisType W>
constexpr auto operator+(T *p, AxisInt<W> y) -> T * {
  return p + *y;
}
template <typename T, AxisType W>
constexpr auto operator-(T *p, AxisInt<W> y) -> T * {
  return p - *y;
}

template <AxisType T>
constexpr auto operator+(AxisInt<T> x, AxisInt<T> y) -> AxisInt<T> {
  return (*x) + (*y);
}
template <AxisType T>
constexpr auto operator-(AxisInt<T> x, AxisInt<T> y) -> AxisInt<T> {
  return (*x) - (*y);
}
template <AxisType T>
constexpr auto operator*(AxisInt<T> x, AxisInt<T> y) -> AxisInt<T> {
  return (*x) * (*y);
}
template <AxisType T>
constexpr auto operator/(AxisInt<T> x, AxisInt<T> y) -> AxisInt<T> {
  return (*x) / (*y);
}
template <AxisType T>
constexpr auto operator%(AxisInt<T> x, AxisInt<T> y) -> AxisInt<T> {
  return (*x) % (*y);
}
template <AxisType S, AxisType T>
constexpr auto operator==(AxisInt<S> x, AxisInt<T> y) -> bool {
  return *x == *y;
}
template <AxisType S, AxisType T>
constexpr auto operator!=(AxisInt<S> x, AxisInt<T> y) -> bool {
  return *x != *y;
}
template <AxisType T>
constexpr auto operator<(AxisInt<T> x, AxisInt<T> y) -> bool {
  return *x < *y;
}
template <AxisType T>
constexpr auto operator<=(AxisInt<T> x, AxisInt<T> y) -> bool {
  return *x <= *y;
}
template <AxisType T>
constexpr auto operator>(AxisInt<T> x, AxisInt<T> y) -> bool {
  return *x > *y;
}
template <AxisType T>
constexpr auto operator>=(AxisInt<T> x, AxisInt<T> y) -> bool {
  return *x >= *y;
}
using Col = AxisInt<AxisType::Column>;
using Row = AxisInt<AxisType::Row>;
using RowStride = AxisInt<AxisType::RowStride>;
using CarInd = std::pair<Row, Col>;

constexpr auto operator*(RowStride x, Row y) -> size_t { return (*x) * (*y); }
constexpr auto operator>=(RowStride x, Col u) -> bool { return (*x) >= (*u); }
constexpr auto operator<(RowStride x, Col u) -> bool { return (*x) < (*u); }

static_assert(sizeof(Row) == sizeof(size_t));
static_assert(sizeof(Col) == sizeof(size_t));
static_assert(sizeof(RowStride) == sizeof(size_t));
constexpr auto operator*(Row r, Col c) -> Row::V { return *r * *c; }

constexpr auto operator<(size_t x, Row y) -> bool { return x < size_t(y); }
constexpr auto operator<(size_t x, Col y) -> bool { return x < size_t(y); }
constexpr auto operator>(size_t x, Row y) -> bool { return x > size_t(y); }
constexpr auto operator>(size_t x, Col y) -> bool { return x > size_t(y); }

constexpr auto operator+(size_t x, Col y) -> Col { return Col{x + size_t(y)}; }
constexpr auto operator-(size_t x, Col y) -> Col { return Col{x - size_t(y)}; }
constexpr auto operator*(size_t x, Col y) -> Col { return Col{x * size_t(y)}; }
constexpr auto operator+(size_t x, Row y) -> Row { return Row{x + size_t(y)}; }
constexpr auto operator-(size_t x, Row y) -> Row { return Row{x - size_t(y)}; }
constexpr auto operator*(size_t x, Row y) -> Row { return Row{x * size_t(y)}; }
constexpr auto operator+(size_t x, RowStride y) -> RowStride {
  return RowStride{x + size_t(y)};
}
constexpr auto operator-(size_t x, RowStride y) -> RowStride {
  return RowStride{x - size_t(y)};
}
constexpr auto operator*(size_t x, RowStride y) -> RowStride {
  return RowStride{x * size_t(y)};
}

constexpr auto max(Col N, RowStride X) -> RowStride {
  return RowStride{std::max(size_t(N), size_t(X))};
}
constexpr auto min(Col N, Col X) -> Col {
  return Col{std::max(Col::V(N), Col::V(X))};
}
constexpr auto min(Row N, Col X) -> size_t {
  return std::min(size_t(N), size_t(X));
}

template <typename T>
concept RowOrCol = std::same_as<T, Row> || std::same_as<T, Col>;

constexpr auto unwrapRow(Row x) -> size_t { return size_t(x); }
constexpr auto unwrapCol(Col x) -> size_t { return size_t(x); }
constexpr auto unwrapRow(auto x) { return x; }
constexpr auto unwrapCol(auto x) { return x; }

} // namespace LinearAlgebra
