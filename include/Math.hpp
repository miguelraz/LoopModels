#pragma once
// We'll follow Julia style, so anything that's not a constructor, destructor,
// nor an operator will be outside of the struct/class.

#include "./TypePromotion.hpp"
#include "BitSets.hpp"
#include "StaticInts.hpp"
#include <algorithm>
#include <bit>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <llvm/ADT/ArrayRef.h>
#include <llvm/ADT/Optional.h>
#include <llvm/ADT/SmallVector.h>
#include <llvm/Support/raw_ostream.h>
#include <numeric>
#include <optional>
#include <sched.h>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
// #ifndef NDEBUG
// #include <memory>
// #include <stacktrace>
// using stacktrace =
//     std::basic_stacktrace<std::allocator<std::stacktrace_entry>>;
// #endif

struct Rational;
namespace LinearAlgebra {

inline auto isZero(auto x) -> bool { return x == 0; }

inline auto allZero(const auto &x) -> bool {
  for (auto &a : x)
    if (!isZero(a))
      return false;
  return true;
}
inline auto allGEZero(const auto &x) -> bool {
  for (auto &a : x)
    if (a < 0)
      return false;
  return true;
}
inline auto allLEZero(const auto &x) -> bool {
  for (auto &a : x)
    if (a > 0)
      return false;
  return true;
}

inline auto countNonZero(const auto &x) -> size_t {
  size_t i = 0;
  for (auto &a : x)
    i += (a != 0);
  return i;
}

template <typename T>
concept AbstractVector =
  HasEltype<T> && requires(T t, size_t i) {
                    { t(i) } -> std::convertible_to<eltype_t<T>>;
                    { t.size() } -> std::convertible_to<size_t>;
                    { t.view() };
                    // {
                    //     std::remove_reference_t<T>::canResize
                    //     } -> std::same_as<const bool &>;
                    // {t.extendOrAssertSize(i)};
                  };

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
template <AxisType T, IntConvertible V = size_t> struct AxisInt {
  static constexpr AxisType Axis = T;
  using I = V;
  [[no_unique_address]] V value;
  // [[no_unique_address]] unsigned int value{0};
  constexpr AxisInt() = default;
  constexpr AxisInt(V value) : value(value) {}
  explicit constexpr operator size_t() const { return value; }
  explicit constexpr operator ptrdiff_t() const { return value; }
  explicit constexpr operator unsigned() const { return value; }
  explicit constexpr operator bool() const { return value; }

  constexpr auto operator+(V i) const -> AxisInt<T, decltype(value + i)> {
    return value + i;
  }
  constexpr auto operator-(V i) const -> AxisInt<T, decltype(value - i)> {
    return value - i;
  }
  constexpr auto operator*(V i) const -> AxisInt<T, decltype(value * i)> {
    return value * i;
  }
  constexpr auto operator/(V i) const {
    auto d = value / i;
    return AxisInt<T, decltype(d)>(d);
  }
  constexpr auto operator%(V i) const {
    auto r = value % i;
    return AxisInt<T, decltype(r)>(r);
  }
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
  template <IntConvertible U> constexpr operator AxisInt<T, U>() const {
    static_assert(std::is_convertible_v<V, U>);
    return AxisInt<T, U>(U(value));
  }
};
template <AxisType T>
inline auto operator<<(llvm::raw_ostream &os, AxisInt<T> x)
  -> llvm::raw_ostream & {
  return os << T << "{" << *x << "}";
}

template <typename T, AxisType W>
constexpr auto operator+(T *p, AxisInt<W> y) -> T * {
  return p + *y;
}
template <typename T, AxisType W>
constexpr auto operator-(T *p, AxisInt<W> y) -> T * {
  return p - *y;
}

template <AxisType T, IntConvertible U, IntConvertible V>
constexpr auto operator+(AxisInt<T, U> x, AxisInt<T, V> y)
  -> AxisInt<T, decltype(x.value + y.value)> {
  return (*x) + (*y);
}
template <AxisType T, IntConvertible U, IntConvertible V>
constexpr auto operator-(AxisInt<T, U> x, AxisInt<T, V> y)
  -> AxisInt<T, decltype(x.value - y.value)> {
  return (*x) - (*y);
}
template <AxisType T, IntConvertible U, IntConvertible V>
constexpr auto operator*(AxisInt<T, U> x, AxisInt<T, V> y)
  -> AxisInt<T, decltype(x.value * y.value)> {
  return (*x) * (*y);
}
template <AxisType T, IntConvertible U, IntConvertible V>
constexpr auto operator/(AxisInt<T, U> x, AxisInt<T, V> y)
  -> AxisInt<T, decltype(x.value / y.value)> {
  return (*x) / (*y);
}
template <AxisType T, IntConvertible U, IntConvertible V>
constexpr auto operator%(AxisInt<T, U> x, AxisInt<T, V> y)
  -> AxisInt<T, decltype(x.value % y.value)> {
  return (*x) % (*y);
}
template <AxisType T>
constexpr auto operator==(AxisInt<T> x, AxisInt<T> y) -> bool {
  return *x == *y;
}
template <AxisType T>
constexpr auto operator!=(AxisInt<T> x, AxisInt<T> y) -> bool {
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
template <IntConvertible I = size_t> using Row = AxisInt<AxisType::Row, I>;
template <IntConvertible I = size_t> using Col = AxisInt<AxisType::Column, I>;
template <IntConvertible I = size_t>
using RowStride = AxisInt<AxisType::RowStride, I>;

static_assert(std::is_convertible_v<Row<Static<size_t, 16>>, Row<>>);
// static_assert(!std::is_convertible_v<Row<>, Row<Static<size_t, 16>>>);
static_assert(std::is_convertible_v<Col<Static<size_t, 16>>, Col<>>);
static_assert(!std::is_convertible_v<Row<>, Col<>>);

template <AxisType T, IntConvertible I>
constexpr auto axis(I i) -> AxisInt<T, I> {
  return AxisInt<T, I>(i);
}
template <std::signed_integral I> constexpr auto toCol(I i) -> Col<ptrdiff_t> {
  return i;
}
template <std::signed_integral I> constexpr auto toRow(I i) -> Row<ptrdiff_t> {
  return i;
}
template <std::signed_integral I>
constexpr auto toRowStride(I i) -> RowStride<ptrdiff_t> {
  return i;
}
template <std::unsigned_integral I> constexpr auto toCol(I i) -> Col<size_t> {
  return i;
}
template <std::unsigned_integral I> constexpr auto toRow(I i) -> Row<size_t> {
  return i;
}
template <std::unsigned_integral I>
constexpr auto toRowStride(I i) -> RowStride<size_t> {
  return i;
}
template <IntConvertible I, I N>
constexpr auto toCol(Static<I, N> i) -> Col<Static<I, N>> {
  return i;
}
template <IntConvertible I, I N>
constexpr auto toRow(Static<I, N> i) -> Row<Static<I, N>> {
  return i;
}
template <IntConvertible I, I N>
constexpr auto toRowStride(Static<I, N> i) -> RowStride<Static<I, N>> {
  return i;
}

template <IntConvertible I, IntConvertible J>
constexpr auto operator*(RowStride<I> x, Row<J> y) {
  return (*x) * (*y);
}
template <IntConvertible I, IntConvertible J>
constexpr auto operator>=(RowStride<I> x, Col<J> u) -> bool {
  return (*x) >= (*u);
}
template <IntConvertible I, IntConvertible J>
constexpr auto operator<(RowStride<I> x, Col<J> u) -> bool {
  return (*x) < (*u);
}

static_assert(sizeof(Row<>) == sizeof(size_t));
static_assert(sizeof(Col<>) == sizeof(size_t));
static_assert(sizeof(RowStride<>) == sizeof(size_t));
template <IntConvertible I, IntConvertible J>
constexpr auto operator*(Row<I> r, Col<J> c) {
  return *r * *c;
}

template <IntConvertible I>
constexpr auto operator<(size_t x, Row<I> y) -> bool {
  return x < size_t(y);
}
template <IntConvertible I>
constexpr auto operator<(size_t x, Col<I> y) -> bool {
  return x < size_t(y);
}
template <IntConvertible I>
constexpr auto operator>(size_t x, Row<I> y) -> bool {
  return x > size_t(y);
}
template <IntConvertible I>
constexpr auto operator>(size_t x, Col<I> y) -> bool {
  return x > size_t(y);
}

template <AxisType T, IntConvertible I>
constexpr auto operator+(size_t x, AxisInt<T, I> y) {
  return axis<T>(x + (*y));
}
template <AxisType T, IntConvertible I>
constexpr auto operator-(size_t x, AxisInt<T, I> y) {
  return axis<T>(x - (*y));
}
template <AxisType T, IntConvertible I>
constexpr auto operator*(size_t x, AxisInt<T, I> y) {
  return axis<T>(x * (*y));
}

template <IntConvertible I, IntConvertible J>
constexpr auto max(Col<I> N, RowStride<J> X) {
  return axis<AxisType::RowStride>(std::max(size_t(N), size_t(X)));
}
template <IntConvertible I, IntConvertible J>
constexpr auto min(Col<I> N, Col<J> X) {
  return axis<AxisType::Column>(std::max(*N, *X));
}
template <typename T>
concept IsRow = std::remove_reference_t<T>::Axis ==
AxisType::Row;
template <typename T>
concept IsCol = std::remove_reference_t<T>::Axis ==
AxisType::Column;
template <typename T>
concept IsRowStride = std::remove_reference_t<T>::Axis ==
AxisType::RowStride;

template <typename T>
concept RowOrCol = IsRow<T> || IsCol<T>;
static_assert(RowOrCol<Row<>>);
static_assert(RowOrCol<Col<>>);
static_assert(!RowOrCol<RowStride<>>);
static_assert(!RowOrCol<int>);

template <typename T>
concept AbstractMatrixCore =
  HasEltype<T> && requires(T t, size_t i) {
                    { t(i, i) } -> std::convertible_to<eltype_t<T>>;
                    { t.numRow() } -> IsRow;
                    { t.numCol() } -> IsCol;
                    { t.size().first } -> IsRow;
                    { t.size().second } -> IsCol;
                    // {
                    //     std::remove_reference_t<T>::canResize
                    //     } -> std::same_as<const bool &>;
                    // {t.extendOrAssertSize(i, i)};
                  };
template <typename T>
concept AbstractMatrix = AbstractMatrixCore<T> && requires(T t, size_t i) {
                                                    {
                                                      t.view()
                                                      } -> AbstractMatrixCore;
                                                  };
template <typename T>
concept AbstractRowMajorMatrix = AbstractMatrix<T> && requires(T t) {
                                                        {
                                                          t.rowStride()
                                                          } -> IsRowStride;
                                                      };

inline auto copyto(AbstractVector auto &y, const AbstractVector auto &x)
  -> auto & {
  const size_t M = x.size();
  y.extendOrAssertSize(M);
  for (size_t i = 0; i < M; ++i)
    y(i) = x(i);
  return y;
}
inline auto copyto(AbstractMatrixCore auto &A, const AbstractMatrixCore auto &B)
  -> auto & {
  const IsRow auto M = B.numRow();
  const IsCol auto N = B.numCol();
  A.extendOrAssertSize(M, N);
  for (size_t r = 0; r < M; ++r)
    for (size_t c = 0; c < N; ++c)
      A(r, c) = B(r, c);
  return A;
}

auto operator==(const AbstractMatrix auto &A, const AbstractMatrix auto &B)
  -> bool {
  const IsRow auto M = B.numRow();
  const IsCol auto N = B.numCol();
  if ((M != A.numRow()) || (N != A.numCol()))
    return false;
  for (size_t r = 0; r < M; ++r)
    for (size_t c = 0; c < N; ++c)
      if (A(r, c) != B(r, c))
        return false;
  return true;
}

struct Add {
  constexpr auto operator()(auto x, auto y) const { return x + y; }
};
struct Sub {
  constexpr auto operator()(auto x) const { return -x; }
  constexpr auto operator()(auto x, auto y) const { return x - y; }
};
struct Mul {
  constexpr auto operator()(auto x, auto y) const { return x * y; }
};
struct Div {
  constexpr auto operator()(auto x, auto y) const { return x / y; }
};

template <typename Op, typename A> struct ElementwiseUnaryOp {
  using eltype = typename A::eltype;
  [[no_unique_address]] const Op op;
  [[no_unique_address]] const A a;
  auto operator()(size_t i) const { return op(a(i)); }
  auto operator()(size_t i, size_t j) const { return op(a(i, j)); }

  [[nodiscard]] constexpr auto size() const { return a.size(); }
  [[nodiscard]] constexpr auto numRow() const { return a.numRow(); }
  [[nodiscard]] constexpr auto numCol() const { return a.numCol(); }
  [[nodiscard]] constexpr auto view() const { return *this; };
};
// scalars broadcast
constexpr auto get(const std::integral auto A, size_t) { return A; }
constexpr auto get(const std::floating_point auto A, size_t) { return A; }
constexpr auto get(const std::integral auto A, size_t, size_t) { return A; }
constexpr auto get(const std::floating_point auto A, size_t, size_t) {
  return A;
}
inline auto get(const AbstractVector auto &A, size_t i) { return A(i); }
inline auto get(const AbstractMatrix auto &A, size_t i, size_t j) {
  return A(i, j);
}

constexpr auto size(const std::integral auto) -> size_t { return 1; }
constexpr auto size(const std::floating_point auto) -> size_t { return 1; }
constexpr auto size(const AbstractVector auto &x) -> size_t { return x.size(); }

template <typename T>
concept Scalar =
  std::integral<T> || std::floating_point<T> || std::same_as<T, Rational>;

template <typename T>
concept VectorOrScalar = AbstractVector<T> || Scalar<T>;
template <typename T>
concept MatrixOrScalar = AbstractMatrix<T> || Scalar<T>;

template <typename Op, VectorOrScalar A, VectorOrScalar B>
struct ElementwiseVectorBinaryOp {
  using eltype = promote_eltype_t<A, B>;
  [[no_unique_address]] Op op;
  [[no_unique_address]] A a;
  [[no_unique_address]] B b;
  ElementwiseVectorBinaryOp(Op op, A a, B b) : op(op), a(a), b(b) {}
  auto operator()(size_t i) const { return op(get(a, i), get(b, i)); }
  [[nodiscard]] constexpr auto size() const -> size_t {
    if constexpr (AbstractVector<A> && AbstractVector<B>) {
      const size_t N = a.size();
      assert(N == b.size());
      return N;
    } else if constexpr (AbstractVector<A>) {
      return a.size();
    } else { // if constexpr (AbstractVector<B>) {
      return b.size();
    }
  }
  [[nodiscard]] constexpr auto view() const -> auto & { return *this; };
};

template <typename Op, MatrixOrScalar A, MatrixOrScalar B>
struct ElementwiseMatrixBinaryOp {
  using eltype = promote_eltype_t<A, B>;
  [[no_unique_address]] Op op;
  [[no_unique_address]] A a;
  [[no_unique_address]] B b;
  ElementwiseMatrixBinaryOp(Op op, A a, B b) : op(op), a(a), b(b) {}
  auto operator()(size_t i, size_t j) const {
    return op(get(a, i, j), get(b, i, j));
  }
  [[nodiscard]] constexpr auto numRow() const {
    static_assert(AbstractMatrix<A> || std::integral<A> ||
                    std::floating_point<A>,
                  "Argument A to elementwise binary op is not a matrix.");
    static_assert(AbstractMatrix<B> || std::integral<B> ||
                    std::floating_point<B>,
                  "Argument B to elementwise binary op is not a matrix.");
    if constexpr (AbstractMatrix<A> && AbstractMatrix<B>) {
      const IsRow auto N = a.numRow();
      assert(N == b.numRow());
      return N;
    } else if constexpr (AbstractMatrix<A>) {
      return a.numRow();
    } else if constexpr (AbstractMatrix<B>) {
      return b.numRow();
    }
  }
  [[nodiscard]] constexpr auto numCol() const {
    static_assert(AbstractMatrix<A> || std::integral<A> ||
                    std::floating_point<A>,
                  "Argument A to elementwise binary op is not a matrix.");
    static_assert(AbstractMatrix<B> || std::integral<B> ||
                    std::floating_point<B>,
                  "Argument B to elementwise binary op is not a matrix.");
    if constexpr (AbstractMatrix<A> && AbstractMatrix<B>) {
      const IsCol auto N = a.numCol();
      assert(N == b.numCol());
      return N;
    } else if constexpr (AbstractMatrix<A>) {
      return a.numCol();
    } else if constexpr (AbstractMatrix<B>) {
      return b.numCol();
    }
  }
  [[nodiscard]] constexpr auto size() const {
    return std::make_pair(numRow(), numCol());
  }
  [[nodiscard]] constexpr auto view() const -> auto & { return *this; };
};

template <typename A> struct Transpose {
  static_assert(AbstractMatrix<A>, "Argument to transpose is not a matrix.");
  static_assert(std::is_trivially_copyable_v<A>,
                "Argument to transpose is not trivially copyable.");

  using eltype = eltype_t<A>;
  [[no_unique_address]] A a;
  auto operator()(size_t i, size_t j) const { return a(j, i); }
  [[nodiscard]] constexpr auto numRow() const {
    return axis<AxisType::Row>(*a.numCol());
  }
  [[nodiscard]] constexpr auto numCol() const {
    return axis<AxisType::Column>(*a.numCol());
  }
  [[nodiscard]] constexpr auto view() const -> auto & { return *this; };
  [[nodiscard]] constexpr auto size() const {
    return std::make_pair(numRow(), numCol());
  }
};
template <AbstractMatrix A, AbstractMatrix B> struct MatMatMul {
  using eltype = promote_eltype_t<A, B>;
  [[no_unique_address]] A a;
  [[no_unique_address]] B b;
  auto operator()(size_t i, size_t j) const {
    static_assert(AbstractMatrix<B>, "B should be an AbstractMatrix");
    auto s = (a(i, 0) * b(0, j)) * 0;
    for (size_t k = 0; k < size_t(a.numCol()); ++k)
      s += a(i, k) * b(k, j);
    return s;
  }
  [[nodiscard]] constexpr auto numRow() const { return a.numRow(); }
  [[nodiscard]] constexpr auto numCol() const { return b.numCol(); }
  [[nodiscard]] constexpr auto size() const {
    return std::make_pair(numRow(), numCol());
  }
  [[nodiscard]] constexpr auto view() const { return *this; };
};
template <AbstractMatrix A, AbstractVector B> struct MatVecMul {
  using eltype = promote_eltype_t<A, B>;
  [[no_unique_address]] A a;
  [[no_unique_address]] B b;
  auto operator()(size_t i) const {
    static_assert(AbstractVector<B>, "B should be an AbstractVector");
    auto s = (a(i, 0) * b(0)) * 0;
    for (size_t k = 0; k < a.numCol(); ++k)
      s += a(i, k) * b(k);
    return s;
  }
  [[nodiscard]] constexpr auto size() const -> size_t {
    return size_t(a.numRow());
  }
  constexpr auto view() const { return *this; };
};

static inline constexpr struct Begin {
  friend inline auto operator<<(llvm::raw_ostream &os, Begin)
    -> llvm::raw_ostream & {
    return os << 0;
  }
} begin;
static inline constexpr struct End {
  friend inline auto operator<<(llvm::raw_ostream &os, End)
    -> llvm::raw_ostream & {
    return os << "end";
  }
} end;

// FIXME: we currently lose strong typing of Row and Col when using relative
// indexing; we should preserve it, perhaps within the OffsetBegin row/struct,
// making them templated?
template <typename T>
concept ScalarValueIndex = std::integral<T> || IsRow<T> || IsCol<T>;

constexpr auto operator+(ScalarValueIndex auto x, Begin) { return x; }
constexpr auto operator+(Begin, ScalarValueIndex auto x) { return x; }
template <IntConvertible I> struct OffsetEnd {
  static constexpr bool IsOffsetEnd = true;
  [[no_unique_address]] I offset;
};
template <typename T>
concept IsOffsetEnd =
  requires(T t) { std::remove_reference_t<T>::IsOffsetEnd; };
static_assert(IsOffsetEnd<OffsetEnd<size_t>>);
static_assert(!IsOffsetEnd<size_t>);

template <IntConvertible I>
inline auto operator<<(llvm::raw_ostream &os, OffsetEnd<I> r)
  -> llvm::raw_ostream & {
  return os << "end - " << r.offset;
}
template <IntConvertible I> OffsetEnd(I) -> OffsetEnd<I>;
constexpr auto operator-(End, ScalarValueIndex auto x) {
  return OffsetEnd(unwrap(x));
}
template <IntConvertible I>
constexpr auto operator-(OffsetEnd<I> y, ScalarValueIndex auto x) {
  return OffsetEnd(y.offset + unwrap(x));
}
template <IntConvertible I>
constexpr auto operator+(OffsetEnd<I> y, ScalarValueIndex auto x) {
  return OffsetEnd(y.offset - unwrap(x));
}

template <typename T>
concept RelativeOffset =
  std::same_as<T, End> || IsOffsetEnd<T> || std::same_as<T, Begin>;

template <typename B, typename E> struct Range {
  [[no_unique_address]] B b;
  [[no_unique_address]] E e;
};
template <IntConvertible B, IntConvertible E> struct Range<B, E> {
  static constexpr bool IsRange = true;

  [[no_unique_address]] B b;
  [[no_unique_address]] E e;
  // wrapper that allows dereferencing
  struct Iterator {
    [[no_unique_address]] size_t i;
    constexpr auto operator==(E e) -> bool { return i == e; }
    auto operator++() -> Iterator & {
      ++i;
      return *this;
    }
    auto operator++(int) -> Iterator { return Iterator{i++}; }
    auto operator--() -> Iterator & {
      --i;
      return *this;
    }
    auto operator--(int) -> Iterator { return Iterator{i--}; }
    auto operator*() -> B { return i; }
  };
  [[nodiscard]] constexpr auto begin() const -> Iterator { return Iterator{b}; }
  [[nodiscard]] constexpr auto end() const -> E { return e; }
  [[nodiscard]] constexpr auto rbegin() const -> Iterator {
    return std::reverse_iterator{end()};
  }
  [[nodiscard]] constexpr auto rend() const -> E {
    return std::reverse_iterator{begin()};
  }
  [[nodiscard]] constexpr auto size() const { return e - b; }
  friend inline auto operator<<(llvm::raw_ostream &os, Range<B, E> r)
    -> llvm::raw_ostream & {
    return os << "[" << r.b << ":" << r.e << ")";
  }
  template <std::integral BB, std::integral EE>
  constexpr operator Range<BB, EE>() const {
    return Range<BB, EE>{BB(b), EE(e)};
  }
};
template <typename T>
concept IsRange = requires(T t) { std::remove_reference_t<T>::IsRange; };

template <typename T> struct StandardizeRangeBound {
  using type = T;
};
template <RowOrCol T> struct StandardizeRangeBound<T> {
  using type = size_t;
};
template <std::unsigned_integral T> struct StandardizeRangeBound<T> {
  using type = size_t;
};
template <std::signed_integral T> struct StandardizeRangeBound<T> {
  using type = ptrdiff_t;
};
template <typename T>
using StandardizeRangeBound_t = typename StandardizeRangeBound<T>::type;

constexpr auto standardizeRangeBound(auto x) { return x; }
constexpr auto standardizeRangeBound(RowOrCol auto x) { return size_t(x); }

constexpr auto standardizeRangeBound(std::unsigned_integral auto x) {
  return size_t(x);
}
constexpr auto standardizeRangeBound(std::signed_integral auto x) {
  return ptrdiff_t(x);
}

template <typename B, typename E>
Range(B b, E e) -> Range<decltype(standardizeRangeBound(b)),
                         decltype(standardizeRangeBound(e))>;

static inline constexpr struct Colon {
  [[nodiscard]] inline constexpr auto operator()(auto B, auto E) const {
    return Range{standardizeRangeBound(B), standardizeRangeBound(E)};
  }
} _; // NOLINT(bugprone-reserved-identifier)

#ifndef NDEBUG
static inline void checkIndex(size_t X, size_t x) { assert(x < X); }
inline void checkIndex(size_t X, End) { assert(X > 0); }
inline void checkIndex(size_t X, Begin) { assert(X > 0); }
inline void checkIndex(size_t X, IsOffsetEnd auto x) { assert(x.offset < X); }
template <typename B> inline void checkIndex(size_t X, Range<B, size_t> x) {
  assert(x.e <= X);
}
template <typename B, typename E> inline void checkIndex(size_t, Range<B, E>) {}
inline void checkIndex(size_t, Colon) {}
#endif

constexpr auto canonicalize(IntConvertible auto e, size_t) -> size_t {
  return e;
}
constexpr auto canonicalize(Begin, size_t) -> size_t { return 0; }
constexpr auto canonicalize(End, IntConvertible auto M) -> size_t {
  return M - 1;
}

constexpr auto canonicalize(IsOffsetEnd auto e, IntConvertible auto M)
  -> size_t {
  return M - 1 - e.offset;
}

constexpr auto canonicalizeForRange(IntConvertible auto e, size_t) { return e; }
constexpr auto canonicalizeForRange(Begin, size_t) -> Static<size_t, 0> {
  return {};
}
constexpr auto canonicalizeForRange(End, IntConvertible auto M) { return M; }
constexpr auto canonicalizeForRange(IsOffsetEnd auto e, IntConvertible auto M) {
  return M - e.offset;
}

template <typename T>
concept ScalarIndex = std::integral<T> || RelativeOffset<T>;

template <typename B, typename E, IntConvertible I>
constexpr auto canonicalizeRange(Range<B, E> r, I M) -> Range<size_t, size_t> {
  return Range<size_t, size_t>{canonicalizeForRange(r.b, M),
                               canonicalizeForRange(r.e, M)};
}
constexpr auto canonicalizeRange(Colon, size_t M) -> Range<size_t, size_t> {
  return Range<size_t, size_t>{0, M};
}

template <typename B, typename E>
constexpr auto operator+(Range<B, E> r, size_t x) {
  return _(r.b + x, r.e + x);
}
template <typename B, typename E>
constexpr auto operator-(Range<B, E> r, size_t x) {
  return _(r.b - x, r.e - x);
}

template <typename T, UnsignedIntConvertible I = size_t,
          SignedIntConvertible J = Static<ptrdiff_t, 1>>
struct PtrVector;
template <typename T, UnsignedIntConvertible I = size_t,
          SignedIntConvertible J = Static<ptrdiff_t, 1>>
struct MutPtrVector;

template <typename T, IntConvertible I> struct StridedIterator {
  T *data;
  I stride;
  constexpr auto operator*() -> T & { return *data; }
  constexpr auto operator->() -> T * { return data; }
  constexpr auto operator*() const -> const T & { return *data; }
  constexpr auto operator->() const -> const T * { return data; }
  constexpr auto operator++() -> StridedIterator & {
    data += stride;
    return *this;
  }
  constexpr auto operator++(int) -> StridedIterator {
    auto tmp = *this;
    ++*this;
    return tmp;
  }
  constexpr auto operator--() -> StridedIterator & {
    data -= stride;
    return *this;
  }
  constexpr auto operator--(int) -> StridedIterator {
    auto tmp = *this;
    --*this;
    return tmp;
  }
  constexpr auto operator+=(size_t x) -> StridedIterator & {
    data += x * stride;
    return *this;
  }
  constexpr auto operator-=(size_t x) -> StridedIterator & {
    data -= x * stride;
    return *this;
  }
  constexpr auto operator==(StridedIterator const &x) const -> bool {
    return data == x.data;
  }
  constexpr auto operator!=(StridedIterator const &x) const -> bool {
    return data != x.data;
  }
};
template <typename T>
constexpr auto default_iterator(T *data, Static<size_t, 1>) -> T * {
  return data;
}
template <typename T, IntConvertible I>
constexpr auto default_iterator(T *data, I stride) -> StridedIterator<T, I> {
  return {data, stride};
}

// CRTP for const vectors
template <typename T, typename V> struct ConstVectorCore {
  static_assert(!std::is_const_v<T>, "const T is redundant");
  using eltype = T;
  [[nodiscard]] constexpr auto rawSize() const {
    return static_cast<const V *>(this)->rawSize();
  }
  [[nodiscard]] constexpr auto size() const -> size_t {
    return static_cast<const V *>(this)->size();
  }
  [[nodiscard]] constexpr auto isEmpty() const -> bool { return size() == 0; }
  [[nodiscard]] constexpr auto data() const -> const T * {
    return static_cast<const V *>(this)->data();
  }
  [[nodiscard]] constexpr auto stride() const {
    return static_cast<const V *>(this)->stride();
  }
  constexpr auto operator[](size_t i) const -> const T & {
    assert(i < size());
    return data()[i * stride()];
  }
  auto operator==(AbstractVector auto &x) -> bool {
    if (size() != x.size())
      return false;
    const T *p = data();
    for (size_t n = 0; n < size(); ++n)
      if (p[n] != x(n))
        return false;
    return true;
  }
  [[nodiscard]] constexpr auto front() const -> const T & { return data()[0]; }
  [[nodiscard]] constexpr auto back() const -> const T & {
    return data()[(size() - 1) * stride()];
  }

  constexpr auto operator[](const ScalarIndex auto i) const -> const T & {
    size_t N = size();
#ifndef NDEBUG
    checkIndex(N, i);
#endif
    return data()[canonicalize(i, N) * stride()];
  }
  constexpr auto operator()(const ScalarIndex auto i) const -> const T & {
    size_t N = size();
#ifndef NDEBUG
    checkIndex(N, i);
#endif
    return data()[canonicalize(i, N) * stride()];
  }
  constexpr auto operator()(Range<size_t, size_t> i) const -> PtrVector<T> {
    assert(i.b <= i.e);
    assert(i.e <= size());
    return PtrVector(data() + i.b * stride(), i.e - i.b, stride());
  }
  template <typename F, typename L>
  constexpr auto operator()(Range<F, L> i) const -> PtrVector<T> {
    return (*this)(canonicalizeRange(i, size()));
  }
  [[nodiscard]] constexpr auto begin() const {
    return default_iterator(data(), stride());
  }
  [[nodiscard]] constexpr auto end() const {
    return default_iterator(data() + size(), stride());
  }
  [[nodiscard]] constexpr auto rbegin() const {
    return std::reverse_iterator(end());
  }
  [[nodiscard]] constexpr auto rend() const {
    return std::reverse_iterator(begin());
  }
  constexpr operator llvm::ArrayRef<T>() const {
    static_assert(std::is_same_v<Static<size_t, 1>, decltype(stride())>);
    return llvm::ArrayRef<T>{data(), size()};
  }
  // llvm::ArrayRef<T> arrayref() const { return llvm::ArrayRef<T>(ptr, M); }
  auto operator==(const PtrVector<T> x) const -> bool {
    return size() == x.size() && std::equal(begin(), end(), x.begin());
  }
  auto operator==(const llvm::ArrayRef<std::remove_const_t<T>> x) const
    -> bool {
    return size() == x.size() && std::equal(begin(), end(), x.begin());
  }
  [[nodiscard]] constexpr auto view() const -> PtrVector<T> { return *this; };
};

template <typename T, UnsignedIntConvertible I, SignedIntConvertible J>
struct PtrVector : ConstVectorCore<T, PtrVector<T, I>> {

  static_assert(!std::is_const_v<T>, "const T is redundant");
  using eltype = T;
  [[no_unique_address]] const T *const mem;
  [[no_unique_address]] const I N;
  [[no_unique_address]] const J X;

  [[nodiscard]] constexpr auto rawSize() const -> I { return N; }
  [[nodiscard]] constexpr auto size() const -> size_t { return N; }
  [[nodiscard]] constexpr auto data() const -> const T * { return mem; }
  [[nodiscard]] static constexpr auto stride() -> Static<size_t, 1> {
    return {};
  }
  constexpr void extendOrAssertSize(size_t M) const { assert(M == size()); }
};

static_assert(
  std::same_as<decltype(std::declval<PtrVector<int64_t, unsigned>>().begin()),
               const int64_t *>);
static_assert(
  std::same_as<decltype(std::declval<PtrVector<int64_t, unsigned>>().end()),
               const int64_t *>);
static_assert(
  sizeof(PtrVector<int64_t, Static<size_t, 16>, Static<ptrdiff_t, 2>>) ==
  sizeof(int64_t *));
static_assert(
  sizeof(PtrVector<int64_t, Static<size_t, 16>, Static<ptrdiff_t, 16>>) ==
  sizeof(int64_t *));

// CRTP for mut vectors
template <typename T, typename V> struct MutVectorCore {
  using eltype = T;
  [[nodiscard]] constexpr auto rawSize() const {
    return static_cast<const V *>(this)->rawSize();
  }
  [[nodiscard]] constexpr auto size() const -> size_t {
    return static_cast<const V *>(this)->size();
  }
  [[nodiscard]] constexpr auto isEmpty() const -> bool { return size() == 0; }
  [[nodiscard]] constexpr auto data() const -> const T * {
    return static_cast<const V *>(this)->data();
  }

  [[nodiscard]] constexpr auto data() -> T * {
    return static_cast<V *>(this)->data();
  }
  [[nodiscard]] constexpr auto stride() const {
    return static_cast<const V *>(this)->stride();
  }

  constexpr auto operator[](const ScalarIndex auto i) -> T & {
    size_t N = size();
#ifndef NDEBUG
    checkIndex(N, i);
#endif
    return data()[canonicalize(i, N) * stride()];
  }
  constexpr auto operator()(const ScalarIndex auto i) -> T & {
    size_t N = size();
#ifndef NDEBUG
    checkIndex(size_t(N), i);
#endif
    return data()[canonicalize(i, N) * stride()];
  }
  constexpr auto operator[](const ScalarIndex auto i) const -> const T & {
    size_t N = size();
#ifndef NDEBUG
    checkIndex(size_t(N), i);
#endif
    return data()[canonicalize(i, N) * stride()];
  }
  constexpr auto operator()(const ScalarIndex auto i) const -> const T & {
    size_t N = size();
#ifndef NDEBUG
    checkIndex(size_t(N), i);
#endif
    return data()[canonicalize(i, N) * stride()];
  }
  [[nodiscard]] auto front() -> T & {
    assert(size() > 0);
    return data()[0];
  }
  [[nodiscard]] auto back() -> T & {
    assert(size() > 0);
    return data()[(size() - 1) * stride()];
  }
  [[nodiscard]] auto front() const -> const T & {
    assert(size() > 0);
    return data()[0];
  }
  [[nodiscard]] auto back() const -> const T & {
    assert(size() > 0);
    return data()[(size() - 1) * stride()];
  }
  // copy constructor
  // MutPtrVector(const MutPtrVector<T> &x) : mem(x.mem), N(x.N) {}
  constexpr auto operator()(Range<size_t, size_t> i) -> MutPtrVector<T> {
    assert(i.b <= i.e);
    assert(i.e <= size());
    return MutPtrVector(data() + i.b * stride(), i.e - i.b, stride());
  }
  constexpr auto operator()(Range<size_t, size_t> i) const -> PtrVector<T> {
    assert(i.b <= i.e);
    assert(i.e <= size());
    return PtrVector(data() + i.b * stride(), i.e - i.b, stride());
  }
  template <typename F, typename L>
  constexpr auto operator()(Range<F, L> i) -> MutPtrVector<T> {
    return (*this)(canonicalizeRange(i, size()));
  }
  template <typename F, typename L>
  constexpr auto operator()(Range<F, L> i) const -> PtrVector<T> {
    return (*this)(canonicalizeRange(i, size()));
  }
  [[nodiscard]] constexpr auto begin() const {
    return default_iterator(data(), stride());
  }
  [[nodiscard]] constexpr auto end() const {
    return default_iterator(data() + size(), stride());
  }
  [[nodiscard]] constexpr auto begin() {
    return default_iterator(data(), stride());
  }
  [[nodiscard]] constexpr auto end() {
    return default_iterator(data() + size(), stride());
  }
  [[nodiscard]] constexpr auto rbegin() const {
    return std::reverse_iterator(end());
  }
  [[nodiscard]] constexpr auto rend() const {
    return std::reverse_iterator(begin());
  }
  [[nodiscard]] constexpr auto rbegin() { return std::reverse_iterator(end()); }
  [[nodiscard]] constexpr auto rend() { return std::reverse_iterator(begin()); }

  constexpr operator PtrVector<T>() const {
    return PtrVector(data(), size(), stride());
  }
  constexpr operator llvm::ArrayRef<T>() const {
    static_assert(std::is_same_v<Static<size_t, 1>, decltype(stride())>);
    return llvm::ArrayRef<T>{data(), size()};
  }
  constexpr operator llvm::MutableArrayRef<T>() {
    static_assert(std::is_same_v<Static<size_t, 1>, decltype(stride())>);
    return llvm::MutableArrayRef<T>{data(), size()};
  }
  // llvm::ArrayRef<T> arrayref() const { return llvm::ArrayRef<T>(ptr, M); }
  auto operator==(const MutPtrVector<T> x) const -> bool {
    return size() == x.size() && std::equal(begin(), end(), x.begin());
  }
  auto operator==(const PtrVector<T> x) const -> bool {
    return size() == x.size() && std::equal(begin(), end(), x.begin());
  }
  auto operator==(const llvm::ArrayRef<T> x) const -> bool {
    return size() == x.size() && std::equal(begin(), end(), x.begin());
  }
  [[nodiscard]] constexpr auto view() const -> PtrVector<T> { return *this; };
  // PtrVector<T> view() const {
  //     return PtrVector<T>{.mem = mem, .N = N};
  // };
  auto operator=(PtrVector<T> x) -> MutPtrVector<T> { return copyto(*this, x); }
  auto operator=(MutPtrVector<T> x) -> MutPtrVector<T> {
    return copyto(*this, x);
  }
  auto operator=(const AbstractVector auto &x) -> MutPtrVector<T> {
    return copyto(*this, x);
  }
  auto operator=(std::integral auto x) -> MutPtrVector<T> {
    for (auto &&y : *this)
      y = x;
    return *this;
  }
  auto operator+=(const AbstractVector auto &x) -> MutPtrVector<T> {
    size_t N = size();
    assert(N == x.size());
    for (size_t i = 0; i < N; ++i)
      data()[i * stride()] += x(i);
    return *this;
  }
  auto operator-=(const AbstractVector auto &x) -> MutPtrVector<T> {
    size_t N = size();
    assert(N == x.size());
    for (size_t i = 0; i < N; ++i)
      data()[i * stride()] -= x(i);
    return *this;
  }
  auto operator*=(const AbstractVector auto &x) -> MutPtrVector<T> {
    size_t N = size();
    assert(N == x.size());
    for (size_t i = 0; i < N; ++i)
      data()[i * stride()] *= x(i);
    return *this;
  }
  auto operator/=(const AbstractVector auto &x) -> MutPtrVector<T> {
    size_t N = size();
    assert(N == x.size());
    for (size_t i = 0; i < N; ++i)
      data()[i * stride()] /= x(i);
    return *this;
  }
  auto operator+=(const std::integral auto x) -> MutPtrVector<T> {
    size_t N = size();
    for (size_t i = 0; i < N; ++i)
      data()[i * stride()] += x;
    return *this;
  }
  auto operator-=(const std::integral auto x) -> MutPtrVector<T> {
    size_t N = size();
    for (size_t i = 0; i < N; ++i)
      data()[i * stride()] -= x;
    return *this;
  }
  auto operator*=(const std::integral auto x) -> MutPtrVector<T> {
    size_t N = size();
    for (size_t i = 0; i < N; ++i)
      data()[i * stride()] *= x;
    return *this;
  }
  auto operator/=(const std::integral auto x) -> MutPtrVector<T> {
    size_t N = size();
    for (size_t i = 0; i < N; ++i)
      data()[i * stride()] /= x;
    return *this;
  }
};
template <typename T, UnsignedIntConvertible I, SignedIntConvertible J>
struct MutPtrVector : MutVectorCore<T, MutPtrVector<T, I, J>> {
  static_assert(!std::is_const_v<T>, "T shouldn't be const");
  using eltype = T;
  // using eltype = std::remove_const_t<T>;
  [[no_unique_address]] T *const mem;
  [[no_unique_address]] const I N;
  [[no_unique_address]] const J X;
  // constexpr MutPtrVector() = default;
  // constexpr MutPtrVector(const MutPtrVector<T> &x) = default;
  // constexpr MutPtrVector(llvm::MutableArrayRef<T> x)
  //   : mem(x.data()), N(x.size()) {}
  // constexpr MutPtrVector(T *mem, size_t N) : mem(mem), N(N) {}
  void extendOrAssertSize(size_t M) const { assert(M == N); }

  [[nodiscard]] constexpr auto data() -> T * { return mem; }
  [[nodiscard]] constexpr auto data() const -> const T * { return mem; }
  [[nodiscard]] constexpr auto size() const -> size_t { return N; }
  [[nodiscard]] constexpr auto rawSize() const -> I { return N; }
  [[nodiscard]] constexpr auto stride() const -> J { return X; }
};

static_assert(std::is_trivially_copyable_v<Static<ptrdiff_t, 1>>);
static_assert(
  std::is_trivially_copyable_v<MutPtrVector<int64_t, unsigned, ptrdiff_t>>);
static_assert(std::is_trivially_copyable_v<MutPtrVector<int64_t, size_t, int>>);
// TODO: should be fixed with Clang 16?
// https://reviews.llvm.org/D140664
// static_assert(std::is_trivially_copyable_v<MutPtrVector<int64_t, size_t,
// Static<ptrdiff_t, 1>>>);
static_assert(
  sizeof(MutPtrVector<int64_t, Static<size_t, 12>, Static<ptrdiff_t, 4>>) ==
  sizeof(int64_t *));
static_assert(
  sizeof(MutPtrVector<int64_t, Static<size_t, 12>, Static<ptrdiff_t, 12>>) ==
  sizeof(int64_t *));

template <typename T, IntConvertible I, IntConvertible J>
PtrVector(T *, I, J) -> PtrVector<T, I, J>;
template <typename T, IntConvertible I, IntConvertible J>
MutPtrVector(T *, I, J) -> MutPtrVector<T, I, J>;

//
// Vectors
//

template <typename T> constexpr auto view(llvm::SmallVectorImpl<T> &x) {
  return MutPtrVector<T>{x.data(), x.size()};
}
template <typename T> constexpr auto view(const llvm::SmallVectorImpl<T> &x) {
  return PtrVector<T>{.mem = x.data(), .N = x.size()};
}
template <typename T> constexpr auto view(llvm::MutableArrayRef<T> x) {
  return MutPtrVector<T>{x.data(), x.size()};
}
template <typename T> constexpr auto view(llvm::ArrayRef<T> x) {
  return PtrVector<T>{.mem = x.data(), .N = x.size()};
}

template <typename T, IntConvertible I = unsigned> struct Vector {
  using eltype = T;
  [[no_unique_address]] llvm::SmallVector<T, 16> data;

  Vector(int N) : data(llvm::SmallVector<T>(N)){};
  Vector(size_t N = 0) : data(llvm::SmallVector<T>(N)){};
  Vector(llvm::SmallVector<T> A) : data(std::move(A)){};

  constexpr auto operator[](const ScalarIndex auto i) -> T & {
    return data[canonicalize(i, data.size())];
  }
  constexpr auto operator()(const ScalarIndex auto i) -> T & {
    return data[canonicalize(i, data.size())];
  }
  constexpr auto operator[](const ScalarIndex auto i) const -> const T & {
    return data[canonicalize(i, data.size())];
  }
  constexpr auto operator()(const ScalarIndex auto i) const -> const T & {
    return data[canonicalize(i, data.size())];
  }
  constexpr auto operator()(Range<size_t, size_t> i) -> MutPtrVector<T> {
    assert(i.b <= i.e);
    assert(i.e <= data.size());
    return MutPtrVector<T>{data.data() + i.b, i.e - i.b};
  }
  constexpr auto operator()(Range<size_t, size_t> i) const -> PtrVector<T> {
    assert(i.b <= i.e);
    assert(i.e <= data.size());
    return PtrVector<T>{.mem = data.data() + i.b, .N = i.e - i.b};
  }
  template <typename F, typename L>
  constexpr auto operator()(Range<F, L> i) -> MutPtrVector<T> {
    return (*this)(canonicalizeRange(i, data.size()));
  }
  template <typename F, typename L>
  constexpr auto operator()(Range<F, L> i) const -> PtrVector<T> {
    return (*this)(canonicalizeRange(i, data.size()));
  }
  constexpr auto operator[](size_t i) -> T & { return data[i]; }
  constexpr auto operator[](size_t i) const -> const T & { return data[i]; }
  // bool operator==(Vector<T, 0> x0) const { return allMatch(*this, x0); }
  constexpr auto begin() { return data.begin(); }
  constexpr auto end() { return data.end(); }
  [[nodiscard]] constexpr auto begin() const { return data.begin(); }
  [[nodiscard]] constexpr auto end() const { return data.end(); }
  [[nodiscard]] constexpr auto size() const -> size_t { return data.size(); }
  // MutPtrVector<T> view() {
  //     return MutPtrVector<T>{.mem = data.data(), .N = data.size()};
  // };
  [[nodiscard]] constexpr auto view() const -> PtrVector<T> {
    return PtrVector<T>{.mem = data.data(), .N = data.size()};
  };
  template <typename A> void push_back(A &&x) {
    data.push_back(std::forward<A>(x));
  }
  template <typename... A> void emplace_back(A &&...x) {
    data.emplace_back(std::forward<A>(x)...);
  }
  Vector(const AbstractVector auto &x) : data(llvm::SmallVector<T>{}) {
    const size_t N = x.size();
    data.resize_for_overwrite(N);
    for (size_t n = 0; n < N; ++n)
      data[n] = x(n);
  }
  void resize(size_t N) { data.resize(N); }
  void resizeForOverwrite(size_t N) { data.resize_for_overwrite(N); }

  constexpr operator MutPtrVector<T>() {
    return MutPtrVector<T>{data.data(), data.size()};
  }
  constexpr operator PtrVector<T>() const {
    return PtrVector<T>{.mem = data.data(), .N = data.size()};
  }
  constexpr operator llvm::MutableArrayRef<T>() {
    return llvm::MutableArrayRef<T>{data.data(), data.size()};
  }
  constexpr operator llvm::ArrayRef<T>() const {
    return llvm::ArrayRef<T>{data.data(), data.size()};
  }
  // MutPtrVector<T> operator=(AbstractVector auto &x) {
  auto operator=(const T &x) -> Vector<T> & {
    MutPtrVector<T> y{*this};
    y = x;
    return *this;
  }
  auto operator=(AbstractVector auto &x) -> Vector<T> & {
    MutPtrVector<T> y{*this};
    y = x;
    return *this;
  }
  auto operator+=(AbstractVector auto &x) -> Vector<T> & {
    MutPtrVector<T> y{*this};
    y += x;
    return *this;
  }
  auto operator-=(AbstractVector auto &x) -> Vector<T> & {
    MutPtrVector<T> y{*this};
    y -= x;
    return *this;
  }
  auto operator*=(AbstractVector auto &x) -> Vector<T> & {
    MutPtrVector<T> y{*this};
    y *= x;
    return *this;
  }
  auto operator/=(AbstractVector auto &x) -> Vector<T> & {
    MutPtrVector<T> y{*this};
    y /= x;
    return *this;
  }
  auto operator+=(const std::integral auto x) -> Vector<T> & {
    for (auto &&y : data)
      y += x;
    return *this;
  }
  auto operator-=(const std::integral auto x) -> Vector<T> & {
    for (auto &&y : data)
      y -= x;
    return *this;
  }
  auto operator*=(const std::integral auto x) -> Vector<T> & {
    for (auto &&y : data)
      y *= x;
    return *this;
  }
  auto operator/=(const std::integral auto x) -> Vector<T> & {
    for (auto &&y : data)
      y /= x;
    return *this;
  }
  template <typename... Ts> Vector(Ts... inputs) : data{inputs...} {}
  void clear() { data.clear(); }
  void extendOrAssertSize(size_t N) const { assert(N == data.size()); }
  void extendOrAssertSize(size_t N) {
    if (N != data.size())
      data.resize_for_overwrite(N);
  }
  auto operator==(const Vector<T> &x) const -> bool {
    return llvm::ArrayRef<T>(*this) == llvm::ArrayRef<T>(x);
  }
  void pushBack(T x) { data.push_back(std::move(x)); }
};

static_assert(std::copyable<Vector<intptr_t>>);
static_assert(AbstractVector<Vector<int64_t>>);
static_assert(!AbstractVector<int64_t>);

template <typename T>
concept DerivedMatrix =
  requires(T t, const T ct) {
    {
      t.data()
      } -> std::convertible_to<typename std::add_pointer_t<
        typename std::add_const_t<typename T::eltype>>>;
    {
      ct.data()
      } -> std::same_as<typename std::add_pointer_t<
        typename std::add_const_t<typename T::eltype>>>;
    { t.numRow() } -> IsRow;
    { t.numCol() } -> IsCol;
    { t.rowStride() } -> IsRowStride;
  };

template <typename T, IntConvertible I = unsigned, IntConvertible J = unsigned,
          IntConvertible K = int>
struct PtrMatrix;
template <typename T, IntConvertible I = unsigned, IntConvertible J = unsigned,
          IntConvertible K = int>
struct MutPtrMatrix;

template <typename T>
concept ScalarRowIndex = ScalarIndex<T> || IsRow<T>;
template <typename T>
concept ScalarColIndex = ScalarIndex<T> || IsCol<T>;

constexpr auto unwrapRow(IsRow auto x) -> size_t { return size_t(x); }
constexpr auto unwrapCol(IsCol auto x) -> size_t { return size_t(x); }
constexpr auto unwrapRow(auto x) { return x; }
constexpr auto unwrapCol(auto x) { return x; }
constexpr auto unwrap(RowOrCol auto x) { return *x; }
constexpr auto unwrap(auto x) { return x; }

template <typename T>
constexpr inline auto
matrixGet(T *ptr, IsRow auto M, IsCol auto N, IsRowStride auto X,
          const ScalarRowIndex auto mm, const ScalarColIndex auto nn) -> T & {
  auto m = unwrapRow(mm);
  auto n = unwrapCol(nn);
#ifndef NDEBUG
  checkIndex(size_t(M), m);
  checkIndex(size_t(N), n);
#endif
  return *(ptr + (canonicalize(n, size_t(N)) + X * canonicalize(m, size_t(M))));
}
template <typename T>
constexpr inline auto matrixGet(const T *ptr, IsRow auto M, IsCol auto N,
                                IsRowStride auto X,
                                const ScalarRowIndex auto mm,
                                const ScalarColIndex auto nn) -> const T & {
  auto m = unwrapRow(mm);
  auto n = unwrapCol(nn);
#ifndef NDEBUG
  checkIndex(size_t(M), m);
  checkIndex(size_t(N), n);
#endif
  return *(ptr + (canonicalize(n, size_t(N)) + X * canonicalize(m, size_t(M))));
}

template <typename T>
concept AbstractSlice = requires(T t, size_t M) {
                          {
                            canonicalizeRange(t, M)
                            } -> std::same_as<Range<size_t, size_t>>;
                        };

template <typename T>
inline constexpr auto matrixGet(const T *ptr, IsRow auto M, IsCol auto N,
                                IsRowStride auto X, const AbstractSlice auto m,
                                const AbstractSlice auto n) -> PtrMatrix<T> {
#ifndef NDEBUG
  checkIndex(size_t(M), m);
  checkIndex(size_t(N), n);
#endif
  IsRange auto mr = canonicalizeRange(m, M);
  IsRange auto nr = canonicalizeRange(n, N);
  return PtrMatrix<T>{ptr + nr.b + X * mr.b, mr.e - mr.b, nr.e - nr.b, X};
}
template <typename T>
inline constexpr auto matrixGet(T *ptr, IsRow auto M, IsCol auto N,
                                IsRowStride auto X, const AbstractSlice auto m,
                                const AbstractSlice auto n) -> MutPtrMatrix<T> {
#ifndef NDEBUG
  checkIndex(size_t(M), m);
  checkIndex(size_t(N), n);
#endif
  IsRange auto mr = canonicalizeRange(m, M);
  IsRange auto nr = canonicalizeRange(n, N);
  return MutPtrMatrix<T>{ptr + nr.b + X * mr.b,
                         axis<AxisType::Row>(mr.e - mr.b),
                         axis<AxisType::Column>(nr.e - nr.b), X};
}

template <typename T>
inline constexpr auto matrixGet(const T *ptr, IsRow auto M, IsCol auto N,
                                IsRowStride auto X,
                                const ScalarRowIndex auto mm,
                                const AbstractSlice auto n) -> PtrVector<T> {
  auto m = unwrapRow(mm);
#ifndef NDEBUG
  checkIndex(size_t(M), m);
  checkIndex(size_t(N), n);
#endif
  size_t mi = canonicalize(m, size_t(M));
  IsRange auto nr = canonicalizeRange(n, size_t(N));
  return PtrVector(ptr + nr.b + X * mi, nr.e - nr.b);
}
template <typename T>
inline constexpr auto
matrixGet(T *ptr, IsRow auto M, IsCol auto N, IsRowStride auto X,
          const ScalarRowIndex auto mm, const AbstractSlice auto n) {
  auto m = unwrapRow(mm);
#ifndef NDEBUG
  checkIndex(size_t(M), m);
  checkIndex(size_t(N), n);
#endif
  size_t mi = canonicalize(m, size_t(M));
  Range<size_t, size_t> nr = canonicalizeRange(n, size_t(N));
  return MutPtrVector(ptr + nr.b + X * mi, nr.e - nr.b);
}

template <typename T>
inline constexpr auto matrixGet(const T *ptr, IsRow auto M, IsCol auto N,
                                IsRowStride auto X, const AbstractSlice auto m,
                                const ScalarColIndex auto nn) {
  auto n = unwrapCol(nn);
#ifndef NDEBUG
  checkIndex(size_t(M), m);
  checkIndex(size_t(N), n);
#endif
  IsRange auto mr = canonicalizeRange(m, size_t(M));
  size_t ni = canonicalize(n, size_t(N));
  return PtrVector{ptr + ni + X * mr.b, mr.e - mr.b, X};
}
template <typename T>
inline constexpr auto matrixGet(T *ptr, IsRow auto M, IsCol auto N,
                                IsRowStride auto X, const AbstractSlice auto m,
                                const ScalarColIndex auto nn) {
  auto n = unwrapCol(nn);
#ifndef NDEBUG
  checkIndex(size_t(M), m);
  checkIndex(size_t(N), n);
#endif
  Range<size_t, size_t> mr = canonicalizeRange(m, size_t(M));
  size_t ni = canonicalize(n, size_t(N));
  return MutPtrVector<T>{ptr + ni + X * mr.b, mr.e - mr.b, X};
}

constexpr auto isSquare(const AbstractMatrix auto &A) -> bool {
  return A.numRow() == A.numCol();
}

template <typename T> constexpr auto diag(MutPtrMatrix<T> A) {
  return MutPtrVector{A.data(), A.minRowCol(),
                      A.rowStride() + Static<ptrdiff_t, 1>()};
}
template <typename T> constexpr auto diag(PtrMatrix<T> A) {
  return PtrVector{A.data(), A.minRowCol(),
                   A.rowStride() + Static<ptrdiff_t, 1>()};
}
template <typename T> constexpr auto antiDiag(MutPtrMatrix<T> A) {
  return MutPtrVector{A.data() + size_t(A.numCol()) - 1, A.minRowCol(),
                      (A.rowStride() - Static<ptrdiff_t, 1>())};
}
template <typename T> constexpr auto antiDiag(PtrMatrix<T> A) {
  return PtrVector{A.data() + size_t(A.numCol()) - 1, A.minRowCol(),
                   (A.rowStride() - Static<ptrdiff_t, 1>())};
}

/// A CRTP type defining const methods for matrices.
template <typename T, typename A> struct ConstMatrixCore {
  [[nodiscard]] constexpr auto data() const -> const T * {
    return static_cast<const A *>(this)->data();
  }
  [[nodiscard]] constexpr auto numRow() const {
    return static_cast<const A *>(this)->numRow();
  }
  [[nodiscard]] constexpr auto numCol() const {
    return static_cast<const A *>(this)->numCol();
  }
  [[nodiscard]] constexpr auto rowStride() const {
    return static_cast<const A *>(this)->rowStride();
  }

  constexpr auto operator()(const ScalarRowIndex auto m,
                            const ScalarColIndex auto n) const -> const T & {
    return matrixGet(data(), numRow(), numCol(), rowStride(), m, n);
  }
  constexpr auto operator()(auto m, auto n) const {
    return matrixGet(data(), numRow(), numCol(), rowStride(), m, n);
  }
  [[nodiscard]] constexpr auto size() const {
    return std::make_pair(numRow(), numCol());
  }
  [[nodiscard]] constexpr auto diag() const {
    return LinearAlgebra::diag(PtrMatrix<T>(*this));
  }
  [[nodiscard]] constexpr auto antiDiag() const {
    return LinearAlgebra::antiDiag(PtrMatrix<T>(*this));
  }
  [[nodiscard]] constexpr auto minRowCol() const -> size_t {
    return std::min(size_t(numRow()), size_t(numCol()));
  }
  [[nodiscard]] constexpr auto isSquare() const -> bool {
    return size_t(numRow()) == size_t(numCol());
  }
  template <IntConvertible I, IntConvertible J, IntConvertible K>
  constexpr operator PtrMatrix<T, I, J, K>() const {
    const T *ptr = data();
    return PtrMatrix<T>(ptr, numRow(), numCol(), rowStride());
  }
  [[nodiscard]] constexpr auto view() const -> PtrMatrix<T> {
    const T *ptr = data();
    return PtrMatrix<T>(ptr, numRow(), numCol(), rowStride());
  }
  [[nodiscard]] constexpr auto transpose() const -> Transpose<PtrMatrix<T>> {
    return Transpose<PtrMatrix<eltype_t<A>>>{view()};
  }
};
template <typename T, typename A> struct MutMatrixCore {

  [[nodiscard]] constexpr auto numRow() const {
    return static_cast<const A *>(this)->numRow();
  }
  [[nodiscard]] constexpr auto numCol() const {
    return static_cast<const A *>(this)->numCol();
  }
  [[nodiscard]] constexpr auto rowStride() const {
    return static_cast<const A *>(this)->rowStride();
  }
  [[nodiscard]] constexpr auto data() -> T * {
    return static_cast<A *>(this)->data();
  }

  constexpr auto operator()(const ScalarRowIndex auto m,
                            const ScalarColIndex auto n) -> T & {
    return matrixGet(data(), numRow(), numCol(), rowStride(), m, n);
  }
  constexpr auto operator()(auto m, auto n) {
    return matrixGet(data(), numRow(), numCol(), rowStride(), m, n);
  }
  constexpr auto diag() {
    return LinearAlgebra::diag(MutPtrMatrix<T>(*static_cast<A *>(this)));
  }
  constexpr auto antiDiag() {
    return LinearAlgebra::antiDiag(MutPtrMatrix<T>(*static_cast<A *>(this)));
  }
  template <IntConvertible I, IntConvertible J, IntConvertible K>
  constexpr operator MutPtrMatrix<T, I, J, K>() {
    return MutPtrMatrix<T>(data(), numRow(), numCol(), rowStride());
  }
  // [[nodiscard]] constexpr auto view() -> MutPtrMatrix<T> {
  //   T *ptr = data();
  //   return MutPtrMatrix<T>{ptr, numRow(), numCol(), rowStride()};
  // }
};

template <typename T, IntConvertible I = unsigned, IntConvertible J = unsigned,
          size_t S = 64>
struct Matrix;

template <typename T, IntConvertible I = size_t> struct SmallSparseMatrix;
template <typename T, IntConvertible I, IntConvertible J, IntConvertible K>
struct PtrMatrix : ConstMatrixCore<T, PtrMatrix<T, I, J, K>> {
  using CBase = ConstMatrixCore<T, PtrMatrix<T, I, J, K>>;
  using CBase::size, CBase::diag, CBase::antiDiag, CBase::operator();

  using eltype = std::remove_reference_t<T>;
  static_assert(!std::is_const_v<T>, "const T is redundant");

  [[no_unique_address]] const T *const mem;
  [[no_unique_address]] I M;
  [[no_unique_address]] J N;
  [[no_unique_address]] K X;
  // [[no_unique_address]] Row M;
  // [[no_unique_address]] Col N;
  // [[no_unique_address]] RowStride X;

  [[nodiscard]] constexpr auto data() const -> const T * { return mem; }
  [[nodiscard]] constexpr auto numRow() const { return toRow(M); }
  [[nodiscard]] constexpr auto numCol() const { return toCol(N); }
  [[nodiscard]] constexpr auto rowStride() const { return toRowStride(X); }
  [[nodiscard]] constexpr auto isSquare() const -> bool {
    return size_t(M) == size_t(N);
  }
  [[nodiscard]] constexpr inline auto view() const -> PtrMatrix<T> {
    return *this;
  };
  [[nodiscard]] constexpr auto transpose() const -> Transpose<PtrMatrix<T>> {
    return Transpose<PtrMatrix<T>>{*this};
  }
  constexpr PtrMatrix(const T *const mem, const IsRow auto M,
                      const IsCol auto N, const IsRowStride auto X)
    : mem(mem), M(M), N(N), X(X) {}
  constexpr void extendOrAssertSize(IsRow auto MM, IsCol auto NN) const {
    assert(MM == M);
    assert(NN == N);
  }
};
static_assert(std::same_as<PtrMatrix<int64_t>::eltype, int64_t>);
static_assert(HasEltype<PtrMatrix<int64_t>>);
template <typename T, IntConvertible I, IntConvertible J, IntConvertible K>
struct MutPtrMatrix : ConstMatrixCore<T, MutPtrMatrix<T, I, J, K>>,
                      MutMatrixCore<T, MutPtrMatrix<T, I, J, K>> {
  using CBase = ConstMatrixCore<T, MutPtrMatrix<T, I, J, K>>;
  using MBase = MutMatrixCore<T, MutPtrMatrix<T, I, J, K>>;

  using CBase::antiDiag, MBase::antiDiag;
  using CBase::diag, MBase::diag;
  using CBase::operator(), MBase::operator();
  using CBase::size, CBase::view, CBase::isSquare, CBase::transpose,
    CBase::operator ::LinearAlgebra::PtrMatrix<T, I, J, J>,
    MBase::operator ::LinearAlgebra::MutPtrMatrix<T, I, J, J>;

  using eltype = std::remove_reference_t<T>;
  static_assert(!std::is_const_v<T>, "MutPtrMatrix should never have const T");
  [[no_unique_address]] T *const mem;
  [[no_unique_address]] I M;
  [[no_unique_address]] J N;
  [[no_unique_address]] K X;
  // [[no_unique_address]] Col N;
  // [[no_unique_address]] RowStride X;

  [[nodiscard]] constexpr auto data() -> T * { return mem; }
  [[nodiscard]] constexpr auto data() const -> const T * { return mem; }
  [[nodiscard]] constexpr auto numRow() const { return toRow(M); }
  [[nodiscard]] constexpr auto numCol() const { return toCol(N); }
  [[nodiscard]] constexpr auto rowStride() const { return toRowStride(X); }
  [[nodiscard]] constexpr auto view() const -> PtrMatrix<T> {
    return PtrMatrix<T>(data(), M, N, X);
  };
  constexpr operator PtrMatrix<T>() const {
    return PtrMatrix<T>(data(), M, N, X);
  }
  template <size_t S>
  MutPtrMatrix(Matrix<T, I, J, S> &mat)
    : mem(*mat.data()), M(*mat.numRow()), N(*mat.numCol()),
      X(*mat.rowStride()) {}
  MutPtrMatrix(Matrix<T, I, J, 64> &mat)
    : mem(*mat.data()), M(*mat.numRow()), N(*mat.numCol()),
      X(*mat.rowStride()) {}

  auto operator=(const SmallSparseMatrix<T> &A) -> MutPtrMatrix<T> {
    assert(M == A.numRow());
    assert(N == A.numCol());
    size_t k = 0;
    for (size_t i = 0; i < M; ++i) {
      uint32_t m = A.rows[i] & 0x00ffffff;
      size_t j = 0;
      while (m) {
        uint32_t tz = std::countr_zero(m);
        m >>= tz + 1;
        j += tz;
        mem[X * i + (j++)] = A.nonZeros[k++];
      }
    }
    assert(k == A.nonZeros.size());
    return *this;
  }
  auto operator=(MutPtrMatrix<T> A) -> MutPtrMatrix<T> {
    return copyto(*this, PtrMatrix<T>(A));
  }
  // rule of 5 requires...
  // constexpr MutPtrMatrix(const MutPtrMatrix<T> &A) = default;
  constexpr MutPtrMatrix(T *mem, IsRow auto MM, IsCol auto NN)
    : mem(mem), M(*MM), N(*NN), X(*NN){};
  constexpr MutPtrMatrix(T *mem, IsRow auto MM, IsCol auto NN,
                         IsRowStride auto XX)
    : mem(mem), M(*MM), N(*NN), X(*XX){};
  // template <typename ARM>
  // constexpr MutPtrMatrix(ARM &A)
  //   : mem(A.data()), M(A.numRow()), N(A.numCol()), X(A.rowStride()) {}

  auto operator=(const AbstractMatrix auto &B) -> MutPtrMatrix<T> {
    return copyto(*this, B);
  }
  auto operator=(const std::integral auto b) -> MutPtrMatrix<T> {
    for (size_t r = 0; r < M; ++r)
      for (size_t c = 0; c < N; ++c)
        (*this)(r, c) = b;
    return *this;
  }
  auto operator+=(const AbstractMatrix auto &B) -> MutPtrMatrix<T> {
    assert(M == B.numRow());
    assert(N == B.numCol());
    for (size_t r = 0; r < M; ++r)
      for (size_t c = 0; c < N; ++c)
        (*this)(r, c) += B(r, c);
    return *this;
  }
  auto operator-=(const AbstractMatrix auto &B) -> MutPtrMatrix<T> {
    assert(M == B.numRow());
    assert(N == B.numCol());
    for (size_t r = 0; r < M; ++r)
      for (size_t c = 0; c < N; ++c)
        (*this)(r, c) -= B(r, c);
    return *this;
  }
  auto operator*=(const std::integral auto b) -> MutPtrMatrix<T> {
    for (size_t r = 0; r < M; ++r)
      for (size_t c = 0; c < N; ++c)
        (*this)(r, c) *= b;
    return *this;
  }
  auto operator/=(const std::integral auto b) -> MutPtrMatrix<T> {
    const size_t M = numRow();
    const size_t N = numCol();
    for (size_t r = 0; r < M; ++r)
      for (size_t c = 0; c < N; ++c)
        (*this)(r, c) /= b;
    return *this;
  }
  [[nodiscard]] constexpr auto isSquare() const -> bool {
    return size_t(M) == size_t(N);
  }
  [[nodiscard]] constexpr auto transpose() const -> Transpose<PtrMatrix<T>> {
    return Transpose<PtrMatrix<T>>{view()};
  }
  constexpr void extendOrAssertSize(IsRow auto M, IsCol auto N) const {
    assert(numRow() == M);
    assert(numCol() == N);
  }
};
template <typename T> constexpr auto ptrVector(T *p, size_t M) {
  if constexpr (std::is_const_v<T>)
    return PtrVector<std::remove_const_t<T>>{.mem = p, .N = M};
  else
    return MutPtrVector<T>{p, M};
}

template <typename T> PtrMatrix(T *, size_t, size_t) -> PtrMatrix<T>;
template <typename T> MutPtrMatrix(T *, size_t, size_t) -> MutPtrMatrix<T>;
template <typename T> PtrMatrix(T *, size_t, size_t, size_t) -> PtrMatrix<T>;
template <typename T>
MutPtrMatrix(T *, size_t, size_t, size_t) -> MutPtrMatrix<T>;

template <AbstractRowMajorMatrix T> PtrMatrix(T &A) -> PtrMatrix<eltype_t<T>>;
template <AbstractRowMajorMatrix T>
MutPtrMatrix(T &A) -> MutPtrMatrix<eltype_t<T>>;

// template <typename T>
// constexpr auto ptrmat(T *ptr, size_t numRow, size_t numCol, size_t stride) {
//     if constexpr (std::is_const_v<T>) {
//         return PtrMatrix<std::remove_const_t<T>>{
//             .mem = ptr, .M = numRow, .N = numCol, .X = stride};
//     } else {
//         return MutPtrMatrix<T>{
//             .mem = ptr, .M = numRow, .N = numCol, .X = stride};
//     }
// }
static_assert(sizeof(PtrMatrix<int64_t>) <=
              4 * sizeof(unsigned int) + sizeof(int64_t *));
static_assert(sizeof(MutPtrMatrix<int64_t>) <=
              4 * sizeof(unsigned int) + sizeof(int64_t *));
static_assert(std::is_trivially_copyable_v<Row<>>);
static_assert(std::is_trivially_copyable_v<Col<Static<size_t, 3>>>);
static_assert(std::is_trivially_copyable_v<RowStride<ptrdiff_t>>);
static_assert(std::is_trivially_copyable_v<const Row<>>);
static_assert(std::is_trivially_copyable_v<const Col<Static<size_t, 77>>>);
static_assert(std::is_trivially_copyable_v<const RowStride<ptrdiff_t>>);
static_assert(std::is_trivially_copyable_v<PtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> is not trivially copyable!");
// Should be fixed with Clang 16:
// https://reviews.llvm.org/D140664
// static_assert(std::is_trivially_copyable_v<PtrVector<int64_t>>,
//               "PtrVector<int64_t,0> is not trivially copyable!");
// static_assert(std::is_trivially_copyable_v<MutPtrMatrix<int64_t>>,
//               "MutPtrMatrix<int64_t> is not trivially copyable!");

static_assert(!AbstractVector<PtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> isa AbstractVector succeeded");
static_assert(!AbstractVector<MutPtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> isa AbstractVector succeeded");
static_assert(!AbstractVector<const PtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> isa AbstractVector succeeded");

static_assert(AbstractMatrix<PtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> isa AbstractMatrix failed");
static_assert(
  std::same_as<std::remove_reference_t<decltype(MutPtrMatrix<int64_t>(
                 nullptr, Row<>{0}, Col<>{0})(size_t(0), size_t(0)))>,
               int64_t>);

static_assert(AbstractMatrix<MutPtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> isa AbstractMatrix failed");
static_assert(AbstractMatrix<const PtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> isa AbstractMatrix failed");
static_assert(AbstractMatrix<const MutPtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> isa AbstractMatrix failed");

static_assert(AbstractVector<MutPtrVector<int64_t>>,
              "PtrVector<int64_t> isa AbstractVector failed");
static_assert(AbstractVector<PtrVector<int64_t>>,
              "PtrVector<const int64_t> isa AbstractVector failed");
static_assert(AbstractVector<const PtrVector<int64_t>>,
              "PtrVector<const int64_t> isa AbstractVector failed");
static_assert(AbstractVector<const MutPtrVector<int64_t>>,
              "PtrVector<const int64_t> isa AbstractVector failed");

static_assert(AbstractVector<Vector<int64_t>>,
              "PtrVector<int64_t> isa AbstractVector failed");

static_assert(!AbstractMatrix<MutPtrVector<int64_t>>,
              "PtrVector<int64_t> isa AbstractMatrix succeeded");
static_assert(!AbstractMatrix<PtrVector<int64_t>>,
              "PtrVector<const int64_t> isa AbstractMatrix succeeded");
static_assert(!AbstractMatrix<const PtrVector<int64_t>>,
              "PtrVector<const int64_t> isa AbstractMatrix succeeded");
static_assert(!AbstractMatrix<const MutPtrVector<int64_t>>,
              "PtrVector<const int64_t> isa AbstractMatrix succeeded");

static_assert(
  AbstractMatrix<ElementwiseMatrixBinaryOp<Mul, PtrMatrix<int64_t>, int>>,
  "ElementwiseBinaryOp isa AbstractMatrix failed");

static_assert(
  !AbstractVector<MatMatMul<PtrMatrix<int64_t>, PtrMatrix<int64_t>>>,
  "MatMul should not be an AbstractVector!");
static_assert(AbstractMatrix<MatMatMul<PtrMatrix<int64_t>, PtrMatrix<int64_t>>>,
              "MatMul is not an AbstractMatrix!");

template <typename T>
concept IntVector = requires(T t, int64_t y) {
                      { t.size() } -> std::convertible_to<size_t>;
                      { t[y] } -> std::convertible_to<int64_t>;
                    };

template <typename T, IntConvertible I>
struct SquarePtrMatrix : ConstMatrixCore<T, SquarePtrMatrix<T, I>> {
  using Base = ConstMatrixCore<T, SquarePtrMatrix<T, I>>;
  using Base::diag, Base::antiDiag,
    Base::operator(), Base::size, Base::view, Base::isSquare, Base::transpose,
    Base::operator ::LinearAlgebra::PtrMatrix<T>;
  using eltype = std::remove_reference_t<T>;
  static_assert(!std::is_const_v<T>, "const T is redundant");
  [[no_unique_address]] const T *const mem;
  [[no_unique_address]] const I M;
  constexpr SquarePtrMatrix(const T *const mem, size_t m) : mem(mem), M(m){};

  [[nodiscard]] constexpr auto numRow() const { return toRow(M); }
  [[nodiscard]] constexpr auto numCol() const { return toCol(M); }
  [[nodiscard]] constexpr auto rowStride() const { return toRowStride(M); }
  constexpr auto data() -> const T * { return mem; }
  [[nodiscard]] constexpr auto data() const -> const T * { return mem; }
};

template <typename T, IntConvertible I>
struct MutSquarePtrMatrix : ConstMatrixCore<T, MutSquarePtrMatrix<T, I>>,
                            MutMatrixCore<T, MutSquarePtrMatrix<T, I>> {
  using CBase = ConstMatrixCore<T, MutSquarePtrMatrix<T, I>>;
  using MBase = MutMatrixCore<T, MutSquarePtrMatrix<T, I>>;
  using CBase::diag, CBase::antiDiag, MBase::diag, MBase::antiDiag,
    CBase::operator(), MBase::operator(), CBase::size, CBase::view,
    CBase::isSquare, CBase::transpose,
    CBase::operator ::LinearAlgebra::PtrMatrix<T, I, I, I>,
    MBase::operator ::LinearAlgebra::MutPtrMatrix<T, I, I, I>;
  using eltype = std::remove_reference_t<T>;
  static_assert(!std::is_const_v<T>, "T should not be const");
  [[no_unique_address]] T *const mem;
  [[no_unique_address]] const I M;

  [[nodiscard]] constexpr auto numRow() const { return toRow(M); }
  [[nodiscard]] constexpr auto numCol() const { return toCol(M); }
  [[nodiscard]] constexpr auto rowStride() const { return toRowStride(M); }
  MutSquarePtrMatrix(T *mem, size_t m) : mem(mem), M(m){};
  constexpr auto data() -> T * { return mem; }
  [[nodiscard]] constexpr auto data() const -> const T * { return mem; }
  constexpr operator SquarePtrMatrix<T, I>() const { return {mem, M}; }
  auto operator=(const AbstractMatrix auto &B) -> MutSquarePtrMatrix<T, I> {
    return copyto(*this, B);
  }
};

template <typename T, IntConvertible I = size_t, unsigned STORAGE = 8>
struct SquareMatrix : ConstMatrixCore<T, SquareMatrix<T, I, STORAGE>>,
                      MutMatrixCore<T, SquareMatrix<T, I, STORAGE>> {
  using CBase = ConstMatrixCore<T, SquareMatrix<T, I, STORAGE>>;
  using MBase = MutMatrixCore<T, SquareMatrix<T, I, STORAGE>>;
  using CBase::diag, CBase::antiDiag, MBase::diag, MBase::antiDiag,
    CBase::operator(), MBase::operator(), CBase::size, CBase::view,
    CBase::isSquare, CBase::transpose,
    CBase::operator ::LinearAlgebra::PtrMatrix<T, I, I, I>,
    MBase::operator ::LinearAlgebra::MutPtrMatrix<T, I, I, I>;
  using eltype = std::remove_reference_t<T>;
  static constexpr unsigned TOTALSTORAGE = STORAGE * STORAGE;
  [[no_unique_address]] llvm::SmallVector<T, TOTALSTORAGE> mem;
  [[no_unique_address]] I M;

  SquareMatrix(size_t m)
    : mem(llvm::SmallVector<T, TOTALSTORAGE>(m * m)), M(m){};

  [[nodiscard]] constexpr auto numRow() const { return toRow(M); }
  [[nodiscard]] constexpr auto numCol() const { return toCol(M); }
  [[nodiscard]] constexpr auto rowStride() const { return toRowStride(M); }

  constexpr auto data() -> T * { return mem.data(); }
  [[nodiscard]] constexpr auto data() const -> const T * { return mem.data(); }

  constexpr auto begin() -> T * { return data(); }
  constexpr auto end() -> T * { return data() + M * M; }
  [[nodiscard]] constexpr auto begin() const -> const T * { return data(); }
  [[nodiscard]] constexpr auto end() const -> const T * {
    return data() + M * M;
  }
  auto operator[](size_t i) -> T & { return mem[i]; }
  auto operator[](size_t i) const -> const T & { return mem[i]; }

  static auto identity(size_t N) -> SquareMatrix<T, I> {
    SquareMatrix<T, I> A(N);
    for (size_t r = 0; r < N; ++r)
      A(r, r) = 1;
    return A;
  }
  inline static auto identity(Row<I> N) -> SquareMatrix<T, I> {
    return identity(size_t(N));
  }
  inline static auto identity(Col<I> N) -> SquareMatrix<T, I> {
    return identity(size_t(N));
  }
  constexpr operator MutSquarePtrMatrix<T, I>() { return {mem.data(), M}; }
  constexpr operator SquarePtrMatrix<T, I>() const { return {mem.data(), M}; }
};

template <typename T, IntConvertible I, IntConvertible J, size_t S>
struct Matrix : ConstMatrixCore<T, Matrix<T, I, J, S>>,
                MutMatrixCore<T, Matrix<T, I, J, S>> {
  using CBase = ConstMatrixCore<T, Matrix<T, I, J, S>>;
  using MBase = MutMatrixCore<T, Matrix<T, I, J, S>>;
  using CBase::antiDiag, MBase::antiDiag;
  using CBase::diag, MBase::diag;
  using CBase::operator(), MBase::operator();
  using CBase::size, CBase::view, CBase::isSquare, CBase::transpose,
    CBase::operator ::LinearAlgebra::PtrMatrix<T, I, J, J>,
    MBase::operator ::LinearAlgebra::MutPtrMatrix<T, I, J, J>;

  using eltype = std::remove_reference_t<T>;
  [[no_unique_address]] llvm::SmallVector<T, S> mem;

  [[no_unique_address]] I M = 0;
  [[no_unique_address]] J N = 0;
  [[no_unique_address]] J X = 0;

  constexpr auto data() -> T * { return mem.data(); }
  [[nodiscard]] constexpr auto data() const -> const T * { return mem.data(); }
  Matrix(llvm::SmallVector<T, S> content, IsRow auto M, IsCol auto N)
    : mem(std::move(content)), M(*M), N(*N), X(rowStride(*N)){};

  Matrix(IsRow auto M, IsCol auto N)
    : mem(llvm::SmallVector<T, S>(M * N)), M(*M), N(*N), X(rowStride(*N)){};

  Matrix() = default;
  Matrix(SquareMatrix<T, I> &&A)
    : mem(std::move(A.mem)), M(A.M), N(A.M), X(A.M){};
  Matrix(const SquareMatrix<T, I> &A)
    : mem(A.begin(), A.end()), M(A.M), N(A.M), X(A.M){};
  Matrix(const AbstractMatrix auto &A)
    : mem(llvm::SmallVector<T, S>{}), M(A.numRow()), N(A.numCol()),
      X(A.numCol()) {
    mem.resize_for_overwrite(M * N);
    for (size_t m = 0; m < M; ++m)
      for (size_t n = 0; n < N; ++n)
        mem[size_t(X * m + n)] = A(m, n);
  }
  constexpr auto begin() { return mem.begin(); }
  constexpr auto end() { return mem.begin() + rowStride() * M; }
  [[nodiscard]] constexpr auto begin() const { return mem.begin(); }
  [[nodiscard]] constexpr auto end() const {
    return mem.begin() + rowStride() * M;
  }
  [[nodiscard]] constexpr auto numRow() const { return toRow(M); }
  [[nodiscard]] constexpr auto numCol() const { return toCol(N); }
  [[nodiscard]] constexpr auto rowStride() const { return toRowStride(X); }

  operator MutPtrMatrix<T, I, J, J>() { return {mem.data(), M, N, X}; }

  static auto uninitialized(Row<I> MM, Col<J> NN) -> Matrix<T, I, J, S> {
    Matrix<T, I, J, S> A(0, 0);
    A.M = MM;
    A.X = A.N = NN;
    A.mem.resize_for_overwrite(MM * NN);
    return A;
  }
  static auto identity(size_t MM) -> Matrix<T, I, J, S> {
    Matrix<T, I, J, S> A(MM, MM);
    for (size_t i = 0; i < MM; ++i)
      A(i, i) = 1;
    return A;
  }
  inline static auto identity(Row<I> N) -> Matrix<T, I, J, S> {
    return identity((*N));
  }
  inline static auto identity(Col<J> N) -> Matrix<T, I, J, S> {
    return identity((*N));
  }
  void clear() {
    M = N = X = 0;
    mem.clear();
  }

  void resize(Row<I> MM, Col<J> NN, RowStride<J> XX) {
    mem.resize(XX * MM);
    size_t minMMM = std::min(size_t(M), size_t(MM));
    if ((XX > X) && M && N)
      // need to copy
      for (size_t m = minMMM - 1; m > 0; --m)
        for (auto n = size_t(N); n-- > 0;)
          mem[size_t(XX * m + n)] = mem[size_t(X * m + n)];
    // zero
    for (size_t m = 0; m < minMMM; ++m)
      for (auto n = size_t(N); n < size_t(NN); ++n)
        mem[size_t(XX * m + n)] = 0;
    for (size_t m = minMMM; m < size_t(MM); ++m)
      for (size_t n = 0; n < size_t(NN); ++n)
        mem[size_t(XX * m + n)] = 0;
    X = *XX;
    M = *MM;
    N = *NN;
  }
  void insertZero(Col<J> i) {
    Col<J> NN = N + 1;
    RowStride<J> XX = rowStride(max(rowStride(X, NN)));
    mem.resize(XX * M);
    size_t nLower = (XX > X) ? size_t(0) : size_t(i);
    if (M && N)
      // need to copy
      for (auto m = size_t(M); m-- > 0;)
        for (auto n = size_t(N); n-- > nLower;)
          mem[size_t(XX * m + n) + (n >= size_t(i))] = mem[size_t(X * m + n)];
    // zero
    for (size_t m = 0; m < M; ++m)
      mem[size_t(XX * m) + size_t(i)] = 0;
    X = *XX;
    N = *NN;
  }
  void resize(Row<I> MM, Col<J> NN) { resize(MM, NN, max(NN, X)); }
  void reserve(Row<I> MM, Col<J> NN) { reserve(MM, max(NN, X)); }
  void reserve(Row<I> MM, RowStride<J> NN) { mem.reserve(NN * MM); }
  void clearReserve(Row<I> MM, Col<J> NN) { clearReserve(MM, rowStride(*NN)); }
  void clearReserve(Row<I> MM, RowStride<J> XX) {
    clear();
    mem.reserve(XX * MM);
  }
  void resizeForOverwrite(Row<I> MM, Col<J> NN, RowStride<J> XX) {
    assert(XX >= NN);
    M = *MM;
    N = *NN;
    X = *XX;
    if (X * M > mem.size())
      mem.resize_for_overwrite(X * M);
  }
  void resizeForOverwrite(Row<I> MM, Col<J> NN) {
    M = *MM;
    N = X = *NN;
    if (X * M > mem.size())
      mem.resize_for_overwrite(X * M);
  }

  void resize(Row<I> MM) {
    Row<I> Mold = M;
    M = *MM;
    if (rowStride() * M > mem.size())
      mem.resize(X * M);
    if (M > Mold)
      (*this)(_(Mold, M), _) = 0;
  }
  void resizeForOverwrite(Row<I> MM) {
    if (rowStride() * MM > mem.size())
      mem.resize_for_overwrite(X * M);
    M = *MM;
  }
  void resize(Col<J> NN) { resize(M, NN); }
  void resizeForOverwrite(Col<J> NN) {
    if (X < NN) {
      X = *NN;
      mem.resize_for_overwrite(X * M);
    }
    N = *NN;
  }
  void erase(Col<J> i) {
    assert(i < N);
    for (size_t m = 0; m < M; ++m)
      for (auto n = size_t(i); n < N - 1; ++n)
        (*this)(m, n) = (*this)(m, n + 1);
    --N;
  }
  void erase(Row<J> i) {
    assert(i < M);
    auto it = mem.begin() + X * i;
    mem.erase(it, it + X);
    --M;
  }
  constexpr void truncate(Col<J> NN) {
    assert(NN <= N);
    N = *NN;
  }
  constexpr void truncate(Row<I> MM) {
    assert(MM <= M);
    M = *MM;
  }
  auto operator=(T x) -> Matrix<T, I, J, S> & {
    const Row<I> M = numRow();
    const Col<J> N = numCol();
    for (size_t r = 0; r < M; ++r)
      for (size_t c = 0; c < N; ++c)
        (*this)(r, c) = x;
    return *this;
  }
  void moveLast(Col<I> j) {
    if (j == N)
      return;
    for (size_t m = 0; m < M; ++m) {
      auto x = (*this)(m, j);
      for (auto n = size_t(j); n < N - 1;) {
        size_t o = n++;
        (*this)(m, o) = (*this)(m, n);
      }
      (*this)(m, N - 1) = x;
    }
  }
  [[nodiscard]] auto deleteCol(size_t c) const -> Matrix<T, I, J, S> {
    Matrix<T, I, J, S> A(M, N - 1);
    for (size_t m = 0; m < M; ++m) {
      A(m, _(0, c)) = (*this)(m, _(0, c));
      A(m, _(c, LinearAlgebra::end)) = (*this)(m, _(c + 1, LinearAlgebra::end));
    }
    return A;
  }
};
template <typename T, IntConvertible I, IntConvertible J, size_t S>
PtrMatrix(Matrix<T, I, J, S> const &) -> PtrMatrix<T, I, J, J>;
template <typename T, IntConvertible I, IntConvertible J, size_t S>
PtrMatrix(Matrix<T, I, J, S> &) -> PtrMatrix<T, I, J, J>;
template <typename T, IntConvertible I, IntConvertible J, size_t S>
MutPtrMatrix(Matrix<T, I, J, S> &) -> MutPtrMatrix<T, I, J, J>;

// template <typename T, IntConvertible I, IntConvertible J>
// PtrMatrix(Matrix<T, I, J, 64> const &) -> PtrMatrix<T, I, J, J>;
// template <typename T, IntConvertible I, IntConvertible J>
// PtrMatrix(Matrix<T, I, J, 64> &) -> PtrMatrix<T, I, J, J>;
// template <typename T, IntConvertible I, IntConvertible J>
// MutPtrMatrix(Matrix<T, I, J, 64> &) -> MutPtrMatrix<T, I, J, J>;

using IntMatrix = Matrix<int64_t>;
static_assert(std::same_as<IntMatrix::eltype, int64_t>);
static_assert(AbstractMatrix<IntMatrix>);
static_assert(AbstractMatrix<Transpose<PtrMatrix<int64_t>>>);
static_assert(std::same_as<eltype_t<Matrix<int64_t>>, int64_t>);

// MutPtrMatrix(IntMatrix &)->MutPtrMatrix<int64_t, unsigned, unsigned,
// unsigned>;

//
// Matrix
//
template <typename T, size_t M, size_t N>
struct Matrix<T, Static<size_t, M>, Static<size_t, N>, 64>
  : ConstMatrixCore<T, Matrix<T, Static<size_t, M>, Static<size_t, N>, 64>>,
    MutMatrixCore<T, Matrix<T, Static<size_t, M>, Static<size_t, N>, 64>> {
  using CBase =
    ConstMatrixCore<T, Matrix<T, Static<size_t, M>, Static<size_t, N>, 64>>;
  using MBase =
    MutMatrixCore<T, Matrix<T, Static<size_t, M>, Static<size_t, N>, 64>>;
  using CBase::antiDiag, MBase::antiDiag;
  using CBase::diag, MBase::diag;
  using CBase::operator(), MBase::operator();
  using CBase::size, CBase::view, CBase::isSquare, CBase::transpose,
    CBase::operator ::LinearAlgebra::PtrMatrix<
      T, Static<size_t, M>, Static<size_t, N>, Static<size_t, N>>,
    MBase::operator ::LinearAlgebra::MutPtrMatrix<
      T, Static<size_t, M>, Static<size_t, N>, Static<size_t, N>>;
  // using eltype = std::remove_cv_t<T>;
  using eltype = std::remove_reference_t<T>;
  // static_assert(M * N == S,
  //               "if specifying non-zero M and N, we should have M*N == S");
  std::array<T, M * N> mem;
  static constexpr auto numRow() -> Row<Static<size_t, M>> { return {}; }
  static constexpr auto numCol() -> Col<Static<size_t, N>> { return {}; }
  static constexpr auto rowStride() -> RowStride<Static<size_t, N>> {
    return {};
  }

  [[nodiscard]] constexpr auto data() -> T * { return mem.data(); }
  [[nodiscard]] constexpr auto data() const -> const T * { return mem.data(); }
};

auto printVectorImpl(llvm::raw_ostream &os, const AbstractVector auto &a)
  -> llvm::raw_ostream & {
  os << "[ ";
  if (size_t M = a.size()) {
    os << a[0];
    for (size_t m = 1; m < M; m++)
      os << ", " << a[m];
  }
  os << " ]";
  return os;
}
template <typename T, IntConvertible I, IntConvertible J>
auto printVector(llvm::raw_ostream &os, PtrVector<T, I, J> a)
  -> llvm::raw_ostream & {
  return printVectorImpl(os, a);
}
template <typename T>
auto printVector(llvm::raw_ostream &os, const llvm::SmallVectorImpl<T> &a)
  -> llvm::raw_ostream & {
  return printVector(os, PtrVector<T>{a.data(), a.size()});
}

template <typename T>
inline auto operator<<(llvm::raw_ostream &os, PtrVector<T> const &A)
  -> llvm::raw_ostream & {
  return printVector(os, A);
}
inline auto operator<<(llvm::raw_ostream &os, const AbstractVector auto &A)
  -> llvm::raw_ostream & {
  return printVector(os, A.view());
}

auto allMatch(const AbstractVector auto &x0, const AbstractVector auto &x1)
  -> bool {
  size_t N = x0.size();
  if (N != x1.size())
    return false;
  for (size_t n = 0; n < N; ++n)
    if (x0(n) != x1(n))
      return false;
  return true;
}

template <typename T, IntConvertible I, IntConvertible J, IntConvertible K>
inline void swap(MutPtrMatrix<T, I, J, K> A, Row<> i, Row<> j) {
  if (i == j)
    return;
  Col<J> N = A.numCol();
  assert((i < A.numRow()) && (j < A.numRow()));
  for (size_t n = 0; n < N; ++n)
    std::swap(A(i, n), A(j, n));
}
template <typename T, IntConvertible I, IntConvertible J, IntConvertible K>
inline void swap(MutPtrMatrix<T, I, J, K> A, Col<> i, Col<> j) {
  if (i == j)
    return;
  Row<I> M = A.numRow();
  assert((i < A.numCol()) && (j < A.numCol()));
  for (size_t m = 0; m < M; ++m)
    std::swap(A(m, i), A(m, j));
}
template <typename T>
inline void swap(llvm::SmallVectorImpl<T> &A, Col<> i, Col<> j) {
  std::swap(A[i], A[j]);
}
template <typename T>
inline void swap(llvm::SmallVectorImpl<T> &A, Row<> i, Row<> j) {
  std::swap(A[i], A[j]);
}

template <int Bits, class T>
constexpr bool is_uint_v =
  sizeof(T) == (Bits / 8) && std::is_integral_v<T> && !std::is_signed_v<T>;

template <class T>
constexpr auto zeroUpper(T x) -> T
requires is_uint_v<16, T>
{
  return x & 0x00ff;
}
template <class T>
constexpr auto zeroLower(T x) -> T
requires is_uint_v<16, T>
{
  return x & 0xff00;
}
template <class T>
constexpr auto upperHalf(T x) -> T
requires is_uint_v<16, T>
{
  return x >> 8;
}

template <class T>
constexpr auto zeroUpper(T x) -> T
requires is_uint_v<32, T>
{
  return x & 0x0000ffff;
}
template <class T>
constexpr auto zeroLower(T x) -> T
requires is_uint_v<32, T>
{
  return x & 0xffff0000;
}
template <class T>
constexpr auto upperHalf(T x) -> T
requires is_uint_v<32, T>
{
  return x >> 16;
}
template <class T>
constexpr auto zeroUpper(T x) -> T
requires is_uint_v<64, T>
{
  return x & 0x00000000ffffffff;
}
template <class T>
constexpr auto zeroLower(T x) -> T
requires is_uint_v<64, T>
{
  return x & 0xffffffff00000000;
}
template <class T>
constexpr auto upperHalf(T x) -> T
requires is_uint_v<64, T>
{
  return x >> 32;
}

template <typename T>
inline auto findMax(llvm::ArrayRef<T> x) -> std::pair<size_t, T> {
  size_t i = 0;
  T max = std::numeric_limits<T>::min();
  for (size_t j = 0; j < x.size(); ++j) {
    T xj = x[j];
    if (max < xj) {
      max = xj;
      i = j;
    }
  }
  return std::make_pair(i, max);
}

template <typename T>
concept TriviallyCopyable = std::is_trivially_copyable_v<T>;

template <typename T>
concept TriviallyCopyableVectorOrScalar =
  std::is_trivially_copyable_v<T> && VectorOrScalar<T>;
template <typename T>
concept TriviallyCopyableMatrixOrScalar =
  std::is_trivially_copyable_v<T> && MatrixOrScalar<T>;

static_assert(std::copy_constructible<PtrMatrix<int64_t>>);
// static_assert(std::is_trivially_copyable_v<MutPtrMatrix<int64_t>>);
static_assert(std::is_trivially_copyable_v<PtrMatrix<int64_t>>);
static_assert(TriviallyCopyableMatrixOrScalar<PtrMatrix<int64_t>>);
static_assert(TriviallyCopyableMatrixOrScalar<int>);
static_assert(TriviallyCopyable<Mul>);
static_assert(TriviallyCopyableMatrixOrScalar<
              ElementwiseMatrixBinaryOp<Mul, PtrMatrix<int64_t>, int>>);
static_assert(TriviallyCopyableMatrixOrScalar<
              MatMatMul<PtrMatrix<int64_t>, PtrMatrix<int64_t>>>);

template <TriviallyCopyable OP, TriviallyCopyableVectorOrScalar A,
          TriviallyCopyableVectorOrScalar B>
ElementwiseVectorBinaryOp(OP, A, B) -> ElementwiseVectorBinaryOp<OP, A, B>;
template <TriviallyCopyable OP, TriviallyCopyableMatrixOrScalar A,
          TriviallyCopyableMatrixOrScalar B>
ElementwiseMatrixBinaryOp(OP, A, B) -> ElementwiseMatrixBinaryOp<OP, A, B>;

inline constexpr auto view(const Scalar auto &x) { return x; }
inline constexpr auto view(const AbstractVector auto &x) { return x.view(); }
inline constexpr auto view(const AbstractMatrixCore auto &x) {
  return x.view();
}

constexpr auto bin2(std::integral auto x) { return (x * (x - 1)) >> 1; }

template <typename T>
auto printMatrix(llvm::raw_ostream &os, PtrMatrix<T> A) -> llvm::raw_ostream & {
  // llvm::raw_ostream &printMatrix(llvm::raw_ostream &os, T const &A) {
  auto [m, n] = A.size();
  if (!m)
    return os << "[ ]";
  for (size_t i = 0; i < m; i++) {
    if (i)
      os << "  ";
    else
      os << "\n[ ";
    if (size_t(n)) {
      for (size_t j = 0; j < n - 1; j++) {
        auto Aij = A(i, j);
        if (Aij >= 0)
          os << " ";
        os << Aij << " ";
      }
    }
    if (n) {
      auto Aij = A(i, n - 1);
      if (Aij >= 0)
        os << " ";
      os << Aij;
    }
    if (i != m - 1)
      os << "\n";
  }
  os << " ]";
  return os;
}

template <typename T, IntConvertible I> struct SmallSparseMatrix {
  // non-zeros
  [[no_unique_address]] llvm::SmallVector<T> nonZeros{};
  // masks, the upper 8 bits give the number of elements in previous rows
  // the remaining 24 bits are a mask indicating non-zeros within this row
  static constexpr size_t maxElemPerRow = 24;
  [[no_unique_address]] llvm::SmallVector<uint32_t> rows;
  [[no_unique_address]] I col;
  [[nodiscard]] constexpr auto numRow() const -> Row<> { return rows.size(); }
  [[nodiscard]] constexpr auto numCol() const -> Col<I> { return col; }
  SmallSparseMatrix(Row<> numRows, Col<I> numCols)
    : rows{llvm::SmallVector<uint32_t>(size_t(numRows))}, col{numCols} {
    assert(size_t(col) <= maxElemPerRow);
  }
  auto get(Row<> i, Col<> j) const -> T {
    assert(j < col);
    uint32_t r(rows[size_t(i)]);
    uint32_t jshift = uint32_t(1) << size_t(j);
    if (r & (jshift)) {
      // offset from previous rows
      uint32_t prevRowOffset = r >> maxElemPerRow;
      uint32_t rowOffset = std::popcount(r & (jshift - 1));
      return nonZeros[rowOffset + prevRowOffset];
    } else {
      return 0;
    }
  }
  constexpr auto operator()(size_t i, size_t j) const -> T {
    return get(toRow(i), col(j));
  }
  void insert(T x, Row<> i, Col<> j) {
    assert(j < col);
    llvm::errs() << "inserting " << x << " at " << size_t(i) << ", "
                 << size_t(j) << "; rows.size() = " << rows.size() << "\n";
    uint32_t r{rows[size_t(i)]};
    uint32_t jshift = uint32_t(1) << size_t(j);
    // offset from previous rows
    uint32_t prevRowOffset = r >> maxElemPerRow;
    uint32_t rowOffset = std::popcount(r & (jshift - 1));
    size_t k = rowOffset + prevRowOffset;
    if (r & jshift) {
      nonZeros[k] = std::move(x);
    } else {
      nonZeros.insert(nonZeros.begin() + k, std::move(x));
      rows[size_t(i)] = r | jshift;
      for (size_t k = size_t(i) + 1; k < rows.size(); ++k)
        rows[k] += uint32_t(1) << maxElemPerRow;
    }
  }

  struct Reference {
    [[no_unique_address]] SmallSparseMatrix<T> *A;
    [[no_unique_address]] size_t i, j;
    operator T() const { return A->get(Row<>(i), Col<>(j)); }
    void operator=(T x) {
      A->insert(std::move(x), Row<>(i), Col<>(j));
      return;
    }
  };
  auto operator()(size_t i, size_t j) -> Reference {
    return Reference{this, i, j};
  }
  operator Matrix<T>() {
    Matrix<T> A(Row<>{numRow()}, Col<>{numCol()});
    assert(numRow() == A.numRow());
    assert(numCol() == A.numCol());
    size_t k = 0;
    for (size_t i = 0; i < numRow(); ++i) {
      uint32_t m = rows[i] & 0x00ffffff;
      size_t j = 0;
      while (m) {
        uint32_t tz = std::countr_zero(m);
        m >>= tz + 1;
        j += tz;
        A(i, j++) = nonZeros[k++];
      }
    }
    assert(k == nonZeros.size());
    return A;
  }
};

template <typename T>
inline auto operator<<(llvm::raw_ostream &os, SmallSparseMatrix<T> const &A)
  -> llvm::raw_ostream & {
  size_t k = 0;
  os << "[ ";
  for (size_t i = 0; i < A.numRow(); ++i) {
    if (i)
      os << "  ";
    uint32_t m = A.rows[i] & 0x00ffffff;
    size_t j = 0;
    while (m) {
      if (j)
        os << " ";
      uint32_t tz = std::countr_zero(m);
      m >>= (tz + 1);
      j += (tz + 1);
      while (tz--)
        os << " 0 ";
      const T &x = A.nonZeros[k++];
      if (x >= 0)
        os << " ";
      os << x;
    }
    for (; j < A.numCol(); ++j)
      os << "  0";
    os << "\n";
  }
  os << " ]";
  assert(k == A.nonZeros.size());
  return os;
}
template <typename T>
inline auto operator<<(llvm::raw_ostream &os, PtrMatrix<T> A)
  -> llvm::raw_ostream & {
  return printMatrix(os, A);
}
template <AbstractMatrix T>
inline auto operator<<(llvm::raw_ostream &os, const T &A)
  -> llvm::raw_ostream & {
  Matrix<std::remove_const_t<typename T::eltype>> B{A};
  return printMatrix(os, PtrMatrix<typename T::eltype>(B));
}

constexpr auto operator-(const AbstractVector auto &a) {
  auto AA{a.view()};
  return ElementwiseUnaryOp<Sub, decltype(AA)>{.op = Sub{}, .a = AA};
}
constexpr auto operator-(const AbstractMatrix auto &a) {
  auto AA{a.view()};
  return ElementwiseUnaryOp<Sub, decltype(AA)>{.op = Sub{}, .a = AA};
}
static_assert(AbstractMatrix<ElementwiseUnaryOp<Sub, PtrMatrix<int64_t>>>);
static_assert(AbstractMatrix<SquareMatrix<int64_t>>);

constexpr auto operator+(const AbstractMatrix auto &a, const auto &b) {
  return ElementwiseMatrixBinaryOp(Add{}, view(a), view(b));
}
constexpr auto operator+(const AbstractVector auto &a, const auto &b) {
  return ElementwiseVectorBinaryOp(Add{}, view(a), view(b));
}
constexpr auto operator+(Scalar auto a, const AbstractMatrix auto &b) {
  return ElementwiseMatrixBinaryOp(Add{}, view(a), view(b));
}
constexpr auto operator+(Scalar auto a, const AbstractVector auto &b) {
  return ElementwiseVectorBinaryOp(Add{}, view(a), view(b));
}
constexpr auto operator-(const AbstractMatrix auto &a, const auto &b) {
  return ElementwiseMatrixBinaryOp(Sub{}, view(a), view(b));
}
constexpr auto operator-(const AbstractVector auto &a, const auto &b) {
  return ElementwiseVectorBinaryOp(Sub{}, view(a), view(b));
}
constexpr auto operator-(Scalar auto a, const AbstractMatrix auto &b) {
  return ElementwiseMatrixBinaryOp(Sub{}, view(a), view(b));
}
constexpr auto operator-(Scalar auto a, const AbstractVector auto &b) {
  return ElementwiseVectorBinaryOp(Sub{}, view(a), view(b));
}
constexpr auto operator/(const AbstractMatrix auto &a, const auto &b) {
  return ElementwiseMatrixBinaryOp(Div{}, view(a), view(b));
}
constexpr auto operator/(const AbstractVector auto &a, const auto &b) {
  return ElementwiseVectorBinaryOp(Div{}, view(a), view(b));
}
constexpr auto operator/(Scalar auto a, const AbstractMatrix auto &b) {
  return ElementwiseMatrixBinaryOp(Div{}, view(a), view(b));
}
constexpr auto operator/(Scalar auto a, const AbstractVector auto &b) {
  return ElementwiseVectorBinaryOp(Div{}, view(a), view(b));
}
constexpr auto operator*(const AbstractMatrix auto &a,
                         const AbstractMatrix auto &b) {
  auto AA{a.view()};
  auto BB{b.view()};
  assert(size_t(AA.numCol()) == size_t(BB.numRow()));
  return MatMatMul<decltype(AA), decltype(BB)>{.a = AA, .b = BB};
}
constexpr auto operator*(const AbstractMatrix auto &a,
                         const AbstractVector auto &b) {
  auto AA{a.view()};
  auto BB{b.view()};
  assert(size_t(AA.numCol()) == BB.size());
  return MatVecMul<decltype(AA), decltype(BB)>{.a = AA, .b = BB};
}
constexpr auto operator*(const AbstractMatrix auto &a, std::integral auto b) {
  return ElementwiseMatrixBinaryOp(Mul{}, view(a), view(b));
}
constexpr auto operator*(const AbstractVector auto &a,
                         const AbstractVector auto &b) {
  return ElementwiseVectorBinaryOp(Mul{}, view(a), view(b));
}
constexpr auto operator*(const AbstractVector auto &a, std::integral auto b) {
  return ElementwiseVectorBinaryOp(Mul{}, view(a), view(b));
}
constexpr auto operator*(Scalar auto a, const AbstractMatrix auto &b) {
  return ElementwiseMatrixBinaryOp(Mul{}, view(a), view(b));
}
constexpr auto operator*(Scalar auto a, const AbstractVector auto &b) {
  return ElementwiseVectorBinaryOp(Mul{}, view(a), view(b));
}

// constexpr auto operator*(AbstractMatrix auto &A, AbstractVector auto &x) {
//     auto AA{A.view()};
//     auto xx{x.view()};
//     return MatMul<decltype(AA), decltype(xx)>{.a = AA, .b = xx};
// }

template <AbstractVector V>
constexpr auto operator*(const Transpose<V> &a, const AbstractVector auto &b) {
  typename V::eltype s = 0;
  for (size_t i = 0; i < b.size(); ++i)
    s += a.a(i) * b(i);
  return s;
}

static_assert(AbstractVector<Vector<int64_t>>);
static_assert(AbstractVector<const Vector<int64_t>>);
static_assert(AbstractVector<Vector<int64_t> &>);
static_assert(AbstractMatrix<IntMatrix>);
static_assert(AbstractMatrix<IntMatrix &>);

static_assert(
  std::copyable<Matrix<int64_t, Static<size_t, 4>, Static<size_t, 4>>>);
static_assert(std::copyable<Matrix<int64_t, Static<size_t, 4>, size_t, 64>>);
static_assert(std::copyable<Matrix<int64_t, Static<size_t, 4>, size_t>>);
static_assert(std::copyable<Matrix<int64_t, size_t, Static<size_t, 4>>>);
static_assert(std::copyable<Matrix<int64_t, size_t, size_t>>);
static_assert(std::copyable<SquareMatrix<int64_t>>);

static_assert(
  DerivedMatrix<Matrix<int64_t, Static<size_t, 4>, Static<size_t, 4>>>);
static_assert(DerivedMatrix<Matrix<int64_t, Static<size_t, 4>, size_t>>);
static_assert(DerivedMatrix<Matrix<int64_t, size_t, Static<size_t, 4>>>);
static_assert(DerivedMatrix<Matrix<int64_t, size_t, size_t>>);
static_assert(DerivedMatrix<IntMatrix>);
static_assert(DerivedMatrix<IntMatrix>);
static_assert(DerivedMatrix<IntMatrix>);

static_assert(std::is_same_v<SquareMatrix<int64_t>::eltype, int64_t>);
static_assert(std::is_same_v<IntMatrix::eltype, int64_t>);

template <typename T, typename I> struct SliceView {
  using eltype = T;
  [[no_unique_address]] MutPtrVector<T> a;
  [[no_unique_address]] llvm::ArrayRef<I> i;
  struct Iterator {
    [[no_unique_address]] MutPtrVector<T> a;
    [[no_unique_address]] llvm::ArrayRef<I> i;
    [[no_unique_address]] size_t j;
    auto operator==(const Iterator &k) const -> bool { return j == k.j; }
    auto operator++() -> Iterator & {
      ++j;
      return *this;
    }
    auto operator*() -> T & { return a[i[j]]; }
    auto operator*() const -> const T & { return a[i[j]]; }
    auto operator->() -> T * { return &a[i[j]]; }
    auto operator->() const -> const T * { return &a[i[j]]; }
  };
  constexpr auto begin() -> Iterator { return Iterator{a, i, 0}; }
  constexpr auto end() -> Iterator { return Iterator{a, i, i.size()}; }
  auto operator()(size_t j) -> T & { return a[i[j]]; }
  auto operator()(size_t j) const -> const T & { return a[i[j]]; }
  [[nodiscard]] constexpr auto size() const -> size_t { return i.size(); }
  constexpr auto view() -> SliceView<T, I> { return *this; }
};

static_assert(AbstractVector<SliceView<int64_t, unsigned>>);
} // namespace LinearAlgebra

// exports:
// NOLINTNEXTLINE(bugprone-reserved-identifier)
using LinearAlgebra::_;
using LinearAlgebra::AbstractVector, LinearAlgebra::AbstractMatrix,
  LinearAlgebra::PtrVector, LinearAlgebra::MutPtrVector, LinearAlgebra::Vector,
  LinearAlgebra::Matrix, LinearAlgebra::SquareMatrix, LinearAlgebra::IntMatrix,
  LinearAlgebra::PtrMatrix, LinearAlgebra::MutPtrMatrix, LinearAlgebra::AxisInt,
  LinearAlgebra::AxisInt, LinearAlgebra::SmallSparseMatrix,
  LinearAlgebra::SquareMatrix, LinearAlgebra::MutSquarePtrMatrix,
  LinearAlgebra::Range, LinearAlgebra::begin, LinearAlgebra::end,
  LinearAlgebra::swap, LinearAlgebra::SquarePtrMatrix, LinearAlgebra::Row,
  LinearAlgebra::RowStride, LinearAlgebra::Col;
