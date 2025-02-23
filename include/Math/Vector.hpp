#pragma once
#include "Math/Indexing.hpp"
#include "TypePromotion.hpp"
#include "Utilities/Valid.hpp"
#include <llvm/ADT/ArrayRef.h>

namespace LinearAlgebra {
template <typename T>
concept AbstractVector =
  HasEltype<T> && requires(T t, size_t i) {
                    { t[i] } -> std::convertible_to<eltype_t<T>>;
                    { t.size() } -> std::convertible_to<size_t>;
                    { t.view() };
                    // {
                    //     std::remove_reference_t<T>::canResize
                    //     } -> std::same_as<const bool &>;
                    // {t.extendOrAssertSize(i)};
                  };

template <typename T> struct PtrVector {
  static_assert(!std::is_const_v<T>, "const T is redundant");
  using eltype = T;
  [[no_unique_address]] NotNull<const T> mem;
  [[no_unique_address]] const size_t N;
  auto operator==(AbstractVector auto &x) -> bool {
    if (N != x.size()) return false;
    for (size_t n = 0; n < N; ++n)
      if (mem[n] != x(n)) return false;
    return true;
  }
  [[nodiscard]] constexpr auto front() const -> const T & { return mem[0]; }
  [[gnu::flatten]] constexpr auto operator[](const ScalarIndex auto i) const
    -> const T & {
#ifndef NDEBUG
    checkIndex(size_t(N), i);
#endif
    return mem[canonicalize(i, N)];
  }
  [[gnu::flatten]] constexpr auto operator[](Range<size_t, size_t> i) const
    -> PtrVector<T> {
    assert(i.b <= i.e);
    assert(i.e <= N);
    return PtrVector<T>{mem + i.b, i.e - i.b};
  }
  template <typename F, typename L>
  [[gnu::flatten]] constexpr auto operator[](Range<F, L> i) const
    -> PtrVector<T> {
    return (*this)[canonicalizeRange(i, N)];
  }
  [[gnu::returns_nonnull]] [[nodiscard]] constexpr auto begin() const
    -> const T * {
    return mem;
  }
  [[gnu::returns_nonnull]] [[nodiscard]] constexpr auto end() const
    -> const T * {
    return mem + N;
  }
  [[nodiscard]] constexpr auto rbegin() const {
    return std::reverse_iterator(end());
  }
  [[nodiscard]] constexpr auto rend() const {
    return std::reverse_iterator(begin());
  }
  [[nodiscard]] constexpr auto size() const -> size_t { return N; }
  constexpr operator llvm::ArrayRef<T>() const {
    return llvm::ArrayRef<T>{mem, N};
  }
  // llvm::ArrayRef<T> arrayref() const { return llvm::ArrayRef<T>(ptr, M); }
  auto operator==(const PtrVector<T> x) const -> bool {
    return llvm::ArrayRef<T>(*this) == llvm::ArrayRef<T>(x);
  }
  auto operator==(const llvm::ArrayRef<std::remove_const_t<T>> x) const
    -> bool {
    return llvm::ArrayRef<std::remove_const_t<T>>(*this) == x;
  }
  [[nodiscard]] constexpr auto view() const -> PtrVector<T> { return *this; };
#ifndef NDEBUG
  void extendOrAssertSize(size_t M) const { assert(M == N); }
#else
  // two defs to avoid unused parameter compiler warning in release builds
  static constexpr void extendOrAssertSize(size_t) {}
#endif
  constexpr PtrVector(NotNull<const T> pt, size_t NN) : mem(pt), N(NN) {}
  PtrVector(llvm::ArrayRef<T> x) : mem(x.data()), N(x.size()) {}
};
template <typename T> struct MutPtrVector {
  static_assert(!std::is_const_v<T>, "T shouldn't be const");
  using eltype = T;
  // using eltype = std::remove_const_t<T>;
  [[no_unique_address]] NotNull<T> mem;
  [[no_unique_address]] const size_t N;
  [[gnu::flatten]] constexpr auto operator[](const ScalarIndex auto i) -> T & {
#ifndef NDEBUG
    checkIndex(size_t(N), i);
#endif
    return mem[canonicalize(i, N)];
  }
  [[gnu::flatten]] constexpr auto operator[](const ScalarIndex auto i) const
    -> const T & {
#ifndef NDEBUG
    checkIndex(size_t(N), i);
#endif
    return mem[canonicalize(i, N)];
  }
  [[nodiscard]] auto front() -> T & {
    assert(N > 0);
    return mem[0];
  }
  [[nodiscard]] auto back() -> T & {
    assert(N > 0);
    return mem[N - 1];
  }
  [[nodiscard]] auto front() const -> const T & {
    assert(N > 0);
    return mem[0];
  }
  [[nodiscard]] auto back() const -> const T & {
    assert(N > 0);
    return mem[N - 1];
  }
  [[nodiscard]] constexpr auto isEmpty() const -> bool { return N == 0; }
  // copy constructor
  constexpr MutPtrVector() = default;
  constexpr MutPtrVector(const MutPtrVector<T> &x) = default;
  constexpr MutPtrVector(llvm::MutableArrayRef<T> x)
    : mem(x.data()), N(x.size()) {}
  constexpr MutPtrVector(T *pt, size_t NN) : mem(pt), N(NN) {}
  constexpr auto operator[](Range<size_t, size_t> i) -> MutPtrVector<T> {
    assert(i.b <= i.e);
    assert(i.e <= N);
    return MutPtrVector<T>{mem + i.b, i.e - i.b};
  }
  constexpr auto operator[](Range<size_t, size_t> i) const -> PtrVector<T> {
    assert(i.b <= i.e);
    assert(i.e <= N);
    return PtrVector<T>{mem + i.b, i.e - i.b};
  }
  template <typename F, typename L>
  constexpr auto operator[](Range<F, L> i) -> MutPtrVector<T> {
    return (*this)[canonicalizeRange(i, N)];
  }
  template <typename F, typename L>
  constexpr auto operator[](Range<F, L> i) const -> PtrVector<T> {
    return (*this)[canonicalizeRange(i, N)];
  }
  [[gnu::returns_nonnull]] constexpr auto begin() -> T * { return mem; }
  [[gnu::returns_nonnull]] constexpr auto end() -> T * { return mem + N; }
  [[gnu::returns_nonnull]] [[nodiscard]] constexpr auto begin() const
    -> const T * {
    return mem;
  }
  [[gnu::returns_nonnull]] [[nodiscard]] constexpr auto end() const
    -> const T * {
    return mem + N;
  }
  [[nodiscard]] constexpr auto size() const -> size_t { return N; }
  constexpr operator PtrVector<T>() const { return PtrVector<T>{mem, N}; }
  constexpr operator llvm::ArrayRef<T>() const {
    return llvm::ArrayRef<T>{mem, N};
  }
  constexpr operator llvm::MutableArrayRef<T>() {
    return llvm::MutableArrayRef<T>{mem, N};
  }
  // llvm::ArrayRef<T> arrayref() const { return llvm::ArrayRef<T>(ptr, M); }
  auto operator==(const MutPtrVector<T> x) const -> bool {
    return llvm::ArrayRef<T>(*this) == llvm::ArrayRef<T>(x);
  }
  auto operator==(const PtrVector<T> x) const -> bool {
    return llvm::ArrayRef<T>(*this) == llvm::ArrayRef<T>(x);
  }
  auto operator==(const llvm::ArrayRef<T> x) const -> bool {
    return llvm::ArrayRef<T>(*this) == x;
  }
  [[nodiscard]] constexpr auto view() const -> PtrVector<T> { return *this; };
  [[gnu::flatten]] auto operator=(PtrVector<T> x) -> MutPtrVector<T> {
    return copyto(*this, x);
  }
  [[gnu::flatten]] auto operator=(MutPtrVector<T> x) -> MutPtrVector<T> {
    return copyto(*this, x);
  }
  [[gnu::flatten]] auto operator=(const AbstractVector auto &x)
    -> MutPtrVector<T> {
    return copyto(*this, x);
  }
  [[gnu::flatten]] auto operator=(std::integral auto x) -> MutPtrVector<T> {
    for (auto &&y : *this) y = x;
    return *this;
  }
  [[gnu::flatten]] auto operator+=(const AbstractVector auto &x)
    -> MutPtrVector<T> {
    assert(N == x.size());
    for (size_t i = 0; i < N; ++i) mem[i] += x[i];
    return *this;
  }
  [[gnu::flatten]] auto operator-=(const AbstractVector auto &x)
    -> MutPtrVector<T> {
    assert(N == x.size());
    for (size_t i = 0; i < N; ++i) mem[i] -= x[i];
    return *this;
  }
  [[gnu::flatten]] auto operator*=(const AbstractVector auto &x)
    -> MutPtrVector<T> {
    assert(N == x.size());
    for (size_t i = 0; i < N; ++i) mem[i] *= x[i];
    return *this;
  }
  [[gnu::flatten]] auto operator/=(const AbstractVector auto &x)
    -> MutPtrVector<T> {
    assert(N == x.size());
    for (size_t i = 0; i < N; ++i) mem[i] /= x[i];
    return *this;
  }
  [[gnu::flatten]] auto operator+=(const std::integral auto x)
    -> MutPtrVector<T> {
    for (size_t i = 0; i < N; ++i) mem[i] += x;
    return *this;
  }
  [[gnu::flatten]] auto operator-=(const std::integral auto x)
    -> MutPtrVector<T> {
    for (size_t i = 0; i < N; ++i) mem[i] -= x;
    return *this;
  }
  [[gnu::flatten]] auto operator*=(const std::integral auto x)
    -> MutPtrVector<T> {
    for (size_t i = 0; i < N; ++i) mem[i] *= x;
    return *this;
  }
  [[gnu::flatten]] auto operator/=(const std::integral auto x)
    -> MutPtrVector<T> {
    for (size_t i = 0; i < N; ++i) mem[i] /= x;
    return *this;
  }
#ifndef NDEBUG
  void extendOrAssertSize(size_t M) const { assert(M == N); }
#else
  static constexpr void extendOrAssertSize(size_t) {}
#endif
};

template <typename T> struct Vector {
  using eltype = T;
  [[no_unique_address]] llvm::SmallVector<T, 16> data;

  Vector(int N) : data(llvm::SmallVector<T, 16>(N)){};
  Vector(size_t N = 0) : data(llvm::SmallVector<T, 16>(N)){};
  Vector(llvm::SmallVector<T, 16> A) : data(std::move(A)){};

  [[gnu::flatten]] constexpr auto operator[](const ScalarIndex auto i) -> T & {
    return data[canonicalize(i, data.size())];
  }
  [[gnu::flatten]] constexpr auto operator[](const ScalarIndex auto i) const
    -> const T & {
    return data[canonicalize(i, data.size())];
  }
  [[gnu::flatten]] constexpr auto operator[](Range<size_t, size_t> i)
    -> MutPtrVector<T> {
    assert(i.b <= i.e);
    assert(i.e <= data.size());
    return MutPtrVector<T>{data.data() + i.b, i.e - i.b};
  }
  [[gnu::flatten]] constexpr auto operator[](Range<size_t, size_t> i) const
    -> PtrVector<T> {
    assert(i.b <= i.e);
    assert(i.e <= data.size());
    return PtrVector<T>{data.data() + i.b, i.e - i.b};
  }
  template <typename F, typename L>
  [[gnu::flatten]] constexpr auto operator[](Range<F, L> i) -> MutPtrVector<T> {
    return (*this)[canonicalizeRange(i, data.size())];
  }
  template <typename F, typename L>
  [[gnu::flatten]] constexpr auto operator[](Range<F, L> i) const
    -> PtrVector<T> {
    return (*this)[canonicalizeRange(i, data.size())];
  }
  [[gnu::flatten]] constexpr auto operator[](size_t i) -> T & {
    return data[i];
  }
  [[gnu::flatten]] constexpr auto operator[](size_t i) const -> const T & {
    return data[i];
  }
  // bool operator==(Vector<T, 0> x0) const { return allMatch(*this, x0); }
  constexpr auto begin() { return data.begin(); }
  constexpr auto end() { return data.end(); }
  [[nodiscard]] constexpr auto begin() const { return data.begin(); }
  [[nodiscard]] constexpr auto end() const { return data.end(); }
  [[nodiscard]] constexpr auto size() const -> size_t { return data.size(); }
  [[nodiscard]] constexpr auto view() const -> PtrVector<T> {
    return PtrVector<T>{data.data(), data.size()};
  };
  template <typename A> void push_back(A &&x) {
    data.push_back(std::forward<A>(x));
  }
  template <typename... A> void emplace_back(A &&...x) {
    data.emplace_back(std::forward<A>(x)...);
  }
  Vector(const AbstractVector auto &x) : data(llvm::SmallVector<T, 16>{}) {
    const size_t N = x.size();
    data.resize_for_overwrite(N);
    for (size_t n = 0; n < N; ++n) data[n] = x[n];
  }
  void resize(size_t N) { data.resize(N); }
  void resizeForOverwrite(size_t N) { data.resize_for_overwrite(N); }

  constexpr operator MutPtrVector<T>() {
    return MutPtrVector<T>{data.data(), data.size()};
  }
  constexpr operator PtrVector<T>() const {
    return PtrVector<T>{data.data(), data.size()};
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
    for (auto &&y : data) y += x;
    return *this;
  }
  auto operator-=(const std::integral auto x) -> Vector<T> & {
    for (auto &&y : data) y -= x;
    return *this;
  }
  auto operator*=(const std::integral auto x) -> Vector<T> & {
    for (auto &&y : data) y *= x;
    return *this;
  }
  auto operator/=(const std::integral auto x) -> Vector<T> & {
    for (auto &&y : data) y /= x;
    return *this;
  }
  template <typename... Ts> Vector(Ts... inputs) : data{inputs...} {}
  void clear() { data.clear(); }
#ifndef NDEBUG
  void extendOrAssertSize(size_t N) const { assert(N == data.size()); }
#else
  constexpr void extendOrAssertSize(size_t) const {}
#endif
  void extendOrAssertSize(size_t N) {
    if (N != data.size()) data.resize_for_overwrite(N);
  }
  auto operator==(const Vector<T> &x) const -> bool {
    return llvm::ArrayRef<T>(*this) == llvm::ArrayRef<T>(x);
  }
  void pushBack(T x) { data.push_back(std::move(x)); }
};

static_assert(std::copyable<Vector<intptr_t>>);
static_assert(AbstractVector<Vector<int64_t>>);
static_assert(!AbstractVector<int64_t>);

template <typename T> struct StridedVector {
  static_assert(!std::is_const_v<T>, "const T is redundant");
  using eltype = T;
  [[no_unique_address]] const T *d;
  [[no_unique_address]] size_t N;
  [[no_unique_address]] RowStride x;
  struct StridedIterator {
    using value_type = const T;
    [[no_unique_address]] const T *d;
    [[no_unique_address]] size_t xx;
    constexpr auto operator++() -> StridedIterator & {
      d += xx;
      return *this;
    }
    constexpr auto operator--() -> StridedIterator & {
      d -= xx;
      return *this;
    }
    constexpr auto operator++(int) {
      auto tmp = *this;
      d += xx;
      return tmp;
    }
    constexpr auto operator--(int) {
      auto tmp = *this;
      d -= xx;
      return tmp;
    }
    auto operator[](ptrdiff_t y) const -> const T & { return d[y * xx]; }
    auto operator-(StridedIterator y) const -> ptrdiff_t {
      return (d - y.d) / xx;
    }
    auto operator+(ptrdiff_t y) const -> StridedIterator {
      return {d + y * xx, xx};
    }
    auto operator-(ptrdiff_t y) const -> StridedIterator {
      return {d + y * xx, xx};
    }
    auto operator+=(ptrdiff_t y) -> StridedIterator & {
      d += y * xx;
      return *this;
    }
    auto operator-=(ptrdiff_t y) -> StridedIterator & {
      d -= y * xx;
      return *this;
    }
    constexpr auto operator*() const -> const T & { return *d; }
    constexpr auto operator->() const -> const T * { return d; }
    constexpr auto operator==(const StridedIterator y) const -> bool {
      return d == y.d;
    }
    constexpr auto operator!=(const StridedIterator y) const -> bool {
      return d != y.d;
    }
    constexpr auto operator>(const StridedIterator y) const -> bool {
      return d > y.d;
    }
    constexpr auto operator<(const StridedIterator y) const -> bool {
      return d < y.d;
    }
    constexpr auto operator>=(const StridedIterator y) const -> bool {
      return d >= y.d;
    }
    constexpr auto operator<=(const StridedIterator y) const -> bool {
      return d <= y.d;
    }
    friend auto operator+(ptrdiff_t y,
                          typename StridedVector<T>::StridedIterator a) ->
      typename StridedVector<T>::StridedIterator {
      return {a.d + y * a.xx, a.xx};
    }
  };
  [[nodiscard]] constexpr auto begin() const {
    return StridedIterator{d, size_t(x)};
  }
  [[nodiscard]] constexpr auto end() const {
    return StridedIterator{d + x * N, size_t(x)};
  }
  [[nodiscard]] constexpr auto rbegin() const {
    return std::reverse_iterator(end());
  }
  [[nodiscard]] constexpr auto rend() const {
    return std::reverse_iterator(begin());
  }
  // [[nodiscard]] constexpr auto begin() const {
  //   return std::ranges::stride_view{llvm::ArrayRef<T>{d, N}, x};
  // }
  constexpr auto operator[](size_t i) const -> const T & {
    return d[size_t(x * i)];
  }

  constexpr auto operator[](Range<size_t, size_t> i) const -> StridedVector<T> {
    return StridedVector<T>{.d = d + x * i.b, .N = i.e - i.b, .x = x};
  }
  template <typename F, typename L>
  constexpr auto operator[](Range<F, L> i) const -> StridedVector<T> {
    return (*this)[canonicalizeRange(i, N)];
  }

  [[nodiscard]] constexpr auto size() const -> size_t { return N; }
  auto operator==(StridedVector<T> a) const -> bool {
    if (size() != a.size()) return false;
    for (size_t i = 0; i < size(); ++i)
      if ((*this)[i] != a[i]) return false;
    return true;
  }
  [[nodiscard]] constexpr auto view() const -> StridedVector<T> {
    return *this;
  }
#ifndef NDEBUG
  void extendOrAssertSize(size_t M) const { assert(N == M); }
#else
  static constexpr void extendOrAssertSize(size_t) {}
#endif
};
template <typename T> struct MutStridedVector {
  static_assert(!std::is_const_v<T>, "T should not be const");
  using eltype = T;
  [[no_unique_address]] T *const d;
  [[no_unique_address]] const size_t N;
  [[no_unique_address]] RowStride x;
  struct StridedIterator {
    using value_type = T;

    [[no_unique_address]] T *d;
    [[no_unique_address]] size_t xx;
    constexpr auto operator++() -> StridedIterator & {
      d += xx;
      return *this;
    }
    constexpr auto operator--() -> StridedIterator & {
      d -= xx;
      return *this;
    }
    constexpr auto operator++(int) {
      auto tmp = *this;
      d += xx;
      return tmp;
    }
    constexpr auto operator--(int) {
      auto tmp = *this;
      d -= xx;
      return tmp;
    }
    auto operator[](ptrdiff_t y) const -> T & { return d[y * xx]; }
    auto operator-(StridedIterator y) const -> ptrdiff_t {
      return (d - y.d) / xx;
    }
    auto operator+(ptrdiff_t y) const -> StridedIterator {
      return {d + y * xx, xx};
    }
    auto operator-(ptrdiff_t y) const -> StridedIterator {
      return {d + y * xx, xx};
    }
    auto operator+=(ptrdiff_t y) -> StridedIterator & {
      d += y * xx;
      return *this;
    }
    auto operator-=(ptrdiff_t y) -> StridedIterator & {
      d -= y * xx;
      return *this;
    }
    constexpr auto operator->() const -> T * { return d; }
    constexpr auto operator*() const -> T & { return *d; }
    // constexpr auto operator->() -> T * { return d; }
    // constexpr auto operator*() -> T & { return *d; }
    // constexpr auto operator->() const -> const T * { return d; }
    // constexpr auto operator*() const -> const T & { return *d; }
    constexpr auto operator==(const StridedIterator y) const -> bool {
      return d == y.d;
    }
    constexpr auto operator!=(const StridedIterator y) const -> bool {
      return d != y.d;
    }
    constexpr auto operator>(const StridedIterator y) const -> bool {
      return d > y.d;
    }
    constexpr auto operator<(const StridedIterator y) const -> bool {
      return d < y.d;
    }
    constexpr auto operator>=(const StridedIterator y) const -> bool {
      return d >= y.d;
    }
    constexpr auto operator<=(const StridedIterator y) const -> bool {
      return d <= y.d;
    }
    friend auto operator+(ptrdiff_t y,
                          typename MutStridedVector<T>::StridedIterator a) ->
      typename MutStridedVector<T>::StridedIterator {
      return {a.d + y * a.xx, a.xx};
    }
  };
  // FIXME: if `x` == 0, then it will not iterate!
  constexpr auto begin() { return StridedIterator{d, size_t(x)}; }
  constexpr auto end() { return StridedIterator{d + x * N, size_t(x)}; }
  [[nodiscard]] constexpr auto begin() const {
    return StridedIterator{d, size_t(x)};
  }
  [[nodiscard]] constexpr auto end() const {
    return StridedIterator{d + x * N, size_t(x)};
  }
  constexpr auto rbegin() { return std::reverse_iterator(end()); }
  constexpr auto rend() { return std::reverse_iterator(begin()); }
  [[nodiscard]] constexpr auto rbegin() const {
    return std::reverse_iterator(end());
  }
  [[nodiscard]] constexpr auto rend() const {
    return std::reverse_iterator(begin());
  }
  constexpr auto operator[](size_t i) -> T & { return d[size_t(x * i)]; }
  constexpr auto operator[](size_t i) const -> const T & {
    return d[size_t(x * i)];
  }
  constexpr auto operator[](Range<size_t, size_t> i) -> MutStridedVector<T> {
    return MutStridedVector<T>{.d = d + x * i.b, .N = i.e - i.b, .x = x};
  }
  constexpr auto operator[](Range<size_t, size_t> i) const -> StridedVector<T> {
    return StridedVector<T>{.d = d + x * i.b, .N = i.e - i.b, .x = x};
  }
  template <typename F, typename L>
  constexpr auto operator[](Range<F, L> i) -> MutStridedVector<T> {
    return (*this)[canonicalizeRange(i, N)];
  }
  template <typename F, typename L>
  constexpr auto operator[](Range<F, L> i) const -> StridedVector<T> {
    return (*this)[canonicalizeRange(i, N)];
  }

  [[nodiscard]] constexpr auto size() const -> size_t { return N; }
  // bool operator==(StridedVector<T> x) const {
  //     if (size() != x.size())
  //         return false;
  //     for (size_t i = 0; i < size(); ++i) {
  //         if ((*this)[i] != x[i])
  //             return false;
  //     }
  //     return true;
  // }
  constexpr operator StridedVector<T>() {
    const T *const p = d;
    return StridedVector<T>{.d = p, .N = N, .x = x};
  }
  [[nodiscard]] constexpr auto view() const -> StridedVector<T> {
    return StridedVector<T>{.d = d, .N = N, .x = x};
  }
  [[gnu::flatten]] auto operator=(const T &y) -> MutStridedVector<T> & {
    for (size_t i = 0; i < N; ++i) d[size_t(x * i)] = y;
    return *this;
  }
  [[gnu::flatten]] auto operator=(const AbstractVector auto &a)
    -> MutStridedVector<T> & {
    return copyto(*this, a);
  }
  [[gnu::flatten]] auto operator=(const MutStridedVector<T> &a)
    -> MutStridedVector<T> & {
    if (this == &a) return *this;
    return copyto(*this, a);
  }
  [[gnu::flatten]] auto operator+=(T a) -> MutStridedVector<T> & {
    MutStridedVector<T> &self = *this;
    for (size_t i = 0; i < N; ++i) self[i] += a;
    return self;
  }
  [[gnu::flatten]] auto operator+=(const AbstractVector auto &a)
    -> MutStridedVector<T> & {
    const size_t M = a.size();
    MutStridedVector<T> &self = *this;
    assert(M == N);
    for (size_t i = 0; i < M; ++i) self[i] += a[i];
    return self;
  }
  [[gnu::flatten]] auto operator-=(const AbstractVector auto &a)
    -> MutStridedVector<T> & {
    const size_t M = a.size();
    MutStridedVector<T> &self = *this;
    assert(M == N);
    for (size_t i = 0; i < M; ++i) self[i] -= a[i];
    return self;
  }
  [[gnu::flatten]] auto operator*=(const AbstractVector auto &a)
    -> MutStridedVector<T> & {
    const size_t M = a.size();
    MutStridedVector<T> &self = *this;
    assert(M == N);
    for (size_t i = 0; i < M; ++i) self[i] *= a[i];
    return self;
  }
  [[gnu::flatten]] auto operator/=(const AbstractVector auto &a)
    -> MutStridedVector<T> & {
    const size_t M = a.size();
    MutStridedVector<T> &self = *this;
    assert(M == N);
    for (size_t i = 0; i < M; ++i) self[i] /= a[i];
    return self;
  }
#ifndef NDEBUG
  void extendOrAssertSize(size_t M) const { assert(N == M); }
#else
  static constexpr void extendOrAssertSize(size_t) {}
#endif
};
static_assert(
  std::weakly_incrementable<StridedVector<int64_t>::StridedIterator>);
static_assert(
  std::input_or_output_iterator<StridedVector<int64_t>::StridedIterator>);

static_assert(std::indirectly_readable<StridedVector<int64_t>::StridedIterator>,
              "failed indirectly readable");
static_assert(
  std::indirectly_readable<MutStridedVector<int64_t>::StridedIterator>,
  "failed indirectly readable");
static_assert(
  std::output_iterator<MutStridedVector<int64_t>::StridedIterator, int>,
  "failed output iterator");
static_assert(std::forward_iterator<StridedVector<int64_t>::StridedIterator>,
              "failed forward iterator");
static_assert(std::input_iterator<StridedVector<int64_t>::StridedIterator>,
              "failed input iterator");
static_assert(
  std::bidirectional_iterator<StridedVector<int64_t>::StridedIterator>,
  "failed bidirectional iterator");

static_assert(std::totally_ordered<StridedVector<int64_t>::StridedIterator>,
              "failed random access iterator");
static_assert(
  std::random_access_iterator<StridedVector<int64_t>::StridedIterator>,
  "failed random access iterator");
static_assert(AbstractVector<StridedVector<int64_t>>);
static_assert(std::is_trivially_copyable_v<StridedVector<int64_t>>);

} // namespace LinearAlgebra
