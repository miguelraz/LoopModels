#pragma once
#include "./Math.hpp"
#include <cstddef>
#include <cstdint>
#include <llvm/ADT/SmallVector.h>

template <typename T> struct EmptyMatrix {
  using eltype = T;
  static constexpr auto getLinearElement(size_t) -> T { return 0; }
  static constexpr auto begin() -> T * { return nullptr; }
  static constexpr auto end() -> T * { return nullptr; }

  static constexpr auto numRow() { return Row<Static<size_t, 0>>{}; }
  // static constexpr auto numRow() -> Row<Static<size_t, 0>> { return {}; }

  static constexpr auto numCol() -> Col<Static<size_t, 0>> { return {}; }
  static constexpr auto rowStride()
    -> LinearAlgebra::RowStride<Static<size_t, 0>> {
    return {};
  }
  static constexpr auto getConstCol() -> size_t { return 0; }

  static constexpr auto data() -> T * { return nullptr; }
  constexpr auto operator()(size_t, size_t) -> T { return 0; }
  static constexpr auto size()
    -> std::pair<Row<Static<size_t, 0>>, Col<Static<size_t, 0>>> {
    return {};
  }
  static constexpr auto view() -> EmptyMatrix<T> { return EmptyMatrix<T>{}; }
};

static_assert(AbstractMatrix<EmptyMatrix<int64_t>>);

template <typename T>
constexpr auto matmul(EmptyMatrix<T>, PtrMatrix<const T>) -> EmptyMatrix<T> {
  return EmptyMatrix<T>{};
}
template <typename T>
constexpr auto matmul(PtrMatrix<const T>, EmptyMatrix<T>) -> EmptyMatrix<T> {
  return EmptyMatrix<T>{};
}

template <typename T, typename S>
concept MaybeMatrix =
  std::is_same_v<T, Matrix<S>> || std::is_same_v<T, EmptyMatrix<S>>;

template <typename T> struct EmptyVector {
  static constexpr auto size() -> size_t { return 0; };
  static constexpr auto begin() -> T * { return nullptr; }
  static constexpr auto end() -> T * { return nullptr; }
};
