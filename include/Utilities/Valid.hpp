#pragma once
#include "Utilities/Invariant.hpp"
#include <llvm/Support/Casting.h>

// TODO: communicate not-null to the compiler somehow?
template <typename T> struct NotNull {
  NotNull() = delete;
  constexpr NotNull(T &v) : value(&v) {}
  constexpr NotNull(T *v) : value(v) { invariant(value != nullptr); }
  constexpr explicit operator bool() const {
    invariant(value != nullptr);
    return true;
  }
  constexpr operator NotNull<const T>() const {
    invariant(value != nullptr);
    return NotNull<const T>(*value);
  }
  [[gnu::returns_nonnull]] constexpr operator T *() {
    invariant(value != nullptr);
    return value;
  }
  [[gnu::returns_nonnull]] constexpr operator T *() const {
    invariant(value != nullptr);
    return value;
  }
  [[gnu::returns_nonnull]] constexpr auto operator->() -> T * {
    invariant(value != nullptr);
    return value;
  }
  constexpr auto operator*() -> T & {
    invariant(value != nullptr);
    return *value;
  }
  [[gnu::returns_nonnull]] constexpr auto operator->() const -> const T * {
    invariant(value != nullptr);
    return value;
  }
  constexpr auto operator*() const -> const T & {
    invariant(value != nullptr);
    return *value;
  }
  constexpr auto operator[](size_t index) -> T & {
    invariant(value != nullptr);
    return value[index];
  }
  constexpr auto operator[](size_t index) const -> const T & {
    invariant(value != nullptr);
    return value[index];
  }
  template <typename C> [[nodiscard]] auto dyn_cast() -> C * {
    return llvm::dyn_cast<C>(value);
  }
  template <typename C> [[nodiscard]] auto cast() -> NotNull<C> {
    return NotNull<C>{llvm::cast<C>(value)};
  }
  template <typename C> [[nodiscard]] auto dyn_cast() const -> C * {
    return llvm::dyn_cast<C>(value);
  }
  template <typename C> [[nodiscard]] auto cast() const -> NotNull<C> {
    return NotNull<C>{llvm::cast<C>(value)};
  }
  template <typename C> [[nodiscard]] auto isa() const -> bool {
    return llvm::isa<C>(value);
  }
  constexpr auto operator+(size_t offset) -> NotNull<T> {
    invariant(value != nullptr);
    return value + offset;
  }
  constexpr auto operator-(size_t offset) -> NotNull<T> {
    invariant(value != nullptr);
    return value - offset;
  }
  constexpr auto operator+(size_t offset) const -> NotNull<T> {
    invariant(value != nullptr);
    return value + offset;
  }
  constexpr auto operator-(size_t offset) const -> NotNull<T> {
    invariant(value != nullptr);
    return value - offset;
  }
  constexpr auto operator++() -> NotNull<T> & {
    invariant(value != nullptr);
    ++value;
    return *this;
  }
  constexpr auto operator++(int) -> NotNull<T> {
    invariant(value != nullptr);
    return value++;
  }
  constexpr auto operator--() -> NotNull<T> & {
    invariant(value != nullptr);
    --value;
    return *this;
  }
  constexpr auto operator--(int) -> NotNull<T> {
    invariant(value != nullptr);
    return value--;
  }

private:
  [[no_unique_address]] T *value;
};
template <typename T> NotNull(T &) -> NotNull<T>;
template <typename T> NotNull(T *) -> NotNull<T *>;
