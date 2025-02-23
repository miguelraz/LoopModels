
#pragma once
#include "../include/Loops.hpp"
#include "../include/MemoryAccess.hpp"
#include "Math/Math.hpp"
#include <cstdint>
#include <llvm/Support/Allocator.h>

struct ArrayReference {
  const llvm::SCEVUnknown *basePointer;
  AffineLoopNest<> *loop;
  IntMatrix indMat;
  IntMatrix offMat;
  llvm::SmallVector<const llvm::SCEV *, 3> sizes;
  ArrayReference(const llvm::SCEVUnknown *p, AffineLoopNest<> &l, size_t dim)
    : basePointer(p), loop(&l), indMat(loop->getNumLoops(), dim),
      offMat(dim, 1), sizes(dim) {}
  ArrayReference(const ArrayReference &other, AffineLoopNest<> *al,
                 PtrMatrix<int64_t> iM)
    : basePointer(other.basePointer), loop(al), indMat(iM),
      offMat(other.offMat), sizes(other.sizes) {}
  auto indexMatrix() -> MutPtrMatrix<int64_t> { return indMat; }
  auto offsetMatrix() -> MutPtrMatrix<int64_t> { return offMat; }
  [[nodiscard]] auto getArrayDim() const -> size_t {
    return size_t(offMat.numRow());
  }
};
inline auto createMemAccess(llvm::BumpPtrAllocator &alloc, ArrayReference &ar,
                            llvm::Instruction *I,
                            llvm::ArrayRef<unsigned> omegas)
  -> NotNull<MemoryAccess> {

  IntMatrix indMatT(ar.indMat.transpose());
  return MemoryAccess::construct(alloc, ar.basePointer, *ar.loop, I, indMatT,
                                 {ar.sizes, {}}, ar.offsetMatrix(), omegas);
}
