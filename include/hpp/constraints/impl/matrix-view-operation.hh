// Copyright (c) 2017, Joseph Mirabel
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//
// This file is part of hpp-constraints.
// hpp-constraints is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-constraints is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-constraints. If not, see <http://www.gnu.org/licenses/>.

#ifndef HPP_CONSTRAINTS_IMPL_MATRIX_VIEW_OPERATION_HH
#define HPP_CONSTRAINTS_IMPL_MATRIX_VIEW_OPERATION_HH

namespace Eigen {
  /** Support for CwiseBinaryOp
   *  3 possible cases:
   *  - view op view
   *  - matrix op view
   *  - view op matrix */

#define HPP_EIGEN_DECLARE_TEMPLATE_ARGS_MATRIX_BLOCK_VIEW                      \
  typename _ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols
#define HPP_EIGEN_MATRIX_BLOCK_VIEW                                            \
  MatrixBlockView<_ArgType, _Rows, _Cols, _allRows, _allCols>

#define HPP_EIGEN_SPECIALIZE_CwiseBinaryOpImpl(                                \
    LHS_TPL, LHS_TYPE, RHS_TPL, RHS_TYPE)                                      \
  template <typename BinaryOp, LHS_TPL, RHS_TPL>                               \
  class CwiseBinaryOpImpl <BinaryOp, LHS_TYPE, RHS_TYPE, Dense>                \
    : public internal::dense_xpr_base<                                         \
      CwiseBinaryOp<BinaryOp, LHS_TYPE, RHS_TYPE > >::type                     \
  {                                                                            \
      typedef LHS_TYPE Lhs_t;                                                  \
      typedef RHS_TYPE Rhs_t;                                                  \
      typedef CwiseBinaryOp<BinaryOp, LHS_TYPE, RHS_TYPE > Derived;            \
    public:                                                                    \
                                                                               \
      typedef typename internal::dense_xpr_base<Derived >::type Base;          \
      EIGEN_DENSE_PUBLIC_INTERFACE( Derived )                                  \
                                                                               \
      template <typename OtherDerived>                                         \
      void evalTo (MatrixBase<OtherDerived>& other) const;                     \
  };

#define HPP_EIGEN_DEFINE_CwiseBinaryOpImpl_evalTo(                             \
    LHS_TPL, LHS_TYPE, RHS_TPL, RHS_TYPE)                                      \
template <typename BinaryOp, LHS_TPL, RHS_TPL>                                 \
template <typename OtherDerived>                                               \
void CwiseBinaryOpImpl<BinaryOp, LHS_TYPE, RHS_TYPE, Dense>::evalTo            \
(MatrixBase<OtherDerived>& other) const

# if HPP_EIGEN_USE_EVALUATOR

#define HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR(                                  \
    LHS_TPL, LHS_TYPE, RHS_TPL, RHS_TYPE)                                      \
    template<typename Derived, typename BinaryOp, LHS_TPL, RHS_TPL,            \
             typename Functor, typename Scalar>                                \
    struct Assignment<Derived,                                                 \
                      CwiseBinaryOp<BinaryOp, LHS_TYPE, RHS_TYPE >,            \
                      Functor, Dense2Dense, Scalar> {                          \
      typedef CwiseBinaryOp<BinaryOp, LHS_TYPE, RHS_TYPE> CwiseDerived;        \
      static EIGEN_STRONG_INLINE void run                                      \
        (Derived& dst, const CwiseDerived& o, const Functor&)                  \
      { o.evalTo(dst); }                                                       \
    };

# else // HPP_EIGEN_USE_EVALUATOR

#define HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR_IMPL(                             \
    LHS_TPL, LHS_TYPE, RHS_TPL, RHS_TYPE,                                      \
    need_to_transpose, EVAL_TO_BODY)                                           \
    template<typename Derived, typename BinaryOp, LHS_TPL, RHS_TPL>            \
    struct assign_selector<Derived,                                            \
                           CwiseBinaryOp<BinaryOp, LHS_TYPE, RHS_TYPE >,       \
                           false,need_to_transpose> {                          \
      typedef CwiseBinaryOp<BinaryOp, LHS_TYPE, RHS_TYPE> CwiseDerived;        \
      static EIGEN_STRONG_INLINE Derived& run                                  \
        (Derived& dst, const CwiseDerived& o)                                  \
      { dst.resize(o.rows(), o.cols()); o.evalTo(dst); return dst; }           \
      template<typename ActualDerived, typename ActualOtherDerived>            \
      static EIGEN_STRONG_INLINE Derived& evalTo                               \
        (ActualDerived& dst, const ActualOtherDerived& other)                  \
      { EVAL_TO_BODY return dst; }                                             \
    };

#define HPP_EIGEN_EVAL_TO_BODY_NORMAL other.evalTo(dst);
#define HPP_EIGEN_EVAL_TO_BODY_TRANSPOSE Transpose<ActualDerived> dstTrans(dst); other.evalTo(dstTrans);

#define HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR(                                  \
    LHS_TPL, LHS_TYPE, RHS_TPL, RHS_TYPE)                                      \
  HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR_IMPL(                                   \
      HPP_EIGEN_LHS_TPL, HPP_EIGEN_LHS_TYPE,                                   \
      HPP_EIGEN_RHS_TPL, HPP_EIGEN_RHS_TYPE,                                   \
      false, HPP_EIGEN_EVAL_TO_BODY_NORMAL)                                    \
    HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR_IMPL(                                 \
      HPP_EIGEN_LHS_TPL, HPP_EIGEN_LHS_TYPE,                                   \
      HPP_EIGEN_RHS_TPL, HPP_EIGEN_RHS_TYPE,                                   \
      true, HPP_EIGEN_EVAL_TO_BODY_TRANSPOSE)

# endif // HPP_EIGEN_USE_EVALUATOR

  // --- matrix op view -- //
#define HPP_EIGEN_LHS_TPL typename Lhs
#define HPP_EIGEN_LHS_TYPE Lhs
#define HPP_EIGEN_RHS_TPL HPP_EIGEN_DECLARE_TEMPLATE_ARGS_MATRIX_BLOCK_VIEW
#define HPP_EIGEN_RHS_TYPE const HPP_EIGEN_MATRIX_BLOCK_VIEW

  HPP_EIGEN_SPECIALIZE_CwiseBinaryOpImpl(
      HPP_EIGEN_LHS_TPL, HPP_EIGEN_LHS_TYPE,
      HPP_EIGEN_RHS_TPL, HPP_EIGEN_RHS_TYPE)
  HPP_EIGEN_DEFINE_CwiseBinaryOpImpl_evalTo(
      HPP_EIGEN_LHS_TPL, HPP_EIGEN_LHS_TYPE,
      HPP_EIGEN_RHS_TPL, HPP_EIGEN_RHS_TYPE)
      {
        typedef const Block<Lhs_t> BlockLhs;
        typedef const typename Rhs_t::
          template block_t< typename Rhs_t::ArgType >::type BlockRhs;
        typedef CwiseBinaryOp < BinaryOp, BlockLhs, BlockRhs > BlockCwiseBOp;

        const Derived& d = derived();
        for(typename Rhs_t::block_iterator block (d.rhs()); block.valid(); ++block) {
          BlockRhs rhs = d.rhs()._block(block);
          BlockLhs lhs = d.lhs().block(block.ro(), block.co(), block.rs(), block.cs());
          other.derived().block(block.ro(), block.co(), block.rs(), block.cs())
            = BlockCwiseBOp (lhs, rhs, d.functor());
        }
      }

  namespace internal {
    HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR(
      HPP_EIGEN_LHS_TPL, HPP_EIGEN_LHS_TYPE,
      HPP_EIGEN_RHS_TPL, HPP_EIGEN_RHS_TYPE)
  }
#undef HPP_EIGEN_LHS_TPL
#undef HPP_EIGEN_LHS_TYPE
#undef HPP_EIGEN_RHS_TPL
#undef HPP_EIGEN_RHS_TYPE

  // --- view op matrix -- //
#define HPP_EIGEN_LHS_TPL HPP_EIGEN_DECLARE_TEMPLATE_ARGS_MATRIX_BLOCK_VIEW
#define HPP_EIGEN_LHS_TYPE const HPP_EIGEN_MATRIX_BLOCK_VIEW
#define HPP_EIGEN_RHS_TPL typename Rhs
#define HPP_EIGEN_RHS_TYPE Rhs
  HPP_EIGEN_SPECIALIZE_CwiseBinaryOpImpl(
      HPP_EIGEN_LHS_TPL, HPP_EIGEN_LHS_TYPE,
      HPP_EIGEN_RHS_TPL, HPP_EIGEN_RHS_TYPE)
  HPP_EIGEN_DEFINE_CwiseBinaryOpImpl_evalTo(
      HPP_EIGEN_LHS_TPL, HPP_EIGEN_LHS_TYPE,
      HPP_EIGEN_RHS_TPL, HPP_EIGEN_RHS_TYPE)
      {
        typedef const typename Lhs_t::
          template block_t< typename Lhs_t::ArgType >::type BlockLhs;
        typedef const Block<Rhs_t> BlockRhs;
        typedef CwiseBinaryOp < BinaryOp, BlockLhs, BlockRhs > BlockCwiseBOp;

        const Derived& d = derived();
        for(typename Lhs_t::block_iterator block (d.lhs()); block.valid(); ++block) {
          BlockLhs lhs = d.lhs()._block(block);
          BlockRhs rhs = d.rhs().block(block.ro(), block.co(), block.rs(), block.cs());
          other.derived().block(block.ro(), block.co(), block.rs(), block.cs())
            = BlockCwiseBOp (lhs, rhs, d.functor());
        }
      }

  namespace internal {
    HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR(
      HPP_EIGEN_LHS_TPL, HPP_EIGEN_LHS_TYPE,
      HPP_EIGEN_RHS_TPL, HPP_EIGEN_RHS_TYPE)
  }
#undef HPP_EIGEN_LHS_TPL
#undef HPP_EIGEN_LHS_TYPE
#undef HPP_EIGEN_RHS_TPL
#undef HPP_EIGEN_RHS_TYPE

  // --- view op view -- //
#define HPP_EIGEN_LHS_TPL HPP_EIGEN_DECLARE_TEMPLATE_ARGS_MATRIX_BLOCK_VIEW
#define HPP_EIGEN_LHS_TYPE const HPP_EIGEN_MATRIX_BLOCK_VIEW
#define HPP_EIGEN_RHS_TPL typename _ArgType2, int _Rows2, int _Cols2, bool _allRows2, bool _allCols2
#define HPP_EIGEN_RHS_TYPE const MatrixBlockView<_ArgType2, _Rows2, _Cols2, _allRows2, _allCols2>

  HPP_EIGEN_SPECIALIZE_CwiseBinaryOpImpl(
      HPP_EIGEN_LHS_TPL, HPP_EIGEN_LHS_TYPE,
      HPP_EIGEN_RHS_TPL, HPP_EIGEN_RHS_TYPE)
  HPP_EIGEN_DEFINE_CwiseBinaryOpImpl_evalTo(
      HPP_EIGEN_LHS_TPL, HPP_EIGEN_LHS_TYPE,
      HPP_EIGEN_RHS_TPL, HPP_EIGEN_RHS_TYPE)
      {
        typedef const typename Lhs_t::
          template block_t< typename Lhs_t::ArgType >::type BlockLhs;
        typedef const typename Rhs_t::
          template block_t< typename Rhs_t::ArgType >::type BlockRhs;
        typedef CwiseBinaryOp < BinaryOp, BlockLhs, BlockRhs > BlockCwiseBOp;

        const Derived& d = derived();
        assert (d.lhs()._blocks() == d.rhs()._blocks());
        typename Lhs_t::block_iterator lblock (d.lhs());
        typename Rhs_t::block_iterator rblock (d.rhs());
        while (lblock.valid()) {
          BlockLhs lhs = d.lhs()._block(lblock);
          BlockRhs rhs = d.rhs()._block(lblock);
          assert (lblock.rs() == rblock.rs() && lblock.cs() == rblock.cs());
          assert (lblock.ro() == rblock.ro() && lblock.co() == rblock.co());
          other.derived().block(rblock.ro(), rblock.co(), rblock.rs(), rblock.cs())
            = BlockCwiseBOp (lhs, rhs, d.functor());
          ++lblock; ++rblock;
        }
        assert (!lblock.valid() && !rblock.valid());
      }

  namespace internal {
    HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR(
      HPP_EIGEN_LHS_TPL, HPP_EIGEN_LHS_TYPE,
      HPP_EIGEN_RHS_TPL, HPP_EIGEN_RHS_TYPE)
  }
#undef HPP_EIGEN_LHS_TPL
#undef HPP_EIGEN_LHS_TYPE
#undef HPP_EIGEN_RHS_TPL
#undef HPP_EIGEN_RHS_TYPE

# if !HPP_EIGEN_USE_EVALUATOR
#  undef HPP_EIGEN_EVAL_TO_BODY_NORMAL
#  undef HPP_EIGEN_EVAL_TO_BODY_TRANSPOSE
#  undef HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR_IMPL
# endif // !HPP_EIGEN_USE_EVALUATOR

#undef HPP_EIGEN_SPECIALIZE_CwiseBinaryOpImpl
#undef HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR
#undef HPP_EIGEN_DEFINE_CwiseBinaryOpImpl_evalTo
#undef HPP_EIGEN_DECLARE_TEMPLATE_ARGS_MATRIX_BLOCK_VIEW
#undef HPP_EIGEN_MATRIX_BLOCK_VIEW

} // namespace Eigen

#endif // HPP_CONSTRAINTS_MATRIX_VIEW_OPERATION_HH
