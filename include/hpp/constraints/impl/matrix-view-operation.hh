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

#define HPP_EIGEN_DECLARE_TEMPLATE_ARGS_MATRIX_BLOCK_VIEW \
  typename _ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols
#define HPP_EIGEN_MATRIX_BLOCK_VIEW \
  MatrixBlockView<_ArgType, _Rows, _Cols, _allRows, _allCols>

  /** matrix op view */
  template <typename BinaryOp, typename Lhs, HPP_EIGEN_DECLARE_TEMPLATE_ARGS_MATRIX_BLOCK_VIEW>
  class CwiseBinaryOpImpl <BinaryOp, Lhs, const HPP_EIGEN_MATRIX_BLOCK_VIEW, Dense>
    : public internal::dense_xpr_base< CwiseBinaryOp<BinaryOp, Lhs, const HPP_EIGEN_MATRIX_BLOCK_VIEW > >::type
  {
      typedef const HPP_EIGEN_MATRIX_BLOCK_VIEW View;
      typedef CwiseBinaryOp<BinaryOp, Lhs, View > Derived;
    public:

      typedef typename internal::dense_xpr_base<Derived >::type Base;
      EIGEN_DENSE_PUBLIC_INTERFACE( Derived )

      template <typename OtherDerived>
      void evalTo (MatrixBase<OtherDerived>& other) const
      {
        typedef Block<Lhs> BlockLhs;
        typedef CwiseBinaryOp < BinaryOp, BlockLhs,
          typename View::template block_t< typename View::ArgType >::type
            > BlockCwiseBinaryOp;

        typedef typename Derived::Index Index;
        const Derived& d = derived();
        Index r = 0, c = 0;
        for(typename Derived::Index k = 0; k < d.rhs()._blocks(); ++k) {
          typename View::template block_t< typename View::ArgType >::type
            rhs = d.rhs()._block(k);
          BlockLhs lhs = d.lhs().block(r, c, rhs.rows(), rhs.cols());
          other.derived().block(r, c, rhs.rows(), rhs.cols())
            = BlockCwiseBinaryOp (lhs, rhs, d.functor());
          // Iteration over blocks is rowwise.
          c += rhs.cols();
          if (c >= other.cols()) {
            c = 0;
            r += rhs.rows();
          }
        }
      }
  };

// #define HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR(eval_before_assign, need_to_transpose)
#define HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR(need_to_transpose) \
    template<typename Derived, typename BinaryOp, HPP_EIGEN_LHS_TPL, HPP_EIGEN_RHS_TPL> \
    struct assign_selector<Derived, CwiseBinaryOp<BinaryOp, HPP_EIGEN_LHS_TYPE, HPP_EIGEN_RHS_TYPE >,false,need_to_transpose> { \
      typedef CwiseBinaryOp<BinaryOp, HPP_EIGEN_LHS_TYPE, HPP_EIGEN_RHS_TYPE> CwiseDerived; \
      static EIGEN_STRONG_INLINE Derived& run(Derived& dst, const CwiseDerived& other) { dst.resize(other.rows(), other.cols()); other.evalTo(dst); return dst; } \
      template<typename ActualDerived, typename ActualOtherDerived> \
        static EIGEN_STRONG_INLINE Derived& evalTo(ActualDerived& dst, const ActualOtherDerived& other) { HPP_EIGEN_EVAL_TO_BODY return dst; } \
    };

  namespace internal {
#define HPP_EIGEN_EVAL_TO_BODY other.evalTo(dst);
#define HPP_EIGEN_LHS_TPL typename OtherDerived
#define HPP_EIGEN_LHS_TYPE OtherDerived
#define HPP_EIGEN_RHS_TPL HPP_EIGEN_DECLARE_TEMPLATE_ARGS_MATRIX_BLOCK_VIEW
#define HPP_EIGEN_RHS_TYPE const HPP_EIGEN_MATRIX_BLOCK_VIEW
    HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR(false)
#undef HPP_EIGEN_LHS_TPL
#undef HPP_EIGEN_LHS_TYPE
#undef HPP_EIGEN_RHS_TPL
#undef HPP_EIGEN_RHS_TYPE

// #define HPP_EIGEN_LHS_TPL HPP_EIGEN_DECLARE_TEMPLATE_ARGS_MATRIX_BLOCK_VIEW
// #define HPP_EIGEN_LHS_TYPE const HPP_EIGEN_MATRIX_BLOCK_VIEW
// #define HPP_EIGEN_RHS_TPL typename OtherDerived
// #define HPP_EIGEN_RHS_TYPE OtherDerived
    // HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR(false)
// #undef HPP_EIGEN_LHS_TPL
// #undef HPP_EIGEN_LHS_TYPE
// #undef HPP_EIGEN_RHS_TPL
// #undef HPP_EIGEN_RHS_TYPE

#undef HPP_EIGEN_EVAL_TO_BODY

#define HPP_EIGEN_EVAL_TO_BODY Transpose<ActualDerived> dstTrans(dst); other.evalTo(dstTrans);
#define HPP_EIGEN_LHS_TPL typename OtherDerived
#define HPP_EIGEN_LHS_TYPE OtherDerived
#define HPP_EIGEN_RHS_TPL HPP_EIGEN_DECLARE_TEMPLATE_ARGS_MATRIX_BLOCK_VIEW
#define HPP_EIGEN_RHS_TYPE const HPP_EIGEN_MATRIX_BLOCK_VIEW
    HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR(true)
#undef HPP_EIGEN_LHS_TPL
#undef HPP_EIGEN_LHS_TYPE
#undef HPP_EIGEN_RHS_TPL
#undef HPP_EIGEN_RHS_TYPE

// #define HPP_EIGEN_LHS_TPL HPP_EIGEN_DECLARE_TEMPLATE_ARGS_MATRIX_BLOCK_VIEW
// #define HPP_EIGEN_LHS_TYPE const HPP_EIGEN_MATRIX_BLOCK_VIEW
// #define HPP_EIGEN_RHS_TPL typename OtherDerived
// #define HPP_EIGEN_RHS_TYPE OtherDerived
    // HPP_EIGEN_SPECIALIZE_ASSIGN_SELECTOR(true)
// #undef HPP_EIGEN_LHS_TPL
// #undef HPP_EIGEN_LHS_TYPE
// #undef HPP_EIGEN_RHS_TPL
// #undef HPP_EIGEN_RHS_TYPE

#undef HPP_EIGEN_EVAL_TO_BODY

  } // namespace internal

#undef HPP_EIGEN_DECLARE_TEMPLATE_ARGS_MATRIX_BLOCK_VIEW
#undef HPP_EIGEN_MATRIX_BLOCK_VIEW

} // namespace Eigen

#endif // HPP_CONSTRAINTS_MATRIX_VIEW_OPERATION_HH
