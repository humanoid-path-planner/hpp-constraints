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

#ifndef HPP_CONSTRAINTS_MATRIX_VIEW_HH
#define HPP_CONSTRAINTS_MATRIX_VIEW_HH

#include <Eigen/Core>
#include <vector>
#include <iostream>

namespace Eigen {
  // template <typename Index> struct interval {
    // Index m_start, m_length;
    // interval (const Index& i, const Index length) : m_start (start), m_length (length) {}
  // };
  template <typename ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols> class MatrixView;

  namespace internal {
    template <typename ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols>
      struct traits< MatrixView <ArgType, _Rows, _Cols, _allRows, _allCols> >
    {
      typedef typename ArgType::Index Index;
      typedef Eigen::Dense StorageKind;
      typedef Eigen::MatrixXpr XprKind;
      // typedef typename ArgType::StorageIndex StorageIndex;
      typedef typename ArgType::Scalar Scalar;
      enum { 
        CoeffReadCost = ArgType::CoeffReadCost,
        Flags = ~AlignedBit & ~DirectAccessBit & ~ActualPacketAccessBit & ~LinearAccessBit & ArgType::Flags,
        RowsAtCompileTime = (_allRows ? ArgType::RowsAtCompileTime : _Rows),
        ColsAtCompileTime = (_allCols ? ArgType::ColsAtCompileTime : _Cols),
        MaxRowsAtCompileTime = ArgType::MaxRowsAtCompileTime,
        MaxColsAtCompileTime = ArgType::MaxColsAtCompileTime
      };
    };

      template <bool row> struct return_first {
        template <typename First, typename Second>
        static inline First& run (First& f, Second&) { return f; }
      };
      template <> struct return_first <false> {
        template <typename First, typename Second>
        static inline Second& run (First&, Second& s) { return s; }
      };

      struct empty_struct {
        typedef MatrixXd::Index Index;
        empty_struct () {}
        template <typename In_t> empty_struct (In_t) {}
        static inline Index size() { return 0; }
        inline const Index& operator[](const Index& i) const { return i; }
      };
  } // namespace internal

  template <int _Rows, int _Cols, bool _allRows, bool _allCols>
  class MatrixIndexes
  {
    public:
      typedef MatrixXd::Index Index;
      typedef std::vector<Index> Indexes_t;
      typedef typename internal::conditional<_allRows, internal::empty_struct, Indexes_t>::type RowIndexes_t;
      typedef typename internal::conditional<_allCols, internal::empty_struct, Indexes_t>::type ColIndexes_t;

      template <typename Derived> struct View {
        typedef MatrixView<Derived, _Rows, _Cols, _allRows, _allCols> type;
      };

      MatrixIndexes () : m_rows(), m_cols() {}

      MatrixIndexes (const Index& rows, const Index& cols)
        : m_rows(rows), m_cols(cols)
      {
        // for (Index i = 0; i < m_rows.size(); ++i) m_rows[i] = i;
        // for (Index i = 0; i < m_cols.size(); ++i) m_cols[i] = i;
      }

      MatrixIndexes (const Index& size)
        : m_rows(size), m_cols(size)
      {
        // for (Index i = 0; i < indexes().size(); ++i) indexes()[i] = i;
      }

      MatrixIndexes (const Indexes_t& rows, const Indexes_t& cols)
        : m_rows(rows), m_cols(cols) {}

      /// Valid only when _allRows or _allCols is true
      MatrixIndexes (const Indexes_t& indexes)
        : m_rows(indexes), m_cols(indexes)
      {}
      
      template <typename Derived>
      EIGEN_STRONG_INLINE typename View<Derived>::type view(MatrixBase<Derived>& other) const {
        return View<Derived>::type (other.derived(), m_rows, m_cols);
      }

      inline const Indexes_t& indexes() const
      {
        EIGEN_STATIC_ASSERT(_allRows && _allCols, internal::YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX)
        return internal::return_first<_allRows>::run(m_rows, m_cols);
      }

      RowIndexes_t m_rows;
      ColIndexes_t m_cols;
  };

  template <typename ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols>
  class MatrixView : public MatrixBase< MatrixView<ArgType, _Rows, _Cols, _allRows, _allCols> >
  {
    public:
      typedef MatrixBase< MatrixView<ArgType, _Rows, _Cols, _allRows, _allCols> > Base;
      EIGEN_GENERIC_PUBLIC_INTERFACE(MatrixView)

      typedef Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> PlainObject;
      // typedef typename internal::ref_selector<MatrixView>::type Nested; 
      typedef typename internal::ref_selector<ArgType>::type ArgTypeNested;
      // typedef typename Base::CoeffReturnType CoeffReturnType;
      // typedef typename Base::Scalar Scalar;

      typedef MatrixIndexes<_Rows, _Cols, _allRows, _allCols> MatrixIndexes_t;
      typedef typename MatrixIndexes_t::Indexes_t Indexes_t;
      typedef typename MatrixIndexes_t::RowIndexes_t RowIndexes_t;
      typedef typename MatrixIndexes_t::ColIndexes_t ColIndexes_t;

      using Base::operator=;

      MatrixView (ArgType& arg, const Indexes_t& rows, const Indexes_t& cols)
        : m_arg (arg), m_rows(rows), m_cols(cols) {}

      /// Valid only when _allRows or _allCols is true
      MatrixView (ArgType& arg, const Indexes_t& indexes)
        : m_arg (arg), m_rows(indexes), m_cols(indexes)
      {}
      
      EIGEN_STRONG_INLINE Index rows() const { if (_allRows) return m_arg.rows(); else return m_rows.size(); }
      EIGEN_STRONG_INLINE Index cols() const { if (_allCols) return m_arg.cols(); else return m_cols.size(); }

      EIGEN_STRONG_INLINE const Index& argIndex(const Index& index) const {
        // EIGEN_STATIC_ASSERT_VECTOR_ONLY(PlainObject)
        if (rows() == 1) return argCol(index);
        else             return argRow(index);
      }
      EIGEN_STRONG_INLINE const Index& argRow(const Index& row) const { if (_allRows) return row; else return m_rows[row]; }
      EIGEN_STRONG_INLINE const Index& argCol(const Index& col) const { if (_allCols) return col; else return m_cols[col]; }

      EIGEN_STRONG_INLINE CoeffReturnType coeff(Index index) const { return m_arg.coeff(argIndex(index)); };
      EIGEN_STRONG_INLINE CoeffReturnType coeff(Index row, Index col) const { return m_arg.coeff(argRow(row), argCol(col)); };
      EIGEN_STRONG_INLINE Scalar& coeffRef(Index index) { return m_arg.coeffRef(argIndex(index)); };
      EIGEN_STRONG_INLINE Scalar& coeffRef(Index row, const Index& col) { return m_arg.coeffRef(argRow(row), argCol(col)); };

      ArgType& m_arg;
      // Indexes_t m_rows, m_cols;
      const RowIndexes_t& m_rows;
      const ColIndexes_t& m_cols;
  };

  template <typename ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols> class MatrixBlockView;

  namespace internal {
    template <typename ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols>
      struct traits< MatrixBlockView <ArgType, _Rows, _Cols, _allRows, _allCols> >
    {
      typedef typename ArgType::Index Index;
      typedef Eigen::Dense StorageKind;
      typedef Eigen::MatrixXpr XprKind;
      // typedef typename ArgType::StorageIndex StorageIndex;
      typedef typename ArgType::Scalar Scalar;
      enum { 
        CoeffReadCost = ArgType::CoeffReadCost,
        Flags = ~AlignedBit & ~DirectAccessBit & ~ActualPacketAccessBit & ~LinearAccessBit & ArgType::Flags,
        RowsAtCompileTime = (_allRows ? ArgType::RowsAtCompileTime : _Rows),
        ColsAtCompileTime = (_allCols ? ArgType::ColsAtCompileTime : _Cols),
        MaxRowsAtCompileTime = ArgType::MaxRowsAtCompileTime,
        MaxColsAtCompileTime = ArgType::MaxColsAtCompileTime
      };
    };

    template<typename Derived, typename ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols>
    struct assign_selector<Derived, MatrixBlockView <ArgType, _Rows, _Cols, _allRows, _allCols>,false,false> {
      typedef MatrixBlockView <ArgType, _Rows, _Cols, _allRows, _allCols> OtherDerived;
      static EIGEN_STRONG_INLINE Derived& run(Derived& dst, const OtherDerived& other) { other.writeTo(dst); return dst; }
      template<typename ActualDerived, typename ActualOtherDerived>
        static EIGEN_STRONG_INLINE Derived& evalTo(ActualDerived& dst, const ActualOtherDerived& other) { other.evalTo(dst); return dst; }
    };
    template<typename Derived, typename ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols>
    struct assign_selector<Derived, MatrixBlockView <ArgType, _Rows, _Cols, _allRows, _allCols>,false,true> {
      typedef MatrixBlockView <ArgType, _Rows, _Cols, _allRows, _allCols> OtherDerived;
      static EIGEN_STRONG_INLINE Derived& run(Derived& dst, const OtherDerived& other) { other.writeTo(dst.transpose()); return dst; }
      template<typename ActualDerived, typename ActualOtherDerived>
        static EIGEN_STRONG_INLINE Derived& evalTo(ActualDerived& dst, const ActualOtherDerived& other) { Transpose<ActualDerived> dstTrans(dst); other.evalTo(dstTrans); return dst; }
    };

    template <typename Dst, typename Src, bool _allRows, bool _allCols> struct evalCols {
      static inline void run (Dst& dst, const Src& src, const typename Dst::Index& row)
      {
        std::size_t col = 0;
        for (std::size_t j = 0; j < src.m_cols.size(); ++j) {
          dst.derived().middleCols(col, src.m_cols[j].second) =
            src.m_arg.block(row,        src.m_cols[j].first,
                            dst.rows(), src.m_cols[j].second);
          col += src.m_cols[j].second;
        }
      }
      static inline void write (const Dst& dst, Src& src, const typename Dst::Index& row)
      {
        std::size_t col = 0;
        for (std::size_t j = 0; j < src.m_cols.size(); ++j) {
          src.m_arg.block(row,        src.m_cols[j].first,
                          dst.rows(), src.m_cols[j].second)
            = dst.middleCols(col, src.m_cols[j].second);
          col += src.m_cols[j].second;
        }
      }
    };
    template <typename Dst, typename Src> struct evalCols<Dst, Src, true , false> {
      static inline void run (Dst& dst, const Src& src)
      {
        std::size_t col = 0;
        for (std::size_t j = 0; j < src.m_cols.size(); ++j) {
          dst.derived().middleCols(col, src.m_cols[j].second) =
            src.m_arg.middleCols(src.m_cols[j].first, src.m_cols[j].second);
          col += src.m_cols[j].second;
        }
      }
      static inline void write (Dst& dst, const Src& src)
      {
        std::size_t col = 0;
        for (std::size_t j = 0; j < src.m_cols.size(); ++j) {
          src.m_arg.middleCols(src.m_cols[j].first, src.m_cols[j].second)
            = dst.derived().middleCols(col, src.m_cols[j].second);
          col += src.m_cols[j].second;
        }
      }
    };
    template <typename Dst, typename Src> struct evalCols<Dst, Src, false, true > {
      static inline void run (Dst& dst, const Src& src, const typename Dst::Index& row)
      {
        dst.derived() = src.m_arg.middleRows(row, dst.rows());
      }
      static inline void write (Dst& dst, Src& src, const typename Dst::Index& row)
      {
        src.m_arg.middleRows(row, dst.rows()) = dst;
      }
    };
    template <typename Dst, typename Src> struct evalCols<Dst, Src, true , true > {};
    template <typename Dst, typename Src, bool _allRows, bool _allCols> struct evalRows {
      static inline void run (Dst& dst, const Src& src)
      {
        evalCols<Dst, Src, _allRows, _allCols>::run(dst.derived(), src);
      }
      static inline void write (Dst& dst, Src& src)
      {
        evalCols<Dst, Src, _allRows, _allCols>::write(dst, src);
      }
    };
    template <typename Dst, typename Src, bool _allCols> struct evalRows <Dst, Src, false, _allCols> {
      static inline void run (Dst& dst, const Src& src)
      {
        std::size_t row = 0;
        for (std::size_t i = 0; i < src.m_rows.size(); ++i) {
          typename Dst::RowsBlockXpr rows = dst.middleRows(row, src.m_rows[i].second);
          evalCols<typename Dst::RowsBlockXpr, Src, false, _allCols>::run(rows,
              src, src.m_rows[i].first);
          row += src.m_rows[i].second;
        }
      }
      static inline void write (Dst& dst, Src& src)
      {
        std::size_t row = 0;
        for (std::size_t i = 0; i < src.m_rows.size(); ++i) {
          typename Dst::ConstRowsBlockXpr rows = dst.middleRows(row, src.m_rows[i].second);
          evalCols<typename Dst::ConstRowsBlockXpr, Src, false, _allCols>::write(rows,
              src, src.m_rows[i].first);
          row += src.m_rows[i].second;
        }
      }
    };
    template <bool print> struct print_indexes { template <typename BlockIndexType> static void run (std::ostream&, const BlockIndexType&) {} };
    template <> struct print_indexes <true> {
      template <typename BlockIndexType>
      static void run (std::ostream& os, const BlockIndexType& bi) {
        for (std::size_t i = 0; i < bi.size(); ++i)
          os << "[ " << bi[i].first << ", " << bi[i].second << "], ";
      }
    };
  } // namespace internal

  template <typename IndexType>
  struct BlockIndex {
    typedef IndexType Index;
    typedef std::pair<Index, Index> type;
    typedef std::vector<type> vector_t;

    static IndexType cardinal (const vector_t& a);

    template <typename Derived>
    static vector_t fromLogicalExpression (const Eigen::ArrayBase<Derived>& array);

    static void sort   (vector_t& a);
    /// Assumes a is sorted
    static void shrink (vector_t& a);

    static bool overlap (const type& a, const type& b);
    /// The sum is the union
    static vector_t sum (const type& a, const type& b);
    /// Assumes b is sorted
    static vector_t difference (const type    & a, const type    & b);
    /// Assumes b is sorted
    static vector_t difference (const vector_t& a, const type    & b);
    /// Assumes b is sorted
    static vector_t difference (const type    & a, const vector_t& b);
  };

  template <bool _allRows, bool _allCols>
  class MatrixBlockIndexes
  {
    public:
      typedef MatrixXd::Index Index;
      typedef BlockIndex<Index>           BlockIndex_t;
      typedef BlockIndex_t::type     BlockIndexType;
      typedef BlockIndex_t::vector_t BlockIndexesType;
      typedef typename internal::conditional<_allRows, internal::empty_struct, BlockIndexesType>::type RowIndexes_t;
      typedef typename internal::conditional<_allCols, internal::empty_struct, BlockIndexesType>::type ColIndexes_t;

      template <typename Derived, int _Rows, int _Cols> struct View {
        typedef MatrixBlockView<Derived, _Rows, _Cols, _allRows, _allCols> type;
        typedef MatrixBlockView<Derived, _Rows, _Cols, _allCols, _allRows> transpose_type;
      };

      MatrixBlockIndexes () : m_nbRows(0), m_nbCols(0), m_rows(), m_cols() {}

      /// \warning rows and cols must be sorted
      MatrixBlockIndexes (const BlockIndexesType& rows, const BlockIndexesType& cols)
        : m_nbRows(BlockIndex<Index>::cardinal(rows)), m_nbCols(BlockIndex<Index>::cardinal(rows)), m_rows(rows), m_cols(cols)
      {}

      /// \warning idx must be sorted and shrinked
      MatrixBlockIndexes (const BlockIndexesType& idx)
        : m_nbRows(_allRows ? 0 : BlockIndex<Index>::cardinal(idx))
        , m_nbCols(_allCols ? 0 : BlockIndex<Index>::cardinal(idx))
        , m_rows(idx), m_cols(idx)
      {}

      /// \warning idx must be sorted and shrinked
      MatrixBlockIndexes (const BlockIndexType& idx)
        : m_nbRows(_allRows ? 0 : idx.second)
        , m_nbCols(_allCols ? 0 : idx.second)
        , m_rows(BlockIndexesType(1,idx)), m_cols(BlockIndexesType(1,idx))
      {}

      /*
      MatrixBlockIndexes (const Index& rows, const Index& cols)
        : m_rows(rows), m_cols(cols)
      {
        // for (Index i = 0; i < m_rows.size(); ++i) m_rows[i] = i;
        // for (Index i = 0; i < m_cols.size(); ++i) m_cols[i] = i;
      }

      MatrixBlockIndexes (const Index& size)
        : m_rows(size), m_cols(size)
      {
        // for (Index i = 0; i < indexes().size(); ++i) indexes()[i] = i;
      }

      /// Valid only when _allRows or _allCols is true
      MatrixBlockIndexes (const Indexes_t& indexes)
        : m_rows(indexes), m_cols(indexes)
      {}
      */

      inline void addRow (const Index& row, const Index size)
      {
        m_rows.push_back(BlockIndexType(row, size));
        m_nbRows += size;
      }

      inline void addCol (const Index& col, const Index size)
      {
        m_cols.push_back(BlockIndexType(col, size));
        m_nbCols += size;
      }

      template<bool Sort, bool Shrink, bool Cardinal>
      inline void updateRows() {
        update<Sort, Shrink, Cardinal> (m_rows, m_nbRows);
      }

      template<bool Sort, bool Shrink, bool Cardinal>
      inline void updateCols() {
        update<Sort, Shrink, Cardinal> (m_cols, m_nbCols);
      }

      template <typename Derived>
      EIGEN_STRONG_INLINE typename View<Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::type lview(const MatrixBase<Derived>& other) const {
        Derived& o = const_cast<MatrixBase<Derived>&>(other).derived();
        if (_allCols || _allRows)
          return typename View<Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::type (o, nbIndexes(), indexes());
        else
          return typename View<Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::type (o, m_nbRows, m_rows, m_nbCols, m_cols);
      }

      template <typename Derived>
      EIGEN_STRONG_INLINE typename View<Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::transpose_type lviewTranspose(const MatrixBase<Derived>& other) const {
        Derived& o = const_cast<MatrixBase<Derived>&>(other).derived();
        if (_allCols || _allRows)
          return typename View<Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::transpose_type (o, nbIndexes(), indexes());
        else
          return typename View<Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::transpose_type (o, m_nbCols, m_cols, m_nbRows, m_rows);
      }

      template <typename Derived>
      EIGEN_STRONG_INLINE typename View<const Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::type rview(const MatrixBase<Derived>& other) const {
        if (_allCols || _allRows)
          return typename View<const Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::type (other.derived(), nbIndexes(), indexes());
        else
          return typename View<const Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::type (other.derived(), m_nbRows, m_rows, m_nbCols, m_cols);
      }

      template <typename Derived>
      EIGEN_STRONG_INLINE typename View<const Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::transpose_type rviewTranspose(const MatrixBase<Derived>& other) const {
        if (_allCols || _allRows)
          return typename View<const Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::transpose_type (other.derived(), nbIndexes(), indexes());
        else
          return typename View<const Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::transpose_type (other.derived(), m_nbCols, m_cols, m_nbRows, m_rows);
      }

      inline const BlockIndexesType& indexes() const
      {
        // EIGEN_STATIC_ASSERT(_allRows && _allCols, internal::YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX)
        return internal::return_first<_allRows>::run(m_cols, m_rows);
      }

      inline const Index& nbIndexes() const
      {
        // EIGEN_STATIC_ASSERT(_allRows && _allCols, internal::YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX)
        return internal::return_first<_allRows>::run(m_nbCols, m_nbRows);
      }

      template<bool Sort, bool Shrink, bool Cardinal>
      inline void updateIndexes() {
        update<Sort, Shrink, Cardinal> (
            internal::return_first<_allRows>::run(m_cols  , m_rows  ), 
            internal::return_first<_allRows>::run(m_nbCols, m_nbRows));
      }

      Index m_nbRows, m_nbCols;
      RowIndexes_t m_rows;
      ColIndexes_t m_cols;

    private:
      template<bool Sort, bool Shrink, bool Cardinal>
      static inline void update(BlockIndexesType& b, Index& idx) {
        if (Sort)     BlockIndex<Index>::sort(b);
        if (Shrink)   BlockIndex<Index>::shrink(b);
        if (Cardinal) idx = BlockIndex<Index>::cardinal(b);
      }
  };

  template <bool _allRows, bool _allCols>
  std::ostream& operator<< (std::ostream& os, MatrixBlockIndexes<_allRows, _allCols> mbi)
  {
    if (!_allRows) {
      os << "Rows: ";
      internal::print_indexes<!_allRows>::run (os, mbi.m_rows);
      os << '\n';
    }
    if (!_allCols) {
      os << "Cols: ";
      internal::print_indexes<!_allCols>::run (os, mbi.m_cols);
      os << '\n';
    }
    return os;
  }

  typedef Eigen::MatrixBlockIndexes<false, true> RowBlockIndexes;
  typedef Eigen::MatrixBlockIndexes<true, false> ColBlockIndexes;

  template <typename _ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols>
  class MatrixBlockView : public MatrixBase< MatrixBlockView<_ArgType, _Rows, _Cols, _allRows, _allCols> >
  {
    public:
      typedef MatrixBase< MatrixBlockView<_ArgType, _Rows, _Cols, _allRows, _allCols> > Base;
      EIGEN_GENERIC_PUBLIC_INTERFACE(MatrixBlockView)

      typedef Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> PlainObject;
      // typedef typename internal::ref_selector<MatrixBlockView>::type Nested; 
      typedef _ArgType ArgType;
      typedef typename internal::ref_selector<ArgType>::type ArgTypeNested;
      // typedef typename Base::CoeffReturnType CoeffReturnType;
      // typedef typename Base::Scalar Scalar;

      typedef MatrixBlockIndexes<_allRows, _allCols> MatrixIndexes_t;
      typedef typename MatrixIndexes_t::BlockIndexesType Indexes_t;
      typedef typename internal::conditional<_allRows, const internal::empty_struct, const Indexes_t& >::type RowIndexes_t;
      typedef typename internal::conditional<_allCols, const internal::empty_struct, const Indexes_t& >::type ColIndexes_t;

      // using Base::operator=;

      MatrixBlockView (ArgType& arg, const Index& nbRows, const RowIndexes_t rows, const Index& nbCols, const ColIndexes_t cols)
        : m_arg (arg), m_nbRows(nbRows), m_rows(rows), m_nbCols(nbCols), m_cols(cols) {}

      /// Valid only when _allRows or _allCols is true
      MatrixBlockView (ArgType& arg, const Index& nbIndexes, const Indexes_t& indexes)
        : m_arg (arg), m_nbRows(_allRows ? arg.rows() : nbIndexes), m_rows(indexes), m_nbCols(_allCols ? arg.cols() : nbIndexes), m_cols(indexes)
      {}
      
      EIGEN_STRONG_INLINE Index rows() const { return m_nbRows; }
      EIGEN_STRONG_INLINE Index cols() const { return m_nbCols; }

      EIGEN_STRONG_INLINE CoeffReturnType coeff(Index index) const {
        assert(false && "It is not possible to access the coefficients of MatrixBlockView this way.");
      }
      EIGEN_STRONG_INLINE CoeffReturnType coeff(Index row, Index col) const {
        assert(false && "It is not possible to access the coefficients of MatrixBlockView this way.");
      }
      EIGEN_STRONG_INLINE Scalar& coeffRef(Index index) {
        assert(false && "It is not possible to access the coefficients of MatrixBlockView this way.");
      }
      EIGEN_STRONG_INLINE Scalar& coeffRef(Index row, const Index& col) {
        assert(false && "It is not possible to access the coefficients of MatrixBlockView this way.");
      }
      /*
      EIGEN_STRONG_INLINE const Index& argIndex(const Index& index) const {
        // EIGEN_STATIC_ASSERT_VECTOR_ONLY(PlainObject)
        if (rows() == 1) return argCol(index);
        else             return argRow(index);
      }
      EIGEN_STRONG_INLINE const Index& argRow(const Index& row) const { if (_allRows) return row; else return m_rows[row]; }
      EIGEN_STRONG_INLINE const Index& argCol(const Index& col) const { if (_allCols) return col; else return m_cols[col]; }
      */

      template <typename Dest>
      EIGEN_STRONG_INLINE void evalTo (Dest& dst) const {
        internal::evalRows<Dest, MatrixBlockView, _allRows, _allCols>::run(dst, *this);
      }

      template <typename Dest>
      EIGEN_STRONG_INLINE void writeTo (Dest& dst) const {
        dst.resize(rows(), cols());
        // dst._resize_to_match(*this);
        evalTo(dst.derived());
      }

      EIGEN_STRONG_INLINE PlainObject eval () const {
        PlainObject dst;
        writeTo(dst);
        return dst;
      }

      template <typename OtherDerived>
      EIGEN_STRONG_INLINE MatrixBlockView& operator= (const EigenBase<OtherDerived>& other) {
        EIGEN_STATIC_ASSERT_LVALUE(ArgType);
        internal::evalRows<const OtherDerived, MatrixBlockView, _allRows, _allCols>::write(other.derived(), *this);
        return *this;
      }

      ArgType& m_arg;
      Index m_nbRows;
      RowIndexes_t m_rows;
      Index m_nbCols;
      ColIndexes_t m_cols;
  };

  // Eigen 3.3.3
  /*
  namespace internal {
    template <typename ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols>
      struct evaluator< MatrixView <ArgType, _Rows, _Cols, _allRows, _allCols> >
      : evaluator_base< MatrixView <ArgType, _Rows, _Cols, _allRows, _allCols> >
      {
        typedef MatrixView <ArgType, _Rows, _Cols, _allRows, _allCols> XprType;
        typedef typename nested_eval<ArgType, XprType::ColsAtCompileTime>::type ArgTypeNested;
        typedef typename remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
        typedef typename XprType::CoeffReturnType CoeffReturnType;
        enum { 
          CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost,
          Flags = ArgType::Flags
        };

        evaluator(const XprType& xpr)
          : m_xpr (xpr), m_argImpl(xpr.m_arg)
        { }
        CoeffReturnType coeff(Index row, Index col) const
        {
          return m_argImpl.coeff(m_xpr.argRow(row), m_xpr.argCol(col));
        }
        ArgTypeNested& m_xpr;
        evaluator<ArgTypeNestedCleaned> m_argImpl;
      };
  }
  */
} // namespace Eigen

#include <hpp/constraints/impl/matrix-view.hh>

#endif // HPP_CONSTRAINTS_MATRIX_VIEW_HH
