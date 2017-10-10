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
#include <hpp/constraints/fwd.hh>

namespace Eigen {

  template <typename ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols> class MatrixBlockView;

  namespace internal {
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
        template <typename In0_t, typename In1_t> empty_struct (In0_t, In1_t) {}
        static inline Index size() { return 0; }
        inline const Index& operator[](const Index& i) const { return i; }
      };

    template <typename ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols>
      struct traits< MatrixBlockView <ArgType, _Rows, _Cols, _allRows, _allCols> >
    {
      typedef typename ArgType::Index Index;
      typedef Eigen::Dense StorageKind;
      typedef Eigen::MatrixXpr XprKind;
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

    template <typename Other, typename View, bool AllCols = View::AllCols> struct evalCols {
      static inline void run (Other& dst, const View& src, const typename Other::Index& row)
      {
        dst.derived() = src.derived().m_arg.middleRows(row, dst.rows());
      }
      static inline void write (const Other& src, View& dst, const typename Other::Index& row)
      {
        dst.m_arg.middleRows(row, src.rows()) = src;
      }
    };
    template <typename Other, typename View> struct evalCols <Other, View, false> {
      static inline void run (Other& dst, const View& src, const typename Other::Index& row)
      {
        std::size_t col = 0;
        for (std::size_t j = 0; j < src.m_cols.size(); ++j) {
          dst.derived().middleCols(col, src.m_cols[j].second) =
            src.m_arg.block(row,        src.m_cols[j].first,
                            dst.rows(), src.m_cols[j].second);
          col += src.m_cols[j].second;
        }
      }
      static inline void write (const Other& src, View& dst, const typename Other::Index& row)
      {
        std::size_t col = 0;
        for (std::size_t j = 0; j < dst.m_cols.size(); ++j) {
          dst.m_arg.block(row,        dst.m_cols[j].first,
              src.rows(), dst.m_cols[j].second)
            = src.derived().middleCols(col, dst.m_cols[j].second);
          col += dst.m_cols[j].second;
        }
      }
    };
    template <typename Other, typename View, bool AllRows = View::AllRows> struct evalRows {
      static inline void run (Other& dst, const View& src)
      {
        evalCols<Other, View>::run(dst.derived(), src, 0);
      }
      static inline void write (const Other& src, View& dst)
      {
        evalCols<Other, View>::write(src, dst, 0);
      }
    };
    template <typename Other, typename View> struct evalRows <Other, View, false> {
      static inline void run (Other& dst, const View& src)
      {
        std::size_t row = 0;
        for (std::size_t i = 0; i < src.m_rows.size(); ++i) {
          typedef typename Other::RowsBlockXpr Rows_t;
          Rows_t rows = dst.middleRows(row, src.m_rows[i].second);
          evalCols<Rows_t, View>::run(rows, src, src.m_rows[i].first);
          row += src.m_rows[i].second;
        }
      }
      static inline void write (const Other& src, View& dst)
      {
        std::size_t row = 0;
        for (std::size_t i = 0; i < dst.m_rows.size(); ++i) {
          typedef typename Other::ConstRowsBlockXpr ConstRows_t;
          ConstRows_t rows = src.middleRows(row, dst.m_rows[i].second);
          evalCols<ConstRows_t, View>::write(rows, dst, dst.m_rows[i].first);
          row += dst.m_rows[i].second;
        }
      }
    };

    template <bool print> struct print_indices { template <typename BlockIndexType> static void run (std::ostream&, const BlockIndexType&) {} };
    template <> struct print_indices <true> {
      template <typename BlockIndexType>
      static void run (std::ostream& os, const BlockIndexType& bi) {
        for (std::size_t i = 0; i < bi.size(); ++i)
          os << "[ " << bi[i].first << ", " << bi[i].second << "], ";
      }
    };
  } // namespace internal

  /// \addtogroup hpp_constraints_tools
  /// \{

  /// List of integer intervals
  ///
  /// Used to select blocks in a vector or in a matrix.
  struct BlockIndex {
    /// Index of vector or matrix
    typedef hpp::constraints::size_type size_type;
    /// Interval of indices [first, first + second - 1]
    typedef std::pair<size_type, size_type> segment_t;
    /// vector of segments
    typedef std::vector<segment_t> segments_t;

    /// Return the number of indices in the vector of segments.
    /// \param a vector of segments
    static size_type cardinal (const segments_t& a);

    /// Build a vector of segments from an array of Boolean.
    /// \param array array of Boolean values
    /// \return the vector of segments corresponding to true values in the
    ///         input.
    template <typename Derived>
    static segments_t fromLogicalExpression
    (const Eigen::ArrayBase<Derived>& array);

    /// Sort segments in increasing order.
    /// Compare lower bounds of intervals and lengths if lower bounds are equal.
    static void sort   (segments_t& a);

    /// Build a sequence of non overlapping segments.
    /// \param a a vector of segments
    /// \note assumes a is sorted
    static void shrink (segments_t& a);

    /// Whether two segments overlap.
    static bool overlap (const segment_t& a, const segment_t& b);

    /// Compute the union of tws segments.
    static segments_t sum (const segment_t& a, const segment_t& b);

    /// Compute the set difference between two segments.
    static segments_t difference (const segment_t& a, const segment_t& b);

    /// Compute the set difference between a vector of segments and a segment.
    /// \note assumes a is sorted
    static segments_t difference (const segments_t& a, const segment_t& b);

    /// Compute the set difference between a segment and a vector of segments.
    /// \note assume b is sorted
    static segments_t difference (const segment_t& a, const segments_t& b);

    /// Compute the set difference between two vectors of segments.
    /// \note assume a and b are sorted
    static segments_t difference (const segments_t& a, const segments_t& b);

    /// Split a set of segment into two sets of segments
    /// \param segments input set of segments,
    /// \param cardinal cardinal of the first set of segments,
    /// \return the first set of segments.
    ///
    /// The second set is stored in the input set of segments.
    static segments_t split (segments_t& segments, const size_type& cardinal);

    /// Extract a subset of a set of segments
    /// \param segments input set of segments
    /// \param start beginning of extracted set of segments (cardinal of subset
    ///        left behind in input set of segments)
    /// \param cardinal cardinal of extracted set of segments,
    /// \return subset of segments.
    static segments_t extract (const segments_t& segments, size_type start,
                               size_type cardinal);
  }; // struct BlockIndex

  /// Collection of indices of matrix blocks
  /// \param _allRows whether the collection is composed of full columns
  /// \param _allCols whether the collection is composed of full rows
  ///
  /// This class enables a user to virtually create a matrix that concatenates
  /// blocks of a larger matrix.
  ///
  /// The smaller matrix is built by methods lview and rview
  /// \li lview returns a smaller matrix that can be written in,
  /// \li rview returns a smaller matrix that cannot be written in.
  template <bool _allRows, bool _allCols>
  class MatrixBlocks
  {
    public:
      /// Index of vector or matrix
      typedef hpp::constraints::size_type size_type;
      /// Interval of indices [first, first + second - 1]
      typedef BlockIndex::segment_t segment_t;
      /// vector of segments
      typedef BlockIndex::segments_t segments_t;
      typedef typename internal::conditional<_allRows, internal::empty_struct, segments_t>::type RowIndices_t;
      typedef typename internal::conditional<_allCols, internal::empty_struct, segments_t>::type ColIndices_t;

      /// Smaller matrix composed by concatenation of the blocks
      template <typename Derived, int _Rows, int _Cols> struct View {
        typedef MatrixBlockView<Derived, _Rows, _Cols, _allRows, _allCols> type;
        typedef MatrixBlockView<Derived, _Rows, _Cols, _allCols, _allRows> transpose_type;
      }; // struct View

      /// Empty constructor
      MatrixBlocks () : m_nbRows(0), m_nbCols(0), m_rows(), m_cols() {}

      /// Constructor by vectors of segments
      /// \param rows set of row indices,
      /// \param cols set of column indices,
      /// \warning rows and cols must be sorted
      MatrixBlocks (const segments_t& rows,
                          const segments_t& cols) :
        m_nbRows(BlockIndex::cardinal(rows)),
        m_nbCols(BlockIndex::cardinal(cols)), m_rows(rows), m_cols(cols)
      {
# ifndef NDEBUG
        // test that input is sorted
        segments_t r (rows); BlockIndex::sort (r);
        assert (r == rows);
        segments_t c (cols); BlockIndex::sort (c);
        assert (c == cols);
#endif
      }

      /// Constructor of single block
      /// \param first indice for row and column
      /// \param size number of indices in the block (row and column)
      /// \note if all rows or all columns are selected (template parameter)
      ///       the block will contain all rows, respectively all columns.
      MatrixBlocks (size_type start, size_type size)
        : m_nbRows(_allRows ? 0 : size)
        , m_nbCols(_allCols ? 0 : size)
        , m_rows(1, BlockIndex::segment_t (start, size))
        , m_cols(1, BlockIndex::segment_t (start, size))
      {}

      /// Constructor by a collection of indices
      /// \param idx collections of indices (for rows and columns)
      /// \warning idx must be sorted and shrinked
      /// \note if all rows or all columns are selected (template parameter)
      ///       the block will contain all rows, respectively all columns.
      MatrixBlocks (const segments_t& idx)
        : m_nbRows(_allRows ? 0 : BlockIndex::cardinal(idx))
        , m_nbCols(_allCols ? 0 : BlockIndex::cardinal(idx))
        , m_rows(idx), m_cols(idx)
      {}

      /// Constructor of a single block
      /// \param idx segment of row and column indices
      /// \note if all rows or all columns are selected (template parameter)
      ///       the block will contain all rows, respectively all columns.
      MatrixBlocks (const segment_t& idx)
        : m_nbRows(_allRows ? 0 : idx.second)
        , m_nbCols(_allCols ? 0 : idx.second)
        , m_rows(segments_t(1,idx)), m_cols(segments_t(1,idx))
      {}

      /// Constructor from other MatrixBlocks
      /// \note This constructor will only be called when
      /// \code
      /// MatrixBlocks<true, false> ( MatrixBlocks<false, true> (...));
      /// MatrixBlocks<false, true> ( MatrixBlocks<true, false> (...));
      /// \endcode
      template <bool _otherAllRows, bool _otherAllCols>
      MatrixBlocks (const MatrixBlocks<_otherAllRows,_otherAllCols>& other)
        : m_nbRows(other.m_nbCols)
        , m_nbCols(other.m_nbRows)
        , m_rows(other.m_cols), m_cols(other.m_rows)
      {
        assert((_allRows != _allCols)
            && (_otherAllRows == _allCols)
            && (_otherAllCols == _allRows));
      }

      /// Copy constructor
      MatrixBlocks (const MatrixBlocks<_allRows,_allCols>& other)
        : m_nbRows(other.m_nbRows)
        , m_nbCols(other.m_nbCols)
        , m_rows(other.m_rows), m_cols(other.m_cols)
      {}

      /// Add consecutive rows
      /// \param row first row to add
      /// \param size number of rows to add
      inline void addRow (const size_type& row, const size_type size)
      {
        m_rows.push_back(segment_t (row, size));
        m_nbRows += size;
      }

      /// Add consecutive columns
      /// \param col first column to add
      /// \param size number of columns to add
      inline void addCol (const size_type& col, const size_type size)
      {
        m_cols.push_back(segment_t (col, size));
        m_nbCols += size;
      }

      /// Selectively recompute set of rows
      /// \param Sort whether set of rows should be sorted,
      /// \param Shrink whether set of rows should be shrunk,
      /// \param Cardinal whether number of rows should be recomputed
      template<bool Sort, bool Shrink, bool Cardinal>
      inline void updateRows() {
        update<Sort, Shrink, Cardinal> (m_rows, m_nbRows);
      }

      /// Selectively recompute set of columns
      /// \param Sort whether set of columns should be sorted,
      /// \param Shrink whether set of columns should be shrunk,
      /// \param Cardinal whether number of columns should be recomputed
      template<bool Sort, bool Shrink, bool Cardinal>
      inline void updateCols() {
        update<Sort, Shrink, Cardinal> (m_cols, m_nbCols);
      }

      /// Writable view of the smaller matrix
      /// \param other matrix to whick block are extracted
      /// \return writable view of the smaller matrix composed by concatenation
      ///         of blocks.
      template <typename Derived>
      EIGEN_STRONG_INLINE typename View<Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::type lview(const MatrixBase<Derived>& other) const {
        Derived& o = const_cast<MatrixBase<Derived>&>(other).derived();
        if (_allCols || _allRows)
          return typename View<Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::type (o, this->nbIndices(), this->indices());
        else
          return typename View<Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::type (o, m_nbRows, m_rows, m_nbCols, m_cols);
      }

      /// Writable view of the smaller matrix transposed
      /// \param other matrix to whick block are extracted
      /// \return writable view of the smaller matrix composed by concatenation
      ///         of blocks and transposed.
      template <typename Derived>
      EIGEN_STRONG_INLINE typename View<Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::transpose_type lviewTranspose(const MatrixBase<Derived>& other) const {
        Derived& o = const_cast<MatrixBase<Derived>&>(other).derived();
        if (_allCols || _allRows)
          return typename View<Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::transpose_type (o, nbIndices(), indices());
        else
          return typename View<Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::transpose_type (o, m_nbCols, m_cols, m_nbRows, m_rows);
      }

      /// Non-writable view of the smaller matrix
      /// \param other matrix to whick block are extracted
      /// \return non-writable view of the smaller matrix composed by
      ///         concatenation of blocks.
      template <typename Derived>
      EIGEN_STRONG_INLINE typename View<const Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::type rview(const MatrixBase<Derived>& other) const {
        if (_allCols || _allRows)
          return typename View<const Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::type (other.derived(), nbIndices(), indices());
        else
          return typename View<const Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::type (other.derived(), m_nbRows, m_rows, m_nbCols, m_cols);
      }

      /// Non-writable view of the smaller matrix transposed
      /// \param other matrix to whick block are extracted
      /// \return non-writable view of the smaller matrix composed by
      ///         concatenation of blocks and transposed.
      template <typename Derived>
      EIGEN_STRONG_INLINE typename View<const Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::transpose_type rviewTranspose(const MatrixBase<Derived>& other) const {
        if (_allCols || _allRows)
          return typename View<const Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::transpose_type (other.derived(), nbIndices(), indices());
        else
          return typename View<const Derived, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>::transpose_type (other.derived(), m_nbCols, m_cols, m_nbRows, m_rows);
      }

      /// Return row or column indices as a vector of segments
      ///
      /// \return rows indices if not all rows are selected
      ///         (see template parameter _allRows),
      ///         column indices if all rows are selected.
      inline const segments_t& indices() const
      {
        return internal::return_first<_allRows>::run(m_cols, m_rows);
      }

      /// Return row indices
      /// \assertion _allRows should be false
      inline const RowIndices_t& rows() const
      {
        assert (!_allRows);
        return m_rows;
      }

      /// Return column indices
      /// \assertion _allCols should be false
      inline const ColIndices_t& cols() const
      {
        assert (!_allCols);
        return m_cols;
      }

      /// Return number of row or column indices
      ///
      /// \return number of rows indices if not all rows are selected
      ///         (see template parameter _allRows),
      ///         number of column indices if all rows are selected.
      inline const size_type& nbIndices() const
      {
        return internal::return_first<_allRows>::run(m_nbCols, m_nbRows);
      }

      /// Return number of row indices
      /// \assertion _allRows should be false
      inline const size_type& nbRows() const
      {
        assert (_allRows);
        return m_nbRows;
      }

      /// Return number of column indices
      /// \assertion _allCols should be false
      inline const size_type& nbCols() const
      {
        assert (_allCols);
        return m_nbCols;
      }

      /// Extract a block
      /// \param i, j, ni, nj upper left corner and lengths of the block
      /// \return new instance
      block (size_type i, size_type j, size_type ni, size_type nj) const
      {
        return MatrixBlock (rows ().extract (i, ni), cols ().extract (j, nj));
      }

      template<bool Sort, bool Shrink, bool Cardinal>
      inline void updateIndices() {
        update<Sort, Shrink, Cardinal> (
            internal::return_first<_allRows>::run(m_cols  , m_rows  ),
            internal::return_first<_allRows>::run(m_nbCols, m_nbRows));
      }

      size_type m_nbRows, m_nbCols;
      RowIndices_t m_rows;
      ColIndices_t m_cols;

    private:
      template<bool Sort, bool Shrink, bool Cardinal>
      static inline void update(segments_t& b, size_type& idx) {
        if (Sort)     BlockIndex::sort(b);
        if (Shrink)   BlockIndex::shrink(b);
        if (Cardinal) idx = BlockIndex::cardinal(b);
      }
  }; // class MatrixBlocks

  template <bool _allRows, bool _allCols>
  std::ostream& operator<< (std::ostream& os, MatrixBlocks<_allRows, _allCols> mbi)
  {
    if (!_allRows) {
      os << "Rows: ";
      internal::print_indices<!_allRows>::run (os, mbi.m_rows);
      if (!_allCols) os << '\n';
    }
    if (!_allCols) {
      os << "Cols: ";
      internal::print_indices<!_allCols>::run (os, mbi.m_cols);
    }
    return os;
  }

  typedef Eigen::MatrixBlocks<false, true> RowBlockIndices;
  typedef Eigen::MatrixBlocks<true, false> ColBlockIndices;

  template <typename _ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols>
  class MatrixBlockView : public MatrixBase< MatrixBlockView<_ArgType, _Rows, _Cols, _allRows, _allCols> >
  {
    public:
    typedef hpp::constraints::size_type size_type;
      enum {
        Rows = _Rows,
        Cols = _Cols,
        AllRows = _allRows,
        AllCols = _allCols
      };
      typedef MatrixBase< MatrixBlockView<_ArgType, _Rows, _Cols, _allRows, _allCols> > Base;
      EIGEN_GENERIC_PUBLIC_INTERFACE(MatrixBlockView)

      typedef Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> PlainObject;
      // typedef typename internal::ref_selector<MatrixBlockView>::type Nested;
      typedef _ArgType ArgType;
      typedef typename internal::ref_selector<ArgType>::type ArgTypeNested;
      // typedef typename Base::CoeffReturnType CoeffReturnType;
      // typedef typename Base::Scalar Scalar;

      typedef MatrixBlocks<_allRows, _allCols> MatrixIndices_t;
      typedef typename MatrixIndices_t::segments_t Indices_t;
      typedef typename internal::conditional<_allRows, const internal::empty_struct, const Indices_t& >::type RowIndices_t;
      typedef typename internal::conditional<_allCols, const internal::empty_struct, const Indices_t& >::type ColIndices_t;

      // using Base::operator=;

      MatrixBlockView (ArgType& arg, const size_type& nbRows,
                       const RowIndices_t rows, const size_type& nbCols,
                       const ColIndices_t cols) :
        m_arg (arg), m_nbRows(nbRows), m_rows(rows), m_nbCols(nbCols),
        m_cols(cols)
      {
      }

      /// Valid only when _allRows or _allCols is true
      MatrixBlockView (ArgType& arg, const size_type& nbIndices,
                       const Indices_t& indices) :
        m_arg (arg), m_nbRows(_allRows ? arg.rows() : nbIndices),
        m_rows(indices), m_nbCols(_allCols ? arg.cols() : nbIndices),
        m_cols(indices)
      {}

      EIGEN_STRONG_INLINE size_type rows() const { return m_nbRows; }
      EIGEN_STRONG_INLINE size_type cols() const { return m_nbCols; }

      EIGEN_STRONG_INLINE CoeffReturnType coeff (size_type index) const
      {
        assert(false && "It is not possible to access the coefficients of "
               "MatrixBlockView this way.");
      }
      EIGEN_STRONG_INLINE CoeffReturnType coeff (size_type row, size_type col)
        const
      {
        assert(false && "It is not possible to access the coefficients of "
               "MatrixBlockView this way.");
      }
      EIGEN_STRONG_INLINE Scalar& coeffRef (size_type index)
      {
        assert(false && "It is not possible to access the coefficients of "
               "MatrixBlockView this way.");
      }
      EIGEN_STRONG_INLINE Scalar& coeffRef (size_type row, const size_type& col)
      {
        assert(false && "It is not possible to access the coefficients of "
               "MatrixBlockView this way.");
      }
      template <typename Dest>
      EIGEN_STRONG_INLINE void evalTo (Dest& dst) const {
        internal::evalRows<Dest, MatrixBlockView>::run(dst, *this);
      }

      template <typename Dest>
      EIGEN_STRONG_INLINE void writeTo (Dest& dst) const {
        dst.resize(rows(), cols());
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
        internal::evalRows<const OtherDerived, MatrixBlockView>::write(other.derived(), *this);
        return *this;
      }

      ArgType& m_arg;
      size_type m_nbRows;
      RowIndices_t m_rows;
      size_type m_nbCols;
      ColIndices_t m_cols;
  }; // MatrixBlockView

  ///\}

} // namespace Eigen

#include <hpp/constraints/impl/matrix-view.hh>

#endif // HPP_CONSTRAINTS_MATRIX_VIEW_HH
