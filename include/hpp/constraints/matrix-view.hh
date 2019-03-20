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
#include <hpp/util/indent.hh>
#include <hpp/pinocchio/util.hh>
#include <hpp/constraints/fwd.hh>

# define HPP_EIGEN_USE_EVALUATOR EIGEN_VERSION_AT_LEAST(3,2,92)

namespace Eigen {

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

    /// In place addition of a segment_t to segments_t.
    static void add (segments_t& a, const segment_t& b);

    /// In place addition of segments_t to segments_t.
    static void add (segments_t& a, const segments_t& b);

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

  template <typename ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols> class MatrixBlockView;

  /// Collection of indices of matrix blocks
  /// \param _allRows whether the collection is composed of full columns
  /// \param _allCols whether the collection is composed of full rows
  ///
  /// This class enables a user to virtually create a matrix that concatenates
  /// blocks of a larger matrix.
  ///
  /// The smaller matrix is built by methods \ref lview and \ref rview
  /// \li \ref lview returns a smaller matrix that can be written in,
  /// \li \ref rview returns a smaller matrix that cannot be written in.
  template <bool _allRows = false, bool _allCols = false> class MatrixBlocks;

  /// \}

  template <bool _allRows = false, bool _allCols = false> class MatrixBlocksRef;

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
        static inline Index size() { return 1; }
        inline const Index& operator[](const Index& i) const { return i; }
      };
      template <bool If> struct get_if { template <typename T1, typename T2> static EIGEN_STRONG_INLINE T1 run (T1 then, T2 Else) { (void)Else; return then; } };
      template <> struct get_if<false> { template <typename T1, typename T2> static EIGEN_STRONG_INLINE T2 run (T1 then, T2 Else) { (void)then; return Else; } };

    template <bool _allRows, bool _allCols>
      struct traits< MatrixBlocks <_allRows, _allCols> >
    {
      enum {
        AllRows = _allRows,
        AllCols = _allCols
      };
      typedef typename internal::conditional<_allRows, internal::empty_struct, BlockIndex::segments_t>::type RowIndices_t;
      typedef typename internal::conditional<_allCols, internal::empty_struct, BlockIndex::segments_t>::type ColIndices_t;
    };

    template <bool _allRows, bool _allCols>
      struct traits< MatrixBlocksRef <_allRows, _allCols> >
    {
      enum {
        AllRows = _allRows,
        AllCols = _allCols
      };
      typedef typename internal::conditional<_allRows, internal::empty_struct, const BlockIndex::segments_t&>::type RowIndices_t;
      typedef typename internal::conditional<_allCols, internal::empty_struct, const BlockIndex::segments_t&>::type ColIndices_t;
    };

    template <typename ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols>
      struct traits< MatrixBlockView <ArgType, _Rows, _Cols, _allRows, _allCols> >
    {
# if HPP_EIGEN_USE_EVALUATOR
      typedef typename ArgType::StorageIndex StorageIndex;
# else // HPP_EIGEN_USE_EVALUATOR
      typedef typename ArgType::Index Index;
# endif // HPP_EIGEN_USE_EVALUATOR
      typedef typename traits<ArgType>::StorageKind StorageKind;
      typedef typename traits<ArgType>::XprKind XprKind;
      typedef typename ArgType::Scalar Scalar;
      enum {
# if !HPP_EIGEN_USE_EVALUATOR
        CoeffReadCost = ArgType::CoeffReadCost,
# endif // !HPP_EIGEN_USE_EVALUATOR
        Flags = ~PacketAccessBit & ~DirectAccessBit & ~ActualPacketAccessBit & ~LinearAccessBit & ArgType::Flags,
        RowsAtCompileTime = (_allRows ? ArgType::RowsAtCompileTime : _Rows),
        ColsAtCompileTime = (_allCols ? ArgType::ColsAtCompileTime : _Cols),
        MaxRowsAtCompileTime = ArgType::MaxRowsAtCompileTime,
        MaxColsAtCompileTime = ArgType::MaxColsAtCompileTime
      };
    };

# if HPP_EIGEN_USE_EVALUATOR
    template<typename Derived, typename ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols, typename Functor, typename Scalar>
    struct Assignment<Derived, MatrixBlockView <ArgType, _Rows, _Cols, _allRows, _allCols>, Functor, Dense2Dense, Scalar> {
      typedef MatrixBlockView <ArgType, _Rows, _Cols, _allRows, _allCols> OtherDerived;
      static EIGEN_STRONG_INLINE void run(Derived& dst, const OtherDerived& src, const Functor& func) {
        dst.resize(src.rows(), src.cols());
        typedef Block<Derived> BlockDerived;
        typedef Assignment<BlockDerived, typename OtherDerived::BlockConstXprType, Functor> AssignmentType;
        for (typename OtherDerived::block_iterator b (src); b.valid(); ++b) {
          BlockDerived bdst (dst.block(b.ro(), b.co(), b.rs(), b.cs()));
          AssignmentType::run(bdst, src._block(b), func);
        }
      }
    };
# else // HPP_EIGEN_USE_EVALUATOR
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
# endif // HPP_EIGEN_USE_EVALUATOR

    template <typename Src, typename Dst> struct eval_matrix_block_view_to {};
    template <typename Src, typename _ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols>
    struct eval_matrix_block_view_to <Src, MatrixBlockView<_ArgType, _Rows, _Cols, _allRows, _allCols> > {
      // MatrixBlockView <- matrix
      typedef MatrixBlockView<_ArgType, _Rows, _Cols, _allRows, _allCols> Dst;
      static void run (const Src& src, Dst& dst) {
        for (typename Dst::block_iterator b (dst); b.valid(); ++b)
          dst._block(b) = src.block(b.ro(), b.co(), b.rs(), b.cs());
      }
    };
    template <typename _ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols, typename Dst>
    struct eval_matrix_block_view_to <MatrixBlockView<_ArgType, _Rows, _Cols, _allRows, _allCols>, Dst > {
      // matrix <- MatrixBlockView
      typedef MatrixBlockView<_ArgType, _Rows, _Cols, _allRows, _allCols> Src;
      static void run (const Src& src, Dst& dst) {
        for (typename Src::block_iterator b (src); b.valid(); ++b)
          dst.block(b.ro(), b.co(), b.rs(), b.cs()) = src._block(b);
      }
    };
    template <typename _ArgType , int _Rows , int _Cols , bool _allRows , bool _allCols ,
              typename _ArgType2, int _Rows2, int _Cols2, bool _allRows2, bool _allCols2>
    struct eval_matrix_block_view_to <
      MatrixBlockView<_ArgType , _Rows , _Cols , _allRows , _allCols >,
      MatrixBlockView<_ArgType2, _Rows2, _Cols2, _allRows2, _allCols2> > {
      // MatrixBlockView <- MatrixBlockView
      typedef MatrixBlockView<_ArgType , _Rows , _Cols , _allRows , _allCols > Src;
      typedef MatrixBlockView<_ArgType2, _Rows2, _Cols2, _allRows2, _allCols2> Dst;
      static void run (const Src& src, Dst& dst) {
        typename Dst::block_iterator db (dst);
        for (typename Src::block_iterator sb (src); sb.valid(); ++sb) {
          assert (db.valid());
          dst._block(db) = src._block(sb);
          ++db;
        }
        assert (!db.valid());
      }
    };

    template <typename ReturnType, typename View, bool AllRows = View::AllRows, bool AllCols = View::AllCols>
    struct access_block_from_matrix_block_view
    {
      typedef typename View::size_type size_type;
      template <typename Derived>
      static ReturnType run (Derived& d, size_type r, size_type c, size_type rs, size_type cs)
      {
        return ReturnType (d, r, c, rs, cs);
      }
    };
    template <typename ReturnType, typename View>
    struct access_block_from_matrix_block_view <ReturnType, View, false, true>
    {
      typedef typename View::size_type size_type;
      template <typename Derived>
      static ReturnType run (Derived& d, size_type r, size_type, size_type rs, size_type)
      {
        return d.middleRows (r, rs);
      }
    };
    template <typename ReturnType, typename View>
    struct access_block_from_matrix_block_view <ReturnType, View, true, false>
    {
      typedef typename View::size_type size_type;
      template <typename Derived>
      static ReturnType run (Derived& d, size_type, size_type c, size_type, size_type cs)
      {
        return d.middleCols (c, cs);
      }
    };

    struct dont_print_indices { template <typename BlockIndexType> static void run (std::ostream&, const BlockIndexType&) {} };
    struct print_indices {
      template <typename BlockIndexType>
      static void run (std::ostream& os, const BlockIndexType& bi) {
        for (std::size_t i = 0; i < bi.size(); ++i)
          os << "[ " << bi[i].first << ", " << bi[i].second << "], ";
      }
    };

# if HPP_EIGEN_USE_EVALUATOR
    template <typename ArgType, int _Rows, int _Cols, bool _allRows, bool _allCols>
    struct unary_evaluator <MatrixBlockView <ArgType, _Rows, _Cols, _allRows, _allCols> >
    : evaluator_base <MatrixBlockView <ArgType, _Rows, _Cols, _allRows, _allCols> >
    {
      typedef MatrixBlockView <ArgType, _Rows, _Cols, _allRows, _allCols> XprType;

      enum {
        CoeffReadCost = evaluator<ArgType>::CoeffReadCost,
        Flags = ~PacketAccessBit & ~DirectAccessBit & ~ActualPacketAccessBit & ~LinearAccessBit & ArgType::Flags,
        Alignment = 0
      };
      EIGEN_DEVICE_FUNC explicit unary_evaluator (const XprType& view)
        : m_view (view)
      {}

      const XprType& m_view;
    };
# endif // HPP_EIGEN_USE_EVALUATOR
  } // namespace internal

#define EIGEN_MATRIX_BLOCKS_PUBLIC_INTERFACE(Derived) \
      enum { \
        AllRows = _allRows, \
        AllCols = _allCols \
      }; \
      typedef MatrixBlocksBase<Derived> Base; \
      typedef typename Base::size_type    size_type; \
      typedef typename Base::segments_t   segments_t; \
      typedef typename Base::segment_t    segment_t; \
      typedef typename Base::RowIndices_t RowIndices_t; \
      typedef typename Base::ColIndices_t ColIndices_t;

  /// \addtogroup hpp_constraints_tools
  /// \{

  template <typename Derived>
  class MatrixBlocksBase
  {
    public:
      enum {
        AllRows = internal::traits<Derived>::AllRows,
        AllCols = internal::traits<Derived>::AllCols,
        OneDimension = bool(AllRows) || bool(AllCols)
      };
      /// Index of vector or matrix
      typedef hpp::constraints::size_type size_type;
      /// Interval of indices [first, first + second - 1]
      typedef BlockIndex::segment_t segment_t;
      /// vector of segments
      typedef BlockIndex::segments_t segments_t;
      typedef typename internal::traits<Derived>::RowIndices_t RowIndices_t;
      typedef typename internal::traits<Derived>::ColIndices_t ColIndices_t;

      /// Smaller matrix composed by concatenation of the blocks
      template <typename MatrixType, int _Rows = MatrixType::RowsAtCompileTime, int _Cols = MatrixType::ColsAtCompileTime> struct View {
        typedef MatrixBlockView<MatrixType, _Rows, _Cols, AllRows, AllCols> type;
      }; // struct View

      Derived const& derived () const { return static_cast<Derived const&> (*this); }
      Derived      & derived ()       { return static_cast<Derived      &> (*this); }

      /// Writable view of the smaller matrix
      /// \param other matrix to whick block are extracted
      /// \return writable view of the smaller matrix composed by concatenation
      ///         of blocks.
      template <typename MatrixType>
      EIGEN_STRONG_INLINE typename View<MatrixType>::type lview(const MatrixBase<MatrixType>& other) const {
        MatrixType& o = const_cast<MatrixBase<MatrixType>&>(other).derived();
        if (Derived::OneDimension)
          return typename View<MatrixType>::type (o, nbIndices(), indices());
        else
          return typename View<MatrixType>::type (o, nbRows(), rows(), nbCols(), cols());
      }

      /// Non-writable view of the smaller matrix
      /// \param other matrix to whick block are extracted
      /// \return non-writable view of the smaller matrix composed by
      ///         concatenation of blocks.
      template <typename MatrixType>
      EIGEN_STRONG_INLINE typename View<const MatrixType>::type rview(const MatrixBase<MatrixType>& other) const {
        if (Derived::OneDimension)
          return typename View<const MatrixType>::type (other.derived(), nbIndices(), indices());
        else
          return typename View<const MatrixType>::type (other.derived(), nbRows(), rows(), nbCols(), cols());
      }

      MatrixBlocksRef<AllCols, AllRows> transpose() const
      {
        return MatrixBlocksRef<AllCols, AllRows> (nbCols(), cols(), nbRows(), rows());
      }

      MatrixBlocksRef<AllRows, true> keepRows() const
      {
        assert (!AllRows);
        return MatrixBlocksRef<AllRows, true> (nbRows(), rows());
      }

      MatrixBlocksRef<true, AllCols> keepCols() const
      {
        assert (!AllCols);
        return MatrixBlocksRef<true, AllCols> (nbCols(), cols());
      }

      /// Return row or column indices as a vector of segments
      ///
      /// \return rows indices if not all rows are selected
      ///         (see template parameter _allRows),
      ///         column indices if all rows are selected.
      inline const segments_t& indices() const
      {
        return internal::return_first<AllRows>::run(cols(), rows());
      }

      /// Return row indices
      /// \warning _allRows should be false
      inline const RowIndices_t& rows() const
      {
        return derived().rows();
      }

      /// Return column indices
      /// \warning _allCols should be false
      inline const ColIndices_t& cols() const
      {
        return derived().cols();
      }

      /// Return number of row or column indices
      ///
      /// \return number of rows indices if not all rows are selected
      ///         (see template parameter _allRows),
      ///         number of column indices if all rows are selected.
      inline const size_type& nbIndices() const
      {
        return internal::return_first<AllRows>::run(nbCols(), nbRows());
      }

      /// Return number of row indices
      /// \warning _allRows should be false
      inline const size_type& nbRows() const
      {
        return derived().nbRows();
      }

      /// Return number of column indices
      /// \warning _allCols should be false
      inline const size_type& nbCols() const
      {
        return derived().nbCols();
      }

      /// Extract a block
      /// \param i, j, ni, nj upper left corner and lengths of the block
      /// \return new instance
      MatrixBlocks <AllRows, AllCols> block
      (size_type i, size_type j, size_type ni, size_type nj) const
      {
        return MatrixBlocks <AllRows, AllCols> (BlockIndex::extract (rows (), i, ni),
                                                BlockIndex::extract (cols (), j, nj));
      }

      /// Extract a set of rows
      /// \param i, ni start and length of the set of rows
      /// \return new instance
      MatrixBlocks <AllRows, AllCols> middleRows
      (size_type i, size_type ni) const
      {
        return MatrixBlocks <AllRows, AllCols> (BlockIndex::extract (rows (), i, ni),
                                                cols ());
      }

      /// Extract a set of cols
      /// \param j, nj start and length of the set of rows
      /// \return new instance
      MatrixBlocks <AllRows, AllCols> middleCols
      (size_type j, size_type nj) const
      {
        return MatrixBlocks <AllRows, AllCols> (rows (),
                                                BlockIndex::extract (cols (), j, nj));
      }

    protected:
      /// Empty constructor
      MatrixBlocksBase () {}

      /// Copy constructor
      MatrixBlocksBase (const MatrixBlocksBase&) {}
  }; // class MatrixBlocks

  template <bool _allRows, bool _allCols>
  class MatrixBlocks : public MatrixBlocksBase < MatrixBlocks <_allRows, _allCols> >
  {
    public:
      EIGEN_MATRIX_BLOCKS_PUBLIC_INTERFACE(MatrixBlocks)

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

      /// Constructor by vectors of segments
      /// \param nbRows number of rows,
      /// \param nbCols number of columns,
      /// \param rows set of row indices,
      /// \param cols set of column indices,
      /// \warning rows and cols must be sorted
      MatrixBlocks (const size_type& nbRows, const RowIndices_t& rows,
                    const size_type& nbCols, const ColIndices_t& cols) :
        m_nbRows(nbRows), m_nbCols(nbCols), m_rows(rows), m_cols(cols)
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
      /// \param start indice for row and column
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

      /// Copy constructor
      template <typename MBDerived>
      MatrixBlocks (const MatrixBlocksBase<MBDerived>& other)
        : m_nbRows(other.nbRows())
        , m_nbCols(other.nbCols())
        , m_rows(other.rows()), m_cols(other.cols())
      {
        EIGEN_STATIC_ASSERT(
            (bool(AllRows) == bool(MBDerived::AllRows)) && (bool(AllCols) == bool(MBDerived::AllCols)),
            YOU_MIXED_MATRICES_OF_DIFFERENT_SIZES);
      }

      /// Clear rows
      inline void clearRows ()
      {
        m_rows.clear();
        m_nbRows = 0;
      }

      /// Clear cols
      inline void clearCols ()
      {
        m_cols.clear();
        m_nbCols = 0;
      }

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
      /// \tparam Sort whether set of rows should be sorted,
      /// \tparam Shrink whether set of rows should be shrunk,
      /// \tparam Cardinal whether number of rows should be recomputed
      template<bool Sort, bool Shrink, bool Cardinal>
      inline void updateRows() {
        update<Sort, Shrink, Cardinal> (m_rows, m_nbRows);
      }

      /// Selectively recompute set of columns
      /// \tparam Sort whether set of columns should be sorted,
      /// \tparam Shrink whether set of columns should be shrunk,
      /// \tparam Cardinal whether number of columns should be recomputed
      template<bool Sort, bool Shrink, bool Cardinal>
      inline void updateCols() {
        update<Sort, Shrink, Cardinal> (m_cols, m_nbCols);
      }

      /// Return row indices
      /// \warning _allRows should be false
      inline const RowIndices_t& rows() const
      {
        return m_rows;
      }

      /// Return column indices
      /// \warning _allCols should be false
      inline const ColIndices_t& cols() const
      {
        return m_cols;
      }

      /// Return number of row indices
      /// \warning _allRows should be false
      inline const size_type& nbRows() const
      {
        return m_nbRows;
      }

      /// Return number of column indices
      /// \warning _allCols should be false
      inline const size_type& nbCols() const
      {
        return m_nbCols;
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

  /// \cond
  template <bool _allRows, bool _allCols>
  class MatrixBlocksRef : public MatrixBlocksBase < MatrixBlocksRef <_allRows, _allCols> >
  {
    public:
      EIGEN_MATRIX_BLOCKS_PUBLIC_INTERFACE(MatrixBlocksRef)

      /// Constructor by vectors of segments
      /// \param rows set of row indices,
      /// \param cols set of column indices,
      /// \warning rows and cols must be sorted
      MatrixBlocksRef (const size_type& nbRows, const RowIndices_t& rows,
                       const size_type& nbCols, const ColIndices_t& cols) :
        m_nbRows(nbRows), m_nbCols(nbCols),
        m_rows  (rows  ), m_cols  (cols  )
      {}

      /// Constructor by two (row or col) MatrixBlocks
      template <typename Derived1, typename Derived2>
      MatrixBlocksRef (const MatrixBlocksBase<Derived1>& rows,
                       const MatrixBlocksBase<Derived2>& cols) :
        m_nbRows(rows.nbIndices()), m_nbCols(cols.nbIndices()),
        m_rows  (rows.  indices()), m_cols  (cols.  indices())
      {
        EIGEN_STATIC_ASSERT( bool(Derived1::OneDimension) && bool (Derived2::OneDimension),
            YOU_MIXED_MATRICES_OF_DIFFERENT_SIZES);
      }

      /// Constructor by a collection of indices
      /// \param idx collections of indices (for rows and columns)
      /// \warning idx must be sorted and shrinked
      /// \note if all rows or all columns are selected (template parameter)
      ///       the block will contain all rows, respectively all columns.
      MatrixBlocksRef (const size_type& nidx, const segments_t& idx)
        : m_nbRows(_allRows ? 0 : nidx)
        , m_nbCols(_allCols ? 0 : nidx)
        , m_rows(idx), m_cols(idx)
      {}

      /// Copy constructor
      MatrixBlocksRef (const MatrixBlocksRef& other)
        : Base (other)
        , m_nbRows(other.nbRows())
        , m_nbCols(other.nbCols())
        , m_rows(other.rows())
        , m_cols(other.cols())
      {}

      /// Return row indices
      /// \warning _allRows should be false
      inline const RowIndices_t& rows() const
      {
        return m_rows;
      }

      /// Return column indices
      /// \warning _allCols should be false
      inline const ColIndices_t& cols() const
      {
        return m_cols;
      }

      /// Return number of row indices
      /// \warning _allRows should be false
      inline const size_type& nbRows() const
      {
        return m_nbRows;
      }

      /// Return number of column indices
      /// \warning _allCols should be false
      inline const size_type& nbCols() const
      {
        return m_nbCols;
      }

      const size_type m_nbRows, m_nbCols;
      RowIndices_t m_rows;
      ColIndices_t m_cols;
  }; // class MatrixBlocksRef
  /// \endcond

  template <typename Derived>
  std::ostream& operator<< (std::ostream& os, const MatrixBlocksBase<Derived>& mbi)
  {
    typedef typename internal::conditional<Derived::AllRows, internal::dont_print_indices, internal::print_indices>::type row_printer;
    typedef typename internal::conditional<Derived::AllCols, internal::dont_print_indices, internal::print_indices>::type col_printer;
    if (!Derived::AllRows) {
      os << "Rows: ";
      row_printer::run (os, mbi.rows());
      if (!Derived::AllCols) os << hpp::iendl;
    }
    if (!Derived::AllCols) {
      os << "Cols: ";
      col_printer::run (os, mbi.cols());
    }
    return os;
  }

  typedef Eigen::MatrixBlocks<false, true> RowBlockIndices;
  typedef Eigen::MatrixBlocks<true, false> ColBlockIndices;

  /// A view of an Eigen matrix.
  ///
  /// Instances of MatrixBlockView are easily built from a MatrixBlocks object.
  ///
  /// As of the date of this documentation, this class does not support all
  /// operations on matrices.
  /// See tests/matrix-view.cc for a list of supported features.
  ///
  /// Although it is usually not useful to iterate over the blocks, it is
  /// possible to do it as follows:
  /// \code
  /// MatrixBlockView view = ...;
  /// for (MatrixBlockView::block_iterator block (view); block.valid(); ++block)
  /// {
  ///   // Access the current block
  ///   view._block(block); // Returns an Eigen::Block object.
  ///   // For the current block
  ///   block.ri() // row, input
  ///   block.ci() // col, input
  ///   block.ro() // row, output
  ///   block.co() // col, output
  ///   block.rs() // number of rows
  ///   block.cs() // number of cols
  /// }
  /// \endcode
  /// \sa MatrixBlocks, MatrixBlockView::block_iterator
  template <typename _ArgType, int _Rows = _ArgType::RowsAtCompileTime, int _Cols = _ArgType::ColsAtCompileTime, bool _allRows = false, bool _allCols = false>
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
      struct block_iterator {
        const MatrixBlockView& view;
        size_type  row,  col;
        internal::variable_if_dynamic<size_type, (_allRows ? 0 : Dynamic) > _ro;
        internal::variable_if_dynamic<size_type, (_allCols ? 0 : Dynamic) > _co;
        block_iterator (const MatrixBlockView& v) : view(v), row (0), col (0), _ro(0), _co(0) {}
        /// <b>R</b>ow in the <b>O</b>utput matrix
        size_type ro() const { return _ro.value(); }
        /// <b>C</b>ol in the <b>O</b>utput matrix
        size_type co() const { return _co.value(); }
        /// <b>R</b>ow in the <b>I</b>nput matrix
        size_type ri() const { return internal::get_if<AllRows>::run(std::make_pair(0,view.m_nbRows), view.m_rows[row]).first; }
        /// <b>C</b>ol in the <b>I</b>nput matrix
        size_type ci() const { return internal::get_if<AllCols>::run(std::make_pair(0,view.m_nbCols), view.m_cols[col]).first; }
        /// number of <b>R</b>ow<b>S</b>
        size_type rs() const { return internal::get_if<AllRows>::run(std::make_pair(0,view.m_nbRows), view.m_rows[row]).second; }
        /// number of <b>C</b>ol<b>S</b>
        size_type cs() const { return internal::get_if<AllCols>::run(std::make_pair(0,view.m_nbCols), view.m_cols[col]).second; }
        // ++it
        block_iterator& operator++()
        {
          _ro.setValue(_ro.value() + rs());
          ++row;
          if (row == (size_type)view.m_rows.size()) {
            row = 0;
            _ro.setValue(0);
            _co.setValue(_co.value() + cs());
            ++col;
            // if (col < (size_type)view.m_cols.size()) _co.setValue(0);
          }
          return *this;
        };
        // it++
        block_iterator operator++(int) { block_iterator copy(*this); operator++(); return copy; };
        bool valid () const {
          return (internal::get_if<AllRows>::run(true, view.m_rows.size()>0))
            && (col < (size_type)internal::get_if<AllCols>::run(1, view.m_cols.size()));
        }
      };
      typedef MatrixBase< MatrixBlockView<_ArgType, _Rows, _Cols, _allRows, _allCols> > Base;
      EIGEN_GENERIC_PUBLIC_INTERFACE(MatrixBlockView)

      typedef Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> PlainObject;
      // typedef typename internal::ref_selector<MatrixBlockView>::type Nested;
      typedef _ArgType ArgType;
      typedef typename internal::ref_selector<ArgType>::type ArgTypeNested;
      typedef typename internal::remove_all<ArgType>::type NestedExpression;
      // typedef typename Base::CoeffReturnType CoeffReturnType;
      // typedef typename Base::Scalar Scalar;

      template <typename Derived>
      struct block_t {
        typedef Block<Derived,
                (AllRows ? Derived::RowsAtCompileTime : Eigen::Dynamic),
                (AllCols ? Derived::ColsAtCompileTime : Eigen::Dynamic),
                (AllCols ? (bool)Derived::IsRowMajor
                 : (AllRows ? (bool)!Derived::IsRowMajor : false)
                )> type ;
      };
      typedef typename block_t<ArgType>::type       BlockXprType;
      typedef typename block_t<const ArgType>::type BlockConstXprType;

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
        internal::eval_matrix_block_view_to<MatrixBlockView, Dest>::run (*this, dst);
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
        internal::eval_matrix_block_view_to<OtherDerived, MatrixBlockView>::run (other.derived(), *this);
        return *this;
      }

      EIGEN_STRONG_INLINE size_type _blocks() const { return m_rows.size() * m_cols.size(); }
      EIGEN_STRONG_INLINE BlockXprType _block(const block_iterator& b)
      {
        return internal::access_block_from_matrix_block_view< BlockXprType, MatrixBlockView >
          ::template run <ArgType> (m_arg, b.ri(), b.ci(), b.rs(), b.cs());
      }
      EIGEN_STRONG_INLINE const BlockConstXprType _block(const block_iterator& b) const
      {
        return internal::access_block_from_matrix_block_view< const BlockConstXprType, MatrixBlockView >
          ::template run <const ArgType> (m_arg, b.ri(), b.ci(), b.rs(), b.cs());
      }
      EIGEN_STRONG_INLINE block_iterator _block_iterator() const { return block_iterator(*this); }

      EIGEN_STRONG_INLINE bool isZero (const RealScalar& prec = NumTraits<Scalar>::dummy_precision()) const
      {
        for (block_iterator block (*this); block.valid(); ++block)
          if (!m_arg.block(
                block.ri(), block.ci(),
                block.rs(), block.cs())
              .isZero(prec))
            return false;
        return true;
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
#include <hpp/constraints/impl/matrix-view-operation.hh>

# undef HPP_EIGEN_USE_EVALUATOR

namespace hpp {
  template <int Option>
  struct prettyPrint<constraints::segment_t, Option> {
    static std::ostream& run (std::ostream& os, const constraints::segment_t& s) {
      return os << "[ " << s.first << ", " << s.first + s.second << " ]";
    }
  };
}

#endif // HPP_CONSTRAINTS_MATRIX_VIEW_HH
