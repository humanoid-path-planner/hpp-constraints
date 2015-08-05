// Copyright (c) 2014, LAAS-CNRS
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

#ifndef HPP_CONSTRAINTS_TOOL_HH
#define HPP_CONSTRAINTS_TOOL_HH

#include <Eigen/SVD>

#include "hpp/constraints/fwd.hh"

#include <hpp/model/joint.hh>
#include <hpp/model/center-of-mass-computation.hh>

namespace hpp {
  namespace constraints {
    void computeLog (vectorOut_t result, double& theta, const fcl::Matrix3f& Rerror);

    void computeJlog (const double& theta, vectorIn_t r, eigen::matrix3_t& Jlog);

    template < typename VectorType, typename MatrixType >
    static void computeCrossMatrix (const VectorType& v, MatrixType& m)
    {
      m.diagonal ().setZero ();
      m (0,1) = -v [2]; m (1,0) = v [2];
      m (0,2) = v [1]; m (2,0) = -v [1];
      m (1,2) = -v [0]; m (2,1) = v [0];
    }

    template <typename InType, typename OutType>
    static OutType convert (const InType& in, const std::size_t s) {
      OutType out(s);
      for (size_t i = 0; i < s; ++i) out[i] = in[i];
      return out;
    }

    /// \defgroup symbolic_calculus Symbolic calculus

    /// \addtogroup symbolic_calculus
    /// \{

    template <typename ValueType, typename JacobianType> class CalculusBaseAbstract;

    template <typename LhsValue, typename RhsValue> class Expression;
    template <typename LhsValue, typename RhsValue> class CrossProduct;
    template <typename LhsValue, typename RhsValue> class ScalarProduct;
    template <typename LhsValue, typename RhsValue> class Difference;
    template <typename LhsValue, typename RhsValue> class Sum;
    template <typename RhsValue> class ScalarMultiply;
    template <typename RhsValue> class RotationMultiply;
    typedef eigen::matrix3_t CrossMatrix;
    typedef Eigen::Matrix <value_type, 1, Eigen::Dynamic> RowJacobianMatrix;
    typedef Eigen::Matrix <value_type, 3, Eigen::Dynamic> JacobianMatrix;

    /// Abstract class defining a basic common interface.
    ///
    /// The purpose of this class is to allow the user to define an expression
    /// without requiring to explicitly write its type. The type will be
    /// automatically deduced by the compiler.
    ///
    /// \code
    ///   // First define a, b and c with basic elements.
    ///   CalculusBaseAbstract<OutValueType, OutJacobianType>::Ptr_t myExpression_ptr =
    ///     CalculusBaseAbstract<OutValueType, OutJacobianType>::create (a + b * c)
    /// \endcode
    template <class ValueType = eigen::vector3_t,
             class JacobianType = JacobianMatrix >
    class CalculusBaseAbstract
    {
      public:
        typedef boost::shared_ptr <CalculusBaseAbstract> Ptr_t;
        virtual const ValueType& value () const = 0;
        virtual const JacobianType& jacobian () const = 0;
        virtual void computeValue () = 0;
        virtual void computeJacobian () = 0;
        virtual void invalidate () = 0;

        template <typename Type>
        static boost::shared_ptr <Type> create (const Type& copy) {
          Type* ptr = new Type (copy);
          return boost::shared_ptr <Type> (ptr);
        }
    };

    /// Main abstract class.
    ///
    /// The framework is using CRTP to virtual function calls overload.
    /// Keep in mind, so far, no mathematical simplification is done on the
    /// expression you provide. It means that you may have a more efficient
    /// expression by calculating and simplifying the jacobian yourself.
    /// These classes provide:
    /// \li a fast way of defining new constraints,
    /// \li a robust way of calculation of jacobians -
    ///     hand calculation error proof
    /// \li a fast way to check the results of your calculations.
    template <class T,
             class ValueType = eigen::vector3_t,
             class JacobianType = JacobianMatrix,
             class CrossType = CrossMatrix >
    class CalculusBase : public CalculusBaseAbstract <ValueType, JacobianType>
    {
      public:
        CalculusBase () : cross_ (CrossMatrix::Zero()) {}

        CalculusBase (const ValueType& value, const JacobianType& jacobian) :
          value_ (value), jacobian_ (jacobian),
          cross_ (CrossMatrix::Zero()) {}

        CalculusBase (const CalculusBase& o) :
          value_ (o.value()), jacobian_ (o.jacobian_),
          cross_ (o.cross())
        {
        }

        inline const ValueType& value () const {
          return value_;
        }
        inline const JacobianType& jacobian () const {
          return jacobian_;
        }
        void computeValue () {
          if (vValid_) return;
          static_cast<T*>(this)->impl_value ();
          vValid_ = true;
        }
        void computeJacobian () {
          if (jValid_) return;
          static_cast<T*>(this)->impl_jacobian ();
          jValid_ = true;
        }
        void invalidate () {
          vValid_ = false;
          jValid_ = false;
          cValid_ = false;
        }
        inline const CrossType& cross () const {
          return cross_;
        }
        void computeCrossValue () {
          if (cValid_) return;
          computeValue ();
          computeCrossMatrix (value_, cross_);
          cValid_ = true;
        }

        template < typename RhsType >
        CrossProduct < T, RhsType > operator^ (const RhsType& rhs) const {
          return CrossProduct < T, RhsType> (static_cast<const T>(*this), rhs);
        }

        template < typename RhsType >
        Difference < T, RhsType > operator- (const RhsType& rhs) const {
          return Difference < T, RhsType> (static_cast<const T>(*this), rhs);
        }

        template < typename RhsType >
        Sum < T, RhsType > operator+ (const RhsType& rhs) const {
          return Sum < T, RhsType> (static_cast<const T>(*this), rhs);
        }

        template < typename RhsType >
        ScalarProduct < T, RhsType > operator* (const RhsType& rhs) const {
          return ScalarProduct < T, RhsType> (static_cast<const T>(*this), rhs);
        }

        RotationMultiply < T > rotate (const JointPtr_t& joint, bool transpose = false) const {
          return RotationMultiply < T > (joint, static_cast<const T>(*this), transpose);
        }

        ScalarMultiply < T > operator* (const value_type& scalar) const {
          return ScalarMultiply < T > (scalar, static_cast<const T>(*this));
        }

      protected:
        ValueType value_;
        JacobianType jacobian_;
        CrossType cross_;

        bool vValid_, jValid_, cValid_;

      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    /// Base class for classes representing an operation.
    template <typename LhsValue, typename RhsValue>
    class Expression
    {
      public:
        typedef boost::shared_ptr <
          Expression < LhsValue, RhsValue >
          > Ptr_t;
        typedef boost::weak_ptr <
          Expression < LhsValue, RhsValue >
          > WkPtr_t;

        static Ptr_t create () {
          Ptr_t p (new Expression ());
          p->init (p);
          return p;
        }

        static Ptr_t create (const LhsValue& lhs, const RhsValue& rhs) {
          Ptr_t p (new Expression (lhs, rhs));
          p->init (p);
          return p;
        }

        const LhsValue& lhs () const {
          return lhs_;
        }

        const RhsValue& rhs () const {
          return rhs_;
        }

        Expression () {}

        Expression (const Expression& other):
          rhs_ (other.rhs()), lhs_ (other.lhs())
        {}

        Expression (const LhsValue& lhs, const RhsValue& rhs):
          rhs_ (rhs), lhs_ (lhs)
        {}

        inline void init (Ptr_t self) {
          self_ = self;
        }

        RhsValue rhs_;
        LhsValue lhs_;
        WkPtr_t self_;

      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    /// Cross product of two expressions.
    template <typename LhsValue, typename RhsValue>
    class CrossProduct :
      public CalculusBase < CrossProduct < LhsValue, RhsValue > >
    {
      public:
        typedef CalculusBase < CrossProduct < LhsValue, RhsValue > >
          Parent_t;

        CrossProduct () {}

        CrossProduct (const CalculusBase <CrossProduct>& other) :
          Parent_t (other),
          e_ (static_cast <const CrossProduct&>(other).e_)
        {}

        CrossProduct (const LhsValue& lhs, const RhsValue& rhs):
          e_ (Expression < LhsValue, RhsValue >::create (lhs, rhs))
        {}

        void impl_value () {
          e_->lhs_.computeCrossValue ();
          e_->rhs_.computeValue ();
          this->value_ = e_->lhs_.cross () * e_->rhs_.value ();
        }
        void impl_jacobian () {
          e_->lhs_.computeCrossValue ();
          e_->rhs_.computeCrossValue ();
          e_->lhs_.computeJacobian ();
          e_->rhs_.computeJacobian ();
          this->jacobian_ = e_->lhs_.cross () * e_->rhs_.jacobian ()
                          - e_->rhs_.cross () * e_->lhs_.jacobian ();
        }
        void invalidate () {
          Parent_t::invalidate ();
          e_->rhs_.invalidate ();
          e_->lhs_.invalidate ();
        }

      protected:
        typename Expression < LhsValue, RhsValue >::Ptr_t e_;

        friend class Expression <LhsValue, RhsValue>;
    };

    /// Scalar product of two expressions.
    template <typename LhsValue, typename RhsValue>
    class ScalarProduct :
      public CalculusBase < ScalarProduct < LhsValue, RhsValue >, value_type, RowJacobianMatrix >
    {
      public:
        typedef CalculusBase < ScalarProduct < LhsValue, RhsValue >, value_type, RowJacobianMatrix >
          Parent_t;

        ScalarProduct () {}

        ScalarProduct (const CalculusBase <ScalarProduct>& other) :
          CalculusBase <ScalarProduct> (other),
          e_ (static_cast <const ScalarProduct&>(other).e_)
        {}

        ScalarProduct (const LhsValue& lhs, const RhsValue& rhs):
          e_ (Expression < LhsValue, RhsValue >::create (lhs, rhs))
        {}

        void impl_value () {
          e_->lhs_.computeValue ();
          e_->rhs_.computeValue ();
          this->value_ = e_->lhs_.value ().dot (e_->rhs_.value ());
        }
        void impl_jacobian () {
          e_->lhs_.computeValue ();
          e_->rhs_.computeValue ();
          e_->lhs_.computeJacobian ();
          e_->rhs_.computeJacobian ();
          this->jacobian_ = e_->lhs_.value ().transpose () * e_->rhs_.jacobian ()
                          + e_->rhs_.value ().transpose () * e_->lhs_.jacobian ();
        }
        void invalidate () {
          Parent_t::invalidate ();
          e_->rhs_.invalidate ();
          e_->lhs_.invalidate ();
        }

      protected:
        typename Expression < LhsValue, RhsValue >::Ptr_t e_;

        friend class Expression <LhsValue, RhsValue>;
    };

    /// Difference of two expressions.
    template <typename LhsValue, typename RhsValue>
    class Difference :
      public CalculusBase < Difference < LhsValue, RhsValue > >
    {
      public:
        typedef CalculusBase < Difference < LhsValue, RhsValue > >
          Parent_t;

        Difference () {}

        Difference (const CalculusBase <Difference>& other) :
          CalculusBase <Difference> (other),
          e_ (static_cast <const Difference&>(other).e_)
        {}

        Difference (const LhsValue& lhs, const RhsValue& rhs):
          e_ (Expression < LhsValue, RhsValue >::create (lhs, rhs))
        {}

        void impl_value () {
          e_->lhs_.computeValue ();
          e_->rhs_.computeValue ();
          this->value_ = e_->lhs_.value () - e_->rhs_.value ();
        }
        void impl_jacobian () {
          e_->lhs_.computeJacobian ();
          e_->rhs_.computeJacobian ();
          this->jacobian_ = e_->lhs_.jacobian () - e_->rhs_.jacobian ();
        }
        void invalidate () {
          Parent_t::invalidate ();
          e_->rhs_.invalidate ();
          e_->lhs_.invalidate ();
        }

      protected:
        typename Expression < LhsValue, RhsValue >::Ptr_t e_;

        friend class Expression <LhsValue, RhsValue>;
    };

    /// Sum of two expressions.
    template <typename LhsValue, typename RhsValue>
    class Sum :
      public CalculusBase < Sum < LhsValue, RhsValue > >
    {
      public:
        typedef CalculusBase < Sum < LhsValue, RhsValue > >
          Parent_t;

        Sum () {}

        Sum (const CalculusBase < Sum >& other) :
          CalculusBase < Sum > (other),
          e_ (static_cast <const Sum&>(other).e_)
        {}

        Sum (const RhsValue& rhs, const LhsValue& lhs):
          e_ (Expression < LhsValue, RhsValue >::create (lhs, rhs))
        {}

        void impl_value () {
          e_->lhs_.computeValue ();
          e_->rhs_.computeValue ();
          this->value_ = e_->lhs_.value () + e_->rhs_.value ();
        }
        void impl_jacobian () {
          e_->lhs_.computeJacobian ();
          e_->rhs_.computeJacobian ();
          this->jacobian_ = e_->lhs_.jacobian () + e_->rhs_.jacobian ();
        }
        void invalidate () {
          Parent_t::invalidate ();
          e_->rhs_.invalidate ();
          e_->lhs_.invalidate ();
        }

      protected:
        typename Expression < LhsValue, RhsValue >::Ptr_t e_;

        friend class Expression <LhsValue, RhsValue>;
    };

    /// Multiplication of an expression by a scalar.
    template <typename RhsValue>
    class ScalarMultiply :
      public CalculusBase < ScalarMultiply < RhsValue > >
    {
      public:
        typedef CalculusBase < ScalarMultiply < RhsValue > >
          Parent_t;

        ScalarMultiply () {}

        ScalarMultiply (const CalculusBase < ScalarMultiply >& other) :
          CalculusBase < ScalarMultiply > (other),
          e_ (static_cast <const ScalarMultiply&>(other).e_)
        {}

        ScalarMultiply (const value_type& scalar, const RhsValue& rhs):
          e_ (Expression < value_type, RhsValue >::create (scalar, rhs))
        {}

        void impl_value () {
          e_->rhs_.computeValue ();
          this->value_ = e_->lhs_ * e_->rhs_.value ();
        }
        void impl_jacobian () {
          e_->rhs_.computeJacobian ();
          this->jacobian_ = e_->lhs_ * e_->rhs_.jacobian ();
        }
        void invalidate () {
          Parent_t::invalidate ();
          e_->rhs_.invalidate ();
        }

      protected:
        typename Expression < value_type, RhsValue >::Ptr_t e_;

        friend class Expression <value_type, RhsValue>;
    };

    /// Multiplication of an expression by a rotation matrix.
    template <typename RhsValue>
    class RotationMultiply :
      public CalculusBase < RotationMultiply < RhsValue > >
    {
      public:
        typedef CalculusBase < RotationMultiply < RhsValue > >
          Parent_t;

        RotationMultiply () {}

        RotationMultiply (const CalculusBase < RotationMultiply >& other) :
          CalculusBase < RotationMultiply > (other),
          e_ (static_cast <const RotationMultiply&>(other).e_),
          transpose_ (other.transpose_)
        {}

        RotationMultiply (const JointPtr_t& joint, const RhsValue& rhs,
            bool transpose = false):
          e_ (Expression < JointPtr_t, RhsValue >::create (joint, rhs)),
          transpose_ (transpose)
        {}

        void impl_value () {
          e_->rhs_.computeValue ();
          computeRotationMatrix ();
          this->value_ = R * e_->rhs_.value ();
        }
        void impl_jacobian () {
          e_->rhs_.computeJacobian ();
          computeRotationMatrix ();
          e_->rhs_.computeCrossValue ();
          const JointJacobian_t& J = e_->lhs_->jacobian ();
          this->jacobian_ = R * (e_->rhs_.cross () * J.bottomRows (3) + e_->rhs_.jacobian ());
        }
        void invalidate () {
          Parent_t::invalidate ();
          e_->rhs_.invalidate ();
        }

      protected:
        typename Expression < JointPtr_t, RhsValue >::Ptr_t e_;

        friend class Expression <JointPtr_t, RhsValue>;

      private:
        void computeRotationMatrix () {
          const fcl::Matrix3f& Rfcl =
            e_->lhs_->currentTransformation ().getRotation ();
          if (transpose_) {
            R (0,0) = Rfcl (0,0); R (1,0) = Rfcl (0,1); R (2,0) = Rfcl (0,2);
            R (0,1) = Rfcl (1,0); R (1,1) = Rfcl (1,1); R (2,1) = Rfcl (1,2);
            R (0,2) = Rfcl (2,0); R (1,2) = Rfcl (2,1); R (2,2) = Rfcl (2,2);
          } else {
            R (0,0) = Rfcl (0,0); R (0,1) = Rfcl (0,1); R (0,2) = Rfcl (0,2);
            R (1,0) = Rfcl (1,0); R (1,1) = Rfcl (1,1); R (1,2) = Rfcl (1,2);
            R (2,0) = Rfcl (2,0); R (2,1) = Rfcl (2,1); R (2,2) = Rfcl (2,2);
          }
        }

        bool transpose_;
        eigen::matrix3_t R;
    };

    /// Basic expression representing a point in a joint frame.
    class PointInJoint : public CalculusBase <PointInJoint>
    {
      public:
        PointInJoint () {}

        PointInJoint (const CalculusBase<PointInJoint>& other) :
          CalculusBase<PointInJoint> (other)
        {
            const PointInJoint& o = static_cast <const PointInJoint&>(other);
          joint_ = o.joint ();
          local_ = o.local ();
          center_= local_.isZero ();
        }

        PointInJoint (const PointInJoint& pointInJoint) :
          CalculusBase<PointInJoint> (pointInJoint),
          joint_ (pointInJoint.joint ()), local_ (pointInJoint.local ()),
          center_ (local_.isZero ())
        {}

        PointInJoint (const JointPtr_t& joint,
            const vector3_t& pointInLocalFrame) :
          joint_ (joint), local_ (pointInLocalFrame),
          center_ (pointInLocalFrame.isZero ())
        {}
        const JointPtr_t& joint () const {
          return joint_;
        }
        const vector3_t& local () const {
          return local_;
        }
        void impl_value () {
          g_ = joint_->currentTransformation ().transform (local_);
          for (int i = 0; i < 3; ++i) this->value_[i] = g_[i];
        }
        void impl_jacobian () {
          const JointJacobian_t& j (joint_->jacobian ());
          if (!center_) {
            computeCrossRXl ();
            this->jacobian_ = - this->cross_ * j.bottomRows (3) + j.topRows (3);
          } else {
            this->jacobian_ = j.topRows (3);
          }
        }
        void computeCrossRXl () {
          if (center_) {
            this->cross_.setZero ();
            return;
          }
          computeCrossMatrix (
              joint_->currentTransformation ().getRotation () * (local_),
              this->cross_);
        }

      protected:
        JointPtr_t joint_;
        vector3_t local_;
        bool center_;

        vector3_t g_;

      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    /// Basic expression representing a vector in a joint frame.
    class VectorInJoint : public CalculusBase <VectorInJoint>
    {
      public:
        VectorInJoint () {}

        VectorInJoint (const CalculusBase<VectorInJoint>& other) :
          CalculusBase<VectorInJoint> (other),
          joint_ (static_cast <const VectorInJoint&>(other).joint()),
          vector_ (static_cast <const VectorInJoint&>(other).vector())
        {}

        VectorInJoint (const VectorInJoint& vectorInJoint) :
          CalculusBase<VectorInJoint> (vectorInJoint),
          joint_ (vectorInJoint.joint ()),
          vector_ (vectorInJoint.vector())
        {}

        VectorInJoint (const JointPtr_t& joint,
            const vector3_t& vectorInLocalFrame) :
          joint_ (joint), vector_ (vectorInLocalFrame)
        {}

        const JointPtr_t& joint () const {
          return joint_;
        }
        const vector3_t& vector () const {
          return vector_;
        }
        void impl_value () {
          g_ = joint_->currentTransformation ().getRotation () * vector_;
          for (int i = 0; i < 3; ++i) this->value_[i] = g_[i];
        }
        void impl_jacobian () {
          const JointJacobian_t& j (joint_->jacobian ());
          computeCrossRXl ();
          this->jacobian_ = - this->cross_ * j.bottomRows (3);
        }
        void computeCrossRXl () {
          computeCrossMatrix (
              joint_->currentTransformation ().getRotation () * vector_,
              this->cross_);
        }

      protected:
        JointPtr_t joint_;
        vector3_t vector_;

        vector3_t g_;

      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    /// Basic expression representing a static point
    ///
    /// Its value is constant and its jacobian is a zero matrix.
    class Point : public CalculusBase <Point, eigen::vector3_t, JacobianMatrix>
    {
      public:
        Point () {}

        Point (const CalculusBase<Point, eigen::vector3_t, JacobianMatrix>& other) :
          CalculusBase <Point, eigen::vector3_t, JacobianMatrix> (other)
        {
        }

        Point (const Point& point) :
          CalculusBase <Point, eigen::vector3_t, JacobianMatrix> (point)
        {
        }

        /// Constructor
        ///
        /// \param point the static point
        /// \param jacobianNbCols number of column of the jacobian
        Point (const vector3_t& point, size_t jacobianNbCols) :
          CalculusBase <Point, eigen::vector3_t, JacobianMatrix>
          (convert <vector3_t, eigen::vector3_t> (point, 3), JacobianMatrix::Zero (3, jacobianNbCols))
        {
        }

        void impl_value () {}
        void impl_jacobian () {}
    };

    /// Basic expression representing a COM.
    class PointCom : public CalculusBase <PointCom>
    {
      public:
        PointCom () {}

        PointCom (const CalculusBase<PointCom>& other):
          CalculusBase <PointCom> (other),
          comc_ (static_cast <const PointCom&>(other).centerOfMassComputation ())
        {
        }

        PointCom (const CenterOfMassComputationPtr_t& comc): comc_ (comc)
        {}

        inline const JacobianMatrix& jacobian () const {
          return comc_->jacobian();
        }

        const CenterOfMassComputationPtr_t& centerOfMassComputation () const {
          return comc_;
        }
        void impl_value () {
          comc_->compute (Device::COM);
          for (int i = 0; i < 3; ++i) this->value_[i] = comc_->com ()[i];
        }
        void impl_jacobian () {
          comc_->compute (Device::JACOBIAN);
          // TODO: there is memory and time to be saved here as this copy is
          // not important.
          //this->jacobian_ = comc_->jacobian ();
        }

      protected:
        CenterOfMassComputationPtr_t comc_;
    };

    /// Matrix having Expression elements
    template <typename ValueType = eigen::vector3_t, typename JacobianType = JacobianMatrix>
    class MatrixOfExpressions :
      public CalculusBase <MatrixOfExpressions <ValueType, JacobianType > ,
                           Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic >,
                           Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic > >
    {
      public:
        typedef Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic >
          Value_t;
        typedef Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic >
          Jacobian_t;
        typedef Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic >
          PseudoInv_t;
        typedef Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic >
          PseudoInvJacobian_t;
        typedef CalculusBase <MatrixOfExpressions, Value_t, Jacobian_t > Parent_t;
        typedef CalculusBaseAbstract <ValueType, JacobianType> Element_t;
        typedef typename Element_t::Ptr_t ElementPtr_t;

        MatrixOfExpressions (const Eigen::Ref<const Value_t>& value,
            const Eigen::Ref<const Jacobian_t>& jacobian) :
          Parent_t (value, jacobian),
          nRows_ (0), nCols_ (0),
          svd_ (value.rows(), value.cols(), Eigen::ComputeThinU | Eigen::ComputeThinV)
        {}

        MatrixOfExpressions (const Parent_t& other) :
          Parent_t (other),
          nRows_ (static_cast <const MatrixOfExpressions&>(other).nRows_),
          nCols_ (static_cast <const MatrixOfExpressions&>(other).nCols_),
          elements_ (static_cast <const MatrixOfExpressions&>(other).elements_),
          svd_ (static_cast <const MatrixOfExpressions&>(other).svd ())
        {
        }

        MatrixOfExpressions (const MatrixOfExpressions& matrix) :
          Parent_t (matrix),
          nRows_ (matrix.nRows_), nCols_ (matrix.nCols_),
          elements_ (matrix.elements_),
          svd_ (matrix.svd())
        {
        }

        void setSize (std::size_t nRows, std::size_t nCols) {
          nRows_ = nRows;
          nCols_ = nCols;
          elements_.resize (nRows_);
          for (std::size_t i = 0; i < nRows; ++i)
            elements_[i].resize(nCols);
        }

        ElementPtr_t& operator() (std::size_t i, std::size_t j) {
          return elements_[i][j];
        }

        void set (std::size_t i, std::size_t j, const ElementPtr_t ptr) {
          elements_[i][j] = ptr;
        }

        void impl_value () {
          size_type r = 0, c = 0, nr = 0, nc = 0;
          for (std::size_t i = 0; i < nRows_; ++i) {
            c = 0;
            nr = elements_[i][0]->value().rows();
            for (std::size_t j = 0; j < nCols_; ++j) {
              elements_[i][j]->computeValue ();
              assert (nr == elements_[i][j]->value().rows());
              nc = elements_[i][j]->value().cols();
              this->value_.block (r, c, nr, nc)
                = elements_[i][j]->value();
              c += nc;
            }
            r += nr;
          }
        }
        void impl_jacobian () {
          size_type r = 0, c = 0, nr = 0, nc = 0;
          for (std::size_t i = 0; i < nRows_; ++i) {
            c = 0;
            nr = elements_[i][0]->jacobian().rows();
            for (std::size_t j = 0; j < nCols_; ++j) {
              elements_[i][j]->computeJacobian ();
              assert (nr == elements_[i][j]->jacobian().rows());
              nc = elements_[i][j]->jacobian().cols();
              this->jacobian_.block (r, c, nr, nc)
                = elements_[i][j]->jacobian();
              c += nc;
            }
            r += nr;
          }
        }

        inline const PseudoInv_t& pinv () const {
          return pi_;
        }
        inline const PseudoInvJacobian_t& pinvJacobian () const {
          return pij_;
        }
        void computePseudoInverse () {
          this->computeValue ();
          svd_.compute (this->value_);
          Eigen::VectorXd singularValues_inv = svd_.singularValues ();
          for (typename Value_t::Index i=0; i < singularValues_inv.cols(); ++i) {
            if (i < svd_.rank ())
              singularValues_inv(i)=1.0/singularValues_inv[i];
            else
              singularValues_inv(i)=0;
          }
          pi_ = svd_.matrixV () * singularValues_inv.asDiagonal () * svd_.matrixU ().transpose();
        }
        void computePseudoInverseJacobian (const Eigen::Ref <const Eigen::Matrix<value_type, Eigen::Dynamic, 1> >& rhs) {
          this->computeJacobian ();
          computePseudoInverse ();
          Jacobian_t cache (this->jacobian_.rows(), elements_[0][0]->jacobian().cols());
          jacobianTimes (pi_ * rhs, cache);
          pij_ = - pi_ * cache;
          cache.resize (this->value_.cols(), elements_[0][0]->jacobian().cols());
          jacobianTransposeTimes (rhs - this->value_ * pi_ * rhs, cache);
          pij_ += (pi_ * pi_.transpose()) * cache;
          jacobianTransposeTimes (pi_.transpose() * pi_ * rhs , cache);
          pij_ += (Jacobian_t::Identity (pi_.rows(), this->value_.cols()) - pi_ * this->value_) * cache;
        }

        void jacobianTimes (const Eigen::Ref <const Eigen::Matrix<value_type, Eigen::Dynamic, 1> >& rhs, Eigen::Ref<Jacobian_t> cache) const {
          size_type r = 0, c = 0, nr = 0, nc = 0;
          cache.setZero();
          for (std::size_t i = 0; i < nRows_; ++i) {
            c = 0;
            nr = elements_[i][0]->jacobian().rows();
            for (std::size_t j = 0; j < nCols_; ++j) {
              elements_[i][j]->computeJacobian ();
              assert (nr == elements_[i][j]->jacobian().rows());
              nc = elements_[i][j]->jacobian().cols();
              cache.middleRows (r,nr) += this->jacobian_.block (r, c, nr, nc) * rhs[j];
              c += nc;
            }
            r += nr;
          }
        }

        void jacobianTransposeTimes (const Eigen::Ref <const Eigen::Matrix<value_type, Eigen::Dynamic, 1> >& rhs, Eigen::Ref<Jacobian_t> cache) const {
          size_type r = 0, c = 0, nr = 0, nc = 0;
          cache.setZero();
          for (std::size_t i = 0; i < nRows_; ++i) {
            c = 0;
            nr = elements_[i][0]->jacobian().rows();
            for (std::size_t j = 0; j < nCols_; ++j) {
              elements_[i][j]->computeJacobian ();
              assert (nr == elements_[i][j]->jacobian().rows());
              nc = elements_[i][j]->jacobian().cols();
              cache.row (j) += rhs.segment (r, nr).transpose() * this->jacobian_.block (r, c, nr, nc);
              c += nc;
            }
            r += nr;
          }
        }

        Eigen::JacobiSVD <Value_t>& svd () { return svd_; }

        void invalidate () {
          Parent_t::invalidate ();
          for (std::size_t i = 0; i < nRows_; ++i)
            for (std::size_t j = 0; j < nCols_; ++j)
              elements_[i][j]->invalidate ();
        }

        std::size_t nRows_, nCols_;
        std::vector <std::vector <ElementPtr_t> > elements_;

      private:
        Eigen::JacobiSVD <Value_t> svd_;
        PseudoInv_t pi_;
        PseudoInvJacobian_t pij_;

      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_TOOL_HH
