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

#ifndef HPP_CONSTRAINTS_SYMBOLIC_CALCULUS_HH
#define HPP_CONSTRAINTS_SYMBOLIC_CALCULUS_HH

#define HPP_CONSTRAINTS_CB_REF boost::shared_ptr
#define HPP_CONSTRAINTS_CB_WKREF boost::shared_ptr

#define HPP_CONSTRAINTS_CB_DEFINE_OPERATOR1(op, InType, OutType) \
    template < typename RhsType > \
        typename Traits <OutType < RhsType > >::Ptr_t op ( \
            const InType& lhs, \
            const HPP_CONSTRAINTS_CB_REF <RhsType>& rhs) { \
          typedef OutType < RhsType> Op_t; \
          return Op_t::create (lhs, rhs); \
        }
#define HPP_CONSTRAINTS_CB_FRIEND_OPERATOR1(op, InType, OutType) \
  template < typename RhsType > \
        friend typename Traits <OutType < RhsType > >::Ptr_t op ( \
            const InType& lhs, \
            const HPP_CONSTRAINTS_CB_REF <RhsType>& rhs);

#define HPP_CONSTRAINTS_CB_DEFINE_OPERATOR2(op, OutType) \
    template < typename LhsType, typename RhsType > \
        typename Traits <OutType < LhsType, RhsType > >::Ptr_t op ( \
            const HPP_CONSTRAINTS_CB_REF <LhsType>& lhs, \
            const HPP_CONSTRAINTS_CB_REF <RhsType>& rhs) { \
          typedef OutType < LhsType, RhsType> Op_t; \
          return Op_t::create (lhs, rhs); \
        }

#define HPP_CONSTRAINTS_CB_FRIEND_OPERATOR2(op, OutType) \
    template < typename LhsType, typename RhsType > \
        friend typename Traits <OutType < LhsType, RhsType > >::Ptr_t op ( \
            const HPP_CONSTRAINTS_CB_REF <LhsType>& lhs, \
            const HPP_CONSTRAINTS_CB_REF <RhsType>& rhs);

#define HPP_CONSTRAINTS_CB_CREATE1(Class, Arg0Type) \
  static typename Traits <Class>::Ptr_t create (Arg0Type arg0) { \
    typename Traits <Class>::Ptr_t ptr (new Class (arg0)); \
    ptr->init (ptr); \
    return ptr; \
  }

#define HPP_CONSTRAINTS_CB_CREATE2(Class, Arg0Type, Arg1Type) \
  static typename Traits <Class>::Ptr_t create (Arg0Type arg0, Arg1Type arg1) { \
    typename Traits <Class>::Ptr_t ptr (new Class (arg0, arg1)); \
    ptr->init (ptr); \
    return ptr; \
  }

#define HPP_CONSTRAINTS_CB_CREATE3(Class, Arg0Type, Arg1Type, Arg2Type) \
  static typename Traits <Class>::Ptr_t create (Arg0Type arg0, Arg1Type arg1, Arg2Type arg2) { \
    typename Traits <Class>::Ptr_t ptr (new Class (arg0, arg1, arg2)); \
    ptr->init (ptr); \
    return ptr; \
  }

#include <Eigen/SVD>

#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/center-of-mass-computation.hh>
#include <hpp/pinocchio/liegroup-element.hh>

#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/svd.hh>
#include <hpp/constraints/tools.hh>
#include <hpp/constraints/macros.hh>

namespace hpp {
  namespace constraints {
    /// \defgroup symbolic_calculus Symbolic calculus

    /// \addtogroup symbolic_calculus
    /// \{

    template <typename ValueType, typename JacobianType> class CalculusBaseAbstract;
    template <typename T> class Traits;

    template <typename LhsValue, typename RhsValue> class Expression;
    template <typename LhsValue, typename RhsValue> class CrossProduct;
    template <typename LhsValue, typename RhsValue> class ScalarProduct;
    template <typename LhsValue, typename RhsValue> class Difference;
    template <typename LhsValue, typename RhsValue> class Sum;
    template <typename RhsValue> class ScalarMultiply;
    template <typename RhsValue> class RotationMultiply;
    typedef eigen::matrix3_t CrossMatrix;
    typedef Eigen::Matrix <value_type, 1, Eigen::Dynamic, Eigen::RowMajor> RowJacobianMatrix;
    typedef Eigen::Matrix <value_type, 3, Eigen::Dynamic, Eigen::RowMajor> JacobianMatrix;

    template <typename Class>
    struct Traits {
      typedef HPP_CONSTRAINTS_CB_REF <Class> Ptr_t;
      typedef HPP_CONSTRAINTS_CB_WKREF <Class> WkPtr_t;
    };
    template <> struct Traits<value_type> {
      typedef value_type Ptr_t;
      typedef value_type WkPtr_t;
    };
    template <> struct Traits<pinocchio::Joint> {
      typedef JointPtr_t Ptr_t;
      typedef JointPtr_t WkPtr_t;
    };

    struct JointTranspose {
      JointTranspose (const JointPtr_t& joint) :
        j_ (joint) {}
      const JointPtr_t j_;
    };
    template <> struct Traits<JointTranspose> {
      typedef JointTranspose Ptr_t;
      typedef JointTranspose WkPtr_t;
    };

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
        typedef ValueType ValueType_t;
        typedef JacobianType JacobianType_t;

        virtual const ValueType& value () const = 0;
        virtual const JacobianType& jacobian () const = 0;
        virtual void computeValue (const ConfigurationIn_t arg) = 0;
        virtual void computeJacobian (const ConfigurationIn_t arg) = 0;
        virtual void invalidate () = 0;
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
             class ValueType = vector3_t,
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
        void computeValue (const ConfigurationIn_t arg) {
          if (vValid_) return;
          static_cast<T*>(this)->impl_value (arg);
          vValid_ = true;
        }
        void computeJacobian (const ConfigurationIn_t arg) {
          if (jValid_) return;
          static_cast<T*>(this)->impl_jacobian (arg);
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
        void computeCrossValue (const ConfigurationIn_t arg) {
          if (cValid_) return;
          computeValue (arg);
          computeCrossMatrix (value_, cross_);
          cValid_ = true;
        }

      protected:
        ValueType value_;
        JacobianType jacobian_;
        CrossType cross_;

        bool vValid_, jValid_, cValid_;

        void init (const typename Traits <T>::Ptr_t& ptr) {
          wkPtr_ = ptr;
        }

      private:
        typename Traits <T>::WkPtr_t  wkPtr_;

      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    HPP_CONSTRAINTS_CB_DEFINE_OPERATOR2(operator-, Difference)
    HPP_CONSTRAINTS_CB_DEFINE_OPERATOR2(operator+, Sum)
    HPP_CONSTRAINTS_CB_DEFINE_OPERATOR2(operator*, ScalarProduct)
    HPP_CONSTRAINTS_CB_DEFINE_OPERATOR2(operator^, CrossProduct)

    HPP_CONSTRAINTS_CB_DEFINE_OPERATOR1(operator*, value_type, ScalarMultiply)
    HPP_CONSTRAINTS_CB_DEFINE_OPERATOR1(operator*, JointPtr_t, RotationMultiply)
    HPP_CONSTRAINTS_CB_DEFINE_OPERATOR1(operator*, JointTranspose, RotationMultiply)

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

        static Ptr_t create (const typename Traits<LhsValue>::Ptr_t& lhs,
            const typename Traits<RhsValue>::Ptr_t& rhs) {
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

        Expression (const typename Traits<LhsValue>::Ptr_t& lhs,
            const typename Traits<RhsValue>::Ptr_t& rhs):
          rhs_ (rhs), lhs_ (lhs)
        {}

        inline void init (Ptr_t self) {
          self_ = self;
        }

        typename Traits<RhsValue>::Ptr_t rhs_;
        typename Traits<LhsValue>::Ptr_t lhs_;
        typename Traits<Expression>::WkPtr_t self_;

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

        HPP_CONSTRAINTS_CB_CREATE2 (CrossProduct, const typename Traits<LhsValue>::Ptr_t&, const typename Traits<RhsValue>::Ptr_t&)

        CrossProduct () {}

        CrossProduct (const CalculusBase <CrossProduct>& other) :
          Parent_t (other),
          e_ (static_cast <const CrossProduct&>(other).e_)
        {}

        CrossProduct (const typename Traits<LhsValue>::Ptr_t& lhs, const typename Traits<RhsValue>::Ptr_t& rhs):
          e_ (Expression < LhsValue, RhsValue >::create (lhs, rhs))
        {}

        void impl_value (const ConfigurationIn_t arg) {
          e_->lhs_->computeCrossValue (arg);
          e_->rhs_->computeValue (arg);
          this->value_ = e_->lhs_->cross () * e_->rhs_->value ();
        }
        void impl_jacobian (const ConfigurationIn_t arg) {
          e_->lhs_->computeCrossValue (arg);
          e_->rhs_->computeCrossValue (arg);
          e_->lhs_->computeJacobian (arg);
          e_->rhs_->computeJacobian (arg);
          this->jacobian_ = e_->lhs_->cross () * e_->rhs_->jacobian ()
                          - e_->rhs_->cross () * e_->lhs_->jacobian ();
        }
        void invalidate () {
          Parent_t::invalidate ();
          e_->rhs_->invalidate ();
          e_->lhs_->invalidate ();
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

        HPP_CONSTRAINTS_CB_CREATE2 (ScalarProduct, const typename Traits<LhsValue>::Ptr_t&, const typename Traits<RhsValue>::Ptr_t&)

        ScalarProduct () {}

        ScalarProduct (const CalculusBase <ScalarProduct>& other) :
          CalculusBase <ScalarProduct> (other),
          e_ (static_cast <const ScalarProduct&>(other).e_)
        {}

        ScalarProduct (const typename Traits<LhsValue>::Ptr_t& lhs, const typename Traits<RhsValue>::Ptr_t& rhs):
          e_ (Expression < LhsValue, RhsValue >::create (lhs, rhs))
        {}

        void impl_value (const ConfigurationIn_t arg) {
          e_->lhs_->computeValue (arg);
          e_->rhs_->computeValue (arg);
          this->value_ = e_->lhs_->value ().dot (e_->rhs_->value ());
        }
        void impl_jacobian (const ConfigurationIn_t arg) {
          e_->lhs_->computeValue (arg);
          e_->rhs_->computeValue (arg);
          e_->lhs_->computeJacobian (arg);
          e_->rhs_->computeJacobian (arg);
          this->jacobian_ = e_->lhs_->value ().transpose () * e_->rhs_->jacobian ()
                          + e_->rhs_->value ().transpose () * e_->lhs_->jacobian ();
        }
        void invalidate () {
          Parent_t::invalidate ();
          e_->rhs_->invalidate ();
          e_->lhs_->invalidate ();
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

        HPP_CONSTRAINTS_CB_CREATE2 (Difference, const typename Traits<LhsValue>::Ptr_t&, const typename Traits<RhsValue>::Ptr_t&)

        Difference () {}

        Difference (const CalculusBase <Difference>& other) :
          CalculusBase <Difference> (other),
          e_ (static_cast <const Difference&>(other).e_)
        {}

        Difference (const typename Traits<LhsValue>::Ptr_t& lhs, const typename Traits<RhsValue>::Ptr_t& rhs):
          e_ (Expression < LhsValue, RhsValue >::create (lhs, rhs))
        {}

        void impl_value (const ConfigurationIn_t arg) {
          e_->lhs_->computeValue (arg);
          e_->rhs_->computeValue (arg);
          this->value_ = e_->lhs_->value () - e_->rhs_->value ();
        }
        void impl_jacobian (const ConfigurationIn_t arg) {
          e_->lhs_->computeJacobian (arg);
          e_->rhs_->computeJacobian (arg);
          this->jacobian_ = e_->lhs_->jacobian () - e_->rhs_->jacobian ();
        }
        void invalidate () {
          Parent_t::invalidate ();
          e_->rhs_->invalidate ();
          e_->lhs_->invalidate ();
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

        HPP_CONSTRAINTS_CB_CREATE2 (Sum, const typename Traits<LhsValue>::Ptr_t&, const typename Traits<RhsValue>::Ptr_t&)

        Sum () {}

        Sum (const CalculusBase < Sum >& other) :
          CalculusBase < Sum > (other),
          e_ (static_cast <const Sum&>(other).e_)
        {}

        Sum (const typename Traits<RhsValue>::Ptr_t& rhs, const typename Traits<LhsValue>::Ptr_t& lhs):
          e_ (Expression < LhsValue, RhsValue >::create (lhs, rhs))
        {}

        void impl_value (const ConfigurationIn_t arg) {
          e_->lhs_->computeValue (arg);
          e_->rhs_->computeValue (arg);
          this->value_ = e_->lhs_->value () + e_->rhs_->value ();
        }
        void impl_jacobian (const ConfigurationIn_t arg) {
          e_->lhs_->computeJacobian (arg);
          e_->rhs_->computeJacobian (arg);
          this->jacobian_ = e_->lhs_->jacobian () + e_->rhs_->jacobian ();
        }
        void invalidate () {
          Parent_t::invalidate ();
          e_->rhs_->invalidate ();
          e_->lhs_->invalidate ();
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

        HPP_CONSTRAINTS_CB_CREATE2 (ScalarMultiply, const typename Traits<value_type>::Ptr_t&, const typename Traits<RhsValue>::Ptr_t&)

        ScalarMultiply () {}

        ScalarMultiply (const CalculusBase < ScalarMultiply >& other) :
          CalculusBase < ScalarMultiply > (other),
          e_ (static_cast <const ScalarMultiply&>(other).e_)
        {}

        ScalarMultiply (const typename Traits<value_type>::Ptr_t& scalar, const typename Traits<RhsValue>::Ptr_t& rhs):
          e_ (Expression < value_type, RhsValue >::create (scalar, rhs))
        {}

        void impl_value (const ConfigurationIn_t arg) {
          e_->rhs_->computeValue (arg);
          this->value_ = e_->lhs_ * e_->rhs_->value ();
        }
        void impl_jacobian (const ConfigurationIn_t arg) {
          e_->rhs_->computeJacobian (arg);
          this->jacobian_ = e_->lhs_ * e_->rhs_->jacobian ();
        }
        void invalidate () {
          Parent_t::invalidate ();
          e_->rhs_->invalidate ();
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

        HPP_CONSTRAINTS_CB_CREATE2 (RotationMultiply, const typename Traits<pinocchio::Joint>::Ptr_t&, const typename Traits<RhsValue>::Ptr_t&)

        HPP_CONSTRAINTS_CB_CREATE2 (RotationMultiply, const typename Traits<JointTranspose>::Ptr_t&, const typename Traits<RhsValue>::Ptr_t&)

        RotationMultiply (const typename Traits<JointTranspose>::Ptr_t& joint, const typename Traits<RhsValue>::Ptr_t& rhs):
          e_ (Expression < pinocchio::Joint, RhsValue >::create (joint.j_, rhs)),
          transpose_ (true)
        {}

        RotationMultiply (const typename Traits<pinocchio::Joint>::Ptr_t& joint, const typename Traits<RhsValue>::Ptr_t& rhs,
            bool transpose = false):
          e_ (Expression < pinocchio::Joint, RhsValue >::create (joint, rhs)),
          transpose_ (transpose)
        {}

        void impl_value (const ConfigurationIn_t arg) {
          e_->rhs_->computeValue (arg);
          const matrix3_t& R = e_->lhs_->currentTransformation ().rotation ();
          if (transpose_)
            this->value_ = R.transpose() * e_->rhs_->value ();
          else
            this->value_ = R             * e_->rhs_->value ();
        }
        void impl_jacobian (const ConfigurationIn_t arg) {
          e_->rhs_->computeJacobian (arg);
          e_->rhs_->computeCrossValue (arg);
          const JointJacobian_t& J = e_->lhs_->jacobian ();
          const matrix3_t& R = e_->lhs_->currentTransformation ().rotation ();
          if (transpose_)
            this->jacobian_ = R.transpose()
              * ((e_->rhs_->cross () * R) * J.bottomRows<3>() + e_->rhs_->jacobian ());
          else
            this->jacobian_ = R
              * ((e_->rhs_->cross () * R) * J.bottomRows<3>() + e_->rhs_->jacobian ());
        }
        void invalidate () {
          Parent_t::invalidate ();
          e_->rhs_->invalidate ();
        }

      protected:
        typename Expression < pinocchio::Joint, RhsValue >::Ptr_t e_;

        friend class Expression <pinocchio::Joint, RhsValue>;

      private:
        bool transpose_;
    };

    /// Basic expression representing a point in a joint frame.
    class PointInJoint : public CalculusBase <PointInJoint>
    {
      public:
        typedef CalculusBase <PointInJoint> Parent_t;

        HPP_CONSTRAINTS_CB_CREATE2 (PointInJoint, const JointPtr_t&, const vector3_t&)
        HPP_CONSTRAINTS_CB_CREATE3 (PointInJoint, const JointPtr_t&, const vector3_t&, const size_type&)

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
        {
          assert (joint_ != NULL);
        }

        PointInJoint (const JointPtr_t& joint,
            const vector3_t& pointInLocalFrame,
            size_type nbDof) :
          joint_ (joint), local_ (pointInLocalFrame),
          center_ (pointInLocalFrame.isZero ())
        {
          if (joint_ == NULL) {
            for (int i = 0; i < 3; ++i) this->value_[i] = local_[i];
            this->jacobian_.resize (3, nbDof);
            this->jacobian_.setZero ();
            this->cross_.setZero ();
          }
        }


        const JointPtr_t& joint () const {
          return joint_;
        }
        const vector3_t& local () const {
          return local_;
        }
        void impl_value (const ConfigurationIn_t ) {
          if (joint_ == NULL) return;
          this->value_ = joint_->currentTransformation ().act (local_);
        }
        void impl_jacobian (const ConfigurationIn_t ) {
          if (joint_ == NULL) return;
          const JointJacobian_t& J (joint_->jacobian ());
          const matrix3_t& R = joint_->currentTransformation ().rotation ();
          this->jacobian_.noalias() = R * J.topRows<3>();
          if (!center_) {
            computeCrossRXl ();
            this->jacobian_.noalias() -= (this->cross_ * R) * J.bottomRows<3>();
          }
        }
        void computeCrossRXl () {
          if (joint_ == NULL) return;
          if (center_) {
            this->cross_.setZero ();
            return;
          }
          computeCrossMatrix (
              joint_->currentTransformation ().rotation () * local_,
              this->cross_);
        }

      protected:
        JointPtr_t joint_;
        vector3_t local_;
        bool center_;

      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    /// Basic expression representing a vector in a joint frame.
    class VectorInJoint : public CalculusBase <VectorInJoint>
    {
      public:

        HPP_CONSTRAINTS_CB_CREATE2 (VectorInJoint, const JointPtr_t&, const vector3_t&)
        HPP_CONSTRAINTS_CB_CREATE3 (VectorInJoint, const JointPtr_t&, const vector3_t&, const size_type&)

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

        VectorInJoint (const JointPtr_t& joint,
            const vector3_t& vectorInLocalFrame,
            const size_type& nbDof) :
          joint_ (joint), vector_ (vectorInLocalFrame)
        {
          if (joint_ == NULL) {
            for (int i = 0; i < 3; ++i) this->value_[i] = vector_[i];
            this->jacobian_.resize (3, nbDof);
            this->jacobian_.setZero ();
            this->cross_.setZero ();
          }
        }

        const JointPtr_t& joint () const {
          return joint_;
        }
        const vector3_t& vector () const {
          return vector_;
        }
        void impl_value (const ConfigurationIn_t ) {
          if (joint_ == NULL) return;
          this->value_ = joint_->currentTransformation ().rotation () * vector_;
        }
        void impl_jacobian (const ConfigurationIn_t ) {
          if (joint_ == NULL) return;
          const JointJacobian_t& J (joint_->jacobian ());
          const matrix3_t& R = joint_->currentTransformation ().rotation ();
          computeCrossRXl ();
          this->jacobian_.noalias() = (- this->cross_ * R ) * J.bottomRows<3>();
        }
        void computeCrossRXl () {
          if (joint_ == NULL) return;
          computeCrossMatrix (
              joint_->currentTransformation ().rotation () * vector_,
              this->cross_);
        }

      protected:
        JointPtr_t joint_;
        vector3_t vector_;

      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    /// Basic expression mapping a function as an expression
    template <typename FunctionType>
    class FunctionExp : public CalculusBase <FunctionExp<FunctionType>, vector_t, matrix_t>
    {
      public:
        typedef CalculusBase <FunctionExp<FunctionType>, vector_t, matrix_t> Parent_t;
        typedef boost::shared_ptr<FunctionType> FunctionTypePtr_t;

        HPP_CONSTRAINTS_CB_CREATE1 (FunctionExp<FunctionType>, const FunctionTypePtr_t&)

        FunctionExp () {}

        FunctionExp (const Parent_t& other) :
          Parent_t (other), f_ (other.f_), lge_ (other.lge_)
        {}

        FunctionExp (const FunctionExp& other) :
          Parent_t (other), f_ (other.f_), lge_ (other.lge_)
        {}

        /// Constructor
        ///
        /// \param func the inner function
        FunctionExp (const FunctionTypePtr_t& func) :
          Parent_t (vector_t::Zero(func->outputSize()),
                    matrix_t::Zero(func->outputDerivativeSize(), func->inputDerivativeSize())),
          f_ (func), lge_ (func->outputSpace())
        {}

        void impl_value (const ConfigurationIn_t arg)
        {
          f_->value (lge_, arg);
          this->value_ = lge_.vector();
        }
        void impl_jacobian (const ConfigurationIn_t arg) { f_->jacobian (this->jacobian_, arg); }

      private:
        FunctionTypePtr_t f_;
        LiegroupElement lge_;
    };

    /// Basic expression representing a static point
    ///
    /// Its value is constant and its jacobian is a zero matrix.
    class Point : public CalculusBase <Point, vector3_t, JacobianMatrix>
    {
      public:
        HPP_CONSTRAINTS_CB_CREATE2 (Point, const vector3_t&, const size_t&)

        Point () {}

        Point (const CalculusBase<Point, vector3_t, JacobianMatrix>& other) :
          CalculusBase <Point, eigen::vector3_t, JacobianMatrix> (other)
        {
        }

        Point (const Point& point) :
          CalculusBase <Point, vector3_t, JacobianMatrix> (point)
        {
        }

        /// Constructor
        ///
        /// \param point the static point
        /// \param jacobianNbCols number of column of the jacobian
        Point (const vector3_t& point, size_t jacobianNbCols) :
          CalculusBase <Point, vector3_t, JacobianMatrix>
          (point, JacobianMatrix::Zero (3, jacobianNbCols))
        {
        }

        void impl_value (const ConfigurationIn_t ) {}
        void impl_jacobian (const ConfigurationIn_t ) {}
    };

    /// Basic expression representing a COM.
    class PointCom : public CalculusBase <PointCom, vector3_t, ComJacobian_t>
    {
      public:
        typedef CalculusBase <PointCom, vector3_t, ComJacobian_t > Parent_t;
        HPP_CONSTRAINTS_CB_CREATE1 (PointCom, const CenterOfMassComputationPtr_t&)

        PointCom () {}

        PointCom (const Parent_t& other):
          Parent_t (other),
          comc_ (static_cast <const PointCom&>(other).centerOfMassComputation ())
        {
        }

        PointCom (const CenterOfMassComputationPtr_t& comc): comc_ (comc)
        {}

        inline const vector3_t& value () const {
          return comc_->com();
        }
        inline const ComJacobian_t& jacobian () const {
          return comc_->jacobian();
        }

        const CenterOfMassComputationPtr_t& centerOfMassComputation () const {
          return comc_;
        }
        void impl_value (const ConfigurationIn_t ) {
          comc_->compute (hpp::pinocchio::COM);
        }
        void impl_jacobian (const ConfigurationIn_t ) {
          comc_->compute (hpp::pinocchio::COMPUTE_ALL);
          // TODO: there is memory and time to be saved here as this copy is
          // not important.
          //this->jacobian_ = comc_->jacobian ();
        }

      protected:
        CenterOfMassComputationPtr_t comc_;
    };

    class JointFrame : public CalculusBase <JointFrame, Eigen::Matrix<value_type, 6, 1>, Eigen::Matrix<value_type, 6, Eigen::Dynamic> >
    {
      public:
        typedef CalculusBase <JointFrame, ValueType_t, JacobianType_t > Parent_t;

        HPP_CONSTRAINTS_CB_CREATE1 (JointFrame, const JointPtr_t&)

        JointFrame () {}

        JointFrame (const Parent_t& other) :
          Parent_t (other),
          joint_ (static_cast <const JointFrame&>(other).joint_)
        {}

        JointFrame (const JointFrame& jf) :
          Parent_t (jf), joint_ (jf.joint ())
        {}

        JointFrame (const JointPtr_t& joint) :
          joint_ (joint)
        {
          assert (joint_ != NULL);
          this->jacobian_.resize(6,joint->robot()->numberDof());
        }

        const JointPtr_t& joint () const {
          return joint_;
        }
        void impl_value (const ConfigurationIn_t ) {
          const Transform3f& M = joint_->currentTransformation ();
          this->value_.head<3>() = M.translation ();
          logSO3 (M.rotation(), theta_, this->value_.tail<3>());
        }
        void impl_jacobian (const ConfigurationIn_t arg) {
          computeValue (arg);
          const JointJacobian_t& J (joint_->jacobian ());
          const matrix3_t& R (joint_->currentTransformation ().rotation ());
          // Compute vector r
          eigen::matrix3_t Jlog;
          assert (theta_ >= 0);
          JlogSO3 (theta_, this->value_.tail<3>(), Jlog);
          this->jacobian_.topRows<3>().noalias() = R * J.topRows<3>();
          this->jacobian_.bottomRows<3>().noalias() = Jlog * J.bottomRows<3>();
        }

      protected:
        JointPtr_t joint_;
        double theta_;

      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
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
        typedef typename Traits<Element_t>::Ptr_t ElementPtr_t;
        typedef Eigen::JacobiSVD <Value_t> SVD_t;

        HPP_CONSTRAINTS_CB_CREATE2 (MatrixOfExpressions, const Eigen::Ref<const Value_t>&, const Eigen::Ref<const Jacobian_t>&)

        MatrixOfExpressions (const Eigen::Ref<const Value_t>& value,
            const Eigen::Ref<const Jacobian_t>& jacobian) :
          Parent_t (value, jacobian),
          nRows_ (0), nCols_ (0),
          svd_ (value.rows(), value.cols(), Eigen::ComputeFullU | Eigen::ComputeFullV),
          piValid_ (false), svdValid_ (false)
        {}

        MatrixOfExpressions (const Parent_t& other) :
          Parent_t (other),
          nRows_ (static_cast <const MatrixOfExpressions&>(other).nRows_),
          nCols_ (static_cast <const MatrixOfExpressions&>(other).nCols_),
          elements_ (static_cast <const MatrixOfExpressions&>(other).elements_),
          svd_ (static_cast <const MatrixOfExpressions&>(other).svd_),
          piValid_ (static_cast <const MatrixOfExpressions&>(other).piValid_),
          svdValid_ (static_cast <const MatrixOfExpressions&>(other).svdValid_)
        {
        }

        MatrixOfExpressions (const MatrixOfExpressions& matrix) :
          Parent_t (matrix),
          nRows_ (matrix.nRows_), nCols_ (matrix.nCols_),
          elements_ (matrix.elements_),
          svd_ (matrix.svd_),
          piValid_ (matrix.piValid_),
          svdValid_ (matrix.svdValid_)
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

        void impl_value (const ConfigurationIn_t arg) {
          size_type r = 0, c = 0, nr = 0, nc = 0;
          for (std::size_t i = 0; i < nRows_; ++i) {
            c = 0;
            nr = elements_[i][0]->value().rows();
            for (std::size_t j = 0; j < nCols_; ++j) {
              elements_[i][j]->computeValue (arg);
              assert (nr == elements_[i][j]->value().rows());
              nc = elements_[i][j]->value().cols();
              this->value_.block (r, c, nr, nc)
                = elements_[i][j]->value();
              c += nc;
            }
            r += nr;
          }
        }
        void impl_jacobian (const ConfigurationIn_t arg) {
          size_type r = 0, c = 0, nr = 0, nc = 0;
          for (std::size_t i = 0; i < nRows_; ++i) {
            c = 0;
            nr = elements_[i][0]->jacobian().rows();
            for (std::size_t j = 0; j < nCols_; ++j) {
              elements_[i][j]->computeJacobian (arg);
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
        void computeSVD (const ConfigurationIn_t arg) {
          if (svdValid_) return;
          this->computeValue (arg);
          svd_.compute (this->value_);
          HPP_DEBUG_SVDCHECK(svd_);
          svdValid_ = true;
        }
        void computePseudoInverse (const ConfigurationIn_t arg) {
          if (piValid_) return;
          this->computeValue (arg);
          this->computeSVD(arg);
          pi_.resize (this->value_.cols(), this->value_.rows());
          pseudoInverse <SVD_t> (svd_, pi_);
          piValid_ = true;
        }
        void computePseudoInverseJacobian (const ConfigurationIn_t arg, const Eigen::Ref <const Eigen::Matrix<value_type, Eigen::Dynamic, 1> >& rhs) {
          this->computeJacobian (arg);
          computePseudoInverse (arg);
          const std::size_t nbDof = elements_[0][0]->jacobian().cols();
          const std::size_t inSize = this->value_.cols();
          const vector_t piTrhs = svd_.solve (rhs);
          assert (pi_.rows () == (size_type)inSize);

          Jacobian_t cache (this->jacobian_.rows(), nbDof);
          jacobianTimes (arg, piTrhs, cache);
          pij_.noalias() = - pi_ * cache;
          cache.resize (inSize, nbDof);

          pkInv_.resize (pi_.cols(), pi_.cols());
          projectorOnKernelOfInv <SVD_t> (svd_, pkInv_, true);
          jacobianTransposeTimes (arg, pkInv_ * rhs, cache);
          pij_.noalias() += (pi_ * pi_.transpose()) * cache;

          jacobianTransposeTimes (arg, pi_.transpose() * piTrhs , cache);
          pk_.resize (inSize, inSize);
          projectorOnKernel <SVD_t> (svd_, pk_, true);
          pij_.noalias() += pk_ * cache;
        }

        void jacobianTimes (const ConfigurationIn_t arg, const Eigen::Ref <const Eigen::Matrix<value_type, Eigen::Dynamic, 1> >& rhs, Eigen::Ref<Jacobian_t> cache) const {
          size_type r = 0, c = 0, nr = 0, nc = 0;
          cache.setZero();
          for (std::size_t i = 0; i < nRows_; ++i) {
            c = 0;
            nr = elements_[i][0]->jacobian().rows();
            for (std::size_t j = 0; j < nCols_; ++j) {
              elements_[i][j]->computeJacobian (arg);
              assert (nr == elements_[i][j]->jacobian().rows());
              nc = elements_[i][j]->jacobian().cols();
              cache.middleRows (r,nr).noalias() +=
                this->jacobian_.block (r, c, nr, nc) * rhs[j];
              c += nc;
            }
            r += nr;
          }
        }

        void jacobianTransposeTimes (const ConfigurationIn_t arg, const Eigen::Ref <const Eigen::Matrix<value_type, Eigen::Dynamic, 1> >& rhs, Eigen::Ref<Jacobian_t> cache) const {
          size_type r = 0, c = 0, nr = 0, nc = 0;
          cache.setZero();
          for (std::size_t i = 0; i < nRows_; ++i) {
            c = 0;
            nr = elements_[i][0]->jacobian().rows();
            for (std::size_t j = 0; j < nCols_; ++j) {
              elements_[i][j]->computeJacobian (arg);
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
          piValid_ = false;
          svdValid_ = false;
        }

        std::size_t nRows_, nCols_;
        std::vector <std::vector <ElementPtr_t> > elements_;

      private:
        SVD_t svd_;
        matrix_t pkInv_, pk_;
        PseudoInv_t pi_;
        PseudoInvJacobian_t pij_;
        bool piValid_, svdValid_;

      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_SYMBOLIC_CALCULUS_HH
