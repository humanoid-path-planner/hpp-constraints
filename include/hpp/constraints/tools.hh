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
      m (0,1) = -v [2]; m (1,0) = v [2];
      m (0,2) = v [1]; m (2,0) = -v [1];
      m (1,2) = -v [0]; m (2,1) = v [0];
    }

    /// \defgroup symbolic_calculus Symbolic calculus

    /// \addtogroup symbolic_calculus
    /// \{

    class CalculusBaseAbstract;
    typedef boost::shared_ptr <CalculusBaseAbstract> CalculusPtr_t;

    template <typename LhsValue, typename RhsValue> class Expression;
    template <typename LhsValue, typename RhsValue> class CrossProduct;
    template <typename LhsValue, typename RhsValue> class Difference;
    template <typename LhsValue, typename RhsValue> class Sum;
    template <typename RhsValue> class ScalarMultiply;
    template <typename RhsValue> class RotationMultiply;
    typedef eigen::matrix3_t CrossMatrix;
    typedef Eigen::Matrix <value_type, 3, Eigen::Dynamic> JacobianMatrix;

    /// Abstract class defining a basic common interface.
    ///
    /// The purpose of this class is to allow the user to define an expression
    /// without requiring to explicitly write its type. The type will be
    /// automatically deduced by the compiler.
    ///
    /// \code
    ///   // First define a, b and c with basic elements.
    ///   CalculusPtr_t myExpression_ptr = CalculusBaseAbstract::create (a + b * c)
    /// \endcode
    class CalculusBaseAbstract
    {
      public:
        virtual const eigen::vector3_t& value () const = 0;
        virtual const JacobianMatrix& jacobian () const = 0;
        virtual void computeValue () = 0;
        virtual void computeJacobian () = 0;

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
    template <class T>
    class CalculusBase : public CalculusBaseAbstract
    {
      public:
        CalculusBase () : cross_ (CrossMatrix::Zero()) {}

        inline const eigen::vector3_t& value () const {
          return value_;
        }
        inline const JacobianMatrix& jacobian () const {
          return jacobian_;
        }
        inline const CrossMatrix& cross () const {
          return cross_;
        }
        void computeCrossValue () {
          T& derived = static_cast <T&> (*this);
          derived.computeValue ();
          computeCrossMatrix (derived.value (), cross_);
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

        ScalarMultiply < T > operator* (const value_type& scalar) const {
          return ScalarMultiply < T > (scalar, static_cast<const T>(*this));
        }

      protected:
        eigen::vector3_t value_;
        JacobianMatrix jacobian_;
        CrossMatrix cross_;
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
    };

    /// Cross product of two expressions.
    template <typename LhsValue, typename RhsValue>
    class CrossProduct :
      public CalculusBase < CrossProduct < LhsValue, RhsValue > >
    {
      public:
        CrossProduct () {}

        CrossProduct (const CalculusBase <CrossProduct>& other) :
          e_ (static_cast <const CrossProduct&>(other).e_)
        {}

        CrossProduct (const LhsValue& lhs, const RhsValue& rhs):
          e_ (Expression < LhsValue, RhsValue >::create (lhs, rhs))
        {}

        void computeValue () {
          e_->lhs_.computeCrossValue ();
          e_->rhs_.computeValue ();
          this->value_ = e_->lhs_.cross () * e_->rhs_.value ();
        }
        void computeJacobian () {
          e_->lhs_.computeCrossValue ();
          e_->rhs_.computeCrossValue ();
          e_->lhs_.computeJacobian ();
          e_->rhs_.computeJacobian ();
          this->jacobian_ = e_->lhs_.cross () * e_->rhs_.jacobian ()
                          - e_->rhs_.cross () * e_->lhs_.jacobian ();
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
        Difference () {}

        Difference (const CalculusBase <Difference>& other) :
          e_ (static_cast <const Difference&>(other).e_)
        {}

        Difference (const LhsValue& lhs, const RhsValue& rhs):
          e_ (Expression < LhsValue, RhsValue >::create (lhs, rhs))
        {}

        void computeValue () {
          e_->lhs_.computeValue ();
          e_->rhs_.computeValue ();
          this->value_ = e_->lhs_.value () - e_->rhs_.value ();
        }
        void computeJacobian () {
          e_->lhs_.computeJacobian ();
          e_->rhs_.computeJacobian ();
          this->jacobian_ = e_->lhs_.jacobian () - e_->rhs_.jacobian ();
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
        Sum () {}

        Sum (const CalculusBase < Sum >& other) :
          e_ (static_cast <const Sum&>(other).e_)
        {}

        Sum (const RhsValue& rhs, const LhsValue& lhs):
          e_ (Expression < LhsValue, RhsValue >::create (lhs, rhs))
        {}

        void computeValue () {
          e_->lhs_.computeValue ();
          e_->rhs_.computeValue ();
          this->value_ = e_->lhs_.value () + e_->rhs_.value ();
        }
        void computeJacobian () {
          e_->lhs_.computeJacobian ();
          e_->rhs_.computeJacobian ();
          this->jacobian_ = e_->lhs_.jacobian () + e_->rhs_.jacobian ();
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
        ScalarMultiply () {}

        ScalarMultiply (const CalculusBase < ScalarMultiply >& other) :
          e_ (static_cast <const ScalarMultiply&>(other).e_)
        {}

        ScalarMultiply (const value_type& scalar, const RhsValue& rhs):
          e_ (Expression < value_type, RhsValue >::create (scalar, rhs))
        {}

        void computeValue () {
          e_->rhs_.computeValue ();
          this->value_ = e_->lhs_ * e_->rhs_.value ();
        }
        void computeJacobian () {
          e_->rhs_.computeJacobian ();
          this->jacobian_ = e_->lhs_ * e_->rhs_.jacobian ();
        }

      protected:
        typename Expression < value_type, RhsValue >::Ptr_t e_;

        friend class Expression <value_type, RhsValue>;
    };

    /// Multiplication of an expression by a rotation matrix.
    template <typename RhsValue>
    class RotationMultiply :
      public CalculusBase < ScalarMultiply < RhsValue > >
    {
      public:
        RotationMultiply () {}

        RotationMultiply (const CalculusBase < RotationMultiply >& other) :
          e_ (static_cast <const RotationMultiply&>(other).e_),
          transpose_ (other.transpose_)
        {}

        RotationMultiply (const JointPtr_t& joint, const RhsValue& rhs,
            bool transpose = false):
          e_ (Expression < JointPtr_t, RhsValue >::create (joint, rhs)),
          transpose_ (transpose)
        {}

        void computeValue () {
          e_->rhs_.computeValue ();
          computeRotationMatrix ();
          this->value_ = R * e_->rhs_.value ();
        }
        void computeJacobian () {
          e_->rhs_.computeJacobian ();
          computeRotationMatrix ();
          e_->rhs_.computeCrossValue ();
          const JointJacobian_t& J = e_->lhs_->jacobian ();
          this->jacobian_ = R * (e_->rhs_.cross () * J.bottomRows (3) + e_->rhs_.jacobian ());
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

        PointInJoint (const CalculusBase<PointInJoint>& other) {
          const PointInJoint& o = static_cast <const PointInJoint&>(other);
          joint_ = o.joint ();
          local_ = o.local ();
          center_= local_.isZero ();
        }

        PointInJoint (const PointInJoint& pointInJoint) :
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
        void computeValue () {
          g_ = joint_->currentTransformation ().transform (local_);
          for (int i = 0; i < 3; ++i) this->value_[i] = g_[i];
        }
        void computeJacobian () {
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
    };

    /// Basic expression representing a static point
    ///
    /// Its value is constant and its jacobian is a zero matrix.
    class Point : public CalculusBase <Point>
    {
      public:
        Point () {}

        Point (const CalculusBase<Point>& other) {
          const Point& o = static_cast <const Point&>(other);
          this->value_ = o.value ();
          this->jacobian_ = o.jacobian ();
        }

        Point (const Point& point) {
          this->value_ = point.value ();
          this->jacobian_ = point.jacobian ();
        }

        /// Constructor
        ///
        /// \param point the static point
        /// \param jacobianNbCols number of column of the jacobian
        Point (const vector3_t& point, size_t jacobianNbCols) {
          for (int i = 0; i < 3; ++i) this->value_[i] = point[i];
          this->jacobian_ = JacobianMatrix::Zero (3, jacobianNbCols);
        }

        void computeValue () {}
        void computeJacobian () {}
    };

    /// Basic expression representing a COM.
    class PointCom : public CalculusBase <PointCom>
    {
      public:
        PointCom () {}

        PointCom (const CalculusBase<PointCom>& other) {
          comc_ = static_cast <const PointCom&>(other).centerOfMassComputation ();
        }

        PointCom (const CenterOfMassComputationPtr_t& comc): comc_ (comc)
        {}

        const CenterOfMassComputationPtr_t& centerOfMassComputation () const {
          return comc_;
        }
        void computeValue () {
          comc_->compute (Device::COM);
          for (int i = 0; i < 3; ++i) this->value_[i] = comc_->com ()[i];
        }
        void computeJacobian () {
          comc_->compute (Device::JACOBIAN);
          // TODO: there is memory and time to be saved here as this copy is
          // not important.
          this->jacobian_ = comc_->jacobian ();
        }

      protected:
        CenterOfMassComputationPtr_t comc_;
    };

    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_TOOL_HH
