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

    template <typename LhsValue, typename RhsValue> class CrossProduct;
    template <typename LhsValue, typename RhsValue> class Difference;
    template <typename LhsValue, typename RhsValue> class Sum;
    template <typename RhsValue> class ScalarMultiply;
    template <typename RhsValue> class RotationMultiply;
    typedef eigen::matrix3_t CrossMatrix;
    typedef Eigen::Matrix <value_type, 3, Eigen::Dynamic> JacobianMatrix;


    template <class T>
    class CalculusBase
    {
      public:
        CalculusBase () : cross_ (CrossMatrix::Zero()) {}

        const eigen::vector3_t& value () const {
          return value_;
        }
        const JacobianMatrix& jacobian () const {
          return jacobian_;
        }
        const CrossMatrix& cross () const {
          return cross_;
        }
        // void computeValue () {
          // return static_cast <T*> (this)->computeValue ();
        // }
        // void computeJacobian () {
          // return static_cast <T*> (this)->computeJacobian ();
        // }
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

    template <typename LhsValue, typename RhsValue>
    class Expression
    {
      public:
        Expression () {}

        Expression (const Expression& other):
          rhs_ (other.rhs()), lhs_ (other.lhs())
        {}

        Expression (const LhsValue& lhs, const RhsValue& rhs):
          rhs_ (rhs), lhs_ (lhs)
        {}

        const LhsValue& lhs () const {
          return lhs_;
        }

        const RhsValue& rhs () const {
          return rhs_;
        }

      protected:
        RhsValue rhs_;
        LhsValue lhs_;
    };

    template <typename LhsValue, typename RhsValue>
    class CrossProduct :
      public CalculusBase < CrossProduct < LhsValue, RhsValue > >,
      public Expression < LhsValue, RhsValue >
    {
      public:
        CrossProduct () {}

        CrossProduct (const CalculusBase <CrossProduct>& other) {
          const CrossProduct& o = static_cast <const CrossProduct&> (other);
          this->rhs_ = o.rhs ();
          this->lhs_ = o.lhs ();
        }

        CrossProduct (const LhsValue& lhs, const RhsValue& rhs):
          Expression < LhsValue, RhsValue > (lhs, rhs)
        {}

        void computeValue () {
          this->lhs_.computeCrossValue ();
          this->rhs_.computeValue ();
          this->value_ = this->lhs_.cross () * this->rhs_.value ();
        }
        void computeJacobian () {
          this->lhs_.computeCrossValue ();
          this->rhs_.computeCrossValue ();
          this->lhs_.computeJacobian ();
          this->rhs_.computeJacobian ();
          this->jacobian_ = this->lhs_.cross () * this->rhs_.jacobian ()
                          - this->rhs_.cross () * this->lhs_.jacobian ();
        }
    };

    template <typename LhsValue, typename RhsValue>
    class Difference :
      public CalculusBase < Difference < LhsValue, RhsValue > >,
      public Expression < LhsValue, RhsValue >
    {
      public:
        Difference () {}

        Difference (const CalculusBase <Difference>& other) {
          const Difference& o = static_cast <const Difference&> (other);
          this->rhs_ = o.rhs ();
          this->lhs_ = o.lhs ();
        }

        Difference (const LhsValue& lhs, const RhsValue& rhs):
          Expression < LhsValue, RhsValue > (lhs, rhs)
        {}

        void computeValue () {
          this->lhs_.computeValue ();
          this->rhs_.computeValue ();
          this->value_ = this->lhs_.value () - this->rhs_.value ();
        }
        void computeJacobian () {
          this->lhs_.computeJacobian ();
          this->rhs_.computeJacobian ();
          this->jacobian_ = this->lhs_.jacobian () - this->rhs_.jacobian ();
        }
    };

    template <typename LhsValue, typename RhsValue>
    class Sum :
      public CalculusBase < Sum < LhsValue, RhsValue > >,
      public Expression < LhsValue, RhsValue >
    {
      public:
        Sum () {}

        Sum (const CalculusBase < Sum >& other) {
          const Sum& o = static_cast <const Sum&> (other);
          this->rhs_ = o.rhs ();
          this->lhs_ = o.lhs ();
        }

        Sum (const RhsValue& rhs, const LhsValue& lhs):
          Expression < LhsValue, RhsValue > (lhs, rhs)
        {}

        void computeValue () {
          this->lhs_.computeValue ();
          this->rhs_.computeValue ();
          this->value_ = this->lhs_.value () + this->rhs_.value ();
        }
        void computeJacobian () {
          this->lhs_.computeJacobian ();
          this->rhs_.computeJacobian ();
          this->jacobian_ = this->lhs_.jacobian () + this->rhs_.jacobian ();
        }
    };

    template <typename RhsValue>
    class ScalarMultiply :
      public CalculusBase < ScalarMultiply < RhsValue > >,
      public Expression < value_type, RhsValue >
    {
      public:
        ScalarMultiply () {}

        ScalarMultiply (const CalculusBase < ScalarMultiply >& other) {
          const ScalarMultiply& o = static_cast <const ScalarMultiply&> (other);
          this->rhs_ = o.rhs ();
          this->lhs_ = o.lhs ();
        }

        ScalarMultiply (const value_type& scalar, const RhsValue& rhs):
          Expression < value_type, RhsValue > (scalar, rhs)
        {}

        void computeValue () {
          this->rhs_.computeValue ();
          this->value_ = this->lhs_ * this->rhs_.value ();
        }
        void computeJacobian () {
          this->rhs_.computeJacobian ();
          this->jacobian_ = this->lhs_ * this->rhs_.jacobian ();
        }
    };

    template <typename RhsValue>
    class RotationMultiply :
      public CalculusBase < ScalarMultiply < RhsValue > >,
      public Expression < JointPtr_t, RhsValue >
    {
      public:
        RotationMultiply () {}

        RotationMultiply (const CalculusBase < RotationMultiply >& other) {
          const RotationMultiply& o = static_cast <const RotationMultiply&> (other);
          this->rhs_ = o.rhs ();
          this->lhs_ = o.lhs ();
        }

        RotationMultiply (const JointPtr_t& joint, const RhsValue& rhs,
            bool transpose = false):
          Expression < JointPtr_t, RhsValue > (joint, rhs),
          transpose_ (transpose)
        {}

        void computeValue () {
          this->rhs_.computeValue ();
          computeRotationMatrix ();
          this->value_ = R * this->rhs_.value ();
        }
        void computeJacobian () {
          this->rhs_.computeJacobian ();
          computeRotationMatrix ();
          const JointJacobian_t& J = this->lhs_->jacobian ();
          this->jacobian_ = R * (this->rhs_.cross () * J.bottomRows (3) * this->rhs_.jacobian ());
        }

      private:
        void computeRotationMatrix () {
          const fcl::Matrix3f& Rfcl =
            this->lhs_->currentTransformation ().getRotation ();
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
              joint_->currentTransformation ().getRotation () * (- local_),
              this->cross_);
        }

      protected:
        JointPtr_t joint_;
        vector3_t local_;
        bool center_;

        vector3_t g_;
    };

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
          // TODO: there is memory and time to be save here as this copy is
          // not important.
          this->jacobian_ = comc_->jacobian ();
        }

      protected:
        CenterOfMassComputationPtr_t comc_;
    };
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_TOOL_HH
