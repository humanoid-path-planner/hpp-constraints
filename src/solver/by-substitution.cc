// Copyright (c) 2017, 2018 CNRS
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

#include <hpp/constraints/solver/by-substitution.hh>
#include <hpp/constraints/solver/impl/by-substitution.hh>
#include <hpp/constraints/solver/impl/hierarchical-iterative.hh>
#include <hpp/constraints/active-set-differentiable-function.hh>

#include <hpp/pinocchio/util.hh>
#include <hpp/pinocchio/configuration.hh>

#include <hpp/constraints/svd.hh>
#include <hpp/constraints/macros.hh>

namespace hpp {
  namespace constraints {
    namespace solver {
      namespace lineSearch {
        template bool Constant::operator()
          (const BySubstitution& solver, vectorOut_t arg, vectorOut_t darg);

        template bool Backtracking::operator()
          (const BySubstitution& solver, vectorOut_t arg, vectorOut_t darg);

        template bool FixedSequence::operator()
          (const BySubstitution& solver, vectorOut_t arg, vectorOut_t darg);

        template bool ErrorNormBased::operator()
          (const BySubstitution& solver, vectorOut_t arg, vectorOut_t darg);
      } // namespace lineSearch

      DifferentiableFunctionPtr_t activeSetFunction
      (const DifferentiableFunctionPtr_t& function,
       const segments_t& pdofs)
      {
        if (pdofs.empty()) return function;
        return ActiveSetDifferentiableFunctionPtr_t
          (new ActiveSetDifferentiableFunction(function, pdofs));
      }

      BySubstitution::BySubstitution (const LiegroupSpacePtr_t& configSpace) :
        HierarchicalIterative(configSpace),
        explicit_ (configSpace->nq (), configSpace->nv ()),
        JeExpanded_ (configSpace->nv (), configSpace->nv ())
      {}

      BySubstitution::BySubstitution (const BySubstitution& other) :
        HierarchicalIterative (other), explicit_ (other.explicit_),
        Je_ (other.Je_), JeExpanded_ (other.JeExpanded_)
      {
        // TODO remove me
        for (LockedJoints_t::const_iterator it = lockedJoints_.begin ();
             it != lockedJoints_.end (); ++it) {
          LockedJointPtr_t lj = HPP_STATIC_PTR_CAST
            (LockedJoint, (*it)->copy ());
          if (!explicitConstraintSet().replace
              ((*it)->explicitFunction(), lj->explicitFunction()))
            throw std::runtime_error
              ("Could not replace lockedJoint function");
        }
      }

      bool BySubstitution::add (const ImplicitPtr_t& nm,
                                const segments_t& passiveDofs,
                                const std::size_t priority)
      {
        if (contains (nm)) {
          hppDout (error, "Constraint " << nm->functionPtr()->name ()
                   << " already in " << this->name () << "." << std::endl);
          return false;
        }
        ComparisonTypes_t types = nm->comparisonType();

        LockedJointPtr_t lj = HPP_DYNAMIC_PTR_CAST (LockedJoint, nm);
        assert (!lj);

        bool addedAsExplicit = false;
        ExplicitPtr_t enm (HPP_DYNAMIC_PTR_CAST (Explicit, nm));
        if (enm) {
          addedAsExplicit = explicitConstraintSet().add
            (enm->explicitFunction(),
             Eigen::RowBlockIndices(enm->inputConf()),
             Eigen::RowBlockIndices(enm->outputConf()),
             Eigen::ColBlockIndices(enm->inputVelocity()),
             Eigen::RowBlockIndices(enm->outputVelocity()),
             types) >= 0;
          if (addedAsExplicit && enm->outputFunction() &&
              enm->outputFunctionInverse()) {
            bool ok = explicitConstraintSet().setG
              (enm->explicitFunction(),
               enm->outputFunction(), enm->outputFunctionInverse());
            assert (ok);
          }
          if (!addedAsExplicit) {
            hppDout (info, "Could not treat " <<
                     enm->explicitFunction()->name()
                     << " as an explicit function.");
          }
        }

        if (!addedAsExplicit) {
          HierarchicalIterative::add (activeSetFunction(nm->functionPtr(),
                                                        passiveDofs), priority,
                                      types);
          // add (Implicit::create
          //      (activeSetFunction(nm->functionPtr(), passiveDofs), types),
          //      segments_t (0), priority);
        } else {
          hppDout (info, "Numerical constraint added as explicit function: "
                   << enm->explicitFunction()->name() << "with "
                   << "input conf " << Eigen::RowBlockIndices(enm->inputConf())
                   << "input vel" << Eigen::RowBlockIndices
                   (enm->inputVelocity())
                   << "output conf " << Eigen::RowBlockIndices
                   (enm->outputConf())
                   << "output vel " << Eigen::RowBlockIndices
                   (enm->outputVelocity()));
          explicitConstraintSetHasChanged();
        }
        hppDout (info, "Constraints " << name() << " has dimension "
                 << dimension());

        functions_.push_back (nm);
        return true;
      }

      void BySubstitution::add (const LockedJointPtr_t& lockedJoint)
      {
        if (lockedJoint->numberDof () == 0) return;
        // If the same dof is already locked, replace by new value
        for (LockedJoints_t::iterator itLock = lockedJoints_.begin ();
             itLock != lockedJoints_.end (); ++itLock) {
          if (lockedJoint->rankInVelocity () == (*itLock)->rankInVelocity ()) {
            if (!explicitConstraintSet().replace
                ((*itLock)->explicitFunction(),
                 lockedJoint->explicitFunction ()))
              {
                throw std::runtime_error
                  ("Could not replace lockedJoint function " +
                   lockedJoint->jointName ());
              }
            *itLock = lockedJoint;
            return;
          }
        }

        ComparisonTypes_t types = lockedJoint->comparisonType();

        bool added = explicitConstraintSet().add
          (lockedJoint->explicitFunction(),
           Eigen::RowBlockIndices(lockedJoint->inputConf()),
           Eigen::RowBlockIndices(lockedJoint->outputConf()),
           Eigen::ColBlockIndices(lockedJoint->inputVelocity()),
           Eigen::RowBlockIndices(lockedJoint->outputVelocity()),
           types) >= 0;

        if (!added) {
          throw std::runtime_error("Could not add lockedJoint function " +
                                   lockedJoint->jointName ());
        }
        if (added) {
          explicitConstraintSet().rightHandSide
            (lockedJoint->explicitFunction(), lockedJoint->rightHandSide());
        }
        explicitConstraintSetHasChanged();

        lockedJoints_.push_back (lockedJoint);
        hppDout (info, "add locked joint " << lockedJoint->jointName ()
                 << " rank in velocity: " << lockedJoint->rankInVelocity ()
                 << ", size: " << lockedJoint->numberDof ());
        hppDout (info, "Intervals: "
                 << explicitConstraintSet().outDers());
        hppDout (info, "Constraints " << name() << " has dimension "
                 << dimension());
      }

      void BySubstitution::explicitConstraintSetHasChanged()
      {
        // set free variables to indices that are not output of the explicit
        // constraint.
        freeVariables (explicit_.notOutDers ().transpose ());
      }

      segments_t BySubstitution::implicitDof () const
      {
        const Eigen::MatrixXi& ioDep = explicit_.inOutDependencies();
        const Eigen::VectorXi& derF = explicit_.derFunction();
        ArrayXb adp (activeDerivativeParameters());
        Eigen::VectorXi out (Eigen::VectorXi::Zero(adp.size()));

        for (size_type i = 0; i < adp.size(); ++i) {
          if (adp(i)) {
            if (derF[i] >= 0) {
              out += ioDep.row(derF[i]);
              out(i) = 0;
            } else
              out(i) += 1;
          }
        }
        return BlockIndex::fromLogicalExpression(out.array().cast<bool>());
      }

      // Note that the jacobian of the implicit constraints have already
      // been computed by computeValue <true>
      // The Jacobian of the implicit constraint of priority i is stored in
      // datas_ [i].jacobian
      void BySubstitution::updateJacobian (vectorIn_t arg) const
      {
        if (explicit_.inDers().nbCols() == 0) return;
        /*                                ------
                         /   in          in u out \
                         |                        |
                   Je_ = |   df                   |
                         |  ---- (qin)      0     |
                         \  dqin                  /
        */
        explicit_.jacobian(JeExpanded_, arg);
        Je_ = explicit_.jacobianNotOutToOut (JeExpanded_);

        hppDnum (info, "Jacobian of explicit system is" << iendl <<
                 setpyformat << pretty_print(Je_));

        for (std::size_t i = 0; i < stacks_.size (); ++i) {
          Data& d = datas_[i];
          hppDnum (info, "Jacobian of stack " << i << " before update:" << iendl
                   << pretty_print(d.reducedJ) << iendl
                   << "Jacobian of explicit variable of stack " << i << ":" << iendl
                   << pretty_print(explicit_.outDers().transpose().rview(d.jacobian).
                                   eval()));
          d.reducedJ.noalias() += Eigen::MatrixBlocksRef<>
            (d.activeRowsOfJ.keepRows(), explicit_.outDers())
            .rview(d.jacobian).eval()
            * Je_;
          hppDnum (info, "Jacobian of stack " << i << " after update:" << iendl
                   << pretty_print(d.reducedJ) << unsetpyformat);
        }
      }

      void BySubstitution::computeActiveRowsOfJ (std::size_t iStack)
      {
        Data& d = datas_[iStack];
        const DifferentiableFunctionStack& f = stacks_[iStack];
        const DifferentiableFunctionStack::Functions_t& fs = f.functions();
        std::size_t row = 0;

        /// ADP: Active Derivative Param
        Eigen::MatrixXi explicitIOdep = explicit_.inOutDofDependencies();
        assert ((explicitIOdep.array() >= 0).all());

        typedef Eigen::MatrixBlocks<false, false> BlockIndices;

        ArrayXb adpF, adpC;
        BlockIndices::segments_t rows;
        for (std::size_t i = 0; i < fs.size (); ++i) {
          bool active;

          // Test on the variable left free by the explicit solver.
          adpF = freeVariables_.rview
            (fs[i]->activeDerivativeParameters().matrix()).eval().array();
          active = adpF.any();
          if (!active && explicitIOdep.size() > 0) {
            // Test on the variable constrained by the explicit solver.
            adpC = explicit_.outDers().rview
              (fs[i]->activeDerivativeParameters().matrix()).eval().array();
            adpF = (explicitIOdep.transpose() * adpC.cast<int>().matrix()).
              array().cast<bool>();
            active = adpF.any();
          }
          if (active) // If at least one element of adp is true
            rows.push_back (BlockIndices::segment_t
                            (row, fs[i]->outputDerivativeSize()));
          row += fs[i]->outputDerivativeSize();
        }
        d.activeRowsOfJ = Eigen::MatrixBlocks<false,false>
          (rows, freeVariables_.m_rows);
        d.activeRowsOfJ.updateRows<true, true, true>();
      }

      void BySubstitution::projectVectorOnKernel
      (ConfigurationIn_t arg, vectorIn_t darg, ConfigurationOut_t result) const
      {
        if (functions_.empty ()) {
          result = darg;
          return;
        }
        computeValue<true> (arg);
        updateJacobian(arg);
        getReducedJacobian (reducedJ_);

        svd_.compute (reducedJ_);

        dqSmall_ = freeVariables_.rview(darg);

        size_type rank = svd_.rank();
        vector_t tmp (getV1(svd_, rank).adjoint() * dqSmall_);
        dqSmall_.noalias() -= getV1(svd_, rank) * tmp;

        freeVariables_.lview(result) = dqSmall_;
      }

      void BySubstitution::projectOnKernel (ConfigurationIn_t from,
                                            ConfigurationIn_t to,
                                            ConfigurationOut_t result)
      {
        // TODO equivalent
        if (functions_.empty ()) {
          result = to;
          return;
        }
        typedef pinocchio::LiegroupElement Lge_t;
        typedef pinocchio::LiegroupConstElementRef LgeConstRef_t;
        LgeConstRef_t O (from, configSpace_);
        LgeConstRef_t M (to, configSpace_);
        OM_ = M - O;

        projectVectorOnKernel (from, OM_, OP_);

        Lge_t P (O + OP_);
        saturate_ (P.vector (), result, saturation_);
      }

      std::ostream& BySubstitution::print (std::ostream& os) const
      {
        os << "BySubstitution" << incendl;
        HierarchicalIterative::print (os) << iendl;
        explicit_.print (os) << decindent;
        return os;
      }

      template BySubstitution::Status BySubstitution::impl_solve
      (vectorOut_t arg, lineSearch::Constant       lineSearch) const;
      template BySubstitution::Status BySubstitution::impl_solve
      (vectorOut_t arg, lineSearch::Backtracking   lineSearch) const;
      template BySubstitution::Status BySubstitution::impl_solve
      (vectorOut_t arg, lineSearch::FixedSequence  lineSearch) const;
      template BySubstitution::Status BySubstitution::impl_solve
      (vectorOut_t arg, lineSearch::ErrorNormBased lineSearch) const;
    } // namespace solver
  } // namespace constraints
} // namespace hpp
