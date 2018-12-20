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
        explicit_ (configSpace),
        JeExpanded_ (configSpace->nv (), configSpace->nv ())
      {}

      BySubstitution::BySubstitution (const BySubstitution& other) :
        HierarchicalIterative (other), explicit_ (other.explicit_),
        Je_ (other.Je_), JeExpanded_ (other.JeExpanded_)
      {
      }

      bool BySubstitution::add (const ImplicitPtr_t& nm,
                                const segments_t& passiveDofs,
                                const std::size_t priority)
      {
        if (contains (nm)) {
          hppDout (error, "Constraint " << nm->functionPtr()->name ()
                   << " already in BySubstitution solver." << std::endl);
          return false;
        }
        ComparisonTypes_t types = nm->comparisonType();

        bool addedAsExplicit = false;
        ExplicitPtr_t enm (HPP_DYNAMIC_PTR_CAST (Explicit, nm));
        if (enm) {
          addedAsExplicit = explicitConstraintSet().add (enm) >= 0;
          if (!addedAsExplicit) {
            hppDout (info, "Could not treat " <<
                     enm->explicitFunction()->name()
                     << " as an explicit function.");
          }
        }

        if (addedAsExplicit) {
          // If added as explicit, add to the list of constraint of Hierarchical
          // iterative
          constraints_.push_back (nm);
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
        } else {
          ImplicitPtr_t constraint
            (Implicit::create (activeSetFunction(nm->functionPtr(),
                                                 passiveDofs), types));
          HierarchicalIterative::add (constraint, priority);
          // add (Implicit::create
          //      (activeSetFunction(nm->functionPtr(), passiveDofs), types),
          //      segments_t (0), priority);
        }
        hppDout (info, "Constraint has dimension "
                 << dimension());

        return true;
      }

      void BySubstitution::add (const DifferentiableFunctionPtr_t& f,
                                const std::size_t& priority,
                                const ComparisonTypes_t& comp)
      {
        HierarchicalIterative::add (f, priority, comp);
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
        const ImplicitConstraintSet::Implicits_t constraints
          (stacks_ [iStack].constraints ());
        std::size_t row = 0;

        /// ADP: Active Derivative Param
        Eigen::MatrixXi explicitIOdep = explicit_.inOutDofDependencies();
        assert ((explicitIOdep.array() >= 0).all());

        typedef Eigen::MatrixBlocks<false, false> BlockIndices;

        ArrayXb adpF, adpC;
        BlockIndices::segments_t rows;
        for (std::size_t i = 0; i < constraints.size (); ++i) {
          bool active;

          // Test on the variable left free by the explicit solver.
          adpF = freeVariables_.rview
            (constraints [i]->function ().activeDerivativeParameters().
             matrix()).eval().array();
          active = adpF.any();
          if (!active && explicitIOdep.size() > 0) {
            // Test on the variable constrained by the explicit solver.
            adpC = explicit_.outDers().rview
              (constraints [i]->function ().activeDerivativeParameters().
               matrix()).eval().array();
            adpF = (explicitIOdep.transpose() * adpC.cast<int>().matrix()).
              array().cast<bool>();
            active = adpF.any();
          }
          if (active) // If at least one element of adp is true
            rows.push_back (BlockIndices::segment_t
                            (row, constraints [i]->function ().
                             outputDerivativeSize()));
          row += constraints [i]->function ().outputDerivativeSize();
        }
        d.activeRowsOfJ = Eigen::MatrixBlocks<false,false>
          (rows, freeVariables_.m_rows);
        d.activeRowsOfJ.updateRows<true, true, true>();
      }

      void BySubstitution::projectVectorOnKernel
      (ConfigurationIn_t arg, vectorIn_t darg, ConfigurationOut_t result) const
      {
        if (constraints_.empty ()) {
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
        if (constraints_.empty ()) {
          result = to;
          return;
        }
        typedef pinocchio::LiegroupElement Lge_t;
        typedef pinocchio::LiegroupElementConstRef LgeConstRef_t;
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

      vector_t BySubstitution::rightHandSideFromConfig
      (ConfigurationIn_t config)
      {
        const size_type top = parent_t::rightHandSideSize();
        const size_type bot = explicit_.rightHandSideSize();
        vector_t rhs (top + bot);
        rhs.head(top) = parent_t::rightHandSideFromConfig (config);
        rhs.tail(bot) = explicit_.rightHandSideFromInput (config);
        return rhs;
      }

      bool BySubstitution::rightHandSideFromConfig
      (const ImplicitPtr_t& constraint, ConfigurationIn_t config)
      {
        if (parent_t::rightHandSideFromConfig (constraint, config))
          return true;
        ExplicitPtr_t exp (HPP_DYNAMIC_PTR_CAST (Explicit, constraint));
        if (exp) {
          return explicit_.rightHandSideFromInput (exp, config);
        }
        return false;
      }

      bool BySubstitution::rightHandSide (const ImplicitPtr_t& constraint,
                                          vectorIn_t rhs)
      {
        if (parent_t::rightHandSide (constraint, rhs))
          return true;
        ExplicitPtr_t exp (HPP_DYNAMIC_PTR_CAST (Explicit, constraint));
        if (exp) {
          return explicit_.rightHandSide (exp, rhs);
        }
        return false;
      }

      bool BySubstitution::getRightHandSide (const ImplicitPtr_t& constraint,vectorOut_t rhs)
      {
	if (parent_t::getRightHandSide ( constraint, rhs))
          return true;
        ExplicitPtr_t exp (HPP_DYNAMIC_PTR_CAST (Explicit, constraint));
        if (exp) {
          return explicit_.getRightHandSide ( exp, rhs);
        }
        return false;
      }

      void BySubstitution::rightHandSide (vectorIn_t rhs)
      {
        const size_type top = parent_t::rightHandSideSize();
        const size_type bot = explicit_.rightHandSideSize();
        parent_t::rightHandSide (rhs.head(top));
        explicit_.rightHandSide (rhs.head(bot));
      }

      vector_t BySubstitution::rightHandSide () const
      {
        const size_type top = parent_t::rightHandSideSize();
        const size_type bot = explicit_.rightHandSideSize();
        vector_t rhs (top + bot);
        rhs.head(top) = parent_t::rightHandSide ();
        rhs.tail(bot) = explicit_.rightHandSide ();
        return rhs;
      }

      size_type BySubstitution::rightHandSideSize () const
      {
        const size_type top = parent_t::rightHandSideSize();
        const size_type bot = explicit_.rightHandSideSize();
        return top + bot;
      }

        /// \}

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
