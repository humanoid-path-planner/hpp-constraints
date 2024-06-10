// Copyright (c) 2017, 2018 CNRS
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#include <boost/serialization/nvp.hpp>
#include <hpp/constraints/active-set-differentiable-function.hh>
#include <hpp/constraints/macros.hh>
#include <hpp/constraints/solver/by-substitution.hh>
#include <hpp/constraints/solver/impl/by-substitution.hh>
#include <hpp/constraints/solver/impl/hierarchical-iterative.hh>
#include <hpp/constraints/svd.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/liegroup-element.hh>
#include <hpp/pinocchio/util.hh>
#include <hpp/util/serialization.hh>

namespace hpp {
namespace constraints {
namespace solver {
namespace lineSearch {
template bool Constant::operator()(const BySubstitution& solver,
                                   vectorOut_t arg, vectorOut_t darg);

template bool Backtracking::operator()(const BySubstitution& solver,
                                       vectorOut_t arg, vectorOut_t darg);

template bool FixedSequence::operator()(const BySubstitution& solver,
                                        vectorOut_t arg, vectorOut_t darg);

template bool ErrorNormBased::operator()(const BySubstitution& solver,
                                         vectorOut_t arg, vectorOut_t darg);
}  // namespace lineSearch

BySubstitution::BySubstitution(const LiegroupSpacePtr_t& configSpace)
    : HierarchicalIterative(configSpace),
      explicit_(configSpace),
      JeExpanded_(configSpace->nv(), configSpace->nv()) {}

BySubstitution::BySubstitution(const BySubstitution& other)
    : HierarchicalIterative(other),
      explicit_(other.explicit_),
      Je_(other.Je_),
      JeExpanded_(other.JeExpanded_) {
  for (NumericalConstraints_t::iterator it(constraints_.begin());
       it != constraints_.end(); ++it) {
    // assert (contains (*it));
  }
}

bool BySubstitution::add(const ImplicitPtr_t& nm, const std::size_t& priority) {
  if (contains(nm)) {
    hppDout(error, "Constraint " << nm->functionPtr()->name()
                                 << " already in BySubstitution solver."
                                 << std::endl);
    return false;
  }
  ComparisonTypes_t types = nm->comparisonType();

  bool addedAsExplicit = false;
  ExplicitPtr_t enm(HPP_DYNAMIC_PTR_CAST(Explicit, nm));
  if (enm) {
    addedAsExplicit = explicitConstraintSet().add(enm) >= 0;
    if (!addedAsExplicit) {
      hppDout(info, "Could not treat " << enm->explicitFunction()->name()
                                       << " as an explicit function.");
    }
  }

  if (addedAsExplicit) {
    // If added as explicit, add to the list of constraint of Hierarchical
    // iterative
    constraints_.push_back(nm);
    hppDout(info,
            "Numerical constraint added as explicit function: "
                << enm->explicitFunction()->name() << "with "
                << "input conf " << Eigen::RowBlockIndices(enm->inputConf())
                << "input vel" << Eigen::RowBlockIndices(enm->inputVelocity())
                << "output conf " << Eigen::RowBlockIndices(enm->outputConf())
                << "output vel "
                << Eigen::RowBlockIndices(enm->outputVelocity()));
    explicitConstraintSetHasChanged();
  } else
    HierarchicalIterative::add(nm, priority);
  hppDout(info, "Constraint has dimension " << dimension());

  assert(contains(nm));
  return true;
}

void BySubstitution::explicitConstraintSetHasChanged() {
  // set free variables to indices that are not output of the explicit
  // constraint.
  freeVariables(explicit_.notOutDers().transpose());
}

segments_t BySubstitution::implicitDof() const {
  const Eigen::MatrixXi& ioDep = explicit_.inOutDependencies();
  const Eigen::VectorXi& derF = explicit_.derFunction();
  ArrayXb adp(activeDerivativeParameters());
  Eigen::VectorXi out(Eigen::VectorXi::Zero(adp.size()));

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
void BySubstitution::updateJacobian(vectorIn_t arg) const {
  if (explicit_.inDers().nbCols() == 0) return;
  /*                                ------
                   /   in          in u out \
                   |                        |
             Je_ = |   df                   |
                   |  ---- (qin)      0     |
                   \  dqin                  /
  */
  explicit_.jacobian(JeExpanded_, arg);
  Je_ = explicit_.jacobianNotOutToOut(JeExpanded_);

  hppDnum(info, "Jacobian of explicit system is" << iendl << setpyformat
                                                 << pretty_print(Je_));

  for (std::size_t i = 0; i < stacks_.size(); ++i) {
    Data& d = datas_[i];
    hppDnum(
        info,
        "Jacobian of stack "
            << i << " before update:" << iendl << pretty_print(d.reducedJ)
            << iendl << "Jacobian of explicit variable of stack " << i << ":"
            << iendl
            << pretty_print(
                   explicit_.outDers().transpose().rview(d.jacobian).eval()));
    d.reducedJ.noalias() += Eigen::MatrixBlocksRef<>(d.activeRowsOfJ.keepRows(),
                                                     explicit_.outDers())
                                .rview(d.jacobian)
                                .eval() *
                            Je_;
    hppDnum(info, "Jacobian of stack " << i << " after update:" << iendl
                                       << pretty_print(d.reducedJ)
                                       << unsetpyformat);
  }
}

void BySubstitution::computeActiveRowsOfJ(std::size_t iStack) {
  Data& d = datas_[iStack];
  const ImplicitConstraintSet::Implicits_t constraints(
      stacks_[iStack].constraints());
  std::size_t row = 0;

  /// ADP: Active Derivative Param
  Eigen::MatrixXi explicitIOdep = explicit_.inOutDofDependencies();
  assert((explicitIOdep.array() >= 0).all());

  typedef Eigen::MatrixBlocks<false, false> BlockIndices;

  ArrayXb adpF, adpC;
  BlockIndices::segments_t rows;
  for (std::size_t i = 0; i < constraints.size(); ++i) {
    bool active;

    // Test on the variable left free by the explicit solver.
    adpF = freeVariables_
               .rview(constraints[i]
                          ->function()
                          .activeDerivativeParameters()
                          .matrix())
               .eval()
               .array();
    active = adpF.any();
    if (!active && explicitIOdep.size() > 0) {
      // Test on the variable constrained by the explicit solver.
      adpC = explicit_.outDers()
                 .rview(constraints[i]
                            ->function()
                            .activeDerivativeParameters()
                            .matrix())
                 .eval()
                 .array();
      adpF = (explicitIOdep.transpose() * adpC.cast<int>().matrix())
                 .array()
                 .cast<bool>();
      active = adpF.any();
    }
    if (active) {  // If at least one element of adp is true
      for (const segment_t s : constraints[i]->activeRows()) {
        rows.emplace_back(s.first + row, s.second);
      }
    }
    row += constraints[i]->function().outputDerivativeSize();
  }
  d.activeRowsOfJ =
      Eigen::MatrixBlocks<false, false>(rows, freeVariables_.m_rows);
  d.activeRowsOfJ.updateRows<true, true, true>();
}

void BySubstitution::projectVectorOnKernel(ConfigurationIn_t arg,
                                           vectorIn_t darg,
                                           vectorOut_t result) const {
  if (constraints_.empty() || reducedDimension() == 0) {
    result = darg;
    return;
  }
  computeValue<true>(arg);
  updateJacobian(arg);
  getReducedJacobian(reducedJ_);

  svd_.compute(reducedJ_);

  // TODO the output of explicit solver should be set to zero ?
  dqSmall_ = freeVariables_.rview(darg);

  size_type rank = svd_.rank();
  vector_t tmp(getV1(svd_, rank).adjoint() * dqSmall_);
  dqSmall_.noalias() -= getV1(svd_, rank) * tmp;

  // Otherwise two uninitialized values may sum up to NaN
  result.setZero();
  freeVariables_.lview(result) = dqSmall_;
}

void BySubstitution::projectOnKernel(ConfigurationIn_t from,
                                     ConfigurationIn_t to,
                                     ConfigurationOut_t result) {
  // TODO equivalent
  if (constraints_.empty()) {
    result = to;
    return;
  }
  typedef pinocchio::LiegroupElement Lge_t;
  typedef pinocchio::LiegroupElementConstRef LgeConstRef_t;
  LgeConstRef_t O(from, configSpace_);
  LgeConstRef_t M(to, configSpace_);
  OM_ = M - O;

  projectVectorOnKernel(from, OM_, OP_);

  assert(hpp::pinocchio::checkNormalized(M));
  Lge_t P(O + OP_);
  assert(hpp::pinocchio::checkNormalized(P));
  saturate_->saturate(P.vector(), result, saturation_);
}

std::ostream& BySubstitution::print(std::ostream& os) const {
  os << "BySubstitution" << incendl;
  HierarchicalIterative::print(os) << iendl;
  explicit_.print(os) << decindent;
  return os;
}

vector_t BySubstitution::rightHandSideFromConfig(ConfigurationIn_t config) {
  const size_type top = parent_t::rightHandSideSize();
  const size_type bot = explicit_.rightHandSideSize();
  vector_t rhs(top + bot);
  rhs.head(top) = parent_t::rightHandSideFromConfig(config);
  rhs.tail(bot) = explicit_.rightHandSideFromInput(config);
  return rhs;
}

bool BySubstitution::rightHandSideFromConfig(const ImplicitPtr_t& constraint,
                                             ConfigurationIn_t config) {
  if (parent_t::rightHandSideFromConfig(constraint, config)) return true;
  ExplicitPtr_t exp(HPP_DYNAMIC_PTR_CAST(Explicit, constraint));
  if (exp) {
    return explicit_.rightHandSideFromInput(exp, config);
  }
  return false;
}

bool BySubstitution::rightHandSide(const ImplicitPtr_t& constraint,
                                   vectorIn_t rhs) {
  if (parent_t::rightHandSide(constraint, rhs)) return true;
  ExplicitPtr_t exp(HPP_DYNAMIC_PTR_CAST(Explicit, constraint));
  if (exp) {
    return explicit_.rightHandSide(exp, rhs);
  }
  return false;
}

bool BySubstitution::getRightHandSide(const ImplicitPtr_t& constraint,
                                      vectorOut_t rhs) const {
  if (parent_t::getRightHandSide(constraint, rhs)) return true;
  ExplicitPtr_t exp(HPP_DYNAMIC_PTR_CAST(Explicit, constraint));
  if (exp) {
    return explicit_.getRightHandSide(exp, rhs);
  }
  return false;
}

void BySubstitution::rightHandSide(vectorIn_t rhs) {
  const size_type top = parent_t::rightHandSideSize();
  const size_type bot = explicit_.rightHandSideSize();
  parent_t::rightHandSide(rhs.head(top));
  explicit_.rightHandSide(rhs.head(bot));
}

vector_t BySubstitution::rightHandSide() const {
  const size_type top = parent_t::rightHandSideSize();
  const size_type bot = explicit_.rightHandSideSize();
  vector_t rhs(top + bot);
  rhs.head(top) = parent_t::rightHandSide();
  rhs.tail(bot) = explicit_.rightHandSide();
  return rhs;
}

size_type BySubstitution::rightHandSideSize() const {
  const size_type top = parent_t::rightHandSideSize();
  const size_type bot = explicit_.rightHandSideSize();
  return top + bot;
}

bool BySubstitution::isConstraintSatisfied(const ImplicitPtr_t& constraint,
                                           vectorIn_t arg, vectorOut_t error,
                                           bool& constraintFound) const {
  constraintFound = false;
  bool satisfied(
      parent_t::isConstraintSatisfied(constraint, arg, error, constraintFound));
  if (constraintFound) return satisfied;
  return explicit_.isConstraintSatisfied(constraint, arg, error,
                                         constraintFound);
}

template <class Archive>
void BySubstitution::load(Archive& ar, const unsigned int version) {
  using namespace boost::serialization;
  (void)version;
  LiegroupSpacePtr_t space;
  ar& BOOST_SERIALIZATION_NVP(space);
  explicit_.init(space);
  ar& make_nvp("base", base_object<HierarchicalIterative>(*this));
}

template <class Archive>
void BySubstitution::save(Archive& ar, const unsigned int version) const {
  using namespace boost::serialization;
  (void)version;
  LiegroupSpacePtr_t space(explicit_.configSpace());
  ar& BOOST_SERIALIZATION_NVP(space);
  ar& make_nvp("base", base_object<HierarchicalIterative>(*this));
}

HPP_SERIALIZATION_SPLIT_IMPLEMENT(BySubstitution);

template BySubstitution::Status BySubstitution::impl_solve(
    vectorOut_t arg, bool optimize, lineSearch::Constant lineSearch) const;
template BySubstitution::Status BySubstitution::impl_solve(
    vectorOut_t arg, bool optimize, lineSearch::Backtracking lineSearch) const;
template BySubstitution::Status BySubstitution::impl_solve(
    vectorOut_t arg, bool optimize, lineSearch::FixedSequence lineSearch) const;
template BySubstitution::Status BySubstitution::impl_solve(
    vectorOut_t arg, bool optimize,
    lineSearch::ErrorNormBased lineSearch) const;
}  // namespace solver
}  // namespace constraints
}  // namespace hpp

BOOST_CLASS_EXPORT(hpp::constraints::solver::BySubstitution)
