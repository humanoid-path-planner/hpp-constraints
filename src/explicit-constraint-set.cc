// Copyright (c) 2017, Joseph Mirabel
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

#include <hpp/constraints/explicit-constraint-set.hh>
#include <hpp/constraints/explicit.hh>
#include <hpp/constraints/matrix-view.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/liegroup.hh>
#include <hpp/pinocchio/util.hh>
#include <hpp/util/indent.hh>
#include <queue>

namespace hpp {
namespace constraints {
typedef Eigen::MatrixBlocksRef<false, false> MatrixBlocksRef;

namespace {
/// Append all indices of a set of segments to a queue of indices
void append(const Eigen::RowBlockIndices& rbi, std::queue<size_type>& q) {
  for (std::size_t i = 0; i < rbi.indices().size(); ++i)
    for (size_type j = 0; j < rbi.indices()[i].second; ++j)
      q.push(rbi.indices()[i].first + j);
}
}  // namespace

Eigen::ColBlockIndices ExplicitConstraintSet::activeParameters() const {
  return inArgs_.transpose();
}

const Eigen::ColBlockIndices&
ExplicitConstraintSet::activeDerivativeParameters() const {
  return inDers_;
}

bool ExplicitConstraintSet::solve(vectorOut_t arg) const {
  for (std::size_t i = 0; i < data_.size(); ++i) {
    solveExplicitConstraint(computationOrder_[i], arg);
  }
  return true;
}

bool ExplicitConstraintSet::isSatisfied(vectorIn_t arg, vectorOut_t error,
                                        value_type errorThreshold) const {
  // Recover default value
  if (errorThreshold == -1) errorThreshold = errorThreshold_;
  value_type squaredNorm = 0;

  size_type row = 0;
  for (std::size_t i = 0; i < data_.size(); ++i) {
    const Data& d(data_[i]);
    const DifferentiableFunction& h(d.constraint->function());
    h.value(d.h_value, arg);
    size_type nRows(h.outputSpace()->nv());
    assert(*(d.h_value.space()) == *(d.rhs_implicit.space()));
    error.segment(row, nRows) = d.h_value - d.rhs_implicit;
    squaredNorm =
        std::max(squaredNorm, error.segment(row, nRows).squaredNorm());
    row += nRows;
  }
  assert(row == error.size());
  hppDout(info, "Max squared error norm is " << squaredNorm);
  return squaredNorm < errorThreshold * errorThreshold;
}

bool ExplicitConstraintSet::isSatisfied(vectorIn_t arg,
                                        value_type errorThreshold) const {
  // Recover default value
  if (errorThreshold == -1) errorThreshold = errorThreshold_;
  diffSmall_.resize(errorSize_);
  return isSatisfied(arg, diffSmall_, errorThreshold);
}

bool ExplicitConstraintSet::isConstraintSatisfied(
    const ImplicitPtr_t& constraint, vectorIn_t arg, vectorOut_t error,
    bool& constraintFound) const {
  value_type squaredNorm = 0;
  constraintFound = false;
  for (std::size_t i = 0; i < data_.size(); ++i) {
    const Data& d(data_[i]);
    if (d.constraint->functionPtr() == constraint->functionPtr() ||
        d.constraint->function() == constraint->function()) {
      const DifferentiableFunction& h(d.constraint->function());
      h.value(d.h_value, arg);
      assert(error.size() == h.outputSpace()->nv());
      assert(*(d.h_value.space()) == *(d.rhs_implicit.space()));
      error = d.h_value - d.rhs_implicit;
      squaredNorm = error.squaredNorm();
      constraintFound = true;
      return squaredNorm < errorThreshold_ * errorThreshold_;
    }
  }
  return false;
}

size_type size(const segments_t& intervals) {
  size_type res = 0;
  for (std::size_t i = 0; i < intervals.size(); ++i) {
    res += intervals[i].second;
  }
  return res;
}

ExplicitConstraintSet::Data::Data(const ExplicitPtr_t& _constraint)
    : constraint(_constraint),
      rhs_implicit(_constraint->functionPtr()->outputSpace()->neutral()),
      h_value(_constraint->functionPtr()->outputSpace()),
      f_value(_constraint->explicitFunction()->outputSpace()),
      res_qout(_constraint->explicitFunction()->outputSpace()) {
  jacobian.resize(_constraint->explicitFunction()->outputDerivativeSize(),
                  _constraint->explicitFunction()->inputDerivativeSize());
  for (std::size_t i = 0; i < constraint->comparisonType().size(); ++i) {
    if (constraint->comparisonType()[i] == Equality) {
      equalityIndices.addRow(i, 1);
    }
  }
  equalityIndices.updateRows<true, true, true>();
}

size_type ExplicitConstraintSet::add(const ExplicitPtr_t& constraint) {
  assert(constraint->outputConf().size() == 1 &&
         "Only contiguous function output is supported.");
  assert(constraint->outputVelocity().size() == 1 &&
         "Only contiguous function output is supported.");
  typedef Eigen::RowBlockIndices RowBlockIndices;
  typedef Eigen::ColBlockIndices ColBlockIndices;
  const RowBlockIndices::segment_t& outIdx(constraint->outputConf()[0]);
  const RowBlockIndices::segment_t& outDerIdx(constraint->outputVelocity()[0]);

  // TODO This sanity check should not be necessary
  // It is done in the while loop below
  // Sanity check: is it explicit ?
  for (std::size_t i = 0; i < constraint->inputConf().size(); ++i)
    if (BlockIndex::overlap(constraint->inputConf()[i], outIdx)) {
      size_type f0(constraint->inputConf()[i].first);
      size_type s0(f0 + constraint->inputConf()[i].second - 1);
      size_type f1(outIdx.first);
      size_type s1(f1 + outIdx.second - 1);
      std::ostringstream oss;
      oss << "Explicit constraint \"" << constraint->function().name()
          << "\" is malformed.";
      oss << " input [" << f0 << "," << s0 << "] and output [" << f1 << ","
          << s1 << "] segments overlap.";
      throw std::logic_error(oss.str().c_str());
    }
  // Sanity check: Comparison type must be either EqualToZero or Equality
  const ComparisonTypes_t& comp(constraint->comparisonType());
  assert(comp.size() ==
         (std::size_t)constraint->function().outputDerivativeSize());
  for (std::size_t i = 0; i < comp.size(); ++i)
    if (comp[i] != EqualToZero && comp[i] != Equality) return -1;
  // Check that no other function already computes its outputs.
  RowBlockIndices outArg(constraint->outputConf());
  if ((outArg.rview(argFunction_).eval().array() >= 0).any()) return -1;
  // Check that it does not insert a loop.
  std::queue<size_type> idxArg;
  // Put all indices of inArg into idxArg
  RowBlockIndices inArg(constraint->inputConf());
  append(inArg, idxArg);
  while (!idxArg.empty()) {
    // iArg must be computed before f
    size_type iArg = idxArg.front();
    idxArg.pop();
    // iArg is an output of the new function -> new function should be
    // computed earlier -> incompatible.
    if (iArg >= outIdx.first && iArg < outIdx.first + outIdx.second) return -1;
    // iArg is not computed by any function
    if (argFunction_[iArg] < 0) continue;
    const Data& d = data_[argFunction_[iArg]];
    append(d.constraint->inputConf(), idxArg);
  }
  size_type nq(configSpace_->nq());
  size_type nv(configSpace_->nv());
  // Add the function
  int idx = int(data_.size());
  outArg.lview(argFunction_).setConstant(idx);
  RowBlockIndices(constraint->outputVelocity())
      .lview(derFunction_)
      .setConstant(idx);
  data_.push_back(Data(constraint));
  errorSize_ += data_.back().rhs_implicit.space()->nv();

  // Update the free dofs
  outArgs_.addRow(outIdx.first, outIdx.second);
  outArgs_.updateIndices<true, true, true>();
  notOutArgs_ = RowBlockIndices(
      BlockIndex::difference(BlockIndex::segment_t(0, nq), outArgs_.indices()));

  BlockIndex::add(inArgs_.m_rows, inArg.rows());
  inArgs_ =
      RowBlockIndices(BlockIndex::difference(inArgs_.rows(), outArgs_.rows()));
  // should be sorted already
  inArgs_.updateIndices<false, true, true>();

  outDers_.addRow(outDerIdx.first, outDerIdx.second);
  outDers_.updateIndices<true, true, true>();
  notOutDers_ = ColBlockIndices(
      BlockIndex::difference(BlockIndex::segment_t(0, nv), outDers_.indices()));

  BlockIndex::add(inDers_.m_cols,
                  ColBlockIndices(constraint->inputVelocity()).cols());
  inDers_ =
      ColBlockIndices(BlockIndex::difference(inDers_.cols(), outDers_.rows()));
  // should be sorted already
  inDers_.updateIndices<false, true, true>();

  /// Computation order
  std::size_t order = 0;
  computationOrder_.resize(data_.size());
  inOutDependencies_ = Eigen::MatrixXi::Zero(data_.size(), nv);
  Computed_t computed(data_.size(), false);
  for (std::size_t i = 0; i < data_.size(); ++i)
    computeOrder(i, order, computed);
  assert(order == data_.size());
  return data_.size() - 1;
}

bool ExplicitConstraintSet::contains(
    const ExplicitPtr_t& numericalConstraint) const {
  for (std::vector<Data>::const_iterator it = data_.begin(); it != data_.end();
       ++it) {
    if ((it->constraint == numericalConstraint) ||
        (*(it->constraint) == *numericalConstraint) ||
        (*(it->constraint->functionPtr()) ==
         *numericalConstraint->functionPtr()))
      return true;
  }
  return false;
}

void ExplicitConstraintSet::solveExplicitConstraint(const std::size_t& iF,
                                                    vectorOut_t arg) const {
  const Data& d = data_[iF];
  // Compute this function
  d.qin = RowBlockIndices(d.constraint->inputConf()).rview(arg);
  d.constraint->outputValue(d.res_qout, d.qin, d.rhs_implicit);
  RowBlockIndices(d.constraint->outputConf()).lview(arg) = d.res_qout.vector();
}

void ExplicitConstraintSet::jacobian(matrixOut_t jacobian,
                                     vectorIn_t arg) const {
  // TODO this could be done only on the complement of inDers_
  jacobian.setZero();
  MatrixBlocksRef(notOutDers_, notOutDers_).lview(jacobian).setIdentity();
  // Compute the function jacobians
  for (std::size_t i = 0; i < data_.size(); ++i) {
    const Data& d = data_[i];
    d.qin = RowBlockIndices(d.constraint->inputConf()).rview(arg);
    // Compute Jacobian of f(qin) + rhs
    // with respect to qin.
    d.constraint->jacobianOutputValue(d.qin, d.f_value, d.rhs_implicit,
                                      d.jacobian);
  }
  for (std::size_t i = 0; i < data_.size(); ++i) {
    computeJacobian(computationOrder_[i], jacobian);
  }
}

void ExplicitConstraintSet::computeJacobian(const std::size_t& iE,
                                            matrixOut_t J) const {
  const Data& d = data_[iE];
  ColBlockIndices inDer(d.constraint->inputVelocity());
  matrix_t Jin(MatrixBlocksRef(inDer, inDers_).rview(J));
  // Jout = d.jacobian * Jin
  MatrixBlocksRef(RowBlockIndices(d.constraint->outputVelocity()), inDers_)
      .lview(J) = d.jacobian * Jin;
}

void ExplicitConstraintSet::computeOrder(const std::size_t& iE,
                                         std::size_t& iOrder,
                                         Computed_t& computed) {
  if (computed[iE]) return;
  const Data& d = data_[iE];
  for (std::size_t i = 0; i < d.constraint->inputVelocity().size(); ++i) {
    const BlockIndex::segment_t& segment(d.constraint->inputVelocity()[i]);
    for (size_type j = 0; j < segment.second; ++j) {
      if (derFunction_[segment.first + j] < 0) {
        inOutDependencies_(iE, segment.first + j) += 1;
      } else {
        assert((std::size_t)derFunction_[segment.first + j] < data_.size());
        computeOrder(derFunction_[segment.first + j], iOrder, computed);
        inOutDependencies_.row(iE) +=
            inOutDependencies_.row(derFunction_[segment.first + j]);
      }
    }
  }
  computationOrder_[iOrder] = iE;
  ++iOrder;
  computed[iE] = true;
}

vector_t ExplicitConstraintSet::rightHandSideFromInput(vectorIn_t arg) {
  // Compute right hand side of explicit formulation
  for (std::size_t i = 0; i < data_.size(); ++i) {
    rightHandSideFromInput(i, arg);
  }
  return rightHandSide();
}

bool ExplicitConstraintSet::rightHandSideFromInput(
    const ExplicitPtr_t& constraint, vectorIn_t arg) {
  for (std::size_t i = 0; i < data_.size(); ++i) {
    Data& d = data_[i];
    if ((d.constraint == constraint) || (*d.constraint == *constraint) ||
        d.constraint->function() == constraint->function()) {
      rightHandSideFromInput(i, arg);
      return true;
    }
  }
  return false;
}

void ExplicitConstraintSet::rightHandSideFromInput(const size_type& i,
                                                   vectorIn_t arg) {
  Data& d = data_[i];
  // compute right hand side of implicit formulation that might be
  // different (RelativePose)
  d.constraint->function().value(d.h_value, arg);
  vector_t logRhs(log(d.h_value));
  // Equality indices apply on the log of the right hand side
  // This is necessary for constraints built with
  // RelativeTransformationR3xSO3.
  vector_t logRhsImplicit(vector_t::Zero(d.rhs_implicit.space()->nv()));
  d.equalityIndices.lview(logRhsImplicit) = d.equalityIndices.rview(logRhs);
  d.rhs_implicit = d.rhs_implicit.space()->exp(logRhsImplicit);
}

void ExplicitConstraintSet::rightHandSide(vectorIn_t rhs) {
  size_type row = 0;
  for (std::size_t i = 0; i < data_.size(); ++i) {
    Data& d = data_[i];
    // comparison types are applied to the log of the right hand side.
    // Initialize a zero vector of velocity size.
    vector_t logRhsImplicit(vector_t::Zero(d.rhs_implicit.space()->nv()));
    // Build liegroupElement from value extracted from input rhs
    LiegroupElementConstRef inputRhs(
        rhs.segment(row, d.rhs_implicit.space()->nq()), d.rhs_implicit.space());
    d.equalityIndices.lview(logRhsImplicit) =
        d.equalityIndices.rview(log(inputRhs));
    ComparisonTypes_t ct(d.constraint->comparisonType());
    for (std::size_t i = 0; i < ct.size(); ++i) {
      assert(ct[i] == Equality || logRhsImplicit[i] == 0);
    }
    d.rhs_implicit = d.rhs_implicit.space()->exp(logRhsImplicit);
    row += d.rhs_implicit.space()->nq();
  }
  assert(row == rhs.size());
}

bool ExplicitConstraintSet::rightHandSide(const ExplicitPtr_t& constraint,
                                          vectorIn_t rhs) {
  for (std::size_t i = 0; i < data_.size(); ++i) {
    Data& d = data_[i];
    if ((d.constraint == constraint) || (*d.constraint == *constraint) ||
        d.constraint->function() == constraint->function()) {
      rightHandSide(i, rhs);
      return true;
    }
  }
  return false;
}

bool ExplicitConstraintSet::getRightHandSide(const ExplicitPtr_t& constraint,
                                             vectorOut_t rhs) const {
  for (std::size_t i = 0; i < data_.size(); ++i) {
    const Data& d = data_[i];
    if ((d.constraint == constraint) || (*d.constraint == *constraint) ||
        d.constraint->function() == constraint->function()) {
      rhs = d.rhs_implicit.vector();
      return true;
    }
  }
  return false;
}

void ExplicitConstraintSet::rightHandSide(const size_type& i, vectorIn_t rhs) {
  assert(i < (size_type)data_.size());
  Data& d = data_[i];
  LiegroupElementConstRef rhs_lge(rhs, d.rhs_implicit.space());
  vector_t logRhsImplicit(vector_t::Zero(d.rhs_implicit.space()->nv()));
  vector_t logRhsInput(log(rhs_lge));
  d.equalityIndices.lview(logRhsImplicit) =
      d.equalityIndices.rview(logRhsInput);
  d.rhs_implicit = d.rhs_implicit.space()->exp(logRhsImplicit);
  ComparisonTypes_t ct(d.constraint->comparisonType());
  for (std::size_t i = 0; i < ct.size(); ++i) {
    assert(ct[i] == Equality ||
           logRhsInput[i] * logRhsInput[i] < errorThreshold_ * errorThreshold_);
  }
}

vector_t ExplicitConstraintSet::rightHandSide() const {
  vector_t rhs(rightHandSideSize());
  size_type row = 0;
  for (std::size_t i = 0; i < data_.size(); ++i) {
    const Data& d = data_[i];
    const size_type nRows = d.rhs_implicit.space()->nq();
    vector_t::SegmentReturnType seg = rhs.segment(row, nRows);
    seg = d.rhs_implicit.vector();
    row += nRows;
  }
  assert(row == rhs.size());
  return rhs;
}

size_type ExplicitConstraintSet::rightHandSideSize() const {
  size_type rhsSize = 0;
  for (std::size_t i = 0; i < data_.size(); ++i)
    rhsSize += data_[i].rhs_implicit.space()->nq();
  return rhsSize;
}

std::ostream& ExplicitConstraintSet::print(std::ostream& os) const {
  os << "ExplicitConstraintSet, " << data_.size() << " functions." << incendl
     << "Other args: " << notOutArgs_ << iendl << "Params: " << inArgs_
     << " -> " << outArgs_ << iendl << "Dofs: " << inDers_ << " -> " << outDers_
     << iendl << "Functions" << incindent;
  for (std::size_t i = 0; i < data_.size(); ++i) {
    const Data& d = data_[computationOrder_[i]];
    os << iendl << i << ": " << RowBlockIndices(d.constraint->inputConf())
       << " -> " << RowBlockIndices(d.constraint->outputConf()) << incendl
       << *d.constraint->explicitFunction() << iendl
       << "Rhs implicit: " << condensed(d.rhs_implicit.vector()) << iendl
       << "Equality: " << d.equalityIndices << decindent;
  }
  return os << decindent << decindent;
}

void ExplicitConstraintSet::replace(
    const std::map<ImplicitPtr_t, ImplicitPtr_t>& oldToNew) {
  for (auto& d : data_) {
    std::map<ImplicitPtr_t, ImplicitPtr_t>::const_iterator it =
        oldToNew.find(d.constraint);
    assert(it != oldToNew.end());
    assert(HPP_DYNAMIC_PTR_CAST(Explicit, it->second));
    ExplicitPtr_t expl(HPP_STATIC_PTR_CAST(Explicit, it->second));
    d.constraint = expl;
  }
}

Eigen::MatrixXi ExplicitConstraintSet::inOutDofDependencies() const {
  Eigen::MatrixXi iod(nv(), inDers_.nbCols());
  if (inDers_.nbCols() == 0) return iod;
  Eigen::RowVectorXi tmp(inDers_.nbCols());
  for (std::size_t i = 0; i < data_.size(); ++i) {
    const Data& d = data_[i];
    tmp = inDers_.rview(inOutDependencies_.row(i));
    RowBlockIndices(d.constraint->outputVelocity()).lview(iod).rowwise() = tmp;
  }
  return outDers_.rview(iod);
}
}  // namespace constraints
}  // namespace hpp
