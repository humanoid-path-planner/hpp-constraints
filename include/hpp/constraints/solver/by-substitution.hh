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

#ifndef HPP_CONSTRAINTS_SOLVER_BY_SUBSTITUTION_HH
#define HPP_CONSTRAINTS_SOLVER_BY_SUBSTITUTION_HH

#include <hpp/constraints/config.hh>
#include <hpp/constraints/explicit-constraint-set.hh>
#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/locked-joint.hh>
#include <hpp/constraints/solver/hierarchical-iterative.hh>
#include <vector>

namespace hpp {
namespace constraints {
namespace solver {
/// \addtogroup solvers
/// \{

/// Solve a non-linear system equations with explicit and implicit constraints
///
/// This solver is defined in paper
/// https://hal.archives-ouvertes.fr/hal-01804774/file/paper.pdf. We
/// give here only a brief description
///
/// The unknows (denoted by \f$\mathbf{q}\f$) of the system of equations
/// is a Lie group. It is usually a robot configuration space or
/// the Cartesian product of robot configuration spaces.
///
/// The solver stores a set of implicit numerical constraints:
/// \f$g_1 (\mathbf{q}) = 0, g_2 (\mathbf{q}) = 0, \cdots\f$. These implicit
/// constraints are added using method HierarchicalIterative::add.
///
/// The solver also stores explicit numerical constraints (constraints where
/// some configuration variables depend on others) in an instance of class
/// ExplicitConstraintSet. This instance is accessible via method
/// BySubstitution::explicitConstraintSet.
///
/// When an explicit constraint is added using method
/// ExplicitConstraintSet::add, this method checks that the explicit
/// constraint is compatible with the previously added ones. If so,
/// the constraint is stored in the explicit constraint set. Otherwise,
/// it has to be added as an implicit constraint.
///
/// See Section III of the above mentioned paper for the description of
/// the constraint resolution.
class HPP_CONSTRAINTS_DLLAPI BySubstitution
    : public solver::HierarchicalIterative {
 public:
  BySubstitution(const LiegroupSpacePtr_t& configSpace);
  BySubstitution(const BySubstitution& other);

  virtual ~BySubstitution() {}

  /// \copydoc HierarchicalIterative::add
  ///
  /// If the constraint is explicit and compatible with previously
  /// inserted constraints, it is added as explicit. Otherwise, it is
  /// added as implicit.
  bool add(const ImplicitPtr_t& numericalConstraint,
           const std::size_t& priority = 0);

  /// \deprecated use add(const ImplicitPtr_t&, const std::size_t)
  HPP_CONSTRAINTS_DEPRECATED bool add(const ImplicitPtr_t& numericalConstraint,
                                      const segments_t& passiveDofs,
                                      const std::size_t priority = 0) {
    if (passiveDofs.size() > 0)
      throw std::invalid_argument(
          "Passive dof in the solver are not "
          "supported anymore. You must build an "
          "ActiveSetDifferentiableFunction yourself.");
    return add(numericalConstraint, priority);
  }

  /// Get the numerical constraints implicit and explicit
  const NumericalConstraints_t& numericalConstraints() const {
    return constraints_;
  }

  /// Get explicit constraint set
  ExplicitConstraintSet& explicitConstraintSet() { return explicit_; }

  /// Set explicit constraint set
  const ExplicitConstraintSet& explicitConstraintSet() const {
    return explicit_;
  }

  /// Return the number of free variables
  size_type numberFreeVariables() const {
    return explicitConstraintSet().notOutDers().nbIndices();
  }

  /// Should be called whenever explicit solver is modified
  void explicitConstraintSetHasChanged();

  template <typename LineSearchType>
  Status solve(vectorOut_t arg, LineSearchType ls = LineSearchType()) const {
    return solve<LineSearchType>(arg, false, ls);
  }

  template <typename LineSearchType>
  Status solve(vectorOut_t arg, bool optimize,
               LineSearchType ls = LineSearchType()) const {
    // TODO when there are only locked joint explicit constraints,
    // there is no need for this intricated loop.
    // if (explicit_.isConstant()) {
    // explicit_.solve(arg);
    // iterative_.solve(arg, ls);
    // } else {
    return impl_solve(arg, optimize, ls);
    // }
  }

  /// Project velocity on constraint tangent space in "from"
  ///
  /// \param from configuration,
  /// \param velocity velocity to project
  ///
  /// \f[
  /// \textbf{q}_{res} = \left(I_n -
  /// J^{+}J(\textbf{q}_{from})\right) (\textbf{v})
  /// \f]
  void projectVectorOnKernel(ConfigurationIn_t from, vectorIn_t velocity,
                             vectorOut_t result) const;

  /// Project configuration "to" on constraint tangent space in "from"
  ///
  /// \param from configuration,
  /// \param to configuration to project
  ///
  /// \f[
  /// \textbf{q}_{res} = \textbf{q}_{from} + \left(I_n -
  /// J^{+}J(\textbf{q}_{from})\right) (\textbf{q}_{to} -
  ///                                   \textbf{q}_{from})
  /// \f]
  virtual void projectOnKernel(ConfigurationIn_t from, ConfigurationIn_t to,
                               ConfigurationOut_t result);

  inline Status solve(vectorOut_t arg) const {
    return solve(arg, DefaultLineSearch());
  }

  /// \name Right hand side accessors
  /// \{

  /// \copydoc HierarchicalIterative::rightHandSideFromConfig(ConfigurationIn_t)
  vector_t rightHandSideFromConfig(ConfigurationIn_t config);

  /// \copydoc HierarchicalIterative::rightHandSideFromConfig(const
  /// ImplicitPtr_t&, ConfigurationIn_t)
  bool rightHandSideFromConfig(const ImplicitPtr_t& constraint,
                               ConfigurationIn_t config);

  /// \copydoc HierarchicalIterative::rightHandSide(const ImplicitPtr_t&,
  /// vectorIn_t)
  bool rightHandSide(const ImplicitPtr_t& constraint, vectorIn_t rhs);

  /// \copydoc HierarchicalIterative::getRightHandSide(const ImplicitPtr_t&,
  /// vectorOut_t)
  bool getRightHandSide(const ImplicitPtr_t& constraint, vectorOut_t rhs) const;

  /// \copydoc HierarchicalIterative::rightHandSide(vectorIn_t)
  void rightHandSide(vectorIn_t rhs);

  /// Get the right hand side
  /// \return the right hand side
  /// \note size of result is equal to total dimension of parameterizable
  ///       constraints (type Equality).
  vector_t rightHandSide() const;

  /// \copydoc HierarchicalIterative::rightHandSideSize()
  size_type rightHandSideSize() const;

  /// \}

  /// Return the size of the error as computed by isSatisfied
  ///
  /// Size of error is different from size of right hand side
  /// (rightHandSideSize) if the solver contains constraints with
  /// values in a Lie group.
  size_type errorSize() const {
    return dimension() + explicit_.outDers().nbIndices();
  }

  /// Whether input vector satisfies the constraints of the solver
  /// \param arg input vector.
  /// Compares to internal error threshold.
  bool isSatisfied(vectorIn_t arg) const {
    return solver::HierarchicalIterative::isSatisfied(arg) &&
           explicit_.isSatisfied(arg);
  }

  /// Whether input vector satisfies the constraints of the solver
  /// \param arg input vector
  /// \param errorThreshold threshold to use instead of the value
  ///        stored in the solver.
  bool isSatisfied(vectorIn_t arg, value_type errorThreshold) const {
    return solver::HierarchicalIterative::isSatisfied(arg, errorThreshold) &&
           explicit_.isSatisfied(arg, errorThreshold);
  }
  /// Whether input vector satisfies the constraints of the solver
  /// \param arg input vector
  /// \retval error the constraint errors dispatched in a vector,
  ///         the head of the vector corresponds to implicit constraints,
  ///         the tail of the vector corresponds to explicit constraints.
  bool isSatisfied(vectorIn_t arg, vectorOut_t error) const {
    assert(error.size() == dimension() + explicit_.errorSize());
    bool iterative = solver::HierarchicalIterative::isSatisfied(arg);
    residualError(error.head(dimension()));
    bool _explicit =
        explicit_.isSatisfied(arg, error.tail(explicit_.errorSize()));
    return iterative && _explicit;
  }

  /// Whether a constraint is satisfied for an input vector
  ///
  /// \param constraint, the constraint in the solver,
  /// \param arg the input vector
  /// \retval error the error of the constraint.
  /// \retval constraintFound whether the constraint belongs to the
  ///         solver,
  /// \return true if constraint belongs to solver and error is below
  ///         the threshold, false otherwise.
  bool isConstraintSatisfied(const ImplicitPtr_t& constraint, vectorIn_t arg,
                             vectorOut_t error, bool& constraintFound) const;

  template <typename LineSearchType>
  bool oneStep(vectorOut_t arg, LineSearchType& lineSearch) const {
    computeValue<true>(arg);
    updateJacobian(arg);
    computeDescentDirection();
    lineSearch(*this, arg, dq_);
    explicit_.solve(arg);
    return solver::HierarchicalIterative::isSatisfied(arg);
  }

  /// Computes the jacobian of the explicit functions and
  /// updates the jacobian of the problem using the chain rule.
  void updateJacobian(vectorIn_t arg) const;

  /// Set error threshold
  void errorThreshold(const value_type& threshold) {
    solver::HierarchicalIterative::errorThreshold(threshold);
    explicit_.errorThreshold(threshold);
  }
  /// Get error threshold
  value_type errorThreshold() const {
    return solver::HierarchicalIterative::errorThreshold();
  }

  /// Return the indices in the input vector which are solved implicitely.
  ///
  /// The other dof which are modified are solved explicitely.
  segments_t implicitDof() const;

  virtual std::ostream& print(std::ostream& os) const;

  bool integrate(vectorIn_t from, vectorIn_t velocity,
                 vectorOut_t result) const {
    bool res = solver::HierarchicalIterative::integrate(from, velocity, result);
    explicit_.solve(result);
    return res;
  }

 protected:
  void computeActiveRowsOfJ(std::size_t iStack);

 private:
  typedef solver::HierarchicalIterative parent_t;

  template <typename LineSearchType>
  Status impl_solve(vectorOut_t arg, bool optimize, LineSearchType ls) const;

  ExplicitConstraintSet explicit_;
  mutable matrix_t Je_, JeExpanded_;

  BySubstitution() {}
  HPP_SERIALIZABLE_SPLIT();
};  // class BySubstitution
/// \}

}  // namespace solver
}  // namespace constraints
}  // namespace hpp

#endif  // HPP_CONSTRAINTS_SOLVER_BY_SUBSTITUTION_HH
