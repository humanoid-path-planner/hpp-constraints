// Copyright (c) 2017, 2018
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

#ifndef HPP_CONSTRAINTS_SOLVER_HIERARCHICAL_ITERATIVE_HH
#define HPP_CONSTRAINTS_SOLVER_HIERARCHICAL_ITERATIVE_HH

#include <functional>
#include <hpp/constraints/config.hh>
#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/implicit-constraint-set.hh>
#include <hpp/constraints/matrix-view.hh>
#include <hpp/util/serialization-fwd.hh>
#include <map>

namespace hpp {
namespace constraints {
namespace solver {
/// \addtogroup solvers
/// \{
namespace lineSearch {
/// No line search. Use \f$\alpha_i = 1\f$
struct Constant {
  template <typename SolverType>
  bool operator()(const SolverType& solver, vectorOut_t arg, vectorOut_t darg);
};

/// Implements the backtracking line search algorithm.
/// See https://en.wikipedia.org/wiki/Backtracking_line_search.
struct Backtracking {
  Backtracking();

  template <typename SolverType>
  bool operator()(const SolverType& solver, vectorOut_t arg, vectorOut_t darg);

  template <typename SolverType>
  inline value_type computeLocalSlope(const SolverType& solver) const;

  value_type c, tau, smallAlpha;  // 0.8 ^ 7 = 0.209, 0.8 ^ 8 = 0.1677
  mutable vector_t arg_darg, df, darg;
};

/// The step size is computed using the recursion
/// \f$ \alpha_{i+1} = \alpha - K \times (\alpha_{max} - \alpha_i) \f$
/// where \f$K\f$ and \f$\alpha_{max}\f$ are some constant values.
struct FixedSequence {
  FixedSequence();

  template <typename SolverType>
  bool operator()(const SolverType& solver, vectorOut_t arg, vectorOut_t darg);

  value_type alpha;
  value_type alphaMax, K;
};

/// The step size is computed using the formula
/// \f$ \alpha_{i} = C - K \times \text{tanh}(a
/// \frac{\|f(\mathbf{q}_i)\|}{\epsilon^2} + b) \f$, where \li \f$\epsilon\f$ is
/// the error threshold: if \f$\|f(\mathbf{q}_i)\|<\epsilon\f$,
/// \f$\mathbf{q}_i\f$ is considered to satisfy the constraint.
struct ErrorNormBased {
  ErrorNormBased(value_type alphaMin, value_type _a, value_type _b);
  ErrorNormBased(value_type alphaMin = 0.2);

  template <typename SolverType>
  bool operator()(const SolverType& solver, vectorOut_t arg, vectorOut_t darg);

  value_type C, K, a, b;
};
}  // namespace lineSearch

namespace saturation {
/// \brief Base class for box constraints.
/// To prevent configuration variables to get out of joint limits during
/// Newton Raphson iterations, the user may provide an object derived
/// from this type to HierarchicalIterative using setter and getter
/// HierarchicalIterative::saturation.
///
/// This class implements no saturation.
struct Base {
  /// This function checks which degrees of freedom are saturated.
  ///
  /// \param q a configuration,
  /// \retval qSat configuration after saturing values out of bounds
  /// \retval saturation vector: for each degree of freedom, saturation
  ///         is set to
  /// \li -1 if the lower bound is reached,
  /// \li  1 if the upper bound is reached,
  /// \li  0 otherwise.
  /// \return true if and only if at least one degree of freedom has been
  ///         saturated
  virtual bool saturate(vectorIn_t q, vectorOut_t qSat,
                        Eigen::VectorXi& saturation);
  virtual ~Base() {}
};
/// \brief saturation from a std::function.
struct Function : Base {
  typedef std::function<bool(vectorIn_t, vectorOut_t, Eigen::VectorXi&)>
      function_t;
  bool saturate(vectorIn_t q, vectorOut_t qSat, Eigen::VectorXi& saturation) {
    return function(q, qSat, saturation);
  }
  Function() {}
  Function(const function_t& function) : function(function) {}
  function_t function;
};
/// \brief simple box constraints
struct Bounds : Base {
  bool saturate(vectorIn_t q, vectorOut_t qSat, Eigen::VectorXi& saturation);
  Bounds() {}
  Bounds(const vector_t& lb, const vector_t& ub) : lb(lb), ub(ub), iq2iv_() {
    iq2iv_.resize(ub.size());
    for (size_type i = 0; i < ub.size(); ++i) {
      iq2iv_[i] = (Eigen::VectorXi::Scalar)i;
    }
  }
  Bounds(const vector_t& lb, const vector_t& ub, const Eigen::VectorXi& iq2iv)
      : lb(lb), ub(ub), iq2iv_(iq2iv) {}
  vector_t lb, ub;
  Eigen::VectorXi iq2iv_;
};
/// \brief Box constraints use a Device joint limits.
struct Device : Base {
  /// \todo write a visitor.
  bool saturate(vectorIn_t q, vectorOut_t qSat, Eigen::VectorXi& saturation);
  Device() {}
  Device(const DevicePtr_t& device) : device(device) {}
  DevicePtr_t device;
};
}  // namespace saturation

/// Solve a system of non-linear equations on a robot configuration
///
/// The non-linear system of equations is built by adding equations with
/// method HierarchicalIterative::add.
///
/// Note that a hierarchy between the equations can be
/// provided. In this case, the solver will try to solve the
/// highest priority equations first, and then to solve the lower priority
/// equations. Note that priorities are in decreasing order: 0 has higher
/// priority than 1.
///
/// The algorithm used is a Newton-Raphson like algorithm that works as
/// follows: for a single level of priority, let \f$f (\mathbf{q}) = 0\f$
/// be the system of equations where
/// \f$f\f$ is a \f$C^1\f$ mapping from the robot configuration space to
/// a Lie group space \f$\mathcal{L}\f$.
///
/// Starting from initial guess \f$\mathbf{q}_0\f$, the method
/// HierarchicalIterative::solve builds a sequence of configurations
/// \f$\mathbf{q}_i\f$ as follows:
/// \f{eqnarray*}
/// \mathbf{q}_{i+1} = \mathbf{q}_i -
///    \alpha_i \frac{\partial f}{\partial \mathbf{q}}(\mathbf{q}_i)^{+}
///    f (\mathbf{q}_i)
/// \f}
/// where
/// \li \f$\frac{\partial f}{\partial \mathbf{q}}(\mathbf{q}_i)^{+}\f$ is
///     the Moore-Penrose pseudo-inverse of the system Jacobian,
/// \li \f$\alpha_i\f$ is a sequence of real numbers depending on the
///     line search strategy. Possible line-search strategies are
///     lineSearch::Constant, lineSearch::Backtracking,
///     lineSearch::FixedSequence, lineSearch::ErrorNormBased.
/// until
/// \li the residual \f$\|f(\mathbf{q})\|\f$ is below an error threshold, or
/// \li the maximal number of iterations has been reached.
///
/// The computation of the direction of descent in the case of multiple
/// level of hierarchy is described in
/// <a href="https://hal.archives-ouvertes.fr/lirmm-00796736">this paper</a>.
///
/// The error threshold can be accessed by methods
/// HierarchicalIterative::errorThreshold. The maximal number of
/// iterations can be accessed by methods
/// HierarchicalIterative::maxIterations.
///
/// \note Solving equations one after the other
///
/// For some applications, it can be more efficient to solve a set of
/// equations one after the other. In other words, an equation is ignored
/// until the previous one is solved (norm below the threshold). To do so,
/// introduce the equations using method HierarchicalIterative::add with
/// increasing value of \c priority, and call method
/// \c HierarchicalIterative::solveLevelByLevel(true).
/// \note Lie group
///
/// The unknowns \f$\mathbf{q}\f$ may take values in a more general set
/// than the configuration space of a robot. This set should be a
/// Cartesian product of Lie groups: hpp::pinocchio::LiegroupSpace.
///
/// \note Saturation
///
/// \copydoc saturate::Base
///
/// \note Right hand side and comparison types
///
/// Instead of \f$f(\mathbf{q}) = 0\f$, other constraints can be defined.
/// Several comparison types are available:
/// \li Equality: \f$f(\mathbf{q}) = rhs\f$, where \f$rhs\f$ is a
/// parameterizable right hand side,
/// \li EqualToZero: \f$f(\mathbf{q}) = 0\f$,
/// \li Superior: \f$f(\mathbf{q}) > 0\f$
/// \li Inferior: \f$f(\mathbf{q}) < 0\f$
/// If several constraint are of type equality, the right hand side of the
/// system of equations can be modified by methods
/// HierarchicalIterative::rightHandSideFromInput,
/// HierarchicalIterative::rightHandSide.
///
/// \note Free variables
///
/// Some variables can be locked, or computed explicitely. In this case,
/// the iterative resolution will only change the other variables called
/// free variables. \sa methods
/// \li \ref freeVariables (const Indices_t& indices) and
/// \li \ref freeVariables (const Indices_t& indices).
class HPP_CONSTRAINTS_DLLAPI HierarchicalIterative {
 public:
  typedef Eigen::RowBlockIndices Indices_t;
  typedef lineSearch::FixedSequence DefaultLineSearch;

  enum Status { ERROR_INCREASED, MAX_ITERATION_REACHED, INFEASIBLE, SUCCESS };
  typedef shared_ptr<saturation::Base> Saturation_t;

  HierarchicalIterative(const LiegroupSpacePtr_t& configSpace);

  HierarchicalIterative(const HierarchicalIterative& other);

  virtual ~HierarchicalIterative() {}

  /// \name Problem definition
  /// \{

  /// Get configuration space on which constraints are defined
  const LiegroupSpacePtr_t& configSpace() const { return configSpace_; }
  /// Check whether a numerical constraint has been added
  /// \param numericalConstraint numerical constraint
  /// \return true if numerical constraint is already in the solver
  virtual bool contains(const ImplicitPtr_t& numericalConstraint) const;

  /// Add an implicit constraint
  ///
  /// \param constraint implicit constraint
  /// \param priority level of priority of the constraint: priority are
  ///        in decreasing order: 0 is the highest priority level,
  virtual bool add(const ImplicitPtr_t& constraint,
                   const std::size_t& priority);

  /// add constraints of another solver
  /// \param other other solver
  ///
  /// Add constraints of other to this solver.
  /// \note right hand side of other is not copied.
  virtual void merge(const HierarchicalIterative& other);

  /// Set the saturation function
  void saturation(const Saturation_t& saturate) { saturate_ = saturate; }

  /// Get the saturation function
  const Saturation_t& saturation() const { return saturate_; }

  /// \}

  /// \name Problem resolution
  /// \{

  /// Solve the system of non linear equations
  ///
  /// \param arg initial guess,
  /// \param ls line search method used.
  ///
  /// Use Newton Rhapson like iterative method until the error is below
  /// the threshold, or until the maximal number of iterations has been
  /// reached.
  ///
  /// \note Explicit constraints are expressed in their implicit
  ///       form: \f$\mathbf{q}_2 = f (\mathbf{q}_1)\f$ is replaced by
  ///       \f$\mathbf{q}_2 - f (\mathbf{q}_1) = 0\f$,
  ///  where \f$\mathbf{q}_1\f$ and \f$\mathbf{q}_2\f$ are vectors
  ///  composed of the components of \f$\mathbf{q}\f$.
  template <typename LineSearchType>
  Status solve(vectorOut_t arg, LineSearchType ls = LineSearchType()) const;

  /// Solve the system of non linear equations
  ///
  /// \param arg initial guess,
  ///
  /// Use Newton Rhapson like iterative method until the error is below
  /// the threshold, or until the maximal number of iterations has been
  /// reached. Use the default line search method (fixed sequence of
  /// \f$\alpha_i\f$).
  ///
  /// \note Explicit constraints are expressed in their implicit
  ///       form: \f$\mathbf{q}_2 = f (\mathbf{q}_1)\f$ is replaced by
  ///       \f$\mathbf{q}_2 - f (\mathbf{q}_1) = 0\f$,
  ///  where \f$\mathbf{q}_1\f$ and \f$\mathbf{q}_2\f$ are vectors
  ///  composed of the components of \f$\mathbf{q}\f$.
  inline Status solve(vectorOut_t arg) const {
    return solve(arg, DefaultLineSearch());
  }

  /// Whether input vector satisfies the constraints of the solver
  /// \param arg input vector
  /// Compares to internal error threshold.
  bool isSatisfied(vectorIn_t arg) const {
    computeValue<false>(arg);
    computeError();
    return squaredNorm_ < squaredErrorThreshold_;
  }

  /// Whether input vector satisfies the constraints of the solver
  /// \param arg input vector
  /// \param errorThreshold threshold to use instead of the value
  ///        stored in the solver.
  bool isSatisfied(vectorIn_t arg, value_type errorThreshold) const {
    computeValue<false>(arg);
    computeError();
    return squaredNorm_ < errorThreshold * errorThreshold;
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

  /// Returns the lowest singular value.
  /// If the jacobian has maximum rank r, then it corresponds to r-th
  /// greatest singular value. This value is zero when the jacobian is
  /// singular.
  const value_type& sigma() const { return sigma_; }

  /// \}

  /// \name Parameters
  /// \{

  /// Set free velocity variables
  ///
  /// The other variables will be left unchanged by the iterative
  /// resolution.
  ///
  /// \param intervals set of index intervals
  void freeVariables(const segments_t intervals) {
    freeVariables_ = Indices_t();
    for (std::size_t i = 0; i < intervals.size(); ++i)
      freeVariables_.addRow(intervals[i].first, intervals[i].second);
    freeVariables_.updateIndices<true, true, true>();
    update();
  }

  /// Set free velocity variables
  ///
  /// The other variables will be left unchanged by the iterative
  /// resolution.
  void freeVariables(const Indices_t& indices) {
    freeVariables_ = indices;
    update();
  }

  /// Get free velocity variables
  const Indices_t& freeVariables() const { return freeVariables_; }

  /// Set maximal number of iterations
  void maxIterations(size_type iterations) { maxIterations_ = iterations; }
  /// Get maximal number of iterations in config projector
  size_type maxIterations() const { return maxIterations_; }

  /// Set error threshold
  void errorThreshold(const value_type& threshold) {
    squaredErrorThreshold_ = threshold * threshold;
  }
  /// Get error threshold
  value_type errorThreshold() const { return sqrt(squaredErrorThreshold_); }
  /// Get error threshold
  value_type squaredErrorThreshold() const { return squaredErrorThreshold_; }

  /// Get the inequality threshold
  value_type inequalityThreshold() const { return inequalityThreshold_; }
  /// set the inequality threshold
  void inequalityThreshold(const value_type& it) { inequalityThreshold_ = it; }

  void lastIsOptional(bool optional) { lastIsOptional_ = optional; }

  bool lastIsOptional() const { return lastIsOptional_; }

  // Set whether the constraints are solved level by level
  //
  // If true, the solver will consider only the higher level of
  // priority that is not satisfied when computing the direction
  // of descent.
  void solveLevelByLevel(bool solveLevelByLevel) {
    solveLevelByLevel_ = solveLevelByLevel;
  }

  // Get whether the constraints are solved level by level
  bool solveLevelByLevel() const { return solveLevelByLevel_; }
  /// \}

  /// \name Stack
  /// \{

  /// Get set of constraints for a give priority level.
  const ImplicitConstraintSet& constraints(const std::size_t priority) {
    assert(priority < stacks_.size());
    return stacks_[priority];
  }

  /// Get constraints (implicit and explicit)
  const NumericalConstraints_t& constraints() const { return constraints_; }

  std::size_t numberStacks() const { return stacks_.size(); }

  const size_type& dimension() const { return dimension_; }

  /// Dimension of the problem after removing the rows of the jacobian
  /// which do not influence the error (only zeros along those lines).
  const size_type& reducedDimension() const { return reducedDimension_; }

  /// Configuration parameters involved in the constraint resolution.
  ArrayXb activeParameters() const;

  /// Velocity parameters involved in the constraint resolution.
  ArrayXb activeDerivativeParameters() const;

  /// \}

  /// Returns the squared norm of the error vector
  value_type residualError() const { return squaredNorm_; }

  /// Returns the error vector
  void residualError(vectorOut_t error) const;

  /// Inclusion of manifolds
  ///
  /// \param solver another solver
  ///
  /// This function returns true if the solution manifold of the solver
  /// given as input is a submanifold of the solution manifold of this
  /// solver. The function tests whether constraints of the input solver
  /// are also in this solver.
  bool definesSubmanifoldOf(const HierarchicalIterative& solver) const;

  /// \name Right hand side accessors
  /// \{

  /// Compute right hand side of equality constraints from a configuration
  /// \param config a configuration.
  ///
  /// for each constraint of type Equality, set right hand side as
  /// \f$rhs = f(\mathbf{q})\f$.
  /// \note Only parameterizable constraints (type Equality) are set
  vector_t rightHandSideFromConfig(ConfigurationIn_t config);

  /// Compute right hand side of a constraint from a configuration
  /// \param constraint the constraint,
  /// \param config a configuration.
  ///
  /// Set right hand side as \f$rhs = f(\mathbf{q})\f$.
  /// \note Only parameterizable constraints (type Equality) are set
  virtual bool rightHandSideFromConfig(const ImplicitPtr_t& constraint,
                                       ConfigurationIn_t config);
  /// Set right hand side of a constraints
  /// \param constraint the constraint,
  /// \param rhs right hand side.
  /// \note Size of rhs should be equal to the total dimension of
  ///       parameterizable constraints (type Equality) .
  virtual bool rightHandSide(const ImplicitPtr_t& constraint, vectorIn_t rhs);

  /// Get right hand side of a constraints
  virtual bool getRightHandSide(const ImplicitPtr_t& constraint,
                                vectorOut_t rhs) const;

  /// Set the right hand side
  /// \param rhs the right hand side
  /// \note Size of rhs should be equal to the total dimension of
  ///       parameterizable constraints (type Equality).
  virtual void rightHandSide(vectorIn_t rhs);

  /// Set the right hand side at a given parameter.
  /// \param s parameter passed to Implicit::rightHandSideAt
  void rightHandSideAt(const value_type& s);

  /// Get the right hand side
  /// \return the right hand side
  /// \note size of result is equal to total dimension of parameterizable
  ///       constraints (type Equality).
  vector_t rightHandSide() const;

  /// Get size of the right hand side
  /// \return sum of dimensions of parameterizable constraints
  ///         (type Equality)
  size_type rightHandSideSize() const;

  /// \}

  /// \name Access to internal data
  /// You should know what you do when you call these functions
  /// \{

  /// Compute the value of each level, and the jacobian if ComputeJac is true.
  template <bool ComputeJac>
  void computeValue(vectorIn_t arg) const;
  void computeSaturation(vectorIn_t arg) const;
  void getValue(vectorOut_t v) const;
  void getReducedJacobian(matrixOut_t J) const;
  /// If lastIsOptional() is true, then the last level is ignored.
  /// \warning computeValue must have been called first.
  void computeError() const;

  /// Accessor to the last step done
  const vector_t& lastStep() const { return dq_; }

  virtual bool integrate(vectorIn_t from, vectorIn_t velocity,
                         vectorOut_t result) const;
  /// \}

  virtual std::ostream& print(std::ostream& os) const;

 protected:
  typedef Eigen::JacobiSVD<matrix_t> SVD_t;

  struct Data {
    /// \cond
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// \endcond
    LiegroupElement output, rightHandSide;
    vector_t error;
    matrix_t jacobian, reducedJ;

    SVD_t svd;
    matrix_t PK;

    mutable size_type maxRank;

    ComparisonTypes_t comparison;
    std::vector<std::size_t> inequalityIndices;
    Eigen::RowBlockIndices equalityIndices;
    Eigen::MatrixBlocks<false, false> activeRowsOfJ;
  };

  /// Allocate datas and update sizes of the problem
  /// Should be called whenever the stack is modified.
  void update();

  /// Compute which rows of the jacobian of stack_[iStack]
  /// are not zero, using the activeDerivativeParameters of the functions.
  /// The result is stored in datas_[i].activeRowsOfJ
  virtual void computeActiveRowsOfJ(std::size_t iStack);

  /// Compute a SVD decomposition of each level and find the best descent
  /// direction at the first order.
  /// Linearization of the system of equations
  /// rhs - v_{i} = J (q_i) (dq_{i+1} - q_{i})
  /// q_{i+1} - q_{i} = J(q_i)^{+} ( rhs - v_{i} )
  /// dq = J(q_i)^{+} ( rhs - v_{i} )
  /// \warning computeValue<true> must have been called first.
  void computeDescentDirection() const;
  void expandDqSmall() const;
  void saturate(vectorOut_t arg) const;

  value_type squaredErrorThreshold_, inequalityThreshold_;
  size_type maxIterations_;

  std::vector<ImplicitConstraintSet> stacks_;
  LiegroupSpacePtr_t configSpace_;
  size_type dimension_, reducedDimension_;
  bool lastIsOptional_;
  bool solveLevelByLevel_;
  /// Unknown of the set of implicit constraints
  Indices_t freeVariables_;
  Saturation_t saturate_;
  /// Members moved from core::ConfigProjector
  NumericalConstraints_t constraints_;
  /// Value rank of constraint in its priority level
  std::map<DifferentiableFunctionPtr_t, size_type> iq_;
  /// Derivative rank of constraint in its priority level
  std::map<DifferentiableFunctionPtr_t, size_type> iv_;
  /// Priority level of constraint
  std::map<DifferentiableFunctionPtr_t, std::size_t> priority_;

  /// The smallest non-zero singular value
  mutable value_type sigma_;

  mutable vector_t dq_, dqSmall_;
  mutable matrix_t reducedJ_;
  mutable Eigen::VectorXi saturation_, reducedSaturation_;
  mutable Configuration_t qSat_;
  mutable ArrayXb tmpSat_;
  mutable value_type squaredNorm_;
  mutable std::vector<Data> datas_;
  mutable SVD_t svd_;
  mutable vector_t OM_;
  mutable vector_t OP_;

  friend struct lineSearch::Backtracking;

 protected:
  HierarchicalIterative() {}

 private:
  HPP_SERIALIZABLE_SPLIT();
};  // class HierarchicalIterative
/// \}
inline std::ostream& operator<<(std::ostream& os,
                                const HierarchicalIterative& hs) {
  return hs.print(os);
}
}  // namespace solver
}  // namespace constraints
}  // namespace hpp

#endif  // HPP_CONSTRAINTS_SOLVER_HIERARCHICAL_ITERATIVE_HH
