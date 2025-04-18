// Copyright (c) 2015, LAAS-CNRS
// Authors: Florent Lamiraux
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

#ifndef HPP_CONSTRAINTS_EXPLICIT_HH
#define HPP_CONSTRAINTS_EXPLICIT_HH

#include <exception>
#include <hpp/constraints/implicit.hh>

namespace hpp {
namespace constraints {
/// \addtogroup constraints
/// \{

/// Exception thrown when a function is evaluated outside its definition domain
class FunctionNotDefinedForThisValue : public std::exception {};

/** Explicit numerical constraint

    \paragraph hpp_constraints_explicit_definition Definition

    An explicit numerical constraint is a constraint such that some
    configuration variables called \b output are function of the
    others called \b input.

    Let
     \li \f$(ic_{1}, \cdots, ic_{n_{ic}})\f$ be the list of indices
         corresponding to ordered input configuration variables,
     \li \f$(oc_{1}, \cdots, oc_{n_{oc}})\f$ be the list of indices
         corresponding to ordered output configuration variables,
     \li \f$(iv_{1}, \cdots, iv_{n_{iv}})\f$ be the list of indices
         corresponding to ordered input degrees of freedom,
     \li \f$(ov_{1}, \cdots, ov_{n_{ov}})\f$ be the list of indices
         corresponding to ordered output degrees of freedom.

    Recall that degrees of freedom refer to velocity vectors.

    Let us notice that \f$n_{ic} + n_{oc}\f$ is less than the robot
    configuration size, and \f$n_{iv} + n_{ov}\f$ is less than the velocity
    size. Some degrees of freedom may indeed be neither input nor output.

    Then the differential function is of the form
    \f{eqnarray*}{
    \mathbf{q}_{out} - g \left(\mathbf{q}_{in}\right)
    \ \ &\mbox{with}&
    \mathbf{q}_{out} = \left(q_{oc_{1}} \cdots q_{oc_{n_{oc}}}\right)^T,
    \ \ \ \mathbf{q}_{in} = (q_{ic_{1}} \cdots q_{ic_{n_{ic}}})^T
    \f}
    It is straightforward that an equality constraint with this function can be
    solved explicitely:
    \f{align*}{
    \mathbf{q}_{out} &- g \left(\mathbf{q}_{in}\right) = rhs \\
    & \mbox {if and only if}\\
    \mathbf{q}_{out} &= g \left(\mathbf{q}_{in}\right) + rhs \\
    \f}

    If function \f$f\f$ takes values in a Lie group (SO(2), SO(3)),
    the above "+" between a Lie group element and a tangent vector
    has to be undestood as the integration of the constant velocity from
    the Lie group element:
    \f{equation*}{
    \mathbf{q} + \mathbf{v} = \mathbf{q}.\exp (\mathbf{v})
    \f}
    where \f$\mathbf{q}\f$ is a Lie group element and \f$\mathbf{v}\f$ is a
    tangent vector.

    Considered as an Implicit instance, the expression of the Jacobian of
    the DifferentiableFunction above depends on the output space of function
    \f$f\f$. The rows corresponding to values in a vector space are
    expressed as follows.

    for any index \f$i\f$ between 0 and the size of velocity vectors, either
    \li \f$\dot{q}_i\f$ is an input degree of freedom:
    \f$\exists j\f$ integer, \f$1 \leq j \leq n_{iv}\f$ such that
    \f$i=iv_{j}\f$,
    \li \f$\dot{q}_i\f$ is an output degree of freedom:
    \f$\exists j\f$ integer, \f$1\leq j \leq n_{ov}\f$ such that
    \f$i=ov_{j}\f$, or
    \li \f$\dot{q}_i\f$ neither input nor output. In this case, the
    corresponding column is equal to 0.
    \f{equation*}{
    J = \left(\begin{array}{cccccccccccc}
    \cdots & ov_1 & \cdots & iv_{1} & \cdots & ov_2 & \cdots & iv_2 & \cdots &
ov_{n_{ov}} & \cdots \\
           &  1   &        &        &        &  0   &        &      &        &
&        \\
           &  0   &        &        &        &  1   &        &      &        &
&        \\
           &       &        &  -\frac{\partial g}{q_1} & & &   & -\frac{\partial
g}{q_2} \\
      &&&&&\\
           & 0    &        &       &         &  0   &        &      &        & 1
    \end{array}\right)
    \f}
    The rows corresponding to values in SO(3) have the following expression.
    \f{equation*}{
    J = \left(\begin{array}{cccccccccccc}
    ov_1 \ ov_2 \ ov_3 & iv_1 \cdots  iv_{n_{iv}} \\
    J_{log}(R_{g}^T R_{out}) & -J_{log}(R_{g}^T R_{out})R_{out}^T R_{g}
\frac{\partial g}{\partial q_{in}} \end{array}\right) \f} where \li
\f$R_{out}\f$ is the rotation matrix corresponding to unit quaternion
\f$(q_{oc1},q_{oc2},q_{oc3},q_{oc4})\f$, \li \f$R_{g}\f$ is the rotation matrix
corresponding to the part of the output value of \f$f\f$ corresponding to SO(3),
    \li \f$J_{log}\f$ is the Jacobian matrix of function that associates
    to a rotation matrix \f$R\f$ the vector \f$\omega\f$ such that
    \f{equation*}{
    R = \exp (\left[\omega\right]_{\times})
    \f}

    \paragraph hpp_constraints_explicit_domain_of_definition "Domain of
definition"

    Some explicit constraints might be defined over only a subspace of the input
space. If the input value is not in the definition subspace, the explicit
constraint will throw an exception of type FunctionNotDefinedForThisValue.

**/
class HPP_CONSTRAINTS_DLLAPI Explicit : public Implicit {
 public:
  /// Copy object and return shared pointer to copy
  virtual ImplicitPtr_t copy() const;

  /// Create instance and return shared pointer
  ///
  /// \param configSpace Configuration space on which the constraint is
  ///        defined,
  /// \param function relation between input configuration variables and
  ///        output configuration variables,
  /// \param inputConf set of integer intervals defining indices
  ///            \f$(ic_{1}, \cdots, ic_{n_{ic}})\f$,
  /// \param outputConf set of integer intervals defining indices
  ///            \f$(oc_{1}, \cdots, oc_{n_{oc}})\f$,
  /// \param inputVelocity set of integer defining indices
  ///            \f$(iv_{1}, \cdots, iv_{n_{iv}})\f$.
  /// \param outputVelocity set of integer defining indices
  ///            \f$(ov_{1}, \cdots, ov_{n_{ov}})\f$.
  /// \param mask mask defining which components of the error are
  ///        taken into account to determine whether the constraint
  ///        is satisfied (See parent class for details).
  static ExplicitPtr_t create(
      const LiegroupSpacePtr_t& configSpace,
      const DifferentiableFunctionPtr_t& function, const segments_t& inputConf,
      const segments_t& outputConf, const segments_t& inputVelocity,
      const segments_t& outputVelocity,
      const ComparisonTypes_t& comp = ComparisonTypes_t(),
      std::vector<bool> mask = std::vector<bool>());

  /// Create a copy and return shared pointer
  static ExplicitPtr_t createCopy(const ExplicitPtr_t& other);

  /// Function that maps input to output
  /// \return function \f$f\f$.
  DifferentiableFunctionPtr_t explicitFunction() const {
    return inputToOutput_;
  }

  /// Get output configuration variables
  const segments_t& outputConf() const { return outputConf_; }
  /// Get output degrees of freedom
  const segments_t& outputVelocity() const { return outputVelocity_; }
  /// Get input configuration variables
  const segments_t& inputConf() const { return inputConf_; }
  /// Get input degrees of freedom
  const segments_t& inputVelocity() const { return inputVelocity_; }
  /// Compute the value of the output configuration variables
  /// \param qin input configuration variables,
  /// \param rhs right hand side of constraint
  ///
  /// The default implementation computes
  /// \f{equation}
  /// g \left((q_{ic_{1}} \cdots q_{ic_{n_{ic}}})^T\right) + rhs
  /// \f}
  virtual void outputValue(LiegroupElementRef result, vectorIn_t qin,
                           LiegroupElementConstRef rhs) const;

  /// Compute Jacobian of output value
  ///
  /// \f{eqnarray*}
  /// J &=& \frac{\partial}{\partial\mathbf{q}_{in}}\left(g(\mathbf{q}_{in})
  ///       + rhs\right).
  /// \f}
  /// \param qin vector of input variables,
  /// \param g_value \f$f(\mathbf{q}_{in})\f$ provided to avoid
  ///        recomputation,
  /// \param rhs right hand side (of implicit formulation).
  virtual void jacobianOutputValue(vectorIn_t qin,
                                   LiegroupElementConstRef g_value,
                                   LiegroupElementConstRef rhs,
                                   matrixOut_t jacobian) const;

 protected:
  /// Constructor
  ///
  /// \param configSpace Configuration space on which the constraint is
  ///        defined,
  /// \param function relation between input configuration variables and
  ///        output configuration variables,
  /// \param inputConf set of integer intervals defining indices
  ///            \f$(ic_{1}, \cdots, ic_{n_{ic}})\f$,
  /// \param outputConf set of integer intervals defining indices
  ///            \f$(oc_{1}, \cdots, oc_{n_{oc}})\f$,
  /// \param inputVelocity set of integer defining indices
  ///            \f$(iv_{1}, \cdots, iv_{n_{iv}})\f$.
  /// \param outputVelocity set of integer defining indices
  ///            \f$(ov_{1}, \cdots, ov_{n_{ov}})\f$.
  /// \param mask mask defining which components of the error are
  ///        taken into account to determine whether the constraint
  ///        is satisfied (See parent class for details).
  Explicit(const LiegroupSpacePtr_t& configSpace,
           const DifferentiableFunctionPtr_t& function,
           const segments_t& inputConf, const segments_t& outputConf,
           const segments_t& inputVelocity, const segments_t& outputVelocity,
           const ComparisonTypes_t& comp, std::vector<bool> mask);

  /// \copydoc Explicit (const LiegroupSpacePtr_t&, const
  /// DifferentiableFunctionPtr_t&, const segments_t& inputConf, const
  /// segments_t& outputConf, const segments_t& inputVelocity, const
  /// segments_t&, const ComparisonTypes_t&); \param implicitFunction
  /// differentiable function of the implicit
  ///        formulation if different from default expression
  ///        \f{equation}{
  ///        \mathbf{q}_{out} - g \left(\mathbf{q}_{in}\right),
  ///        \f}
  Explicit(const DifferentiableFunctionPtr_t& implicitFunction,
           const DifferentiableFunctionPtr_t& function,
           const segments_t& inputConf, const segments_t& outputConf,
           const segments_t& inputVelocity, const segments_t& outputVelocity,
           const ComparisonTypes_t& comp, std::vector<bool> mask);

  /// Copy constructor
  Explicit(const Explicit& other);

  bool isEqual(const Implicit& other, bool swapAndTest) const;

  // Store weak pointer to itself
  void init(const ExplicitWkPtr_t& weak);

 protected:
  // Relation between input and output configuration variables
  DifferentiableFunctionPtr_t inputToOutput_;
  segments_t inputConf_;
  segments_t outputConf_;
  segments_t inputVelocity_;
  segments_t outputVelocity_;

  Explicit() {}

 private:
  ExplicitWkPtr_t weak_;

  HPP_SERIALIZABLE();
};  // class Explicit
/// \}
}  // namespace constraints
}  // namespace hpp

BOOST_CLASS_EXPORT_KEY(hpp::constraints::Explicit)

#endif  // HPP_CONSTRAINTS_EXPLICIT_HH
