// Copyright (c) 2015, LAAS-CNRS
// Authors: Florent Lamiraux
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

#ifndef HPP_CONSTRAINTS_EXPLICIT_HH
# define HPP_CONSTRAINTS_EXPLICIT_HH

# include <hpp/constraints/implicit.hh>

namespace hpp {
  namespace constraints {
    /// \addtogroup constraints
    /// \{

    /** Explicit numerical constraint

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
        \f{equation*}{
        \left(\begin{array}{c}
        q_{oc_{1}} \\ \vdots \\ q_{oc_{n_{oc}}}
        \end{array}\right) -
        f \left((q_{ic_{1}} \cdots q_{ic_{n_{ic}}})^T\right)
        \f}
        It is straightforward that an equality constraint with this function can
        solved explicitely:
        \f{align*}{
        \left(\begin{array}{c}
        q_{oc_{1}} \\ \vdots \\ q_{oc_{n_{oc}}}
        \end{array}\right) &-
        f \left((q_{ic_{1}} \cdots q_{ic_{n_{ic}}})^T\right) = rhs \\
        & \mbox {if and only if}\\
        \left(\begin{array}{c}
        q_{oc_{1}} \\ \vdots \\ q_{oc_{n_{oc}}}
        \end{array}\right) &=
        f \left((q_{ic_{1}} \cdots q_{ic_{n_{ic}}})^T\right) + rhs \\
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

        Considered as a Implicit instance, the expression of the Jacobian of
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
        \cdots & ov_1 & \cdots & iv_{1} & \cdots & ov_2 & \cdots & iv_2 & \cdots & ov_{n_{ov}} & \cdots \\
               &  1   &        &        &        &  0   &        &      &        &             &        \\
               &  0   &        &        &        &  1   &        &      &        &             &        \\
               &       &        &  -\frac{\partial f}{q_1} & & &   & -\frac{\partial f}{q_2} \\
          &&&&&\\
               & 0    &        &       &         &  0   &        &      &        &  1
        \end{array}\right)
        \f}
        The rows corresponding to values in SO(3) have the following expression.
        \f{equation*}{
        J = \left(\begin{array}{cccccccccccc}
        ov_1 \ ov_2 \ ov_3 & iv_1 \cdots  iv_{n_{iv}} \\
        J_{log}(R_{f}^T R_{out}) & -J_{log}(R_{f}^T R_{out})R_{out}^T R_{f} \frac{\partial f}{\partial q_{in}}
        \end{array}\right)
        \f}
        where
        \li \f$R_{out}\f$ is the rotation matrix corresponding to unit
        quaternion \f$(q_{oc1},q_{oc2},q_{oc3},q_{oc4})\f$,
        \li \f$R_{f}\f$ is the rotation matrix corresponding to the part of the
        output value of \f$f\f$ corresponding to SO(3),
        \li \f$J_{log}\f$ is the Jacobian matrix of function that associates
        to a rotation matrix \f$R\f$ the vector \f$\omega\f$ such that
        \f{equation*}{
        R = \exp (\left[\omega\right]_{\times})
        \f}



    **/
    class HPP_CONSTRAINTS_DLLAPI Explicit : public virtual Implicit
    {
    public:
      /// Copy object and return shared pointer to copy
      virtual ImplicitPtr_t copy () const;

      /// Create instance and return shared pointer
      ///
      /// \param robot Robot for which the constraint is defined.
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
      /// \note comparison type for this constraint is always equality
      /// \deprecated Call method that takes LiegroupSpacePtr_t instead of
      ///             DevicePtr_t as input and used robot->configSpace () as
      ///             argument.
      static ExplicitPtr_t create
        (const DevicePtr_t& robot,
         const DifferentiableFunctionPtr_t& function,
	 const segments_t& inputConf,
	 const segments_t& outputConf,
	 const segments_t& inputVelocity,
	 const segments_t& outputVelocity,
         const ComparisonTypes_t& comp = ComparisonTypes_t())
        HPP_CONSTRAINTS_DEPRECATED;

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
      /// \note comparison type for this constraint is always equality
      static ExplicitPtr_t create
        (const LiegroupSpacePtr_t& configSpace,
         const DifferentiableFunctionPtr_t& function,
	 const segments_t& inputConf,
	 const segments_t& outputConf,
	 const segments_t& inputVelocity,
	 const segments_t& outputVelocity,
         const ComparisonTypes_t& comp = ComparisonTypes_t());

      /// Create a copy and return shared pointer
      static ExplicitPtr_t createCopy (const ExplicitPtr_t& other);

      /// Function that maps input to output
      /// \return function \f$f\f$.
      virtual DifferentiableFunctionPtr_t explicitFunction() const
      {
        return inputToOutput_;
      }

      /// Get output configuration variables
      const segments_t& outputConf () const
      {
	return outputConf_;
      }
      /// Get output degrees of freedom
      const segments_t& outputVelocity () const
      {
	return outputVelocity_;
      }
      /// Get input configuration variables
      const segments_t& inputConf () const
      {
	return inputConf_;
      }
      /// Get input degrees of freedom
      const segments_t& inputVelocity () const
      {
	return inputVelocity_;
      }
      /// Convert right hand side
      ///
      /// \param implicitRhs right hand side of implicit formulation,
      /// \retval explicitRhs right hand side of explicit formulation.
      ///
      /// When implicit formulation is different from explicit formulation,
      ///\sa hpp::constraints::explicit_::RelativePose, right hand side are
      /// also different. This method converts right hand side from implicit
      /// to explicit formulation.
      ///
      /// When implicit formulation derives from explicit one, this method
      /// copies the first argument to the second one.
      virtual void implicitToExplicitRhs (vectorIn_t implicitRhs,
                                          vectorOut_t explicitRhs)
      {
        explicitRhs = implicitRhs;
      }

    protected:
      /// Constructor
      ///
      /// \param robot Robot for which the constraint is defined.
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
      /// \note comparison type for this constraint is always equality
      /// \deprecated Use constructor that takes LiegroupSpacePtr_t instead of
      ///             DevicePtr_t as input and used robot->configSpace () as
      ///             argument.
      Explicit
	(const DevicePtr_t& robot, const DifferentiableFunctionPtr_t& function,
	 const segments_t& inputConf,
	 const segments_t& outputConf,
	 const segments_t& inputVelocity,
	 const segments_t& outputVelocity,
         const ComparisonTypes_t& comp) HPP_CONSTRAINTS_DEPRECATED;

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
      /// \note comparison type for this constraint is always equality
      Explicit
	(const LiegroupSpacePtr_t& configSpace,
         const DifferentiableFunctionPtr_t& function,
	 const segments_t& inputConf,
	 const segments_t& outputConf,
	 const segments_t& inputVelocity,
	 const segments_t& outputVelocity,
         const ComparisonTypes_t& comp);

      /// Copy constructor
      Explicit (const Explicit& other);

      // Store weak pointer to itself
      void init (const ExplicitWkPtr_t& weak);
    protected:
      // Relation between input and output configuration variables
      DifferentiableFunctionPtr_t inputToOutput_;
      segments_t inputConf_;
      segments_t outputConf_;
      segments_t inputVelocity_;
      segments_t outputVelocity_;
    private:
      ExplicitWkPtr_t weak_;
    }; // class Explicit
    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_EXPLICIT_HH
