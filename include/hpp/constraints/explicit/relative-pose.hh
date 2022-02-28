// Copyright (c) 2018, LAAS-CNRS
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

#ifndef HPP_CONSTRAINTS_EXPLICIT_RELATIVE_POSE_HH
# define HPP_CONSTRAINTS_EXPLICIT_RELATIVE_POSE_HH

# include <hpp/constraints/explicit.hh>
# include <pinocchio/spatial/se3.hpp>

namespace hpp {
  namespace constraints {
    namespace explicit_ {
      /// Constraint of relative pose between two frames on a kinematic chain
      ///
      /// Given an input configuration \f$\mathbf{q}\f$, solving this constraint
      /// consists in computing output variables with respect to input
      /// variables:
      /// \f[
      /// \mathbf{q}_{out} = g(\mathbf{q}_{in})
      ///\f]
      /// where \f$\mathbf{q}_{out}\f$ are the configuration variables of
      /// input joint2, \f${q}_{in}\f$ are the configuration variables
      /// of input joint1 and parent joints, and \f$g\f$ is a mapping
      /// with values is SE(3). Note that joint2 should be a
      /// freeflyer joint.
      ///
      /// \note According to the documentation of class Explicit, the
      /// implicit formulation should be
      /// \f[
      /// f(\mathbf{q}) = \mathbf{q}_{out} - g(\mathbf{q}_{in}).
      /// \f]
      /// As function \f$g\f$ takes values in SE(3), the above expression is
      /// equivalent to
      /// \f[
      /// f(\mathbf{q}) = \log_{SE(3)}\left(g(\mathbf{q}_{in})^{-1}
      ///                 \mathbf{q}_{out}\right)
      /// \f]
      /// that represents the log of the error of input joint2 pose
      /// (\f$\mathbf{q}_{out}\f$) with respect to its desired value
      /// (\f$g(\mathbf{q}_{in}\f$). The problem with this expression is
      /// that it is different from the corresponding implicit formulation
      /// hpp::constraints::RelativeTransformationR3xSO3 that compares the poses of input
      /// joint1 and joint2. For manipulation planning applications where pairs
      /// of constraints and complements are replaced by an explicit constraint,
      /// this difference of formulation results in inconsistencies such as a
      /// configuration satisfying one formulation (the error being below the
      /// threshold) but not the other one. To cope with this issue, the default
      /// implicit formulation is replaced by the one defined by class
      /// hpp::constraints::RelativeTransformationR3xSO3.
      class HPP_CONSTRAINTS_DLLAPI RelativePose : public Explicit
      {
      public:
        /// Copy object and return shared pointer to copy
        virtual ImplicitPtr_t copy () const;
        /// Create instance and return shared pointer
        ///
        /// \param name the name of the constraints,
        /// \param robot the robot the constraints is applied to,
        /// \param joint1 the first joint the transformation of which is
        ///               constrained,
        /// \param joint2 the second joint the transformation of which is
        ///               constrained,
        /// \param frame1 position of a fixed frame in joint 1,
        /// \param frame2 position of a fixed frame in joint 2,
        /// \param comp vector of comparison types
	/// \param mask mask defining which components of the error are
	///        taken into account to determine whether the constraint
	///        is satisfied.
        /// \note if joint1 is 0x0, joint 1 frame is considered to be the global
        ///       frame.
        static RelativePosePtr_t create
          (const std::string& name, const DevicePtr_t& robot,
           const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
           const Transform3f& frame1, const Transform3f& frame2,
           ComparisonTypes_t comp, std::vector<bool> mask=std::vector<bool>());

        static RelativePosePtr_t createCopy (const RelativePosePtr_t& other);

        /// Compute the value of the output configuration variables
        /// \param qin input configuration variables,
        /// \param rhs right hand side of constraint
        ///
        /// \f{equation}
        /// f \left(\mathbf{q}_{in}\right) + rhs_{expl}
        /// \f}
        /// where \f$rhs_{expl}\f$ is the explicit right hand side converted
        /// using the following expression:
        /// \f{equation}
        /// rhs_{expl} = \log_{SE(3)}\left( F_{2/J_2} rhs_{impl} F_{2/J_2}^{-1}
	///              \right)
        /// \f}
        virtual void outputValue(LiegroupElementRef result, vectorIn_t qin,
                                 LiegroupElementConstRef rhs) const;

        /// Compute Jacobian of output value
        ///
        /// \f{eqnarray*}
        /// J &=& \frac{\partial}{\partial\mathbf{q}_{in}}\left(
        ///       f(\mathbf{q}_{in}) + rhs\right).
        /// \f}
        /// \param qin vector of input variables,
        /// \param f_value \f$f(\mathbf{q}_{in})\f$ to avoid recomputation,
        /// \param rhs right hand side (of implicit formulation).
        virtual void jacobianOutputValue(vectorIn_t qin, LiegroupElementConstRef
                                         f_value, LiegroupElementConstRef rhs,
                                         matrixOut_t jacobian) const;
      protected:
        /// Constructor
        ///
        /// \param name the name of the constraints,
        /// \param robot the robot the constraints is applied to,
        /// \param joint1 the first joint the transformation of which is
        ///               constrained,
        /// \param joint2 the second joint the transformation of which is
        ///               constrained,
        /// \param frame1 position of a fixed frame in joint 1,
        /// \param frame2 position of a fixed frame in joint 2,
        /// \param comp vector of comparison types
	/// \param mask mask defining which components of the error are
	///        taken into account to determine whether the constraint
	///        is satisfied.
        /// \note if joint1 is 0x0, joint 1 frame is considered to be the global
        ///       frame.
        RelativePose (const std::string& name, const DevicePtr_t& robot,
                      const JointConstPtr_t& joint1,
                      const JointConstPtr_t& joint2,
                      const Transform3f& frame1, const Transform3f& frame2,
                      ComparisonTypes_t comp = ComparisonTypes_t(),
		      std::vector<bool> mask = std::vector<bool>(6, true));

        /// Copy constructor
        RelativePose (const RelativePose& other);

        /// Store weak pointer to itself
        void init (RelativePoseWkPtr_t weak);
      private:
        /** Convert right hand side

            \param implicitRhs right hand side of implicit formulation, this
	           is an element of $SE(3)$
            \retval explicitRhs right hand side of explicit formulation, this
                    is an element of $\mathbf{R}^6$.

            For this constraint, the implicit formulation does not derive
            from  the explicit formulation. The explicit form writes

            \f{eqnarray}
            rhs_{expl} &=& \log_{SE(3)} \left(F_{2/J_2} F_{1/J_1}^{-1} J_1^{-1}
            J_2\right)\\
            rhs_{impl} &=& F_{1/J_1}^{-1} J_1^{-1}J_2 F_{2/J_2}
            \f}
            Thus
            \f{equation}
            rhs_{expl} = \log_{SE(3)}\left( F_{2/J_2}rhs_{impl}
            F_{2/J_2}^{-1}\right)
            \f}
        */
        void implicitToExplicitRhs (LiegroupElementConstRef implicitRhs,
                                    vectorOut_t explicitRhs) const;
        /** Convert right hand side

            \param explicitRhs right hand side of explicit formulation,
            \retval implicitRhs right hand side of implicit formulation.

            For this constraint, the implicit formulation does not derive
            from  the explicit formulation. The explicit form writes

            \f{eqnarray}
            rhs_{expl} &=& \log_{SE(3)} \left(F_{2/J_2} F_{1/J_1}^{-1} J_1^{-1}
            J_2\right)\\
            rhs_{impl} &=& F_{1/J_1}^{-1} J_1^{-1}J_2 F_{2/J_2}
            \f}
            Thus
            \f{equation}
            rhs_{impl} = F_{2/J_2}^{-1} \exp_{SE(3)}(rhs_{expl}) F_{2/J_2}
            \f}
        */
        void explicitToImplicitRhs (vectorIn_t explicitRhs,
                                    LiegroupElementRef implicitRhs) const;
        // Create LiegroupSpace instances to avoid useless allocation.
        static LiegroupSpacePtr_t SE3;
        static LiegroupSpacePtr_t R3xSO3;
        JointConstPtr_t joint1_, joint2_;
        Transform3f frame1_;
        Transform3f frame2_;
        RelativePoseWkPtr_t weak_;

        RelativePose() {}
        HPP_SERIALIZABLE();
      }; // class RelativePose
    } // namespace explicit_
  } // namespace constraints
} // namespace hpp

BOOST_CLASS_EXPORT_KEY(hpp::constraints::explicit_::RelativePose)

#endif //HPP_CONSTRAINTS_EXPLICIT_RELATIVE_POSE_HH
