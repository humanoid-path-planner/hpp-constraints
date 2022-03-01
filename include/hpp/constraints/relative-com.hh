//
// Copyright (c) 2014 CNRS
// Authors: Florent Lamiraux
//
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

#ifndef HPP_CONSTRAINTS_RELATIVE_COM_HH
# define HPP_CONSTRAINTS_RELATIVE_COM_HH

# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/config.hh>
# include <hpp/constraints/fwd.hh>

namespace hpp {
  namespace constraints {

    /// \addtogroup constraints
    /// \{

    /**
     *  Constraint on the relative position of the center of mass
     *
     *  The value of the function is defined as the position of the center
     *  of mass in the reference frame of a joint.
     *
     *  \f{eqnarray*}
     *  \mathbf{f}(\mathbf{q}) &=& R^T \left(\mathbf{x} - \mathbf{t}\right)
     *  - \mathbf{x}^{*}\\
     *  \mathbf{\dot{f}} &=& R^T \left(
     *  J_{com} + [\mathbf{x}-\mathbf{t}]_{\times}J_{joint}^{\omega}
     *  - J_{joint}^{\mathbf{v}}\right)\mathbf{\dot{q}}
     *  \f}
     *
     *  where
     *  \li \f$
     *      \left(\begin{array}{cc} R & \mathbf{t} \\ 0 & 1\end{array}\right)
     *      \f$
     *  is the position of the joint,
     *  \li \f$\mathbf{x}\f$ is the position of the center of mass,
     *  \li \f$\mathbf{x}^{*}\f$ is the desired position of the center of mass
     *      expressed in joint frame.
    **/
    class HPP_CONSTRAINTS_DLLAPI RelativeCom : public DifferentiableFunction
    {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      /// Return a shared pointer to a new instance
      static RelativeComPtr_t create (const DevicePtr_t& robot,
				      const JointPtr_t& joint,
				      const vector3_t reference,
                                      std::vector <bool> mask = std::vector<bool>(3, true))
      { return create ("RelativeCom", robot, joint, reference, mask); }
      static RelativeComPtr_t create (const std::string& name,
                                      const DevicePtr_t& robot,
				      const JointPtr_t& joint,
				      const vector3_t reference,
                                      std::vector <bool> mask = std::vector<bool>(3, true));
      static RelativeComPtr_t create (const DevicePtr_t& robot,
                                      const CenterOfMassComputationPtr_t& comc,
				      const JointPtr_t& joint,
				      const vector3_t reference,
                                      std::vector <bool> mask = std::vector<bool>(3, true));
      static RelativeComPtr_t create (const std::string& name,
                                      const DevicePtr_t& robot,
                                      const CenterOfMassComputationPtr_t& comc,
				      const JointPtr_t& joint,
				      const vector3_t reference,
                                      std::vector <bool> mask = std::vector<bool> (3, true));
      virtual ~RelativeCom () {}
      RelativeCom (const DevicePtr_t& robot,
          const CenterOfMassComputationPtr_t& comc,
          const JointPtr_t& joint, const vector3_t reference,
          std::vector <bool> mask,
          const std::string& name);

      virtual std::ostream& print (std::ostream& o) const;
    protected:
      /// Compute value of error
      ///
      /// \param argument configuration of the robot,
      /// \retval result error vector
      virtual void impl_compute	(LiegroupElementRef result,
				 ConfigurationIn_t argument)
	const;
      virtual void impl_jacobian (matrixOut_t jacobian,
				  ConfigurationIn_t arg) const;

      bool isEqual(const DifferentiableFunction& other) const {
        const RelativeCom& castother = dynamic_cast<const RelativeCom&>(other);
        if (!DifferentiableFunction::isEqual(other))
          return false;
        
        if (robot_ != castother.robot_)
          return false;
        if (comc_ != castother.comc_)
          return false;
        if (joint_ != castother.joint_)
          return false;
        if (mask_ != castother.mask_)
          return false;
        if (nominalCase_ != castother.nominalCase_)
          return false;
        
        return true;
      }
    private:
      DevicePtr_t robot_;
      CenterOfMassComputationPtr_t comc_;
      JointPtr_t joint_;
      vector3_t reference_;
      std::vector <bool> mask_;
      bool nominalCase_;
      mutable ComJacobian_t jacobian_;

      RelativeCom() {}
      HPP_SERIALIZABLE();
    }; // class RelativeCom
    /// \}
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_RELATIVE_COM_HH
