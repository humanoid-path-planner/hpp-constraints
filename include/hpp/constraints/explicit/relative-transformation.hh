// Copyright (c) 2015 - 2018, CNRS
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr), Florent Lamiraux
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

#ifndef HPP_CONSTRAINTS_EXPLICIT_RELATIVE_TRANSFORMATION_HH
# define HPP_CONSTRAINTS_EXPLICIT_RELATIVE_TRANSFORMATION_HH

# include <hpp/constraints/matrix-view.hh>
# include <hpp/constraints/generic-transformation.hh>
# include <hpp/constraints/explicit.hh>

namespace hpp {
  namespace constraints {
    using constraints::RelativeTransformation;
    using constraints::RelativeTransformationPtr_t;

    namespace explicit_ {
    /// \addtogroup constraints
    /// \{

    /// Relative transformation as a mapping between configuration variables
    ///
    /// When the positions of two joints are constrained by a full
    /// transformation constraint, if the second joint is hold by a freeflyer
    /// joint, the position of this latter joint can be
    /// explicitely expressed with respect to the position of the first joint.
    ///
    /// This class provides this expression. The input configuration variables
    /// are the joint values of all joints except the above mentioned freeflyer
    /// joint. The output configuration variables are the 7 configuration
    /// variables of the freeflyer joint.
    ///
    class HPP_CONSTRAINTS_DLLAPI RelativeTransformation :
      public DifferentiableFunction
    {
    public:
      /// Return a shared pointer to a new instance
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraints is applied to,
      /// \param joint1 input joint
      /// \param joint2 output joint: position of this joint is computed with
      ///        respect to joint1 position
      /// \param frame1 position of a fixed frame in joint 1,
      /// \param frame2 position of a fixed frame in joint 2,
      static RelativeTransformationPtr_t create
	(const std::string& name      , const DevicePtr_t& robot,
	 const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
	 const Transform3f& frame1    , const Transform3f& frame2);

      /// Get joint 1
      const JointConstPtr_t& joint1 () const
      {
        return joint1_;
      }

      /// Get joint 2
      const JointConstPtr_t& joint2 () const
      {
        return joint2_;
      }

    protected:
      typedef Eigen::BlockIndex BlockIndex;
      typedef Eigen::RowBlockIndices RowBlockIndices;
      typedef Eigen::ColBlockIndices ColBlockIndices;

      RelativeTransformation (
          const std::string& name      , const DevicePtr_t& robot,
          const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
          const Transform3f& frame1    , const Transform3f& frame2,
          const segments_t inConf , const segments_t outConf,
          const segments_t inVel  , const segments_t outVel ,
          std::vector <bool> mask = std::vector<bool>(6,true));

      RelativeTransformation
	(const RelativeTransformation& other) :
	DifferentiableFunction (other),
	robot_ (other.robot_),
        parentJoint_ (other.parentJoint_),
        inConf_ (other.inConf_),   inVel_  (other.inVel_),
        outConf_ (other.outConf_), outVel_ (other.outVel_),
        F1inJ1_invF2inJ2_ (other.F1inJ1_invF2inJ2_)
	{
	}

      // Store weak pointer to itself
      void init (const RelativeTransformationWkPtr_t& weak)
      {
        weak_ = weak;
      }

      /// Compute the value (dimension 7) of the freeflyer joint 2
      ///
      /// \param argument vector of input configuration variables (all joints
      ///        except freeflyer joint)
      /// \retval result vector of output configuration variables corresponding
      ///         to the freeflyer value.
      void impl_compute (LiegroupElementRef result, vectorIn_t argument) const;

      void impl_jacobian (matrixOut_t jacobian, vectorIn_t arg) const;

    private:
      void forwardKinematics (vectorIn_t arg) const;

      DevicePtr_t robot_;
      // Parent of the R3 joint.
      JointConstPtr_t parentJoint_;
      JointConstPtr_t joint1_, joint2_;
      Transform3f frame1_, frame2_;
      RowBlockIndices inConf_;
      ColBlockIndices inVel_;
      RowBlockIndices outConf_ , outVel_;
      Transform3f F1inJ1_invF2inJ2_;

      RelativeTransformationWkPtr_t weak_;

      // Tmp variables
      mutable vector_t qsmall_, q_;
      mutable matrix_t tmpJac_, J2_parent_minus_J1_;
    }; // class RelativeTransformation
    /// \}

    } // namespace explicit_
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_EXPLICIT_RELATIVE_TRANSFORMATION_HH
