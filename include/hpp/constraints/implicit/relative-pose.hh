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

#ifndef HPP_CONSTRAINTS_IMPLICIT_RELATIVE_POSE_HH
# define HPP_CONSTRAINTS_IMPLICIT_RELATIVE_POSE_HH

# include <hpp/constraints/implicit.hh>

namespace hpp {
  namespace constraints {
    namespace implicit {
      /// Constraint of relative pose between two frames on a kinematic chain
      class HPP_CONSTRAINTS_DLLAPI RelativePose : public virtual Implicit
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
        /// \param mask vector of 6 boolean defining which coordinates of the
        ///        error vector to take into account.
        /// \note if joint1 is 0x0, joint 1 frame is considered to be the global
        ///       frame.
        static RelativePosePtr_t create
          (const std::string& name, const DevicePtr_t& robot,
           const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
           const Transform3f& frame1, const Transform3f& frame2,
           std::vector <bool> mask = std::vector<bool>(6,true),
           ComparisonTypes_t comp = std::vector <ComparisonType> (),
           vectorIn_t rhs = vector_t ());

        static RelativePosePtr_t createCopy (const RelativePosePtr_t& other);

        /// Get joint1
        const JointConstPtr_t& joint1 () const
        {
          return joint1_;
        }

        /// Get joint2
        const JointConstPtr_t& joint2 () const
        {
          return joint2_;
        }

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
        /// \param mask vector of 6 boolean defining which coordinates of the
        ///        error vector to take into account.
        /// \note if joint1 is 0x0, joint 1 frame is considered to be the global
        ///       frame.
        RelativePose (const std::string& name, const DevicePtr_t& robot,
                      const JointConstPtr_t& joint1,
                      const JointConstPtr_t& joint2,
                      const Transform3f& frame1, const Transform3f& frame2,
                      std::vector <bool> mask = std::vector<bool>(6,true),
                      ComparisonTypes_t comp = std::vector <ComparisonType> (),
                      vectorIn_t rhs = vector_t ());

        /// Copy constructor
        RelativePose (const RelativePose& other);
        /// Store shared pointer to itself
        void init (RelativePoseWkPtr_t weak);
      private:
        JointConstPtr_t joint1_, joint2_;
        RelativePoseWkPtr_t weak_;
      }; // class RelativePose
    } // namespace implicit
  } // namespace constraints
} // namespace hpp
#endif //HPP_CONSTRAINTS_IMPLICIT_RELATIVE_POSE_HH
