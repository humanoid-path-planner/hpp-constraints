// Copyright (c) 2018, LAAS-CNRS
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

#include <hpp/pinocchio/device.hh>
#include <hpp/constraints/explicit/relative-pose.hh>
#include <hpp/constraints/explicit/relative-transformation.hh>
#include <hpp/constraints/matrix-view.hh>

namespace hpp {
  namespace constraints {
    namespace explicit_ {
      LiegroupSpacePtr_t RelativePose::SE3 (LiegroupSpace::SE3 ());
      LiegroupSpacePtr_t RelativePose::R3xSO3 (LiegroupSpace::R3xSO3 ());

      namespace {
        using Eigen::BlockIndex;
        BlockIndex::segments_t vectorOfBoolToIntervals (std::vector<bool>& v)
        {
          BlockIndex::segments_t ret;
          for (std::size_t i = 0; i < v.size(); ++i)
            if (v[i]) ret.push_back(BlockIndex::segment_t (i, 1));
          BlockIndex::shrink (ret);
          return ret;
        }

        // This function computes the configuration variables of
        //   - joint1 and its parents up to the root joint,
        //   - joint2 parent and its parent up to root joint.
        // Note that configuration variables of joint2 do not belong to
        // the resulting set.
        BlockIndex::segments_t inputConfVariables
        (const DevicePtr_t& robot, JointConstPtr_t joint1,
         JointConstPtr_t joint2)
        {
          std::vector<bool> conf (robot->configSize(), false);
          while (joint1 && joint1->index() != 0) {
            for (size_type i = 0; i < joint1->configSize(); ++i)
              conf[joint1->rankInConfiguration() + i] =
                !conf[joint1->rankInConfiguration() + i];
            hppDout (info, "Adding joint1 " << joint1->name ()
                     << " as input variable.");
            joint1 = joint1->parentJoint();
          }
          joint2 = joint2->parentJoint ();
          while (joint2 && joint2->index() != 0) {
            for (size_type i = 0; i < joint2->configSize(); ++i)
              conf[joint2->rankInConfiguration() + i] =
                !conf[joint2->rankInConfiguration() + i];
            hppDout (info, "Adding joint2 " << joint2->name ()
                     << " as input variable.");
            joint2 = joint2->parentJoint();
          }
          return vectorOfBoolToIntervals (conf);
        }

        // This function computes the velocity variables of
        //   - joint1 and its parents up to the root joint,
        //   - joint2 parent and its parent up to root joint.
        // Note that velocity variables of joint2 do not belong to
        // the resulting set.
        BlockIndex::segments_t inputVelocityVariables
        (const DevicePtr_t& robot, JointConstPtr_t joint1,
         JointConstPtr_t joint2)
        {
          std::vector<bool> vel (robot->numberDof(), false);
          while (joint1 && joint1->index() != 0) {
            for (size_type i = 0; i < joint1->numberDof(); ++i)
              vel[joint1->rankInVelocity() + i] =
                !vel[joint1->rankInVelocity() + i];
            hppDout (info, "Adding joint1 " << joint1->name ()
                     << " as input variable.");
            joint1 = joint1->parentJoint();
          }
          joint2 = joint2->parentJoint ();
          while (joint2 && joint2->index() != 0) {
            for (size_type i = 0; i < joint2->numberDof(); ++i)
              vel[joint2->rankInVelocity() + i] =
                !vel[joint2->rankInVelocity() + i];
            hppDout (info, "Adding joint2 " << joint2->name ()
                     << " as input variable.");
            joint2 = joint2->parentJoint();
          }
          return vectorOfBoolToIntervals (vel);
        }

        // Return configuration variables of a joint
        BlockIndex::segments_t jointConfInterval (JointConstPtr_t j) {
          return BlockIndex::segments_t(1, BlockIndex::segment_t
                                        (j->rankInConfiguration(),
                                         j->configSize()));
        }
        // Return velocity variables of a joint
        BlockIndex::segments_t jointVelInterval (JointConstPtr_t j) {
          return BlockIndex::segments_t(1, BlockIndex::segment_t
                                        (j->rankInVelocity(), j->numberDof()));
        }
      }

      RelativePosePtr_t RelativePose::create
        (const std::string& name, const DevicePtr_t& robot,
         const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
         const Transform3f& frame1, const Transform3f& frame2,
         std::vector <bool> mask, ComparisonTypes_t comp, vectorIn_t rhs)
      {
        RelativePose* ptr (new RelativePose
                           (name, robot, joint1, joint2, frame1, frame2,
                            mask, comp, rhs));
        RelativePosePtr_t shPtr (ptr);
        RelativePoseWkPtr_t wkPtr (shPtr);
        ptr->init (wkPtr);
        return shPtr;
      }

      RelativePosePtr_t RelativePose::createCopy
      (const RelativePosePtr_t& other)
      {
        RelativePose* ptr (new RelativePose (*other));
        RelativePosePtr_t shPtr (ptr);
        RelativePoseWkPtr_t wkPtr (shPtr);
        ptr->init (wkPtr);
        return shPtr;
      }

      ImplicitPtr_t RelativePose::copy () const
      {
        return createCopy (weak_.lock ());
      }

      void RelativePose::implicitToExplicitRhs (vectorIn_t implicitRhs,
                                                vectorOut_t explicitRhs) const
      {
        assert (implicitRhs.size () == 6);
        assert (explicitRhs.size () == 6);
        // p1 = exp_{R^3xSO(3)} (implicitRhs)
        LiegroupElement p1 (R3xSO3->exp (implicitRhs));
        // convert p1 to Transform3f M1
        Transform3f M1 ((Quaternion_t (p1.vector ().tail <4> ()))
                        .toRotationMatrix (), p1.vector ().head <3> ());
        // M2 = F_{2/J_2} M1 F_{2/J_2}^{-1}
        Transform3f M2 (frame2_ * M1 * frame2_.inverse ());
        // convert M2 to LiegroupElement p2
        vector7_t v2;
        v2.head <3> () = M2.translation ();
        v2.tail <4> () = Quaternion_t (M2.rotation ()).coeffs ();
        LiegroupElement p2 (v2, SE3);
        explicitRhs = p2 - SE3->neutral ();
      }

      void RelativePose::explicitToImplicitRhs (vectorIn_t explicitRhs,
                                                vectorOut_t implicitRhs) const
      {
        assert (implicitRhs.size () == 6);
        assert (explicitRhs.size () == 6);
        // p1 = exp_{SE(3)} (explicitRhs)
        LiegroupElement p1 (SE3->exp (explicitRhs));
        // convert p1 to Transform3f M1
        Transform3f M1 ((Quaternion_t (p1.vector ().tail <4> ()))
                        .toRotationMatrix (), p1.vector ().head <3> ());
        // M2 = F_{2/J_2}^{-1} M1 F_{2/J_2}
        Transform3f M2 (frame2_.inverse () * M1 * frame2_);
        // convert M2 to LiegroupElement p2
        vector7_t v2;
        v2.head <3> () = M2.translation ();
        v2.tail <4> () = Quaternion_t (M2.rotation ()).coeffs ();
        LiegroupElement p2 (v2, R3xSO3);
        implicitRhs = p2 - R3xSO3->neutral ();
      }

      RelativePose::RelativePose
      (const std::string& name, const DevicePtr_t& robot,
       const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
       const Transform3f& frame1, const Transform3f& frame2,
       std::vector <bool> mask, ComparisonTypes_t comp, vectorIn_t rhs) :
        Implicit (GenericTransformation
                  <RelativeBit | PositionBit | OrientationBit>::create
                  (name, robot, joint1, joint2, frame1, frame2, mask),
                  comp, rhs),
        Explicit (robot->configSpace (), RelativeTransformation::create
                  (name, robot, joint1, joint2, frame1, frame2),
                  inputConfVariables (robot, joint1, joint2),
                  jointConfInterval (joint2),
                  inputVelocityVariables (robot, joint1, joint2),
                  jointVelInterval (joint2), comp),
        implicit::RelativePose (name, robot, joint1, joint2, frame1, frame2,
                                mask, comp, rhs),
        frame1_ (frame1), frame2_ (frame2)
        {
          assert ((*joint2->configurationSpace () ==
                   *pinocchio::LiegroupSpace::SE3 ()) ||
                  (*joint2->configurationSpace () ==
                   *pinocchio::LiegroupSpace::R3xSO3 ()));
        }

      RelativePose::RelativePose (const RelativePose& other) :
        Implicit (other), Explicit (other), implicit::RelativePose (other),
        frame1_ (other.frame1_), frame2_ (other.frame2_)
      {
      }

      void RelativePose::init (RelativePoseWkPtr_t weak)
      {
        Implicit::init (weak);
        Explicit::init (weak);
        implicit::RelativePose::init (weak);
        weak_ = weak;
      }
  } // namespace explicit_
  } // namespace constraints
} // namespace hpp
