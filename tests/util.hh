// Copyright (c) 2017, Joseph Mirabel
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
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

#ifndef TEST_UTIL_HH
# define TEST_UTIL_HH

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/extra-config-space.hh>
#include <hpp/pinocchio/joint.hh>

bool saturate (const hpp::pinocchio::DevicePtr_t& robot,
    hpp::pinocchio::vectorIn_t q,
    Eigen::VectorXi& sat)
{
  using hpp::pinocchio::size_type;
  bool ret = false;
  const se3::Model& model = robot->model();

  for (std::size_t i = 1; i < model.joints.size(); ++i) {
    const size_type nq = model.joints[i].nq();
    const size_type nv = model.joints[i].nv();
    const size_type idx_q = model.joints[i].idx_q();
    const size_type idx_v = model.joints[i].idx_v();
    for (size_type j = 0; j < nq; ++j) {
      const size_type iq = idx_q + j;
      const size_type iv = idx_v + std::min(j,nv-1);
      if        (q[iq] >= model.upperPositionLimit[iq]) {
        sat[iv] =  1;
        ret = true;
      } else if (q[iq] <= model.lowerPositionLimit[iq]) {
        sat[iv] = -1;
        ret = true;
      } else
        sat[iv] =  0;
    }
  }

  const hpp::pinocchio::ExtraConfigSpace& ecs = robot->extraConfigSpace();
  const size_type& d = ecs.dimension();

  for (size_type k = 0; k < d; ++k) {
    const size_type iq = model.nq + k;
    const size_type iv = model.nv + k;
    if        (q[iq] >= ecs.upper(k)) {
      sat[iv] =  1;
      ret = true;
    } else if (q[iq] <= ecs.lower(k)) {
      sat[iv] = -1;
      ret = true;
    } else
      sat[iv] =  0;
  }
  return ret;
}

#endif // TEST_UTIL_HH
