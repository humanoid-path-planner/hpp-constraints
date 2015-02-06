// Copyright (c) 2015, LAAS-CNRS
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

#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>
#include <hpp/model/configuration.hh>
#include <hpp/model/object-factory.hh>

#include <hpp/constraints/tools.hh>

#define BOOST_TEST_MODULE SymbolicCalculus
#include <boost/test/included/unit_test.hpp>

#include <stdlib.h>
#include <limits>
#include <math.h>

using namespace hpp::constraints;

struct DataWrapper {
  eigen::vector3_t value;
  JacobianMatrix jacobian;
};

class PointTester : public CalculusBase <PointTester>
{
  public:
    PointTester () : datas (NULL) {}

    PointTester (const CalculusBase<PointTester>& other) {
      const PointTester& o = static_cast <const PointTester&>(other);
      datas = o.datas;
    }

    PointTester (DataWrapper* d): datas (d)
    {}

    void computeValue () {
      this->value_ = datas->value;
    }
    void computeJacobian () {
      this->jacobian_ = datas->jacobian;
    }

    DataWrapper* datas;
};

namespace crossProduct {
  typedef Eigen::Matrix <value_type, 1, 6> Config;
  typedef Eigen::Matrix <value_type, 3, 1> Value;
  typedef Eigen::Matrix <value_type, 3, 6> Jacobian;

  void setWrappers (Config cfg, DataWrapper* d1, DataWrapper* d2)
  {
    d1->value = cfg.leftCols (3);
    d2->value = cfg.rightCols (3);
    d1->jacobian = Jacobian::Zero ();
    d2->jacobian = Jacobian::Zero ();
    d1->jacobian.leftCols (3).setIdentity ();
    d2->jacobian.rightCols (3).setIdentity ();
  }

  void value (Config c, Value& v)
  {
    v[0] =   c[1] * c[5] - c[2] * c[4];
    v[1] = - c[0] * c[5] + c[2] * c[3];
    v[2] =   c[0] * c[4] - c[1] * c[3];
  }

  void jacobian (Config cfg, Jacobian& j)
  {
    j (0,0) =      0 ; j (0,1) =  cfg[5]; j (0,2) = -cfg[4]; j (0,3) =      0 ; j (0,4) = -cfg[2]; j (0,5) =  cfg[1];
    j (1,0) = -cfg[5]; j (1,1) =      0 ; j (1,2) =  cfg[3]; j (1,3) =  cfg[2]; j (1,4) =      0 ; j (1,5) = -cfg[0];
    j (2,0) =  cfg[4]; j (2,1) = -cfg[3]; j (2,2) =      0 ; j (2,3) = -cfg[1]; j (2,4) =  cfg[0]; j (2,5) =      0 ;
  }
}

BOOST_AUTO_TEST_CASE (CrossProductTest) {
  DataWrapper* d1 = new DataWrapper ();
  DataWrapper* d2 = new DataWrapper ();
  PointTester p1 (d1), p2 (d2);
  typedef CrossProduct <PointTester, PointTester> CP_t;
  CP_t cp = p1 ^ p2;
  crossProduct::Value v;
  crossProduct::Jacobian j;
  for (size_t i = 0; i < 100; i++) {
    crossProduct::Config cfg = crossProduct::Config::Random ();
    crossProduct::setWrappers (cfg, d1, d2);
    cp.computeValue ();
    crossProduct::value (cfg, v);
    crossProduct::jacobian (cfg, j);
    BOOST_CHECK (cp.value ().isApprox (v));
    cp.computeJacobian ();
    BOOST_CHECK (cp.jacobian ().isApprox (j));
  }
  delete d1;
  delete d2;
}

namespace difference {
  typedef Eigen::Matrix <value_type, 1, 6> Config;
  typedef Eigen::Matrix <value_type, 3, 1> Value;
  typedef Eigen::Matrix <value_type, 3, 6> Jacobian;

  void setWrappers (Config cfg, DataWrapper* d1, DataWrapper* d2)
  {
    d1->value = cfg.leftCols (3);
    d2->value = cfg.rightCols (3);
    d1->jacobian = Jacobian::Identity ();
    d2->jacobian = Jacobian::Zero ();
    d2->jacobian.rightCols (3).setIdentity ();
  }

  void value (Config c, Value& v)
  {
    v = c.leftCols (3) - c.rightCols (3);
  }

  void jacobian (Config , Jacobian& j)
  {
    j.rightCols (3).setIdentity ();
    j = -j;
    j.leftCols (3).setIdentity ();
  }
}

BOOST_AUTO_TEST_CASE (DifferenceTest) {
  DataWrapper* d1 = new DataWrapper ();
  DataWrapper* d2 = new DataWrapper ();
  PointTester p1 (d1), p2 (d2);
  typedef Difference <PointTester, PointTester> D_t;
  D_t cp = p1 - p2;
  difference::Value v;
  difference::Jacobian j;
  for (size_t i = 0; i < 100; i++) {
    difference::Config cfg = difference::Config::Random ();
    difference::setWrappers (cfg, d1, d2);
    cp.computeValue ();
    difference::value (cfg, v);
    difference::jacobian (cfg, j);
    BOOST_CHECK (cp.value ().isApprox (v));
    cp.computeJacobian ();
    BOOST_CHECK (cp.jacobian ().isApprox (j));
  }
  delete d1;
  delete d2;
}

namespace sum {
  typedef Eigen::Matrix <value_type, 1, 6> Config;
  typedef Eigen::Matrix <value_type, 3, 1> Value;
  typedef Eigen::Matrix <value_type, 3, 6> Jacobian;

  void setWrappers (Config cfg, DataWrapper* d1, DataWrapper* d2)
  {
    d1->value = cfg.leftCols (3);
    d2->value = cfg.rightCols (3);
    d1->jacobian = Jacobian::Identity ();
    d2->jacobian = Jacobian::Zero ();
    d2->jacobian.rightCols (3).setIdentity ();
  }

  void value (Config c, Value& v)
  {
    v = c.leftCols (3) + c.rightCols (3);
  }

  void jacobian (Config , Jacobian& j)
  {
    j.rightCols (3).setIdentity ();
    j.leftCols (3).setIdentity ();
  }
}

BOOST_AUTO_TEST_CASE (SumTest) {
  DataWrapper* d1 = new DataWrapper ();
  DataWrapper* d2 = new DataWrapper ();
  PointTester p1 (d1), p2 (d2);
  typedef Sum <PointTester, PointTester> D_t;
  D_t cp = p1 + p2;
  sum::Value v;
  sum::Jacobian j;
  for (size_t i = 0; i < 100; i++) {
    sum::Config cfg = sum::Config::Random ();
    sum::setWrappers (cfg, d1, d2);
    cp.computeValue ();
    sum::value (cfg, v);
    sum::jacobian (cfg, j);
    BOOST_CHECK (cp.value ().isApprox (v));
    cp.computeJacobian ();
    BOOST_CHECK (cp.jacobian ().isApprox (j));
  }
  delete d1;
  delete d2;
}

namespace scalarMultiply {
  typedef Eigen::Matrix <value_type, 1, 4> Config;
  typedef Eigen::Matrix <value_type, 3, 1> Value;
  typedef Eigen::Matrix <value_type, 3, 3> Jacobian;

  void setWrappers (Config cfg, DataWrapper* d1, DataWrapper* )
  {
    d1->value = cfg.leftCols (3);
    d1->jacobian = Jacobian::Identity ();
  }

  void value (Config c, Value& v)
  {
    v = c[3] * c.leftCols (3);
  }

  void jacobian (Config c, Jacobian& j)
  {
    j.setIdentity ();
    j = c[3] * j;
  }
}

BOOST_AUTO_TEST_CASE (ScalarMultiplyTest) {
  DataWrapper* d1 = new DataWrapper ();
  DataWrapper* d2 = new DataWrapper ();
  PointTester p1 (d1), p2 (d2);
  typedef ScalarMultiply <PointTester> D_t;
  value_type scalar = 3;
  D_t cp = p1 * scalar;
  scalarMultiply::Value v;
  scalarMultiply::Jacobian j;
  for (size_t i = 0; i < 100; i++) {
    scalarMultiply::Config cfg = scalarMultiply::Config::Random ();
    cfg[3] = scalar;
    scalarMultiply::setWrappers (cfg, d1, d2);
    cp.computeValue ();
    scalarMultiply::value (cfg, v);
    scalarMultiply::jacobian (cfg, j);
    BOOST_CHECK (cp.value ().isApprox (v));
    cp.computeJacobian ();
    BOOST_CHECK (cp.jacobian ().isApprox (j));
  }
  delete d1;
  delete d2;
}
