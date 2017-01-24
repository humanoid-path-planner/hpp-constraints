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

#define BOOST_TEST_MODULE SymbolicCalculus
#include <boost/test/included/unit_test.hpp>

#include <stdlib.h>
#include <limits>
#include <math.h>

#include <Eigen/Geometry>

#include <hpp/constraints/symbolic-calculus.hh>
#include <hpp/constraints/symbolic-function.hh>

using namespace hpp::constraints;

typedef SymbolicFunction<JointFrame> JointFrameFunction;

template <class ValueType = eigen::vector3_t,
          class JacobianType = JacobianMatrix >
class PointTesterT : public CalculusBase <PointTesterT<ValueType, JacobianType>, ValueType, JacobianType>
{
  public:
    typedef CalculusBase <PointTesterT<ValueType, JacobianType>, ValueType, JacobianType> Parent_t;
    struct DataWrapper {
      ValueType value;
      JacobianType jacobian;
      
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    HPP_CONSTRAINTS_CB_CREATE1 (PointTesterT, DataWrapper*)

    PointTesterT (DataWrapper* d = NULL) :
      datas (d)
    {}

    PointTesterT (const Parent_t& other) :
      Parent_t (other),
      datas (static_cast <const PointTesterT&>(other).datas)
    {
    }

    void impl_value () {
      this->value_ = datas->value;
    }
    void impl_jacobian () {
      this->jacobian_ = datas->jacobian;
    }

    DataWrapper* datas;
};

namespace crossProduct {
  typedef Eigen::Matrix <value_type, 1, 6> Config;
  typedef Eigen::Matrix <value_type, 3, 1> Value;
  typedef Eigen::Matrix <value_type, 3, 6, Eigen::RowMajor> Jacobian;
  typedef PointTesterT <Value, Jacobian> PointTester;
  typedef PointTester::DataWrapper DataWrapper;

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

  void jacobian (Config cfg, Eigen::Ref<Jacobian> j)
  {
    j (0,0) =      0 ; j (0,1) =  cfg[5]; j (0,2) = -cfg[4]; j (0,3) =      0 ; j (0,4) = -cfg[2]; j (0,5) =  cfg[1];
    j (1,0) = -cfg[5]; j (1,1) =      0 ; j (1,2) =  cfg[3]; j (1,3) =  cfg[2]; j (1,4) =      0 ; j (1,5) = -cfg[0];
    j (2,0) =  cfg[4]; j (2,1) = -cfg[3]; j (2,2) =      0 ; j (2,3) = -cfg[1]; j (2,4) =  cfg[0]; j (2,5) =      0 ;
  }
}

BOOST_AUTO_TEST_CASE (CrossProductTest) {
  using namespace crossProduct;
  DataWrapper* d1 = new DataWrapper ();
  DataWrapper* d2 = new DataWrapper ();
  Traits<PointTester>::Ptr_t p1 = PointTester::create (d1),
    p2 = PointTester::create (d2);
  typedef CrossProduct <PointTester, PointTester> CP_t;
  Traits<CP_t>::Ptr_t cp = p1 ^ p2;
  Value v;
  Jacobian j;
  for (size_t i = 0; i < 100; i++) {
    Config cfg = Config::Random ();
    cp->invalidate ();
    setWrappers (cfg, d1, d2);
    cp->computeValue ();
    value (cfg, v);
    jacobian (cfg, j);
    BOOST_CHECK (cp->value ().isApprox (v));
    cp->computeJacobian ();
    BOOST_CHECK (cp->jacobian ().isApprox (j));
  }
  delete d1;
  delete d2;
}

namespace difference {
  typedef Eigen::Matrix <value_type, 1, 6> Config;
  typedef Eigen::Matrix <value_type, 3, 1> Value;
  typedef Eigen::Matrix <value_type, 3, 6, Eigen::RowMajor> Jacobian;
  typedef PointTesterT <Value, Jacobian> PointTester;
  typedef PointTester::DataWrapper DataWrapper;

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
  using namespace difference;
  DataWrapper* d1 = new DataWrapper ();
  DataWrapper* d2 = new DataWrapper ();
  Traits<PointTester>::Ptr_t p1 = PointTester::create (d1), p2 = PointTester::create (d2);
  typedef Difference <PointTester, PointTester> D_t;
  Traits<D_t>::Ptr_t cp = p1 - p2;
  Value v;
  Jacobian j;
  for (size_t i = 0; i < 100; i++) {
    Config cfg = Config::Random ();
    cp->invalidate ();
    setWrappers (cfg, d1, d2);
    cp->computeValue ();
    value (cfg, v);
    jacobian (cfg, j);
    BOOST_CHECK (cp->value ().isApprox (v));
    cp->computeJacobian ();
    BOOST_CHECK (cp->jacobian ().isApprox (j));
  }
  delete d1;
  delete d2;
}

namespace sum {
  typedef Eigen::Matrix <value_type, 1, 6> Config;
  typedef Eigen::Matrix <value_type, 3, 1> Value;
  typedef Eigen::Matrix <value_type, 3, 6, Eigen::RowMajor> Jacobian;
  typedef PointTesterT <Value, Jacobian> PointTester;
  typedef PointTester::DataWrapper DataWrapper;

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
  using namespace sum;
  DataWrapper* d1 = new DataWrapper ();
  DataWrapper* d2 = new DataWrapper ();
  Traits<PointTester>::Ptr_t p1 = PointTester::create (d1), p2 = PointTester::create (d2);
  typedef Sum <PointTester, PointTester> D_t;
  Traits<D_t>::Ptr_t cp = p1 + p2;
  Value v;
  Jacobian j;
  for (size_t i = 0; i < 100; i++) {
    Config cfg = Config::Random ();
    cp->invalidate ();
    setWrappers (cfg, d1, d2);
    cp->computeValue ();
    value (cfg, v);
    jacobian (cfg, j);
    BOOST_CHECK (cp->value ().isApprox (v));
    cp->computeJacobian ();
    BOOST_CHECK (cp->jacobian ().isApprox (j));
  }
  delete d1;
  delete d2;
}

namespace scalarMultiply {
  typedef Eigen::Matrix <value_type, 1, 4> Config;
  typedef Eigen::Matrix <value_type, 3, 1> Value;
  typedef Eigen::Matrix <value_type, 3, 3, Eigen::RowMajor> Jacobian;
  typedef PointTesterT <Value, Jacobian> PointTester;
  typedef PointTester::DataWrapper DataWrapper;

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
  using namespace scalarMultiply;
  DataWrapper* d1 = new DataWrapper ();
  DataWrapper* d2 = new DataWrapper ();
  Traits<PointTester>::Ptr_t p1 = PointTester::create (d1), p2 = PointTester::create (d2);
  typedef ScalarMultiply <PointTester> D_t;
  value_type scalar = 3;
  Traits<D_t>::Ptr_t cp = scalar * p1;
  Value v;
  Jacobian j;
  for (size_t i = 0; i < 100; i++) {
    Config cfg = Config::Random ();
    cfg[3] = scalar;
    cp->invalidate ();
    setWrappers (cfg, d1, d2);
    cp->computeValue ();
    value (cfg, v);
    jacobian (cfg, j);
    BOOST_CHECK (cp->value ().isApprox (v));
    cp->computeJacobian ();
    BOOST_CHECK (cp->jacobian ().isApprox (j));
  }
  delete d1;
  delete d2;
}

namespace scalarProduct {
  typedef Eigen::Matrix <value_type, 1, 6> Config;
  typedef Eigen::Matrix <value_type, 1, 1> OutputValue;
  typedef Eigen::Matrix <value_type, 1, 6, Eigen::RowMajor> OutputJacobian;
  typedef Eigen::Matrix <value_type, 3, 1> Value;
  typedef Eigen::Matrix <value_type, 3, 6, Eigen::RowMajor> Jacobian;
  typedef PointTesterT <Value, Jacobian> PointTester;
  typedef PointTester::DataWrapper DataWrapper;

  void setWrappers (Config cfg, DataWrapper* d1, DataWrapper* d2)
  {
    d1->value = cfg.leftCols (3);
    d2->value = cfg.rightCols (3);
    d1->jacobian = Jacobian::Zero ();
    d2->jacobian = Jacobian::Zero ();
    d1->jacobian.leftCols (3).setIdentity ();
    d2->jacobian.rightCols (3).setIdentity ();
  }

  void value (Config c, OutputValue& v)
  {
    v[0] = c.segment (0,3).dot (c.segment (3,3));
  }

  void jacobian (Config cfg, OutputJacobian& j)
  {
    j (0,0) = cfg[3]; j (0,1) = cfg[4]; j (0,2) = cfg[5]; j (0,3) = cfg[0]; j (0,4) = cfg[1]; j (0,5) = cfg[2];
  }
}

BOOST_AUTO_TEST_CASE (ScalarProductTest) {
  using namespace scalarProduct;
  DataWrapper* d1 = new DataWrapper ();
  DataWrapper* d2 = new DataWrapper ();
  Traits<PointTester>::Ptr_t p1 = PointTester::create (d1), p2 = PointTester::create (d2);
  typedef ScalarProduct <PointTester, PointTester> SP_t;
  Traits<SP_t>::Ptr_t sp = p1 * p2;
  OutputValue v;
  OutputJacobian j;
  for (size_t i = 0; i < 100; i++) {
    Config cfg = Config::Random ();
    sp->invalidate ();
    setWrappers (cfg, d1, d2);
    sp->computeValue ();
    value (cfg, v);
    jacobian (cfg, j);
    BOOST_CHECK (v.isApproxToConstant (sp->value ()));
    sp->computeJacobian ();
    BOOST_CHECK (sp->jacobian ().isApprox (j));
  }
  delete d1;
  delete d2;
}

namespace matrixOfExp {
  typedef Eigen::Matrix <value_type, 1, 6> Config;
  typedef Eigen::Matrix <value_type, 3, 1> Value;
  typedef Eigen::Matrix <value_type, 3, Eigen::Dynamic, Eigen::RowMajor> Jacobian;
  typedef Eigen::Matrix <value_type, 6, 2> OutputValue;
  typedef Eigen::Matrix <value_type, 6, 2*6, Eigen::RowMajor> OutputJacobian;
  typedef Eigen::Matrix <value_type, 2, 1> RValue;
  typedef Eigen::Matrix <value_type, 6, 6> OutputJacobianTimesRValue;
  typedef PointTesterT <Value, Jacobian> PointTester;
  typedef PointTester::DataWrapper DataWrapper;

  void setWrappers (Config cfg, DataWrapper* d1, DataWrapper* d2)
  {
    d1->value = cfg.leftCols (3);
    d2->value = cfg.rightCols (3);
    d1->jacobian = Jacobian::Zero (3,6);
    d2->jacobian = Jacobian::Zero (3,6);
    d1->jacobian.leftCols (3).setIdentity ();
    d2->jacobian.rightCols (3).setIdentity ();
  }

  void value (Config c, OutputValue& v)
  {
    v.block (0,0,3,1).transpose() = c.segment<3> (0);
    v.block (3,0,3,1).transpose() = c.segment<3> (0).cross (c.segment<3> (3));
    v.block (0,1,3,1).transpose() = c.segment (0,3) - c.segment (3,3);
    v.block (3,1,3,1).transpose() = c.segment (0,3) + c.segment (3,3);
  }

  void jacobian (Config cfg, OutputJacobian& j)
  {
    j.block<3,3> (0,0).setIdentity();
    j.block<3,3> (0,3).setZero();

    crossProduct::jacobian (cfg, j.block<3,6> (3,0));
    //j.block (0,3,3,3).setZero();

    j.block<3,3> (0,6).setIdentity();
    j.block<3,3> (0,9).setIdentity();
    j.block<3,3> (0,9) *= -1;

    j.block<3,3> (3,6).setIdentity();
    j.block<3,3> (3,9).setIdentity();
  }

  void jacobianTimes (Eigen::Ref <const OutputJacobian> jin,
      Eigen::Ref <const RValue > rvalue, Eigen::Ref <OutputJacobianTimesRValue> jout)
  {
    jout.block (0,0,3,6) = rvalue[0] * jin.block (0,0,3,6)
                         + rvalue[1] * jin.block (0,6,3,6);
    jout.block (3,0,3,6) = rvalue[0] * jin.block (3,0,3,6)
                         + rvalue[1] * jin.block (3,6,3,6);
  }
}

BOOST_AUTO_TEST_CASE (MatrixOfExpTest) {
  using namespace matrixOfExp;
  DataWrapper* d1 = new DataWrapper ();
  DataWrapper* d2 = new DataWrapper ();
  Traits<PointTester>::Ptr_t p1 = PointTester::create (d1), p2 = PointTester::create (d2);
  typedef MatrixOfExpressions <Value, Jacobian> MoE_t;
  OutputValue v;
  OutputJacobian j;
  RValue rvalue;
  OutputJacobianTimesRValue jout, cache;
  MoE_t moe(v,j);
  moe.setSize (2,2);
  moe.set(0,0,p1);
  moe.set(1,0,p1 ^ p2);
  moe.set(0,1,p1 - p2);
  moe.set(1,1,p1 + p2);
  for (size_t i = 0; i < 100; i++) {
    Config cfg = Config::Random ();
    moe.invalidate ();
    setWrappers (cfg, d1, d2);
    moe.computeValue ();
    value (cfg, v);
    jacobian (cfg, j);
    BOOST_CHECK (moe.value().isApprox (v));
    moe.computeJacobian ();
    BOOST_CHECK (moe.jacobian ().isApprox (j));
    for (size_t l = 0; l < 100; l++) {
      rvalue = RValue::Random ();
      moe.jacobianTimes (rvalue, cache);
      jacobianTimes (j, rvalue, jout);
      BOOST_CHECK (cache.isApprox (jout));
    }
  }
  delete d1;
  delete d2;
}
