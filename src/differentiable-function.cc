// Copyright (c) 2015, Joseph Mirabel
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

#include <hpp/constraints/differentiable-function.hh>

#include <pinocchio/multibody/liegroup/liegroup.hpp>
#include <pinocchio/algorithm/finite-differences.hpp>

#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/joint-collection.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/liegroup.hh>

namespace hpp {
  namespace constraints {
    namespace {
      using hpp::pinocchio::DefaultLieGroupMap;
      typedef std::vector<pinocchio::JointIndex> JointIndexVector;

      struct FiniteDiffRobotOp
      {
        FiniteDiffRobotOp (const DevicePtr_t& r, const value_type& epsilon)
          : robot(r), model(robot->model()), 
          epsilon(epsilon),
          v(robot->numberDof())
        {}

        inline value_type step (const size_type& i, const vector_t& x) const
        {
          assert(i >= 0);
          value_type r = std::abs(x[i]);

          if (r == 0) return epsilon;
          else        return epsilon * r;
        }

        template <bool forward>
        inline void integrate (const vector_t& x, const vector_t& h, const size_type& /*i*/, vector_t& result) const
        {
          // Use only the joint corresponding to velocity index i
          if (forward)
            hpp::pinocchio::integrate<false, DefaultLieGroupMap> (robot, x,  h, result);
          else
            hpp::pinocchio::integrate<false, DefaultLieGroupMap> (robot, x, -h, result);
        }

        inline value_type difference (const vector_t& x0, const vector_t& x1, const size_type& i) const
        {
          hpp::pinocchio::difference <DefaultLieGroupMap> (robot, x0, x1, v);
          return v[i];
        }

        inline void reset (const vector_t& x, const size_type& /*i*/, vector_t& result) const
        {
          // Use only the joint corresponding to velocity index i
          result = x;
        }

        const DevicePtr_t& robot;
        const pinocchio::Model& model;
        const value_type& epsilon;
        mutable vector_t v;
      };

      struct FiniteDiffVectorSpaceOp
      {
        FiniteDiffVectorSpaceOp (const value_type& epsilon) : epsilon(epsilon) {}

        inline value_type step (const size_type i, const vector_t& x) const
        {
          const value_type r = std::abs(x[i]);

          if (r == 0) return epsilon;
          else        return epsilon * r;
        }

        template <bool forward>
        inline void integrate (const vector_t& x, const vector_t& h, const size_type& i, vector_t& result) const
        {
          result[i] = x[i] + (forward ? h[i] : -h[i]);
        }

        inline value_type difference (const vector_t& x0, const vector_t& x1, const size_type& i) const
        {
          return x0[i] - x1[i];
        }

        inline void reset (const vector_t& x, const size_type& i, vector_t& result) const
        {
          result[i] = x[i];
        }

        const value_type& epsilon;
      };

      template <typename FiniteDiffOp, typename Function>
        void finiteDiffCentral(matrixOut_t jacobian, vectorIn_t x,
            const FiniteDiffOp& op, const Function& f)
        {
          size_type n = jacobian.cols();
          vector_t x_pdx = x;
          vector_t x_mdx = x;
          vector_t h = vector_t::Zero (jacobian.cols());
          LiegroupElement f_x_mdx (f.outputSpace ()),
            f_x_pdx (f.outputSpace ());
          const ArrayXb& adp = f.activeDerivativeParameters();

          for (size_type j = 0; j < n; ++j) {
            if (!adp[j]) {
              jacobian.col (j).setZero();
              continue;
            }

            h[j] = op.step(j, x);

            op.template integrate<false>(x, h, j, x_mdx);
            f.value (f_x_mdx, x_mdx);

            op.template integrate<true >(x, h, j, x_pdx);
            f.value (f_x_pdx, x_pdx);

            jacobian.col (j) = ((f_x_pdx - f_x_mdx) / h[j]) / 2;

            op.reset(x, j, x_mdx);
            op.reset(x, j, x_pdx);
            h[j] = 0;
          }
          if (jacobian.hasNaN ()) {
            hppDout (error, "Central finite difference: NaN");
          }
        }

      template <typename FiniteDiffOp, typename Function>
        void finiteDiffForward(matrixOut_t jacobian, vectorIn_t x,
            const FiniteDiffOp& op, const Function& f)
        {
          size_type n = jacobian.cols();
          vector_t x_dx = x;
          vector_t h = vector_t::Zero (jacobian.cols());
          LiegroupElement f_x (f.outputSpace ()), f_x_pdx (f.outputSpace ());
          const ArrayXb& adp = f.activeDerivativeParameters();

          f.value (f_x, x);

          for (size_type j = 0; j < n; ++j) {
            if (!adp[j]) {
              jacobian.col (j).setZero();
              continue;
            }

            h[j] = op.step(j, x);

            op.template integrate<true >(x, h, j, x_dx);
            f.value (f_x_pdx, x_dx);

            jacobian.col (j) = (f_x_pdx - f_x) / h[j];

            op.reset(x, j, x_dx);
            h[j] = 0;
          }
          if (jacobian.hasNaN ()) {
            hppDout (warning, "Finite difference of \"" << f.name() << "\" has NaN values.");
          }
        }
    }

    void DifferentiableFunction::finiteDifferenceForward
      (matrixOut_t jacobian, vectorIn_t x,
       DevicePtr_t robot, value_type eps) const
      {
        if (robot)
          finiteDiffForward(jacobian, x, FiniteDiffRobotOp(robot, eps), *this);
        else
          finiteDiffForward(jacobian, x, FiniteDiffVectorSpaceOp(eps), *this);
      }

    void DifferentiableFunction::finiteDifferenceCentral
      (matrixOut_t jacobian, vectorIn_t x,
       DevicePtr_t robot, value_type eps) const
      {
        if (robot)
          finiteDiffCentral(jacobian, x, FiniteDiffRobotOp(robot, eps), *this);
        else
          finiteDiffCentral(jacobian, x, FiniteDiffVectorSpaceOp(eps), *this);
      }

    DifferentiableFunction::DifferentiableFunction
    (size_type sizeInput, size_type sizeInputDerivative,
     size_type sizeOutput, std::string name) :
      inputSize_ (sizeInput), inputDerivativeSize_ (sizeInputDerivative),
      outputSpace_ (LiegroupSpace::Rn (sizeOutput)),
      activeParameters_ (ArrayXb::Constant (sizeInput, true)),
      activeDerivativeParameters_
      (ArrayXb::Constant (sizeInputDerivative, true)),
      name_ (name)
      {
      }

    DifferentiableFunction::DifferentiableFunction
    (size_type sizeInput, size_type sizeInputDerivative,
     const LiegroupSpacePtr_t& outputSpace, std::string name) :
      inputSize_ (sizeInput), inputDerivativeSize_ (sizeInputDerivative),
      outputSpace_ (outputSpace), activeParameters_
      (ArrayXb::Constant (sizeInput, true)),
      activeDerivativeParameters_
      (ArrayXb::Constant (sizeInputDerivative, true)),
      name_ (name), context_ ()
    {
    }

    std::ostream& DifferentiableFunction::print (std::ostream& o) const
    {
      return o << "Differentiable function: " << name ();
    }
  } // namespace constraints
} // namespace hpp
