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

#ifndef HPP_CONSTRAINTS_EXPLICIT_SOLVER_HH
#define HPP_CONSTRAINTS_EXPLICIT_SOLVER_HH

#include <vector>

#include <boost/function.hpp>

#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/config.hh>

#include <hpp/constraints/matrix-view.hh>
#include <hpp/constraints/differentiable-function-stack.hh>

namespace hpp {
  namespace constraints {
    void difference (const DevicePtr_t& robot,
        const Eigen::BlockIndex<size_type>::vector_t& indexes,
        vectorIn_t arg0,
        vectorIn_t arg1,
        vectorOut_t result);

    /// \addtogroup solvers
    /// \{

    /// Solve system of explicit functions
    ///
    /// The solver works on a given set of variables \f$ X = (x_i) \f$.
    /// It contains a set of functions \f$ f_j \f$ that takes as input a subset of \f$ X \f$ and
    /// outputs values of corresponding to another subset of \f$ X \f$.
    /// There can be no cycles in the dependencies. Moreover, two functions must have
    /// non-intersecting output subset.
    /// For instance, \f$ (x_0, x_2) = f_0( x_1, x_3 ) \f$ and
    /// \f$ (x_3) = f_1( x_4 ) \f$ is a valid input. It would not be possible to
    /// add \f$ (x_0) = f_2( x_2 ) \f$ because it would introduce a cycle,
    /// or \f$ (x_3) = f_3( x_1 ) \f$ because \f$ x_3 \f$ would be computed by
    /// two different function.
    ///
    /// The resolution consists in modyfing the output values of each function,
    /// while respecting the dependendy order. Considering \f$ f_0, f_1 \f$
    /// above, \f$ f_1 \f$ must be computed before.
    class HPP_CONSTRAINTS_DLLAPI ExplicitSolver
    {
      public:
        typedef Eigen::RowBlockIndexes RowBlockIndexes;
        typedef Eigen::ColBlockIndexes ColBlockIndexes;
        // typedef Eigen::MatrixBlockIndexes<false, false> MatrixBlockIndexes;
        typedef Eigen::MatrixBlockView<matrix_t, Eigen::Dynamic, Eigen::Dynamic, false, false> MatrixBlockView;
        /// This function sets \f{ result = arg1 - arg0 \f}
        /// \note result may be of a different size than arg0 and arg1
        typedef boost::function<void (vectorIn_t arg0, vectorIn_t arg1, vectorOut_t result)> Difference_t;

        bool solve (vectorOut_t arg) const;

        bool isSatisfied (vectorIn_t arg) const;

        bool isSatisfied (vectorIn_t arg, vectorOut_t error) const;

        /// Returns true if the function was added, false otherwise
        /// A function can be added iif its outputs do not interfere with the
        /// output of another function.
        bool add (const DifferentiableFunctionPtr_t& f,
            const RowBlockIndexes& inArg,
            const RowBlockIndexes& outArg,
            const ColBlockIndexes& inDer,
            const RowBlockIndexes& outDer);

        /// \warning the two functions must have the same input and output
        /// indexes.
        bool replace (const DifferentiableFunctionPtr_t& oldf,
                      const DifferentiableFunctionPtr_t& newd);

        ExplicitSolver (const std::size_t& argSize, const std::size_t derSize)
          : argSize_ (argSize), derSize_ (derSize)
          ,  inArgs_ (),  inDers_ ()
          , outArgs_ (), outDers_ ()
          , argFunction_ (Eigen::VectorXi::Constant(argSize, -1))
          , derFunction_ (Eigen::VectorXi::Constant(derSize, -1))
          , squaredErrorThreshold_ (Eigen::NumTraits<value_type>::epsilon())
          // , Jg (derSize, derSize)
          , arg_ (argSize), diff_(derSize), diffSmall_()
        {
          inArgs_.addRow(0, argSize);
          inDers_.addCol(0, derSize);
        }

        /// Set error threshold
        void errorThreshold (const value_type& threshold)
        {
          squaredErrorThreshold_ = threshold * threshold;
        }
        /// Get error threshold
        value_type errorThreshold () const
        {
          return sqrt (squaredErrorThreshold_);
        }
        /// Get error threshold
        value_type squaredErrorThreshold () const
        {
          return squaredErrorThreshold_;
        }

        /// The set of variable indexes which are not affected by the
        /// resolution.
        const RowBlockIndexes& inArgs () const
        {
          return inArgs_;
        }

        /// The set of derivative variable indexes which are not affected by the
        /// resolution.
        const ColBlockIndexes& inDers () const
        {
          return inDers_;
        }

        /// Configuration parameters involved in the constraint resolution.
        ColBlockIndexes activeParameters () const;

        /// Velocity parameters involved in the constraint resolution.
        ColBlockIndexes activeDerivativeParameters () const;

        /// The set of variable indexes which are computed.
        const RowBlockIndexes& outArgs () const
        {
          return outArgs_;
        }

        /// The set of derivative variable indexes which are computed.
        const RowBlockIndexes& outDers () const
        {
          return outDers_;
        }

        /// The number of variables
        const std::size_t& argSize () const
        {
          return argSize_;
        }

        /// The number of derivative variables
        const std::size_t& derSize () const
        {
          return derSize_;
        }

        inline MatrixBlockView viewJacobian(matrix_t& jacobian) const
        {
          return MatrixBlockView(jacobian,
              outDers_.nbIndexes(), outDers_.indexes(),
              inDers_.nbIndexes(), inDers_.indexes());
        }

        /// Set the integration function
        void difference (const Difference_t& difference)
        {
          difference_ = difference;
        }

        /// Get the integration function
        const Difference_t& difference () const
        {
          return difference_;
        }

        // /// \param jacobian must be of dimensions (derSize - freeDers().nbIndexes(), freeDers().nbIndexes())
        /// \param jacobian must be of dimensions (derSize, derSize) but only a subsegment will be used.
        /// \warning it is assumed solve(arg) has been called before.
        void jacobian(matrixOut_t jacobian, vectorIn_t arg) const;

      private:
        typedef std::vector<bool> Computed_t;

        void computeFunction(const std::size_t& i, vectorOut_t arg) const;
        void computeJacobian(const std::size_t& i, matrixOut_t J) const;
        void computeOrder(const std::size_t& iF, std::size_t& iOrder, Computed_t& computed);

        const std::size_t argSize_, derSize_;

        struct Function {
          Function (DifferentiableFunctionPtr_t _f, RowBlockIndexes ia, RowBlockIndexes oa, ColBlockIndexes id, RowBlockIndexes od)
            : f (_f), inArg (ia), outArg (oa), inDer (id), outDer (od)
          {
            value   .resize(_f->outputSize());
            jacobian.resize(_f->outputDerivativeSize(), _f->inputDerivativeSize());
          }
          DifferentiableFunctionPtr_t f;
          RowBlockIndexes inArg, outArg;
          ColBlockIndexes inDer;
          RowBlockIndexes outDer;

          mutable vector_t value;
          mutable matrix_t jacobian;
        };

        RowBlockIndexes inArgs_;
        ColBlockIndexes inDers_;
        RowBlockIndexes outArgs_, outDers_;

        std::vector<Function> functions_;
        std::vector<std::size_t> computationOrder_;
        /// For dof i, dofFunction_[i] is the index of the function that computes it.
        /// -1 means it is the output of no function.
        Eigen::VectorXi argFunction_, derFunction_;
        Difference_t difference_;
        value_type squaredErrorThreshold_;
        // mutable matrix_t Jg;
        mutable vector_t arg_, diff_, diffSmall_;
    }; // class ExplicitSolver
    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_EXPLICIT_SOLVER_HH
