/*
    Simple Simplex Solver Class
    Copyright (C) 2012  Tamas Bolner
	For more information, visit: http://blog.bolner.hu/2012/08/22/simplex/
	
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    updated by Joseph Mirabel on August 2015.
*/

#ifndef HPP_CONSTRAINTS_SIMPLEX_SOLVER_HH
#define HPP_CONSTRAINTS_SIMPLEX_SOLVER_HH
#include <Eigen/Dense>

namespace hpp {
  namespace constraints {
    class SimplexSolver {
      public:
        enum Status {
          NOT_SOLVED,
          INFEASIBLE,
          UNBOUNDED,
          SOLVED
        };

        /// Multiple query solver.
        SimplexSolver(const size_t nVariables,
            const size_t nConstraints);

        Eigen::Ref <Eigen::MatrixXd> constraints ();
        Eigen::Ref <Eigen::VectorXd> rhs ();
        Status checkFeasability ();
        /// Build the solution found by checkFeasability.
        /// \note checkFeasability must have been called first.
        Status buildSolution ();

        /// Single-query solver.
        SimplexSolver(const Eigen::MatrixXd &constraints);

        Status status ();
        bool hasSolution();
        double getOptimum();
        Eigen::VectorXd getSolution();
        Eigen::VectorXd getSlack();

        static void buildSignMatrix (const Eigen::Ref <Eigen::VectorXd> rhs,
            Eigen::Ref<Eigen::MatrixXd> sign);

      protected:
        void init ();

      private:
        Eigen::MatrixXd tableau;
        Eigen::VectorXi selectedVariables_;
        Status foundSolution;
        double optimum;
        Eigen::VectorXd solution;
        int numberOfVariables,
            numberOfConstraints;

        int findPivot_min(int column);
        Status simplexAlgorithm(int variableNum);
        int getPivotRow(int column);
    };

  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_SIMPLEX_SOLVER_HH
