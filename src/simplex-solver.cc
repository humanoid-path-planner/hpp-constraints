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

#include "hpp/constraints/simplex-solver.hh"

#include <stdexcept>
#include <Eigen/Dense>

namespace hpp {
  namespace constraints {
    /**
     * Constructor
     *
     * @param const size_t nVariables
     * @param const size_t nConstraints
     * @returns SimplexSolver
     */
    SimplexSolver::SimplexSolver(const size_t nVariables,
            const size_t nConstraints) :
      foundSolution (NOT_SOLVED), optimum (0),
      numberOfVariables (nVariables), numberOfConstraints (nConstraints)
    {
      init ();
    }

    Eigen::Ref <Eigen::MatrixXd> SimplexSolver::constraints ()
    {
      return tableau.block (1, 0, numberOfConstraints, numberOfVariables);
    }

    Eigen::Ref <Eigen::VectorXd> SimplexSolver::rhs ()
    {
      return tableau.col (numberOfConstraints + numberOfVariables)
        .segment (1, numberOfConstraints);
    }

    SimplexSolver::Status SimplexSolver::checkFeasability ()
    {
      foundSolution = simplexAlgorithm(numberOfVariables);
      return foundSolution;
    }

    SimplexSolver::Status SimplexSolver::buildSolution ()
    {
      Eigen::Ref <Eigen::VectorXd> b = tableau.col (tableau.cols()-1);

      // Maximize
      solution.setZero ();
      for (int i = 0; i < numberOfConstraints; i++) {
        int k = selectedVariables_ (i);
        int temp = getPivotRow(k);

        if (temp > 0) {
          // Basic variable
          solution(k) = b (temp);
        } else {
          // Non-basic variable
          // This should not happen
          throw std::logic_error ("Basic variables should have a valid pivot row");
        }
      }
      if (!solution.segment (numberOfVariables, numberOfConstraints).isZero()) {
        // Problem is unfeasible.
        foundSolution = INFEASIBLE;
        //if (foundSolution != INFEASIBLE)
          //throw std::logic_error ("A slack variable is not zero and the status "
              //"is not INFEASIBLE");
      }
      return foundSolution;
    }

    /**
     * Constructor
     * 
     * @param const MatrixXd &constraints Full matrix for the constraints.
     *                                    Contains also the righthand-side values.
     * @returns SimplexSolver
     */
    SimplexSolver::SimplexSolver(const Eigen::MatrixXd &constraints) :
      foundSolution (NOT_SOLVED), optimum (0),
      numberOfVariables (constraints.cols () - 1),
      numberOfConstraints (constraints.rows())
    {
      init ();

      // The bottom rows represent the constraints
      this->constraints() = constraints.leftCols(numberOfVariables);
      this->rhs ()        = constraints.col (numberOfVariables);

      /*
         Simplex algorithm
         */
      checkFeasability ();

      /*
         Fetch solution
         */
      buildSolution ();
    }

    /**
     * Returns Status corresponding to the current status
     *
     * @returns Status
     */
    SimplexSolver::Status SimplexSolver::status () {
      return foundSolution;
    }

    /**
     * Returns true if a solution has been found.
     * Return false otherwise.
     *
     * @returns boolean
     */
    bool SimplexSolver::hasSolution() {
      return foundSolution == SOLVED || foundSolution == UNBOUNDED;
    }

    /**
     * Returns the maximum/minimum value of
     * the objective function.
     * 
     * @returns double
     */
    double SimplexSolver::getOptimum() {
      return this->optimum;
    }

    /**
     * Returns a vector with the variable
     * values for the solution.
     *
     * return VectorXd
     */
    Eigen::VectorXd SimplexSolver::getSolution() {
      return solution.segment (0, numberOfVariables);
    }

    /**
     * Returns a vector with the slack variable values.
     *
     * return VectorXd
     */
    Eigen::VectorXd SimplexSolver::getSlack() {
      return solution.segment (numberOfVariables, numberOfConstraints);
    }

    void SimplexSolver::init () {
      /*
         Validate input parameters
         */
      if (numberOfConstraints < 1) {
        throw(std::invalid_argument("SimplexSolver: The constraint "
              "matrix must contain at least one row."));
      }

      if (numberOfVariables < 1) {
        throw(std::invalid_argument("SimplexSolver: The constraint "
              "matrix must contain at least two columns."));
      }

      selectedVariables_.resize (numberOfConstraints);
      for (int i = 0; i < numberOfConstraints; ++i)
        selectedVariables_(i) = i + numberOfVariables;

      /*
         Initialize the tableau
         */
      // Maximize
      tableau.resize(numberOfConstraints + 1,
          numberOfVariables + numberOfConstraints + 1);

      // The first row is the objective function
      tableau.row (0).rightCols (numberOfConstraints + 1).setZero ();
      tableau.row (0).leftCols (numberOfVariables).setConstant (-1);

      // The bottom rows represent the constraints
      tableau.bottomRows (numberOfConstraints)
        .leftCols (numberOfVariables).setZero ();
      tableau.bottomRows (numberOfConstraints)
        .middleCols (numberOfVariables, numberOfConstraints).setIdentity ();
      tableau.bottomRows (numberOfConstraints)
        .col (numberOfVariables + numberOfConstraints).setZero ();

      solution.resize(numberOfVariables + numberOfConstraints);
      solution.setZero ();
    }

    /**
     * Searches for the pivot row in the given column, by calculating the ratios.
     * Tries to find smallest non-negative ratio.
     * Returns -1 if all possible pivots are 0 or if the ratios are negative.
     * Deals with cases like this:  0/negative < 0/positive
     * 
     * @param int column
     * @returns int Returns the number of the pivot row, or -1 if found none.
     */
    int SimplexSolver::findPivot_min(int column) {
      int minIndex = -1;
      int constantColumn = this->tableau.cols() - 1;
      double minRatio = 0;
      double minConstant = 0;	// For the "0/negative < 0/positive" difference the constants have to be tracked also.
      double ratio;
      int rowNum = this->tableau.rows();

      for (int i = 1; i < rowNum; i++) {
        if (std::abs (this->tableau(i, column))
            <= Eigen::NumTraits<double>::dummy_precision()) {
          continue;
        }

        ratio = this->tableau(i, constantColumn) / this->tableau(i, column);
        if (ratio < 0) {
          //The ratio must be non-negative
          continue;
        }

        if (minIndex == -1) {
          // First pivot candidate
          minIndex = i;
          minRatio = ratio;
          minConstant = this->tableau(i, constantColumn);
        } else {
          if (ratio == 0 && ratio == minRatio) {
            // 0/negative < 0/positive
            if (this->tableau(i, constantColumn) < minConstant) {
              minIndex = i;
              minRatio = ratio;
              minConstant = this->tableau(i, constantColumn);
            }
          } else if (ratio < minRatio) {
            minIndex = i;
            minRatio = ratio;
            minConstant = this->tableau(i, constantColumn);
          }
        }
      }

      return minIndex;
    }

    /**
     * Iterates through the this->tableau matrix to solve the problem.
     * 
     * @param int variableNum The number of variables (dimensions). (different for the minimization problem)
     * @returns bool Returns true if a solution has been found. Returns false otherwise.
     */
    SimplexSolver::Status SimplexSolver::simplexAlgorithm(int variableNum) {
      Eigen::MatrixXd::Index pivotColumn;
      int pivotRow;
      double minC;
      bool hasNegativeElmt;

      /// Check the tableau:
      for (Eigen::MatrixXd::Index i = 0; i < numberOfVariables; ++i) {
        // When all the coefficients of a variable is 0, we consider the
        // variable inactive.
        // TODO: if the tableau(i,0) is non-negative, then the problem is unbounded.
        // Otherwise, the corresponding variable is zero.
        if (tableau.bottomRows (numberOfConstraints).col (i).isZero ()) {
          tableau.col (i).setZero ();
        }
      }

      while (true) {
        pivotColumn = -1;
        minC = 0;
        hasNegativeElmt = false;
        /*
           Find pivot column, check for halt condition
           */
        for (Eigen::MatrixXd::Index i = 0; i < variableNum; ++i) {
          if (tableau (0, i) < minC) {
            hasNegativeElmt = true;
            for (Eigen::MatrixXd::Index j = 1; j < numberOfConstraints + 1; ++j)
              if (tableau(j,i) > Eigen::NumTraits<double>::dummy_precision()) {
                minC = tableau (0, i);
                pivotColumn = i;
                break;
              }
          }
        }
        if (pivotColumn == -1) {
          /// TODO Check if the last variables are 0.
          /// If not, the problem is INFEASIBLE
          if (hasNegativeElmt) {
            /// Unbounded problem.
            return UNBOUNDED;
          }
          break;
        }
        /* Find pivot row */
        pivotRow = this->findPivot_min(pivotColumn);
        if (pivotRow == -1) {
          // No solution
          return INFEASIBLE;
        }
        /*
        this->tableau.row(0).leftCols(variableNum).minCoeff(&pivotColumn);
        if (this->tableau(0, pivotColumn) >= 0) {
          //Found no negative coefficient
          break;
        }
        pivotRow = this->findPivot_min(pivotColumn);
        if (pivotRow == -1) {
          // No solution
          return false;
        }// */

        /*
           Do pivot operation
           */
        selectedVariables_(pivotRow - 1) = pivotColumn;
        this->tableau.row(pivotRow) /= this->tableau(pivotRow, pivotColumn);
        this->tableau(pivotRow, pivotColumn) = 1; // For possible precision issues
        for (int i = 0; i < this->tableau.rows(); i++) {
          if (i == pivotRow) continue;

          this->tableau.row(i) -= this->tableau.row(pivotRow) * this->tableau(i, pivotColumn);
          this->tableau(i, pivotColumn) = 0; // For possible precision issues
        }
      }

      return SOLVED;
    }

    /**
     * If the given column has only one coefficient with value 1 (except in topmost row), and all other
     * coefficients are zero, then returns the row of the non-zero value.
     * Otherwise return -1.
     * This method is used in the final step of maximization, when we read
     * the solution from the tableau.
     * 
     * @param int column
     * @returns int
     */
    int SimplexSolver::getPivotRow(int column) {
      int one_row = -1;

      for (int i = 1; i < this->tableau.rows(); i++) {
        if (this->tableau(i, column) == 1) {
          if (one_row >= 0) {
            return -1;
          } else {
            one_row = i;
            continue;
          }
        } else if (this->tableau(i, column) != 0) {
          return -1;
        }
      }

      return one_row;
    }

    void SimplexSolver::buildSignMatrix (const Eigen::Ref <Eigen::VectorXd> rhs,
        Eigen::Ref<Eigen::MatrixXd> sign)
    {
      if (sign.cols () != rhs.rows ()) {
        throw(std::invalid_argument("SimplexSolver: Sign matrix must be a "
              "diagonal matrix with as many elements as in rhs."));
      }
      sign.diagonal () = (2 * (rhs.array() > Eigen::ArrayXd::Zero (rhs.rows())).cast <int>() - 1).cast <double>();
    }
  } // namespace constraints
} // namespace hpp
