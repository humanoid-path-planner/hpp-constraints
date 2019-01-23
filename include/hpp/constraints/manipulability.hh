// Copyright (c) 2018 CNRS
// Authors: Joseph Mirabel
//
// This file is part of hpp-constraints
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
// hpp-constraints  If not, see
// <http://www.gnu.org/licenses/>.

#ifndef HPP_CONSTRAINTS_MANIPULABILITY_HH
# define HPP_CONSTRAINTS_MANIPULABILITY_HH

# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/config.hh>
# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/matrix-view.hh>

namespace hpp {
  namespace constraints {
    HPP_PREDEF_CLASS (Manipulability);
    typedef boost::shared_ptr <Manipulability> ManipulabilityPtr_t;

    /// \addtogroup constraints
    /// \{

    /// Differentiable function
    class HPP_CONSTRAINTS_DLLAPI Manipulability : public DifferentiableFunction
    {
    public:
      virtual ~Manipulability () {}

      static ManipulabilityPtr_t create (DifferentiableFunctionPtr_t function,
          DevicePtr_t robot, std::string name)
      {
        return ManipulabilityPtr_t (new Manipulability (function, robot, name));
      }

    protected:
      /// \brief Concrete class constructor should call this constructor.
      ///
      /// \param function the function which must be analysed
      /// \param name function's name
      Manipulability (DifferentiableFunctionPtr_t function,
          DevicePtr_t robot, std::string name);

      void impl_compute (LiegroupElementRef result, vectorIn_t argument) const;

      void impl_jacobian (matrixOut_t jacobian, vectorIn_t arg) const;

    private:
      DifferentiableFunctionPtr_t function_;
      DevicePtr_t robot_;

      Eigen::ColBlockIndices cols_;

      mutable matrix_t J_, J_JT_;
    }; // class Manipulability
    /// \}
  } // namespace constraints
} // namespace hpp


#endif // HPP_CONSTRAINTS_MANIPULABILITY_HH
