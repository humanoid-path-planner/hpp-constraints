//
// Copyright (c) 2016 CNRS
// Authors: Joseph Mirabel
//
// This file is part of hpp-constraints.
// hpp-constraints is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-constraints is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-constraints. If not, see
// <http://www.gnu.org/licenses/>.

#ifndef HPP_CONSTRAINTS_DIFFERENTIABLE_FUNCTION_SET_HH
# define HPP_CONSTRAINTS_DIFFERENTIABLE_FUNCTION_SET_HH

# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/differentiable-function.hh>

namespace hpp {
  namespace constraints {
    /// \addtogroup constraints
    /// \{

    /// Set of differentiable functions
    ///
    /// This class also handles selection of cols of the output matrix.
    class HPP_CONSTRAINTS_DLLAPI DifferentiableFunctionSet :
      public DifferentiableFunction
    {
      public:
        typedef std::vector<DifferentiableFunctionPtr_t> Functions_t;

        /// Return a shared pointer to a new instance
        ///
        /// \param name the name of the constraints,
        static DifferentiableFunctionSetPtr_t create (const std::string& name)
        {
          return DifferentiableFunctionSetPtr_t
            (new DifferentiableFunctionSet(name));
        }

        virtual ~DifferentiableFunctionSet () throw () {}

        /// \name Function stack management
        /// \{

        /// Get the stack of functions
        const Functions_t& functions () const
        {
          return functions_;
        }

        void add (const DifferentiableFunctionPtr_t& func)
        {
          if (functions_.empty()) {
            inputSize_           = func->inputSize();
            inputDerivativeSize_ = func->inputDerivativeSize();
            activeParameters_ = func->activeParameters();
            activeDerivativeParameters_ = func->activeDerivativeParameters();
          } else {
            assert (inputSize_           == func->inputSize());
            assert (inputDerivativeSize_ == func->inputDerivativeSize());

            activeParameters_ =
              activeParameters_ || func->activeParameters();
            activeDerivativeParameters_ =
              activeDerivativeParameters_ || func->activeDerivativeParameters();
          }
          functions_.push_back(func);
          result_.push_back (LiegroupElement (func->outputSpace ()));
          *outputSpace_ *= func->outputSpace ();
        }

        /// The output columns selection of other is not taken into account.
        void merge (const DifferentiableFunctionSetPtr_t& other)
        {
          const Functions_t& functions = other->functions();
          for (Functions_t::const_iterator _f = functions.begin();
              _f != functions.end(); ++_f)
            add (*_f);
        }

        /// \}

        std::ostream& print (std::ostream& os) const;

        /// Constructor
        ///
        /// \param name the name of the constraints,
        DifferentiableFunctionSet (const std::string& name)
          : DifferentiableFunction (0, 0, 0, name)
        {}

        DifferentiableFunctionSet ()
          : DifferentiableFunction (0, 0, 0, "Stack")
        {}

      protected:
        void impl_compute (LiegroupElementRef result, ConfigurationIn_t arg)
          const throw ()
        {
          size_type row = 0;
          std::size_t i = 0;
          for (Functions_t::const_iterator _f = functions_.begin();
              _f != functions_.end(); ++_f) {
            const DifferentiableFunction& f = **_f;
            f.impl_compute(result_ [i], arg);
            result.vector ().segment(row, f.outputSize()) =
              result_ [i].vector ();
            row += f.outputSize(); ++i;
          }
        }
        void impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t arg) const throw ()
        {
          size_type row = 0;
          for (Functions_t::const_iterator _f = functions_.begin();
              _f != functions_.end(); ++_f) {
            const DifferentiableFunction& f = **_f;
            f.impl_jacobian(jacobian.middleRows(row, f.outputDerivativeSize()), arg);
            row += f.outputDerivativeSize();
          }
        }
      private:
        Functions_t functions_;
        mutable std::vector <LiegroupElement> result_;
    }; // class DifferentiableFunctionSet
    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_DIFFERENTIABLE_FUNCTION_SET_HH
