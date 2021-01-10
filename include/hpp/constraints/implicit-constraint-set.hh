//
// Copyright (c) 2016 - 2018 CNRS
// Authors: Joseph Mirabel, Florent Lamiraux
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

#ifndef HPP_CONSTRAINTS_IMPLICIT_CONSTRAINT_SET_HH
# define HPP_CONSTRAINTS_IMPLICIT_CONSTRAINT_SET_HH

# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/implicit.hh>
# include <hpp/constraints/differentiable-function-set.hh>

namespace hpp {
  namespace constraints {
    /// \addtogroup constraints
    /// \{

    /// Set of implicit constraints
    ///
    /// This class also handles selection of cols of the output matrix.
    class HPP_CONSTRAINTS_DLLAPI ImplicitConstraintSet :
      public Implicit
    {
      public:
        typedef std::vector<ImplicitPtr_t> Implicits_t;

        /// Return a shared pointer to a new instance
        ///
        /// \param name the name of the constraints,
        static ImplicitConstraintSetPtr_t create (const std::string& name)
        {
          return ImplicitConstraintSetPtr_t
            (new ImplicitConstraintSet(name));
        }

        virtual ~ImplicitConstraintSet () {}

        /// \name Function stack management
        /// \{

        void add (const ImplicitPtr_t& constraint)
        {
          assert (HPP_DYNAMIC_PTR_CAST (DifferentiableFunctionSet,
                                        function_));
          DifferentiableFunctionSetPtr_t functions
            (HPP_STATIC_PTR_CAST (DifferentiableFunctionSet, function_));
          functions->add (constraint->functionPtr ());
          constraints_.push_back(constraint);
          // Handle comparison types
          const ComparisonTypes_t& comp (constraint->comparisonType ());
          for (std::size_t i = 0; i < comp.size(); ++i) {
            if (comp[i] == Equality)
            {
              equalityIndices_.addRow(comparison_.size(), 1);
            }
            comparison_.push_back (comp[i]);
          }
          equalityIndices_.updateRows<true, true, true>();
	  // Handle mask
	  mask_.insert(mask_.end(), constraint->mask_.begin(),
		       constraint->mask_.end());
	  // Recompute active rows
	  computeActiveRows();
        }

        /// Get indices of constraint coordinates that are equality
        const Eigen::RowBlockIndices& equalityIndices () const
        {
          return equalityIndices_;
        }
        /// Get constraints
        const Implicits_t& constraints () const
        {
          return constraints_;
        }

        /// The output columns selection of other is not taken into account.
        void merge (const ImplicitConstraintSetPtr_t& other)
        {
          const Implicits_t& constraints = other->constraints();
          for (Implicits_t::const_iterator constraint = constraints.begin();
              constraint != constraints.end(); ++constraint)
            add (*constraint);
        }

        /// \}

        std::ostream& print (std::ostream& os) const
        {
          function_->print (os);
          return os;
        }

        /// Constructor
        ///
        /// \param name the name of the constraints,
        ImplicitConstraintSet (const std::string& name)
          : Implicit (DifferentiableFunctionSet::create (name),
                      ComparisonTypes_t (), std::vector<bool>())
        {
        }

        ImplicitConstraintSet ()
          : Implicit (DifferentiableFunctionSet::create ("Stack"),
                      ComparisonTypes_t (), std::vector<bool>())
          {}

        ImplicitConstraintSet (const ImplicitConstraintSet& o)
          : Implicit (DifferentiableFunctionSet::create ("Stack"),
                      ComparisonTypes_t (), std::vector<bool>())
          {
            const Implicits_t& constraints = o.constraints();
            for (Implicits_t::const_iterator constraint = constraints.begin();
                constraint != constraints.end(); ++constraint)
              add (*constraint);
          }
        
      private:
        Implicits_t constraints_;
        Eigen::RowBlockIndices equalityIndices_;

        HPP_SERIALIZABLE();
    }; // class ImplicitConstraintSet
    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_IMPLICIT_CONSTRAINT_SET_HH
