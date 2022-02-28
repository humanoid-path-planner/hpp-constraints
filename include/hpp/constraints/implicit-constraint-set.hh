//
// Copyright (c) 2016 - 2018 CNRS
// Authors: Joseph Mirabel, Florent Lamiraux
//

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

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
            comparison_.push_back (comp[i]);
          }
	  // Handle mask
	  mask_.insert(mask_.end(), constraint->mask_.begin(),
		       constraint->mask_.end());
	  // Recompute active rows
	  computeActiveRows();
	  computeIndices();
          // Resize temporary variables
          output_ = LiegroupElement(functions->outputSpace());
          logOutput_.resize(functions->outputSpace()->nv());
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

        HPP_SERIALIZABLE();
    }; // class ImplicitConstraintSet
    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_IMPLICIT_CONSTRAINT_SET_HH
