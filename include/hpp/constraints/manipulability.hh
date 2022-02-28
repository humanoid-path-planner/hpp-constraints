// Copyright (c) 2018 CNRS
// Authors: Joseph Mirabel
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

#ifndef HPP_CONSTRAINTS_MANIPULABILITY_HH
# define HPP_CONSTRAINTS_MANIPULABILITY_HH

# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/config.hh>
# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/matrix-view.hh>

namespace hpp {
  namespace constraints {
    HPP_PREDEF_CLASS (Manipulability);
    typedef shared_ptr <Manipulability> ManipulabilityPtr_t;

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

      bool isEqual(const DifferentiableFunction& other) const {
        const Manipulability& castother = dynamic_cast<const Manipulability&>(other);
        if (!DifferentiableFunction::isEqual(other))
          return false;
        
        if (function_ != castother.function_)
          return false;
        if (robot_ != castother.robot_)
          return false;
        if (cols_.cols() != castother.cols_.cols())
          return false;
        if (J_ != castother.J_)
          return false;
        if (J_JT_ != castother.J_JT_)
          return false;
        
        return true;
      }

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
