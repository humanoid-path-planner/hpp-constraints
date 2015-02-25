namespace hpp {
  namespace constraints {
    /// \mainpage
    ///
    /// This package provides a representation of the concept of mathematical
    /// functions from a manifold, for instace the configuration space of a robot
    /// (hpp::model::Device) to a finite dimensional vector space.<BR>
    ///
    /// Some robots can be subject to constraints.
    /// \li closed kinematic chains,
    /// \li quasi-static equilibrium for legged robots
    ///
    /// are examples of such constraints. These constraints can be defined
    /// implicitely by an equation the left hand side of which is a \ref
    /// DifferentiableFunction "differentiable function" of the
    /// robot configuration.
    ///
    /// \defgroup constraints Constraints
  } // namespace constraints
} // namespace hpp
