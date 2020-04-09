/// \mainpage
/// \anchor hpp_constraints_documentation
///
/// \section hpp_constraints_concepts Mathematical Concepts
///
/// This package provides a representation of the following mathematical
/// concepts:
/// \li \link hpp::constraints::DifferentiableFunction differentiable
/// functions\endlink from a \link hpp::pinocchio::LiegroupSpace Lie
/// group\endlink, for instance the configuration space of a robot
/// (hpp::pinocchio::Device) to a another \link
/// hpp::pinocchio::LiegroupSpace Lie group\endlink;
/// \li \link hpp::constraints::Implicit implicit constraints\endlink
/// are composed of a left hand side that is a differentiable
/// function, a \link hpp::constraints::ComparisonTypes_t comparison
/// vector\endlink, and a right hand side;
/// \li \link hpp::constraints::Explicit explicit constraints\endlink are
/// numerical constraints that can be solved by computing some input
/// variables with respect to others.
///
/// \section hpp_constraints_solvers Numerical Solvers
///
/// The constraints defined in the previous section can be solved by some
/// numerical solvers.
/// \li \link hpp::constraints::solver::HierarchicalIterative
/// HierarchicalIterative \endlink solves a set of implicit constraints with
/// priority order following a Newton like iterative procedure.
/// \li \link hpp::constraints::solver::BySubstitution BySubstitution \endlink
/// derives from the former and substitutes output variables of explicit
/// constraints with their explicit expression, solving numerically over the
/// remaining variables.
/// \li Compatible explicit constraints are gathered into \link
/// hpp::constraints::ExplicitConstraintSet ExplicitConstraintSet\endlink and
/// used by BySubstitution solver.
///
/// One of the most common use cases of explicit constraint is the
/// \link hpp::constraints::explicit_::RelativePose relative pose\endlink
/// constraint where a free floating object is hold by a robot gripper. In
/// this case, the configuration variables of the object can be explicitely
/// expressed with respect to the robot configuration variables.
///
/// \defgroup constraints Constraints
/// \defgroup solvers Constraint solvers
/// \defgroup hpp_constraints_tools Tools
/// \ingroup constraints
