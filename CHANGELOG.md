# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

## [9.0.2] - 2026-07-24



## [9.0.0] - 2026-07-07



## [7.0.0] - 2026-03-06



## [6.1.0] - 2025-10-23



## [6.0.0] - 2024-12-07

Changes in v6.0.0
- hpp-fcl dependency has been replaced by coal
- updates for coal v3


## [5.2.0] - 2024-10-09

Changes in v5.2.0:
- nix: move package to nixpkgs
- ci: use https
- setup mergify


## [5.1.0] - 2024-07-02

Changes in v5.1.0:
- fix assert in Eigen
- remove use of deprecated symbols
- Nix: initial support
- update tooling


## [5.0.0] - 2024-03-31

Changes in v5.0.0:
- :warning: removed deprecated symbols
- update to hpp-pinocchio v5.0.0
- make RelativeCom use thread-safe device
- fix compilation warnings
- update packaging
- update tooling


## [4.15.1] - 2023-01-20



## [4.14.0] - 2022-11-02



## [4.13.0] - 2022-05-31



## [4.12.0] - 2021-10-06

Changes in v4.12.0:
- Fix ConvexShapeContactComplement


## [4.11.0] - 2021-05-04

Changes in v4.11.0:
- Solvers now handle constraints with right hand sides in Lie groups.
- A mask has been added to class Implicit to select which lines of the
  constraint should be taken into account (active rows).
- Use hpp::shared_ptr instead of boost::shared_ptr
- Remove usage of boost list_of
- [solver] Add methods isSatisfied with input error threshold.


## [4.10.1] - 2020-09-24

Changes since v4.9.0:
* ConvexShapeContact classes have been improved.
  - stable position of objects is now unique for any right hand side value of
    the complement constraint.
* Use cmake to handle dependencies (cmake submodule).
* Architecture of constraint classes has been simplified.
* Use Romeo instead of Baxter in tests.
* Enable users to use GenericTransformation with a joint2 equal to 0x0.
* Add serialization functions.

## [4.10.0] - 2020-08-17



## [4.9.1] - 2020-05-14

Fix unit tests.

## [4.9.0] - 2020-04-29

Changes in v4.9.0:
- Computation of constraint right hand side have been fixed
  - in ExplicitConstraintSet, the computation was wrong, it has been fixed,
  - in Implicit, the computation has been made similar to the one in
    HierarchicalIterative.
- Template flag definitions in GenericTransformation have been fixed
  - to comply with C++11 standard.
- throw declaration have been removed to comply with C++11 standard.
- Configuration variable RUN_TESTS has benn replaced by BUILD_TESTING
  - for homogeneity with other hpp packages.
- CMake Exports

## [4.8.0] - 2019-11-28

Changes since v4.7.0:
- declare pinocchio dependency
- update to changes in pinocchio
- Fix indentation in output streams
- update CMake

## [4.7.0] - 2019-10-04

Changes since v4.6.0:
- include pinocchio before boost
- Update symbolic calculus class ScalarProduct + add GJK executable.
- Do not create an implicit constraint when not necessary.
- explicit cast template arguments
- [HierarchicalIterative] Add a method testing inclusion of manifolds
- [explicit::RelativePose] Fix bug in implicit to explicit rhs conversion.


## [4.5.0] - 2019-04-24

Changes since v4.4.0:
- fix MatrixView for eigen 3.3.4
- fix compilation warnings
- Add a test on right hand side of solvers.
- s/BOOST_MESSAGE/BOOST_TEST_MESSAGE


## [4.4.0] - 2019-03-19

Changes since v4.3.0:
- Update Eigen minimum version.
- Update to Pinocchio v2 + varying right hand side.


## [4.3.0] - 2019-01-31

- Fix implicit constraint set copy constructor
- Add HierarchicalIterative and BySubstitution::getRightHandSide
- [CMake] add required dependency on romeo_description if we run tests
- [HierarchicalIterative] Add accessor to free variables
- Minor fixes


## [4.2.0] - 2018-10-11

Changes since v4.1:
- * Refactor ExplicitConstraintSet: rename members and methods to better fit RSS paper notation and pinocchio convention:
  - argSize -> nq,
  - derSize -> nv,
  - freeArgs -> notOutArgs_,
  - freeDers -> notOutDers_,
  - viewJacobian -> jacobianNotOutToOut.
* Constructor and create methods of Explicit take a LiegroupSpace instead of a robot,
* Modify prototype of BySubstitution::projectVectorOnKernel
  - vectorIn_t -> ConfigurationIn_t,
  - vectorOut_t -> ConfigurationOut_t.
* In class HierarchicalIterative, rename "reduction" -> "free variables".
- [Documentation] Add documentation for class HierarchicalIterativeSolver.

## [3.2] - 2017-03-17



## [4.0] - 2018-03-14

* Update dependency versions.
* Update to changes in hpp-pinocchio
* Fix compilation warning
* Clean test-jacobians
* Fix test
* Export whether qpOASES is used in .pc file.
* Add a cache variable to avoid running unit tests.
* Make dependency to qpOASES optional.
* Add GenericTransformation::print
* Make GenericTransformation::joint[12] const method
* Update to changes in hpp-pinocchio and pinocchio
* Add possibility to check the jacobians.
* Fix finite difference (do not saturate when integrating)
* Fix compilation with numerical output activated
* Fix test-jacobians.cc
* Rewrite finite difference algorithm
* Make finite difference algorithm public function of DifferentiableFunction
* Add diagnostic for a test that fails
* Clean GenericTransformation
* [WIP] Update GenericTransformation
* Remove dep from tools.hh to macros.hh
* [WIP] Check computeLog and computeJlog with impl in GenericTransformation
* Fix computeLog and computeJlog
* Use JointConstPtr_t instead of JointPtr_t when relevant
* Add missing header
* Fix old hpp/_constraints/orientation.hh
* Add option to CMakeLists and update doc
* Fix compilation due to changes in pinocchio
* Fix DiffentiableFunction
* Add a constructor to class DistanceBetweenBodies
* Clean code
* Add StaticStability and QPStaticStability
* Fix compilation warning
* Fix DistanceBetweenBodies
* Fix DistanceBetweenBodies and update to API changes of Device
* Fix compilation warning
* Update to change in hpp-pinocchio
* Add ConfigurationConstraint, ConvexShape and ConvexShapeContact
* Add comment in p/m comparison of com related functions
* Add DistanceBetweenPointsInBodies
* Add SymbolicFunction and update COM related p/m comparison
* Add p/m comparison for DistanceBetweenBodies
* Clean pmdiff/tvalue.cc
* Add DistanceBetweenBodies
* Uncommented p/m comparison tests
* Add ComBetweenFeet
* Fix symbolic calculus
* Clean pinocchio/model comparison test
* Add RelativeCom
* Add pinocchio version of GenericTransformation
* Move function using hpp-model to _constraints folder (hh and cc files)


[Unreleased]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v9.0.2...HEAD
[9.0.2]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v9.0.0...v9.0.2
[9.0.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v7.0.0...v9.0.0
[7.0.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v6.1.0...v7.0.0
[6.1.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v6.0.0...v6.1.0
[6.0.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v5.2.0...v6.0.0
[5.2.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v5.1.0...v5.2.0
[5.1.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v5.0.0...v5.1.0
[5.0.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.15.1...v5.0.0
[4.15.1]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.14.0...v4.15.1
[4.14.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.13.0...v4.14.0
[4.13.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.12.0...v4.13.0
[4.12.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.11.0...v4.12.0
[4.11.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.10.1...v4.11.0
[4.10.1]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.10.0...v4.10.1
[4.10.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.9.1...v4.10.0
[4.9.1]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.9.0...v4.9.1
[4.9.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.8.0...v4.9.0
[4.8.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.7.0...v4.8.0
[4.7.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.5.0...v4.7.0
[4.5.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.4.0...v4.5.0
[4.4.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.3.0...v4.4.0
[4.3.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.2.0...v4.3.0
[4.2.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v4.0...v4.2.0
[4.0]: https://github.com/humanoid-path-planner/hpp-constraints/compare/v3.2...v4.0
[3.2]: https://github.com/humanoid-path-planner/hpp-constraints/releases/tag/v3.2
