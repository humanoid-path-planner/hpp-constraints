# hpp-constraints

[![Pipeline status](https://gitlab.laas.fr/humanoid-path-planner/hpp-constraints/badges/master/pipeline.svg)](https://gitlab.laas.fr/humanoid-path-planner/hpp-constraints/commits/master)
[![Coverage report](https://gitlab.laas.fr/humanoid-path-planner/hpp-constraints/badges/master/coverage.svg?job=doc-coverage)](https://gepettoweb.laas.fr/doc/humanoid-path-planner/hpp-constraints/master/coverage/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/humanoid-path-planner/hpp-constraints/master.svg)](https://results.pre-commit.ci/latest/github/humanoid-path-planner/hpp-constraints)

Numerical constraints and solvers of the
[Humanoid Path Planner](https://humanoid-path-planner.github.io/hpp-doc/).

- Differentiable functions from a Lie group, for instance the configuration
  space of a robot, to another Lie group: position and orientation of frames
  (`GenericTransformation`), center of mass (`RelativeCom`), distance between
  bodies, contact between convex shapes, static stability...
- Implicit constraints, defined by a differentiable function, comparison types
  (equality or inequality) and a right hand side.
- Explicit constraints, where some configuration variables are computed from
  the others, for instance the pose of an object held by a gripper
  (`explicit_::RelativePose`) or locked joints.
- Solvers: `solver::HierarchicalIterative` solves implicit constraints by
  order of priority with a Newton-like method, and `solver::BySubstitution`
  also substitutes the variables computed by explicit constraints.

The Python bindings are in module `pyhpp.constraints` of
[hpp-python](https://github.com/humanoid-path-planner/hpp-python).

## Documentation

- [C++ API](https://gepetto.github.io/doc/hpp-constraints/doxygen-html/index.html)
- [Python API](https://humanoid-path-planner.github.io/hpp-doc/reference/hpp-python/api/constraints.html)
- [Installation](https://humanoid-path-planner.github.io/hpp-doc/installation/installation.html)
