# RoME.jl NEWS

RoME.jl follows semver, with only a few case specific exceptions.  Please see repo's [Milestones](https://github.com/JuliaRobotics/RoME.jl/milestones?state=closed) page for a more complete list of changes.  This NEWS file lists select changes like to produce breaking changes downstream.  Note serious efforts are taken to have both breaking and smaller changes go through a proper deprecation and warning printout cycle, consistent with JuliaLang convention.

## v0.24

- Manifolds based inertial odometry (preintegration).  Replaces previous 2016-2016 generation `InertialPose3` variables and factors.
- Testing enhancements and fixes to restore tests for RoMEPlotting and Caesar Docs.
- Minor code updates for standardized usage of new PyCaesar.jl.

## v0.23

- Bug fixes and maintenance.
- Transfer work to using Manopt.jl as new solver, slowly replacing previous Optim.jl approach (see IncrementalInference.jl).
- Extensions (weakdeps) replacement for legacy `Requires.jl`.

## v0.21

- Add `SnoopPrecompile` and Julia v1.8 min compat.
- Remove `FactorMetadata` and `ConvPerThread` usage as per IIF v0.31.
- Move code and files to subfolders `services`, `entities`, and `legacy`.
- Added new features `homograph_to_coordinates` and `coordinates_to_homography`.

## v0.20

- Multiple factor upgrades and fixes to Manifolds.jl v0.8.11 and above.
- Restore parametric batch solving.
- Convert all Manifolds usage to `ArrayPartition` while removing `ProductRepr`.
- Various bug fixes.

## v0.19

- Various standardizations and code quality improvements.
- Expand parametric solving support for `SpecialEuclidean(2)` or `(3)`.

## v0.18

- Various enhancements and maintenance fixes.
- Better consolidation of factor serialization and clearing deprecations.
- Maintenance on FluxPose2Pose2 factors to restore tests.
- Fix various factor bugs in serialization and constructors.
- Simplify CI testing.

## v0.17

- Graph generator API changing to `generateGraph_ABC`.
- Factors that can default to field `.Z` for easier/better use of dispatch (#538).

## v0.15

- A new canonical generator's name is changed to `generateCanonicalFG_Honeycomb!` (#445), and instead keeping the previous but recent function name (#440) `Beehive` available for a different upcoming canonical graph generator.
- Adding new `generateCanonicalFG_Helix2D!` plus convenience wrappers `Slew` and `Spiral`.
- Deprecating `generateCanonicalFG_ZeroPose2` and replaced by `generateCanonicalFG_ZeroPose`, already defaulting to keyword `varType=Pose2`.

