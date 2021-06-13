# RoME.jl NEWS

RoME.jl follows semver, with only a few case specific exceptions.  Please see repo's [Milestones](https://github.com/JuliaRobotics/RoME.jl/milestones?state=closed) page for a more complete list of changes.  This NEWS file lists select changes like to produce breaking changes downstream.  Note serious efforts are taken to have both breaking and smaller changes go through a proper deprecation and warning printout cycle, consistent with JuliaLang convention.

## v0.15.1 -> v0.15.2 (New graph generators)

- A new canonical generator's name is changed to `generateCanonicalFG_Honeycomb!` (#445), and instead keeping the previous but recent function name (#440) `Beehive` available for a different upcoming canonical graph generator.
- Adding new `generateCanonicalFG_Helix2D!` plus convenience wrappers `Slew` and `Spiral`.
