# RoME.jl

| Stable | Dev | Coverage | Docs |
|--------|-----|----------|------|
| [![version][rjl-ver-img]][rjl-releases] | [![CI][rjl-ci-dev-img]][rjl-ci-dev-url] | [![codecov.io][rjl-cov-img]][rjl-cov-url] | [![docs][cjl-docs-img]][cjl-docs-url] <br> [![][caesar-slack-badge]][caesar-slack] |


## Introduction

Robot Motion Estimate (RoME.jl) is part of the overall [Caesar.jl][cjl-url] and provides a set of graph variables, factors, and utility features for robotics-related navigation, tracking, and mapping (i.e. SLAM).  RoME.jl helps build front-ends using the [Multi-modal iSAM] backend solver which is implemented over at [IncrementalInference.jl][iif-url].  See [the related references of interest here](http://www.juliarobotics.org/Caesar.jl/latest/refs/literature/).  Most notably, this package provides common navigation-type variables and factors to be included in more general [DistributedFactorGraphs.jl](https://github.com/JuliaRobotics/DistributedFactorGraphs.jl) graph objects.

[NavAbility.io](http://www.navability.io) helps the with administration and support of the Caesar.jl community, please reach out for any additional information (info@navability.io) or via the caesarjl Slack badge-link above.
## Installation

You can directly install with:

```julia
using Pkg
Pkg.add("RoME")
```

If you are interested in a broader toolkit, which includes a visualizer and database interaction, please see [Caesar.jl][cjl-url].

## Examples

See the common Caesar.jl documenation for more details [![cjl-docs-img]][cjl-docs-url].  Further examples can be found in the examples and test folders.
## Consider Citing

Consider citing our work using the common reference at [Caesar.jl Citation with IncrementalInference.jl DOI](https://github.com/JuliaRobotics/Caesar.jl#contributors).  We are grateful for many, many contributions within the Julia package ecosystem -- see the [Juliahub.com](https://juliahub.com/ui/Packages/RoME/VVxXB) page for dependencies.

## Comments and Issues Welcome

Please don't hesitate to open issues or suggestions in line with [JuliaRobotics code of conduct](https://github.com/JuliaRobotics/administration/blob/master/code_of_conduct.md).  Find [the Gist here](https://gist.github.com/dehann/5f943d833f5fb06f4e00a2f4fb9f945a).

<!-- md variables duplicated in Caesar.jl README -->
[rjl-url]: http://www.github.com/JuliaRobotics/RoME.jl
[rjl-cov-img]: https://codecov.io/github/JuliaRobotics/RoME.jl/coverage.svg?branch=master
[rjl-cov-url]: https://codecov.io/github/JuliaRobotics/RoME.jl?branch=master
[rjl-ci-dev-img]: https://github.com/JuliaRobotics/RoME.jl/actions/workflows/ci.yml/badge.svg
[rjl-ci-dev-url]: https://github.com/JuliaRobotics/RoME.jl/actions/workflows/ci.yml
[rjl-ver-img]: https://juliahub.com/docs/RoME/version.svg
[rjl-milestones]: https://github.com/JuliaRobotics/RoME.jl/milestones
[rjl-releases]: https://github.com/JuliaRobotics/RoME.jl/releases
[rjl-juliahub]: https://juliahub.com/ui/Packages/RoME/VVxXB

[iif-url]: https://github.com/JuliaRobotics/IncrementalInference.jl

[cjl-url]: https://github.com/JuliaRobotics/Caesar.jl
[cjl-docs-img]: https://img.shields.io/badge/docs-latest-blue.svg
[cjl-docs-url]: http://juliarobotics.github.io/Caesar.jl/latest/
[caesar-slack-badge]: https://img.shields.io/badge/Caesarjl-Slack-green.svg?style=popout
[caesar-slack]: https://join.slack.com/t/caesarjl/shared_invite/zt-ucs06bwg-y2tEbddwX1vR18MASnOLsw
