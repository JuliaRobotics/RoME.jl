# RoME.jl

| Stable v0.14 | Stable v0.15 | Dev | Coverage | Docs |
|--------------|--------------|-----|----------|------|
| [![Build Status](https://travis-ci.org/JuliaRobotics/RoME.jl.svg?branch=release%2Fv0.14)](https://travis-ci.org/JuliaRobotics/RoME.jl) | [![Build Status](https://travis-ci.org/JuliaRobotics/RoME.jl.svg?branch=release%2Fv0.15)](https://travis-ci.org/JuliaRobotics/RoME.jl) | [![CI](https://github.com/JuliaRobotics/RoME.jl/workflows/CI/badge.svg?branch=master)](https://github.com/JuliaRobotics/RoME.jl/actions?query=workflow%3ACI) | [![codecov.io][r-cov-img]][r-cov-url] | [![docs][docs-shield]][caesar-docs] <br> [![][caesar-slack-badge]][caesar-slack] |

## Introduction

Robot Motion Estimate: A set of functions for developing robotics-related navigation, tracking, and mapping (i.e. SLAM) front-ends using the [Multi-modal iSAM] backend solver.  The backend is implemented more abstractly in the [IncrementalInference.jl package](https://github.com/JuliaRobotics/IncrementalInference.jl).  See [the related references of interest here](http://www.juliarobotics.org/Caesar.jl/latest/refs/literature/).  Most notably, this package provides common navigation-type variables and factors to be included in more general [DistributedFactorGraphs.jl](https://github.com/JuliaRobotics/DistributedFactorGraphs.jl) graph objects.

Please contact info@navability.io for further support on this package, [NavAbility](https://www.navability.io).

## Installation

You can directly install with:

```julia
using Pkg
Pkg.add("RoME")
```

If you are interested in a broader toolkit, which includes a visualizer and database interaction, please see [Caesar.jl](https://github.com/dehann/Caesar.jl).

## Consider Citing

We are grateful for many, many contributions within the Julia package ecosystem -- see the [`Caesar.jl/Project.toml`](https://github.com/JuliaRobotics/Caesar.jl/blob/master/Project.toml) (and upstream package) files for a far reaching list of contributions.

Consider citing our work using the common reference at [Caesar.jl Citation with IncrementalInference.jl DOI](https://github.com/JuliaRobotics/Caesar.jl#contributors).

## Examples

[![docs][docs-shield]][caesar-docs]
See project wide Caesar.jl documentation for more details (click on badge).

<a href="https://vimeo.com/190052649" target="_blank"><img src="https://raw.githubusercontent.com/JuliaRobotics/IncrementalInference.jl/master/doc/images/mmisamvid01.gif" alt="IMAGE ALT TEXT HERE" width="480" border="10" /></a>

## Comments and Issues Welcome

Please don't hesitate to open issues or suggestions in line with [JuliaRobotics code of conduct](https://github.com/JuliaRobotics/administration/blob/master/code_of_conduct.md).  Find [the Gist here](https://gist.github.com/dehann/5f943d833f5fb06f4e00a2f4fb9f945a).


[r-cov-img]: https://codecov.io/github/JuliaRobotics/RoME.jl/coverage.svg?branch=master
[r-cov-url]: https://codecov.io/github/JuliaRobotics/RoME.jl?branch=master
[r-build-img]: https://travis-ci.org/JuliaRobotics/RoME.jl.svg?branch=master
[r-build-stbl]: https://travis-ci.org/JuliaRobotics/RoME.jl.svg?branch=release%2Fv0.15
[r-build-url]: https://travis-ci.org/JuliaRobotics/RoME.jl

[docs-shield]: https://img.shields.io/badge/docs-latest-blue.svg
[caesar-docs]: http://juliarobotics.github.io/Caesar.jl/latest/
[caesar-slack-badge]: https://img.shields.io/badge/Caesarjl-Slack-green.svg?style=popout
[caesar-slack]: https://join.slack.com/t/caesarjl/shared_invite/zt-ucs06bwg-y2tEbddwX1vR18MASnOLsw
