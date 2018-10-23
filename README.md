# RoME.jl

[![Build Status](https://travis-ci.org/JuliaRobotics/RoME.jl.svg?branch=master)](https://travis-ci.org/JuliaRobotics/RoME.jl)
[![codecov.io](https://codecov.io/github/JuliaRobotics/RoME.jl/coverage.svg?branch=master)](https://codecov.io/github/JuliaRobotics/RoME.jl?branch=master)

Robot Motion Estimate: A set of functions for developing front-ends for SLAM in [Julia](www.julialang.org) which adds transform, visualization and convenience functions to the [Multi-modal iSAM](http://frc.ri.cmu.edu/~kaess/pub/Fourie16iros.pdf) backend solver. The back-end solver is implemented in [IncrementalInference.jl](https://github.com/JuliaRobotics/IncrementalInference.jl).

## Introduction

This package forms part of the [Caesar.jl](https://github.com/JuliaRobotics/Caesar.jl) robot state estimate toolkit, towards robust solutions in robot navigation and mapping, which includes visualization and database interaction features. Robot style wrapper function and front-end factor graph generation functions are provided. Plot based visualization of robot belief based navigation variables is provided.

## Example

<a href="https://vimeo.com/190052649" target="_blank"><img src="https://raw.githubusercontent.com/JuliaRobotics/IncrementalInference.jl/master/doc/images/mmisamvid01.gif" alt="IMAGE ALT TEXT HERE" width="480" border="10" /></a>

## Installation

You can directly install all the RoME and Multi-modal iSAM functionality in Julia 0.6+ with:

    julia> Pkg.add("RoME")

If you are interested in a broader toolkit, which includes a visualizer and database interaction, please see [Caesar.jl](https://github.com/dehann/Caesar.jl).

## Examples

There are a few use cases in the examples folder, including the Victoria Park dataset. The code was recently refactored and several examples are due to appear.

Feel free to create and issue to resolve problems or if something is unclear.

## Comments and Issues Welcome

Please don't hesitate to open issues or suggestions in line with [JuliaRobotics code of conduct](https://github.com/JuliaRobotics/administration/blob/master/code_of_conduct.md).  Find [the Gist here](https://gist.github.com/dehann/5f943d833f5fb06f4e00a2f4fb9f945a).
