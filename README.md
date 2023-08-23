# Ephemerides.jl

_A Modern Binary Ephemeris Reader for Julia, in Julia_

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaspacemissiondesign.github.io/Ephemerides.jl/stable/) 
[![Dev Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaspacemissiondesign.github.io/Ephemerides.jl/dev/) 
[![Build Status](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaSpaceMissionDesign/Ephemerides.jl/branch/main/graph/badge.svg?token=3SJCV229XX)](https://codecov.io/gh/JuliaSpaceMissionDesign/Ephemerides.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

Ephemerides.jl is a Julia library that provides fast, thread-safe and allocation-free access to binary JPL
[SPK](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html) and [PCK](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html) ephemeris files. Completely written in Julia, it enables 
Automatic-Differentiation (AD) via [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) 
across all of its function calls.

This package is meant to be used in combination with [FrameTransformations.jl](https://github.com/JuliaSpaceMissionDesign/FrameTransformations.jl), which enables transformations between different point and axes. Indeed, differently from traditional ephemeris readers such as [CALCEPH](https://www.imcce.fr/inpop/calceph) and [SPICE](https://naif.jpl.nasa.gov/naif/toolkit.html), this package is only meant to read the data stored in the binary kernels. It does not perform transformations between reference frames nor concatanations of state vectors. For example, if ephemeris data for point 399 (Earth) is defined with respect to point 3 (Earth-Moon Barycenter) in the ICRF axes, with this package we will only be able to compute the state vector from 399 to 3 or viceversa. 

## Installation

## Quickstart

## Current Limitations
- The supported JPL binary SPK/PCK segments are limited to: type 1, 2, 3, 8, 12, 21. 
- Binary INPOP kernels are not supported. 

## Documentation 
For further information on this package please refer to the [stable documentation](https://juliaspacemissiondesign.github.io/Ephemerides.jl/stable/)

## Support
If you found this package useful, please consider starring the repository. We also encourage 
you to take a look at other astrodynamical packages of the [JSMD](https://github.com/JuliaSpaceMissionDesign/) organisation.