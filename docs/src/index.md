# Welcome to Ephemerides.jl!

Ephemerides.jl is a Julia library that provides fast, thread-safe and allocation-free access to binary JPL
[SPK](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html) and [PCK](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html) ephemeris files. Completely written in Julia, it enables 
Automatic-Differentiation (AD) via [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) 
across all of its function calls. 

It outperforms both [SPICE.jl](https://github.com/JuliaAstro/SPICE.jl) and [CALCEPH.jl](https://github.com/JuliaAstro/CALCEPH.jl) calls for most types of SPK segments and supports state vector and orientation angles computation up to order 3 (jerk).

This package is meant to be used in combination with [FrameTransformations.jl](https://github.com/JuliaSpaceMissionDesign/FrameTransformations.jl), which enables transformations between different point and axes. Indeed, differently from traditional ephemeris readers such as [CALCEPH](https://www.imcce.fr/inpop/calceph) and [SPICE](https://naif.jpl.nasa.gov/naif/toolkit.html), this package is only meant to read the data stored in the binary kernels. It does not perform transformations between reference frames nor concatenations of state vectors. For example, if ephemeris data for point 399 (Earth) is defined with respect to point 3 (Earth-Moon Barycenter) in the ICRF axes, with this package we will only be able to compute the state vector from 399 to 3 or viceversa. 


## Installation

## Quickstart

Load SPK and PCK ephemeris kernels: 

```julia
using Ephemerides 

# Load a single SPK kernel 
eph_spk = EphemerisProvider("de440.bsp")

# Load a single PCK kernel
eph_pck = EphemerisProvider("pa440.bsp")

# Load multiple SPK and PCK kernels
eph = EphemerisProvider(["de440.bsp", "pa440.bsp"])
```

Inspect the kernels properties:
```julia
# Retrieve the list of NAIF ID for all the available points 
points = ephem_available_points(eph)

# Retrieve the list of NAIF ID for all the available axes
axes = ephem_available_axes(eph)
```

Retrieve state and orientation data:
```julia

# TDB seconds at 2000-01-01T12:00:00 (J2000)
time = 0.0

# Compute the position of point 399 with respect to 3 at J2000
pos = ephem_vector3(eph, 3, 399, time)

# Compute the position and its derivatives for point 299 with respect to 2
pvaj = ephem_vector12(eph, 2, 299, time)

# Compute the orientation of axes 31006 (PA440) with respect to 1 (ICRF) at J2000
angles = ephem_rotation3(eph, 1, 31006, time)
```

## Current Limitations
- The supported JPL binary SPK/PCK segments are limited to: type 1, 2, 3, 8, 12, 21. 
- Binary INPOP kernels are not supported. 
- Acceleration and jerk computations are unavailable for SPK segments of type 1 and 21.