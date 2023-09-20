# Loading Kernels

This tutorial will walk you through the basic features and interfaces that allow you to load binary ephemeris kernels.

The supported sources of ephemerides are currently limited to binary PCK and SPK segments of type: 1, 2, 3, 8, 9, 12, 13, 18, 19, 20 and 21. 

!!! note 
    Support for IMCCE INPOP ephemerides is yet to be implemented.

Before retrieving position and orientation data of celestial objects, the user is first required to load the ephemerides files into an [`EphemerisProvider`](@ref) object. 

```julia
using Ephemerides 

# Load a single ephemeris file 
eph1 = EphemerisProvider("kernel1.bsp")

# Load multiple ephemeris files simultaneously
eph2 = EphemerisProvider(["kernel1.bsp", "kernel2.bsp"])
```

You must specify the relative or absolute path(s) of the file(s) to load. Either one or multiple ephemeris files can be simultaneously loaded into a single `EphemerisProvider` object. However, 
once a provider has been created, no more kernels can be loaded inside it. 