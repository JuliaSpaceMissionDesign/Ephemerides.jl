
module Ephemerides 

using LazyArtifacts
using Mmap
using PreallocationTools
using PrecompileTools: PrecompileTools
using StaticArrays

import JSMDInterfaces.Ephemeris as jEph

# Utilities
include("utils.jl")

# Interpolation algorithms 
include("interp/cache.jl")
include("interp/hermite.jl")
include("interp/lagrange.jl")
include("interp/chebyshev.jl")

# Twobody utilities 
include("twobody.jl")

# SPK segment types definitions
include("spk/spktypes.jl")

# DAF definitions
include("daf.jl")

# SPK segments implementations
include("spk/spk1.jl")
include("spk/spk2.jl")
include("spk/spk5.jl")
include("spk/spk8.jl")
include("spk/spk9.jl")
include("spk/spk14.jl")
include("spk/spk17.jl")
include("spk/spk18.jl")
include("spk/spk19.jl")
include("spk/spk20.jl")

# Descriptors linking and provider 
include("links.jl")
include("provider.jl")

include("properties.jl")
include("transform.jl")

# Provide compatibility with JSMDInterfaces
include("interfaces.jl")

# Package precompilation routines
include("precompile.jl")

end