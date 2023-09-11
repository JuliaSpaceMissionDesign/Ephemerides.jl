
module Ephemerides 

using Mmap
using PreallocationTools
using StaticArrays

import JSMDInterfaces.Ephemeris as jEph

# Utilities
include("utils.jl")

# Interpolation algorithms 
include("interp/cache.jl")
include("interp/hermite.jl")
include("interp/lagrange.jl")

# SPK segment types definitions
include("spk/spktypes.jl")

# DAF definitions
include("daf.jl")

# SPK segments implementations
include("spk/spk1.jl")
include("spk/spk2.jl")
include("spk/spk3.jl")
include("spk/spk8.jl")
include("spk/spk9.jl")
include("spk/spk18.jl")

# Descriptors linking and provider 
include("links.jl")
include("provider.jl")

include("properties.jl")
include("transform.jl")

# Provide compatibility with JSMDInterfaces
include("interfaces.jl")

end