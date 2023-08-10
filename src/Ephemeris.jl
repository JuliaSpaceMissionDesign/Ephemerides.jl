
module Ephemeris 

using ForwardDiff: Dual, partials, value
using Mmap
using PreallocationTools
using StaticArrays

import JSMDInterfaces.Ephemeris as jEph

include("utils.jl")

# SPK segment types definitions
include("spk/spktypes.jl")

# DAF definitions
include("daf.jl")

# SPK segments implementations
include("spk/spk1.jl")
include("spk/spk2.jl")
include("spk/spk3.jl")

# Descriptors linking and provider 
include("links.jl")
include("provider.jl")

include("properties.jl")
include("transform.jl")

# Provide compatibility with JSMDInterfaces
include("interfaces.jl")

end

# TODO: thread-safe version (i.e., tante cache quanti threads)