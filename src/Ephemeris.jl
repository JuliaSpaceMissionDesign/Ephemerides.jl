
using Mmap
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

# Descriptors linking and provider 
include("links.jl")
include("provider.jl")

include("properties.jl")
include("compute.jl")
include("orient.jl")

# Provide compatibility with JSMDInterfaces
# include("interfaces.jl")

# TODO: thread-safe version (i.e., tante cache quanti threads)