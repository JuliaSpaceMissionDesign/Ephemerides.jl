module Ephem
    
    using NodeGraphs: NodeGraph, SimpleGraph, get_nodes, add_edge!, add_vertex!

    export SPK, DAF, position!, velocity!, state!

    const SECONDS_PER_DAY = 86400.0
    
    include("types.jl")
    include("daf.jl")
    include("segment.jl")
    include("spk.jl")
    include("compute.jl")
 
end

    
include("types.jl")
include("daf.jl")
include("segment.jl")
include("spk.jl")