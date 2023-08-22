using Ephemerides
using Test 

using CalcephEphemeris
using JSMDUtils.Math: D¹, D², D³
using LazyArtifacts

import JSMDInterfaces.Ephemeris as jEphem

@testset "Download all artifacts" begin
    @info artifact"testdata"
    @info "All artifacts downloaded"
end;


@testset "Ephemerides" verbose=true begin 
    include("spk/spk.jl")
    include("properties.jl")
    include("interfaces.jl")
end;