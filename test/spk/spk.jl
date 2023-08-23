
@testset "SPK Segments" verbose=true begin 

    struct TestSPK <: Ephemerides.AbstractSPKSegment end

    @test_throws Ephemerides.spk_field(TestSPK()) ErrorException
    @test_throws Ephemerides.header(TestSPK()) ErrorException
    @test_throws Ephemerides.cache(TestSPK()) ErrorException

    include("spk1.jl")
    include("spk2.jl")
    include("spk8.jl")
    include("spk12.jl")
    include("spk21.jl")

end;