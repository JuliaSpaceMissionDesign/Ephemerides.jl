
@testset "SPK Segments" verbose=true begin 

    struct TestSPK <: Ephemerides.AbstractSPKSegment end

    @test_throws ErrorException Ephemerides.spk_field(TestSPK()) 
    @test_throws ErrorException Ephemerides.header(TestSPK())
    @test_throws ErrorException Ephemerides.cache(TestSPK())

    include("spk1.jl")
    include("spk2.jl")
    include("spk3.jl")
    include("spk8.jl")
    include("spk9.jl")
    include("spk12.jl")
    include("spk13.jl")
    include("spk18.jl")
    include("spk19.jl")
    include("spk21.jl")

end;