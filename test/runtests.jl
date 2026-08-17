using Test

@testset "PhyloDistances.jl" begin
    include("test_taxa.jl")
    include("test_splits.jl")
    include("test_interface.jl")
end
