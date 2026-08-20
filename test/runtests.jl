using Test

@testset "PhyloDistances.jl" begin
    include("test_taxa.jl")
    include("test_random.jl")
    include("test_splits.jl")
    include("test_information.jl")
    include("test_clustertable.jl")
    include("test_assignment.jl")
    include("test_splitmatching.jl")
    include("test_interface.jl")
    include("test_robinsonfoulds.jl")
    include("test_quartet.jl")
    include("test_generalizedrf.jl")
    include("test_fixtures.jl")
end
