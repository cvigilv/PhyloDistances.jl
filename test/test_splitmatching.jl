using PhyloDistances
using PhyloDistances: splitmatching
using Test

"""The canonical mask over `index` marking exactly `labels`, oriented away from taxon 1."""
function splitmask(index, labels...)
    m = falses(length(index))
    for label in labels
        m[index[label]] = true
    end
    return first(m) ? .!m : m
end

# A trivial scorer independent of any real metric: the number of shared marked taxa. Its
# optimum is hand-verifiable, which is what makes it useful here rather than reusing a
# metric's own scorer to test the framework that metric depends on.
_overlap(a, b, ntaxa) = count(a .& b)

@testset "identical split sets match every split to itself" begin
    tree = readnw("(T1:1.0,T2:1.0,((T3:1.0,T4:1.0):1.0,T5:1.0):1.0);")
    index = taxonindex(tree)
    s = splits(tree, index)

    result = splitmatching(_overlap, s, s; maximize = true)
    @test result.matching == 1:length(s)
    @test result.score == sum(count, s)
end

@testset "hand-computed rectangular matching" begin
    # Tree 1 has three informative splits: {T5,T6}, {T3,T4}, {T3,T4,T5,T6}. Tree 2 has two:
    # {T3,T4}, {T3,T4,T5}. Every pairing's overlap is small enough to check by hand:
    #   {T5,T6}  vs {T3,T4}=0, {T3,T4,T5}=1
    #   {T3,T4}  vs {T3,T4}=2, {T3,T4,T5}=2
    #   {T3,T4,T5,T6} vs {T3,T4}=2, {T3,T4,T5}=3
    # The only way to place both of tree 2's splits at distinct tree-1 splits and beat a
    # total of 4 is {T3,T4}<->{T3,T4} (2) with {T3,T4,T5,T6}<->{T3,T4,T5} (3), total 5 —
    # which leaves {T5,T6} unmatched even though it has a nonzero overlap available.
    tree = readnw("(T1:1,T2:1,((T3:1,T4:1):1,(T5:1,T6:1):1):1);")
    other = readnw("(T1:1,T2:1,T6:1,((T3:1,T4:1):1,T5:1):1);")
    index = taxonindex(tree, other)
    s1, s2 = splits(tree, index), splits(other, index)
    @test length(s1) == 3
    @test length(s2) == 2

    result = splitmatching(_overlap, s1, s2; maximize = true)
    @test result.score == 5

    i56 = findfirst(==(splitmask(index, "T5", "T6")), collect(s1))
    i34 = findfirst(==(splitmask(index, "T3", "T4")), collect(s1))
    i3456 = findfirst(==(splitmask(index, "T3", "T4", "T5", "T6")), collect(s1))
    j34 = findfirst(==(splitmask(index, "T3", "T4")), collect(s2))
    j345 = findfirst(==(splitmask(index, "T3", "T4", "T5")), collect(s2))

    @test result.matching[i56] == 0
    @test result.matching[i34] == j34
    @test result.matching[i3456] == j345
end

@testset "symmetry: swapping the trees swaps the matching's direction, not its score" begin
    a = readnw("(T1:1,T2:1,((T3:1,T4:1):1,T5:1):1);")
    b = readnw("(T1:1,T2:1,T3:1,(T4:1,T5:1):1);")
    index = taxonindex(a, b)
    sa, sb = splits(a, index), splits(b, index)

    forward = splitmatching(_overlap, sa, sb; maximize = true)
    backward = splitmatching((x, y, n) -> _overlap(y, x, n), sb, sa; maximize = true)
    @test forward.score == backward.score
end

@testset "minimize picks the worst-overlap pairing instead" begin
    tree = readnw("(T1:1,T2:1,((T3:1,T4:1):1,T5:1):1);")
    index = taxonindex(tree)
    s = splits(tree, index)

    maxresult = splitmatching(_overlap, s, s; maximize = true)
    minresult = splitmatching(_overlap, s, s; maximize = false)
    @test maxresult.score >= minresult.score
end

@testset "mismatched taxon orderings are rejected" begin
    a = splits(readnw("(T1:1,T2:1,(T3:1,T4:1):1);"))
    b = splits(readnw("(X1:1,X2:1,(X3:1,X4:1):1);"))
    @test_throws "different taxon orderings" splitmatching(_overlap, a, b)
end
