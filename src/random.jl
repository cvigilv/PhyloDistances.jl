"""
    randomtree([rng], ntaxa; labels, branchlength) -> NewickTree.Node

An unrooted binary tree on `ntaxa` taxa, drawn uniformly from the topologies on that many
taxa.

Growing the tree by repeatedly splitting a uniformly chosen branch gives every unrooted
binary topology equal probability, because each of the `2k - 3` branches of a `k`-taxon
tree is an equally likely place for taxon `k + 1` to attach.

`labels` names the taxa, defaulting to `T1` … `T<ntaxa>`. `branchlength` is called with
the random number generator for each branch and defaults to a uniform draw on `[0, 1)`.
The root has three children, marking the result unrooted.

Passing a seeded generator makes the result reproducible:

```julia
randomtree(Xoshiro(1), 12)
```
"""
function randomtree(
    rng::Random.AbstractRNG,
    ntaxa::Integer;
    labels::AbstractVector = ["T$i" for i in 1:ntaxa],
    branchlength = rand,
)
    ntaxa >= 3 || throw(ArgumentError(
        "an unrooted tree needs at least three taxa, got $ntaxa"
    ))
    length(labels) == ntaxa || throw(DimensionMismatch(
        "got $(length(labels)) labels for $ntaxa taxa"
    ))

    nextid = Ref(0)
    newid() = UInt16(nextid[] += 1)
    leafdata(name) = NewickTree.NewickData(branchlength(rng), 0.0, String(name))

    # Three taxa on a degree-3 root: the only unrooted topology on three taxa, and the
    # seed every larger one grows from.
    root = NewickTree.Node(newid(), NewickTree.NewickData())
    branches = [NewickTree.Node(newid(), leafdata(labels[i]), root) for i in 1:3]

    for i in 4:ntaxa
        # Every non-root node stands for the branch above it, so choosing one uniformly
        # chooses a branch uniformly.
        below = rand(rng, branches)
        parent = below.parent

        # The new node takes the chosen branch's place, and the branch hangs beneath it
        # alongside the new taxon.
        joint = NewickTree.Node(
            newid(), NewickTree.NewickData(branchlength(rng), 0.0, "")
        )
        parent.children[findfirst(child -> child === below, parent.children)] = joint
        joint.parent = parent

        push!(joint, below)
        leaf = NewickTree.Node(newid(), leafdata(labels[i]), joint)

        push!(branches, joint)
        push!(branches, leaf)
    end

    return root
end

randomtree(ntaxa::Integer; kwargs...) = randomtree(Random.default_rng(), ntaxa; kwargs...)

"""
    nni!([rng], tree) -> tree

Apply one nearest-neighbour interchange to `tree` in place, and return it.

An interchange picks an internal branch and swaps a subtree hanging off one end with one
hanging off the other, which is the smallest rearrangement that changes a topology. The
taxon set is unchanged.

Throws for a tree with no internal branch, which has no interchange to make.
"""
function nni!(rng::Random.AbstractRNG, tree)
    # An internal branch is one above a node that is neither the root nor a leaf.
    internal = [
        node for node in AbstractTrees.PostOrderDFS(tree)
        if node !== tree && !isempty(AbstractTrees.children(node))
    ]
    isempty(internal) && throw(ArgumentError(
        "tree has no internal branch to rearrange; " *
        "an unrooted tree needs at least four taxa"
    ))

    below = rand(rng, internal)
    above = below.parent

    # Swapping a subtree beside the branch with one beneath it moves exactly one subtree
    # across, and cannot detach anything, since a sibling is never inside the subtree.
    beside = rand(rng, [c for c in above.children if c !== below])
    beneath = rand(rng, below.children)

    above.children[findfirst(child -> child === beside, above.children)] = beneath
    below.children[findfirst(child -> child === beneath, below.children)] = beside
    beside.parent, beneath.parent = below, above

    return tree
end

nni!(tree) = nni!(Random.default_rng(), tree)

"""
    perturb([rng], tree, nmoves) -> NewickTree.Node

A copy of `tree` with `nmoves` nearest-neighbour interchanges applied.

Distance from `tree` grows with `nmoves` in expectation, which is what makes this a way to
check a metric that has no closed form. It is only an expectation: moves are drawn
independently and a later one can undo an earlier one, so the result is at most `nmoves`
interchanges away and often fewer.

`tree` itself is not modified.
"""
function perturb(rng::Random.AbstractRNG, tree, nmoves::Integer)
    nmoves >= 0 || throw(ArgumentError("nmoves must be non-negative, got $nmoves"))

    perturbed = deepcopy(tree)
    for _ in 1:nmoves
        nni!(rng, perturbed)
    end
    return perturbed
end

perturb(tree, nmoves::Integer) = perturb(Random.default_rng(), tree, nmoves)
