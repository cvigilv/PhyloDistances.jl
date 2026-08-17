# PhyloDistances.jl

Distance and similarity metrics between phylogenetic trees.

Trees are read from Newick with [NewickTree.jl](https://github.com/arzwa/NewickTree.jl) and
traversed through the
[AbstractTrees.jl](https://github.com/JuliaCollections/AbstractTrees.jl) interface, so any
node type implementing that interface is accepted.

## Interface

Metrics implement the [Distances.jl](https://github.com/JuliaStats/Distances.jl) interface.
A metric is a value, and applying it to two trees compares them:

```julia
using PhyloDistances

d = SomeMetric()(tree1, tree2)
d = evaluate(SomeMetric(), tree1, tree2)   # the same thing
```

Everything that selects a variant of the computation is a field of the metric, so a
constructed metric specifies the computation completely and no keyword arguments are needed
when applying it:

```julia
SomeMetric(; convention = :primary, normalize = true)(tree1, tree2)
SomeParameterizedMetric(2; normalize = true)(tree1, tree2)
```

That is what lets the whole Distances.jl toolkit work unchanged:

```julia
D = pairwise(SomeMetric(), trees)          # all-pairs matrix
pairwise!(D, SomeMetric(), trees)          # in place
colwise(SomeMetric(), trees, references)   # elementwise along two collections
```

`pairwise`, `pairwise!`, `colwise`, `colwise!`, `evaluate` and `result_type` are
re-exported from Distances.jl, so loading both packages is safe: the names refer to the
same functions rather than clashing.

Distance matrices produced this way feed directly into ecosystem tools that expect one —
clustering, multidimensional scaling, and so on.

### Distances and similarities

Distances subtype `TreeMetric`, which is a `Distances.SemiMetric`.

Quantities that are *largest* on identical trees — mutual clustering information, shared
phylogenetic information, Nye similarity, maximum agreement subtree size — subtype
`TreeSimilarity` instead, which deliberately sits outside the Distances.jl hierarchy.
`Distances.PreMetric` requires `d(x, x) == 0`, which a similarity does not satisfy, and
declaring one a `SemiMetric` would be actively wrong: `pairwise` takes the zero diagonal on
faith rather than computing it, so the result would report zero self-similarity. Both kinds
are used identically; only the type hierarchy differs.

### Taxa and rooting

Trees are read with `readnw`, re-exported from NewickTree.jl. `taxonindex` gives a tree's
canonical taxon ordering — labels sorted, so positions depend on the taxon set rather than
on the order leaves happened to appear — and every array a metric builds is indexed by
those positions:

```julia
tree = readnw("(((Human:0.1,Chimp:0.1):0.2,Gorilla:0.3):0.1,Orang:0.5,Gibbon:0.6);")
index = taxonindex(tree)
taxa(index)      # ["Chimp", "Gibbon", "Gorilla", "Human", "Orang"]
index["Human"]   # 4
```

The two-argument form additionally requires that both trees span the same taxon set and
names the labels that differ otherwise. Repeated labels within a tree are rejected: they
would make the correspondence between two trees ambiguous.

Newick carries no explicit marker for rooting, so `PhyloDistances.isrooted` reads it from
the root's child count — two means rooted, three or more means the root is an arbitrary
starting point. When that does not match what a metric is defined on:

- **a rooted tree given to an unrooted metric warns and proceeds.** The root's two child
  branches induce the same bipartition, so discarding the root position leaves a
  well-defined unrooted tree, and the warning says so.
- **an unrooted tree given to a rooted metric throws.** There is no canonical way to root
  it — midpoint and outgroup rooting give different answers — so the choice belongs to the
  caller. Root the tree explicitly, for instance with `NewickTree.set_outgroup`.

Node types other than `NewickTree.Node` are accepted once
`PhyloDistances.taxonlabel(leaf)` is defined for them; everything else goes through the
AbstractTrees.jl interface.

### Conventions

Several tree metrics have more than one formulation in the literature, differing in
normalization, in the treatment of trivial splits, or in the base of the logarithm. The
`convention` field selects which is computed, so a value can be traced back to a specific
definition:

```julia
SomeMetric(; convention = :treedist)   # default
SomeMetric(; convention = :primary)
```

`:treedist` reproduces the R package [TreeDist](https://github.com/ms609/TreeDist);
`:primary` follows the source that first defined the metric, where the two differ.

### Normalization

The `normalize` field scales the result onto a range comparable across taxon set sizes. It
takes more than a flag:

```julia
SomeMetric()(t1, t2)                      # raw value
SomeMetric(; normalize = true)(t1, t2)    # the metric's own scheme
SomeMetric(; normalize = max)(t1, t2)     # relative to the more informative tree
SomeMetric(; normalize = min)(t1, t2)     # relative to the less informative tree
SomeMetric(; normalize = 12)(t1, t2)      # divided by a value you supply
```

A function is applied to what each tree carries — its split count for a split-based
metric, its information content for an information-based one — so `max` and `min` express
the result relative to one tree or the other. `true` selects the scheme the metric defines
for itself, which for Robinson-Foulds is the sum, the total splits present in the pair.

### Traits

`PhyloDistances.issimilarity(m)` reports whether larger values mean *more* similar trees;
`PhyloDistances.requiresrooted(m)` reports whether the metric is defined on rooted trees;
`PhyloDistances.convention(m)` and `PhyloDistances.isnormalized(m)` report how a
constructed metric is configured.
