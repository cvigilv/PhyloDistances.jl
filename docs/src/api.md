# API reference

## Tree distances

```@docs
PhyloDistances.RobinsonFoulds
PhyloDistances.WeightedRobinsonFoulds
PhyloDistances.InfoRobinsonFoulds
PhyloDistances.QuartetDistance
PhyloDistances.JaccardRobinsonFoulds
PhyloDistances.ClusteringInfoDistance
PhyloDistances.MatchingSplitInfoDistance
```

## Tree similarities

```@docs
PhyloDistances.NyeSimilarity
PhyloDistances.MutualClusteringInfo
PhyloDistances.SharedPhylogeneticInfo
```

## Comparison interface

```@docs
PhyloDistances.TreeComparison
PhyloDistances.TreeMetric
PhyloDistances.TreeSimilarity
PhyloDistances.Convention
PhyloDistances.TreeDistConvention
PhyloDistances.PrimaryConvention
PhyloDistances.convention
PhyloDistances.normalization
PhyloDistances.isnormalized
PhyloDistances.issimilarity
PhyloDistances.requiresrooted
PhyloDistances.normalizer
PhyloDistances.normalizerinfo
PhyloDistances.pairwise
```

## Tree data

```@docs
PhyloDistances.TaxonIndex
PhyloDistances.taxa
PhyloDistances.taxonindex
PhyloDistances.taxonlabel
PhyloDistances.taxonlabels
PhyloDistances.isrooted
PhyloDistances.branchlength
PhyloDistances.Splits
PhyloDistances.splits
PhyloDistances.istrivial
PhyloDistances.incidencematrix
```

## Information and matching

```@docs
PhyloDistances.SplitInfoTable
PhyloDistances.log2rooted
PhyloDistances.log2unrooted
PhyloDistances.splitinfo
PhyloDistances.clusteringentropy
PhyloDistances.mutualinformation
PhyloDistances.jointentropy
PhyloDistances.splitmatching
```

## Tree generation

```@docs
PhyloDistances.randomtree
PhyloDistances.perturb
PhyloDistances.nni!
```

## Package

```@docs
PhyloDistances
```
