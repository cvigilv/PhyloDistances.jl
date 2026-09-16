```@raw html
---
layout: home

hero:
  name: "PhyloDistances.jl"
  text: "Compare phylogenetic trees in Julia"
  tagline: "Tree distances and similarities for AbstractTrees-compatible trees."
  actions:
    - theme: brand
      text: API reference
      link: /api
    - theme: alt
      text: View on GitHub
      link: https://github.com/cvigilv/PhyloDistances.jl

features:
  - title: Tree comparisons
    details: "Robinson-Foulds, quartet, split-matching, clustering-information, and phylogenetic-information measures."
  - title: Julia interfaces
    details: "Distance types implement Distances.jl and accept NewickTree.jl trees or other AbstractTrees-compatible types."
---
```

## What PhyloDistances.jl does

PhyloDistances.jl computes distances and similarities between phylogenetic trees that span
the same taxa. Metrics cover topology, branch lengths, split matching, and information
content.

Each metric is a callable value that also works with the `pairwise` and `colwise` functions
from Distances.jl. Trees can come from NewickTree.jl or another node type that implements
the AbstractTrees.jl interface.
