# Code availability wording

The following wording describes the core software publicly provided here:

> The source code for IDEA-design and IDEA-GP is available at
> https://github.com/nariew111/IDEA under the MIT License. The repository includes
> installation instructions, dependency specifications, input/output guides,
> example configurations, a small artificial IDEA-design selection example,
> archived fixed scoring weights, and scripts for exporting and rebuilding the
> documented calibration resources. The dated five-sample calibration
> reconstruction is identified separately from historical analysis inputs.

The initial uploaded software snapshot is commit
[`c00553fd8d40ece452dd30b9e98015e62f6d9c7e`](https://github.com/nariew111/IDEA/tree/c00553fd8d40ece452dd30b9e98015e62f6d9c7e).
Later documentation and framework-figure changes can be tracked separately in
the Git history and [changelog](../CHANGELOG.md). They do not change the
computational engines in this snapshot.

The [scoring and calibration bundle](../resources/idea_design/README.md) includes
the unchanged historical weights, actual historical reader priorities and a new
reconstruction dated 2026-11-13. It must not be described as recovering the
original five-sample ABS or reproducing the manuscript's existing panels.

The manuscript should additionally identify the exact uploaded software commit
or release used for the study and the actual access location of the
manuscript-specific analysis scripts and fixed configurations. See the
[reproducibility record](manuscript_reproducibility.md) for the status of those
materials. The software package version alone does not establish that its
contents are identical to every historical manuscript analysis.

The IDEA-design web resource is at https://alphaindex.zju.edu.cn/IDEA. This is the
panel-generation resource; it should not be described as an online IDEA-GP
prediction service.
