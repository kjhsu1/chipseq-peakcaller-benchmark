# ChIP-seq Simulator Comparison

## Compatibility Checklist

A simulator is compatible with this benchmark if it can support, directly or by
a simple documented wrapper:

- treatment and control/input read generation
- explicit control-depth sweeps
- reproducible seed control
- FASTQ output for reuse by the existing alignment and peak-calling workflows
- known binding regions or a defensible truth-proxy handoff for ontology windows
- parameterization from real ChIP-seq data or a clear reason manual parameters
  are scientifically adequate
- Snakemake and Slurm execution without interactive steps
- realistic ChIP-seq noise sources such as shearing, pulldown, PCR duplicates,
  sequencing error, and genome-wide background

## Recommendation

Use **ChIPs** as the first realistic simulation engine.

ChIPs is the best fit for this project because it explicitly models major
ChIP-seq experimental steps, provides a `learn` module for estimating model
parameters from existing aligned reads and peaks, and provides `simreads` for
generating FASTQ from peaks plus model/experiment parameters. The `simreads`
mode also supports whole-cell extract control simulation using `-t wce`, which
maps cleanly onto the control-depth sweep idea.

The key caveat is that ChIPs still requires peak/truth inputs. For this project,
those should be treated as reference/truth-proxy regions derived from real data
or controlled truth regions, not as literal biological truth.

## Candidate Table

| Simulator | Evidence | Strengths | Weaknesses | Install/runtime concerns | Fit for control-depth sweeps | Fit for ontology evaluation | Recommendation |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ChIPs | BMC Bioinformatics 2021; GitHub; Bioconda | Learns parameters from BAM+peaks; simulates shearing, pulldown, PCR, sequencing; emits FASTQ; supports WCE control mode | Requires peak/truth input; number of copies is not identical to real cell count | Best cluster path is a Bioconda environment; local Mac execution is not required for final sweeps | Strong: `--numreads`, `--seed`, and WCE control mode support sweep wrappers | Strong: peak BED input gives a natural truth-proxy window handoff | Primary choice |
| ChIPulate | PLOS Computational Biology 2019; GitHub | Detailed TF binding and experimental-factor model | Focused on TF occupancy; ChIPs paper notes it does not simulate background fragments outside peak regions | Requires a separate custom wrapper and likely more manual parameter fitting | Moderate for TF-only experiments; weaker for broad mark/control background benchmark | Moderate: binding regions exist, but background ontology is less natural | Keep as background comparison, not primary |
| isChIP | Briefings in Bioinformatics 2021 | Command-line; uses BED template of true binding events; generates FASTQ-like outputs | Does not infer model parameters from existing datasets according to ChIPs comparison | Possible fallback, but parameterization burden is higher than ChIPs | Moderate: can likely vary read/cycle settings | Strong if template truth is well defined | Fallback if ChIPs installation fails |
| ChIPsim | Bioconductor | Mature R/Bioconductor package; general ChIP-seq simulation framework | Current package focus is nucleosome positioning; less direct CLI fit for this repo | Would add an R/Bioconductor workflow path and more environment surface area | Moderate but would require R workflow integration | Moderate; possible but heavier wrapper burden | Not primary |
| Direct real pileup resampling | Existing draft idea | Simple conceptually; directly uses selected real studies | Treats one noisy experiment too much like a true generator; copies experiment-specific noise | Easy mechanically, but scientific risk is high | Easy mechanically but scientifically weak | Truth definition is weak | Reject as primary |

## Sources

- ChIPs paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04097-5
- ChIPs GitHub: https://github.com/gymreklab/chips
- ChIPs Bioconda package: https://anaconda.org/bioconda/chips
- ChIPulate paper: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006921
- isChIP paper: https://academic.oup.com/bib/article/doi/10.1093/bib/bbaa352/6035271
- ChIPsim package: https://bioconductor.org/packages/release/bioc/html/ChIPsim.html

## Decision

The realistic benchmark should use ChIPs first, while preserving the selected
real worm broad-mark and ENCODE SKN-1 studies as reference inputs for learning,
calibration, and ontology labels. If ChIPs installation or runtime behavior
fails on the cluster, `isChIP` is the most reasonable fallback because it keeps
the command-line, peak-template, and FASTQ-output shape closest to this project.
