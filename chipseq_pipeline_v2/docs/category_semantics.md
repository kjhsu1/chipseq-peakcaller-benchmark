# Category Semantics

The current controlled simulator baseline is the TF-clean six:

- `flatearth_*`: no GC/accessibility-shaped realism background beyond the flat baseline.
- `wavy_*`: smoother structured background intended to mimic nonuniform but not hotspot-heavy signal.
- `hilly_*`: dense hotspot-like background driven by mappability enrichment; this label is retained for continuity, but it should not be read as "3% of peaks" or sparse hill placement.
- `*_peak_*`: narrower TF-like planted signal with `tf_sigma <= 5`.
- `*_plateau_*`: broader histone-like planted signal with `tf_sigma >= 15`.

The decode-mode contract is now explicit:

- `*_peak_*` categories must use narrow decoding.
- `*_plateau_*` categories must use broad decoding.
- The workflow now validates this invariant from `tf_sigma` at build time so a
  peak/broad or plateau/narrow mismatch fails immediately.

Important interpretation limits:

- "plateau" versus "peak" is directionally useful, but the naming is only partly trustworthy because the benchmark still couples shape labels to current simulator settings rather than an external morphological ontology.
- The current hilly labeling corresponds to a dense hotspot process, not a literal percentage of loci.
