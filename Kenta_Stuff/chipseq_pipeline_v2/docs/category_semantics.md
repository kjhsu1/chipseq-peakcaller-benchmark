# Category Semantics

The current controlled simulator baseline is the TF-clean six:

- `flatearth_*`: no GC/accessibility-shaped realism background beyond the flat baseline.
- `wavy_*`: smoother structured background intended to mimic nonuniform but not hotspot-heavy signal.
- `hilly_*`: dense hotspot-like background driven by mappability enrichment; this label is retained for continuity, but it should not be read as "3% of peaks" or sparse hill placement.
- `*_peak_*`: narrower TF-like planted signal with `tf_sigma <= 5`.
- `*_plateau_*`: broader histone-like planted signal with `tf_sigma >= 15`.

Important interpretation limits:

- "plateau" versus "peak" is directionally useful, but the naming is only partly trustworthy because the benchmark still couples shape labels to current simulator settings rather than an external morphological ontology.
- Some category differences are confounded with MACS2 broad versus narrow calling mode.
- The current hilly labeling corresponds to a dense hotspot process, not a literal percentage of loci.
