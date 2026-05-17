# Scientific Validity Report: ChIPs Real-Study Ontology Analysis

## Purpose

This report evaluates whether downstream ontology analysis is statistically,
logically, and biologically valid for the real-study ChIP-seq/ChIPs workflow
before implementing that step in Snakemake.

## Evidence Reviewed

- ChIPs documentation states that `learn` estimates model parameters from an
  aligned BAM plus peak file and writes a JSON model. It also states that
  `simreads` can simulate treatment reads from BED/HOMER peaks or simulate
  whole-cell-extract control data with `-t wce`.
  Source: https://github.com/gymreklab/chips
- The ChIPs paper describes ChIPs as a flexible ChIP-seq simulation toolkit and
  emphasizes inference of model parameters from real data as a key advantage.
  Source: https://doi.org/10.1186/s12859-021-04097-5
- The Bioconda package page describes ChIPs as a tool for simulating
  ChIP-sequencing experiments and provides package availability for version
  `2.4`.
  Source: https://bioconda.github.io/recipes/chips/README.html
- ENCODE ChIP-seq pipeline documentation separates mapping, signal generation,
  peak calling, replicate treatment, and reproducibility handling. It also
  emphasizes matching input controls and replicate/run-type structure.
  Sources:
  https://github.com/ENCODE-DCC/chip-seq-pipeline and
  https://www.encodeproject.org/chip-seq/histone-encode3/
- MACS2 expects ChIP treatment input and can use a control/input file; broad
  peak mode is explicitly supported with `--broad`.
  Source: https://carta.tech/man-pages/man1/macs2.1.html

## Key Finding

ChIPs is appropriate for realistic read simulation and control-depth sweeps, but
it does not expose a clean closed-form intrinsic PMF for every parameter set in
the way the older synthetic simulator concept did. Therefore, ontology analysis
should not claim exact PMF-known ground truth.

The valid interpretation is:

- The real-study peak BED used by ChIPs is a signal-template truth proxy.
- The learned JSON model is a study-level empirical simulation model.
- The generated treatment/control BAMs provide empirical coverage behavior for
  each simulated run.
- The final MACS2 peak BEDs can be evaluated against the signal-template regions
  and grouped by ontology labels derived from region shape, score, enrichment,
  and background behavior.

## Statistical Validity

The strongest statistically valid analysis is not a binary declaration that a
region is truly bound or unbound. Instead, it is a sensitivity analysis:

- Does peak-caller recovery change as control depth changes?
- Which classes of regions lose recall first?
- Which classes of regions accumulate false positives?
- Are broad domains and narrow peaks affected differently?
- Do high-background regions behave differently from low-background regions?

This is valid because the sweep varies coverage and seed while holding the
study-derived simulation model fixed. It lets us estimate conditional
performance patterns rather than absolute biological truth.

The main statistical limits are:

- ChIPs stochastic output means each seed is one draw from the model, not a full
  analytical distribution.
- The input peak BED is a truth proxy, not a perfect truth set.
- If the real-study peak calls used for training contain false positives or miss
  real sites, the ontology inherits that bias.
- MACS2 output is thresholded and interval-shaped, so per-region overlap metrics
  must be interpreted as operational recovery, not exact molecular binding.

## Logical Validity

The workflow is logically sound if the ontology step answers this question:

For regions used as ChIPs signal templates, how does MACS2 recovery vary by
region class and control depth after realistic simulated read generation?

The workflow is not logically sound if it claims:

- ChIPs gives the exact intrinsic PMF for every locus.
- Every non-template region is guaranteed biological null.
- A called peak outside the template is automatically a simulator error.

To avoid those problems, the ontology workflow should use careful labels:

- `truth_proxy_region`
- `called_overlap`
- `called_region_recovered`
- `template_region_missed`
- `control_depth_sensitive`

## Biological Validity

The design is biologically reasonable because it preserves major ChIP-seq
features:

- treatment signal is modeled from real ChIP-seq reads and real peak calls
- control reads are generated as WCE/background reads
- narrow and broad modes are separated through the selected study and MACS2 mode
- read alignment and peak calling are performed downstream, rather than scoring
  simulator output directly

The main biological limits are:

- The current workflow merges available selected replicates to learn one model
  per study. This is acceptable for a benchmark simulation, but it is not the
  same as ENCODE-style replicate reproducibility or IDR analysis.
- The broad H3K9me2-like study may be sensitive to repetitive sequence and
  multi-mapping behavior; ontology labels should include background-sensitive
  classes where possible.
- The control simulation is WCE-like background, not a real matched input
  library with all lab-specific biases.

## Adjusted Experimental Plan

The ontology downstream step should be implemented as empirical truth-proxy
analysis:

1. For each ChIPs run, use the same real-study peak BED that was supplied to
   `chips simreads` as the template/truth-proxy region set.
2. Read the simulated treatment and control BAMs for that run.
3. Read the called MACS2 peak BED for that run.
4. For each template region, compute:
   - region width
   - template peak score
   - treatment read/fragment overlap count
   - control read/fragment overlap count
   - log2 treatment/control enrichment with a pseudocount
   - whether the region overlaps any called peak
5. Classify each region using the existing ontology helper thresholds.
6. Aggregate precision/recall/F1 by ontology class and control depth.
7. Report the result as truth-proxy recovery, not exact intrinsic-PMF recovery.

## Implementation Requirements

The Snakemake workflow should add:

- a per-run ontology scoring rule downstream of ChIPs MACS2 peak calls
- a per-run ontology classification rule using existing `classify_regions.py`
- a combined ontology table rule across all ChIPs runs
- a final ontology summary rule using existing `evaluate_by_region_ontology.py`

The implementation should avoid introducing new dependencies. The repo already
uses Python, pandas, matplotlib, and pysam in the realstudy scripts, which are
sufficient for per-region BAM overlap counting and summary tables.

## Conclusion

Ontology analysis is scientifically defensible with ChIPs if framed as
empirical truth-proxy analysis. The workflow should not claim access to an exact
intrinsic PMF. It should evaluate called-peak recovery against ChIPs template
regions and generated treatment/control coverage, grouped into biologically
interpretable ontology classes.
