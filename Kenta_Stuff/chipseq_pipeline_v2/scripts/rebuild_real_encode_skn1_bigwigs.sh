#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Rebuild pooled real-data SKN-1 BigWigs from ENCODE ce11 BAMs.

Required arguments:
  --treat1 PATH    Treatment replicate 1 BAM (ENCFF775SIN in the original run)
  --treat2 PATH    Treatment replicate 2 BAM (ENCFF927CRC in the original run)
  --control PATH   Control/input BAM (ENCFF904DXM in the original run)

Optional arguments:
  --outdir PATH        Output directory
  --threads INT        Worker threads for samtools/deepTools (default: 2)
  --bin-size INT       BigWig bin size (default: 10)
  --extend-reads INT   Single-end fragment extension length (default: 125)
  --help               Show this help text

Example:
  conda run -n chipseq_align bash scripts/rebuild_real_encode_skn1_bigwigs.sh \
    --treat1 /path/to/ENCFF775SIN.bam \
    --treat2 /path/to/ENCFF927CRC.bam \
    --control /path/to/ENCFF904DXM.bam \
    --outdir real_data_rebuilds/encode_skn1_ce11_pooled_rebuild_20260503_154535
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

require_tool() {
  command -v "$1" >/dev/null 2>&1 || die "Required executable not found in PATH: $1"
}

copy_input() {
  local src="$1"
  local dest_dir="$2"
  local base
  base="$(basename "$src")"
  cp -f "$src" "$dest_dir/$base"
}

sort_and_index() {
  local src_bam="$1"
  local out_bam="$2"
  local threads="$3"
  samtools sort -@ "$threads" -o "$out_bam" "$src_bam"
  samtools index "$out_bam"
}

main() {
  local treat1=""
  local treat2=""
  local control=""
  local outdir=""
  local threads=2
  local bin_size=10
  local extend_reads=125

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --treat1) treat1="$2"; shift 2 ;;
      --treat2) treat2="$2"; shift 2 ;;
      --control) control="$2"; shift 2 ;;
      --outdir) outdir="$2"; shift 2 ;;
      --threads) threads="$2"; shift 2 ;;
      --bin-size) bin_size="$2"; shift 2 ;;
      --extend-reads) extend_reads="$2"; shift 2 ;;
      --help|-h) usage; exit 0 ;;
      *) die "Unknown argument: $1" ;;
    esac
  done

  [[ -n "$treat1" ]] || die "--treat1 is required"
  [[ -n "$treat2" ]] || die "--treat2 is required"
  [[ -n "$control" ]] || die "--control is required"
  [[ -f "$treat1" ]] || die "Input BAM not found: $treat1"
  [[ -f "$treat2" ]] || die "Input BAM not found: $treat2"
  [[ -f "$control" ]] || die "Input BAM not found: $control"

  if [[ -z "$outdir" ]]; then
    outdir="real_data_rebuilds/encode_skn1_ce11_pooled_rebuild_$(date +%Y%m%d_%H%M%S)"
  fi

  require_tool samtools
  require_tool bamCoverage
  require_tool bamCompare

  export MPLCONFIGDIR="${MPLCONFIGDIR:-$(mktemp -d)}"

  mkdir -p "$outdir"/{downloads,bam,bw,logs,meta}

  copy_input "$treat1" "$outdir/downloads"
  copy_input "$treat2" "$outdir/downloads"
  copy_input "$control" "$outdir/downloads"

  local treat1_name treat2_name control_name pooled_name
  treat1_name="$(basename "${treat1%.bam}")"
  treat2_name="$(basename "${treat2%.bam}")"
  control_name="$(basename "${control%.bam}")"
  pooled_name="${treat1_name}_${treat2_name}.pooled"

  sort_and_index "$outdir/downloads/$(basename "$treat1")" "$outdir/bam/${treat1_name}.sorted.bam" "$threads"
  sort_and_index "$outdir/downloads/$(basename "$treat2")" "$outdir/bam/${treat2_name}.sorted.bam" "$threads"
  sort_and_index "$outdir/downloads/$(basename "$control")" "$outdir/bam/${control_name}.sorted.bam" "$threads"

  samtools merge -f "$outdir/bam/${pooled_name}.bam" \
    "$outdir/bam/${treat1_name}.sorted.bam" \
    "$outdir/bam/${treat2_name}.sorted.bam"
  samtools index "$outdir/bam/${pooled_name}.bam"

  samtools flagstat "$outdir/bam/${treat1_name}.sorted.bam" > "$outdir/logs/${treat1_name}.sorted.flagstat.txt"
  samtools flagstat "$outdir/bam/${treat2_name}.sorted.bam" > "$outdir/logs/${treat2_name}.sorted.flagstat.txt"
  samtools flagstat "$outdir/bam/${control_name}.sorted.bam" > "$outdir/logs/${control_name}.sorted.flagstat.txt"
  samtools flagstat "$outdir/bam/${pooled_name}.bam" > "$outdir/logs/${pooled_name}.flagstat.txt"

  bamCoverage \
    -b "$outdir/bam/${pooled_name}.bam" \
    -o "$outdir/bw/${pooled_name}.ce11.CPM.extend${extend_reads}.bin${bin_size}.bw" \
    --binSize "$bin_size" \
    --normalizeUsing CPM \
    --extendReads "$extend_reads" \
    -p "$threads"

  bamCompare \
    -b1 "$outdir/bam/${pooled_name}.bam" \
    -b2 "$outdir/bam/${control_name}.sorted.bam" \
    -o "$outdir/bw/${pooled_name}_vs_${control_name}.ce11.log2ratio.extend${extend_reads}.bin${bin_size}.bw" \
    --binSize "$bin_size" \
    --extendReads "$extend_reads" \
    --operation log2 \
    --pseudocount 1 1 \
    --scaleFactorsMethod readCount \
    -p "$threads"

  {
    echo "samtools $(samtools --version | sed -n '1p' | awk '{print $2}')"
    echo "bamCoverage $(bamCoverage --version | awk '{print $2}')"
    echo "bamCompare $(bamCompare --version | awk '{print $2}')"
  } > "$outdir/meta/tool_versions.txt"
}

main "$@"
