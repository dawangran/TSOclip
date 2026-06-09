#!/usr/bin/env sh
set -eu

bin="${1:-./tsoclip}"
repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
tmpdir="${TMPDIR:-/tmp}/tsoclip_tests.$$"
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"' EXIT

fail(){
  printf 'FAIL: %s\n' "$1" >&2
  exit 1
}

write_core_fastq(){
  fq="$1"
  printf '@full\nAAACCCACGTTAGCTAAC\n+\nIIIIIIIIIIIIIIIIII\n@partial\nTTTTTTACGTTA\n+\nIIIIIIIIIIII\n@insert\nGGGGGGACGTTAGCTTAAC\n+\nIIIIIIIIIIIIIIIIIII\n@concat\nCCCCCACGTTAGCTAACACGTTAGCTAAC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@internal\nACGTTAGCTAACAAAAA\n+\nIIIIIIIIIIIIIIIII\n@none\nCCCCCCGGGGGG\n+\nIIIIIIIIIIII\n' > "$fq"
}

run_core(){
  fq="$1"
  tsv="$2"
  trim="$3"
  shift 3
  "$bin" \
    --fastq "$fq" \
    --tso ACGTTAGCTAAC \
    --out-tsv "$tsv" \
    --out-trim-fastq "$trim" \
    --plain-out \
    --tail-window 30 \
    --tso-min-overlap 4 \
    --tso-max-mm 1 \
    --tso-max-mmr 0.20 \
    --tso-max-shift 2 \
    --max-tail-gap 2 \
    --max-chain-gap 6 \
    --strict-n \
    --min-keep-len 1 \
    "$@" >/dev/null
}

test_version(){
  "$bin" --version | awk '$1=="TSOclip" && $2 ~ /^[0-9]/ { ok=1 } END { exit ok ? 0 : 1 }' \
    || fail "--version output"
}

test_core_plain(){
  fq="$tmpdir/core.fastq"
  tsv="$tmpdir/core.tsv"
  trim="$tmpdir/core.trim.fastq"
  write_core_fastq "$fq"
  run_core "$fq" "$tsv" "$trim"

  awk -F '\t' '
    $1=="full"     { if ($4 != 6 || $10 != 0 || $14 != 6 || $17 != "FULL") exit 1; found++ }
    $1=="partial"  { if ($4 != 6 || $14 != 6 || $17 != "PARTIAL") exit 1; found++ }
    $1=="insert"   { if ($4 != 6 || $10 != 1 || $14 != 6 || $17 != "FULL") exit 1; found++ }
    $1=="concat"   { if ($4 != 5 || $14 != 5 || $17 != "FULL") exit 1; found++ }
    $1=="internal" { if ($3 == 0 || $17 != "NO_CUT") exit 1; found++ }
    $1=="none"     { if ($3 != 0 || $17 != "NO_HIT") exit 1; found++ }
    END            { if (found != 6) exit 1 }
  ' "$tsv" || fail "core TSV classification"

  awk '
    NR==2  { if ($0 != "AAACCC") exit 1 }
    NR==6  { if ($0 != "TTTTTT") exit 1 }
    NR==10 { if ($0 != "GGGGGG") exit 1 }
    NR==14 { if ($0 != "CCCCC") exit 1 }
    NR==18 { if ($0 != "ACGTTAGCTAACAAAAA") exit 1 }
    NR==22 { if ($0 != "CCCCCCGGGGGG") exit 1 }
    END    { if (NR != 24) exit 1 }
  ' "$trim" || fail "core trimmed FASTQ"
}

test_emit_only_hit(){
  fq="$tmpdir/emit.fastq"
  tsv="$tmpdir/emit.tsv"
  trim="$tmpdir/emit.trim.fastq"
  write_core_fastq "$fq"
  run_core "$fq" "$tsv" "$trim" --emit-only-hit 1
  awk 'BEGIN{n=0} /^@/{n++} END{exit n==4 ? 0 : 1}' "$trim" \
    || fail "--emit-only-hit output count"
}

test_no_json(){
  fq="$tmpdir/nojson.fastq"
  tsv="$tmpdir/nojson.tsv"
  trim="$tmpdir/nojson.trim.fastq"
  write_core_fastq "$fq"
  run_core "$fq" "$tsv" "$trim" --no-json
  awk -F '\t' 'NR>1 && $18!="[]" { exit 1 } END{ exit 0 }' "$tsv" \
    || fail "--no-json all_hits_json"
}

test_gzip_io(){
  fq="$tmpdir/gzip.fastq"
  gz="$tmpdir/gzip.fastq.gz"
  tsv="$tmpdir/gzip.tsv"
  trim_gz="$tmpdir/gzip.trim.fastq.gz"
  out="$tmpdir/gzip.out.fastq"
  write_core_fastq "$fq"
  gzip -c "$fq" > "$gz"

  "$bin" \
    --fastq "$gz" \
    --tso ACGTTAGCTAAC \
    --out-tsv "$tsv" \
    --out-trim-fastq "$trim_gz" \
    --tail-window 30 \
    --tso-min-overlap 4 \
    --tso-max-mm 1 \
    --tso-max-mmr 0.20 \
    --tso-max-shift 2 \
    --max-tail-gap 2 \
    --max-chain-gap 6 \
    --strict-n \
    --min-keep-len 1 >/dev/null

  gzip -cd "$trim_gz" > "$out"
  awk 'NR==2 { exit $0=="AAACCC" ? 0 : 1 }' "$out" \
    || fail "gzip input/output"
}

test_n_modes(){
  fq="$tmpdir/n.fastq"
  tsv_wild="$tmpdir/n_wild.tsv"
  tsv_strict="$tmpdir/n_strict.tsv"
  trim="$tmpdir/n.trim.fastq"
  printf '@n\nAAAACGTTNGCTAAC\n+\nIIIIIIIIIIIIIIII\n' > "$fq"

  "$bin" --fastq "$fq" --tso ACGTTAGCTAAC --out-tsv "$tsv_wild" --out-trim-fastq "$trim" \
    --plain-out --tail-window 20 --tso-min-overlap 4 --tso-max-mm 1 >/dev/null
  "$bin" --fastq "$fq" --tso ACGTTAGCTAAC --out-tsv "$tsv_strict" --out-trim-fastq "$trim" \
    --plain-out --tail-window 20 --tso-min-overlap 4 --tso-max-mm 1 --strict-n >/dev/null

  awk -F '\t' '$1=="n" { exit $10==0 ? 0 : 1 }' "$tsv_wild" \
    || fail "N wildcard mode"
  awk -F '\t' '$1=="n" { exit $10==1 ? 0 : 1 }' "$tsv_strict" \
    || fail "N strict mode"
}

test_validation(){
  fq="$tmpdir/invalid.fastq"
  tsv="$tmpdir/invalid.tsv"
  trim="$tmpdir/invalid.fastq.out"
  printf '@x\nACGT\n+\nIIII\n' > "$fq"

  if "$bin" --fastq "$fq" --tso ACGT --out-tsv "$tsv" --out-trim-fastq "$trim" \
    --plain-out --tso-min-overlap 4 >/dev/null 2>&1; then
    fail "invalid --tso-min-overlap accepted"
  fi

  if "$bin" --fastq "$fq" --tso ACGTAC --out-tsv "$tsv" --out-trim-fastq "$trim" \
    --plain-out --tso-min-overlap 3 --tso-max-mmr 1.1 >/dev/null 2>&1; then
    fail "invalid --tso-max-mmr accepted"
  fi

  if "$bin" --fastq "$fq" --tso "" --out-tsv "$tsv" --out-trim-fastq "$trim" \
    --plain-out --tso-min-overlap 1 >/dev/null 2>&1; then
    fail "empty --tso accepted"
  fi
}

test_simulator(){
  sim_fq="$tmpdir/sim.fastq"
  sim_truth="$tmpdir/sim.truth.tsv"
  sim_fq_2="$tmpdir/sim2.fastq"
  sim_truth_2="$tmpdir/sim2.truth.tsv"
  sim_tsv="$tmpdir/sim.tsoclip.tsv"
  sim_trim="$tmpdir/sim.trim.fastq"

  python3 "$repo_dir/scripts/simulate_tso_reads.py" \
    --out-fastq "$sim_fq" \
    --out-truth "$sim_truth" \
    --tso ACGTTAGCTAAC \
    --reads-per-scenario 3 \
    --seed 11 \
    --min-cdna-len 40 \
    --max-cdna-len 60

  python3 "$repo_dir/scripts/simulate_tso_reads.py" \
    --out-fastq "$sim_fq_2" \
    --out-truth "$sim_truth_2" \
    --tso ACGTTAGCTAAC \
    --reads-per-scenario 3 \
    --seed 11 \
    --min-cdna-len 40 \
    --max-cdna-len 60

  cmp "$sim_fq" "$sim_fq_2" >/dev/null || fail "simulator FASTQ reproducibility"
  cmp "$sim_truth" "$sim_truth_2" >/dev/null || fail "simulator truth reproducibility"

  awk 'BEGIN{n=0} /^@/{n++} END{exit n==24 ? 0 : 1}' "$sim_fq" \
    || fail "simulator FASTQ read count"
  awk -F '\t' '
    NR==1 { next }
    { seen[$2]++; rows++ }
    END {
      split("full partial concat mismatch insertion deletion internal no_hit", s, " ");
      for (i in s) if (seen[s[i]] != 3) exit 1;
      exit rows==24 ? 0 : 1;
    }
  ' "$sim_truth" || fail "simulator truth scenarios"

  "$bin" \
    --fastq "$sim_fq" \
    --tso ACGTTAGCTAAC \
    --out-tsv "$sim_tsv" \
    --out-trim-fastq "$sim_trim" \
    --plain-out \
    --tail-window 80 \
    --tso-min-overlap 4 \
    --tso-max-mm 2 \
    --tso-max-mmr 0.25 \
    --tso-max-shift 2 \
    --max-tail-gap 2 \
    --max-chain-gap 6 \
    --min-keep-len 1 >/dev/null

  awk 'NR>1 { rows++ } END{ exit rows==24 ? 0 : 1 }' "$sim_tsv" \
    || fail "simulator TSOclip row count"
}

test_version
test_core_plain
test_emit_only_hit
test_no_json
test_gzip_io
test_n_modes
test_validation
test_simulator

printf 'tests ok\n'
