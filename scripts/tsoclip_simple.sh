#!/usr/bin/env bash
set -euo pipefail

PLATFORM="ont"      # ont|hifi
MODE="balance"      # recall|balance|precision
WINDOW=150
THREADS=16
BATCH=40000
GZBUF=1024
GZLVL=1
EMIT_ONLY=1
MINKEEP=50
NOJSON=1
NASMATCH=1
CONCAT_FULL=0
PLAIN_OUT=0

usage(){ cat <<'EOF'
Usage:
  tsoclip_simple.sh --fastq <in.fq.gz|-> --tso <SEQ> --out-prefix <prefix>
                    [--platform ont|hifi] [--mode recall|balance|precision]
                    [--window 150] [--threads 16] [--plain-out]

Examples:
  pigz -dc input.fastq.gz | \
    ./scripts/tsoclip_simple.sh --fastq - --tso CCCCTCTGCGTTGATACCACTGCTT --out-prefix run1 \
      --platform ont --mode balance --window 150 --threads 16 \
  | pigz -c -p 16 > run1.trimmed.fastq.gz
EOF
}

FASTQ=""; TSO=""; OUTP=""; PLAIN_OUT=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq) FASTQ="$2"; shift 2;;
    --tso) TSO="$2"; shift 2;;
    --out-prefix) OUTP="$2"; shift 2;;
    --platform) PLATFORM="$2"; shift 2;;
    --mode) MODE="$2"; shift 2;;
    --window) WINDOW="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --plain-out) PLAIN_OUT=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown: $1"; usage; exit 1;;
  esac
done
[[ -z "$FASTQ" || -z "$TSO" || -z "$OUTP" ]] && { usage; exit 1; }

# presets
if [[ "$PLATFORM" == "ont" ]]; then
  case "$MODE" in
    recall)   MM=6; SHIFT=5; OVL=10; MMR=0.40; HITS=100; SPACE=6;;
    balance)  MM=5; SHIFT=4; OVL=12; MMR=0.20; HITS=100; SPACE=6;;
    precision)MM=4; SHIFT=3; OVL=12; MMR=0.20; HITS=80;  SPACE=8;;
    *) echo "bad --mode"; exit 2;;
  esac
  NASMATCH=1
else
  case "$MODE" in
    recall)   MM=5; SHIFT=3; OVL=10; MMR=0.30; HITS=80; SPACE=8;;
    balance)  MM=4; SHIFT=2; OVL=12; MMR=0.25; HITS=64; SPACE=10;;
    precision)MM=3; SHIFT=2; OVL=14; MMR=0.20; HITS=48; SPACE=12;;
    *) echo "bad --mode"; exit 2;;
  esac
  NASMATCH=0
fi

COMMON_ARGS=(
  --tso "$TSO"
  --tail-window "$WINDOW"
  --tso-max-mm "$MM"
  --tso-max-shift "$SHIFT"
  --tso-min-overlap "$OVL"
  --tso-max-mmr "$MMR"
  --tso-max-hits "$HITS"
  --min-spacing "$SPACE"
  --threads "$THREADS"
  --batch-size "$BATCH"
  --gzbuf-kb "$GZBUF"
  --gzip-level "$GZLVL"
  --emit-only-hit $EMIT_ONLY
  --min-keep-len "$MINKEEP"
)
[[ "$NOJSON" == "1" ]] && COMMON_ARGS+=( --no-json )
[[ "$NASMATCH" == "1" ]] && COMMON_ARGS+=( --n-as-match )
[[ "$CONCAT_FULL" == "1" ]] && COMMON_ARGS+=( --concat-prefer-full )
[[ "$PLAIN_OUT" == "1" ]] && COMMON_ARGS+=( --plain-out )

TSV="${OUTP}.hits.tsv"

if [[ "$FASTQ" == "-" ]]; then
  ./tsoclip --fastq - "${COMMON_ARGS[@]}" --out-tsv "$TSV" --out-trim-fastq -
else
  pigz -dc "$FASTQ" | ./tsoclip --fastq - "${COMMON_ARGS[@]}" --out-tsv "$TSV" --out-trim-fastq -
fi
