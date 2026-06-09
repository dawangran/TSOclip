# TSOclip

**Fast TSO detection & trimming for long reads.**
- Tail-window search (default **150 bp**)
- **Partial = TSO prefix** anywhere inside the tail window
- Full and partial hits scored with bounded edit distance
- Tail-anchored concatenation-aware hit chaining
- Optional OpenMP multi-threading, gzip I/O, single-write per read
- Per-read cut coordinates plus edit/mismatch/insertion/deletion evidence
- `--plain-out` to pipe into `pigz` safely

## Build

```bash
make
# OpenMP is enabled automatically when the compiler supports -fopenmp.

make check
# or: cc -O3 -march=native -pipe -DNDEBUG -o tsoclip src/tsoclip.c -lz
```

The binary reports its version with:

```bash
./tsoclip --version
```

Install to `/usr/local/bin` by default:

```bash
make install
# or choose a prefix:
make install PREFIX="$HOME/.local"
```

## Tests

```bash
make check
```

The regression suite covers full/partial TSO hits, a single insertion, concatenated TSO chains, internal TSO-like matches that should not be trimmed, `--emit-only-hit`, `--no-json`, gzip input/output, `N` wildcard versus `--strict-n`, and invalid parameter handling.

## Simulated Benchmark Data

Generate reproducible FASTQ reads plus per-read truth labels:

```bash
python3 scripts/simulate_tso_reads.py \
  --out-fastq sim.reads.fastq \
  --out-truth sim.truth.tsv \
  --tso ACGTTAGCTAAC \
  --reads-per-scenario 100 \
  --seed 1
```

The simulator emits `full`, `partial`, `concat`, `mismatch`, `insertion`, `deletion`, `internal`, and `no_hit` scenarios. The truth TSV contains the expected cut start, expected cut reason, TSO copy count, partial length, and simulated mismatch/insertion/deletion counts.

Example run against simulated reads:

```bash
./tsoclip \
  --fastq sim.reads.fastq \
  --tso ACGTTAGCTAAC \
  --out-tsv sim.hits.tsv \
  --out-trim-fastq sim.trimmed.fastq \
  --plain-out \
  --tail-window 150 \
  --tso-min-overlap 4
```

## Docker

```bash
docker build -t tsoclip:$(cat VERSION) .
docker run --rm tsoclip:$(cat VERSION) --version
```

## Quick start

Only output **hits** to gz file (program writes gzip itself):

```bash
tsoclip \
  --fastq input.fastq.gz \
  --tso NNNNNNNNNNNNNNNNNNNNNN \
  --out-tsv run.hits.tsv \
  --out-trim-fastq run.trimmed.fastq.gz \
  --emit-only-hit 1 --min-keep-len 50
```

Or **pipe** with `pigz` (program outputs **plain FASTQ**):

```bash
pigz -dc input.fastq.gz | \
tsoclip \
  --fastq - \
  --tso NNNNNNNNNNNNNNNNNNNNNN \
  --out-tsv run.hits.tsv \
  --out-trim-fastq - \
  --plain-out \
| pigz -c -p 16 > run.trimmed.fastq.gz
```

## Recommended params 

```
--tail-window 150 \
--tso-min-overlap 12 --tso-max-mmr 0.20 \
--tso-max-mm 5 --tso-max-shift 4 \
--tso-max-hits 100 --min-spacing 6 \
--max-tail-gap 2 --max-chain-gap 6 \
--threads 16 --batch-size 40000 \
--no-json --gzbuf-kb 1024 --gzip-level 1
```

## Output

- `run.hits.tsv` with columns:
  ```
  read_id  read_len  tso_hit_count  pick_start  pick_end  pick_mm  pick_overlap  pick_mmr  pick_partial  pick_edits  pick_ins  pick_del  pick_error_rate  cut_start  cut_end  cut_len  cut_reason  all_hits_json
  ```
- `run.trimmed.fastq.gz`: trimmed FASTQ.

## Notes

- **Do not double-compress**: if you don’t pass `--plain-out`, TSOclip already writes gzip.
- **Partial = prefix**: incomplete TSO at 3′ is matched via TSO **prefix** within tail window.
- Trimming requires a hit chain anchored near the 3′ end. `--max-tail-gap` controls how many bases may remain after the tail-most hit, and `--max-chain-gap` controls how far apart chained TSO hits may be.
- `pick_*` reports the hit used for trimming. `cut_start`, `cut_end`, and `cut_len` report the actual removed suffix.
- Hits that are detected inside the tail window but are not connected to the 3′ end are reported in `all_hits_json` with `cut_reason=NO_CUT`.
- `pick_mmr` is retained as a compatibility column and now reports the same edit error rate as `pick_error_rate`.
- `--n-as-match` is the default and treats `N` in either the read or TSO as a wildcard. Use `--strict-n` to count `N` as an ordinary base.
- `--tso-max-mm` is retained for compatibility but now means maximum edit distance for full-length TSO hits.

## Release

The release checklist is in `docs/RELEASE.md`. For manuscript or App Note work, cite an archived release rather than a moving branch.

## License
MIT
