# TSOclip

**Fast TSO detection & trimming for long reads.**
- Tail-window search (default **150 bp**)
- **Partial = TSO prefix** anywhere inside the tail window
- Full match with ±shift (small indels)
- Concatenation-aware, multi-hit de-noising
- OpenMP multi-threaded, gzip I/O, single-write per read
- `--plain-out` to pipe into `pigz` safely

## Build

```bash
make
# or: gcc -O3 -march=native -pipe -fopenmp -DNDEBUG -o tsoclip src/tsoclip.c -lz
```

## Quick start

Only output **hits** to gz file (program writes gzip itself):

```bash
tsoclip \
  --fastq input.fastq.gz \
  --tso CCCCTCTGCGTTGATACCACTGCTT \
  --out-tsv run.hits.tsv \
  --out-trim-fastq run.trimmed.fastq.gz \
  --emit-only-hit 1 --min-keep-len 50
```

Or **pipe** with `pigz` (program outputs **plain FASTQ**):

```bash
pigz -dc input.fastq.gz | \
tsoclip \
  --fastq - \
  --tso CCCCTCTGCGTTGATACCACTGCTT \
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
--n-as-match \
--threads 16 --batch-size 40000 \
--no-json --gzbuf-kb 1024 --gzip-level 1
```

## Output

- `run.hits.tsv` with columns:
  ```
  read_id  read_len  tso_hit_count  pick_start  pick_end  pick_mm  pick_overlap  pick_mmr  pick_partial  all_hits_json
  ```
- `run.trimmed.fastq.gz`: trimmed FASTQ.

## Notes

- **Do not double-compress**: if you don’t pass `--plain-out`, TSOclip already writes gzip.
- **Partial = prefix**: incomplete TSO at 3′ is matched via TSO **prefix** within tail window.

## License
MIT
