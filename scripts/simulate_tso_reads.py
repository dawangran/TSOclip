#!/usr/bin/env python3
"""Simulate FASTQ reads with benchmark-ready TSO artifact truth labels."""

from __future__ import annotations

import argparse
import random
from dataclasses import dataclass
from pathlib import Path


BASES = "ACGT"


@dataclass
class Artifact:
    sequence: str
    start: int
    end: int
    reason: str
    copies: int
    partial_len: int
    mismatches: int
    insertions: int
    deletions: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate simulated FASTQ reads and TSV truth labels for TSOclip benchmarks."
    )
    parser.add_argument("--out-fastq", required=True, help="Output simulated FASTQ path.")
    parser.add_argument("--out-truth", required=True, help="Output truth TSV path.")
    parser.add_argument("--tso", default="ACGTTAGCTAAC", help="TSO sequence to embed.")
    parser.add_argument("--reads-per-scenario", type=int, default=100, help="Reads per scenario.")
    parser.add_argument("--seed", type=int, default=1, help="Random seed.")
    parser.add_argument("--min-cdna-len", type=int, default=80, help="Minimum cDNA prefix length.")
    parser.add_argument("--max-cdna-len", type=int, default=240, help="Maximum cDNA prefix length.")
    parser.add_argument("--quality-char", default="I", help="Single FASTQ quality character.")
    parser.add_argument(
        "--scenarios",
        default="full,partial,concat,mismatch,insertion,deletion,internal,no_hit",
        help="Comma-separated scenarios: full,partial,concat,mismatch,insertion,deletion,internal,no_hit.",
    )
    return parser.parse_args()


def rand_dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(BASES) for _ in range(length))


def mutate_base(rng: random.Random, base: str) -> str:
    choices = [b for b in BASES if b != base]
    return rng.choice(choices)


def add_mismatch(rng: random.Random, seq: str) -> tuple[str, int]:
    if not seq:
        return seq, 0
    pos = rng.randrange(len(seq))
    return seq[:pos] + mutate_base(rng, seq[pos]) + seq[pos + 1 :], 1


def add_insertion(rng: random.Random, seq: str) -> tuple[str, int]:
    pos = rng.randrange(len(seq) + 1)
    return seq[:pos] + rng.choice(BASES) + seq[pos:], 1


def add_deletion(rng: random.Random, seq: str) -> tuple[str, int]:
    if not seq:
        return seq, 0
    pos = rng.randrange(len(seq))
    return seq[:pos] + seq[pos + 1 :], 1


def make_artifact(rng: random.Random, scenario: str, tso: str, cdna_len: int) -> Artifact:
    if scenario == "full":
        seq = tso
        return Artifact(seq, cdna_len, cdna_len + len(seq), "FULL", 1, 0, 0, 0, 0)

    if scenario == "partial":
        partial_len = rng.randint(max(1, min(8, len(tso) - 1)), len(tso) - 1)
        seq = tso[:partial_len]
        return Artifact(seq, cdna_len, cdna_len + len(seq), "PARTIAL", 1, partial_len, 0, 0, 0)

    if scenario == "concat":
        copies = rng.randint(2, 3)
        seq = tso * copies
        return Artifact(seq, cdna_len, cdna_len + len(seq), "FULL", copies, 0, 0, 0, 0)

    if scenario == "mismatch":
        seq, mismatches = add_mismatch(rng, tso)
        return Artifact(seq, cdna_len, cdna_len + len(seq), "FULL", 1, 0, mismatches, 0, 0)

    if scenario == "insertion":
        seq, insertions = add_insertion(rng, tso)
        return Artifact(seq, cdna_len, cdna_len + len(seq), "FULL", 1, 0, 0, insertions, 0)

    if scenario == "deletion":
        seq, deletions = add_deletion(rng, tso)
        return Artifact(seq, cdna_len, cdna_len + len(seq), "FULL", 1, 0, 0, 0, deletions)

    if scenario in {"internal", "no_hit"}:
        return Artifact("", -1, -1, "NO_CUT" if scenario == "internal" else "NO_HIT", 0, 0, 0, 0, 0)

    raise ValueError(f"unknown scenario: {scenario}")


def make_read(rng: random.Random, scenario: str, tso: str, min_len: int, max_len: int) -> tuple[str, Artifact]:
    cdna_len = rng.randint(min_len, max_len)
    cdna = rand_dna(rng, cdna_len)

    if scenario == "internal":
        insert_at = rng.randint(10, max(10, cdna_len - 10))
        tail = rand_dna(rng, max(12, len(tso) // 2))
        read = cdna[:insert_at] + tso + cdna[insert_at:] + tail
        artifact = Artifact(tso, insert_at, insert_at + len(tso), "NO_CUT", 0, 0, 0, 0, 0)
        return read, artifact

    if scenario == "no_hit":
        return cdna, Artifact("", -1, -1, "NO_HIT", 0, 0, 0, 0, 0)

    artifact = make_artifact(rng, scenario, tso, cdna_len)
    return cdna + artifact.sequence, artifact


def write_outputs(args: argparse.Namespace) -> None:
    rng = random.Random(args.seed)
    tso = args.tso.upper()
    scenarios = [s.strip() for s in args.scenarios.split(",") if s.strip()]
    allowed = {"full", "partial", "concat", "mismatch", "insertion", "deletion", "internal", "no_hit"}
    unknown = sorted(set(scenarios) - allowed)
    if unknown:
        raise SystemExit(f"unknown scenarios: {','.join(unknown)}")
    if not tso:
        raise SystemExit("--tso must not be empty")
    if args.reads_per_scenario <= 0:
        raise SystemExit("--reads-per-scenario must be positive")
    if args.min_cdna_len <= 0 or args.max_cdna_len < args.min_cdna_len:
        raise SystemExit("invalid cDNA length range")
    if len(args.quality_char) != 1:
        raise SystemExit("--quality-char must be a single character")

    fastq_path = Path(args.out_fastq)
    truth_path = Path(args.out_truth)
    fastq_path.parent.mkdir(parents=True, exist_ok=True)
    truth_path.parent.mkdir(parents=True, exist_ok=True)

    truth_header = [
        "read_id",
        "scenario",
        "read_len",
        "artifact_start",
        "artifact_end",
        "expected_cut_start",
        "expected_cut_reason",
        "tso_copies",
        "partial_len",
        "truth_mismatches",
        "truth_insertions",
        "truth_deletions",
        "tso",
    ]

    with fastq_path.open("w", encoding="ascii") as fq, truth_path.open("w", encoding="ascii") as truth:
        truth.write("\t".join(truth_header) + "\n")
        idx = 0
        for scenario in scenarios:
            for rep in range(args.reads_per_scenario):
                idx += 1
                read_id = f"sim_{idx:06d}_{scenario}_{rep + 1:04d}"
                read, artifact = make_read(rng, scenario, tso, args.min_cdna_len, args.max_cdna_len)
                qual = args.quality_char * len(read)
                fq.write(f"@{read_id}\n{read}\n+\n{qual}\n")

                expected_cut_start = artifact.start if artifact.reason in {"FULL", "PARTIAL"} else -1
                row = [
                    read_id,
                    scenario,
                    str(len(read)),
                    str(artifact.start),
                    str(artifact.end),
                    str(expected_cut_start),
                    artifact.reason,
                    str(artifact.copies),
                    str(artifact.partial_len),
                    str(artifact.mismatches),
                    str(artifact.insertions),
                    str(artifact.deletions),
                    tso,
                ]
                truth.write("\t".join(row) + "\n")


def main() -> int:
    write_outputs(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
