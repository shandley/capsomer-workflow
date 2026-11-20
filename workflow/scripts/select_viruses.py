"""Choose the viral sequences of one sample from two classifiers.

Everything geNomad lists as a virus is kept, proviruses included. A contig
that only Jaeger calls viral is added whole, unless geNomad listed it as a
plasmid or found a provirus on it. The table says which tool or tools called
each sequence, so a later reader can keep only the ones both agree on.

Jaeger adds most on contigs under about 5 kb, where geNomad's default rule
asks for a hallmark gene. That is also where Jaeger errs most: in our test it
called 7% of non-viral 2 kb pieces viral. A minimum length for contigs taken
on Jaeger's call alone trades those finds against those errors.
"""

import argparse
import csv
from collections.abc import Iterator
from pathlib import Path

JAEGER_VIRAL = {"phage", "virus"}


def records(path: Path) -> Iterator[tuple[str, str]]:
    name, chunks = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("contigs", type=Path)
    parser.add_argument("genomad_viruses", type=Path, help="geNomad's *_virus.fna")
    parser.add_argument(
        "genomad_plasmids", type=Path, help="geNomad's *_plasmid_summary.tsv"
    )
    parser.add_argument("fasta_out", type=Path)
    parser.add_argument("table_out", type=Path)
    parser.add_argument("--jaeger", type=Path, help="Jaeger's per-sequence table")
    parser.add_argument(
        "--jaeger-min-length",
        type=int,
        default=0,
        help="shortest contig accepted on Jaeger's call alone",
    )
    args = parser.parse_args()

    jaeger: set[str] = set()
    if args.jaeger:
        with args.jaeger.open() as handle:
            jaeger = {
                row["contig_id"]
                for row in csv.DictReader(handle, delimiter="\t")
                if row["prediction"] in JAEGER_VIRAL
            }
    with args.genomad_plasmids.open() as handle:
        plasmids = {row["seq_name"] for row in csv.DictReader(handle, delimiter="\t")}

    used: set[str] = set()
    with args.fasta_out.open("w") as fasta, args.table_out.open("w") as table:
        table.write("name\tcontig\tevidence\n")
        for name, bases in records(args.genomad_viruses):
            contig = name.split("|provirus_")[0]
            used.add(contig)
            whole = contig == name
            evidence = "geNomad and Jaeger" if whole and contig in jaeger else "geNomad"
            fasta.write(f">{name}\n{bases}\n")
            table.write(f"{name}\t{contig}\t{evidence}\n")
        for name, bases in records(args.contigs):
            long_enough = len(bases) >= args.jaeger_min_length
            if (
                name in jaeger
                and name not in used
                and name not in plasmids
                and long_enough
            ):
                fasta.write(f">{name}\n{bases}\n")
                table.write(f"{name}\t{name}\tJaeger\n")


if __name__ == "__main__":
    main()
