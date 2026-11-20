"""Keep contigs of a minimum length and put the sample name in front of each contig name."""

import argparse
from collections.abc import Iterator
from pathlib import Path


def records(path: Path) -> Iterator[tuple[str, str]]:
    header, chunks = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembly", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--min-length", type=int, default=0)
    args = parser.parse_args()
    with args.out.open("w") as out:
        for header, bases in records(args.assembly):
            if len(bases) >= args.min_length:
                out.write(f">{args.sample}__{header}\n{bases}\n")


if __name__ == "__main__":
    main()
