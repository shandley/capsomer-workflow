"""Turn the workflow's outputs into one Capsomer bundle, and check it."""

import csv
import re
import sys
import tomllib
from pathlib import Path
from typing import Any

from capsomer.builder import Builder
from capsomer.importers.checkv import add_checkv
from capsomer.importers.clusters import add_clusters
from capsomer.importers.contigs import add_contigs
from capsomer.importers.genomad import add_genomad
from capsomer.importers.mapping import add_coverm, add_read_mapping
from capsomer.importers.study import add_fastp, init_bundle
from capsomer.validate import conforms, format_report, validate

PIPELINE = "pipeline"
FILE_COLUMNS = ("r1", "r2")


def versions(path: Path) -> dict[str, str]:
    """Tool name -> version, from lines of 'tool<TAB>whatever the tool prints for --version'."""
    found: dict[str, str] = {}
    for line in path.read_text().splitlines():
        tool, _, printed = line.partition("\t")
        if printed.split():
            found[tool] = printed.split()[-1].removeprefix("v")
    return found


def study_table(source: Path, out: Path) -> None:
    """The sample table without its read file columns, which do not belong in a bundle."""
    with source.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    columns = [c for c in rows[0] if c not in FILE_COLUMNS]
    with out.open("w") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def one_answer_per_sequence(builder: Builder) -> None:
    """Flag one geNomad row per sequence as the pipeline's answer.

    Each sample has its own geNomad run. A contig assembled identically in two
    samples is one sequence with a call from each run, and only one can be preferred.
    """
    for table in ("viral_call", "taxon_assignment"):
        builder.con.execute(
            f"UPDATE {table} AS t SET is_preferred = (t.run_id = (SELECT min(o.run_id) FROM {table} o "
            "WHERE o.sequence_id = t.sequence_id AND o.run_id LIKE 'genomad.%')) "
            "WHERE t.run_id LIKE 'genomad.%'"
        )


def build(snakemake: Any) -> bool:
    config = snakemake.config
    results = Path(snakemake.params.results)
    samples: list[str] = list(snakemake.params.samples)
    bundle = results / "bundle"
    bundle.mkdir(parents=True, exist_ok=True)
    for old in [*bundle.glob("*.parquet"), bundle / "manifest.json"]:
        old.unlink(missing_ok=True)
    (results / "work").mkdir(exist_ok=True)
    table = results / "work" / "samples.tsv"
    study_table(Path(snakemake.input.samples), table)
    meta = config.get("bundle") or {}
    init_bundle(bundle, table, bundle_id=meta.get("id"), title=meta.get("title"))

    version = versions(Path(snakemake.input.versions))
    project = tomllib.loads(
        (Path(snakemake.scriptdir).parents[1] / "pyproject.toml").read_text()
    )
    assembler = config["assembler"]
    genomad_db = Path(config["databases"]["genomad"])
    checkv_db = Path(config["databases"]["checkv"])
    genomad_release = (
        (genomad_db / "version.txt").read_text().strip()
        if (genomad_db / "version.txt").exists()
        else None
    )
    checkv_match = re.search(r"v(\d+(?:\.\d+)*)$", checkv_db.name)
    clustering = config["clustering"]
    viruses = results / "clusters" / "viral_sequences.fna"
    representatives = results / "clusters" / "representatives.fna"

    with Builder(bundle) as b:
        b.add_run(
            PIPELINE,
            tool=project["project"]["name"],
            tool_version=project["project"]["version"],
            step="pipeline",
            description="Contig-first virome workflow.",
        )
        for s in samples:
            contigs = results / "contigs" / f"{s}.fasta"
            add_fastp(
                b,
                results / "qc" / f"{s}.fastp.json",
                run_id=f"trim.{s}",
                library_id=s,
                parent=PIPELINE,
            )
            add_contigs(
                b,
                contigs,
                run_id=f"assembly.{s}",
                tool=assembler,
                tool_version=version.get(
                    "megahit" if assembler == "megahit" else "spades.py"
                ),
                libraries=[s],
                inputs=[f"trim.{s}"],
                parent=PIPELINE,
                parameters={
                    "min_contig_length": config["min_contig_length"],
                    "assembler_depth": "k-mer depth from the contig header",
                    "contig_names": "sample name, two underscores, the assembler's name",
                },
            )
            add_genomad(
                b,
                results / "genomad" / s,
                contigs,
                run_id=f"genomad.{s}",
                tool_version=version.get("genomad"),
                database_version=genomad_release,
                inputs=[f"assembly.{s}"],
                parent=PIPELINE,
            )
            add_checkv(
                b,
                results / "checkv" / s,
                results / "genomad" / s / f"{s}_summary" / f"{s}_virus.fna",
                run_id=f"checkv.{s}",
                tool_version=version.get("checkv"),
                database_version=checkv_match.group(1) if checkv_match else None,
                inputs=[f"genomad.{s}"],
                parent=PIPELINE,
            )
        add_clusters(
            b,
            results / "clusters" / "clusters.tsv",
            viruses,
            run_id="cluster",
            tool="vclust",
            tool_version=version.get("vclust"),
            level="species",
            layout="vclust",
            inputs=[f"genomad.{s}" for s in samples],
            parent=PIPELINE,
            parameters={
                "algorithm": clustering["algorithm"],
                "metric": "ani",
                "min_ani": clustering["min_ani"],
                "min_query_coverage": clustering["min_coverage"],
                "members": "sequences geNomad called virus, from every sample",
            },
        )
        for s in samples:
            add_read_mapping(
                b,
                representatives,
                run_id=f"map.{s}",
                tool="minimap2",
                tool_version=version.get("minimap2"),
                library_id=s,
                inputs=[f"trim.{s}", "cluster"],
                parent=PIPELINE,
                parameters={
                    "preset": "sr",
                    "targets": "representatives of the species clusters",
                },
            )
        add_coverm(
            b,
            results / "map" / "coverm.tsv",
            representatives,
            run_id="coverage",
            libraries={s: s for s in samples},
            tool_version=version.get("coverm"),
            inputs=[f"map.{s}" for s in samples],
            parent=PIPELINE,
            parameters={
                "mapped_unit": "reads",
                "methods": "count mean covered_fraction covered_bases",
            },
        )
        one_answer_per_sequence(b)

    findings = validate(bundle)
    report = format_report(findings)
    Path(snakemake.output.report).write_text(report + "\n")
    print(report)
    return conforms(findings)


if __name__ == "__main__":
    job: Any = globals()["snakemake"]
    with open(job.log[0], "w") as log:
        sys.stdout = sys.stderr = log
        ok = build(job)
    sys.exit(0 if ok else 1)
