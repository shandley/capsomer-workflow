# capsomer-workflow

A Snakemake workflow for virome analysis that works on contigs first and writes every result into one [Capsomer](https://github.com/shandley/capsomer) dataset.

Status: an early draft. It has run end to end on two simulated samples.

## What it does

Most viromes are dominated by phage, and phage analysis is best developed for assembled contigs. The workflow follows that route:

| Step | Tool | Per sample or pooled |
|---|---|---|
| Trim adapters and low-quality bases | fastp | per sample |
| Assemble | MEGAHIT or metaSPAdes | per sample |
| Decide which contigs are viral, and name them | geNomad | per sample |
| Estimate completeness and contamination | CheckV | per sample |
| Group viral contigs into species-level clusters | vclust | pooled over all samples |
| Map reads to the cluster representatives | minimap2 | per sample |
| Count reads and depth on each representative | CoverM | all samples |
| Write and check the Capsomer dataset | capsomer | all samples |

The identity and coverage cutoffs for clustering, the minimum contig length and the assembler are settings in `config/config.yaml`. Each is recorded in the dataset with the run that used it.

## What it does not do yet

- Remove host reads before assembly.
- Search reads against protein and nucleotide databases. That module is planned for discovery: reads that map to no contig are the ones worth searching for distant relatives of known viruses.
- Predict hosts, lifestyle or gene function.
- Run on a cluster. The rules carry no resource requests yet.

## Run it

The workflow needs [uv](https://docs.astral.sh/uv/) and the tools in `workflow/envs/tools.yaml`. Either put the tools on your `PATH`, or let Snakemake build the environment with `--use-conda`.

```sh
uv sync

# Databases, once (about 8 GB together).
genomad download-database databases
checkv download_database databases

# Edit config/config.yaml and config/samples.tsv, then:
uv run snakemake --cores 12 -n     # dry run
uv run snakemake --cores 12
```

`config/samples.tsv` is tab-separated with one row per sample: `sample_id`, `r1`, `r2`, and any study metadata columns you want carried into the dataset. The read file paths are not copied into the dataset.

MEGAHIT 1.2.9 crashes with more than one thread on Apple Silicon. On such a machine set `megahit: 1` under `threads` in the configuration.

## What you get

`<output>/bundle` is the Capsomer dataset: a folder of Parquet tables and a manifest that describes them. `<output>/bundle/validation.txt` is the result of checking it against the standard. Open the dataset with DuckDB, Python, R or the `capsomer sql` command.

```sh
uv run capsomer sql results/bundle "SELECT run_id, call, count(*) FROM viral_call GROUP BY ALL"
```

Contig names in the dataset are the sample name, two underscores, and the assembler's own name. A contig assembled identically in two samples is one sequence in the dataset, with one source row per sample.

The intermediate files of each tool stay under `<output>` in a folder named for the step.

## Licence

MIT. See `LICENSE`.
