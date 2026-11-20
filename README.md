# capsomer-workflow

A Snakemake workflow for virome analysis that works on contigs first and writes every result into one [Capsomer](https://github.com/shandley/capsomer) dataset.

Status: an early draft. It has run end to end on two simulated samples.

The workflow is opinionated. Each step uses the tool that did best in a comparison on sequences of known origin, and the reason is stated below. Conventional choices that did worse are left out.

## What it does

Most viromes are dominated by phage, and phage analysis is best developed for assembled contigs. The workflow follows that route:

| Step | Tool | Per sample or pooled |
|---|---|---|
| Trim adapters and low-quality bases | fastp | per sample |
| Assemble | MEGAHIT or metaSPAdes | per sample |
| Decide which contigs are viral, and name them | geNomad and Jaeger | per sample |
| Estimate completeness and contamination | CheckV | per sample |
| Group viral contigs into species-level clusters | vclust | pooled over all samples |
| Map reads to the cluster representatives | minimap2 | per sample |
| Count reads and depth on each representative | CoverM | all samples |
| Write and check the Capsomer dataset | capsomer | all samples |

## How viral contigs are chosen

A contig is treated as viral if geNomad calls it viral, or if Jaeger does and geNomad did not call it a plasmid. Both tools' calls are stored in the dataset, along with the workflow's decision and which tool or tools stood behind it, so you can keep only the contigs both agree on.

The pairing comes from a [comparison of five classifiers](https://github.com/shandley/capsomer/tree/main/examples/classifier_benchmark) on 5 kb pieces, which included 600 viral genomes made public after every tool's database was built:

- geNomad and Jaeger found the most viruses and made the fewest calls on non-viral DNA.
- They work differently (marker genes against a neural network on the sequence alone) and miss different viruses. On pieces unlike anything in RefSeq, geNomad found 79% alone and 90% with Jaeger.
- VirSorter2 found fewer viruses than geNomad at every level of novelty and took six to nine times as long, so it is not included.

On one simulated test sample, adding Jaeger raised the truly viral contigs recovered from 73 to 80 of 112, with no non-viral contig added. That sample contains only RefSeq viruses.

Set `jaeger: enabled: false` in the configuration to use geNomad alone.

The identity and coverage cutoffs for clustering, the minimum contig length and the assembler are settings in `config/config.yaml`. Each is recorded in the dataset with the run that used it.

## What it does not do yet

- Remove host reads before assembly.
- Search reads against protein and nucleotide databases. That module is planned for discovery: reads that map to no contig are the ones worth searching for distant relatives of known viruses.
- Predict hosts, lifestyle or gene function.
- Run on a cluster. The rules carry no resource requests yet.

## Run it

The workflow needs [uv](https://docs.astral.sh/uv/), the tools in `workflow/envs/tools.yaml`, and Jaeger, which has its own environment in `workflow/envs/jaeger.yaml`. Either put the tools on your `PATH` and give the path to `jaeger` in the configuration, or let Snakemake build the environments with `--use-conda`.

```sh
uv sync

# Databases, once (about 8 GB together).
genomad download-database databases
checkv download_database databases
jaeger download --path databases/jaeger --model_name jaeger_448931d3_1.3M_fragment

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
