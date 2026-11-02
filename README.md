# Pprc1

## Repository layout

```
src/
  main.R                      Index of all analysis scripts, in run order
  GEO_submission.R            Builds the GEO metadata and processed-file tables
  analysis/
    Annotations/              Reference annotation built once, reused everywhere
    AIDscreen/                Degron screen: alignment -> insertions -> enrichment
    ChIP/                     ChIP-seq and CUT&RUN: demultiplex -> peaks
    PROseq/                   PRO-seq: demultiplex -> counts -> DESeq2
    Fig*.R, FigS*.R           One script per main / supplementary figure
  functions/                  Shared helpers (sourced, not a package)
  scripts/                    Standalone scripts submitted to the cluster
db/                           Intermediate and processed files (not tracked)
Rdata/                        Serialised metadata and processed objects (not tracked)
results/                      Final figure PDFs and summary tables (not tracked)
```

`db/`, `Rdata/` and `results/` are listed in `.gitignore`: they hold large or
regenerable files and are not part of the repository. The pipelines recreate them.

## Analysis order

[`src/main.R`](src/main.R) is the canonical index — it opens each script in sequence
and documents the intended order. Steps depend on the outputs of earlier steps, so run
them in this order on a fresh checkout.


## Dependencies

### R

Developed under R 4.3. CRAN and Bioconductor packages used across the analysis:

```r
install.packages(c("data.table", "tidyverse", "dplyr", "tidyr", "ggplot2",
                   "ggrastr", "ggpointdensity", "ggseqlogo", "viridis",
                   "circlize", "Hmisc", "here", "metafor"))

BiocManager::install(c("DESeq2", "GenomicRanges", "rtracklayer", "Rsamtools",
                       "ComplexHeatmap", "AnnotationDbi", "org.Mm.eg.db",
                       "org.Hs.eg.db", "GO.db", "BSgenome.Mmusculus.UCSC.mm10"))
```

The pipelines also depend on **`vlite`**, an in-house package that wraps the alignment,
peak-calling and cluster-submission steps. The package can be downloaded from
https://github.com/vloubiere/vlite

Anyone running this will need to point that call at their own
checkout of `vlite`, or replace the `vlite::*` calls with equivalent command lines.

## Running this elsewhere

Two things to be aware of before re-running any of it:

- **vlite package.** Scripts point to the local copy of the vlite R package. The package needs to be downloaded to be accessible.
- **Absolute paths.** Scripts set an absolute working directory near the top and write
  intermediates to cluster scratch. Both need changing for another environment.
- **Cluster submission.** Alignment and peak calling are submitted as batch jobs through
  `vlite::bsub` and assume a Slurm cluster. The downstream parsing, statistics
  and figure scripts run on a single machine.

## License

MIT — see [LICENSE](LICENSE).
