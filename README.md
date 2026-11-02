# PPRC1

Analysis code for the manuscript:

> **_TITLE_**
> _Authors_
> _Journal / preprint, year._ DOI: [_doi_](https://doi.org/_doi_)

This repository contains the code used to process the sequencing data and to generate
every figure panel in the paper. It covers four experiments in mouse embryonic stem
cells: a genome-wide gene-trap degron screen, ChIP-seq, CUT&RUN, and PRO-seq.

## Data availability

Raw and processed sequencing data are deposited in GEO under accession
[**GSE_XXXXXX**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE_XXXXXX).

Per-library sequencing statistics (total reads, uniquely mapping reads, read length)
are reported in [`results/sequencing_stats.txt`](results/sequencing_stats.txt).

## Experiments

| Assay | Genome | Description |
| --- | --- | --- |
| **Degron (AID) gene-trap screen** | mm10 | Haploid mESC gene-trap library carrying an auxin-inducible degron reporter, in two reporter backgrounds (`EF1a`, `Atp5k`) and three reading frames. Insertions are mapped by iPCR and compared between the `GFPlow`, `GFPhigh` and `Input` fractions. |
| **ChIP-seq** | mm10 | Pprc1 and associated factors, paired-end. |
| **CUT&RUN** | mm10 | Paired-end, processed through the same pipeline as the ChIP-seq. |
| **PRO-seq** | mm10 (+ `dm3` spike-in) | Nascent transcription on auxin-mediated Pprc1 depletion, spike-in normalised. |

Sample-level metadata for each assay lives in the three tab-separated tables at the
repository root — `metadata_AIDscreen.txt`, `metadata_ChIP.txt` and
`metadata_PROseq.txt`. These drive the pipelines: every script reads its sample sheet
from one of them rather than hard-coding file lists.

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
them in this order on a fresh checkout:

**1. Annotations** — reference tracks used by all downstream steps.

```
src/analysis/Annotations/gene_trap_compatible_introns.R
src/analysis/Annotations/Retrieve_FANTOM5_CAGE.R
src/analysis/Annotations/EPDnew_promoters.R
src/analysis/Annotations/GRCm38v102_PROseq_annotations.R
```

Gene models come from Ensembl release 102 (GRCm38/mm10); promoters are taken from
EPDnew and corrected against FANTOM5 CAGE.

**2. Degron screen**

```
src/analysis/AIDscreen/01_AIDscreen_pipeline.R    Trimming, Bowtie2, dedup, MAPQ filter
src/analysis/AIDscreen/02_Data_parsing.R          Insertions -> per-gene count tables
src/analysis/AIDscreen/03_Enrichment_CMH.R        CMH / Fisher enrichment per reporter
```

**3. ChIP-seq and CUT&RUN**

```
src/analysis/ChIP/01_ChIP_demultiplexing.R
src/analysis/ChIP/02_ChIP_processing.R            Bowtie2 -> MAPQ30 -> coverage tracks
src/analysis/ChIP/03_ChIP_peakCalling.R           MACS2, paired-end mode
src/analysis/ChIP/04_ChIP_Pprc1confidentPeaks.R   Reproducible peak set across replicates
```

**4. PRO-seq**

```
src/analysis/PROseq/01_Demultiplex.R
src/analysis/PROseq/02_Alignment.R                Bowtie, UMI collapsing, dm3 spike-in
src/analysis/PROseq/03_Counting.R
src/analysis/PROseq/04_DESeq_enrichment.R         DESeq2, spike-in size factors
src/analysis/PROseq/05_PROseq_CNR_combined.R      PRO-seq vs CUT&RUN integration
```

Differentially transcribed genes are called at `padj < 0.05` and
`|log2 fold change| > log2(1.5)`.

**5. Figures** — each script writes its panels as PDFs into `results/`.

```
src/analysis/Fig1.R   Fig2.R   Fig3.R   Fig5.R
src/analysis/FigS1.R  FigS2.R  FigS3.R  FigS4.R  FigS8.R
```

**6. GEO submission**

```
src/GEO_submission.R
```

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
peak-calling and cluster-submission steps. It is not on CRAN and is loaded from a local
path:

```r
devtools::load_all("/groups/stark/vloubiere/vlite")
```

Anyone running this outside the Stark lab will need to point that call at their own
checkout of `vlite`, or replace the `vlite::*` calls with equivalent command lines.

### Command-line tools

`bowtie`, `bowtie2`, `samtools`, `trim_galore`, `macs2`, and `perl` (for the
demultiplexing scripts in `src/functions/`).

## Running this elsewhere

Two things to be aware of before re-running any of it:

- **Absolute paths.** Scripts set an absolute working directory near the top
  (`/groups/stark/vojtek/Bioinfo_project/PPRC1_paper/Pprc1_git/`) and write
  intermediates to cluster scratch (`/scratch-cbe/users/...`). Both need changing for
  another environment.
- **Cluster submission.** Alignment and peak calling are submitted as batch jobs through
  `vlite::bsub` and assume a Slurm cluster (VBC CBE). The downstream parsing, statistics
  and figure scripts run on a single machine.

## License

MIT — see [LICENSE](LICENSE).
