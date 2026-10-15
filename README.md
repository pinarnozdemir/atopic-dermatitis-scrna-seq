# Atopic Dermatitis | Single-Cell RNA-seq Analysis

**An R-based scRNA-seq analysis workflow comparing lesional and non-lesional skin in atopic dermatitis (GSE147424).**

This project presents a bioinformatics workflow developed during an internship, along with Turkish-language analysis notes. It explores cellular heterogeneity, differential gene expression, pathway enrichment, and inferred cell–cell communication in atopic dermatitis.

> **Project status:** Research / internship portfolio project. The original analysis script is provided as written; the repository has not been independently rerun or validated in a clean environment. Interpretation and methods should be reviewed before reuse in a publication.

## Research question

How do cell populations, gene expression programs, and inferred intercellular communication differ between lesional and non-lesional skin samples from individuals with atopic dermatitis?

**Important distinction:** Non-lesional skin is not the same as skin from unaffected healthy controls.

## Workflow

1. Load an **existing** Seurat object (`combined.rds`) representing GSE147424.
2. Inspect quality-control metrics and filter cells by detected feature count.
3. Use `scDblFinder` for doublet prediction and keep predicted singlets.
4. Normalize, select variable features, scale, and run PCA.
5. Apply Harmony sample-level integration and perform clustering, UMAP, and t-SNE.
6. Identify cluster markers and assess representative skin-cell markers.
7. Annotate clusters with SingleR using Human Primary Cell Atlas and Monaco immune references.
8. Run lesional-versus-non-lesional differential expression analysis, both globally and by annotated cell type.
9. Analyze GO Biological Process and KEGG enrichment; export gene lists for external Metascape analysis.
10. Estimate and compare cell–cell communication networks using CellChat.

## Repository structure

```text
atopic-dermatitis-scrna-seq/
├── README.md
├── .gitignore
├── scripts/
│   └── atopic_dermatitis_analysis.R   # Original analysis code, unchanged
└── docs/
    └── Atopic_Dermatitis_Analysis_Notes_TR.docx
```

## Data and prerequisites

- **Dataset accession:** [GSE147424 — NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147424).
- The script **does not** download GEO data, build a Seurat object, or convert raw matrices into a Seurat object.
- To execute it as written, supply a pre-existing Seurat object at `~/Desktop/GSE147424_RAW/combined.rds` and use that directory as the working directory.
- The object is expected to contain an RNA assay and metadata fields including `sample_id` and `condition`, with condition labels `Lesional` and `Non_lesional`. Standard Seurat QC fields (`nFeature_RNA`, `nCount_RNA`) are also used.
- Data preparation steps for producing `combined.rds` are **not included** in this repository, and the repository is **not end-to-end reproducible from raw GEO files**.

### R packages used

`Seurat`, `dplyr`, `ggplot2`, `patchwork`, `Matrix`, `SingleR`, `celldex`, `SingleCellExperiment`, `scDblFinder`, `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `AnnotationDbi`, `CellChat`, `NMF`, `ggalluvial`, `ComplexHeatmap`, `harmony`.

Install versions compatible with your R, Bioconductor and operating system environments. No package lockfile or exact versions were supplied with the original project.

## How to run

After preparing the required Seurat object and installing all R dependencies:

```bash
Rscript scripts/atopic_dermatitis_analysis.R
```

The script includes fixed **absolute paths**, so it currently assumes the `~/Desktop/GSE147424_RAW` directory. Adjust both the `setwd(...)` and `readRDS(...)` paths **locally** if necessary; the committed analysis script remains unchanged.

The workflow writes plots, DEG/enrichment files, CellChat summaries, Metascape input lists, and intermediate/final `.rds` objects to the working directory. These generated files are excluded from Git tracking by `.gitignore`.

## Analysis decisions and limitations

- The script applies `nFeature_RNA > 200` and `nFeature_RNA < 6000` as QC thresholds.
- Mitochondrial-gene percentage is **not used** for filtering because the source script indicates mitochondrial gene features could not be reliably identified in the input object.
- The provided input appears to be processed; the script itself warns that `scDblFinder` predictions should be interpreted cautiously.
- Harmony correction uses `sample_id`; the script uses the first 20 Harmony dimensions and clustering resolution 0.5.
- Differential expression uses Seurat's Wilcoxon-based `FindMarkers` on cells, not donor-level pseudobulk modeling. If samples originate from paired donors, donor dependence should be accounted for before making inferential claims.
- GO/KEGG and cell-communication results are computational predictions and do not on their own establish mechanism or causal relationships.
- No numerical result is claimed in this README unless supported by the committed script or notes. The full Turkish analysis/interpretation document is available in `docs/`.

## Documentation

See [`docs/Atopic_Dermatitis_Analysis_Notes_TR.docx`](docs/Atopic_Dermatitis_Analysis_Notes_TR.docx) for the accompanying Turkish descriptions of the workflow and interpretations. This is a project notebook, **not a peer-reviewed manuscript**.

## Author

**Pınar Nur Özdemir**  
Bioengineering | Single-cell transcriptomics and bioinformatics

## Citation and reuse

Please cite the original GEO dataset and relevant analysis software if you build upon this workflow. Original dataset licensing, participant metadata, and use restrictions should be reviewed at the GEO record before redistributing source data. No license for this repository has been specified.

## Turkish analysis notes

The full Turkish-language methodology and explanations are available as [GitHub-readable Markdown](docs/ANALIZ_NOTLARI_TR.md). The [original Word document](docs/Atopic_Dermatitis_Analysis_Notes_TR.docx) contains the embedded figures.
