# CoTRA

Comprehensive Toolbox for RNA-seq Data Analysis

CoTRA is an interactive Shiny based platform for bulk RNA-seq and single cell RNA-seq analysis. The platform integrates commonly used transcriptomics workflows into a user friendly interface for data preprocessing, visualization, differential expression analysis, enrichment analysis, trajectory inference, cell communication, and downstream biological interpretation.

---

# Features

## Bulk RNA-seq
- Quality control visualization
- Sample grouping and filtering
- Differential expression analysis
- Volcano plots
- Heatmaps
- PCA visualization
- Functional enrichment
- GO enrichment
- KEGG enrichment
- Reactome enrichment
- Hallmark enrichment
- GSEA analysis
- Interactive tables and plots
- PDF and SVG export
- Automated report generation

## Single cell RNA-seq
- Seurat workflow integration
- QC filtering
- Clustering
- Dimensionality reduction
- Marker gene analysis
- Cell type annotation
- Pathway activity analysis
- Trajectory analysis
- Cell communication analysis
- RNA velocity support
- Interactive visualization
- Session export and restoration

---

# System Requirements

## Recommended
- 16 GB RAM minimum
- 32 GB RAM recommended for large scRNA-seq datasets
- Multi-core CPU recommended

## Software
- R >= 4.3
- RStudio recommended
- Git optional
- Docker optional

---

# Installation

## Option 1. Standard Installation

### Step 1. Install R

Download and install R:

- Windows:
  https://cran.r-project.org/bin/windows/base/

- macOS:
  https://cran.r-project.org/bin/macosx/

- Linux:
  https://cran.r-project.org/

---

### Step 2. Install RStudio

Download:
https://posit.co/download/rstudio-desktop/

---

### Step 3. Clone or Download CoTRA

Using Git:

```bash
git clone https://github.com/YOUR_USERNAME/CoTRA.git
