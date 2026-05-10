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

Windows:
https://cran.r-project.org/bin/windows/base/

macOS:
https://cran.r-project.org/bin/macosx/

Linux:
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
```

Or download ZIP from GitHub and extract it.

---

### Step 4. Install Required Packages

Open R or RStudio.

Run:

```r
source("install_cotra_packages.R")
```

The installer:
- skips already installed packages
- installs CRAN packages
- installs Bioconductor packages
- installs GitHub packages
- minimizes user interaction

Installation may take time depending on internet speed and system performance.

---

# Platform Specific Setup

## Windows

Usually no extra setup required.

If package compilation fails:
- install Rtools

Download:
https://cran.r-project.org/bin/windows/Rtools/

Restart R after installation.

---

## macOS

Install Xcode command line tools:

```bash
xcode-select --install
```

Install Homebrew if needed:

https://brew.sh/

Recommended libraries:

```bash
brew install cmake hdf5 harfbuzz fribidi
```

---

## Linux Ubuntu/Debian

Install required system libraries:

```bash
sudo apt update

sudo apt install -y \
  build-essential \
  gfortran \
  cmake \
  git \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libhdf5-dev \
  libfontconfig1-dev \
  libcairo2-dev \
  libxt-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype6-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev
```

---

# Running CoTRA

Open R or RStudio inside the CoTRA folder.

Run:

```r
shiny::runApp()
```

Or:

```r
source("app.R")
```

The application will open automatically in your browser.

---

# Docker Installation

## Install Docker

Download Docker Desktop:

https://www.docker.com/products/docker-desktop/

---

## Build and Run CoTRA

Inside the CoTRA folder:

```bash
docker compose up --build
```

Open browser:

```text
http://localhost:3838
```

---

# Dockerfile Example

```dockerfile
FROM rocker/shiny:latest

WORKDIR /srv/shiny-server/CoTRA

COPY . /srv/shiny-server/CoTRA

RUN Rscript install_cotra_packages.R

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/CoTRA', host='0.0.0.0', port=3838)"]
```

---

# Project Structure

```text
CoTRA/
├── app.R
├── global.R
├── server.R
├── ui.R
├── modules/
├── CoTRA_installer/
├── data/
├── outputs/
├── www/
├── Dockerfile
└── README.md
```

---

# Recommended Input Data

## Bulk RNA-seq
- raw count matrix

## scRNA-seq
- Seurat object
- Cell Ranger output
- expression matrices
- metadata annotations

---

# Output

CoTRA generates:
- publication ready figures
- interactive plots
- downloadable tables
- PDF reports
- SVG figures
- session files

---

# Citation

If you use CoTRA in your research, please cite the related publication when available.

---

# License

MIT License

---

# Support

For issues and feature requests:

https://github.com/YOUR_USERNAME/CoTRA/issues

---

# Acknowledgements

CoTRA integrates multiple open source R and Bioconductor packages including:
- Seurat
- clusterProfiler
- ReactomePA
- monocle3
- slingshot
- CellChat
- NicheNet
- GSVA
- enrichplot
- DESeq2
- limma
- fgsea
- plotly
- Shiny

---

# Future Development

Planned additions:
- Spatial transcriptomics
- Multiomics integration
- GPU acceleration
- Cloud deployment
- Automated workflow export
- AI assisted interpretation
