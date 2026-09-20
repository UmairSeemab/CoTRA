# CoTRA

CoTRA (Comprehensive Toolbox for RNA Sequencing Data Analysis) is an R/Shiny application for bulk RNA-seq and single-cell RNA-seq analysis.

CoTRA can be installed in two ways:

1. **Container installation (recommended)** using Docker on Windows, macOS, and Ubuntu/Linux, or Apptainer on HPC systems.
2. **Native R installation** directly in R/RStudio.

The container route is recommended when you want a predefined environment and want to avoid installing the full R, CRAN, Bioconductor, GitHub, and system-library dependency stack manually.

For native installation, CoTRA requires:

```text
R >= 4.4.0
```

---

# 1. Docker container installation

The CoTRA container includes R, CoTRA, its declared R-package dependencies, HDF5 support, Linux system libraries required by the dependency stack, and Pandoc.

The container image is:

```text
ghcr.io/umairseemab/cotra:latest
```

CoTRA is available in the browser at:

```text
http://localhost:3838
```

The container uses:

```text
/data
/results
```

for persistent input data and results.

> The container provides the CoTRA analysis environment. Upstream command-line tools such as FastQC, STAR, HTSeq, and MultiQC are not part of the current CoTRA application container.

---

## Windows with Docker

### Requirements

Install and start Docker Desktop.

No separate R, RStudio, CRAN, or Bioconductor installation is required.

### Recommended launcher

Clone or download the CoTRA repository and run:

```text
launchers/CoTRA-Windows.bat
```

or from PowerShell:

```powershell
.\launchers\CoTRA-Windows.ps1
```

The launcher pulls the CoTRA image, creates persistent data and result folders, starts CoTRA, and opens it in your browser.

To stop CoTRA:

```powershell
docker stop cotra
```

---

## macOS with Docker

### Requirements

Install and start Docker Desktop.

The container supports Intel/AMD64 and Apple Silicon/ARM64 systems.

After cloning or downloading the CoTRA repository:

```bash
chmod +x launchers/CoTRA-macOS.command
./launchers/CoTRA-macOS.command
```

To stop CoTRA:

```bash
docker stop cotra
```

---

## Ubuntu/Linux with Docker

### Requirements

Install Docker Engine and make sure your user can run Docker.

After cloning the CoTRA repository:

```bash
chmod +x launchers/CoTRA-Ubuntu.sh
./launchers/CoTRA-Ubuntu.sh
```

To stop CoTRA:

```bash
docker stop cotra
```

---

## Run directly with Docker

You can also run CoTRA without the launcher.

Pull the image:

```bash
docker pull ghcr.io/umairseemab/cotra:latest
```

Create persistent folders:

```bash
mkdir -p CoTRA_data CoTRA_results
```

Run CoTRA on macOS or Linux:

```bash
docker run --rm --name cotra \
  -p 3838:3838 \
  -v "$(pwd)/CoTRA_data:/data" \
  -v "$(pwd)/CoTRA_results:/results" \
  ghcr.io/umairseemab/cotra:latest
```

Then open:

```text
http://localhost:3838
```

For Windows, the supplied Windows launcher is recommended because it handles Windows volume paths automatically.

---

## Container data and results

Files placed in:

```text
CoTRA_data/
```

are available inside the container under:

```text
/data
```

CoTRA's default output folder inside the container is:

```text
/results/CoTRA_Results
```

which appears on the host computer as:

```text
CoTRA_results/CoTRA_Results/
```

---

# 2. Native R installation

Use native installation when you want to run CoTRA directly from R or RStudio.

The common installation sequence is:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "UmairSeemab/CoTRA",
  dependencies = FALSE,
  upgrade = "never",
  force = TRUE
)

library(CoTRA)

CoTRA::install_cotra_dependencies(
  ask = FALSE,
  update = FALSE
)

CoTRA::check_cotra_dependencies()

CoTRA::runCoTRA()
```

`CoTRA::install_cotra_dependencies()` installs the CRAN, Bioconductor, and GitHub dependencies required by CoTRA.

---

## Native installation on Windows

### Requirements

Install:

- R >= 4.4.0
- RStudio, optional but recommended

### Using the supplied installer

Clone or download the CoTRA repository.

From the repository folder, run:

```text
install_CoTRA_Windows.bat
```

After installation, open R or RStudio:

```r
library(CoTRA)
CoTRA::runCoTRA()
```

You can also use the common native R installation commands shown above.

---

## Native installation on macOS

CoTRA supports Intel Macs and Apple Silicon Macs.

For Apple Silicon, install the ARM64 version of R.

### Install macOS build tools

Open Terminal:

```bash
xcode-select --install
```

Install Homebrew if it is not already installed:

```text
https://brew.sh
```

Install required system libraries:

```bash
brew install hdf5 pkgconf openssl@3
```

### Run the CoTRA installer

Clone or download the CoTRA repository:

```bash
chmod +x install_CoTRA_MacOS.command
./install_CoTRA_MacOS.command
```

The installer checks the required Apple build tools and Homebrew libraries before installing CoTRA.

After installation:

```r
library(CoTRA)
CoTRA::runCoTRA()
```

---

## Native installation on Ubuntu/Linux

### Requirements

Install R >= 4.4.0.

Install the Linux libraries required by the CoTRA dependency stack:

```bash
sudo apt-get update

sudo apt-get install -y \
  git curl ca-certificates build-essential cmake pkg-config pandoc \
  libcurl4-openssl-dev libssl-dev libxml2-dev libxslt1-dev libgit2-dev \
  libfontconfig1-dev libfreetype6-dev libpng-dev libjpeg-dev libtiff-dev \
  libcairo2-dev libharfbuzz-dev libfribidi-dev libxt-dev \
  libglpk-dev libgmp3-dev libmpfr-dev libgsl-dev libhdf5-dev \
  libudunits2-dev libgdal-dev libgeos-dev libproj-dev \
  libmagick++-dev libpoppler-cpp-dev libv8-dev libicu-dev \
  libbz2-dev liblzma-dev libpcre2-dev zlib1g-dev
```

Clone CoTRA:

```bash
git clone https://github.com/UmairSeemab/CoTRA.git
cd CoTRA
```

Run the installer:

```bash
chmod +x install_CoTRA_Linux.sh
./install_CoTRA_Linux.sh
```

After installation:

```r
library(CoTRA)
CoTRA::runCoTRA()
```

---

# 3. HPC installation with Apptainer

On HPC systems where Docker is unavailable to normal users, use the same CoTRA image through Apptainer.

Pull the image:

```bash
apptainer pull CoTRA.sif docker://ghcr.io/umairseemab/cotra:latest
```

Create input and output folders:

```bash
mkdir -p CoTRA_data CoTRA_results
```

Run CoTRA:

```bash
apptainer run \
  --bind "$PWD/CoTRA_data:/data" \
  --bind "$PWD/CoTRA_results:/results" \
  CoTRA.sif
```

CoTRA listens on port:

```text
3838
```

On a remote HPC system, use the site's supported SSH tunnelling or interactive web-application method to access port 3838 from your local browser.

---

## CSC Roihu

CSC Roihu supports Apptainer containers. The CoTRA Docker/OCI image can therefore be used as an Apptainer SIF image.

The CPU side of Roihu is the normal choice for the current CoTRA workflows.

Move to your CSC project directory, for example:

```bash
cd /scratch/project_XXXXXXX
```

CSC recommends using the node-local temporary directory for the Apptainer cache:

```bash
export APPTAINER_CACHEDIR="$TMPDIR"
export APPTAINER_TMPDIR="$TMPDIR"
```

Pull CoTRA:

```bash
apptainer pull CoTRA.sif docker://ghcr.io/umairseemab/cotra:latest
```

Create persistent folders:

```bash
mkdir -p CoTRA_data CoTRA_results
```

Run:

```bash
apptainer run \
  --bind "$PWD/CoTRA_data:/data" \
  --bind "$PWD/CoTRA_results:/results" \
  CoTRA.sif
```

### Example SLURM job

Adjust the CSC account, time, memory, CPUs, and partition for your project.

```bash
#!/bin/bash
#SBATCH --account=project_XXXXXXX
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

export APPTAINER_CACHEDIR="$TMPDIR"
export APPTAINER_TMPDIR="$TMPDIR"

apptainer run \
  --bind "$PWD/CoTRA_data:/data" \
  --bind "$PWD/CoTRA_results:/results" \
  CoTRA.sif
```

Replace:

```text
project_XXXXXXX
```

with your CSC project number.

Current CSC Roihu systems provide Apptainer directly, so a separate `module load apptainer` step is not normally required.

Use CSC's current SSH tunnelling or interactive web-application instructions to access CoTRA from a compute node.

---

# Output folder

For native installations, choose an output folder from the CoTRA Home page.

If no output folder is selected, CoTRA uses:

```text
~/CoTRA_Results
```

For Docker and Apptainer installations, use the mounted `/results` location so generated files remain available after the container stops.
