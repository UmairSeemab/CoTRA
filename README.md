# CoTRA

CoTRA, Comprehensive Toolbox for RNA Sequencing Data Analysis, is an R/Shiny application for bulk RNA-seq and single-cell RNA-seq analysis.

CoTRA can be used in two ways:

1. **Container installation, recommended**: Docker on Windows, macOS, and Ubuntu/Linux, or Apptainer on HPC systems. R and the CoTRA R-package dependency stack are installed inside the image.
2. **Native R installation**: install CoTRA directly into your local R/RStudio environment.

The container route is recommended when you want the most consistent environment across computers and want to avoid repeated local R/Bioconductor dependency installation.

---

## 1. Recommended: run CoTRA in a container

### What the container provides

The CoTRA image is built from a versioned Rocker R image and contains:

- R 4.5.2
- CoTRA
- CoTRA CRAN dependencies
- CoTRA Bioconductor dependencies
- CoTRA GitHub dependencies declared in `R/dependencies.R`
- HDF5 support through `hdf5r` and system HDF5 libraries
- system libraries required by common CoTRA dependencies
- Pandoc
- persistent `/data` and `/results` locations
- an application health check

The GitHub Actions workflow targets:

```text
linux/amd64
linux/arm64
```

This allows the same containerized CoTRA environment to be used through Docker Desktop on Windows and macOS and through Docker Engine on Ubuntu/Linux. Apple Silicon Macs use the `arm64` image when the multi-platform build succeeds.

The first GitHub Actions build should be treated as the validation step for the complete dependency stack on both architectures.

### Important scope

The container covers the CoTRA application and its R/system-library environment. Upstream command-line preprocessing and alignment tools such as FastQC, STAR, HTSeq, and MultiQC are not bundled in this image because they are not part of the current CoTRA Shiny application runtime.

If a future CoTRA module directly invokes one of these tools, add it to the Docker image and validate it in CI.

---

## 2. Container file layout

Add the following files to the repository:

```text
CoTRA/
├── README.md
├── Dockerfile
├── docker-compose.yml
├── .dockerignore
├── docker/
│   ├── install-dependencies.R
│   ├── start-cotra.sh
│   └── healthcheck.R
├── launchers/
│   ├── CoTRA-Windows.bat
│   ├── CoTRA-Windows.ps1
│   ├── CoTRA-macOS.command
│   └── CoTRA-Ubuntu.sh
├── apptainer/
│   └── CoTRA.def
└── .github/
    └── workflows/
        └── docker-build.yml
```

Do not maintain a separate `CONTAINER_README.md`. Container instructions are included here so users have one installation document.

---

## 3. Publish the CoTRA image with GitHub Actions

After the container files are committed, the workflow at:

```text
.github/workflows/docker-build.yml
```

builds and publishes the image to GitHub Container Registry.

The expected image name is:

```text
ghcr.io/umairseemab/cotra:latest
```

A Git tag such as:

```text
v0.1.0
```

also produces a versioned image tag.

After the first successful workflow run, check the GitHub package settings. If anonymous users cannot pull the image, change the package visibility to **Public**.

The workflow has `packages: write` permission and publishes both `linux/amd64` and `linux/arm64` images.

---

## 4. Windows

### Requirements

Install and start Docker Desktop.

No separate R, RStudio, CRAN, or Bioconductor installation is required for the container route.

### Recommended launcher

Download or clone the CoTRA repository and run:

```text
launchers/CoTRA-Windows.ps1
```

Alternatively, double-click:

```text
launchers/CoTRA-Windows.bat
```

The launcher will:

- check that Docker is available
- pull the latest CoTRA image
- create persistent `CoTRA_data` and `CoTRA_results` folders
- start the container
- wait for the health check
- open `http://localhost:3838`

To stop CoTRA:

```powershell
docker stop cotra
```

---

## 5. macOS

### Requirements

Install and start Docker Desktop.

The GitHub Actions workflow targets both Intel/AMD64 and ARM64 container architectures.

Make the launcher executable once after cloning if needed:

```bash
chmod +x launchers/CoTRA-macOS.command
```

Run:

```bash
./launchers/CoTRA-macOS.command
```

The launcher opens:

```text
http://localhost:3838
```

To stop CoTRA:

```bash
docker stop cotra
```

---

## 6. Ubuntu/Linux

### Requirements

Install Docker Engine and ensure your user can run Docker.

Make the launcher executable once after cloning:

```bash
chmod +x launchers/CoTRA-Ubuntu.sh
```

Run:

```bash
./launchers/CoTRA-Ubuntu.sh
```

Then open:

```text
http://localhost:3838
```

To stop CoTRA:

```bash
docker stop cotra
```

---

## 7. Run directly with Docker

After the GHCR image has been published:

```bash
docker pull ghcr.io/umairseemab/cotra:latest
```

Create local data and result directories:

```bash
mkdir -p CoTRA_data CoTRA_results
```

Run CoTRA:

```bash
docker run --rm --name cotra \
  -p 3838:3838 \
  -v "$(pwd)/CoTRA_data:/data" \
  -v "$(pwd)/CoTRA_results:/results" \
  ghcr.io/umairseemab/cotra:latest
```

Open:

```text
http://localhost:3838
```

On Windows PowerShell, the provided Windows launcher is simpler than manually translating the volume paths.

---

## 8. Docker Compose

From the repository root:

```bash
docker compose pull
docker compose up -d
```

Open:

```text
http://localhost:3838
```

View logs:

```bash
docker compose logs -f cotra
```

Stop CoTRA:

```bash
docker compose down
```

---

## 9. Data and output folders

The container uses two mounted locations:

```text
/data
/results
```

They map to these folders in the repository directory:

```text
CoTRA_data/
CoTRA_results/
```

You can place input datasets in:

```text
CoTRA_data/
```

Inside the container they are available under:

```text
/data
```

The container launcher sets its home directory to `/results`. CoTRA's existing default output path `~/CoTRA_Results` therefore becomes:

```text
/results/CoTRA_Results
```

On the host computer, those default results appear under:

```text
CoTRA_results/CoTRA_Results/
```

You can still select another writable output directory from the CoTRA interface.

---

## 10. Build the image locally

You can build CoTRA before publishing it to GHCR.

From the repository root:

```bash
docker build -t cotra:local .
```

Run the local image:

```bash
docker run --rm --name cotra \
  -p 3838:3838 \
  -v "$(pwd)/CoTRA_data:/data" \
  -v "$(pwd)/CoTRA_results:/results" \
  cotra:local
```

Check its status:

```bash
docker ps
```

View logs:

```bash
docker logs cotra
```

---

## 11. Reproducibility and versioning

The Dockerfile pins the base R version:

```text
R 4.5.2
```

A built container image preserves the exact packages installed in that image. For reproducible releases, use versioned image tags rather than relying only on `latest`.

Example:

```text
ghcr.io/umairseemab/cotra:v0.1.0
```

The repository does not currently contain a complete `renv.lock` for every CoTRA dependency. Therefore, rebuilding the Dockerfile at a later date can still pick up newer package versions from CRAN, Bioconductor, or GitHub even though the R base version is fixed.

For strict source-level reproducibility, add a tested lockfile or pin individual dependency versions/commits in a future release. The published versioned container image itself remains the most direct record of the tested software environment.

---

## 12. Apptainer for HPC systems

The same OCI image can be used on HPC systems that provide Apptainer.

Pull the image:

```bash
apptainer pull CoTRA.sif docker://ghcr.io/umairseemab/cotra:latest
```

Create host folders:

```bash
mkdir -p CoTRA_data CoTRA_results
```

Run with persistent mounts:

```bash
apptainer run \
  --bind "$PWD/CoTRA_data:/data" \
  --bind "$PWD/CoTRA_results:/results" \
  CoTRA.sif
```

For a remote HPC compute node, you may also need SSH tunnelling or the site-specific method for exposing port `3838` to your browser.

### Example SLURM job

```bash
#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

module load apptainer

mkdir -p CoTRA_data CoTRA_results

apptainer run \
  --bind "$PWD/CoTRA_data:/data" \
  --bind "$PWD/CoTRA_results:/results" \
  CoTRA.sif
```

For CSC/Roihu, use the Apptainer module and storage paths appropriate for your CSC project and follow CSC's current port-forwarding guidance for interactive web applications.

---

# Native R installation

Container installation is recommended for users who want a predefined environment. Native installation remains available for developers and users who prefer direct R/RStudio access.

## Requirements

CoTRA currently declares:

```text
R >= 4.4.0
```

Install CoTRA from GitHub:

```r
install.packages("remotes")
remotes::install_github(
  "UmairSeemab/CoTRA",
  dependencies = TRUE,
  upgrade = "never"
)
```

Install the analysis dependency stack:

```r
library(CoTRA)
CoTRA::install_cotra_dependencies(ask = FALSE)
```

For 10x Genomics HDF5 input, ensure `hdf5r` is installed:

```r
install.packages("hdf5r")
```

Run CoTRA:

```r
library(CoTRA)
CoTRA::runCoTRA()
```

`runCoTRA()` creates a temporary writable copy of the Shiny application so CoTRA does not write generated output inside the installed R package directory.

---

## Existing native installation helper files

The repository also contains:

```text
install_cotra_packages.R
install_CoTRA_Windows.bat
install_CoTRA_MacOS.command
install_CoTRA_Linux.sh
```

These are native R installation helpers. They are separate from the Docker launchers under `launchers/`.

Use the container launchers when you want to avoid installing the R dependency stack on the host computer.

---

## Native installation without administrator rights on Windows

A user-writable R library can be used when administrator rights are unavailable.

Example:

```r
r_version <- paste0(
  R.version$major,
  ".",
  strsplit(R.version$minor, "\\.")[[1]][1]
)

local_appdata <- Sys.getenv("LOCALAPPDATA")
if (!nzchar(local_appdata)) {
  local_appdata <- path.expand("~")
}

user_lib <- file.path(
  local_appdata,
  "R",
  "win-library",
  r_version
)

dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(R_LIBS_USER = user_lib)
.libPaths(unique(c(user_lib, .libPaths())))

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  timeout = 2000
)

install.packages(
  c("BiocManager", "remotes", "hdf5r"),
  lib = user_lib,
  dependencies = TRUE
)

remotes::install_github(
  "UmairSeemab/CoTRA",
  lib = user_lib,
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE
)

library(CoTRA, lib.loc = user_lib)
CoTRA::install_cotra_dependencies(ask = FALSE)
CoTRA::runCoTRA()
```

If compiled packages require system libraries or build tools, the container installation is preferable because those system dependencies are installed inside the image.

---

## Native installation on CSC/Roihu

If you prefer a native R library instead of Apptainer, create a project-specific library tied to the active R major/minor version.

```r
project_directory <- "/projappl/project_XXXXXXX"

r_major_minor <- paste(
  R.version$major,
  strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1],
  sep = "."
)

libpath <- file.path(
  project_directory,
  paste0("CoTRA_Rlibs_R", r_major_minor)
)

dir.create(libpath, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(libpath, .libPaths())))
Sys.setenv(R_LIBS_USER = libpath)

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  timeout = 2000
)

install.packages(
  c("BiocManager", "remotes", "hdf5r"),
  lib = libpath,
  dependencies = TRUE
)

remotes::install_github(
  "UmairSeemab/CoTRA",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE,
  lib = libpath
)

library(CoTRA, lib.loc = libpath)
CoTRA::install_cotra_dependencies(ask = FALSE)
```

Run later with:

```r
project_directory <- "/projappl/project_XXXXXXX"

r_major_minor <- paste(
  R.version$major,
  strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1],
  sep = "."
)

libpath <- file.path(
  project_directory,
  paste0("CoTRA_Rlibs_R", r_major_minor)
)

.libPaths(unique(c(libpath, .libPaths())))
Sys.setenv(R_LIBS_USER = libpath)

library(CoTRA, lib.loc = libpath)
CoTRA::runCoTRA()
```

Replace `project_XXXXXXX` with your CSC project directory.

---

# Output folder

For a native installation, open the CoTRA Home page and select an output folder. CoTRA stores generated reports, figures, tables, ZIP files, and session files in that directory.

If no output folder is selected, native CoTRA uses:

```text
~/CoTRA_Results
```

For the supplied Docker container, this default is mapped into the persistent host `CoTRA_results` directory as described above.

---

# Bioconductor preparation

CoTRA keeps reusable package functions under `R/` and the Shiny application under `inst/app/`.

Before a Bioconductor-oriented package submission, run the appropriate package checks, for example:

```r
devtools::check()
BiocCheck::BiocCheck()
```

---

# Container troubleshooting

## Port 3838 is already in use

Run on another host port:

```bash
docker run --rm --name cotra \
  -p 3839:3838 \
  ghcr.io/umairseemab/cotra:latest
```

Then open:

```text
http://localhost:3839
```

## View container logs

```bash
docker logs cotra
```

## Check container health

```bash
docker inspect --format '{{.State.Health.Status}}' cotra
```

## Remove an old container

```bash
docker rm -f cotra
```

## Refresh the image

```bash
docker pull ghcr.io/umairseemab/cotra:latest
```

## Browser-based report rendering

The image includes Pandoc, but it does not currently install Google Chrome or Chromium. If a report path specifically requires `webshot2`, `pagedown`, or another headless-Chrome operation, validate that workflow separately before treating browser-based PDF capture as supported across both AMD64 and ARM64 images.

---

# Development notes

The Dockerfile installs dependencies before launching the application and fails the image build if the dependency verification step reports missing declared CoTRA dependencies. This moves most package installation failures from the end user's computer into CI, where they can be fixed once for all users.

When changing CoTRA dependencies:

1. Update `DESCRIPTION` where appropriate.
2. Update `R/dependencies.R`.
3. Rebuild the container.
4. Check the GitHub Actions build for both architectures.
5. Test representative bulk RNA-seq and scRNA-seq workflows.
6. Publish a versioned container tag for a release.

