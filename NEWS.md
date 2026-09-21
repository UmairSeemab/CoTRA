# CoTRA 1.1.0

## Containerized deployment

- Added Docker-based deployment for CoTRA.
- Added a predefined container environment containing R, CoTRA dependencies, required Linux system libraries, HDF5 support, and Pandoc.
- Added Docker Compose configuration for simplified local deployment.
- Added support for AMD64 and ARM64 container images.
- Added container installation instructions for Windows, macOS, and Linux.
- Added guidance for running the CoTRA OCI/Docker image through Apptainer on compatible HPC systems.
- Added persistent input and output directory mounting for containerized analyses.
- Added automated container build configuration through GitHub Actions.

## Documentation

- Expanded installation documentation for native and containerized deployment.
- Added Docker and HPC/Apptainer usage instructions.
- Clarified storage of input data and analysis results when running CoTRA in containers.

## Software

- Updated CoTRA package version to 1.1.0.
- No changes were made to the core statistical methods underlying the previously reported implementation-validation and computational-benchmarking analyses.
