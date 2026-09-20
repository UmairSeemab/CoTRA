# CoTRA container files: GitHub upload locations

Copy these files into the existing `UmairSeemab/CoTRA` repository using exactly these paths.

```text
CoTRA/
├── README.md                              ← replace current README.md
├── Dockerfile                             ← repository root
├── docker-compose.yml                     ← repository root
├── .dockerignore                          ← repository root
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

Do not create a separate `CONTAINER_README.md`. The container documentation is already merged into the supplied root `README.md`.

## After uploading

1. Commit the files to a branch or `main`.
2. Open the repository's **Actions** tab and run **Build CoTRA container** if it did not start automatically.
3. Confirm that the workflow completes for both `linux/amd64` and `linux/arm64`.
4. Open the repository/package page and verify that `ghcr.io/umairseemab/cotra:latest` exists.
5. If anonymous users cannot pull the image, change the GHCR package visibility to **Public**.
6. Test one launcher on each platform you support.
7. For a stable release, create a Git tag such as `v0.1.0`. The workflow will publish a matching versioned container tag.

## Executable permission

GitHub's web uploader does not reliably preserve executable bits. After cloning on macOS/Linux, users can run:

```bash
chmod +x launchers/CoTRA-macOS.command
chmod +x launchers/CoTRA-Ubuntu.sh
chmod +x docker/start-cotra.sh
```

If you commit from Linux/macOS, set these executable bits before `git add` so Git stores them in the repository.
