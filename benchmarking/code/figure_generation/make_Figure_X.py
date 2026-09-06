#!/usr/bin/env python3
"""
Reproduce CoTRA Main Figure X from the supplied CSV files.

Outputs
-------
Figure_XA.png / .svg  Bulk synthetic runtime
Figure_XB.png / .svg  Bulk synthetic peak RAM
Figure_XC.png / .svg  Synthetic scRNA combined runtime
Figure_XD.png / .svg  Synthetic scRNA peak RAM
Main_Figure_X.png      Combined four-panel figure, 300 dpi
Main_Figure_X.tif      Combined four-panel TIFF, LZW compressed, 300 dpi
Main_Figure_X.svg      Combined SVG wrapper
                       (contains the four high-resolution panels)

PLOS Computational Biology-oriented settings
---------------------------------------------
- 300 dpi
- Combined canvas: 2070 x 1770 pixels
- White RGB background
- TIFF LZW compression
- Text target: 10 pt
- Default font requested in this script: Arial
- Panel labels A-D
- No figure title/caption embedded in the figure

PLOS currently lists Arial, Times, or Symbol as accepted figure fonts.
If you explicitly want Arial Black, run:
    python make_Figure_X.py --font "Arial Black"

Note: Arial/Arial Black must be installed on the machine where this
script is executed. If unavailable, the script warns and uses a
compatible sans-serif fallback.
"""

import argparse
import base64
import csv
import io
import os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager
from PIL import Image


DPI = 300

# Combined target size:
# 6.90 x 5.90 inches at 300 dpi = 2070 x 1770 pixels.
PANEL_W_IN = 3.45
PANEL_H_IN = 2.95
COMBINED_W_PX = 2070
COMBINED_H_PX = 1770


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def choose_font(requested):
    """Return requested font if installed, otherwise an allowed-like fallback."""
    available = {f.name for f in font_manager.fontManager.ttflist}
    if requested in available:
        return requested

    fallbacks = ["Arial", "Arimo", "Liberation Sans", "DejaVu Sans"]
    for name in fallbacks:
        if name in available:
            print(
                f"WARNING: requested font '{requested}' is not installed. "
                f"Using '{name}' instead."
            )
            return name

    print(
        f"WARNING: requested font '{requested}' is not installed. "
        "Using matplotlib's default sans-serif font."
    )
    return "sans-serif"


def set_plot_defaults(font_name, bold=False):
    plt.rcParams.update({
        "font.family": font_name,
        "font.size": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "svg.fonttype": "none",
    })
    return "black" if bold else "normal"


def finish_axis(ax, panel_letter, weight):
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", which="major", length=3, width=0.8)
    ax.margins(x=0.05)

    ax.text(
        -0.13, 1.04, panel_letter,
        transform=ax.transAxes,
        fontsize=10,
        fontweight=weight,
        ha="left",
        va="bottom",
    )

    ax.xaxis.label.set_fontweight(weight)
    ax.yaxis.label.set_fontweight(weight)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight(weight)

    leg = ax.get_legend()
    if leg is not None:
        for txt in leg.get_texts():
            txt.set_fontweight(weight)


def save_panel(fig, out_base):
    fig.tight_layout(pad=0.6)
    fig.savefig(
        out_base.with_suffix(".png"),
        dpi=DPI,
        facecolor="white",
        edgecolor="white",
    )
    fig.savefig(
        out_base.with_suffix(".svg"),
        facecolor="white",
        edgecolor="white",
    )
    plt.close(fig)


def panel_bulk_runtime(rows, out_base, weight):
    fig = plt.figure(figsize=(PANEL_W_IN, PANEL_H_IN), dpi=DPI)
    ax = fig.add_subplot(111)

    for method in ["DESeq2", "edgeR"]:
        rr = sorted(
            [r for r in rows if r["method"] == method],
            key=lambda r: int(r["samples"]),
        )
        x = [int(r["samples"]) for r in rr]
        y = [float(r["runtime_median_s"]) for r in rr]
        lo = [
            float(r["runtime_median_s"]) - float(r["runtime_q1_s"])
            for r in rr
        ]
        hi = [
            float(r["runtime_q3_s"]) - float(r["runtime_median_s"])
            for r in rr
        ]

        ax.errorbar(
            x, y,
            yerr=[lo, hi],
            marker="o",
            capsize=2.5,
            linewidth=1.2,
            label=method,
        )

    ax.set_xlabel("Number of samples")
    ax.set_ylabel("Median runtime (s)")
    ax.legend(frameon=False, loc="upper left")
    finish_axis(ax, "A", weight)
    save_panel(fig, out_base)


def panel_bulk_ram(rows, out_base, weight):
    fig = plt.figure(figsize=(PANEL_W_IN, PANEL_H_IN), dpi=DPI)
    ax = fig.add_subplot(111)

    for method in ["DESeq2", "edgeR"]:
        rr = sorted(
            [r for r in rows if r["method"] == method],
            key=lambda r: int(r["samples"]),
        )
        x = [int(r["samples"]) for r in rr]
        y = [float(r["peak_ram_median_gb"]) for r in rr]
        lo = [
            float(r["peak_ram_median_gb"]) - float(r["peak_ram_q1_gb"])
            for r in rr
        ]
        hi = [
            float(r["peak_ram_q3_gb"]) - float(r["peak_ram_median_gb"])
            for r in rr
        ]

        ax.errorbar(
            x, y,
            yerr=[lo, hi],
            marker="o",
            capsize=2.5,
            linewidth=1.2,
            label=method,
        )

    ax.set_xlabel("Number of samples")
    ax.set_ylabel("Peak RAM (GB)")
    ax.legend(frameon=False, loc="upper left")
    finish_axis(ax, "B", weight)
    save_panel(fig, out_base)


def panel_scrna_runtime(rows, out_base, weight):
    fig = plt.figure(figsize=(PANEL_W_IN, PANEL_H_IN), dpi=DPI)
    ax = fig.add_subplot(111)

    rr = sorted(rows, key=lambda r: int(r["cells"]))
    x = [int(r["cells"]) for r in rr]
    y = [float(r["combined_runtime_median_s"]) for r in rr]
    lo = [
        float(r["combined_runtime_median_s"])
        - float(r["combined_runtime_q1_s"])
        for r in rr
    ]
    hi = [
        float(r["combined_runtime_q3_s"])
        - float(r["combined_runtime_median_s"])
        for r in rr
    ]

    ax.errorbar(
        x, y,
        yerr=[lo, hi],
        marker="o",
        capsize=2.5,
        linewidth=1.2,
    )

    ax.set_xlabel("Number of cells")
    ax.set_ylabel("Median combined runtime (s)")
    finish_axis(ax, "C", weight)
    save_panel(fig, out_base)


def panel_scrna_ram(rows, out_base, weight):
    fig = plt.figure(figsize=(PANEL_W_IN, PANEL_H_IN), dpi=DPI)
    ax = fig.add_subplot(111)

    rr = sorted(rows, key=lambda r: int(r["cells"]))
    x = [int(r["cells"]) for r in rr]
    y = [float(r["peak_ram_median_gb"]) for r in rr]
    lo = [
        float(r["peak_ram_median_gb"]) - float(r["peak_ram_q1_gb"])
        for r in rr
    ]
    hi = [
        float(r["peak_ram_q3_gb"]) - float(r["peak_ram_median_gb"])
        for r in rr
    ]

    ax.errorbar(
        x, y,
        yerr=[lo, hi],
        marker="o",
        capsize=2.5,
        linewidth=1.2,
    )

    ax.set_xlabel("Number of cells")
    ax.set_ylabel("Peak RAM (GB)")
    finish_axis(ax, "D", weight)
    save_panel(fig, out_base)


def combine_png_panels(panel_paths, out_png, out_tif):
    imgs = [Image.open(p).convert("RGB") for p in panel_paths]

    width = max(im.width for im in imgs)
    height = max(im.height for im in imgs)

    canvas = Image.new(
        "RGB",
        (width * 2, height * 2),
        "white",
    )

    for i, im in enumerate(imgs):
        x = (i % 2) * width
        y = (i // 2) * height
        canvas.paste(im, (x, y))

    if canvas.size != (COMBINED_W_PX, COMBINED_H_PX):
        canvas = canvas.resize(
            (COMBINED_W_PX, COMBINED_H_PX),
            Image.Resampling.LANCZOS,
        )

    canvas.save(
        out_png,
        dpi=(DPI, DPI),
        optimize=True,
    )

    canvas.save(
        out_tif,
        dpi=(DPI, DPI),
        compression="tiff_lzw",
    )


def make_svg_wrapper(panel_pngs, out_svg):
    """
    Create a single-page SVG containing the four high-resolution panels.
    The TIFF remains the recommended PLOS submission file.
    """
    images = [Image.open(p) for p in panel_pngs]
    w = max(im.width for im in images)
    h = max(im.height for im in images)

    total_w = w * 2
    total_h = h * 2

    positions = [
        (0, 0),
        (w, 0),
        (0, h),
        (w, h),
    ]

    lines = [
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{total_w}" height="{total_h}" '
            f'viewBox="0 0 {total_w} {total_h}">'
        ),
        '<rect width="100%" height="100%" fill="white"/>',
    ]

    for (x, y), path in zip(positions, panel_pngs):
        encoded = base64.b64encode(path.read_bytes()).decode("ascii")
        uri = "data:image/png;base64," + encoded
        lines.append(
            f'<image x="{x}" y="{y}" '
            f'width="{w}" height="{h}" href="{uri}"/>'
        )

    lines.append("</svg>")
    out_svg.write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bulk",
        default="Figure_X_bulk_synthetic_data.csv",
        help="Bulk synthetic plotting CSV.",
    )
    parser.add_argument(
        "--scrna",
        default="Figure_X_scRNA_synthetic_data.csv",
        help="Synthetic scRNA plotting CSV.",
    )
    parser.add_argument(
        "--outdir",
        default="Figure_X_output",
        help="Output directory.",
    )
    parser.add_argument(
        "--font",
        default="Arial",
        help='Figure font. PLOS-safe default is "Arial".',
    )
    parser.add_argument(
        "--bold",
        action="store_true",
        help=(
            "Use a heavy/bold appearance. "
            'For the requested look, use --font "Arial Black" --bold.'
        ),
    )
    args = parser.parse_args()

    bulk = read_csv(args.bulk)
    scrna = read_csv(args.scrna)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    font = choose_font(args.font)
    weight = set_plot_defaults(font, bold=args.bold)

    bases = [
        outdir / "Figure_XA",
        outdir / "Figure_XB",
        outdir / "Figure_XC",
        outdir / "Figure_XD",
    ]

    panel_bulk_runtime(bulk, bases[0], weight)
    panel_bulk_ram(bulk, bases[1], weight)
    panel_scrna_runtime(scrna, bases[2], weight)
    panel_scrna_ram(scrna, bases[3], weight)

    panel_pngs = [p.with_suffix(".png") for p in bases]

    combine_png_panels(
        panel_pngs,
        outdir / "Main_Figure_X.png",
        outdir / "Main_Figure_X.tif",
    )

    make_svg_wrapper(
        panel_pngs,
        outdir / "Main_Figure_X.svg",
    )

    print("Figure X created successfully.")
    print(f"Font used: {font}")
    print(f"Output directory: {outdir.resolve()}")
    print("Recommended PLOS submission file: Main_Figure_X.tif")


if __name__ == "__main__":
    main()
