"""Shared plotting style and output paths for the figure scripts.

Palette
-------
Categorical hues are taken in fixed slot order and never cycled; three of the
light-mode slots sit below 3:1 against the chart surface, so every chart here
carries a legend and, where there is room, a direct label as well.  The set was
validated for colour-vision deficiency separation on the adjacent-pair list
(worst adjacent CVD dE 9.1, normal-vision dE 19.6).

Ordered quantities (a sweep over ADC bits, over MASH order) use the single-hue
blue ramp instead, starting no lighter than step 250 so the lightest line still
reads against the surface.
"""

from __future__ import annotations

import pathlib

import matplotlib as mpl
import matplotlib.pyplot as plt

# ------------------------------------------------------------------ paths --
ROOT = pathlib.Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

# ------------------------------------------------------------- the palette --
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

#: categorical slots, in fixed order -- index by role, never cycle
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
          "#4a3aa7", "#e34948"]

#: single-hue ordinal ramp (blue 250..700); lightest step still clears 2:1
ORDINAL = ["#86b6ef", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#184f95",
           "#0d366b"]


def ordinal(n: int) -> list[str]:
    """``n`` evenly spaced steps of the ordinal ramp, light to dark."""
    if n <= 1:
        return [ORDINAL[len(ORDINAL) // 2]]
    idx = [round(i * (len(ORDINAL) - 1) / (n - 1)) for i in range(n)]
    return [ORDINAL[i] for i in idx]


def use_style() -> None:
    """Apply the house style.  Call once at the top of a script."""
    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "sans-serif",
        "font.sans-serif": ["Segoe UI", "DejaVu Sans", "sans-serif"],
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.titleweight": "semibold",
        "axes.titlelocation": "left",
        "axes.titlepad": 10,
        "axes.labelsize": 9,
        "axes.labelcolor": INK_2,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.7,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK_2,
        "ytick.labelcolor": INK_2,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "lines.linewidth": 2.0,
        "lines.markersize": 8,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "legend.labelcolor": INK_2,
        "figure.dpi": 130,
        "savefig.dpi": 160,
        "savefig.bbox": "tight",
    })


def title(ax, text: str, subtitle: str | None = None) -> None:
    """Left-aligned title with an optional muted second line.

    The subtitle sits in the gap the extra title pad opens up, so the two never
    land on top of each other.
    """
    if subtitle:
        ax.set_title(text, color=INK, pad=24)
        ax.text(0.0, 1.015, subtitle, transform=ax.transAxes, fontsize=8.5,
                color=MUTED, va="bottom", ha="left")
    else:
        ax.set_title(text, color=INK)


def save(fig, name: str) -> pathlib.Path:
    """Write ``name`` into ``results/`` and report where it went."""
    RESULTS.mkdir(exist_ok=True)
    path = RESULTS / name
    fig.savefig(path)
    plt.close(fig)
    print(f"  wrote {path.relative_to(ROOT)}")
    return path


def annotate(ax, x, y, text, color, dx=6, dy=0, ha="left", va="center"):
    """Direct label in the series colour's own ink weight.

    Values and labels stay in text ink; the line beside them carries identity.
    Colour is used only for the small leader dot.
    """
    ax.annotate(text, xy=(x, y), xytext=(dx, dy), textcoords="offset points",
                color=INK_2, fontsize=8.5, ha=ha, va=va)
    ax.plot([x], [y], marker="o", markersize=5, color=color,
            markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=5,
            linestyle="none")
