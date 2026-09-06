# Generates partition_schemes.pdf: the two partitions shown in the
# Grothendieck exposition of the AI companion paper, in the same visual
# style as the scheme figures of the theory paper (orange / near-black
# regions with slight transparency, gray grid drawn on top, white dashed
# baseline at z2 = 0).
#
# Run with:  uv run --with matplotlib --with numpy python make_partition_schemes.py

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.size": 12,
    "mathtext.fontset": "dejavusans",
})

ORANGE = "#ff7f0e"
DARK = "#262626"
GRID = "#8fa3b8"

x = np.linspace(-4, 4, 1601)
theta = 0.3

panels = [
    (np.zeros_like(x), r"$\mathrm{sgn}(z_2)$"),
    (theta * (x**3 - 3 * x), r"$\mathrm{sgn}\left(z_2 - \vartheta\,(z_1^3 - 3z_1)\right)$"),
]

fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.0))
for ax, (P, title) in zip(axes, panels):
    ax.fill_between(x, P, 4.5, color=ORANGE, alpha=0.88, lw=0, zorder=1)
    ax.fill_between(x, -4.5, P, color=DARK, alpha=0.92, lw=0, zorder=1)
    ax.grid(True, color=GRID, lw=0.9, alpha=0.55, zorder=2)
    ax.set_axisbelow(False)
    ax.axhline(0.0, color="white", lw=2.0, ls=(0, (6, 4)), zorder=3)
    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    ax.set_xticks(range(-4, 5))
    ax.set_yticks(range(-4, 5))
    ax.set_xlabel(r"$z_1$", fontsize=14)
    ax.set_ylabel(r"$z_2$", fontsize=14)
    ax.set_title(title, fontsize=14, pad=10)
    ax.set_aspect("equal")
    ax.text(0, 2.7, r"$+1$", color="white", fontsize=17, ha="center", va="center", zorder=4)
    ax.text(0, -2.7, r"$-1$", color="white", fontsize=17, ha="center", va="center", zorder=4)

fig.tight_layout()
fig.savefig("partition_schemes.pdf", bbox_inches="tight")
fig.savefig("/tmp/partition_schemes.png", dpi=150, bbox_inches="tight")
print("done")
