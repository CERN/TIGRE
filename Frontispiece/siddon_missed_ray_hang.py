# -*- coding: utf-8 -*-
"""The geometry that makes Siddon_projection.cu never return.

Top: the x-y plane in the kernel's voxel coordinates, drawn with y HORIZONTAL
so the ray (which runs exactly along -Y) lies across the page and the two
X-planes are horizontal lines. Real captured numbers for the ray that hangs
(proj 53, row 619, col 313 of the test geometry). Bottom: why an odd
detector dimension supplies the first zero direction component.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
BLUE, BLUE_L = "#2a78d6", "#cde2fb"
CRIT, CRIT_L, GOOD, GREY = "#d03b3b", "#fbeaea", "#0ca30c", "#9a9891"

# captured on the GPU, voxel units. nVoxel is (Z,Y,X)-ordered: X=64, Y=128.
NX, NY = 64.0, 128.0
S_X, S_Y = 280.44769, 497.73584          # source
P_X, P_Y = 280.44794, -225.03015         # detector pixel

fig = plt.figure(figsize=(13.5, 9.6), facecolor=SURF)
gs = fig.add_gridspec(2, 1, height_ratios=[2.35, 1.0], hspace=0.30,
                      left=0.055, right=0.985, top=0.885, bottom=0.04)

# ================================================================ panel A
# horizontal axis = y (voxel), vertical axis = x (voxel)
ax = fig.add_subplot(gs[0])
ax.set_facecolor(SURF)
ax.set_aspect("equal")
ax.set_xlim(-300, 600)
ax.set_ylim(-75, 345)

# volume and the two X-planes Siddon measures distance to
ax.add_patch(Rectangle((0, 0), NY, NX, facecolor=BLUE_L, edgecolor=BLUE,
                       lw=1.4, zorder=2))
ax.text(NY / 2, 49, "volume\nx in [0, 64],  y in [0, 128]", ha="center",
        va="center", fontsize=8.6, color=INK, zorder=3)
for xp, lab in ((0.0, "X-plane  x = 0"), (NX, "X-plane  x = nVoxelX = 64")):
    ax.plot([-300, 600], [xp, xp], color=BLUE, lw=0.9, ls=(0, (4, 3)),
            zorder=1, alpha=0.85)
    ax.text(-296, xp + 4, lab, ha="left", va="bottom", fontsize=8.5,
            color=BLUE)
ax.text(-292, -44, r"Siddon: $a_{xm}, a_{xM}$ = ray parameter where it "
        "crosses each X-plane  =  (plane $-$ source.x) / ray.x",
        ha="left", va="top", fontsize=8.5, color=BLUE)

# --- the ordinary miss: same source, slight tilt -> finite, harmless
ax.plot([S_Y, P_Y], [S_X, 222.0], color=GREY, lw=1.3, ls=(0, (2, 3)),
        zorder=3)
ax.text(-292, 214, "ordinary miss:  ray.x $\\neq$ 0, quotients finite,\n"
        "Np saturates to 0 -- harmless (257 million of these)",
        ha="left", va="top", fontsize=8.5, color=INK2)

# --- the ray that hangs: exactly along -Y at x = 280.45
ax.plot([S_Y, P_Y], [S_X, P_X], color=CRIT, lw=2.4, zorder=5)
ax.plot(S_Y, S_X, marker="o", ms=9, mfc=CRIT, mec=SURF, mew=1.5, zorder=6)
ax.plot(P_Y, P_X, marker="s", ms=7.5, mfc=SURF, mec=CRIT, mew=1.7, zorder=6)
ax.text(S_Y + 10, S_X + 6, "source  (x 280.45, y 497.74)", fontsize=9,
        color=INK, va="bottom", ha="left")
ax.text(P_Y - 4, P_X - 9, "detector pixel\n(x 280.45, y $-$225.03)",
        fontsize=9, color=INK, va="top", ha="left")
ax.text(200, 300,
        "THE RAY THAT HANGS  --  exactly along $-$Y\n"
        "ray.x = 280.44794 $-$ 280.44769 = 0.00025 < eps  $\\to$  0      "
        "ray.z = 63.999996 $-$ 64 = $-$4e$-$6 < eps  $\\to$  0\n"
        "runs 216 voxels above the volume and never enters it",
        fontsize=9, color=CRIT, va="bottom", ha="center",
        bbox=dict(boxstyle="round,pad=0.4", fc=SURF, ec=CRIT, lw=0.9))

# both quotients point the SAME way: both planes are below the source
for yy, xp, lab in ((400, NX, r"$(64-280.4)\,/\,0 = -\infty$"),
                    (330, 0.0, r"$(0-280.4)\,/\,0 = -\infty$")):
    ar = FancyArrowPatch((yy, S_X), (yy, xp), arrowstyle="-|>",
                         mutation_scale=13, color=CRIT, lw=1.4, zorder=4)
    ax.add_patch(ar)
    ax.text(yy + 7, (S_X + xp) / 2 + 20, lab, ha="left", va="center",
            fontsize=9.2, color=CRIT, rotation=90)
ax.text(70, 150,
        "source.x = 280.4 is OUTSIDE [0, 64]: both X-planes lie on the\n"
        "SAME side, so fminf and fmaxf collapse onto one value,\n"
        r"$a_{xm} = a_{xM} = -\infty$  $\Rightarrow$  $a_M = \min(a_{xM}, a_{yM}, a_{zM}) = -\infty$" "\n"
        r"$\Rightarrow$  $a_M \cdot$ ray.z $= (-\infty)\cdot 0 =$ NaN  $\Rightarrow$  Np = NaN" "\n"
        r"$\Rightarrow$  (unsigned long) NaN $= 2^{63}$ loop iterations",
        ha="center", va="center", fontsize=9.2, color=INK,
        bbox=dict(boxstyle="round,pad=0.5", fc=CRIT_L, ec=CRIT, lw=0.9))

# --- the safe contrast: source INSIDE the volume's x-extent
G_X = 32.0
ax.plot([S_Y, P_Y], [G_X, G_X], color=GOOD, lw=1.8, zorder=4)
ax.plot(S_Y, G_X, marker="o", ms=8, mfc=GOOD, mec=SURF, mew=1.5, zorder=6)
ax.text(S_Y + 10, G_X - 6, "same ray, but with\nsource.x = 32, inside [0, 64]",
        fontsize=8.6, color=GOOD, va="top", ha="left")
for yy, xp, lab in ((470, NX, r"$+\infty$"), (470, 0.0, r"$-\infty$")):
    ar = FancyArrowPatch((yy, G_X), (yy, xp), arrowstyle="-|>",
                         mutation_scale=11, color=GOOD, lw=1.3, zorder=5)
    ax.add_patch(ar)
    ax.text(yy + 6, xp + (5 if xp else -5), lab, ha="left",
            va="bottom" if xp else "top", fontsize=9, color=GOOD)
ax.text(-292, -22,
        "planes on OPPOSITE sides: $a_{xm} = -\\infty$, $a_{xM} = +\\infty$ "
        "stay apart, $a_M$ comes from y or z and is finite.  "
        "(And this ray hits the volume anyway.)",
        ha="left", va="top", fontsize=8.6, color=GOOD)

for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(colors=INK2, labelsize=8, length=0)
ax.set_xlabel("y  (voxel units)   --   the ray runs along $-$Y, right to left",
              color=INK2, fontsize=9)
ax.set_ylabel("x  (voxel units)", color=INK2, fontsize=9)
ax.grid(True, color="#e6e5e1", lw=0.6)
ax.set_title("A.  The x-y plane in kernel voxel coordinates  --  real numbers "
             "captured from the ray at proj 53, row 619, col 313",
             loc="left", fontsize=10.5, color=INK, pad=8)

# ================================================================ panel B
axb = fig.add_subplot(gs[1])
axb.set_facecolor(SURF)
axb.set_aspect("equal")
axb.set_xlim(-0.6, 23.4)
axb.set_ylim(-1.2, 7.9)
axb.axis("off")
axb.set_title("B.  Why ODD detector dimensions supply the first zero "
              "component", loc="left", fontsize=10.5, color=INK, pad=8)


def grid(x0, y0, n, label, odd):
    for i in range(n):
        for j in range(n):
            on_row = odd and (j == n // 2)
            on_col = odd and (i == n // 2)
            fc = "#ffffff"
            if on_row or on_col:
                fc = CRIT_L
            if on_row and on_col:
                fc = "#f3bdbd"
            axb.add_patch(Rectangle((x0 + i, y0 + j), 1, 1, facecolor=fc,
                                    edgecolor="#bdbcb6", lw=0.7, zorder=2))
    cx, cy = x0 + n / 2.0, y0 + n / 2.0      # optical axis through centre
    axb.plot([x0 - 0.35, x0 + n + 0.35], [cy, cy], color=BLUE, lw=1.1,
             ls=(0, (4, 3)), zorder=3)
    axb.plot([cx, cx], [y0 - 0.35, y0 + n + 0.35], color=BLUE, lw=1.1,
             ls=(0, (4, 3)), zorder=3)
    axb.plot(cx, cy, marker="+", ms=11, mew=1.6, color=BLUE, zorder=4)
    axb.text(x0 + n / 2.0, 7.55, label, ha="center", va="top", fontsize=10,
             color=INK)
    return cx, cy


# align the two grids on a common axis height
cx1, cy1 = grid(0.0, 0.0, 7, "ODD  (7 x 7)", True)
cx2, cy2 = grid(9.0, 0.5, 6, "EVEN  (6 x 6)", False)

axb.text(8.0, 6.6, "optical axis\n(dashed)", ha="center", va="center",
         fontsize=8.3, color=BLUE)

# the half-pixel offset on the even grid, drawn outside the crosshair
axb.annotate("", xy=(cx2 + 3.55, cy2), xytext=(cx2 + 3.55, cy2 + 0.5),
             arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.0,
                             shrinkA=0, shrinkB=0))
axb.text(cx2 + 3.7, cy2 + 0.25, "$\\frac{1}{2}$ px", fontsize=8.5,
         color=INK2, va="center", ha="left")
axb.text(cx2, 0.15, "nearest pixel is $\\frac{1}{2}$ px off the axis\n"
         "= 0.27 voxels = 270 x eps", ha="center", va="top", fontsize=8.3,
         color=INK2)

axb.text(17.0, 5.9,
         "ODD:  a pixel sits EXACTLY on the axis, so the whole centre\n"
         "row has ray.z = 0. That is 1179 columns x 179 views = 211 000\n"
         "rays on which a second coincidence (ray.x = 0) can land.\n"
         "Both hanging rays were on row 619 = (1239 $-$ 1) / 2.",
         ha="left", va="top", fontsize=8.7, color=CRIT)
axb.text(17.0, 2.5,
         "EVEN:  no pixel on the axis. The nearest is half a pixel\n"
         "off -- 0.21 mm, 0.27 voxels, 270 x above eps = 0.001 --\n"
         "so ray.z is never snapped to 0 and nothing collapses.\n"
         "Adding ONE pixel to either dimension is enough to fix it.",
         ha="left", va="top", fontsize=8.7, color=GOOD)

fig.text(0.055, 0.955,
         "Siddon forward projection never returns: a missed ray falls "
         "through to a loop bound of (unsigned long) NaN = 2$^{63}$",
         fontsize=13, color=INK, weight="bold")
fig.text(0.055, 0.918,
         "All four at once:  two ray components under eps (the ray lies along "
         "a voxel axis)   +   source outside the volume along a zeroed axis   "
         "+   the ray misses   +   no early return.      2 rays in 261 million.",
         fontsize=9.0, color=INK2)

import os
out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "siddon_missed_ray_hang.png")
fig.savefig(out, dpi=170, facecolor=SURF)
print("wrote", out)
