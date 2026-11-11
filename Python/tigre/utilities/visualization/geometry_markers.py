"""Coordinate markers shared by the detector builder and geometry animation."""
import numpy as np


def draw_world_axes(axes, length):
    """Draw labelled +x/+y/+z arrows at world zero, at 30 percent opacity."""
    arrows = []
    for axis, (label, color) in enumerate(zip("xyz", ("red", "green", "blue"))):
        direction = np.eye(3)[axis] * length
        arrows.append(axes.quiver(0, 0, 0, *direction, color=color, linewidth=2,
            alpha=0.3, arrow_length_ratio=0.25, label=f"World +{label}"))
        axes.text(*(direction * 1.12), label, color=color, fontweight="bold", alpha=0.3)
    return arrows
