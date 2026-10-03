"""Compatibility shim: shared publication style now lives in ``publication/style.py``.

``save_panel_crops`` remains here for the historical exploratory ``main.py``.
"""

from pathlib import Path

from publication.style import (  # noqa: F401
    COLORS, EDGE_RETENTION_CMAP, NULL_MODEL_COLORS, SEMANTIC_COLORS, SPECIES_COLORS,
    color_ramp, configure, save_figure,
)


def save_panel_crops(
    fig, axes, output_root, config_dirs, filename, extra_axes=None
) -> None:
    """Save each axis of a horizontal comparison figure as its own PNG."""
    from matplotlib.transforms import Bbox

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    extras = extra_axes if extra_axes is not None else [None] * len(axes)
    for ax, extra_ax, dirname in zip(axes, extras, config_dirs):
        target = Path(output_root) / dirname / filename
        target.parent.mkdir(parents=True, exist_ok=True)
        visible_axes = {ax}
        if extra_ax is not None:
            visible_axes.add(extra_ax)
        original_visibility = {item: item.get_visible() for item in fig.axes}
        for item in fig.axes:
            item.set_visible(item in visible_axes)
        boxes = [ax.get_tightbbox(renderer)]
        if extra_ax is not None:
            boxes.append(extra_ax.get_tightbbox(renderer))
        bbox = Bbox.union(boxes).transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(
            target, dpi=300, bbox_inches=bbox.expanded(1.12, 1.12),
            pad_inches=0.10, facecolor="white",
        )
        for item, visible in original_visibility.items():
            item.set_visible(visible)
