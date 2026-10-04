"""Terminal primary family: Figures 1 amendment-6, S1-S3 and Table 1 from validated caches."""

from __future__ import annotations

from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from publication import style, tables
from publication.captions.terminal import NULL_SCHEMATIC
from publication.figures import (terminal_canonical, terminal_cell_types, terminal_global, terminal_s3,
                                 terminal_subtree_map, terminal_within_type)


def render_panels(data, out: Path) -> None:
    """Every panel asset the terminal wrappers include, written into ``out``."""
    out = Path(out)
    style.configure()
    terminal_canonical.render_figure1_amendment(out)
    plt.close(terminal_global.plot_main(data.global_twr, data.global_nulls, display=data.global_display,
                                        out_dir=out, panel_letter=None))
    plt.close(terminal_global.plot_support_b(data.global_twr, display=data.global_display, out_dir=out))
    plt.close(terminal_global.plot_support_c(data.global_twr, display=data.global_display, out_dir=out))
    shutil.copy2(NULL_SCHEMATIC, out / NULL_SCHEMATIC.name)
    terminal_subtree_map.emit_figure(data.subtree_nodes, data.subtree_edges, data.min_cells,
                                     out / "fig4_ce_subtree_map_panel")
    terminal_canonical.plot_primary_rows(data.canonical_display, out)
    terminal_canonical.plot_summary(data.canonical_display, metric="r",
                                    component_names=terminal_canonical.SUPPLEMENT_COMPONENTS, out_dir=out)
    terminal_cell_types.plot_retention_heatmap(data.cell_type_retention, out_dir=out)
    terminal_cell_types.plot_endpoint_ledger(data.cell_type_keypoints, out_dir=out)
    within = data.within_type
    terminal_within_type.plot_within_type(within["fronts"], within["summary"], within["aggregate"],
                                          within["endpoints"], out_dir=out, display_mode="endpoint")
    terminal_s3.render_publication(data.s3_projections, data.s3_references, out_dir=out)
    tables.write_subtree_statistics(out, data.canonical, data.subtree_summary, data.min_cells)
