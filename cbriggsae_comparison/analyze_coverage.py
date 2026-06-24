"""
Analyze temporal coverage of cells across both datasets to find
comparable timepoints for pairwise correlation analysis.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from parse_data import load_all_data, sbd_cells_to_dataframe

OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)


def analyze_coverage():
    csv_df, sbd2_cells, sbd3_cells = load_all_data()

    active2 = [c for c in sbd2_cells if c["active"]]
    active3 = [c for c in sbd3_cells if c["active"]]

    # ============================================================
    # 1. CSV: cell count per timepoint
    # ============================================================
    csv_coverage = csv_df.groupby("time")["cell"].nunique()

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    ax = axes[0, 0]
    ax.plot(csv_coverage.index, csv_coverage.values, 'b-', alpha=0.7)
    ax.set_xlabel("Time")
    ax.set_ylabel("Number of cells")
    ax.set_title("CSV: Cell count per timepoint")
    ax.axvline(x=200, color='r', linestyle='--', alpha=0.5, label='t=200 (end)')
    ax.legend()

    # ============================================================
    # 2. SBD files: cell count per frame
    # ============================================================
    for idx, (cells, label) in enumerate([(active2, "SBD2"), (active3, "SBD3")]):
        # Build frame -> n_cells mapping
        frame_counts = {}
        for c in cells:
            for f in c["frame"]:
                frame_counts[f] = frame_counts.get(f, 0) + 1

        frames = sorted(frame_counts.keys())
        counts = [frame_counts[f] for f in frames]

        ax = axes[0, 1 + idx]
        ax.plot(frames, counts, 'g-', alpha=0.7, linewidth=0.5)
        ax.set_xlabel("Frame")
        ax.set_ylabel("Number of cells")
        ax.set_title(f"{label}: Cell count per frame")

        # Mark the "knee" where coverage drops rapidly
        # Find the frame where coverage is max
        max_frame = frames[np.argmax(counts)]
        max_count = max(counts)
        ax.axvline(x=max_frame, color='orange', linestyle='--', alpha=0.5,
                   label=f'peak at frame {max_frame}')
        ax.legend()

    # ============================================================
    # 3. SBD: individual cell track spans
    # ============================================================
    for idx, (cells, label) in enumerate([(active2, "SBD2"), (active3, "SBD3")]):
        ax = axes[1, idx]
        # Sort cells by first frame
        sorted_cells = sorted(cells, key=lambda c: c["frame"].min())
        for i, c in enumerate(sorted_cells):
            if i % 10 == 0:  # plot every 10th cell for readability
                ax.plot([c["frame"].min(), c["frame"].max()], [i, i],
                        'b-', alpha=0.3, linewidth=0.5)
        ax.set_xlabel("Frame")
        ax.set_ylabel("Cell index (sorted by start)")
        ax.set_title(f"{label}: Cell track spans")

    # ============================================================
    # 4. Coverage at late frames (zoom in)
    # ============================================================
    ax = axes[1, 2]
    for cells, label, color in [(active2, "SBD2", "blue"), (active3, "SBD3", "red")]:
        frame_counts = {}
        for c in cells:
            for f in c["frame"]:
                frame_counts[f] = frame_counts.get(f, 0) + 1

        frames = sorted(frame_counts.keys())
        counts = [frame_counts[f] for f in frames]

        # Only show 2nd half of frames
        max_f = max(frames)
        mid_f = max_f // 2
        late_frames = [f for f in frames if f >= mid_f]
        late_counts = [frame_counts[f] for f in late_frames]

        ax.plot(late_frames, late_counts, '-', alpha=0.7, linewidth=0.8,
                color=color, label=label)

    ax.set_xlabel("Frame (2nd half)")
    ax.set_ylabel("Number of cells")
    ax.set_title("SBD: Cell count - late frames only")
    ax.legend()

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "coverage_analysis.png", dpi=150)
    plt.close()
    print("Saved coverage_analysis.png")

    # ============================================================
    # 5. Key metrics: where to stop
    # ============================================================
    print("\n=== Coverage Summary ===")
    print(f"CSV: {csv_coverage.idxmax()} cells at peak t={csv_coverage.max()}, "
          f"cells at t=200: {csv_df[csv_df.time==200].cell.nunique()}")

    for cells, label in [(active2, "SBD2"), (active3, "SBD3")]:
        frame_counts = {}
        for c in cells:
            for f in c["frame"]:
                frame_counts[f] = frame_counts.get(f, 0) + 1

        frames = sorted(frame_counts.keys())
        max_count = max(frame_counts.values())
        peak_frame = [f for f, c in frame_counts.items() if c == max_count][0]

        # Find the last frame where coverage is >= 50% of max
        half_max = max_count / 2
        last_good_frame = max(f for f, c in frame_counts.items() if c >= half_max)

        # Find the last frame where coverage is >= 80% of max
        p80_max = max_count * 0.8
        last_good_80 = max(f for f, c in frame_counts.items() if c >= p80_max)

        # Count cells that have observations up to last_good_frame
        cells_at_end = sum(1 for c in cells if c["frame"].max() >= last_good_frame)

        print(f"\n{label}:")
        print(f"  Max coverage: {max_count} cells at frame {peak_frame}")
        print(f"  Last frame with >=50% max: frame {last_good_frame} "
              f"({frame_counts[last_good_frame]} cells)")
        print(f"  Last frame with >=80% max: frame {last_good_80} "
              f"({frame_counts[last_good_80]} cells)")
        print(f"  Cells tracked to frame {last_good_frame}: {cells_at_end}")

    return csv_df, active2, active3


if __name__ == "__main__":
    analyze_coverage()
