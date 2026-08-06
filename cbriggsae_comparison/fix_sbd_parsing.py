"""Fix SBD parsing: rename diameter->z, cell_id->birth_frame, and update 3D side-note."""
import nbformat as nbf

NB = "cbriggsae_comparison/analysis.ipynb"
with open(NB) as f:
    nb = nbf.read(f, as_version=4)

# ---- Fix 3D side-note markdown (cell with "Only 2 of the 5 datasets...") ----
for i, cell in enumerate(nb.cells):
    src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]

    if "Only 2 of the 5 datasets carry a z coordinate" in src:
        print(f"Fixing 3D side-note markdown at cell {i}")

        new_src = src.replace(
            "- **Only 2 of the 5 datasets carry a z coordinate.** The SBD/SIMI·BioCell files store\n"
            "  `frame, x, y, diameter` only, so 6 of the 7 pairwise comparisons *cannot* be 3D.",
            "- The SBD/SIMI·BioCell observation lines are `frame, x, y, z` — so **all 5 datasets carry\n"
            "  3D coordinates**. However, z is on an uncalibrated unit in both SBD files and in Cb-CSV,\n"
            "  making a naive 3D comparison invalid without per-dataset calibration."
        )
        new_src = new_src.replace(
            "The one pair where 3D is possible is the gold-standard **Cb-CSV vs Ce-tracks**. The catch\n"
            "is z calibration:",
            "The catch with 3D is z calibration across datasets:"
        )
        new_src = new_src.replace(
            "| Dataset | z / x extent ratio | z status |\n"
            "|---|---|---|\n"
            "| Ce-tracks | ~0.49 | isotropic voxels; **0.1625 µm/px** (collaborator-provided; gives ~50×25 µm ✓) |\n"
            "| Cb-CSV | ~0.04 | z on a different unit than xy — needs calibration |",
            "| Dataset | z / x extent ratio | z status |\n"
            "|---|---|---|\n"
            "| Ce-tracks | ~0.49 | isotropic voxels; **0.1625 µm/px** (confirmed; ~50×25 µm ✓) |\n"
            "| Cb-CSV | ~0.04 | z present but on a different unit than xy — needs calibration |\n"
            "| Cb-SBD / Ce-SBD | ~0.07 | z present but on a different unit than xy — needs calibration |"
        )
        new_src = new_src.replace(
            "A **naive** raw-3D distance mixes Cb's nearly-flat z (a compressed/different unit) with\n"
            "Ce's full-scale z, which artificially *lowers* the correlation. We instead calibrate\n"
            "Cb from anatomy: *C. briggsae* and *C. elegans* embryos are nearly the same size, so we\n"
            "anchor Cb's length (x) to ~50 µm — square pixels fix y at the same scale — and set the\n"
            "z-diameter equal to the y-diameter (roughly circular cross-section). Only the **relative**\n"
            "z-vs-xy weighting affects the Spearman correlation, and that is exactly what this pins down.",
            "A **naive** raw-3D distance mixes uncalibrated z scales, which artificially *lowers* the\n"
            "correlation. We calibrate Cb from anatomy: both species' embryos are ~50 µm long, so we\n"
            "anchor Cb's length (x) to 50 µm, square pixels fix y at the same scale, and set the\n"
            "z-diameter equal to the y-diameter (roughly circular cross-section). Only the **relative**\n"
            "z-vs-xy weighting affects the Spearman correlation. The same calibration logic applies to\n"
            "the SBD datasets, but since we already validated 2D ≈ 3D on the checkable gold-standard\n"
            "pair, the main workflow stays 2D."
        )
        cell["source"] = new_src

    # Fix takeaway cell
    if "Takeaway." in src and "This confirms the 2D result" in src:
        print(f"Fixing 3D takeaway at cell {i}")
        new_src = src.replace(
            "This confirms the 2D result is not hiding a z-axis\n"
            "discrepancy, and it is why the main workflow stays 2D: it needs no per-dataset z\n"
            "calibration, applies to all five datasets, and gives the same answer where 3D is checkable.",
            "This confirms the 2D result is not hiding a z-axis discrepancy. The main workflow\n"
            "stays 2D because it avoids per-dataset z calibration, applies uniformly to all five\n"
            "datasets, and gives the same answer where 3D is checkable."
        )
        cell["source"] = new_src


# ---- Fix inline parse_sbd in all code cells ----
fixes = 0
for i, cell in enumerate(nb.cells):
    if cell["cell_type"] != "code":
        continue
    src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]

    changed = False
    if "diams" in src and "float(op[3])" in src:
        src = src.replace("diams.append(float(op[3]))", "zs.append(float(op[3]))")
        changed = True
    if "diams = []" in src:
        src = src.replace("diams = []", "zs = []")
        changed = True
    if "frames, xs, ys, diams" in src:
        src = src.replace("frames, xs, ys, diams", "frames, xs, ys, zs")
        changed = True
    if "'diameter': np.array(diams)" in src:
        src = src.replace("'diameter': np.array(diams)", "'z': np.array(zs)")
        changed = True
    if "diameter" in src and "c0['diameter']" in src:
        src = src.replace("c0['diameter']", "c0['z']")
        changed = True
    if "'diameter': c0['diameter']" in src or '"diameter"' in src:
        src = src.replace("'diameter': c0['diameter']", "'z': c0['z']")
        src = src.replace('"diameter"', '"z"')
        changed = True
    # cell_id -> birth_frame
    if "cell_id = int(lines[2].split()[0])" in src:
        src = src.replace("cell_id = int(lines[2].split()[0])",
                          "birth_frame = int(lines[2].split()[0])")
        changed = True
    if "'cell_id': cell_id" in src:
        src = src.replace("'cell_id': cell_id", "'birth_frame': birth_frame")
        changed = True
    if "start_time={c0['start_time']}" in src and "c0['cell_id']" in src:
        src = src.replace("c0['cell_id']", "c0['birth_frame']")
        changed = True

    if changed:
        cell["source"] = src
        fixes += 1
        print(f"Fixed code cell {i}")

print(f"\nFixed {fixes} code cells")

nbf.write(nb, NB)
print("Notebook saved.")
