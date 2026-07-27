"""Marimo notebook: Terminal-only Pareto front analysis."""

import marimo

__generated_with = "0.23.5"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Pareto Front Analysis of Terminal Cell Lineage

    **Comparative Summary: C. elegans vs C. briggsae**

    For terminal cells, does the real lineage sit on the Pareto frontier
    trading off physical proximity vs. expression similarity?  All Pareto
    fronts use **std-based scaling**: costs as z-scores relative to the
    cousin-randomisation null, so (0,0) = null mean and 1 unit = 1σ.
    """)
    return


@app.cell
def _():
    import json
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.optimize import linear_sum_assignment
    from scipy.spatial.distance import cosine
    from matplotlib.lines import Line2D
    from collections import defaultdict
    from joblib import Parallel, delayed
    import matplotlib as mpl
    mpl.rcParams['figure.dpi'] = 150
    sns.set_style("whitegrid")
    return Line2D, Parallel, delayed, np, pd, plt


@app.cell
def _():
    import sys, os
    _repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if _repo_root not in sys.path:
        sys.path.insert(0, _repo_root)
    from terminal_pareto import data_loader as dl
    from terminal_pareto import pareto_engine as pe
    from terminal_pareto import lineage_metrics as lm

    return dl, lm, pe


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Data Loading
    """)
    return


@app.cell
def _(dl):
    lineage_data = dl.load_json(dl._REPO_ROOT + '/data/cell_lineage.json')
    apoptotic_cells = dl.load_apoptotic()
    xyz_ce, valid_ce = dl.load_elegans_tracking(dl.T_CE)
    xyz_ce_early, valid_ce_early = dl.load_elegans_tracking(dl.T_CE_EARLY)
    xyz_cb, valid_cb = dl.load_briggsae_tracking(dl.T_CB)
    xyz_cb_new, valid_cb_new = dl.load_briggsae_tracking(dl.T_CB_NEW, dl.CB_NEW_PATH)
    protein_exp = dl.load_protein_expression()
    prot_sel = dl.load_prot_sel()
    ce_rna = dl.load_ce_rna()
    cb_rna = dl.load_cb_rna()
    rna_sel = dl.load_rna_sel()
    v_prot = [n for n in valid_ce if n in protein_exp.index]
    v_ce_rna = [n for n in valid_ce if n in ce_rna.index]
    v_prot_early = [n for n in valid_ce_early if n in protein_exp.index]
    v_ce_rna_early = [n for n in valid_ce_early if n in ce_rna.index]
    v_cb_rna = [n for n in valid_cb if n in cb_rna.index]
    v_cb_rna_new = [n for n in valid_cb_new if n in cb_rna.index]
    print(f"CE: T=255→{len(valid_ce)}, T=230→{len(valid_ce_early)}")
    print(f"CB old: {len(valid_cb)}, CB new: {len(valid_cb_new)}")
    print(f"Valid: prot={len(v_prot)}, ceRNA={len(v_ce_rna)}, cbRNA(new)={len(v_cb_rna_new)}")
    return (
        cb_rna,
        ce_rna,
        lineage_data,
        prot_sel,
        protein_exp,
        rna_sel,
        v_cb_rna,
        v_cb_rna_new,
        v_ce_rna,
        v_ce_rna_early,
        v_prot,
        v_prot_early,
        xyz_cb,
        xyz_cb_new,
        xyz_ce,
        xyz_ce_early,
    )


@app.cell
def _(lineage_data, pe):
    gp_map = pe.build_grandparent_map(lineage_data)
    print(f"Grandparent map: {len(gp_map)} cells")
    return (gp_map,)


@app.cell
def _(lineage_data, lm):
    tree_index = lm.build_lineage_tree_index(lineage_data)
    max_depth = max(v['depth'] for v in tree_index.values())
    print(f"Lineage tree index: {len(tree_index)} nodes, max depth={max_depth}")
    return (tree_index,)



@app.cell
def _(
    cb_rna,
    ce_rna,
    dl,
    lineage_data,
    pe,
    prot_sel,
    protein_exp,
    rna_sel,
    v_cb_rna,
    v_cb_rna_new,
    v_ce_rna,
    v_ce_rna_early,
    v_prot,
    v_prot_early,
    xyz_cb,
    xyz_cb_new,
    xyz_ce,
    xyz_ce_early,
):
    tn_prot, tp_prot = dl.collect_terminals(lineage_data, v_prot)
    tn_er, tp_er = dl.collect_terminals(lineage_data, v_ce_rna)
    tn_prot_early, tp_prot_early = dl.collect_terminals(lineage_data, v_prot_early)
    tn_er_early, tp_er_early = dl.collect_terminals(lineage_data, v_ce_rna_early)
    tn_cb, tp_cb = dl.collect_terminals(lineage_data, v_cb_rna)
    tn_cb_new, tp_cb_new = dl.collect_terminals(lineage_data, v_cb_rna_new)
    print(f"Terminals: prot={len(tn_prot)}, RNA={len(tn_er)}, CB_new={len(tn_cb_new)}")
    xm_prot, em_prot, fd_prot = pe.build_cost_matrices(tn_prot, tp_prot, xyz_ce, protein_exp, prot_sel)
    xm_er, em_er, fd_er = pe.build_cost_matrices(tn_er, tp_er, xyz_ce, ce_rna, rna_sel)
    xm_prot_early, em_prot_early, _ = pe.build_cost_matrices(tn_prot_early, tp_prot_early, xyz_ce_early, protein_exp, prot_sel)
    xm_er_early, em_er_early, _ = pe.build_cost_matrices(tn_er_early, tp_er_early, xyz_ce_early, ce_rna, rna_sel)
    xm_cb, em_cb, fd_cb = pe.build_cost_matrices(tn_cb, tp_cb, xyz_cb, cb_rna, rna_sel)
    xm_cb_new, em_cb_new, fd_cb_new = pe.build_cost_matrices(tn_cb_new, tp_cb_new, xyz_cb_new, cb_rna, rna_sel)
    print("Cost matrices built.")
    return (
        em_cb,
        em_cb_new,
        em_er,
        em_er_early,
        em_prot,
        em_prot_early,
        fd_cb_new,
        fd_er,
        fd_prot,
        tn_cb,
        tn_cb_new,
        tn_er,
        tn_er_early,
        tn_prot,
        tn_prot_early,
        tp_cb_new,
        tp_er,
        tp_prot,
        xm_cb,
        xm_cb_new,
        xm_er,
        xm_er_early,
        xm_prot,
        xm_prot_early,
    )


@app.cell
def _(
    em_cb,
    em_cb_new,
    em_er,
    em_er_early,
    em_prot,
    em_prot_early,
    gp_map,
    pe,
    tn_cb,
    tn_cb_new,
    tn_er,
    tn_er_early,
    tn_prot,
    tn_prot_early,
    xm_cb,
    xm_cb_new,
    xm_er,
    xm_er_early,
    xm_prot,
    xm_prot_early,
):
    prot_groups = pe.build_cousin_groups(tn_prot, gp_map)
    prot_rs = pe.compute_cousin_random_stats(xm_prot, em_prot, prot_groups, n_random=500, seed=42)
    prot_early_groups = pe.build_cousin_groups(tn_prot_early, gp_map)
    prot_early_rs = pe.compute_cousin_random_stats(xm_prot_early, em_prot_early, prot_early_groups, n_random=500, seed=42)
    er_groups = pe.build_cousin_groups(tn_er, gp_map)
    er_rs = pe.compute_cousin_random_stats(xm_er, em_er, er_groups, n_random=500, seed=42)
    er_early_groups = pe.build_cousin_groups(tn_er_early, gp_map)
    er_early_rs = pe.compute_cousin_random_stats(xm_er_early, em_er_early, er_early_groups, n_random=500, seed=42)
    cb_groups = pe.build_cousin_groups(tn_cb, gp_map)
    cb_rs = pe.compute_cousin_random_stats(xm_cb, em_cb, cb_groups, n_random=500, seed=42)
    cb_new_groups = pe.build_cousin_groups(tn_cb_new, gp_map)
    cb_new_rs = pe.compute_cousin_random_stats(xm_cb_new, em_cb_new, cb_new_groups, n_random=500, seed=42)
    print(f"xyz_std: Prot={prot_rs['xyz_std']:.1f}, RNA={er_rs['xyz_std']:.1f}, CB_new={cb_new_rs['xyz_std']:.1f}")
    markers = pe.MARKERS
    return cb_new_rs, er_rs, markers, prot_rs


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ---
    ## 1. Tracking Replicate Consistency
    """)
    return


@app.cell
def _(dl, gp_map, lineage_data, pe):
    ce_replicates = dl.CE_REPLICATES
    cb_replicates = dl.CB_REPLICATES_OLD + dl.CB_REPLICATES_NEW

    def _run_rep(xyz_map, valid_names, exp_df, sel_features):
        return pe.run_replicate_pareto(xyz_map, valid_names, exp_df, sel_features, lineage_data, gp_map, n_random=300, seed=42)

    RUN_REP = _run_rep
    print(f"Replicates: {len(dl.CB_REPLICATES_OLD)} old + {len(dl.CB_REPLICATES_NEW)} new briggsae")
    return RUN_REP, cb_replicates, ce_replicates


@app.cell
def _(RUN_REP, ce_replicates, dl, prot_sel, protein_exp):
    _reps = {}
    for _lbl, _path, _tcut in ce_replicates:
        _xyz, _val = dl.load_elegans_tracking(_tcut, path=_path)
        _r = RUN_REP(_xyz, _val, protein_exp, prot_sel)
        _reps[_lbl] = dict(res=_r, tcut=_tcut)
        print(f"  {_lbl} (T={_tcut}): n={_r['n']}")
    ce_protein_reps = _reps
    return (ce_protein_reps,)


@app.cell
def _(RUN_REP, ce_replicates, ce_rna, dl, rna_sel):
    _reps = {}
    for _lbl, _path, _tcut in ce_replicates:
        _xyz, _val = dl.load_elegans_tracking(_tcut, path=_path)
        _r = RUN_REP(_xyz, _val, ce_rna, rna_sel)
        _reps[_lbl] = dict(res=_r, tcut=_tcut)
        print(f"  {_lbl} (T={_tcut}): n={_r['n']}")
    ce_rna_reps = _reps
    return (ce_rna_reps,)


@app.cell
def _(RUN_REP, cb_replicates, cb_rna, dl, rna_sel):
    _reps = {}
    for _lbl, _path, _tcut in cb_replicates:
        _xyz, _val = dl.load_briggsae_tracking(_tcut, path=_path)
        _r = RUN_REP(_xyz, _val, cb_rna, rna_sel)
        _reps[_lbl] = dict(res=_r, tcut=_tcut)
        print(f"  {_lbl} (T={_tcut}): n={_r['n']}")
    cb_rna_reps = _reps
    return (cb_rna_reps,)


@app.cell
def _(Line2D, ce_protein_reps, markers, pe, plt):
    def _plot():
        fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(16, 7))
        first = True
        for label, d in ce_protein_reps.items():
            r = d['res']; color = '#1976D2' if 'embryo1' in label else ('#4CAF50' if 'embryo2' in label else '#FF9800')
            xa, ea, edge, kp = r['xyz_arr'], r['exp_arr'], r['edge_arr'], r['kp']
            lx, le = pe.lineage_std_position(r['xm'], r['em'], r['rs'])
            if first:
                nx, ny = pe.get_null_cloud(r['rs'], max_points=250)
                ax_l.scatter(nx, ny, c=color, s=8, alpha=0.10, zorder=0); first = False
            ax_l.plot(xa, ea, color=color, lw=1.5, alpha=0.85, label=f'{label} (T={d["tcut"]}, n={r["n"]})')
            ax_l.scatter(lx, le, color=color, marker='X', s=70, edgecolors='white', lw=0.8, zorder=6)
            ax_r.plot(xa, edge, color=color, lw=2, alpha=0.85, label=f'{label} (T={d["tcut"]}, n={r["n"]})')
            for pk, mk in markers.items():
                ax_r.scatter(xa[kp[pk]], edge[kp[pk]], color=color, marker=mk, s=50, zorder=8, edgecolors='white', lw=0.5)
        ax_l.axhline(0, color='gray', lw=0.8, ls=':', alpha=0.4); ax_l.axvline(0, color='gray', lw=0.8, ls=':', alpha=0.4)
        ax_l.set_xlabel('Spatial Cost (σ from null)'); ax_l.set_ylabel('Expression Cost (σ from null)')
        ax_l.set_title('C. elegans Protein — Pareto Fronts by Replicate', fontsize=12)
        ax_l.legend(fontsize=8, loc='upper right'); ax_l.grid(True, alpha=0.3)
        ax_r.set_xlabel('Spatial Cost (σ from null)'); ax_r.set_ylabel('Edge Retention Ratio')
        ax_r.set_title('C. elegans Protein — Edge Retention by Replicate', fontsize=12)
        ax_r.set_ylim(0, 1.05); ax_r.grid(True, alpha=0.3)
        ax_r.legend(handles=list(ax_r.get_lines()) + [
            Line2D([0],[0], color='black', marker='s', ls='', markersize=6, label='Expr opt'),
            Line2D([0],[0], color='black', marker='D', ls='', markersize=6, label='Max ER'),
            Line2D([0],[0], color='black', marker='o', ls='', markersize=6, label='Spatial opt'),
        ], fontsize=7.5, loc='best')
        fig.suptitle('Tracking Replicate Consistency — C. elegans Protein', fontsize=14, fontweight='bold')
        fig.tight_layout()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'ce_protein_replicates.png'), dpi=150, bbox_inches='tight')
        return fig
    _fig = _plot()
    return _fig


@app.cell
def _(Line2D, ce_rna_reps, markers, pe, plt):
    def _plot():
        fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(16, 7))
        first = True
        for label, d in ce_rna_reps.items():
            r = d['res']; color = '#1976D2' if 'embryo1' in label else ('#4CAF50' if 'embryo2' in label else '#FF9800')
            xa, ea, edge, kp = r['xyz_arr'], r['exp_arr'], r['edge_arr'], r['kp']
            lx, le = pe.lineage_std_position(r['xm'], r['em'], r['rs'])
            if first:
                nx, ny = pe.get_null_cloud(r['rs'], max_points=250)
                ax_l.scatter(nx, ny, c=color, s=8, alpha=0.10, zorder=0); first = False
            ax_l.plot(xa, ea, color=color, lw=1.5, alpha=0.85, label=f'{label} (T={d["tcut"]}, n={r["n"]})')
            ax_l.scatter(lx, le, color=color, marker='X', s=70, edgecolors='white', lw=0.8, zorder=6)
            ax_r.plot(xa, edge, color=color, lw=2, alpha=0.85, label=f'{label} (T={d["tcut"]}, n={r["n"]})')
            for pk, mk in markers.items():
                ax_r.scatter(xa[kp[pk]], edge[kp[pk]], color=color, marker=mk, s=50, zorder=8, edgecolors='white', lw=0.5)
        ax_l.axhline(0, color='gray', lw=0.8, ls=':', alpha=0.4); ax_l.axvline(0, color='gray', lw=0.8, ls=':', alpha=0.4)
        ax_l.set_xlabel('Spatial Cost (σ from null)'); ax_l.set_ylabel('Expression Cost (σ from null)')
        ax_l.set_title('C. elegans RNA — Pareto Fronts by Replicate', fontsize=12)
        ax_l.legend(fontsize=8, loc='upper right'); ax_l.grid(True, alpha=0.3)
        ax_r.set_xlabel('Spatial Cost (σ from null)'); ax_r.set_ylabel('Edge Retention Ratio')
        ax_r.set_title('C. elegans RNA — Edge Retention by Replicate', fontsize=12)
        ax_r.set_ylim(0, 1.05); ax_r.grid(True, alpha=0.3)
        ax_r.legend(handles=list(ax_r.get_lines()) + [
            Line2D([0],[0], color='black', marker='s', ls='', markersize=6, label='Expr opt'),
            Line2D([0],[0], color='black', marker='D', ls='', markersize=6, label='Max ER'),
            Line2D([0],[0], color='black', marker='o', ls='', markersize=6, label='Spatial opt'),
        ], fontsize=7.5, loc='best')
        fig.suptitle('Tracking Replicate Consistency — C. elegans RNA', fontsize=14, fontweight='bold')
        fig.tight_layout()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'ce_rna_replicates.png'), dpi=150, bbox_inches='tight')
        return fig
    _fig = _plot()
    return _fig


@app.cell
def _(Line2D, cb_rna_reps, markers, pe, plt):
    def _plot():
        fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(16, 7))
        first = True
        for label, d in cb_rna_reps.items():
            r = d['res']; tcut = d['tcut']
            is_new = 'she-1' in label or 'AF16' in label
            color = '#1B5E20' if is_new else '#1976D2'
            xa, ea, edge, kp = r['xyz_arr'], r['exp_arr'], r['edge_arr'], r['kp']
            lx, le = pe.lineage_std_position(r['xm'], r['em'], r['rs'])
            if first:
                nx, ny = pe.get_null_cloud(r['rs'], max_points=250)
                ax_l.scatter(nx, ny, c=color, s=8, alpha=0.10, zorder=0); first = False
            ax_l.plot(xa, ea, color=color, lw=1.5, alpha=0.85, label=f'{label} (T={tcut}, n={r["n"]})')
            ax_l.scatter(lx, le, color=color, marker='X', s=70, edgecolors='white', lw=0.8, zorder=6)
            ax_r.plot(xa, edge, color=color, lw=2, alpha=0.85, label=f'{label} (T={tcut}, n={r["n"]})')
            for pk, mk in markers.items():
                ax_r.scatter(xa[kp[pk]], edge[kp[pk]], color=color, marker=mk, s=50, zorder=8, edgecolors='white', lw=0.5)
        ax_l.axhline(0, color='gray', lw=0.8, ls=':', alpha=0.4); ax_l.axvline(0, color='gray', lw=0.8, ls=':', alpha=0.4)
        ax_l.set_xlabel('Spatial Cost (σ from null)'); ax_l.set_ylabel('Expression Cost (σ from null)')
        ax_l.set_title('C. briggsae RNA — Pareto Fronts by Replicate', fontsize=12)
        ax_l.legend(fontsize=6, loc='upper right'); ax_l.grid(True, alpha=0.3)
        ax_r.set_xlabel('Spatial Cost (σ from null)'); ax_r.set_ylabel('Edge Retention Ratio')
        ax_r.set_title('C. briggsae RNA — Edge Retention by Replicate', fontsize=12)
        ax_r.set_ylim(0, 1.05); ax_r.grid(True, alpha=0.3)
        ax_r.legend(handles=list(ax_r.get_lines()) + [
            Line2D([0],[0], color='black', marker='s', ls='', markersize=6, label='Expr opt'),
            Line2D([0],[0], color='black', marker='D', ls='', markersize=6, label='Max ER'),
            Line2D([0],[0], color='black', marker='o', ls='', markersize=6, label='Spatial opt'),
        ], fontsize=7.5, loc='best')
        fig.suptitle('Tracking Replicate Consistency — C. briggsae RNA', fontsize=14, fontweight='bold')
        fig.tight_layout()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'cb_rna_replicates.png'), dpi=150, bbox_inches='tight')
        return fig
    _fig = _plot()
    return _fig


@app.cell
def _(Line2D, cb_rna_reps, markers, pe, plt):
    def _plot():
        _wt_colors = {'210519ZZY0874p1 (AF16)': '#E65100', '210519ZZY0874p5 (AF16)': '#FF9800'}
        fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(16, 7))
        for label, d in cb_rna_reps.items():
            r = d['res']; is_wt = 'AF16' in label
            color = _wt_colors.get(label, '#43A047')
            lw = 2.5 if is_wt else 1.5; al = 0.9 if is_wt else 0.7
            xa, ea, edge, kp = r['xyz_arr'], r['exp_arr'], r['edge_arr'], r['kp']
            lx, le = pe.lineage_std_position(r['xm'], r['em'], r['rs'])
            nx, ny = pe.get_null_cloud(r['rs'], max_points=200)
            ax_l.scatter(nx, ny, c=color, s=8, alpha=0.08, zorder=0)
            ax_l.plot(xa, ea, color=color, lw=lw, alpha=al, label=f'{label} (n={r["n"]})')
            ax_l.scatter(lx, le, color=color, marker='X', s=70, edgecolors='white', lw=0.8, zorder=6)
            ax_r.plot(xa, edge, color=color, lw=lw, alpha=al)
            for pk, mk in markers.items():
                ax_r.scatter(xa[kp[pk]], edge[kp[pk]], color=color, marker=mk, s=40, zorder=8, edgecolors='white', lw=0.5)
        ax_l.axhline(0, color='gray', lw=0.8, ls=':', alpha=0.4); ax_l.axvline(0, color='gray', lw=0.8, ls=':', alpha=0.4)
        ax_l.set_xlabel('Spatial Cost (σ from null)'); ax_l.set_ylabel('Expression Cost (σ from null)')
        ax_l.set_title('C. briggsae RNA — WT vs she-1 Pareto Fronts', fontsize=12)
        ax_l.legend(fontsize=7.5, loc='upper right'); ax_l.grid(True, alpha=0.3)
        ax_r.set_xlabel('Spatial Cost (σ from null)'); ax_r.set_ylabel('Edge Retention Ratio')
        ax_r.set_title('C. briggsae RNA — WT vs she-1 Edge Retention', fontsize=12)
        ax_r.set_ylim(0, 1.05); ax_r.grid(True, alpha=0.3)
        ax_r.legend(handles=[
            Line2D([0],[0], color='#E65100', lw=2.5, label='WT (AF16)'),
            Line2D([0],[0], color='#43A047', lw=1.5, label='she-1 mutant'),
        ], fontsize=7.5, loc='best')
        fig.suptitle('C. briggsae — Wildtype vs she-1 Mutant', fontsize=14, fontweight='bold')
        fig.tight_layout()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'cb_wt_vs_she1.png'), dpi=150, bbox_inches='tight')
        return fig
    _fig = _plot()
    return _fig


@app.cell
def _(cb_rna_reps, ce_protein_reps, ce_rna_reps, pd):
    def _b():
        _rows = []
        for _lb, _d in ce_protein_reps.items():
            _r = _d['res']; _ea = _r['edge_arr']; _kp = _r['kp']
            _rows.append({'Config': 'CE Protein', 'Replicate': _lb, 'T': _d['tcut'], 'N': _r['n'],
                'Expr Opt ER': _ea[_kp['expr_opt_idx']], 'Max ER': _ea[_kp['max_er_idx']], 'Spatial Opt ER': _ea[_kp['spatial_opt_idx']]})
        for _lb, _d in ce_rna_reps.items():
            _r = _d['res']; _ea = _r['edge_arr']; _kp = _r['kp']
            _rows.append({'Config': 'CE RNA', 'Replicate': _lb, 'T': _d['tcut'], 'N': _r['n'],
                'Expr Opt ER': _ea[_kp['expr_opt_idx']], 'Max ER': _ea[_kp['max_er_idx']], 'Spatial Opt ER': _ea[_kp['spatial_opt_idx']]})
        for _lb, _d in cb_rna_reps.items():
            _r = _d['res']; _ea = _r['edge_arr']; _kp = _r['kp']
            _rows.append({'Config': 'CB RNA', 'Replicate': _lb, 'T': _d['tcut'], 'N': _r['n'],
                'Expr Opt ER': _ea[_kp['expr_opt_idx']], 'Max ER': _ea[_kp['max_er_idx']], 'Spatial Opt ER': _ea[_kp['spatial_opt_idx']]})
        return pd.DataFrame(_rows).set_index(['Config', 'Replicate']).round(3)
    _b()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ---
    ## 2. Full Tree — Cross-Dataset Comparison
    """)
    return


@app.cell
def _(
    em_cb_new,
    em_er,
    em_prot,
    pd,
    pe,
    tn_cb_new,
    tn_er,
    tn_prot,
    tp_cb_new,
    tp_er,
    tp_prot,
    xm_cb_new,
    xm_er,
    xm_prot,
):
    res_prot_xs = pe.pareto_sweep(xm_prot, em_prot, tp_prot)
    res_er_xs   = pe.pareto_sweep(xm_er,   em_er,   tp_er)
    res_cb_new  = pe.pareto_sweep(xm_cb_new, em_cb_new, tp_cb_new)
    full_tree_df = pd.DataFrame({
        'Dataset': ['C. elegans\nProtein (T=255)', 'C. elegans\nRNA (T=255)', 'C. briggsae\nRNA (T=143)'],
        'N Cells': [len(tn_prot), len(tn_er), len(tn_cb_new)],
        'Expr Opt ER': [res_prot_xs['expr_opt_er'], res_er_xs['expr_opt_er'], res_cb_new['expr_opt_er']],
        'Max ER': [res_prot_xs['max_er'], res_er_xs['max_er'], res_cb_new['max_er']],
        'Spatial Opt ER': [res_prot_xs['spatial_opt_er'], res_er_xs['spatial_opt_er'], res_cb_new['spatial_opt_er']],
        'Closest Dist': [res_prot_xs['closest_dist'], res_er_xs['closest_dist'], res_cb_new['closest_dist']],
    }).set_index('Dataset')
    full_tree_df
    return res_cb_new, res_er_xs, res_prot_xs


@app.cell
def _(
    Line2D,
    cb_new_rs,
    em_cb_new,
    em_er,
    em_prot,
    er_rs,
    markers,
    pe,
    plt,
    prot_rs,
    tp_cb_new,
    tp_er,
    tp_prot,
    xm_cb_new,
    xm_er,
    xm_prot,
):
    def _plot():
        datasets = [
            ('C. elegans Protein (T=255)', xm_prot, em_prot, tp_prot, '#1976D2', prot_rs),
            ('C. elegans RNA (T=255)',     xm_er,   em_er,   tp_er,   '#4CAF50', er_rs),
            ('C. briggsae RNA (T=143)',    xm_cb_new, em_cb_new, tp_cb_new, '#1B5E20', cb_new_rs),
        ]
        fig = plt.figure(figsize=(22, 14))
        gs_pareto = fig.add_gridspec(3, 1, left=0.03, right=0.38, hspace=0.35)
        for idx, (name, xm, em, tp, color, rs) in enumerate(datasets):
            ax = fig.add_subplot(gs_pareto[idx])
            xa, ea, _, kp = pe.compute_std_scaled_pareto(xm, em, tp, rs)
            lx, le = pe.lineage_std_position(xm, em, rs); n = len(tp)
            nx, ny = pe.get_null_cloud(rs, max_points=300)
            ax.scatter(nx, ny, c=color, s=8, alpha=0.18, zorder=0)
            ax.axhline(0, color='gray', lw=0.8, ls=':', alpha=0.4); ax.axvline(0, color='gray', lw=0.8, ls=':', alpha=0.4)
            ax.plot(xa, ea, color=color, lw=1.8, alpha=0.9, label='Pareto front')
            ax.scatter(lx, le, color=color, marker='X', s=80, edgecolors='white', lw=0.8, zorder=6)
            ax.scatter(lx, le, color='black', marker='X', s=80, edgecolors='white', lw=0.8, zorder=10,
                       label=f'Lineage ({lx:.0f}σ, {le:.0f}σ)')
            ax.annotate(f'Lineage\n({lx:.0f}σ, {le:.0f}σ)', (lx, le), xytext=(10, 10),
                        textcoords='offset points', fontsize=8, color='black', fontweight='bold')
            ax.set_xlabel('Spatial Cost (σ from null)'); ax.set_ylabel('Expression Cost (σ from null)')
            ax.set_title(f'{name} (n={n})', fontsize=11, color=color, fontweight='bold')
            ax.legend(fontsize=7, loc='upper right'); ax.grid(True, alpha=0.3)
        ax_er = fig.add_subplot(fig.add_gridspec(1, 1, left=0.42, right=0.98, top=0.92, bottom=0.08)[0])
        for name, xm, em, tp, color, rs in datasets:
            xa, ___, edge, kp = pe.compute_std_scaled_pareto(xm, em, tp, rs); n = len(tp)
            ax_er.plot(xa, edge, color=color, lw=2, alpha=0.85, label=f'{name} (n={n})')
            for pk, mk in markers.items():
                _i = kp[pk]
                ax_er.scatter(xa[_i], edge[_i], color=color, marker=mk, s=50, zorder=8, edgecolors='white', lw=0.5)
                ax_er.annotate(f'{edge[_i]:.3f}', (xa[_i], edge[_i]+0.04), color=color, fontsize=7, ha='center', fontweight='bold')
        ax_er.legend(handles=list(ax_er.get_lines()) + [
            Line2D([0],[0], color='black', marker='s', ls='', markersize=6, label='Expr opt'),
            Line2D([0],[0], color='black', marker='D', ls='', markersize=6, label='Max ER'),
            Line2D([0],[0], color='black', marker='o', ls='', markersize=6, label='Spatial opt'),
        ], fontsize=7.5, loc='best')
        ax_er.set_xlabel('Spatial Cost (σ from null)'); ax_er.set_ylabel('Edge Retention Ratio')
        ax_er.set_title('Edge Retention along Pareto Front', fontsize=12)
        ax_er.set_ylim(0, 1.05); ax_er.grid(True, alpha=0.3)
        fig.suptitle('Full Tree Cross-Species Comparison — Std-Based Scaling\n(centred at null mean, all at full cutoffs)',
                     fontsize=14, fontweight='bold', y=0.98)
        plt.show()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'cross_species_comparison.png'), dpi=150, bbox_inches='tight')
        return fig
    _fig = _plot()
    return _fig


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ---
    ### 2a. Cross-Species — 2D (XY-plane only)

    C. briggsae tracking has lower z-axis accuracy. Restricting to XY removes this confound.
    C. elegans (isotropic) serves as a control.
    """)
    return


@app.cell
def _(
    cb_rna,
    ce_rna,
    gp_map,
    pd,
    pe,
    prot_sel,
    protein_exp,
    res_cb_new,
    res_er_xs,
    res_prot_xs,
    rna_sel,
    tn_cb_new,
    tn_er,
    tn_prot,
    tp_cb_new,
    tp_er,
    tp_prot,
    xyz_cb_new,
    xyz_ce,
):
    xyz_ce_2d = {k: v[1:] for k, v in xyz_ce.items()}
    xyz_cb_new_2d = {k: v[1:] for k, v in xyz_cb_new.items()}
    xm_prot_2d, em_prot_2d, _ = pe.build_cost_matrices(tn_prot, tp_prot, xyz_ce_2d, protein_exp, prot_sel)
    xm_er_2d, em_er_2d, _ = pe.build_cost_matrices(tn_er, tp_er, xyz_ce_2d, ce_rna, rna_sel)
    xm_cb_2d, em_cb_2d, _ = pe.build_cost_matrices(tn_cb_new, tp_cb_new, xyz_cb_new_2d, cb_rna, rna_sel)
    prot_2d_groups = pe.build_cousin_groups(tn_prot, gp_map)
    prot_2d_rs = pe.compute_cousin_random_stats(xm_prot_2d, em_prot_2d, prot_2d_groups, n_random=500, seed=42)
    er_2d_groups = pe.build_cousin_groups(tn_er, gp_map)
    er_2d_rs = pe.compute_cousin_random_stats(xm_er_2d, em_er_2d, er_2d_groups, n_random=500, seed=42)
    cb_2d_groups = pe.build_cousin_groups(tn_cb_new, gp_map)
    cb_2d_rs = pe.compute_cousin_random_stats(xm_cb_2d, em_cb_2d, cb_2d_groups, n_random=500, seed=42)
    print(f"2D xyz_std: Prot={prot_2d_rs['xyz_std']:.1f}, RNA={er_2d_rs['xyz_std']:.1f}, CB={cb_2d_rs['xyz_std']:.1f}")
    res_2d_prot = pe.pareto_sweep(xm_prot_2d, em_prot_2d, tp_prot)
    res_2d_er   = pe.pareto_sweep(xm_er_2d,   em_er_2d,   tp_er)
    res_2d_cb   = pe.pareto_sweep(xm_cb_2d,   em_cb_2d,   tp_cb_new)
    _2d_rows = [
        dict(Config='CE Prot 3D', N=len(tn_prot), ExprER=res_prot_xs['expr_opt_er'], MaxER=res_prot_xs['max_er'], SpatER=res_prot_xs['spatial_opt_er']),
        dict(Config='CE Prot 2D', N=len(tn_prot), ExprER=res_2d_prot['expr_opt_er'], MaxER=res_2d_prot['max_er'], SpatER=res_2d_prot['spatial_opt_er']),
        dict(Config='CE RNA 3D',  N=len(tn_er),   ExprER=res_er_xs['expr_opt_er'], MaxER=res_er_xs['max_er'], SpatER=res_er_xs['spatial_opt_er']),
        dict(Config='CE RNA 2D',  N=len(tn_er),   ExprER=res_2d_er['expr_opt_er'], MaxER=res_2d_er['max_er'], SpatER=res_2d_er['spatial_opt_er']),
        dict(Config='CB RNA 3D',  N=len(tn_cb_new), ExprER=res_cb_new['expr_opt_er'], MaxER=res_cb_new['max_er'], SpatER=res_cb_new['spatial_opt_er']),
        dict(Config='CB RNA 2D',  N=len(tn_cb_new), ExprER=res_2d_cb['expr_opt_er'], MaxER=res_2d_cb['max_er'], SpatER=res_2d_cb['spatial_opt_er']),
    ]
    pd.DataFrame(_2d_rows).set_index('Config').round(3)
    return (
        cb_2d_rs,
        em_cb_2d,
        em_er_2d,
        em_prot_2d,
        er_2d_rs,
        prot_2d_rs,
        res_2d_cb,
        res_2d_er,
        res_2d_prot,
        xm_cb_2d,
        xm_er_2d,
        xm_prot_2d,
    )


@app.cell
def _(
    Line2D,
    cb_2d_rs,
    em_cb_2d,
    em_er_2d,
    em_prot_2d,
    er_2d_rs,
    markers,
    pe,
    plt,
    prot_2d_rs,
    tp_cb_new,
    tp_er,
    tp_prot,
    xm_cb_2d,
    xm_er_2d,
    xm_prot_2d,
):
    def _plot():
        ds = [
            ('C. elegans Protein (2D)', xm_prot_2d, em_prot_2d, tp_prot, '#1976D2', prot_2d_rs),
            ('C. elegans RNA (2D)',     xm_er_2d,   em_er_2d,   tp_er,   '#4CAF50', er_2d_rs),
            ('C. briggsae RNA (2D)',    xm_cb_2d,   em_cb_2d,   tp_cb_new,'#1B5E20', cb_2d_rs),
        ]
        fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(16, 7))
        for name, xm, em, tp, color, rs in ds:
            xa, ea, edge, kp = pe.compute_std_scaled_pareto(xm, em, tp, rs)
            lx, le = pe.lineage_std_position(xm, em, rs); n = len(tp)
            nx, ny = pe.get_null_cloud(rs, max_points=250)
            ax_l.scatter(nx, ny, c=color, s=8, alpha=0.10, zorder=0)
            ax_l.plot(xa, ea, color=color, lw=1.8, alpha=0.9, label=f'{name} (n={n})')
            ax_l.scatter(lx, le, color=color, marker='X', s=80, edgecolors='white', lw=0.8, zorder=6)
            ax_r.plot(xa, edge, color=color, lw=2, alpha=0.85, label=f'{name} (n={n})')
            for pk, mk in markers.items():
                _i = kp[pk]
                ax_r.scatter(xa[_i], edge[_i], color=color, marker=mk, s=50, zorder=8, edgecolors='white', lw=0.5)
                ax_r.annotate(f'{edge[_i]:.3f}', (xa[_i], edge[_i]+0.04), color=color, fontsize=7, ha='center', fontweight='bold')
        ax_l.axhline(0, color='gray', lw=0.8, ls=':', alpha=0.4); ax_l.axvline(0, color='gray', lw=0.8, ls=':', alpha=0.4)
        ax_l.set_xlabel('Spatial Cost 2D (σ from null)'); ax_l.set_ylabel('Expression Cost (σ from null)')
        ax_l.set_title('Cross-Species Pareto Fronts — 2D', fontsize=12)
        ax_l.legend(fontsize=8, loc='upper right'); ax_l.grid(True, alpha=0.3)
        ax_r.set_xlabel('Spatial Cost 2D (σ from null)'); ax_r.set_ylabel('Edge Retention Ratio')
        ax_r.set_title('Edge Retention — 2D', fontsize=12)
        ax_r.set_ylim(0, 1.05); ax_r.grid(True, alpha=0.3)
        ax_r.legend(handles=list(ax_r.get_lines()) + [
            Line2D([0],[0], color='black', marker='s', ls='', markersize=6, label='Expr opt'),
            Line2D([0],[0], color='black', marker='D', ls='', markersize=6, label='Max ER'),
            Line2D([0],[0], color='black', marker='o', ls='', markersize=6, label='Spatial opt'),
        ], fontsize=7.5, loc='best')
        fig.suptitle('Cross-Species — 2D XY-plane', fontsize=14, fontweight='bold')
        fig.tight_layout()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'cross_species_2d.png'), dpi=150, bbox_inches='tight')
        return fig
    _fig = _plot()
    return _fig


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ---
    ### 2b. Z-Axis Noise Null Model

    Replace each real z-coordinate with Gaussian noise calibrated to xy spatial scale.
    If real 3D ~ 2D+z_noise then the z-axis is noise; if 3D >> 2D+z_noise it carries signal.
    """)
    return


@app.cell
def _(
    em_cb_new,
    em_er,
    em_prot,
    gp_map,
    lineage_data,
    pd,
    pe,
    res_2d_cb,
    res_2d_er,
    res_2d_prot,
    res_cb_new,
    res_er_xs,
    res_prot_xs,
    tn_cb_new,
    tn_er,
    tn_prot,
    tp_cb_new,
    tp_er,
    tp_prot,
    xyz_cb_new,
    xyz_ce,
):
    print("Running z-noise null model (10 draws each)...")
    znoise_prot = pe.run_z_noise_pareto(xyz_ce, tn_prot, tp_prot, em_prot, lineage_data, gp_map, n_draws=10)
    znoise_er   = pe.run_z_noise_pareto(xyz_ce, tn_er,   tp_er,   em_er,   lineage_data, gp_map, n_draws=10)
    znoise_cb   = pe.run_z_noise_pareto(xyz_cb_new, tn_cb_new, tp_cb_new, em_cb_new, lineage_data, gp_map, n_draws=10)
    print("Done.")
    _z_rows = [
        dict(Config='CE Prot 3D', MaxER=res_prot_xs['max_er'], ExprER=res_prot_xs['expr_opt_er'], SpatER=res_prot_xs['spatial_opt_er']),
        dict(Config='CE Prot 2D', MaxER=res_2d_prot['max_er'], ExprER=res_2d_prot['expr_opt_er'], SpatER=res_2d_prot['spatial_opt_er']),
        dict(Config='CE Prot 2D+zN', MaxER=f"{znoise_prot['max_er'][0]:.3f}+-{znoise_prot['max_er'][1]:.3f}", ExprER=f"{znoise_prot['expr_er'][0]:.3f}+-{znoise_prot['expr_er'][1]:.3f}", SpatER=f"{znoise_prot['spatial_er'][0]:.3f}+-{znoise_prot['spatial_er'][1]:.3f}"),
        dict(Config='CE RNA 3D',  MaxER=res_er_xs['max_er'], ExprER=res_er_xs['expr_opt_er'], SpatER=res_er_xs['spatial_opt_er']),
        dict(Config='CE RNA 2D',  MaxER=res_2d_er['max_er'], ExprER=res_2d_er['expr_opt_er'], SpatER=res_2d_er['spatial_opt_er']),
        dict(Config='CE RNA 2D+zN', MaxER=f"{znoise_er['max_er'][0]:.3f}+-{znoise_er['max_er'][1]:.3f}", ExprER=f"{znoise_er['expr_er'][0]:.3f}+-{znoise_er['expr_er'][1]:.3f}", SpatER=f"{znoise_er['spatial_er'][0]:.3f}+-{znoise_er['spatial_er'][1]:.3f}"),
        dict(Config='CB RNA 3D',  MaxER=res_cb_new['max_er'], ExprER=res_cb_new['expr_opt_er'], SpatER=res_cb_new['spatial_opt_er']),
        dict(Config='CB RNA 2D',  MaxER=res_2d_cb['max_er'], ExprER=res_2d_cb['expr_opt_er'], SpatER=res_2d_cb['spatial_opt_er']),
        dict(Config='CB RNA 2D+zN', MaxER=f"{znoise_cb['max_er'][0]:.3f}+-{znoise_cb['max_er'][1]:.3f}", ExprER=f"{znoise_cb['expr_er'][0]:.3f}+-{znoise_cb['expr_er'][1]:.3f}", SpatER=f"{znoise_cb['spatial_er'][0]:.3f}+-{znoise_cb['spatial_er'][1]:.3f}"),
    ]
    pd.DataFrame(_z_rows).set_index('Config')
    return znoise_cb, znoise_er, znoise_prot


@app.cell
def _(
    cb_2d_rs,
    cb_new_rs,
    em_cb_2d,
    em_cb_new,
    em_er,
    em_er_2d,
    em_prot,
    em_prot_2d,
    er_2d_rs,
    er_rs,
    pe,
    plt,
    prot_2d_rs,
    prot_rs,
    tp_cb_new,
    tp_er,
    tp_prot,
    xm_cb_2d,
    xm_cb_new,
    xm_er,
    xm_er_2d,
    xm_prot,
    xm_prot_2d,
    znoise_cb,
    znoise_er,
    znoise_prot,
):
    def _plot():
        _p3 = pe.compute_std_scaled_pareto(xm_prot, em_prot, tp_prot, prot_rs)
        _e3 = pe.compute_std_scaled_pareto(xm_er, em_er, tp_er, er_rs)
        _c3 = pe.compute_std_scaled_pareto(xm_cb_new, em_cb_new, tp_cb_new, cb_new_rs)
        _p2 = pe.compute_std_scaled_pareto(xm_prot_2d, em_prot_2d, tp_prot, prot_2d_rs)
        _e2 = pe.compute_std_scaled_pareto(xm_er_2d, em_er_2d, tp_er, er_2d_rs)
        _c2 = pe.compute_std_scaled_pareto(xm_cb_2d, em_cb_2d, tp_cb_new, cb_2d_rs)
        cfgs = [
            ('CE Prot (T=255)', _p3, _p2, znoise_prot, '#1976D2'),
            ('CE RNA (T=255)',  _e3, _e2, znoise_er,   '#4CAF50'),
            ('CB RNA (T=143)',  _c3, _c2, znoise_cb,   '#1B5E20'),
        ]
        fig, axes = plt.subplots(2, 3, figsize=(21, 10), sharex='col')
        for col, (title, r3d, r2d, rnz, color) in enumerate(cfgs):
            ax_p, ax_e = axes[0, col], axes[1, col]
            n = len(r3d[2])
            ax_p.plot(r3d[0], r3d[1], color=color, lw=2.5, alpha=0.9, ls='-', label='3D')
            ax_e.plot(r3d[0], r3d[2], color=color, lw=2.5, alpha=0.9, ls='-')
            ax_p.plot(r2d[0], r2d[1], color=color, lw=1.8, alpha=0.7, ls='--', label='2D')
            ax_e.plot(r2d[0], r2d[2], color=color, lw=1.8, alpha=0.7, ls='--')
            ax_p.plot(rnz['median_xyz_arr'], rnz['median_exp_arr'], color='gray', lw=1.5, alpha=0.7, ls=':', label=f'2D+z-noise')
            ax_e.plot(rnz['median_xyz_arr'], rnz['median_edge_arr'], color='gray', lw=1.5, alpha=0.7, ls=':')
            ax_p.axhline(0, color='gray', lw=0.6, ls=':', alpha=0.3); ax_p.axvline(0, color='gray', lw=0.6, ls=':', alpha=0.3)
            ax_p.set_title(f'{title} (n={n})', fontsize=11, color=color, fontweight='bold')
            ax_p.legend(fontsize=7, loc='upper right'); ax_p.grid(True, alpha=0.3)
            ax_p.set_ylabel('Expression Cost (σ from null)')
            ax_e.set_xlabel('Spatial Cost (σ from null)'); ax_e.set_ylabel('Edge Retention')
            ax_e.set_ylim(0, 1.05); ax_e.grid(True, alpha=0.3)
        fig.suptitle('Z-Axis Noise Null — 3D vs 2D vs 2D+Gaussian(z)', fontsize=14, fontweight='bold')
        fig.tight_layout()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'z_noise_null.png'), dpi=150, bbox_inches='tight')
        return fig
    _fig = _plot()
    return _fig


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ---
    ### 2c. Tree-Aware Lineage Distance

    **Tree-Weighted Edge Retention (TWER)** measures how far in tree-space each
    reassigned edge moves. Unlike flat ER (binary preserved/changed), TWER weights
    each change by the lineage-tree path length between the real parent and the
    Pareto-assigned parent. A first-cousin swap (tree dist=2) barely affects TWER;
    a cross-branch reassignment (tree dist≥6) heavily penalises it.

    **Edge Change Profile** tracks *which* edges change at each α step along the
    Pareto sweep, revealing whether deep structural changes happen early or late.
    """)
    return


@app.cell
def _(
    cb_new_rs,
    em_cb_new,
    em_er,
    em_prot,
    er_rs,
    lm,
    prot_rs,
    tn_cb_new,
    tn_er,
    tn_prot,
    tp_cb_new,
    tp_er,
    tp_prot,
    tree_index,
    xm_cb_new,
    xm_er,
    xm_prot,
):
    _ITER = 300
    _twer_prot = lm.combined_lineage_proximity(
        xm_prot, em_prot, tp_prot, tn_prot, tree_index, prot_rs, iteration=_ITER, n_random=100)
    _twer_er   = lm.combined_lineage_proximity(
        xm_er,   em_er,   tp_er,   tn_er,   tree_index, er_rs,   iteration=_ITER, n_random=100)
    _twer_cb   = lm.combined_lineage_proximity(
        xm_cb_new, em_cb_new, tp_cb_new, tn_cb_new, tree_index, cb_new_rs, iteration=_ITER, n_random=100)
    print(f"TWER at max ER: Prot={1-_twer_prot['twer'][_twer_prot['kp']['max_er_idx']]:.3f}  "
          f"RNA={1-_twer_er['twer'][_twer_er['kp']['max_er_idx']]:.3f}  "
          f"CB={1-_twer_cb['twer'][_twer_cb['kp']['max_er_idx']]:.3f}")
    twer_data = _twer_prot, _twer_er, _twer_cb
    return (twer_data,)


@app.cell
def _(Line2D, _save_plot, plt, twer_data):
    def _plot():
        _configs = [
            ('C. elegans Protein', twer_data[0], '#1976D2'),
            ('C. elegans RNA',     twer_data[1], '#4CAF50'),
            ('C. briggsae RNA',    twer_data[2], '#1B5E20'),
        ]
        fig, axes = plt.subplots(1, 3, figsize=(20, 6))

        for col, (title, twr, color) in enumerate(_configs):
            ax = axes[col]
            _er = twr['traditional_er']
            _twer = 1 - twr['twer']
            _xa = twr['xyz_arr']
            _n = twr['n']
            _kp = twr['kp']

            # TWER curve
            ax.plot(_xa, _twer, color=color, lw=2.5, alpha=0.9, label='TWER (tree-weighted)')
            # Traditional ER curve for comparison
            ax.plot(_xa, _er, color='gray', lw=1.2, alpha=0.5, ls='--', label='ER (flat)')

            # Mark key points
            for pk, mk in [('expr_opt_idx', 's'), ('max_er_idx', 'D'), ('spatial_opt_idx', 'o')]:
                _i = _kp[pk]
                ax.scatter(_xa[_i], _twer[_i], color=color, marker=mk, s=50, zorder=8,
                          edgecolors='white', lw=0.5)
                ax.scatter(_xa[_i], _er[_i], color='gray', marker=mk, s=30, zorder=8,
                          edgecolors='white', lw=0.3, alpha=0.7)

            ax.axvline(0, color='gray', lw=0.6, ls=':', alpha=0.3)
            ax.set_xlabel('Spatial Cost (σ from null)')
            ax.set_ylabel('Retention Score')
            ax.set_title(f'{title} (n={_n})', fontsize=11, color=color, fontweight='bold')
            ax.set_ylim(0, 1.05)
            ax.legend(fontsize=7.5, loc='best')
            ax.grid(True, alpha=0.3)

        _handles = [
            Line2D([0],[0], color='black', marker='s', ls='', markersize=6, label='Expr opt'),
            Line2D([0],[0], color='black', marker='D', ls='', markersize=6, label='Max ER'),
            Line2D([0],[0], color='black', marker='o', ls='', markersize=6, label='Spatial opt'),
        ]
        fig.legend(handles=_handles, fontsize=8, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.02))
        fig.suptitle('Tree-Weighted Edge Retention (TWER) vs Traditional ER',
                     fontsize=14, fontweight='bold')
        fig.tight_layout()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'twer_vs_er.png'), dpi=150, bbox_inches='tight')
        return fig
    _fig = _plot()
    return _fig


@app.cell
def _(np, plt, twer_data):
    def _plot():
        _configs = [
            ('C. elegans Protein', twer_data[0], '#1976D2'),
            ('C. elegans RNA',     twer_data[1], '#4CAF50'),
            ('C. briggsae RNA',    twer_data[2], '#1B5E20'),
        ]
        fig, axes = plt.subplots(1, 3, figsize=(20, 6))

        for col, (title, twr, color) in enumerate(_configs):
            ax = axes[col]
            _ep = twr['edge_profile']
            _steps = _ep['steps']        # α values where changes occurred
            _mean_dists = _ep['mean_tree_dists']
            _xa = twr['xyz_arr']         # spatial cost (σ) at each α
            _n = twr['n']

            if len(_steps) > 0:
                # Map α → spatial cost for x-axis
                _n_pts = len(_xa) - 1  # iteration count
                _x_spatial = np.array([_xa[int(s * _n_pts)] for s in _steps])

                # Scatter: spatial cost vs mean tree distance of changes
                ax.scatter(_x_spatial, _mean_dists, color=color,
                           s=30, alpha=0.7, edgecolors='white', lw=0.3)

                # Smooth trend
                if len(_steps) > 5:
                    _window = max(len(_steps) // 10, 3)
                    _smooth = _mean_dists.copy()
                    for _j in range(len(_mean_dists) - _window + 1):
                        _smooth[_j + _window//2] = _mean_dists[_j:_j+_window].mean()
                    ax.plot(_x_spatial, _smooth, color='black', lw=1.5, alpha=0.6, ls='-')

                ax.axhline(2, color='green', lw=0.8, ls=':', alpha=0.5, label='Cousin-level (≤2)')
                ax.axhline(4, color='orange', lw=0.8, ls=':', alpha=0.5, label='Close relative (≤4)')

                # Split by x-position (left = near expression-opt ≈ low σ,
                # right = near spatial-opt ≈ high σ) and annotate at data coords
                _mid_x = (_x_spatial.min() + _x_spatial.max()) / 2
                _left = _x_spatial < _mid_x
                _right = ~_left
                _left_m = _mean_dists[_left].mean() if _left.any() else 0
                _right_m = _mean_dists[_right].mean() if _right.any() else 0
                _ylim = ax.get_ylim()
                _ypos = _ylim[1] * 0.92
                if _left.any():
                    ax.annotate(f'Spatial-side mean: {_left_m:.1f}',
                               xy=(_x_spatial[_left].mean(), _left_m),
                               xytext=(_x_spatial[_left].mean(), _ypos),
                               fontsize=8, color='darkgreen', ha='center',
                               arrowprops=dict(arrowstyle='->', color='darkgreen', lw=1))
                if _right.any():
                    ax.annotate(f'Expr-side mean: {_right_m:.1f}',
                               xy=(_x_spatial[_right].mean(), _right_m),
                               xytext=(_x_spatial[_right].mean(), _ypos * 0.88),
                               fontsize=8, color='darkred', ha='center',
                               arrowprops=dict(arrowstyle='->', color='darkred', lw=1))

            ax.axvline(0, color='gray', lw=0.6, ls=':', alpha=0.3)
            ax.set_xlabel('Spatial Cost (σ from null)')
            ax.set_ylabel('Mean Tree Distance of Changes')
            ax.set_title(f'{title} — Edge Change Profile (n={_n})',
                        fontsize=11, color=color, fontweight='bold')
            ax.legend(fontsize=7, loc='upper right')
            ax.grid(True, alpha=0.3)

        fig.suptitle('Edge Change Profile along Pareto Sweep — Mean Tree Distance of Changed Edges',
                     fontsize=14, fontweight='bold')
        fig.tight_layout()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'edge_change_profile.png'), dpi=150, bbox_inches='tight')
        return fig
    _fig = _plot()
    return _fig


@app.cell
def _(np, plt, twer_data):
    def _plot():
        _configs = [
            ('C. elegans Protein', twer_data[0], '#1976D2'),
            ('C. elegans RNA',     twer_data[1], '#4CAF50'),
            ('C. briggsae RNA',    twer_data[2], '#1B5E20'),
        ]
        fig, axes = plt.subplots(1, 3, figsize=(20, 6))

        for col, (title, twr, color) in enumerate(_configs):
            ax = axes[col]
            _lmd = twr['lineage_mean_dist']
            _xa = twr['xyz_arr']
            _null = twr['lineage_null']
            _n = twr['n']

            # Pareto curve
            ax.plot(_xa, _lmd, color=color, lw=2.5, alpha=0.9, label='Pareto front')
            # Null baselines
            _fm = _null['full_mean']
            _cm = _null['cousin_mean']
            ax.axhline(_fm, color='darkred', lw=1.2, ls='--', alpha=0.7,
                       label=f'Full random ({_fm:.1f})')
            ax.axhline(_cm, color='darkgreen', lw=1.2, ls='--', alpha=0.7,
                       label=f'Cousin shuffle ({_cm:.1f})')

            ax.axvline(0, color='gray', lw=0.6, ls=':', alpha=0.3)
            ax.set_xlabel('Spatial Cost (σ from null)')
            ax.set_ylabel('Mean Tree Distance from Real Lineage')
            ax.set_title(f'{title} (n={_n})', fontsize=11, color=color, fontweight='bold')
            ax.legend(fontsize=7, loc='upper right')
            ax.grid(True, alpha=0.3)

        fig.suptitle('Mean Lineage Tree Distance along Pareto Front\n(null baselines from 100k Monte Carlo samples)',
                     fontsize=14, fontweight='bold')
        fig.tight_layout()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'lineage_distance_along_pareto.png'), dpi=150, bbox_inches='tight')
        return fig
    _fig = _plot()
    return _fig


@app.cell
def _(np, plt, twer_data):
    def _plot():
        _configs = [
            ('C. elegans Protein', twer_data[0], '#1976D2'),
            ('C. elegans RNA',     twer_data[1], '#4CAF50'),
            ('C. briggsae RNA',    twer_data[2], '#1B5E20'),
        ]
        fig, axes = plt.subplots(1, 3, figsize=(20, 6))

        for col, (title, twr, color) in enumerate(_configs):
            ax = axes[col]
            _ct = twr['cost_tree_tradeoff']
            _td = _ct['tree_dists']
            _dx = _ct['delta_xyz']
            _de = _ct['delta_exp']
            _n = twr['n']

            # Cost savings vs tree distance
            ax.plot(_td, _dx, color='#1976D2', lw=2, alpha=0.85, label='Δ xyz cost')
            ax.plot(_td, _de, color='#E65100', lw=2, alpha=0.85, label='Δ exp cost')
            ax.axhline(0, color='gray', lw=0.8, ls=':', alpha=0.5)
            # Mark real lineage at origin
            ax.scatter([0], [0], color='black', marker='X', s=80, zorder=6,
                      label='Real lineage')

            # Annotate direction: which end is which
            ax.annotate('← spatial-opt', xy=(_td[0], _dx[0]), fontsize=7, color='gray',
                       ha='right', va='bottom')
            ax.annotate('expr-opt →', xy=(_td[-1], _dx[-1]), fontsize=7, color='gray',
                       ha='left', va='top')

            ax.set_xlabel('Mean Tree Distance from Real Lineage')
            ax.set_ylabel('Δ Cost from Real Lineage (σ)')
            ax.set_title(f'{title} (n={_n})', fontsize=11, color=color, fontweight='bold')
            ax.legend(fontsize=7, loc='best')
            ax.grid(True, alpha=0.3)

        fig.suptitle('Cost Savings vs Tree Disruption along Pareto Front\n'
                     '(positive Δ = better than real lineage; negative = worse)',
                     fontsize=14, fontweight='bold')
        fig.tight_layout()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'cost_vs_tree_tradeoff.png'), dpi=150, bbox_inches='tight')
        return fig
    _fig = _plot()
    return _fig


@app.cell
def _(np, plt, twer_data):
    def _plot():
        _configs = [
            ('C. elegans Protein', twer_data[0], '#1976D2'),
            ('C. elegans RNA',     twer_data[1], '#4CAF50'),
            ('C. briggsae RNA',    twer_data[2], '#1B5E20'),
        ]
        fig, axes = plt.subplots(1, 3, figsize=(20, 6))

        for col, (title, twr, color) in enumerate(_configs):
            ax = axes[col]
            _ct = twr['cost_tree_tradeoff']
            _xa = twr['xyz_arr']
            _mx = _ct['marginal_xyz']
            _me = _ct['marginal_exp']
            _mc = _ct['marginal_combined']
            _n = twr['n']

            # Marginal efficiency along the front
            # Smooth with a small window for readability
            _w = max(len(_mx) // 20, 3)
            _smooth_mx = np.convolve(_mx, np.ones(_w)/_w, mode='same')
            _smooth_me = np.convolve(_me, np.ones(_w)/_w, mode='same')
            _smooth_mc = np.convolve(_mc, np.ones(_w)/_w, mode='same')

            ax.plot(_xa, _smooth_mx, color='#1976D2', lw=1.5, alpha=0.7, ls='--',
                   label='Marginal xyz')
            ax.plot(_xa, _smooth_me, color='#E65100', lw=1.5, alpha=0.7, ls='--',
                   label='Marginal exp')
            ax.plot(_xa, _smooth_mc, color='black', lw=2, alpha=0.9,
                   label='Combined')

            ax.axhline(0, color='gray', lw=0.6, ls=':', alpha=0.4)
            ax.axvline(0, color='gray', lw=0.6, ls=':', alpha=0.3)
            ax.set_xlabel('Spatial Cost (σ from null)')
            ax.set_ylabel('Δ Savings per Δ Tree Distance (σ/edge)')
            ax.set_title(f'{title} (n={_n})', fontsize=11, color=color, fontweight='bold')
            ax.legend(fontsize=7, loc='best')
            ax.grid(True, alpha=0.3)

        fig.suptitle('Marginal Efficiency: Cost Savings per Unit of Tree Disruption\n'
                     '(higher = more savings per edge of tree change)',
                     fontsize=14, fontweight='bold')
        fig.tight_layout()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'marginal_efficiency.png'), dpi=150, bbox_inches='tight')
        return fig
    _fig = _plot()
    return _fig


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ---
    ## 3. Subtree Analysis
    """)
    return


@app.cell
def _(
    dl,
    gp_map,
    lineage_data,
    pd,
    pe,
    prot_sel,
    protein_exp,
    tn_prot,
    tp_prot,
    v_prot,
    xyz_ce,
):
    SUBTREES = ['AB', 'ABa', 'ABp', 'P1']
    def _c():
        _ps = {}
        for _sr in SUBTREES + ['Full']:
            _lb = _sr if _sr != 'Full' else 'Full tree'
            if _sr == 'Full': _tn, _tp = tn_prot, tp_prot
            else: _tn, _tp = dl.collect_terminals(lineage_data, v_prot, subtree=_sr)
            _xm, _em, _ = pe.build_cost_matrices(_tn, _tp, xyz_ce, protein_exp, prot_sel)
            _res = pe.pareto_sweep(_xm, _em, _tp)
            _sg = pe.build_cousin_groups(_tn, gp_map)
            _srs = pe.compute_cousin_random_stats(_xm, _em, _sg, n_random=300, seed=42)
            _xa, _ea, _, _ = pe.compute_std_scaled_pareto(_xm, _em, _tp, _srs)
            _lx, _le = pe.lineage_std_position(_xm, _em, _srs)
            _rd = pe.relative_pareto_distance(_xa, _ea, _lx, _le)
            _ps[_lb] = dict(tn=_tn, tp=_tp, xm=_xm, em=_em, res=_res, rel_pareto_dist=_rd)
            print(f"  {_lb:12s}: {len(_tn):4d} cells | Expr ER={_res['expr_opt_er']:.3f} | Max ER={_res['max_er']:.3f} | Rel Dist={_rd:.3f}")
        return _ps
    prot_sub = _c()
    prot_sub_df = pd.DataFrame({
        'Subtree': list(prot_sub.keys()), 'N': [len(prot_sub[k]['tn']) for k in prot_sub],
        'Expr ER': [prot_sub[k]['res']['expr_opt_er'] for k in prot_sub],
        'Max ER': [prot_sub[k]['res']['max_er'] for k in prot_sub],
        'Spatial ER': [prot_sub[k]['res']['spatial_opt_er'] for k in prot_sub],
        'Rel Dist': [round(prot_sub[k]['rel_pareto_dist'], 3) for k in prot_sub],
    }).set_index('Subtree')
    prot_sub_df
    return SUBTREES, prot_sub


@app.cell
def _(gp_map, markers, pe, plt, prot_sub):
    def _plot():
        sc = {'Full tree': '#333', 'AB': '#1976D2', 'ABa': '#4CAF50', 'ABp': '#FF9800', 'P1': '#E91E63'}
        fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(16, 7))
        for label in ['Full tree', 'AB', 'ABa', 'ABp', 'P1']:
            d = prot_sub[label]; color = sc[label]; n = len(d['tn'])
            sg = pe.build_cousin_groups(d['tn'], gp_map)
            sr = pe.compute_cousin_random_stats(d['xm'], d['em'], sg, n_random=300, seed=42)
            xa, ea, edge, kp = pe.compute_std_scaled_pareto(d['xm'], d['em'], d['tp'], sr)
            lx, le = pe.lineage_std_position(d['xm'], d['em'], sr)
            nx, ny = pe.get_null_cloud(sr, max_points=250)
            ax_l.scatter(nx, ny, c=color, s=8, alpha=0.10, zorder=0)
            ax_l.plot(xa, ea, color=color, lw=1.5, alpha=0.85, label=f'{label} (n={n})')
            ax_l.scatter(lx, le, color=color, marker='X', s=70, edgecolors='white', lw=0.8, zorder=6)
            ax_r.plot(xa, edge, color=color, lw=2, alpha=0.85, label=f'{label} (n={n})')
            for pk, mk in markers.items():
                ax_r.scatter(xa[kp[pk]], edge[kp[pk]], color=color, marker=mk, s=50, zorder=8, edgecolors='white', lw=0.5)
        ax_l.axhline(0, color='gray', lw=0.8, ls=':', alpha=0.4); ax_l.axvline(0, color='gray', lw=0.8, ls=':', alpha=0.4)
        ax_l.set_xlabel('Spatial Cost (σ from null)'); ax_l.set_ylabel('Expression Cost (σ from null)')
        ax_l.set_title('CE Protein — Subtree Pareto Fronts', fontsize=12)
        ax_l.legend(fontsize=8, loc='upper right'); ax_l.grid(True, alpha=0.3)
        ax_r.set_xlabel('Spatial Cost (σ from null)'); ax_r.set_ylabel('Edge Retention')
        ax_r.set_title('CE Protein — Edge Retention by Subtree', fontsize=12)
        ax_r.set_ylim(0, 1.05); ax_r.grid(True, alpha=0.3); ax_r.legend(fontsize=7.5, loc='best')
        fig.suptitle('Subtree Analysis — C. elegans Protein (top 20)', fontsize=14, fontweight='bold')
        fig.tight_layout()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'ce_protein_subtrees.png'), dpi=150, bbox_inches='tight')
        return fig
    _fig = _plot()
    return _fig


@app.cell
def _(
    SUBTREES,
    ce_rna,
    dl,
    gp_map,
    lineage_data,
    markers,
    pd,
    pe,
    plt,
    rna_sel,
    tn_er,
    tp_er,
    v_ce_rna,
    xyz_ce,
):
    def _c():
        _es = {}
        for _sr in SUBTREES + ['Full']:
            _lb = _sr if _sr != 'Full' else 'Full tree'
            if _sr == 'Full': _tn, _tp = tn_er, tp_er
            else: _tn, _tp = dl.collect_terminals(lineage_data, v_ce_rna, subtree=_sr)
            _xm, _em, _ = pe.build_cost_matrices(_tn, _tp, xyz_ce, ce_rna, rna_sel)
            _res = pe.pareto_sweep(_xm, _em, _tp)
            _sg = pe.build_cousin_groups(_tn, gp_map)
            _srs = pe.compute_cousin_random_stats(_xm, _em, _sg, n_random=300, seed=42)
            _xa, _ea, _, _ = pe.compute_std_scaled_pareto(_xm, _em, _tp, _srs)
            _lx, _le = pe.lineage_std_position(_xm, _em, _srs)
            _rd = pe.relative_pareto_distance(_xa, _ea, _lx, _le)
            _es[_lb] = dict(tn=_tn, tp=_tp, xm=_xm, em=_em, res=_res, rel_pareto_dist=_rd)
            print(f"  {_lb:12s}: {len(_tn):4d} cells | Expr ER={_res['expr_opt_er']:.3f} | Max ER={_res['max_er']:.3f} | Rel Dist={_rd:.3f}")
        return _es
    er_sub = _c()
    er_sub_df = pd.DataFrame({
        'Subtree': list(er_sub.keys()), 'N': [len(er_sub[k]['tn']) for k in er_sub],
        'Expr ER': [er_sub[k]['res']['expr_opt_er'] for k in er_sub],
        'Max ER': [er_sub[k]['res']['max_er'] for k in er_sub],
        'Spatial ER': [er_sub[k]['res']['spatial_opt_er'] for k in er_sub],
        'Rel Dist': [round(er_sub[k]['rel_pareto_dist'], 3) for k in er_sub],
    }).set_index('Subtree')
    er_sub_df

    def _plot():
        sc = {'Full tree': '#333', 'AB': '#1976D2', 'ABa': '#4CAF50', 'ABp': '#FF9800', 'P1': '#E91E63'}
        fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(16, 7))
        for label in ['Full tree', 'AB', 'ABa', 'ABp', 'P1']:
            d = er_sub[label]; color = sc[label]; n = len(d['tn'])
            sg = pe.build_cousin_groups(d['tn'], gp_map)
            sr = pe.compute_cousin_random_stats(d['xm'], d['em'], sg, n_random=300, seed=42)
            xa, ea, edge, kp = pe.compute_std_scaled_pareto(d['xm'], d['em'], d['tp'], sr)
            lx, le = pe.lineage_std_position(d['xm'], d['em'], sr)
            nx, ny = pe.get_null_cloud(sr, max_points=250)
            ax_l.scatter(nx, ny, c=color, s=8, alpha=0.10, zorder=0)
            ax_l.plot(xa, ea, color=color, lw=1.5, alpha=0.85, label=f'{label} (n={n})')
            ax_l.scatter(lx, le, color=color, marker='X', s=70, edgecolors='white', lw=0.8, zorder=6)
            ax_r.plot(xa, edge, color=color, lw=2, alpha=0.85, label=f'{label} (n={n})')
            for pk, mk in markers.items():
                ax_r.scatter(xa[kp[pk]], edge[kp[pk]], color=color, marker=mk, s=50, zorder=8, edgecolors='white', lw=0.5)
        ax_l.axhline(0, color='gray', lw=0.8, ls=':', alpha=0.4); ax_l.axvline(0, color='gray', lw=0.8, ls=':', alpha=0.4)
        ax_l.set_xlabel('Spatial Cost (σ from null)'); ax_l.set_ylabel('Expression Cost (σ from null)')
        ax_l.set_title('CE RNA — Subtree Pareto Fronts', fontsize=12)
        ax_l.legend(fontsize=8, loc='upper right'); ax_l.grid(True, alpha=0.3)
        ax_r.set_xlabel('Spatial Cost (σ from null)'); ax_r.set_ylabel('Edge Retention')
        ax_r.set_title('CE RNA — Edge Retention by Subtree', fontsize=12)
        ax_r.set_ylim(0, 1.05); ax_r.grid(True, alpha=0.3); ax_r.legend(fontsize=7.5, loc='best')
        fig.suptitle('Subtree Analysis — C. elegans RNA', fontsize=14, fontweight='bold')
        fig.tight_layout()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'ce_rna_subtrees.png'), dpi=150, bbox_inches='tight')
        return fig
    _plot()
    return (er_sub,)


@app.cell
def _(
    SUBTREES,
    cb_rna,
    dl,
    gp_map,
    lineage_data,
    markers,
    pd,
    pe,
    plt,
    rna_sel,
    tn_cb_new,
    tp_cb_new,
    v_cb_rna_new,
    xyz_cb_new,
):
    def _c():
        _cs = {}
        for _sr in SUBTREES + ['Full']:
            _lb = _sr if _sr != 'Full' else 'Full tree'
            if _sr == 'Full': _tn, _tp = tn_cb_new, tp_cb_new
            else: _tn, _tp = dl.collect_terminals(lineage_data, v_cb_rna_new, subtree=_sr)
            _xm, _em, _ = pe.build_cost_matrices(_tn, _tp, xyz_cb_new, cb_rna, rna_sel)
            _res = pe.pareto_sweep(_xm, _em, _tp)
            _sg = pe.build_cousin_groups(_tn, gp_map)
            _srs = pe.compute_cousin_random_stats(_xm, _em, _sg, n_random=300, seed=42)
            _xa, _ea, _, _ = pe.compute_std_scaled_pareto(_xm, _em, _tp, _srs)
            _lx, _le = pe.lineage_std_position(_xm, _em, _srs)
            _rd = pe.relative_pareto_distance(_xa, _ea, _lx, _le)
            _cs[_lb] = dict(tn=_tn, tp=_tp, xm=_xm, em=_em, res=_res, rel_pareto_dist=_rd)
            print(f"  {_lb:12s}: {len(_tn):4d} cells | Expr ER={_res['expr_opt_er']:.3f} | Max ER={_res['max_er']:.3f} | Rel Dist={_rd:.3f}")
        return _cs
    cb_sub = _c()
    cb_sub_df = pd.DataFrame({
        'Subtree': list(cb_sub.keys()), 'N': [len(cb_sub[k]['tn']) for k in cb_sub],
        'Expr ER': [cb_sub[k]['res']['expr_opt_er'] for k in cb_sub],
        'Max ER': [cb_sub[k]['res']['max_er'] for k in cb_sub],
        'Spatial ER': [cb_sub[k]['res']['spatial_opt_er'] for k in cb_sub],
        'Rel Dist': [round(cb_sub[k]['rel_pareto_dist'], 3) for k in cb_sub],
    }).set_index('Subtree')
    cb_sub_df

    def _plot():
        sc = {'Full tree': '#333', 'AB': '#1976D2', 'ABa': '#4CAF50', 'ABp': '#FF9800', 'P1': '#E91E63'}
        fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(16, 7))
        for label in ['Full tree', 'AB', 'ABa', 'ABp', 'P1']:
            d = cb_sub[label]; color = sc[label]; n = len(d['tn'])
            sg = pe.build_cousin_groups(d['tn'], gp_map)
            sr = pe.compute_cousin_random_stats(d['xm'], d['em'], sg, n_random=300, seed=42)
            xa, ea, edge, kp = pe.compute_std_scaled_pareto(d['xm'], d['em'], d['tp'], sr)
            lx, le = pe.lineage_std_position(d['xm'], d['em'], sr)
            nx, ny = pe.get_null_cloud(sr, max_points=250)
            ax_l.scatter(nx, ny, c=color, s=8, alpha=0.10, zorder=0)
            ax_l.plot(xa, ea, color=color, lw=1.5, alpha=0.85, label=f'{label} (n={n})')
            ax_l.scatter(lx, le, color=color, marker='X', s=70, edgecolors='white', lw=0.8, zorder=6)
            ax_r.plot(xa, edge, color=color, lw=2, alpha=0.85, label=f'{label} (n={n})')
            for pk, mk in markers.items():
                ax_r.scatter(xa[kp[pk]], edge[kp[pk]], color=color, marker=mk, s=50, zorder=8, edgecolors='white', lw=0.5)
        ax_l.axhline(0, color='gray', lw=0.8, ls=':', alpha=0.4); ax_l.axvline(0, color='gray', lw=0.8, ls=':', alpha=0.4)
        ax_l.set_xlabel('Spatial Cost (σ from null)'); ax_l.set_ylabel('Expression Cost (σ from null)')
        ax_l.set_title('CB RNA — Subtree Pareto Fronts', fontsize=12)
        ax_l.legend(fontsize=8, loc='upper right'); ax_l.grid(True, alpha=0.3)
        ax_r.set_xlabel('Spatial Cost (σ from null)'); ax_r.set_ylabel('Edge Retention')
        ax_r.set_title('CB RNA — Edge Retention by Subtree', fontsize=12)
        ax_r.set_ylim(0, 1.05); ax_r.grid(True, alpha=0.3); ax_r.legend(fontsize=7.5, loc='best')
        fig.suptitle('Subtree Analysis — C. briggsae RNA', fontsize=14, fontweight='bold')
        fig.tight_layout()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'cb_rna_subtrees.png'), dpi=150, bbox_inches='tight')
        return fig
    _plot()
    return (cb_sub,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ---
    ## 4. Edge Perturbation Analysis

    How close is the real lineage to a local optimum?  Enumerate all C(n,2) and C(n,3)
    edge swaps; count what fraction improve spatial, expression, or both costs.
    """)
    return


@app.cell
def _(em_cb_new, em_er, em_prot, pd, pe, prot_sub, xm_cb_new, xm_er, xm_prot):
    pert_results = []
    for label, xm, em in [('CE Prot (T=255)', xm_prot, em_prot), ('CE RNA (T=255)', xm_er, em_er), ('CB RNA (T=143)', xm_cb_new, em_cb_new)]:
        pert_results.append(pe.run_perturbation_test(xm, em, label))
    for _lb in ['Full tree', 'AB', 'ABa', 'ABp', 'P1']:
        _d = prot_sub[_lb]
        pert_results.append(pe.run_perturbation_test(_d['xm'], _d['em'], f'CE Prot Subtree: {_lb}'))
    pert_rows = [{ 'Config': _r['label'], 'N': _r['n'],
        'C2 Total': _r['c2']['total'], 'C2 Save XYZ%': _r['c2']['pct_xyz'],
        'C2 Save Exp%': _r['c2']['pct_exp'], 'C2 Save Both%': _r['c2']['pct_both'],
        'C3 Total': _r['c3']['total'], 'C3 Save XYZ%': _r['c3']['pct_xyz'],
        'C3 Save Exp%': _r['c3']['pct_exp'], 'C3 Save Both%': _r['c3']['pct_both'],
    } for _r in pert_results]
    pert_df = pd.DataFrame(pert_rows).set_index('Config')
    pert_df
    return (pert_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ---
    ## 5. Random-20-Features Null Model

    To test whether our 20 selected proteins/TFs carry real lineage signal, we build
    cost matrices from random sets of 20 expression features and run the full Pareto
    sweep on each.  The real feature set should significantly outperform random.
    """)
    return


@app.cell
def _(
    Parallel,
    delayed,
    fd_cb_new,
    fd_er,
    fd_prot,
    np,
    pd,
    pe,
    plt,
    res_cb_new,
    res_er_xs,
    res_prot_xs,
    tn_cb_new,
    tn_er,
    tn_prot,
    tp_cb_new,
    tp_er,
    tp_prot,
    xm_cb_new,
    xm_er,
    xm_prot,
):
    from tqdm import tqdm as _tqdm
    _NR = 200
    print(f"Building random cost matrices (n={_NR}) for 3 configs...")

    def _rb(xm, tp, tn, fd, label):
        rmats = pe.build_random_cost_matrices(tp, tn, fd, n_features=20, n_random=_NR, seed_base=0)
        rres = Parallel(n_jobs=-1, prefer="processes")(
            delayed(pe.pareto_sweep)(xm, rmat, tp) for rmat in _tqdm(rmats, desc=f"Random: {label}")
        )
        return pd.DataFrame(rres)

    prot_rand_df = _rb(xm_prot, tp_prot, tn_prot, fd_prot, "CE Protein")
    er_rand_df   = _rb(xm_er,   tp_er,   tn_er,   fd_er,   "CE RNA")
    cb_rand_df   = _rb(xm_cb_new, tp_cb_new, tn_cb_new, fd_cb_new, "CB RNA")

    for _nm, _rdf, _real in [("CE Protein", prot_rand_df, res_prot_xs),
                              ("CE RNA",    er_rand_df,   res_er_xs),
                              ("CB RNA",    cb_rand_df,   res_cb_new)]:
        print(f"  {_nm}: Expr ER={_real['expr_opt_er']:.3f} (real) vs {_rdf['expr_opt_er'].mean():.3f}+-{_rdf['expr_opt_er'].std():.3f} (random)")
        print(f"          Max ER={_real['max_er']:.3f} (real) vs {_rdf['max_er'].mean():.3f}+-{_rdf['max_er'].std():.3f} (random)")

    def _plot():
        fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
        ds = ['CE Protein', 'CE RNA', 'CB RNA']
        cols = ['#1976D2', '#4CAF50', '#FF9800']
        real_expr = [res_prot_xs['expr_opt_er'], res_er_xs['expr_opt_er'], res_cb_new['expr_opt_er']]
        rand_expr_m = [prot_rand_df['expr_opt_er'].mean(), er_rand_df['expr_opt_er'].mean(), cb_rand_df['expr_opt_er'].mean()]
        rand_expr_s = [prot_rand_df['expr_opt_er'].std(), er_rand_df['expr_opt_er'].std(), cb_rand_df['expr_opt_er'].std()]
        real_max = [res_prot_xs['max_er'], res_er_xs['max_er'], res_cb_new['max_er']]
        rand_max_m = [prot_rand_df['max_er'].mean(), er_rand_df['max_er'].mean(), cb_rand_df['max_er'].mean()]
        rand_max_s = [prot_rand_df['max_er'].std(), er_rand_df['max_er'].std(), cb_rand_df['max_er'].std()]
        x = np.arange(3); w = 0.35
        for ax, title, rv, rm, rs in [
            (axes[0], 'Expr Opt ER (α=0)', real_expr, rand_expr_m, rand_expr_s),
            (axes[1], 'Max ER', real_max, rand_max_m, rand_max_s)]:
            ax.bar(x-w/2, rv, w, color=cols, alpha=0.9, label='Selected 20')
            ax.bar(x+w/2, rm, w, color='gray', alpha=0.6, label='Random 20', yerr=rs, capsize=4)
            ax.set_xticks(x); ax.set_xticklabels(ds, fontsize=10)
            ax.set_ylabel('Edge Retention'); ax.set_title(title, fontsize=12)
            ax.legend(fontsize=9); ax.grid(True, alpha=0.3, axis='y'); ax.set_ylim(0, 1.05)
        fig.suptitle('Real vs Random 20 Features', fontsize=14, fontweight='bold')
        fig.tight_layout()
        import os as _os
        _out = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'terminal_pareto', 'output')
        _os.makedirs(_out, exist_ok=True)
        fig.savefig(_os.path.join(_out, 'random_features_null.png'), dpi=150, bbox_inches='tight')
        return fig
    _plot()
    return cb_rand_df, er_rand_df, prot_rand_df


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ---
    ## 6. Grand Summary
    """)
    return


@app.cell
def _(
    cb_rand_df,
    cb_sub,
    er_rand_df,
    er_sub,
    pd,
    pert_df,
    prot_rand_df,
    prot_sub,
    res_cb_new,
    res_er_xs,
    res_prot_xs,
    tn_cb_new,
    tn_er,
    tn_prot,
):
    all_rows = []
    for _lb, _res, _n in [('CE Protein (T=255)', res_prot_xs, len(tn_prot)), ('CE RNA (T=255)', res_er_xs, len(tn_er)), ('CB RNA (T=143)', res_cb_new, len(tn_cb_new))]:
        all_rows.append(dict(Experiment=f'Full — {_lb}', Group='Full tree', N=_n, ExprER=_res['expr_opt_er'], MaxER=_res['max_er'], SpatER=_res['spatial_opt_er']))
    for _en, _sd in [('Subtree — CE Prot', prot_sub), ('Subtree — CE RNA', er_sub), ('Subtree — CB RNA', cb_sub)]:
        for _lb in ['AB', 'ABa', 'ABp', 'P1']:
            _d = _sd[_lb]
            all_rows.append(dict(Experiment=_en, Group=_lb, N=len(_d['tn']), ExprER=_d['res']['expr_opt_er'], MaxER=_d['res']['max_er'], SpatER=_d['res']['spatial_opt_er']))
    for _, _row in pert_df.iterrows():
        all_rows.append(dict(Experiment='Perturbation', Group=_row.name, N=int(_row['N']), C2SaveBoth=_row['C2 Save Both%'], C3SaveBoth=_row['C3 Save Both%']))
    # Random-feature summary
    for _nm, _rdf, _real in [('CE Prot Random', prot_rand_df, res_prot_xs), ('CE RNA Random', er_rand_df, res_er_xs), ('CB RNA Random', cb_rand_df, res_cb_new)]:
        all_rows.append(dict(Experiment=_nm, Group='Random 20', N=200, RandExprER=f"{_rdf['expr_opt_er'].mean():.3f}+-{_rdf['expr_opt_er'].std():.3f}", RandMaxER=f"{_rdf['max_er'].mean():.3f}+-{_rdf['max_er'].std():.3f}"))
    grand_df = pd.DataFrame(all_rows)
    grand_df.head(20)
    return


if __name__ == "__main__":
    app.run()
