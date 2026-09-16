"""Render a CE-only S3 amendment from validated tracking/common-assignment caches.

Does not recompute optimization or write accepted publication files. The
standalone caption lives in figS3_ce_tracking_variance_amendment.tex.
"""
from pathlib import Path
import hashlib
import itertools
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend import Legend
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from terminal_pareto.figS3_ce_tracking_robustness import plot_fronts, REPLICATE_STYLES
TRACK = ROOT/'terminal_pareto/output/tracking_geometry_sensitivity'
COMMON = ROOT/'terminal_pareto/output/common_tracking_assignment'
OUT = ROOT/'terminal_pareto/output/tracking_s3_amendment'
REP_COLORS = {'embryo1':'#4477AA','embryo2':'#CC6677','embryo3':'#228833'}
PAIR_COLORS = {('embryo1','embryo2'):'#4477AA',('embryo1','embryo3'):'#CC7722',
               ('embryo2','embryo3'):'#8866AA'}


def load_and_validate():
    provenance=json.loads((TRACK/'provenance.json').read_text())
    for name,digest in provenance['sha256'].items():
        assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest, name
    cache=np.load(TRACK/'assignments_and_costs.npz')
    rows=np.arange(len(cache['children']))
    transfer=pd.read_csv(TRACK/'transfer.csv')
    neighborhoods=pd.read_csv(TRACK/'neighborhoods.csv')
    summary=pd.read_csv(TRACK/'front_summary.csv').set_index('replicate')
    for r in REP_COLORS:
        p=cache[r+'_assignments']; x=cache[r+'_travel_matrix']; e=cache['expression_matrix']
        assert np.array_equal(np.sort(p,axis=1),np.broadcast_to(rows,p.shape))
        for q in REP_COLORS:
            df=transfer[(transfer.source==r)&(transfer.target==q)]
            np.testing.assert_allclose(cache[q+'_travel_matrix'][rows,p].sum(axis=1),df.travel,atol=1e-9)
            np.testing.assert_allclose(e[rows,p].sum(axis=1),df.expression,atol=1e-9)
    common=[]
    parent_index={p:i for i,p in enumerate(cache['natural_parents'])}
    for folder in [COMMON,COMMON/'expression_tradeoff']:
        cp=json.loads((folder/'provenance.json').read_text())
        assert hashlib.sha256(Path(cp['cache']).read_bytes()).hexdigest()==cp['cache_sha256']
        data=pd.read_csv(folder/'common_assignment_summary.csv')
        edges=pd.read_csv(folder/'common_assignments.csv')
        for record in data.itertuples():
            assert record.solver_status==0 and record.mip_gap==0
            a=edges[edges.expression_reduction_required_pct==record.expression_reduction_required_pct]
            assert np.array_equal(a.child.to_numpy(),cache['children'])
            assert a.common_parent.value_counts().to_dict()==a.natural_parent.value_counts().to_dict()
            pr=np.array([parent_index[p] for p in a.common_parent])
            np.testing.assert_allclose(cache['expression_matrix'][pr,rows].sum(),record.expression,atol=1e-9)
            for r in REP_COLORS:
                np.testing.assert_allclose(cache[r+'_travel_matrix'][pr,rows].sum(),getattr(record,r+'_travel'),atol=1e-9)
        common.append(data)
    common=pd.concat(common).sort_values('expression_reduction_required_pct').reset_index(drop=True)
    assert len(common)==3
    endpoints=transfer[(transfer.source!=transfer.target)&(transfer.alpha==1)].copy()
    endpoints['penalty_fraction_available']=[
        (row.travel-summary.loc[row.target,'travel_optimum'])/
        (summary.loc[row.target,'travel_natural']-summary.loc[row.target,'travel_optimum'])
        for row in endpoints.itertuples()]
    endpoints=endpoints.sort_values(['source','target'])
    return transfer,neighborhoods,common,endpoints


def style_axis(ax,letter,title):
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(labelsize=8,length=3)
    ax.set_title(title,loc='left',fontsize=10,pad=11)
    ax.text(-.18,1.075,letter,transform=ax.transAxes,fontsize=12,fontweight='bold')
    ax.grid(axis='y',color='#E8E8E8',linewidth=.5,zorder=0)
    ax.set_axisbelow(True)


def render_common_fronts(common):
    """Display every common solution in each embryo's own measured cost space."""
    base=ROOT/'terminal_pareto/output'
    fronts=pd.read_csv(base/'ce_tracking_replicate_fronts.csv')
    nulls=pd.read_csv(base/'ce_tracking_replicate_null_clouds.csv')
    audit=pd.read_csv(base/'ce_tracking_replicate_audit.csv')
    scales=pd.read_csv(TRACK/'front_summary.csv').set_index('replicate')
    # Confirm the displayed accepted fronts are on exactly the same cost scale.
    transfer=pd.read_csv(TRACK/'transfer.csv')
    for r in REP_COLORS:
        native=transfer[(transfer.source==r)&(transfer.target==r)]
        reference=fronts[fronts.replicate==r].sort_values('sweep_index')
        np.testing.assert_allclose((native.travel-scales.loc[r,'travel_natural'])/
            scales.loc[r,'travel_null_sd'],reference.travel_sigma,atol=1e-10)
        np.testing.assert_allclose((native.expression-common.expression_natural.iloc[0])/
            scales.loc[r,'expression_null_sd'],reference.cell_state_sigma,atol=1e-10)
    fig=plot_fronts(fronts,nulls,audit)
    fig.set_size_inches(8.2,5.55)
    ax=fig.axes[0]
    ax.set_position([.095,.52,.69,.38])
    ax.set_title('Matched fronts and jointly optimized assignments (275 edges)',loc='left',fontsize=10,pad=7)
    ax.set(xlabel='Travel change (first-cousin-null SD)',ylabel='Cell-state change (first-cousin-null SD)')
    ax.xaxis.label.set_size(8);ax.yaxis.label.set_size(8)
    ax.tick_params(labelsize=7)
    ax.grid(color='#E4E4E4',linewidth=.5)
    for legend in [child for child in ax.get_children() if isinstance(child,Legend)]:
        legend.set_frame_on(False)
        for text in legend.get_texts():text.set_fontsize(7)
        legend.get_title().set_fontsize(7.3)
    ax.get_legend().set_bbox_to_anchor((1.015,.43))
    markers=['*','s','^']; sizes=[65,24,30]; joint_color='#773377'
    handles=[Line2D([],[],marker=marker,color=joint_color,ls='',markersize=8 if i==0 else 5,
                    label=f'{budget:g}% required')
             for i,(marker,budget) in enumerate(zip(markers,common.expression_reduction_required_pct))]
    fig.legend(handles=handles,title='Common assignments: required cell-state reduction',
               loc='upper center',bbox_to_anchor=(.52,.438),ncol=3,frameon=False,
               fontsize=7.5,title_fontsize=8,handlelength=1.2,columnspacing=2)
    coordinates=[]
    for i,r in enumerate(REP_COLORS):
        mini=fig.add_axes([.095+i*.30,.10,.235,.215])
        f=fronts[fronts.replicate==r].sort_values('sweep_index')
        mini.axhspan(-2,0,xmin=0,xmax=3.15/3.45,color='#E7F1E8',zorder=0)
        mini.plot(f.travel_sigma,f.cell_state_sigma,color=REPLICATE_STYLES[r]['color'],
                  lw=1.2,marker='.',markersize=2,ls=REPLICATE_STYLES[r]['ls'])
        mini.axhline(0,color='#999999',lw=.6,ls=':')
        mini.axvline(0,color='#999999',lw=.6,ls=':')
        for j,row in enumerate(common.itertuples()):
            x=(getattr(row,r+'_travel')-getattr(row,r+'_natural_travel'))/scales.loc[r,'travel_null_sd']
            y=(row.expression-row.expression_natural)/scales.loc[r,'expression_null_sd']
            assert x<0 and y<0
            mini.scatter([x],[y],marker=markers[j],s=sizes[j],color=joint_color,
                         edgecolors='white',linewidth=.4,zorder=5)
            coordinates.append(dict(replicate=r,required_state_reduction_pct=row.expression_reduction_required_pct,
                                    travel_change_sigma=x,state_change_sigma=y))
        mini.scatter([0],[0],marker='X',s=42,color='#222222',edgecolor='white',linewidth=.5,zorder=6)
        gain=common.iloc[0][r+'_travel_reduction_pct']
        mini.set_title(f'Embryo {i+1}: {gain:.2f}% travel gain (star)',fontsize=7.3,pad=5)
        mini.set(xlim=(-3.15,.30),ylim=(-2,.4),xticks=[-3,-2,-1,0],yticks=[-2,-1,0],
                 xlabel='Travel change (null SD)')
        mini.xaxis.label.set_size(7)
        if i==0:mini.set_ylabel('Cell-state change (null SD)',fontsize=7)
        mini.tick_params(labelsize=7,length=2)
        mini.spines[['top','right']].set_visible(False)
    fig.savefig(OUT/'figS3A_ce_tracking_common_fronts.pdf')
    fig.savefig(OUT/'figS3A_ce_tracking_common_fronts.png',dpi=180)
    plt.close(fig)
    pd.DataFrame(coordinates).to_csv(OUT/'common_assignments_in_front_space.csv',index=False)


def render():
    OUT.mkdir(parents=True,exist_ok=True)
    transfer,neighborhoods,common,endpoints=load_and_validate()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,
        'axes.labelsize':9,'axes.linewidth':.7,'pdf.fonttype':42,'ps.fonttype':42})
    render_common_fronts(common)
    fig=plt.figure(figsize=(8.2,7.5))
    gs=fig.add_gridspec(3,2,left=.095,right=.98,top=.91,bottom=.07,hspace=.83,wspace=.34,
                       height_ratios=[1,1,1.04])
    ax,bx,cx,dx=[fig.add_subplot(gs[i,j]) for i,j in [(0,0),(0,1),(1,0),(1,1)]]
    bottom=gs[2,:].subgridspec(1,2,width_ratios=[1.25,1],wspace=.23)
    ex,tx=fig.add_subplot(bottom[0,0]),fig.add_subplot(bottom[0,1])
    fig.suptitle('Tracking geometry changes assignments; common improvements remain',fontsize=12.5,y=.987)
    fig.text(.5,.952,'C. elegans  |  3 tracking replicates  |  275 matched terminal edges  |  fixed protein costs',
             ha='center',fontsize=9,color='#555555')

    style_axis(ax,'C','Agreement in selected parents')
    style_axis(bx,'D','Overlap of candidate neighborhoods')
    for (r,q),color in PAIR_COLORS.items():
        data=transfer[(transfer.source==r)&(transfer.target==q)]
        ax.plot(data.alpha,100*data.parent_agreement,color=color,lw=1.4,label=f'{r[-1]} vs {q[-1]}')
        n=neighborhoods[(neighborhoods.source==r)&(neighborhoods.target==q)].groupby('k').overlap.mean()
        bx.plot(n.index,100*n.values,color=color,lw=1.4,marker='o',markersize=3.5)
    reference=neighborhoods.groupby('k').chance_overlap.mean()
    bx.plot(reference.index,100*reference.values,color='#999999',ls=':',lw=1.3)
    bx.text(10,17,'Uniform-set reference',fontsize=7.5,color='#777777')
    ax.set(xlabel='Travel weight',ylabel='Parent agreement (%)',xlim=(0,1),ylim=(40,102))
    ax.legend(frameon=False,fontsize=7.5,loc='lower left',ncol=3,handlelength=1.4,columnspacing=1)
    bx.set(xlabel='Nearest distinct candidate parents, k',ylabel='Mean overlap (%)',ylim=(0,102),xticks=[3,5,10,20])

    style_axis(cx,'E','Cost of transfer at the same weight')
    for r,q in itertools.permutations(REP_COLORS,2):
        data=transfer[(transfer.source==r)&(transfer.target==q)]
        pair=tuple(sorted([r,q])); color=PAIR_COLORS[pair]
        cx.plot(data.alpha,data.weighted_regret_pct_optimum,color=color,
                ls='-' if r<q else '--',lw=1.15,label=f'{r[-1]} to {q[-1]}')
    cx.set(xlabel='Travel weight',ylabel='Excess weighted cost (%)',xlim=(0,1),ylim=(-.4,25))
    cx.legend(frameon=False,ncol=2,fontsize=7.2,loc='upper left',handlelength=1.7,columnspacing=1)

    style_axis(dx,'F','Travel-endpoint transfer loss')
    labels=[]
    for i,row in enumerate(endpoints.itertuples()):
        color=PAIR_COLORS[tuple(sorted([row.source,row.target]))]
        dx.plot([0,100*row.penalty_fraction_available],[i,i],color=color,lw=1,alpha=.35)
        dx.scatter(100*row.penalty_fraction_available,i,color=color,s=24,zorder=3,
                   marker='o' if row.source<row.target else 's')
        labels.append(f'{row.source[-1]} to {row.target[-1]}')
    dx.axvline(100,color='#555555',lw=.8,ls='--')
    dx.text(102,-.75,'100%: saving erased',fontsize=7.2,color='#555555')
    dx.set(xlabel='Penalty / available travel saving (%)',yticks=range(6),yticklabels=labels,xlim=(0,182),ylim=(5.6,-1.1))
    dx.grid(False)

    style_axis(ex,'G','One assignment shared across all embryos')
    for r,marker in zip(REP_COLORS,['o','s','^']):
        ex.plot(common.expression_reduction_pct,common[r+'_travel_reduction_pct'],
                color=REP_COLORS[r],marker=marker,markersize=4,lw=.9,ls='--',label=f'Embryo {r[-1]}')
    ex.scatter([0],[0],marker='+',s=60,color='#222222',zorder=5)
    ex.text(.04,.11,'Natural',fontsize=8,color='#333333')
    ex.axhline(0,color='#AAAAAA',lw=.6)
    ex.set(xlabel='Cell-state cost reduction (%)',ylabel='Travel cost reduction (%)',xlim=(-.04,1.12),ylim=(-.13,2.75))
    ex.legend(frameon=False,fontsize=7.2,loc='upper right',ncol=3,handlelength=1,columnspacing=.6,
              bbox_to_anchor=(1.02,1.03))
    tx.axis('off')
    tx.text(0,1.0,'Guaranteed common improvement',fontsize=8.7,fontweight='bold',va='top')
    contents=[]
    for row in common.itertuples():
        contents.append([f'{row.expression_reduction_required_pct:g}%',
                         f'{row.guaranteed_travel_reduction_pct:.2f}%',
                         f'{row.natural_edge_retention*100:.1f}%'])
    table=tx.table(cellText=contents,colLabels=['Required\nstate gain','Worst travel\ngain','Natural edges\nretained'],
                   cellLoc='center',colLoc='center',bbox=[0,.27,1,.55],colWidths=[.29,.32,.39])
    table.auto_set_font_size(False);table.set_fontsize(8)
    for (r,c),cell in table.get_celld().items():
        cell.set_linewidth(.5);cell.set_edgecolor('#DDDDDD')
        if r==0:cell.set_facecolor('#F2F4F6');cell.set_text_props(weight='bold',fontsize=7.5)
    tx.text(0,.13,'One feasible common assignment per row.\nZero reported solver gap for all three.',
            fontsize=7.5,color='#555555',va='top',linespacing=1.4)
    fig.savefig(OUT/'figS3C_G_ce_tracking_variance.pdf')
    fig.savefig(OUT/'figS3C_G_ce_tracking_variance.png',dpi=180)
    plt.close(fig)
    endpoints.to_csv(OUT/'travel_endpoint_transfer_summary.csv',index=False)
    common.to_csv(OUT/'common_assignment_summary.csv',index=False)
    # Add common-assignment views only to the candidate; accepted files stay intact.
    existing=ROOT/'terminal_pareto/output/publication/figS3_ce_tracking_robustness.tex'
    old_text=existing.read_text().replace('\\begin{document}',
        '\\renewcommand{\\thefigure}{S3}\n'
        + '\\graphicspath{{'+str(OUT)+'/}{'+str(existing.parent)+'/}}\n\\begin{document}')
    old_text=old_text.replace('{figS3A_ce_tracking_replicate_fronts.pdf}',
                              '{figS3A_ce_tracking_common_fronts.pdf}')
    old_text=old_text.replace('The three fronts remain closely aligned under the matched analysis.',
        'The three fronts remain closely aligned under the matched analysis. '
        'Below, three zooms place the common assignments from panel G in each '
        "embryo's own cost space. A symbol denotes the same assignment in all "
        'three zooms: star, square, and triangle require 0\\%, 0.5\\%, and '
        '1\\% cell-state reduction, respectively. All lie in the shaded '
        'lower-left region, with both costs below nature. The starred '
        'assignment saves 2.41\\%, 2.32\\%, and 2.28\\% travel. Common '
        'assignments need not lie on any individual-embryo front.')
    (OUT/'figS3_existing_numbered.tex').write_text(old_text)
    sources=[Path(__file__),ROOT/'terminal_pareto/figS3_ce_tracking_variance_amendment.tex',
             TRACK/'assignments_and_costs.npz',TRACK/'transfer.csv',TRACK/'neighborhoods.csv',
             COMMON/'common_assignment_summary.csv',COMMON/'expression_tradeoff/common_assignment_summary.csv',
             COMMON/'common_assignments.csv',COMMON/'expression_tradeoff/common_assignments.csv',
             ROOT/'terminal_pareto/output/publication/figS3_ce_tracking_robustness.pdf',existing,
             ROOT/'terminal_pareto/figS3_ce_tracking_robustness.py',
             ROOT/'terminal_pareto/output/ce_tracking_replicate_fronts.csv',
             ROOT/'terminal_pareto/output/ce_tracking_replicate_null_clouds.csv',
             ROOT/'terminal_pareto/output/ce_tracking_replicate_audit.csv']
    (OUT/'figure_provenance.json').write_text(json.dumps(dict(
        sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        verified='Source hashes, all transferred costs, common assignment costs/capacities, solver statuses',
        scope='CE only; candidate panel A adds common-assignment zooms; accepted S3 assets unchanged'),indent=2)+'\n')
    print('Rendered and validated:',OUT,flush=True)


if __name__=='__main__':
    render()
