#!/usr/bin/env python
"""Figures for the Tait-graph part of the Background: (1) a knot diagram and its Tait graph,
(2) Reidemeister moves as Tait-graph operations, (3) the flype."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tait  # noqa
import numpy as np
import matplotlib
os.makedirs('figures', exist_ok=True)
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

# ---------------------------------------------------------------- (1) example: 6_2 diagram + Tait graph
def fig_example():
    import json
    from spherogram import Link
    from spherogram.links.orthogonal import OrthogonalLinkDiagram
    from tait import knotinfo, graph as tg
    nm = '6_2'
    pd = knotinfo.pd_code(nm)
    L = Link([list(q) for q in pd])
    O = OrthogonalLinkDiagram(L)
    verts, arrows, crossings = O.plink_data()
    verts = [np.array(v, float) for v in verts]
    segs = [(verts[a], verts[b]) for a, b in arrows]
    unit = min(np.linalg.norm(p2-p1) for p1, p2 in segs if np.linalg.norm(p2-p1) > 1e-9)
    gap = 0.30*unit
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.6, 3.6))
    # left: diagram
    cuts = {i: [] for i in range(len(segs))}
    for (iu, io, _, lab) in crossings:
        p1, p2 = segs[iu]; q1, q2 = segs[io]
        d1 = p2-p1; d2 = q2-q1
        den = d1[0]*d2[1]-d1[1]*d2[0]
        t = ((q1-p1)[0]*d2[1]-(q1-p1)[1]*d2[0])/den
        cuts[iu].append(p1+t*d1)
    for i, (p1, p2) in enumerate(segs):
        d = p2-p1; Ls = np.linalg.norm(d)
        if Ls < 1e-9: continue
        u = d/Ls
        cur = 0.0; pieces = []
        for c in sorted(float(np.dot(P-p1, u)) for P in cuts[i]):
            pieces.append((cur, max(cur, c-gap))); cur = min(Ls, c+gap)
        pieces.append((cur, Ls))
        for a, b in pieces:
            if b-a > 1e-6:
                q1, q2 = p1+a*u, p1+b*u
                ax1.plot([q1[0], q2[0]], [-q1[1], -q2[1]], '-', color='black', lw=2.0, solid_capstyle='round')
    ax1.set_title('a diagram of $6_2$', fontsize=10)
    # right: Tait graph from our machinery, spring layout
    g = tg.from_link(Link([list(q) for q in pd]))
    import networkx as nx
    G = nx.MultiGraph()
    for v in g.vertices: G.add_node(v)
    for e,(u,w) in g.edges.items(): G.add_edge(u, w, key=e, sign=g.signs[e])
    pos = nx.spring_layout(G, seed=4, iterations=300)
    for (u, w, k, dat) in G.edges(keys=True, data=True):
        # separate parallel edges by curvature
        par = [kk for (_,_,kk) in G.edges(u, keys=True) if G.has_edge(u, w, kk)]
        idx = sorted(kk for (uu,ww,kk) in G.edges(keys=True) if {uu,ww}=={u,w}).index(k)
        tot = len([1 for (uu,ww,kk) in G.edges(keys=True) if {uu,ww}=={u,w}])
        rad = 0.0 if tot == 1 else (-0.25 + 0.5*idx/(tot-1))
        style = '-' if dat['sign'] > 0 else (0, (4, 2))
        ar = FancyArrowPatch(pos[u], pos[w], connectionstyle=f'arc3,rad={rad}', arrowstyle='-',
                             lw=1.8, color='black', linestyle=style)
        ax2.add_patch(ar)
    for v, p in pos.items():
        ax2.plot(*p, 'o', ms=9, color='black', zorder=3)
    ax2.set_title('its Tait graph (all signs equal: alternating)', fontsize=10)
    for ax in (ax1, ax2):
        ax.set_aspect('equal'); ax.axis('off')
    fig.tight_layout(); fig.savefig('figures/fig_tait_example.pdf', bbox_inches='tight')
    fig.savefig('figures/fig_tait_example.png', dpi=140, bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- helpers for schematic panels
def node(ax, x, y):
    ax.plot(x, y, 'o', ms=8, color='black', zorder=3)

def edge(ax, p, q, sign='+', rad=0.0, label=None, lpos=None):
    style = '-' if sign == '+' else (0, (4, 2))
    ax.add_patch(FancyArrowPatch(p, q, connectionstyle=f'arc3,rad={rad}', arrowstyle='-',
                                 lw=1.8, color='black', linestyle=style))
    if label:
        m = lpos if lpos is not None else ((p[0]+q[0])/2, (p[1]+q[1])/2+0.16+0.35*rad)
        ax.text(*m, label, fontsize=9, ha='center', va='center')

def arrow_lr(ax, x, y):
    ax.annotate('', xy=(x+0.35, y), xytext=(x-0.35, y),
                arrowprops=dict(arrowstyle='<->', lw=1.2))

# ---------------------------------------------------------------- (2) R-moves on the Tait graph
def fig_moves():
    fig, axes = plt.subplots(1, 3, figsize=(8.6, 2.6))
    # R1: loop / pendant edge
    ax = axes[0]
    node(ax, 0.0, 0.5)
    th = np.linspace(0, 2*np.pi, 100)
    ax.plot(0.28*np.cos(th)+0.28, 0.28*np.sin(th)+0.5, 'k-', lw=1.8)
    ax.text(0.62, 0.88, '$\\pm$', fontsize=10)
    node(ax, 1.5, 0.75); node(ax, 1.5, 0.25)
    edge(ax, (1.5, 0.75), (1.5, 0.25), '+', label='$\\pm$', lpos=(1.72, 0.5))
    arrow_lr(ax, 2.35, 0.5)
    node(ax, 3.1, 0.5)
    ax.set_title('R1: loop / pendant edge', fontsize=10)
    # R2: parallel pair and series pair
    ax = axes[1]
    node(ax, 0.0, 0.5); node(ax, 1.0, 0.5)
    edge(ax, (0, 0.5), (1, 0.5), '+', rad=0.45, label='$+$', lpos=(0.5, 0.98))
    edge(ax, (0, 0.5), (1, 0.5), '-', rad=-0.45, label='$-$', lpos=(0.5, 0.03))
    arrow_lr(ax, 1.6, 0.5)
    node(ax, 2.2, 0.5); node(ax, 3.0, 0.5)
    ax.text(2.6, 0.28, '(parallel)', fontsize=8, ha='center')
    node(ax, 0.0, -0.45); node(ax, 0.75, -0.45); node(ax, 1.5, -0.45)
    edge(ax, (0, -0.45), (0.75, -0.45), '+', label='$+$', lpos=(0.37, -0.28))
    edge(ax, (0.75, -0.45), (1.5, -0.45), '-', label='$-$', lpos=(1.12, -0.28))
    arrow_lr(ax, 2.0, -0.45)
    node(ax, 2.7, -0.45)
    ax.text(2.7, -0.75, '(series: ends merge)', fontsize=8, ha='center')
    ax.set_title('R2: parallel / series pair', fontsize=10)
    # R3: Y -- Delta
    ax = axes[2]
    c = (0.5, 0.42)
    outer = [(0.0, 0.0), (1.0, 0.0), (0.5, 1.05)]
    for p in outer: node(ax, *p)
    node(ax, *c)
    for p, sl, lp in zip(outer, ['$s_1$', '$s_2$', '$s_3$'],
                         [(0.06, 0.34), (0.95, 0.34), (0.68, 0.78)]):
        edge(ax, c, p, '+', label=sl, lpos=lp)
    arrow_lr(ax, 1.62, 0.45)
    sh = 2.35
    outer2 = [(sh, 0.0), (sh+1.0, 0.0), (sh+0.5, 1.05)]
    for p in outer2: node(ax, *p)
    edge(ax, outer2[0], outer2[1], '-', label='$-s_3$', lpos=(sh+0.5, -0.2))
    edge(ax, outer2[1], outer2[2], '-', label='$-s_1$', lpos=(sh+1.03, 0.6))
    edge(ax, outer2[0], outer2[2], '-', label='$-s_2$', lpos=(sh-0.05, 0.6))
    ax.set_title('R3: $Y$--$\\Delta$, signs negated', fontsize=10)
    for ax in axes:
        ax.set_aspect('equal'); ax.axis('off')
        ax.set_xlim(-0.6, 3.7); ax.set_ylim(-1.0, 1.35)
    fig.tight_layout(); fig.savefig('figures/fig_tait_moves.pdf', bbox_inches='tight')
    fig.savefig('figures/fig_tait_moves.png', dpi=140, bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- (3) flype
def fig_flype():
    fig, ax = plt.subplots(figsize=(7.2, 2.3))
    def tangle(x, y, lab):
        ax.add_patch(Rectangle((x, y), 1.0, 0.9, fill=False, lw=1.6))
        ax.text(x+0.5, y+0.45, lab, fontsize=12, ha='center', va='center')
    # left: crossing then tangle
    y0 = 0.0
    ax.plot([-1.1, -0.45], [y0+0.7, y0+0.2], 'k-', lw=1.8)
    ax.plot([-1.1, -0.82], [y0+0.2, y0+0.415], 'k-', lw=1.8)
    ax.plot([-0.73, -0.45], [y0+0.485, y0+0.7], 'k-', lw=1.8)
    ax.plot([-0.45, -0.2], [y0+0.7, y0+0.7], 'k-', lw=1.8)
    ax.plot([-0.45, -0.2], [y0+0.2, y0+0.2], 'k-', lw=1.8)
    tangle(-0.2, y0-0.02, '$T$')
    ax.plot([0.8, 1.15], [y0+0.7, y0+0.7], 'k-', lw=1.8)
    ax.plot([0.8, 1.15], [y0+0.2, y0+0.2], 'k-', lw=1.8)
    ax.annotate('', xy=(2.15, y0+0.45), xytext=(1.45, y0+0.45), arrowprops=dict(arrowstyle='<->', lw=1.2))
    # right: tangle rotated then crossing
    x2 = 2.45
    ax.plot([x2, x2+0.35], [y0+0.7, y0+0.7], 'k-', lw=1.8)
    ax.plot([x2, x2+0.35], [y0+0.2, y0+0.2], 'k-', lw=1.8)
    tangle(x2+0.35, y0-0.02, '$T^{\\updownarrow}$')
    xr = x2+1.35
    ax.plot([xr, xr+0.65], [y0+0.7, y0+0.2], 'k-', lw=1.8)
    ax.plot([xr, xr+0.28], [y0+0.2, y0+0.415], 'k-', lw=1.8)
    ax.plot([xr+0.37, xr+0.65], [y0+0.485, y0+0.7], 'k-', lw=1.8)
    ax.set_xlim(-1.35, 4.6); ax.set_ylim(-0.35, 1.15)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_title('the flype: on the Tait graph, a block at a 2-separation is inverted', fontsize=10)
    fig.tight_layout(); fig.savefig('figures/fig_flype.pdf', bbox_inches='tight')
    fig.savefig('figures/fig_flype.png', dpi=140, bbox_inches='tight'); plt.close(fig)

fig_example(); fig_moves(); fig_flype(); print('tait figures done')
