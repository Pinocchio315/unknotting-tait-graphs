#!/usr/bin/env python
"""Draw a knot diagram from a PD code as a clean orthogonal (plink-style) figure, circling given crossings.
Uses spherogram's OrthogonalLinkDiagram (the engine behind PLink's smooth export); pure matplotlib output."""
import json, sys
import os
from pathlib import Path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tait_graphs'))
import tait  # noqa (sqlite shim)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from spherogram import Link
from spherogram.links.orthogonal import OrthogonalLinkDiagram


def seg_intersection(p1, p2, q1, q2):
    p1, p2, q1, q2 = map(lambda t: np.array(t, float), (p1, p2, q1, q2))
    d1, d2 = p2 - p1, q2 - q1
    den = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(den) < 1e-12:
        return None
    t = ((q1 - p1)[0] * d2[1] - (q1 - p1)[1] * d2[0]) / den
    s = ((q1 - p1)[0] * d1[1] - (q1 - p1)[1] * d1[0]) / den
    if -1e-9 <= t <= 1 + 1e-9 and -1e-9 <= s <= 1 + 1e-9:
        return p1 + t * d1
    return None


def draw(pd, marked_rows, out_pdf, title=None):
    out_pdf = str(out_pdf)
    Path(out_pdf).parent.mkdir(parents=True, exist_ok=True)
    # Rows refer to the supplied PD, not to crossings in a relabelled diagram.
    # OrthogonalLinkDiagram carries these labels into plink_data; reject a
    # missing mark instead of silently publishing an incomplete certificate.
    if len(set(marked_rows)) != len(marked_rows) or any(i < 0 or i >= len(pd) for i in marked_rows):
        raise ValueError('Invalid marked PD rows')
    L = Link([list(q) for q in pd])
    O = OrthogonalLinkDiagram(L)
    verts, arrows, crossings = O.plink_data()
    verts = [np.array(v, float) for v in verts]
    segs = [(verts[a], verts[b]) for a, b in arrows]
    unit = min(np.linalg.norm(p2 - p1) for p1, p2 in segs if np.linalg.norm(p2 - p1) > 1e-9)
    gap = 0.30 * unit
    cuts = {i: [] for i in range(len(segs))}
    marks = []
    for (iu, io, _, label) in crossings:
        P = seg_intersection(*segs[iu], *segs[io])
        if P is None:
            continue
        cuts[iu].append(P)
        if label in marked_rows:
            marks.append(P)
    if len(marks) != len(marked_rows):
        raise ValueError('The drawing did not preserve every marked crossing')
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    for i, (p1, p2) in enumerate(segs):
        d = p2 - p1
        Lseg = np.linalg.norm(d)
        if Lseg < 1e-9:
            continue
        u = d / Lseg
        ts = sorted([0.0] + [float(np.dot(P - p1, u)) for P in cuts[i]] + [Lseg])
        pieces = []
        cur = 0.0
        cutpts = sorted(float(np.dot(P - p1, u)) for P in cuts[i])
        for c in cutpts:
            pieces.append((cur, max(cur, c - gap)))
            cur = min(Lseg, c + gap)
        pieces.append((cur, Lseg))
        for a, b in pieces:
            if b - a > 1e-6:
                q1, q2 = p1 + a * u, p1 + b * u
                ax.plot([q1[0], q2[0]], [-q1[1], -q2[1]], '-', color='black', lw=2.0,
                        solid_capstyle='round', zorder=2)
    for P in marks:
        ax.add_patch(plt.Circle((P[0], -P[1]), 0.55 * unit, fill=False, color='red', lw=1.7, zorder=3))
    if title:
        ax.set_title(title, fontsize=11)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches='tight')
    fig.savefig(out_pdf.replace('.pdf', '.png'), dpi=140, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    import argparse
    from verify_presentation_diagrams import load_paper_diagrams, load_additional_diagrams
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('names', nargs='*', help='KnotInfo names; default: all nineteen v1.8 diagrams')
    parser.add_argument('--outdir', type=Path, default=Path('figures'))
    args = parser.parse_args()
    # Both generations use this one renderer; data, marks, and endpoint claims
    # are loaded by the same structural checks as the presentation verifier.
    diagrams = {n: (pd, marks) for n, pd, marks in load_paper_diagrams()}
    diagrams.update({n: (pd, marks) for n, pd, marks, _ in load_additional_diagrams()})
    for nm in args.names or diagrams:
        pd, marks = diagrams[nm]
        draw(pd, marks, args.outdir / f'fig_{nm}.pdf',
             title=nm.replace('_', ''))
        print('drawn', nm, len(pd), 'crossings, marked', marks)
