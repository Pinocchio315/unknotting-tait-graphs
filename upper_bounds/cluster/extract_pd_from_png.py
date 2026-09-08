"""Extract a PD code from a clean knot-diagram PNG (blue strokes, gaps at under-crossings)."""
import sys, numpy as np
from PIL import Image
from scipy import ndimage as ndi

def load_mask(path, closing=1):
    a = np.asarray(Image.open(path).convert('RGB')).astype(int)
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    red = (r > 180) & (g < 120) & (b < 120)
    blue = (b > 120) & (b - r > 40) & (b - g > 40)
    blue &= ~red
    # bridge tiny gaps (red box outline cut the stroke): closing with a 5x5 disk
    if closing > 1:
        blue = ndi.binary_closing(blue, structure=np.ones((closing, closing)), iterations=1)
    return blue, red

def thin(img):
    """Zhang-Suen thinning."""
    img = img.copy().astype(np.uint8)
    def step(img, it):
        P = np.pad(img, 1)
        p2 = P[:-2, 1:-1]; p3 = P[:-2, 2:]; p4 = P[1:-1, 2:]; p5 = P[2:, 2:]
        p6 = P[2:, 1:-1]; p7 = P[2:, :-2]; p8 = P[1:-1, :-2]; p9 = P[:-2, :-2]
        nb = [p2, p3, p4, p5, p6, p7, p8, p9]
        B = sum(n.astype(int) for n in nb)
        A = sum(((nb[i] == 0) & (nb[(i + 1) % 8] == 1)).astype(int) for i in range(8))
        if it == 0:
            c1 = (p2 * p4 * p6 == 0); c2 = (p4 * p6 * p8 == 0)
        else:
            c1 = (p2 * p4 * p8 == 0); c2 = (p2 * p6 * p8 == 0)
        rem = (img == 1) & (B >= 2) & (B <= 6) & (A == 1) & c1 & c2
        img[rem] = 0
        return rem.any()
    while True:
        a = step(img, 0); b = step(img, 1)
        if not (a or b): break
    return img.astype(bool)

def neighbours(sk, y, x):
    H, W = sk.shape; out = []
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy == 0 and dx == 0: continue
            yy, xx = y + dy, x + dx
            if 0 <= yy < H and 0 <= xx < W and sk[yy, xx]: out.append((yy, xx))
    return out

def prune_spurs(sk, maxlen=8):
    sk = sk.copy()
    changed = True
    while changed:
        changed = False
        ys, xs = np.nonzero(sk)
        deg = {(y, x): len(neighbours(sk, y, x)) for y, x in zip(ys, xs)}
        for (y, x), d in deg.items():
            if d == 1:
                # walk until a junction
                path = [(y, x)]; prev = None; cur = (y, x)
                while True:
                    nb = [n for n in neighbours(sk, *cur) if n != prev and n not in path]
                    if len(nb) != 1: break
                    prev, cur = cur, nb[0]; path.append(cur)
                    if len(path) > maxlen: break
                if len(path) <= maxlen and len(neighbours(sk, *cur)) >= 3:
                    for p in path[:-1]: sk[p] = False
                    changed = True
                    break
    return sk

def trace_components(blue):
    """One polyline per connected blue stroke: thin the stroke, then take the longest
    shortest path (graph diameter) through its skeleton, which ignores spurs and tiny loops."""
    import collections
    lab, n = ndi.label(blue, structure=np.ones((3, 3)))
    comps = []
    for i in range(1, n + 1):
        mask = lab == i
        if mask.sum() < 3: continue
        sk = thin(mask)
        pts = list(zip(*np.nonzero(sk)))
        if len(pts) < 2:
            ys, xs = np.nonzero(mask); pts = list(zip(ys, xs))
        pset = set(pts)
        def bfs(src):
            dist = {src: 0}; prev = {src: None}; dq = collections.deque([src])
            while dq:
                c = dq.popleft()
                for q in neighbours(sk, *c):
                    if q in pset and q not in dist:
                        dist[q] = dist[c] + 1; prev[q] = c; dq.append(q)
            far = max(dist, key=dist.get)
            return far, dist, prev
        a, _, _ = bfs(pts[0])
        b, dist, prev = bfs(a)
        path = []; c = b
        while c is not None: path.append(c); c = prev[c]
        path.reverse()
        comps.append(np.array([(x, y) for y, x in path], dtype=float))
    return comps


def tangent(poly, at_end, k=6):
    if at_end == 0:
        a, b = poly[min(k, len(poly) - 1)], poly[0]
    else:
        a, b = poly[max(len(poly) - 1 - k, 0)], poly[-1]
    v = b - a; return v / (np.linalg.norm(v) + 1e-9)

def extract(path, verbose=True, closing=1):
    blue, red = load_mask(path, closing)
    comps = trace_components(blue)
    ends = []   # (comp, which_end, point, outward tangent)
    for ci, poly in enumerate(comps):
        ends.append((ci, 0, poly[0], tangent(poly, 0)))
        ends.append((ci, 1, poly[-1], tangent(poly, 1)))
    # pair endpoints across gaps
    used = set(); gaps = []
    for i, (ci, ei, p, t) in enumerate(ends):
        if i in used: continue
        best = None
        for j, (cj, ej, q, u) in enumerate(ends):
            if j == i or j in used: continue
            d = q - p; dist = np.linalg.norm(d)
            if dist > 48 or dist < 1: continue
            dn = d / dist
            score = float(np.dot(t, dn)) + float(np.dot(u, -dn))
            if np.dot(t, dn) > 0.6 and np.dot(u, -dn) > 0.6:
                if best is None or dist < best[0]: best = (dist, j)
        if best is None:
            print('unpaired endpoint', i, p, file=sys.stderr); continue
        j = best[1]; used |= {i, j}
        gaps.append((i, j))
    if verbose: print(f'{len(comps)} pieces, {len(gaps)} gaps', file=sys.stderr)
    # over strand for each gap: nearest polyline point to the gap midpoint
    crossings = []
    for gi, (i, j) in enumerate(gaps):
        p, q = ends[i][2], ends[j][2]; m = (p + q) / 2
        best = None
        for ci, poly in enumerate(comps):
            d = np.linalg.norm(poly - m, axis=1); k = int(np.argmin(d))
            if best is None or d[k] < best[0]: best = (d[k], ci, k)
        crossings.append({'gap': (i, j), 'mid': m, 'over': (best[1], best[2]), 'join': best[0] > 6})
    joins = [c for c in crossings if c['join']]
    if verbose: print(f'{len(joins)} gaps are plain joins (no over strand)', file=sys.stderr)
    # traversal
    end_index = {(ends[i][0], ends[i][1]): i for i in range(len(ends))}
    gap_of_end = {}
    for gi, (i, j) in enumerate(gaps): gap_of_end[i] = (gi, j); gap_of_end[j] = (gi, i)
    over_at = {}
    for gi, c in enumerate(crossings):
        if not c['join']: over_at.setdefault(c['over'][0], []).append((c['over'][1], gi))
    passages = []   # (crossing id, 'over'/'under', direction vector)
    ci, direction = 0, +1   # start at comp 0 from end 0 to end 1
    start = (ci, direction); visited_pieces = 0
    while True:
        poly = comps[ci]
        idxs = sorted(over_at.get(ci, []), key=lambda t: t[0], reverse=(direction < 0))
        for k, gi in idxs:
            a = poly[max(k - 6, 0)]; b = poly[min(k + 6, len(poly) - 1)]
            v = (b - a) * direction; v = v / (np.linalg.norm(v) + 1e-9)
            passages.append((gi, 'over', v))
        exit_end = 1 if direction > 0 else 0
        ei = end_index[(ci, exit_end)]
        gi, ej = gap_of_end[ei]
        p, q = ends[ei][2], ends[ej][2]; v = (q - p); v = v / np.linalg.norm(v)
        if not crossings[gi]['join']: passages.append((gi, 'under', v))
        cj, ej_which = ends[ej][0], ends[ej][1]
        ci, direction = cj, (+1 if ej_which == 0 else -1)
        visited_pieces += 1
        if (ci, direction) == start: break
        if visited_pieces > 10 * len(comps): raise RuntimeError('traversal did not close')
    real = [gi for gi, c in enumerate(crossings) if not c['join']]
    n = len(real); renum = {gi: k for k, gi in enumerate(real)}
    passages = [(renum[gi], kind, v) for gi, kind, v in passages]
    crossings = [crossings[gi] for gi in real]
    if len(passages) != 2 * n: print('WARNING passages', len(passages), 'crossings', n, file=sys.stderr)
    # PD code
    N = len(passages)
    def rot_ccw(v): return np.array([v[1], -v[0]])   # visual ccw rotation in image coordinates
    info = {}
    for idx, (gi, kind, v) in enumerate(passages):
        info.setdefault(gi, {})[kind] = (idx, v)
    pd = []
    for gi in range(n):
        u_idx, u_v = info[gi]['under']; o_idx, o_v = info[gi]['over']
        a_in, a_out = (u_idx) % N or N, (u_idx + 1) % N or N        # arcs: arc k from passage k to k+1 (1-based)
        o_in, o_out = (o_idx) % N or N, (o_idx + 1) % N or N
        # position of over arcs: incoming over arc lies at -o_v, outgoing at +o_v; j is the one at rot_ccw(-u_v)
        ref = rot_ccw(-u_v)
        if np.dot(ref, -o_v) > np.dot(ref, o_v): j, l = o_in, o_out
        else: j, l = o_out, o_in
        pd.append([a_in, j, a_out, l])
    return pd, crossings, comps, red

if __name__ == '__main__':
    pd, crossings, comps, red = extract(sys.argv[1])
    print(pd)
