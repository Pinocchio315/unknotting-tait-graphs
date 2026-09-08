# Embedded Tait graphs

The shared `tait/` package supports the diagrammatic constructions in Sections 2–4 of *How to Compute the Unknotting Number: Theory and Computations*. A diagram is represented by a signed plane multigraph with a rotation system. Changing one crossing reverses one edge sign; isotopies must preserve the embedded graph's diagrammatic information.

| Module | Purpose |
|---|---|
| `graph.py` | Link, PD-code, and embedded-graph conversions, including duals |
| `moves.py` | Reidemeister moves, flypes, and crossing changes |
| `canonical.py` | Canonical codes of signed planar maps for diagram comparison |
| `invariants.py` | Goeritz forms, determinant, homology and linking pairings, Jones polynomial, knot Floer genus, and double branched covers |
| `knotinfo.py` | Access to KnotInfo diagrams and Jones-polynomial parsing |
| `search.py`, `polygon.py` | Earlier exploratory utilities, outside the paper verification pipeline |

The production upper-bound search and its recorded isotopies live in `../upper_bounds/`. In `invariants.py`, the cyclic linking-pairing generator is computed by the shared exact helper in `../lower_bounds/linking_pairing.py`. Double branched covers are constructed using the meridian of the knot exterior, rather than selecting a Dehn filling because it has the expected homology order. The homology comparison is only a consistency check. Linking-pairing signs and Casson–Walker signs require a common orientation convention before they can be combined.

From `code/`:

```python
import sys
sys.path.insert(0, 'tait_graphs')
from spherogram import Link
from tait import knotinfo, graph, invariants
knot = graph.from_link(Link(knotinfo.pd_code('13n_1587')))
print(invariants.determinant(knot))  # 75
print(invariants.h1_torsion(knot))  # [75]
```

Run `python -m unittest discover -s tests -p test_upper_bounds.py -v` for bounded exact-invariant, graph-move, and branched-cover regressions. Numerical SnapPy identification is used as a computational comparison; these checks do not constitute a formal verification of every stored search result. `make_tait_figs.py` and `figures/` retain preliminary figure material.
