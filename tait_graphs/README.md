# Embedded Tait graphs

The shared `tait/` package supports the diagrammatic constructions in Sections 2–4 of *Computation of Unknotting Numbers: Which Knot Breaks the Bernhard–Jablan Conjecture?*, manuscript v1.7. A diagram is represented by a signed plane multigraph with a rotation system. Changing one crossing reverses one edge sign; isotopies must preserve the embedded graph's diagrammatic information. See the [repository README](../README.md) for installation and the frozen `v1.3` computational profile used by the current manuscript.

| Module | Purpose |
|---|---|
| `graph.py` | Link, PD-code, and embedded-graph conversions, including duals |
| `moves.py` | Reidemeister moves, flypes, and crossing changes |
| `canonical.py` | Canonical codes of signed planar maps for diagram comparison |
| `invariants.py` | Goeritz forms, determinant, homology and linking pairings, Jones polynomial, knot Floer genus, and double branched covers |
| `knotinfo.py` | Access to KnotInfo diagrams and Jones-polynomial parsing |

The production upper-bound search and its recorded isotopies live in `../upper_bounds/`. In `invariants.py`, the cyclic linking-pairing generator is computed by the shared exact helper in `../lower_bounds/linking_pairing.py`. Double branched covers are constructed using the meridian of the knot exterior, rather than selecting a Dehn filling because it has the expected homology order. The homology comparison is only a consistency check. Linking-pairing signs and Casson–Walker signs require a common orientation convention before they can be combined.

This directory provides shared representations and invariants; running a helper
does not update the manuscript's tables. For the deposited crossing-change
formula audit, see [crossing_changes/README.md](../crossing_changes/README.md).
It compares the determinant, signature, dual-resistance, and self-linking
formulas with 69,632 recorded crossing changes from 5,546 parent knots.

From the repository root (`code/` in the manuscript workspace):

```python
import sys
sys.path.insert(0, 'tait_graphs')
from spherogram import Link
from tait import knotinfo, graph, invariants
knot = graph.from_link(Link(knotinfo.pd_code('13n_1587')))
print(invariants.determinant(knot))  # 75
print(invariants.h1_torsion(knot))  # [75]
```

Run `python -m unittest discover -s tests -p test_upper_bounds.py -v` for bounded exact-invariant, graph-move, and branched-cover regressions. Numerical SnapPy identification is used as a computational comparison; these checks do not constitute a formal verification of every stored search result.
