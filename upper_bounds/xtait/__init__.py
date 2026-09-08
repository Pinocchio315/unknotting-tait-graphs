"""xtait — self-contained (stdlib-only) embedded-Tait-graph toolkit for the crossing-additivity campaign.

This package is a frozen, dependency-free snapshot of the project's diagram machinery plus the P0
pass-move module.  It deliberately does NOT import spherogram / snappy / regina / database_knotinfo /
networkx / sympy, so the whole `crossing_additivity/` directory can be copied to any server with a plain
Python >= 3.9 and run as-is.

Modules:
  graph      EmbeddedTaitGraph (signed planar map with rotation system), PD <-> graph, dual, splice sums
  moves      exact isotopy moves R1±, R2± (parallel/series), R3 (Y–Δ), flype  (+ crossing change)
  canonical  canonical code of a diagram (label/colour-class independent) for de-duplication
  jones      exact Kauffman bracket / Jones polynomial from a PD (state sum; for verification only)
  passmove   P0: pass move (pick up an all-over/all-under strand and reroute it) as a recorded move
  reduce     pass-first crossing-minimisation (greedy + best-first search with budgets)
  verify     independent certificate replayer

Conventions are identical to the main project's `tait` package (documented in graph.py); everything is
re-verified internally by `selftest.py` using the Kauffman bracket as an independent oracle.
"""
__version__ = '0.1.0'
