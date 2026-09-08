"""tait — embedding-preserving Tait graphs for knot diagrams (05_unknotting_number project).

Modules
-------
knotinfo    KnotInfo table access (database_knotinfo CSV) and Jones-string parsing
graph       EmbeddedTaitGraph: signed planar map (rotation system), Link <-> graph, PD, JSON/GraphML
moves       Reidemeister moves R1/R2/R3, flype, crossing change on the embedded graph (+ enumerators)
canonical   canonical code of a signed planar map (search de-duplication / transposition tables)
invariants  Goeritz matrix, det, H1(Sigma_2), linking form (Lickorish), Casson–Walker, Regina Jones,
            SnapPy double branched cover, HFK/unknot detection, identification

Run inside the pyenv virtualenv `unknot-venv` (snappy, spherogram, regina, knot_floer_homology,
database_knotinfo).  No Sage needed.
"""
from . import _compat  # noqa: F401  (sqlite3 shim for Pythons built without _sqlite3)
from .graph import EmbeddedTaitGraph, from_link, from_pd, to_pd, to_link, ribbon_faces  # noqa: F401
from .canonical import canonical_code  # noqa: F401
from . import moves, invariants, knotinfo  # noqa: F401

__version__ = '0.1.0'
