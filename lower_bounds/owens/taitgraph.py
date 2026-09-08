"""taitgraph — shim for the shared embedded-Tait-graph module.

The self-contained Slurm bundle shipped a verbatim copy of ``tait/graph.py`` under this
name so that the obstruction package had no external project dependencies.  In this
repository the shared package ``tait_graphs/tait`` is used instead; this module simply
re-exports it (the file it replaces was byte-identical to ``tait/graph.py``).
"""
import os as _os
import sys as _sys

_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                  _os.pardir, _os.pardir, 'tait_graphs'))
import tait.graph as _graph  # noqa: E402

_sys.modules[__name__] = _graph
