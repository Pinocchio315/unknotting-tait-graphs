"""Compatibility shims applied on `import tait`.

sqlite3: spherogram -> snappy_manifolds needs the stdlib sqlite3 module.  Some self-built Pythons (e.g. on a
cluster) lack the `_sqlite3` extension; in that case the manylinux wheel `pysqlite3-binary` is used instead
(`pip install pysqlite3-binary`).  Import this module (or `tait`) before spherogram/snappy.
"""
import sys

try:
    import sqlite3  # noqa: F401
except ImportError:  # pragma: no cover
    try:
        import pysqlite3
        sys.modules['sqlite3'] = pysqlite3
    except ImportError as exc:
        raise ImportError('stdlib sqlite3 is missing and pysqlite3-binary is not installed: '
                          'pip install pysqlite3-binary') from exc
