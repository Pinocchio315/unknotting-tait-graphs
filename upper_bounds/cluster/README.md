# Enumeration of equivalent diagrams

This directory contains the enumeration search used with the upper-bound constructions and candidate analysis in Section 4 of *Computation of Unknotting Numbers: Which Knot Breaks the Bernhard–Jablan Conjecture?*. Within a chosen crossing budget, the search expands diagrams by pass moves, Reidemeister III moves, flypes, and optional Reidemeister II additions. Canonical codes avoid repeated screening. The search can propose a sequence of crossing changes to the unknot or to a partner with a known upper bound.

The search kernels use the Python standard library. Local endpoint recognition and certificate verification additionally require the topology environment described in the top-level README.

| File | Purpose |
|---|---|
| `enum_search.py` | Enumerate reachable diagrams, screen crossings, and record candidate constructions |
| `moves_io.py`, `xtait/` | Record and replay the allowed diagram isotopies |
| `flipdet.py` | Modular determinant screening |
| `lickorish.py` | Exact rational linking-pairing filter using Sherman–Morrison updates |
| `alexander.py` | Alexander-polynomial fingerprint for candidate screening |
| `resolve_candidates.py` | Compare candidate endpoints with named knots or detect the unknot |
| `verify_enum_certificates.py` | Authenticate source diagrams, replay moves, and verify marked crossing changes |
| `build_data.py` | Prepare exploratory targets and partner codes from the installed database |
| `summarize.py` | Summarize recorded search coverage and timings |
| `slurm_enum.sh`, `slurm_single.sh`, `slurm_u23_children.sh` | Archived cluster launch configurations |

## Verification and search scope

From this directory, verify an existing result directory with:

```sh
python summarize.py /path/to/run
python verify_enum_certificates.py /path/to/run --known-bounds ../../results/paper_v1_1_snapshot.json
```

The verifier requires exact agreement of the recorded and replayed PD codes; sorting the four labels at a crossing would discard over/under information and is not a valid comparison. The marked crossing count and partner upper bound are checked independently. A determinant-one endpoint must still be recognized as the unknot. The linking-pairing filter passes determinant-one cases because their torsion pairing supplies no obstruction. Numerical SnapPy isometries are computational identification checks, not formal proof certificates.

A completed finite enumeration concerns only the implemented moves and its recorded crossing budget. Timeouts and unsuccessful searches give no lower bound on the unknotting number. The child-knot campaign seeks an unknotting crossing for possible partners of the open alternating knots; a verified partner with `u = 1` can yield an upper bound of 2 for a parent, while a failed search proves nothing about the parent's unknotting number.

For exploratory runs, inspect a launcher before submitting it, for example `bash slurm_enum.sh --dry-run`. Search-data builders read the installed database and must not be confused with the frozen manuscript snapshot. Historical partner lists, single-knot configurations, and logs remain archived under their original filenames; their presence does not add a claim to the paper. New verification reports should be written outside the deposited `results/` files.
