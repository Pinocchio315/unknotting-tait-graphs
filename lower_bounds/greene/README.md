# lower_bounds/greene — correction terms of the double branched cover from a knot diagram

The obstructions of Section 3 compare the correction terms of `Σ₂(K)` with the definite forms allowed by
an unknotting sequence.  For an alternating knot the correction terms come from the sharp Goeritz form.
This directory removes that restriction: when the reduced Khovanov homology over `F₂` has rank equal to
the determinant, the branched-cover spectral sequence makes `Σ₂(K)` an `L`-space, and Greene's
spanning-tree model computes the correction terms from any diagram.  The comparison with the definite
forms is then unchanged, so the same tests apply to non-alternating knots.

The method was written to decide `u(12n491)`, the knot left open by Brittenham and Hermiller's
counterexample to the Bernhard–Jablan conjecture (arXiv:1705.05985).

## Results obtained here

| knot | previous range | result | consequence |
|---|---|---|---|
| `12n491` | `[1,2]` | `u = 2` | with `u(12n288) = u(12n501) = 2` from `../montesinos/`, identifies `13n3370` as the counterexample |
| `13n3370` | `[1,3]` | `u = 2` | four correction terms remain ambiguous; all nine candidate vectors are obstructed |
| `13n1587` | `[1,2]` | `u = 2` | settles the second example of arXiv:1705.05985, §3 |
| `13n1669` | `[2,4]` | not obstructed | the rank-two test admits the form `(18,18,17)`; the value stays open |

Deposited records are in `../../results/bernhard_jablan/`.

## Files

| file | role |
|---|---|
| `qa_search.py` | quasi-alternating certificate by resolving crossings (`det L = det L₀ + det L₁`), an independent proof of the `L`-space property |
| `greene_dinv.py` | Greene's model: Kauffman states of a marked diagram, their absolute gradings (Theorem 4.1), `Spin^c` structures (§4.5) and the solitary states that generate the `E₁` page (Theorem 6.8) |
| `half_integral.py` | the half-integral surgery pattern of Ni and Wu, which excludes `u = 1` |
| `pin_dinv.py` | intersects the solitary-state pages of every marking and colouring: `d(t)` is one of the solitary gradings in the class `t`, the classes of two markings are matched by the automorphisms of `H² = Z/D` preserving the self-pairing of `c₁`, and `d(t) = d(-t)` |
| `greene_ranks.py` | the rank `n` comparison driven by supplied correction terms; the forms come from `../owens/` |
| `build_greene_targets.py` | freezes the sweep input (`data/greene_targets.json`): diagram, determinant, signature, range, mod-2 Khovanov vector |
| `greene_sweep.py` | the sweep: one subprocess per knot, resume-safe JSONL output, per-knot time limit |
| `slurm_greene.sh` | the array launcher, sixteen tasks by default |
| `make_server_package.py` | collects everything the sweep needs into `greene_server/` and a tar archive |
| `summarize_greene.py` | combines the result files of every task, checks the controls, and lists the bounds |

## The sweep

`build_greene_targets.py` selects the knots to which a test applies.  A knot qualifies when it is
non-alternating, its reduced mod-2 Khovanov homology has rank equal to the determinant, and
`H₁(Σ₂(K))` is cyclic.  The recorded range and the signature then choose the test: a range `[1, b]` gives
the Ni–Wu test for `u = 1`, and a range `[n, b]` with `|σ| = 2n` gives the rank `n` test, which excludes
`n` crossing changes and so raises the lower bound to `n + 1`.

    targets   1967   u1 1521, rank2 385, rank3 54, rank4 7
    controls  1452   u1 574, rank2 640, rank3 210, rank4 28

A control has a known unknotting number equal to `n`, so the obstruction must not fire on it; the sweep
reports an obstructed control as an implementation error and exits nonzero.

Build the server package and copy it across.  `make_server_package.py` collects the modules of this
directory, the form enumerations of `../owens`, and the frozen input file into one flat directory, so
the server needs only `numpy`, `sympy` and the standard library.  The KnotInfo database, SnapPy and the
rest of the repository are not required.

```sh
python make_server_package.py
scp greene_server.tar.gz SERVER:~/
ssh SERVER "tar xzf greene_server.tar.gz && cd greene_server && sbatch slurm_greene.sh"
```

Sixteen array tasks each take one sixteenth of the work.  A task sweeps its share of the controls first
and goes on to the targets only if none of them is obstructed, so an implementation error stops the job
before any target is computed.  Resubmitting the same command resumes: a knot already recorded in the
output file is skipped, so a wall-clock limit costs only the unfinished knots.

The rank-four enumeration is much slower than the rest.  To keep it out of the main job and run it
separately:

```sh
sbatch --export=ALL,SKIP=rank4 slurm_greene.sh
sbatch --export=ALL,TEST=rank4,ARRAY_SIZE=4 --array=0-3 slurm_greene.sh
```

`bash slurm_greene.sh --dry-run` prints the resolved commands without running them.  Each task writes
`runs/<set>/results_<set>_<task>_<size>.jsonl`, one line per knot.  Typical cost is a few seconds per
knot for the `u1`, `rank2` and `rank3` tests; the rank-four enumeration is much slower and has its own
submission line above.

Verdicts: `OBSTRUCTED` raises the lower bound, and `PASS`, `UNDECIDED`, `TIMEOUT` and `ERROR` supply no
bound.  A time limit is not a negative mathematical result.

A run split between the cluster and another machine leaves several result files.  `summarize_greene.py`
reads any number of them, keeps one record per knot (a conclusive verdict replaces an error or a
timeout, and two different conclusive verdicts for one knot are reported as a contradiction), states the
coverage against the frozen list, checks that no control is obstructed, and writes the bounds:

```sh
python summarize_greene.py runs/*/*.jsonl --out ../../results/greene/sweep.json
```

## Validation

`python greene_ranks.py --validate 5_1 7_3 7_5 8_2 9_10 9_13 7_1 9_3 9_6 9_9 10_2 10_46` recomputes the
correction terms of twelve alternating knots with Greene's model, checks that they agree with the sharp
Goeritz form up to relabelling, and compares the verdict with the published rank-two and rank-three
tests of `../owens/`.  All twelve agree.  On the Montesinos knots `8_20`, `9_43` and `9_44` the pinned
correction terms agree with the star plumbings of `../montesinos/`, and the surgery test admits the
knots of unknotting number one (`3_1`, `4_1`, `5_2`, `8_20`, `9_44`) and excludes `5_1`, `7_4`, `9_43`.
The correction terms of `12n491` were recomputed from an independent minimal diagram and agree.
