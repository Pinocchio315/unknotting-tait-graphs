# Independent review of the new Bernhard–Jablan discussion

Review date: 2026-09-08. Scope: additions in `main_v1.3.tex` relative to v1.2,
the new 13n1587 correction terms, and the full crossing-neighbour set of
14a2539. This review did not modify either manuscript or production code.

## Findings requiring manuscript changes

1. **Incorrect grading sentence.** The new background paragraph says that the
   spectral-sequence differentials preserve both the Maslov grading and the
   Spin^c structure. They preserve the Spin^c summand and have Maslov degree
   minus one. The preceding paragraph already says this correctly and supplies
   the necessary argument: each page is graded, and taking homology cannot
   create a grading absent from the previous page. Delete the redundant new
   sentence or retain only its second observation about supplying correction
   terms to the surgery obstruction. This wording error does not invalidate
   the existing solitary-state candidate-set inference.

2. **Correction terms are not automatically uniquely determined.** Equality
   of reduced mod-2 Khovanov rank and determinant certifies the L-space
   hypothesis. The solitary-state procedure then supplies finite candidate
   sets; it need not determine every correction term. The paper's own
   13n3370 example illustrates this distinction. Replace claims that Greene's
   model supplies the complete list for every such knot by claims about
   exact values or certified candidate sets, with obstruction only when
   every surviving possibility is excluded.

3. **The 13n1587 computation does not identify a unique strong counterexample.**
   It proves u(13n1587)=2 and narrows the strong counterexample to that knot
   or a member of its BJ set with unknotting number one. No exhaustive BJ set
   for this nonalternating knot is supplied. Replace “locates the
   counterexample” and “settles the second example” by this precise scope.

4. **Do not omit the composite member for 14a2539.** Its BJ set has six
   members, not five. Scharlemann's theorem gives u(8_6 # 4_1)>=2, but does
   not establish equality three. The conclusion that the u=3 alternative
   moves the question to the five prime knots omits a logical possibility.
   A rigorous statement requiring no new theorem about composite minimal
   diagrams is: if u(14a2539)=3, either 14a2539 fails the strong equality,
   one of the five prime neighbours has unknotting number two and fails it,
   or u(8_6 # 4_1)=2. In the last possibility, the depth lemma still bounds
   the possible further descent. No claim about the composite's strong BJ
   number is necessary.

5. **Priority language.** The conditional phrase “is the first alternating
   knot known to do so” is a claim about the state of the literature rather
   than a consequence of the proposition. “Would give an alternating
   counterexample” states the proved implication without a priority claim.

6. **Stale upper-bound descriptions.** The deposited
   `results/13n_1587_u_le_2_braid_certificate.json` and
   `upper_bounds/verify_braid_certificate.py` still describe the exact value
   as unresolved, and attribute the upper-bound construction to Applebaum.
   The construction was already given in Brittenham–Hermiller §3. The
   certificate and verifier remain valid upper-bound computations, but their
   status and attribution should be updated. Section 4's sentence saying
   the range becomes [1,2] should explicitly say this is the upper-bound
   contribution before combining with the new lower bound.

## Mathematical checks

The new depth lemma is valid. If a knot with u<u_BJ^w satisfies the strong
equality, a minimizing BJ neighbour lowers u by exactly one. The weak
recurrence implies that this neighbour still has u<u_BJ^w. Iteration cannot
reach the unknot with that strict inequality, so it must stop at a strong
counterexample after fewer than the initial u steps. Cite Brittenham–Hermiller
Lemma 2.1 for the recurrence; the argument is a finite refinement of their
Corollary 2.2.

Owens's general Theorem 5 is not restricted to alternating or Montesinos
knots, nor does it require an L-space. It assumes the specified signed
unknotting sequence and uses the actual labelled correction terms of the
branched cover. Applying it to nonalternating knots is therefore valid once
the correction-term inputs and the signature-selected orientation are
justified. For a proposed sequence of n changes with |sigma|=2n, the signature
forces every change to have the relevant sign after mirroring. The L-space
hypothesis is used upstream to turn solitary-state gradings into candidates
for the correction terms, not as a new hypothesis of Owens's obstruction.

The blanket assertion in `greene_ranks.py` that Greene's output always has
d(spin)=-sigma/4 is too broad for arbitrary L-space covers. This equality
should be an explicitly checked consistency condition for the particular
inputs, not the general theorem used to infer orientation. The current guard
returns `ERROR` when the equality fails, so this issue makes the routine more
restrictive; it does not on its own create a false obstruction. Alternating
and quasi-alternating equality is supported by Manolescu–Owens and
Lisca–Owens respectively.

For 14a2539, signature two permits both two-negative and mixed-sign
possibilities. Owens's mixed-sign test alone cannot rule out u=2. In addition,
the mixed-sign candidate (m1,m2,a)=(1,154,0) really survives. The statement
that the obstruction is inconclusive is correct; no stronger lower bound
should be inferred from this test.

## Replayed computations

### 13n1587

Using the PD code and mod-2 Khovanov vector stored in the existing result file,
`pin_dinv.result_record` recomputed all 52 pages and the two intersection
passes. It reproduced all 75 labelled correction terms exactly. The mod-2
rank is 75, equal to the determinant. The independent replay is stored as
`13n1587_replay.json`.

All 2*phi(75)*75=6000 affine comparisons were rerun. For one orientation,
all 3000 fail the nonnegative-even-gap condition. For the other, 2996 fail
that condition and four fail consistency of repeated V indices. There are
no surviving fits. This is a replay of the supplied algorithm, not an
independent implementation of Greene's state-grading formula.

`upper_bounds/verify_braid_certificate.py` also completed successfully:
the two braid words differ at one generator, their closures match 13n1587
and 10_113, the changed Tait diagram matches the partner, the stored reduction
path preserves the partner, and a crossing of 10_113 produces genus-zero
knot Floer homology. These knot identifications use numerical SnapPy exterior
isometries, as the verifier states. Together with the independently published
upper bound, the lower-bound replay supports the exact value two.

### 14a2539

`replay_bjset14a2539.py` constructs the source independently from Spherogram's
named table and changes each of its 14 crossings. The resulting set is:

| Member | Determinant | Number of source crossings |
|---|---:|---:|
| 10_11 | 43 | 2 |
| 12a178 | 139 | 3 |
| 12a183 | 121 | 3 |
| 12a196 | 127 | 2 |
| 12a684 | 135 | 2 |
| 8_6 # 4_1 | 115 | 2 |

All named identities were reproduced by SnapPy. Except for 12a684, the
Reidemeister-simplified prime factors also have exactly the same Regina
canonical diagram signature as the independent named reference. The remaining
12a684 identification was then certified by exact diagram moves: for both
source crossings (zero-based rows 0 and 5), two pass moves reduce the changed
diagram from 14 to 12 crossings, and one flype gives full canonical signed
planar-map equality with the independent named 12a684 reference, up to mirror
and turning over. The deterministic search took 0.238 seconds and visited
12 diagrams for each row. `exact_12a684_identification.py` and
`12a684_exact_identification.json` preserve the input and endpoint PD codes,
all pass/flype parameters, and the matching full canonical codes. Thus none
of the fourteen neighbour identifications depends solely on a numerical
isometry. All changed PD codes and initial identifications are in
`bjset14a2539_replay.json`. The finite replay does not itself prove
completeness over all minimal diagrams; the flyping lemma supplies that
separate mathematical step.

All five prime members have comparison range [2,3] and occur in the existing
dichotomy appendix. Their zero-candidate status is checked against
`git show 6468494:results/open23/priority_u23_final.json`, the committed input
at the start of this review, rather than concurrent edits to the live scan.
The stored comparison values of the composite factors
are u(8_6)=2 and u(4_1)=1. Thus Scharlemann gives the stated lower bound two
for the composite, and subadditivity gives upper bound three.

The positive Goeritz form of rank eight has determinant 307 and spin
correction term -1/2. Its labelled correction terms admit the mixed-sign
candidate (1,154,0); `14a2539_survivor.json` records this replay. The other
checkerboard's spin value is +1/2, as expected on orientation reversal.

## Primary references checked

- Brittenham–Hermiller, *A counterexample to the Bernhard–Jablan unknotting
  conjecture*, arXiv:1705.05985v2: Lemma 2.1 and Corollary 2.2, printed
  pp. 4–5; flype Lemma 2.4, p. 7; the two further pairs, §3, p. 11.
  https://arxiv.org/pdf/1705.05985
- Greene, *A spanning tree model for the Heegaard Floer homology of a
  branched double-cover*, arXiv:0805.1381: Theorem 4.1 and §4.5 for gradings
  and Spin^c labels; §6.1, pp. 31–32, for the filtration; Theorem 6.8,
  p. 36, for solitary states; §7.1, pp. 37–38, for extracting correction
  terms and changing markings. https://arxiv.org/pdf/0805.1381
- Owens, *Unknotting information from Heegaard Floer homology*,
  arXiv:math/0506485: general Theorem 5, p. 4; proof, p. 19.
  https://arxiv.org/pdf/math/0506485
- Lisca–Owens, *Signatures, Heegaard Floer correction terms and
  quasi-alternating links*, arXiv:1302.3190, for the quasi-alternating
  signature identity. https://arxiv.org/abs/1302.3190
