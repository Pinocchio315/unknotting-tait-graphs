# Review-only common-prime connected-sum calculation

This standalone calculation is retained as review evidence. It is not a result
claimed in the revised v1.3 manuscript and contributes no manuscript macros or
prime-knot totals.

For a fixed prime `p`, the rank `r_p(K) = dim H1(Sigma_2(K); F_p)` is additive
under connected sum and satisfies `r_p(K) <= u(K)`. If the same prime reaches
the exact unknotting number of each summand, subadditivity supplies the reverse
inequality and proves `u(A # B) = u(A) + u(B)`.

Taking the maximum of `r_p` over primes does **not** preserve additivity. For
example, the groups `Z/3 + Z/15` and `Z/5 + Z/15` each need two generators,
whereas their direct sum needs three. The code tests the common-prime condition
for every counted pair and assigns a pair with several witnesses to its least
witness prime, avoiding double counting.

`homology_inputs_2026-09-08.json.gz` records the determinant and the cyclic
factors of the double-cover homology for all 12,965 knots in the study. Its
source CSV is the public `database_knotinfo` 2026.8.1 snapshot with SHA-256
`e37cfb46060b0ba5165a730e482275c8a82460c6a462653e662188cfe952598d`.
This agrees with the previously archived public release metadata. The extract
contains tabulated inputs; it is not an independent computation of homology.

`common_prime_pairs_2026-09-08.json` records the 1,972 eligible summands, their
exact input values and their sharp primes. Enumerating unordered pairs with
repetition gives **483,817 certificates**, including **269,032** for which at
least one summand has unknotting number at least two. Mirroring choices are not
counted separately. These are pairs of table labels, not distinct isotopy types
and not additional prime-knot determinations.

From `code/`, reproduce the pair analysis with no database dependency:

```text
python connected_sum_generators.py
```

To reproduce the input extraction, obtain the CSV member of the public archive
identified in `../comparison/releases/manifest.json`, then run:

```text
python connected_sum_generators.py --freeze-csv /path/to/knotinfo_data_complete.csv
```

The freezer rejects a source whose checksum differs. The standalone script reads the deposited extract and recomputes the pair
counts. Neither current nor historical manuscript reporting imports it.
