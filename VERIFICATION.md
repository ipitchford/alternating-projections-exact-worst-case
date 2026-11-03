# Verification record

The theorems rest on the analytic proofs in the paper. The code below checks identities, inequalities, the attaining
construction and many concrete instances in exact or high-precision arithmetic. It can reveal errors; it does not
prove universally quantified statements. No proof-assistant verification was performed.

Run: `python3 -I -B run_checks.py --out results` (standard library; also valid with `-O`). Recorded run: CPython
3.11.14, macOS arm64, 3 November 2026: 88,902 explicit checks, 0 failures, both with and without `-O`. Each suite
exits non-zero on any failed check. Planted single-inequality errors made every suite fail under both modes.

| Suite | What it checks | Arithmetic | Checks |
|---|---|---|---|
| `checks/exact_rational_checks.py` | inverse state (9), halfspace update (12), transition polynomials (16), value comparison (24), radius monotonicity, normal decompositions and metric inequalities of Lemma 3 | exact rationals | 15,717 |
| `checks/second_exact_checks.py` | independently written: state identities incl. endpoints, pointwise value comparisons, quadratic reductions, two-normal metric cases | exact rationals | 26,676 |
| `checks/check_convex_polygons.py` | exact projections onto 160 pairs of polygons/halfspaces containing the unit ball, 8 horizons each, against (3); largest ratio 0.99940 | exact rationals, 80-digit decimals | 4,688 |
| `checks/extension_checks.py` | Section 8 functions for both orders and unequal radii: inverse, identities (41), inequality (42), kappa comparison, monotonicity of V; attainment of Theorems 7-8; randomised search over 720 polygon / disc-halfspace instances, m = 1..5, both orders (largest ratio 0.99862); Zamani-Glineur limit; asymptotic (43); cycle-count tables (Section 9); two-radius vs equal-radius comparison | exact rationals, 60-digit decimals | 41,821 |

Outputs: `results/*_log.txt`, `results/*.json`, `results/run_summary.json` (records the manuscript SHA-256).
Figures: `code/plot_exact_results.py` regenerates `figures/*.png` and their CSV data (binary floating point, for display).

Internal review computations (not shipped) added further adversarial searches: about 1,000 optimised polytope
searches, 4,000 random polygon pairs, 34,425 polygon/polytope evaluations against Theorem 1, and about 22,000,
28,800 and 24,000 comparisons against Theorems 7-8. None exceeded the bounds.
