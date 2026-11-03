# AI index: exact worst-case convergence of alternating projections (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) and the
paper [exact_worst_case_alternating_projections.pdf](exact_worst_case_alternating_projections.pdf).

## Normalisation

A, B closed convex in R^d (Theorem 1) or a real Hilbert space (Theorem 7). Common closed ball B(q,delta) ⊆ A∩B
(Theorem 1), or B(q,a) ⊆ A and B(q,b) ⊆ B with a common centre (Theorems 7-8). ||x0 - q|| <= R. Iterates
y0 = P_A x0, then alternate P_B, P_A, ...; x_N = (P_B P_A)^N x0 = y_{2N-1}. W_N(R,delta) = sup dist(x_N, A∩B)
(unsquared). kappa = 1 - 2 delta^2/R^2 = cos(2 arcsin(delta/R)). Source: Luner-Grimmer arXiv:2410.14811 v1/v2
Prop. 6.1 (lower construction) and Conj. 6.2 (stated for beta-smooth sets, max over c in [0,1], equivalent to [0,kappa]).

## Claims and where they are proved

| ID | Claim | Location in [manuscript.md](manuscript.md) | Evidence type | Priority status |
|---|---|---|---|---|
| C1 | W_N = 0 if R <= sqrt2 delta; else W_N = max_{0<=c<=kappa} c^{2N-1}(sqrt(R^2-delta^2) - delta sqrt((1+c)/(1-c))); unique maximiser; attained by two planar halfspaces; dimension-independent | Theorem 1; Sections 2-6 | proof | formula and lower construction: Luner-Grimmer; upper bound: not found earlier |
| C2 | One-step comparison r_{j+1}^2 <= r^2 - t^2 - 2 delta t, t_{j+2} <= c(r,t) t; final-state admissibility (17) | Lemma 2, (17) | proof + exact checks | new as stated |
| C3 | Sharp error bound dist(y, A∩B) <= E(r, dist(y,B)) for admissible states | Lemma 3 | proof + exact checks | sharp form not found; non-sharp antecedents cited |
| C4 | Monotone value function V_k, contraction under halfspace update (23) | Lemma 4 | proof + exact checks | new |
| C5 | Unique maximiser via quartic (28) | Section 6 | proof | evaluation of known formula |
| C6 | W_N ~ R^4 kappa^{2N}/(8 e delta^2 sqrt(R^2-delta^2) N); also when N delta^2/R^2 -> infinity | Corollary 5 | proof + numerics | consequence of C1 |
| C7 | Critical limit W_N/R -> L(lambda) when N delta^2/R^2 -> lambda; uniform additive profile; expansions (35) | Corollary 6 | proof | consequence of C1 |
| C8 | Unequal radii a, b (common centre), any number k >= 2 of projections, real Hilbert space: exact worst case max c^{k-1} D_AB(R,c), threshold R^2 <= a^2+b^2; k=1 equal radii: R - delta | Theorem 7, Section 8 | proof + exact checks + randomised search | new |
| C9 | Exact worst residual max c^m sqrt(1-c^2) D_AB(R,c); delta -> 0 limit equals Zamani-Glineur (arXiv:2407.15195, Thm 4.4/4.5) | Theorem 8, Section 8 | proof + high-precision check of limit | new; limit is a consistency check |
| C10 | Asymptotic (43) for unequal radii; order dependence of the constant | Section 8 | proof (as Cor. 5) + numerics | new |
| C11 | Classical baseline: per-projection factor sqrt(1 - delta^2/(R+delta)^2) | Proposition 9 | proof | classical-type argument, stated for comparison |
| C12 | Cycle counts (exact vs leading order vs (A4) vs classical) | Section 9 table | computation (checks/extension_checks.py) | illustrative |
| C13 | Fixed instances: N dist(x_N)/kappa^{2N} -> 0 | Corollary A1 | proof | new |
| C14 | Near-worst instances: dist(P_A x0, A∩B) = Theta(1/N) for N >= N0 | Corollary A2 | proof | new |
| C15 | dist(y_m, A∩B) <= kappa^m dist(P_A x0, A∩B)/sqrt(1-kappa^2) | Corollary A3 | proof | new as stated |

## Verification artefacts

See [VERIFICATION.md](VERIFICATION.md). `run_checks.py` runs four suites (88,902 explicit checks, valid under `-O`):
[checks/exact_rational_checks.py](checks/exact_rational_checks.py), [checks/second_exact_checks.py](checks/second_exact_checks.py),
[checks/check_convex_polygons.py](checks/check_convex_polygons.py), [checks/extension_checks.py](checks/extension_checks.py).
Outputs in [results/](results/). Figures: [code/plot_exact_results.py](code/plot_exact_results.py) -> [figures/](figures/).

## Review and provenance

- [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md): internal AI-assisted review, objections and dispositions.
- [PRIOR_ART.md](PRIOR_ART.md): dated prior-art and prior-credit search (3 Nov 2026), including GitHub issues.
- [crossref_records.jsonl](crossref_records.jsonl): Crossref records for the DOI references.
- [refs.bib](refs.bib): bibliography of the reference list, built from those Crossref records and arXiv API records.
- [RESEARCH_GATES.json](RESEARCH_GATES.json): contribution map, sources read, limitation triage and extension scout.
- [CITATION.cff](CITATION.cff) and [.zenodo.json](.zenodo.json): citation metadata; version DOI 10.5281/zenodo.23922446.
- [MANIFEST.sha256](MANIFEST.sha256): SHA-256 of every file; [.github/workflows/replay.yml](.github/workflows/replay.yml) reruns the checks (normal and -O).

## Limits

No human peer review; no formal verification; finite checks are not proofs; bounded priority search. Open:
three or more sets, relaxed projections, balls with different centres, smooth/strongly convex classes.
