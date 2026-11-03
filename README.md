# Exact worst-case convergence of alternating projections with a common interior ball

Version 0.1.0-candidate (3 November 2026). Author: Anonymous. Text and data CC0-1.0; code MIT (see LICENSE.md).
Unrefereed candidate research paper.

## In plain English

Alternating projections look for a point that satisfies two convex constraints by repeatedly jumping to the nearest
point of one set, then of the other. Suppose both sets contain a common ball of radius δ and the starting point is
within distance R of the ball's centre. How far from a solution can the method still be after N rounds, in the worst
case? Luner and Grimmer (arXiv:2410.14811, versions 1-2) found a pair of half-planes that is this slow and conjectured
that nothing is worse. This paper proves their conjecture, for every N and in every dimension, with a short
hand-checkable argument. It also:

- extends the exact answer to sets containing balls of **different** radii about a common centre, to any number
  (at least two) of single projections, and to any real Hilbert space (Theorem 7);
- gives the exact worst-case **residual** (distance to the other set), which reproduces the exact rate of Zamani and
  Glineur (2024) when the ball shrinks to a point (Theorem 8);
- turns the formula into practice: sharp asymptotics, cycle counts (for R/δ = 10 and accuracy 10^-3: 93 cycles in the
  worst case, against 833 guaranteed by the classical interior-point bound), and a computable error bound.

The formula and the extremal half-planes are Luner and Grimmer's; the contribution is the matching upper bound and
the extensions.

## Reading order

1. `exact_worst_case_alternating_projections.pdf` (or `manuscript.md`): the paper. Section 1 gives context and a
   roadmap; Sections 2-5 prove the upper bound; Section 8 the extensions; Section 9 the practical reading.
2. `REVIEW_SUMMARY.md`: what the internal AI-assisted review found and how each point was handled.
3. `PRIOR_ART.md`: the dated prior-art and prior-credit search.
4. `VERIFICATION.md` and `AI_INDEX.md`: what the code checks, and a claim-by-claim map.
5. `RESEARCH_GATES.json`, `refs.bib`, `CITATION.cff`, `MANIFEST.sha256`: machine-readable scrutiny record,
   bibliography, citation metadata and file hashes. `requirements.txt` pins the figure dependencies.

## Reproduce the checks

Standard library only (Python 3.10 or later; recorded with CPython 3.11.14):

~~~bash
python3 -I -B run_checks.py --out results     # also valid with -O
~~~

Four suites run (88,902 explicit checks): exact rational checks of the scalar identities and inequalities, an
independently written second exact suite, exact projections onto 160 polygon pairs, and the extension suite
(two-radius identities and inequalities, attainment, a randomised polygon search, the Zamani-Glineur limit,
asymptotics and the cycle-count tables). Every check is explicit and the runner exits non-zero on any failure, also
under `python -O`. These are finite checks; the theorems rest on the proofs in the paper.

Figures: `python3 -I -B code/plot_exact_results.py` (needs numpy 2.3.5, matplotlib 3.10.8). PDF:
`python3 -I -B code/build_pdf.py` (needs pandoc and XeLaTeX with Latin Modern fonts).

## Status

Internal AI-assisted review only (five reviewer roles, an adversarial review by a second model, a proof audit and a
re-review). No external human peer review, no proof-assistant verification. A bounded search found no earlier proof;
this is not a certificate of priority.

## How to cite

Anonymous (2026). *Exact worst-case convergence of alternating projections with a common interior ball*. Version
0.1.0-candidate. Evidence Press. https://doi.org/10.5281/zenodo.23922446 (concept DOI for all versions: https://doi.org/10.5281/zenodo.23922445).
Release page: https://evidencepress.org/releases/alternating-projections-exact-worst-case/ . Machine-readable
citation metadata: `CITATION.cff`.

## Source

Luner, A., & Grimmer, B. (2024). Performance estimation for smooth and strongly convex sets (arXiv:2410.14811v1,
Section 6, Proposition 6.1 and Conjecture 6.2). https://arxiv.org/abs/2410.14811v1
