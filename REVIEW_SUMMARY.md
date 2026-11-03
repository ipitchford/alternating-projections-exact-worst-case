# Review summary (public, internal AI-assisted review)

Date: 3 November 2026. This summarises an internal review carried out with AI reviewers. It is not external human
peer review or independent reproduction.

## Round 1

Five reviewer roles (editor, methodology, domain, interdisciplinary perspective, devil's advocate) each read the
first version independently. An adversarial review was run with a second AI model. The lead re-ran all shipped
scripts from a clean copy.

**Correctness.** No reviewer found an error or gap in the proof of Theorem 1. Each traced Lemmas 2-4 and the
telescoping step, including the edge cases (one cycle, zero residuals, the threshold R ≤ √2δ). Independent numerical
attacks found no instance exceeding the formula: about 1,000 optimised searches, 4,000 random polygon pairs, and
34,425 polygon and polytope evaluations in two and three dimensions. The best attempts came within 0.1% of the bound,
and the closest was the known two-halfspace example itself.

**Main objections and how they were handled.**

| Objection | Disposition |
|---|---|
| The literature review was far too thin. Classical interior-point rates (Gubin-Polyak-Raik), the subspace angle identity (Kayalar-Weinert), error bounds, and the performance-estimation literature, including Zamani-Glineur's exact alternating-projection residual rate, were all missing. | New Section 1.2 and 18 references, each checked against Crossref or arXiv. |
| There was no classical baseline, and the Figure 2 caption referred to "general upper estimates" that the paper never gave. | Proposition 9 now gives the classical bound in the paper's normalisation, Figure 2 plots it, and a cycle-count table was added. |
| The 1/N factor comes from measuring the start from the ball centre. Near-worst instances start Θ(1/N) from the intersection, and this was hidden in the appendix. | This is now stated in the abstract and Section 9. New Corollary A3 gives geometric contraction relative to the initial distance. |
| Several cheap generalisations were left unstated: any number of projections and real Hilbert space. | Added as part of Theorem 7. |
| Unequal ball radii looked tractable with the same method. | Proved as Theorem 7, with an explicit formula, an order-dependence remark and asymptotics (43). |
| The exact worst residual and its small-ball limit were missing. | Proved as Theorem 8. Its limit as δ → 0 equals Zamani and Glineur's exact rate, with the iteration counts matched. |
| The source maximises over [0,1] and its table reports squared distances. | Both normalisations are now explained in Section 1.2. |
| Version history: v2 of the source keeps the conjecture, and the journal version omits it. | Both are now cited. |
| Terse steps and overloaded notation. | One-line justifications were added and symbols renamed. |
| Check scripts relied on `assert`, so they silently passed under `python -O`. | All checks are now explicit. Planted errors make every suite fail, also under `-O`. |
| Three or more sets. | Stated as open. The method does not extend because inequality (7) fails. |

## Verification round

A proof auditor checked all new and rewritten proofs line by line. A second-model audit checked the new proofs
again, and a re-reviewer checked every round-1 item, every reference, the replay and the PDF. None found a
correctness issue. Their adversarial searches against Theorems 7 and 8 (about 22,000, 28,800 and 24,000
comparisons) found nothing above the bounds.

The remaining minor points were fixed in the final text:

- Wording on projection counts.
- The admissibility hypothesis in the stopping test.
- "About half" now carries the condition δ ≪ r.
- The Θ(1/N) statement now refers to the iterate after the first projection.
- An explicit uniform-limit argument for the Zamani-Glineur limit.
- The citation year of Drusvyatskiy, Li and Wolkowicz.
- Proposition numbering.
- The one-projection bound now uses the positive part.
- A fairer comparison column in the cycle-count table.

## What the review does not establish

There is no human refereeing and no formal (proof-assistant) verification. The prior-art search was bounded and is
not a priority certificate.
