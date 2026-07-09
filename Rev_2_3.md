# Revision — Editor/Referee comments: near-flat likelihood regions & cross-derivative standard errors

## The comments being addressed

**Editor (near-flat regions):** Related to the optimization strategy, the referee also
raises concerns about the potential presence of near-flat regions of the likelihood
surface, where multiple parameter configurations may yield very similar likelihood
values. The manuscript would benefit from a clearer discussion of whether such
behavior has been encountered in practice, how frequently it occurs, and how users
should diagnose or handle such situations. More explicit guidance on the role of
random initialization, parameter bounds, and stopping rules would improve robustness
and reproducibility.

**Referee (near-flat regions):** In addition, given the type of stopping rule used,
could the authors comment on whether they have encountered near-flat regions of the
log-likelihood surface, where multiple parameter configurations yield very similar
likelihood values?

**Editor (cross-derivatives):** Standard errors are currently computed using a
finite-difference approximation to only the diagonal elements of the Hessian matrix,
thereby ignoring cross-derivatives. As noted by the referee, such cross-terms may be
non-negligible in Cox frailty models, potentially affecting the accuracy of standard
errors and confidence intervals. The authors are encouraged to provide stronger
justification for this approximation or to consider alternatives such as computing
the full Hessian at convergence or adopting sandwich or robust variance estimators.
Clarifying the inferential implications of the current approach would be valuable
for applied users.

**Referee (cross-derivatives):** In the manuscript, the standard errors are computed
using a finite-difference approximation to only the diagonal of the Hessian. This
yields a diagonal approximation to the observed information and therefore ignores
cross-derivatives, which may be non-negligible in Cox/frailty models. Consequently,
the reported standard errors and confidence intervals may be inaccurate. Could the
authors provide further justification for neglecting the cross-derivatives, or
consider alternatives such as computing the full Hessian at the final estimate or
using a sandwich/robust variance estimator?

---

## Reasoning

### Near-flat regions

This is directly answerable with the `GradientCheck`/`BoundaryAdjacent` diagnostic
built for the earlier optimality-conditions comment (`Rev_1.md`, `REVISION_NOTES.md`
§7). Today's (`n_extrarun=10`) investigation already found a structured, interpretable
instance of near-flatness rather than an arbitrary one: `mu1` shrinking towards its
lower bound drags `nu` along a ridge (`mu1/nu ~= 0` regardless of `nu`'s value once
`mu1` is negligible — see `REVISION_NOTES.md` §7bis), and several `gamma_k` shrink
towards their lower bound in intervals where the data carry little information about
between-cluster heterogeneity. In all such cases the diagnostic shows the *interior*
parameters have a small gradient (genuine stationary point) while only the
boundary-adjacent dispersion parameters show a large, expected, one-sided gradient —
i.e. near-flatness is confined to specific, identifiable directions, not diffuse
across the whole parameter space.

### Cross-derivatives

The referee's own suggested alternative — a full Hessian computed once at the final
estimate, not at every iteration — is directly buildable by extending the same
finite-difference machinery already used for `GradientCheck` and the existing
diagonal-only `params_se()`. This is a read-only diagnostic computed after
optimization finishes (does not touch the optimizer or any previously reported
number), so it carries the same low-risk profile as today's other additions. A
sandwich/robust estimator is a worse fit for this setting: the model is a fully
parametric likelihood, not a working-independence/GEE-type model that a sandwich
correction is designed to fix, so a correctly computed model-based (inverse-Hessian)
covariance is the more directly applicable option here — noted as such in the draft
response below.

One expected complication, confirmed once the full Hessian was actually computed
(see "Results" below): parameters sitting at/near a declared boundary (`mu1`, `nu`,
several `gamma_k`) can make the *full* Hessian ill-conditioned or non-negative-definite
overall, since a near-zero-curvature boundary direction couples into every other
row/column once the matrix is inverted — unlike the diagonal-only approach, which is
immune to other parameters' degeneracy by construction. Standard asymptotic
likelihood theory (and a naive full-Hessian inversion) does not directly apply at a
boundary solution; this is itself a relevant, honest point to make in the response
rather than to gloss over.

**On whether independence assumptions in the paper already justify the diagonal
shortcut** (author question, 2026-07-08): checked directly against the paper text.
The only independence statements are model-level: `alpha_j` and `eps_jk` (Section 3.3,
Appendix B) are independent *random variables*, and units within a group are
conditionally independent given the frailty. Section 4.4's own stated justification
for the diagonal-only Hessian is purely computational ("...would be computationally
expensive and inefficient... more efficient since the standard error formula
requires only the diagonal elements" — no independence argument is invoked there.
Model components being independent random variables does not imply their *estimators*
(`mu1_hat`, `nu_hat`, `gamma_k_hat`, `beta_hat`) have uncorrelated sampling error —
those are different claims, and today's own Hessian computation contradicts the
latter directly: `beta` and `phi` are correlated through the shared `A_j..` term
(material SE changes below), and the `mu1`/`nu` block isn't even negative-definite.
So if an independence-based justification was intended anywhere, it would not
actually hold up against the paper's own model structure — worth being explicit about
this distinction in the response rather than leaning on it.

---

## Draft responses

### To the near-flat-regions comment

> We thank the editor and referee for this comment, which is closely related to the
> request for a first-order optimality check addressed above. To investigate
> directly, we extended the finite-difference gradient diagnostic introduced in
> response to that comment to separate the model's parameters into those away from
> any declared range boundary ("interior") and those close to one
> ("boundary-adjacent").
>
> Applying this to the paper's worked example, near-flat regions do occur, but in a
> structured, interpretable way rather than at arbitrary points of the parameter
> space. In particular, we observe cases where $\mu_1$ shrinks towards its lower
> bound; since $\mu_1/\nu \approx 0$ regardless of the value of $\nu$ once $\mu_1$ is
> negligible, $\nu$ becomes only weakly identified along a ridge in the
> $(\mu_1,\nu)$-plane — the log-likelihood is nearly flat along that direction,
> consistent with $\nu$'s large standard error already discussed in Section 7.
> Similarly, some $\gamma_k$ associated with time intervals in which few or no
> clusters register any event shrink towards their lower bound, consistent with a
> variance component for which the data carry little information, rather than with a
> failure to converge. In all such cases, the *interior* parameters (regression
> coefficients, baseline log-hazards, and the remaining, better-identified
> $\gamma_k$) show a small finite-difference gradient, consistent with a genuine
> interior stationary point, while only the boundary-adjacent dispersion parameters
> show a large, one-sided gradient — the expected signature of a boundary solution
> rather than evidence against convergence.
>
> We now report this diagnostic by default in the package output (`print()`/
> `summary()` of an `AdPaik` object), together with an explicit flag identifying
> boundary-adjacent parameters, and expand Section 4.3/6 with guidance for users:
> (i) *random initialization* — `n_extrarun` already governs how many additional
> random restarts are attempted; we now recommend inspecting whether the reported
> optimum is stable across restarts, particularly when boundary-adjacent parameters
> are flagged; (ii) *parameter bounds* — tighter, data-informed bounds (as already
> recommended via a preliminary time-independent fit, Section 6.1) reduce the chance
> of a spurious boundary solution for a parameter that is in fact well identified;
> (iii) *stopping rule* — convergence is assessed on the log-likelihood alone, so we
> recommend the new gradient diagnostic as a complementary check, treating
> boundary-adjacent parameters as requiring domain judgement (is a near-vanishing
> dispersion plausible here?) rather than as a diagnostic failure.
>
> To directly address "how frequently" this occurs and the role of random
> initialization, we refit the worked example from 5 different random seeds
> (`Examples/MultiSeedStabilityCheck.R`). The log-likelihood at convergence agreed to
> within $3\times10^{-5}$ across all 5 seeds, the regression coefficients agreed to
> within $10^{-4}$, and — notably — the *same* set of 6 boundary-adjacent parameters
> ($\mu_1$, $\nu$, $\gamma_1$, $\gamma_4$, $\gamma_9$, $\gamma_{10}$) was identified in
> every single run, despite the number of coordinate sweeps needed to converge
> varying (23 to 35). For this dataset, the boundary/near-flat behavior is therefore
> a robust, reproducible feature of the likelihood surface given the data, not an
> artifact of a particular random start; we obtained no evidence of the optimizer
> landing on materially different, competing optima.

Confirmed at the paper's own `n_extrarun=60`/seed-1 setup (see "Results" below): same
boundary-adjacent set (`mu1`, `nu`, `gamma_1`, `gamma_4`, `gamma_9`, `gamma_10`), same
small interior gradient (`‖grad‖_2 = 0.0221`, `max|grad| = 0.0162`) — the pattern found
earlier at the reduced `n_extrarun=10` setting was not an artifact of that reduction.

### To the cross-derivatives comment

> We agree that the diagonal-only Hessian approximation ignores off-diagonal
> curvature that may be non-negligible, and thank the referee for raising this. The
> diagonal approximation was adopted for computational efficiency (Section 4.4): a
> full Hessian recomputed at every iteration of a multi-start, coordinate-wise
> optimization would be considerably more expensive. However, this cost only applies
> during optimization — at the *final* reported optimum, a full Hessian needs to be
> computed only once. We therefore implemented a full numerical Hessian at
> convergence, via the same centered finite-difference approach already used for the
> diagonal (Algorithm 5), extended to off-diagonal pairs, and inverted it (where
> well-conditioned) to obtain a proper, correlation-aware covariance matrix.
>
> For the worked example of Section 6, restricting to the parameters not sitting at a
> declared range boundary (18 of 24; see the discussion of the previous comment), the
> resulting covariance matrix is positive definite, and the ratio of the
> cross-derivative-corrected standard error to the diagonal-only one ranges from 1.00
> to 2.11 across these 18 parameters (median 1.12, mean 1.21). The effect is
> parameter-type specific: the baseline log-hazards $\phi_k$ increase by roughly
> 10-44%; the regression coefficients are the most affected, with the standard error
> of $\widehat{\beta}_{\text{CFUP}}$ increasing by about 45% (from 0.0345 to 0.0502)
> and that of $\widehat{\beta}_{\text{GenderMale}}$ more than doubling (from 0.0518 to
> 0.1093, a 111% increase); the (non-boundary-adjacent) $\gamma_k$ are essentially
> unaffected (at most 6%). Recomputed with the corrected standard error, the 95\% CI
> for $\widehat{\beta}_{\text{GenderMale}}$ widens from $(0.117, 0.320)$ to
> approximately $(0.004, 0.432)$ — the estimated effect remains distinguishable from
> zero, but far less comfortably so than the diagonal-only approximation suggested.
> We therefore agree with the referee that cross-derivatives are not negligible here,
> and are revising the standard errors reported in Section 6.2 (and the corresponding
> discussion in Section 4.4/7) accordingly.
>
> For the parameters directly involved in the near-flat region discussed above
> ($\mu_1$, $\nu$), the raw $2\times2$ Hessian sub-block is not negative definite, so
> no valid correlation/standard error can be reported for this pair from a
> finite-difference Hessian at all — a further, independent illustration that
> standard (Wald-type) asymptotic inference does not apply at this boundary solution,
> rather than an artifact of the diagonal-only shortcut specifically.
>
> We did not pursue a sandwich/robust variance estimator in this revision: the model
> is a fully parametric likelihood-based one, without the working-independence/
> GEE-type structure a sandwich correction is designed to fix, so a correctly
> computed model-based (inverse-Hessian) covariance is the more directly applicable
> improvement here. We note a robust variant as a possible direction for future work.

---

## Results

Run: `Examples/ReplicationCode.R` settings exactly (`set.seed(1)`, default
`n_extrarun=60`, `formula <- time_to_event ~ Gender + CFUP + cluster(group)`,
`categories_range_max <- c(-eps, 0.5, 1-eps, 1, 10)`), `TimeDepFrail` at the state of
branch `revision` (log-space log-likelihood + `gradient_check()` + `hessian_check()`,
this session). Fit time ~474s, full-Hessian diagnostic ~34s on top.

**Fit** (matches the paper closely; small residual differences are the expected
effect of the mathematically-equivalent log-space log-likelihood reformulation,
§7 above): `Loglikelihood = -2175.049` (paper: `-2175.135`), `AIC = 4398.098` (paper:
`4398.27`), `coef = (GenderMale=0.2181392, CFUP=-1.270622)` (paper: `0.217802,
-1.270657`), `NRun = 31` (paper: `31`, exact match).

**GradientCheck**: `GradientNorm (interior) = 0.02208516`, `GradientMaxAbs (interior)
= 0.01616302`. `BoundaryAdjacent = [13, 14, 15, 18, 23, 24]` = `mu1`, `nu`, `gamma_1`,
`gamma_4`, `gamma_9`, `gamma_10` — identical set and near-identical gradient values to
the earlier `n_extrarun=10` exploratory run, confirming the pattern is not an artifact
of that reduced setting.

**hessian_check**: `InteriorIndices = [1..12, 16, 17, 19, 20, 21, 22]` (18 of 24: all
`phi_k`, both `beta_r`, and the 6 `gamma_k` not flagged boundary-adjacent).
`IsPositiveDefinite (interior submatrix) = TRUE`. `Mu1NuBlockPositiveDefinite =
FALSE`, hence `CorrelationMu1Nu = NA` (not computable — see reasoning above).

`SE_full / SE_diag` by parameter (interior only):

| index | parameter | SE_diag | SE_full | ratio |
|---|---|---|---|---|
| 1 | phi_1 | 0.11471 | 0.15611 | 1.361 |
| 2 | phi_2 | 0.16550 | 0.19651 | 1.187 |
| 3 | phi_3 | 0.10010 | 0.14419 | 1.440 |
| 4 | phi_4 | 0.20412 | 0.22886 | 1.121 |
| 5 | phi_5 | 0.14846 | 0.18030 | 1.214 |
| 6 | phi_6 | 0.12382 | 0.16255 | 1.313 |
| 7 | phi_7 | 0.29700 | 0.32979 | 1.110 |
| 8 | phi_8 | 0.15820 | 0.18805 | 1.189 |
| 9 | phi_9 | 0.40829 | 0.42035 | 1.030 |
| 10 | phi_10 | 0.31624 | 0.33162 | 1.049 |
| 11 | beta_GenderMale | 0.05184 | 0.10929 | **2.108** |
| 12 | beta_CFUP | 0.03449 | 0.05016 | 1.454 |
| 16 | gamma_2 | 0.13623 | 0.13703 | 1.006 |
| 17 | gamma_3 | 0.06276 | 0.06327 | 1.008 |
| 19 | gamma_5 | 0.12466 | 0.12594 | 1.010 |
| 20 | gamma_6 | 0.07305 | 0.07498 | 1.026 |
| 21 | gamma_7 | 0.53060 | 0.56240 | 1.060 |
| 22 | gamma_8 | 0.12265 | 0.12307 | 1.003 |

Summary of ratios: min 1.003, 1st Qu. 1.027, median 1.116, mean 1.205, 3rd Qu. 1.288,
max 2.108.

**Recomputed 95% CI for the headline regression coefficients**, diagonal-only vs.
full-Hessian SE:
- `GenderMale`: `0.2181 ± 1.96*0.0518 = (0.117, 0.320)` -> `0.2181 ± 1.96*0.1093 =
  (0.004, 0.432)`. Still excludes 0, materially closer to it.
- `CFUP`: `-1.2706 ± 1.96*0.0345 = (-1.338, -1.203)` -> `-1.2706 ± 1.96*0.0502 =
  (-1.369, -1.172)`. Still clearly negative; smaller relative change.

### Multi-seed stability check

Script: `Examples/MultiSeedStabilityCheck.R` (new). Same `Examples/ReplicationCode.R`
settings, `full_hessian_se = TRUE`, refit from 5 different seeds (`1, 2, 42, 123,
2024`) run in parallel. Purpose: directly answer the editor's "how frequently does
this occur" and "role of random initialization" questions with evidence rather than
a single anecdotal run.

| seed | LL | NRun | beta_GenderMale | beta_CFUP | GradNorm (interior) | BoundaryAdjacent | SEratio_GenderMale | SEratio_CFUP |
|---|---|---|---|---|---|---|---|---|
| 1 | -2175.04880 | 31 | 0.21814 | -1.27062 | 0.0221 | [13,14,15,18,23,24] | 2.1082 | 1.4544 |
| 2 | -2175.04883 | 35 | 0.21817 | -1.27064 | 0.0083 | [13,14,15,18,23,24] | 2.1081 | 1.4544 |
| 42 | -2175.04883 | 30 | 0.21812 | -1.27064 | 0.0163 | [13,14,15,18,23,24] | 2.1082 | 1.4544 |
| 123 | -2175.04883 | 29 | 0.21813 | -1.27062 | 0.0241 | [13,14,15,18,23,24] | 2.1081 | 1.4543 |
| 2024 | -2175.04883 | 23 | 0.21823 | -1.27062 | 0.0066 | [13,14,15,18,23,24] | 2.1082 | 1.4544 |

**Conclusion**: across all 5 seeds, `Loglikelihood` agrees to within `3e-5`, `beta`
estimates agree to within `1e-4`, the interior gradient norm is small every time
(`0.007-0.024`), and the **same** set of 6 boundary-adjacent parameters (`mu1`, `nu`,
`gamma_1`, `gamma_4`, `gamma_9`, `gamma_10`) is flagged in *every single run*, with
essentially identical `SE_ratio` corrections. The number of coordinate sweeps needed
to converge varies (23-35, i.e. the *paths* differ), but the destination doesn't. For
this dataset, the boundary/near-flat behavior documented above is a robust,
reproducible feature of the likelihood surface given the data — not an artifact of a
particular random start — and we found no evidence of the optimizer landing on
materially different, competing optima.

**Caveat/scope**: this shows local stability around this basin across 5 restarts; it
is not a global search and cannot rule out a distant, materially different optimum
existing elsewhere in the 24-dimensional space. It answers "does the optimizer
reliably find the same thing from different starting points" (yes, here), not "is
this provably the global optimum."

### Correction to `gradient_check()` itself (2026-07-09)

While re-examining the boundary-adjacent gradient values above, an inconsistency
surfaced: `gamma_4`/`gamma_9` reported a *positive* gradient while `gamma_1`/`gamma_10`
were negative — but a genuine lower-bound optimum should show a consistently negative
one-directional score (moving away from the bound should make the fit worse). Traced
to the symmetric-shrunk-step formula comparing two points that straddle the optimum at
a resolution finer than Brent's own search tolerance (`tol_optimize`), which measures
curvature noise rather than a real slope. Fixed by using a one-sided difference with
the full `h_grad` step for any parameter where the symmetric step would otherwise be
shrunk below it. All four `gamma_k` are now consistently negative
(`gamma_1: -11.55`, `gamma_4: -1.81`, `gamma_9: -0.459`, `gamma_10: -0.00995` — see
`REVISION_NOTES.md` §9 for full detail and before/after numbers). This does not change
any conclusion above: `BoundaryAdjacent` classification, interior `GradientNorm`, and
every qualitative claim already relied only on sign/relative-magnitude, not the exact
old values.

## Decision: opt-in argument (implemented and verified)

Discussed three options: (a) keep `hessian_check()` standalone/opt-in, (b) always
compute it and add to `AdPaikModel()`'s default output, (c) replace
`StandardErrorParameters` outright wherever the interior submatrix is positive
definite. Chose (a), implemented as an argument rather than a fully separate
function: `full_hessian_se = FALSE` on `AdPaikModel()` (default unchanged).
Reasoning: the extra cost scales `O(n_p^2)` (`~32s` at `n_p=24` here) — fine for this
worked example, but this package's own stated selling point (abstract, Section 4.4)
is computational efficiency for large-scale datasets, and someone with many more time
intervals/regressors could have a much larger `n_p` where that cost stops being
trivial. An opt-in flag preserves default behavior/existing efficiency claims while
making the corrected inference a single-argument-away, well-documented option; it
does not require deciding this trade-off for every user.

**Implementation** (`R/AdPaikModel.R`, `R/summary_and_print.R`):
- New argument `full_hessian_se = FALSE`; when `TRUE`, `hessian_check()` is called
  once, right after the existing diagonal-only `params_se()`, using
  `GradientCheck$BoundaryAdjacent` and the diagonal SEs already computed.
- New `"HessianCheck"` element in the returned `AdPaik` object (`NULL` when
  `full_hessian_se = FALSE`), appended after `"GradientCheck"` (position 25 — safe
  with respect to `check.result()`'s positional validation, same reasoning as
  `GradientCheck` in `Rev_1.md`).
- `summary.AdPaik()`/`print.summary.AdPaik()`: when `HessianCheck` is present and the
  regressor block is fully interior (as it is here), show a
  "Corrected (full-Hessian, cross-derivative-aware) standard error" block alongside
  the existing regressor line.
- `print.AdPaik()`: shows the full `SE_full`/`SE_ratio` vectors, the interior
  positive-definiteness flag, and the `mu1`/`nu` correlation (or the honest
  "not available" message when the raw 2x2 block isn't negative definite).

**Verified end-to-end** (`Examples/ReplicationCode.R` settings, `set.seed(1)`):
- `full_hessian_se = FALSE` (default): `HessianCheck` is `NULL`, fit time unchanged
  (~469s), identical `Loglikelihood`/`NRun` to before this change.
- `full_hessian_se = TRUE`: reaches the exact same optimum as the default run
  (`all.equal(OptimalParameters)` TRUE), adds only ~32s, and reproduces the standalone
  `hessian_check()` numbers above exactly (`SE_ratio` for `GenderMale`/`CFUP` =
  `2.108228`/`1.454368`, matching to full precision). `summary()`/`print()` output
  confirmed correct (corrected regressor SEs shown, boundary-adjacent parameters
  correctly `NA`, `mu1`/`nu` correlation correctly reported as unavailable).

Not yet done: updating the actual paper text (Section 6.2's reported regressor SEs,
Section 4.4's justification discussion) to reflect this — left for the authors to
decide how to phrase, using the draft responses and numbers above.
