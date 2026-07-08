# TimeDepFrail — Revision Working Notes

Branch: `work` (created off `main`, `cb19433`). Purpose: build up a verified,
line-by-line understanding of what `R/AdPaikModel.R` and friends actually compute,
cross-checked against the JSS paper ("TimeDepFrail: Time-Dependent Shared Frailty Cox
Models in R"), before touching anything for the revision.

Rule for this file: only write down something as "confirmed" after tracing it against
the paper's equations/algorithms line by line. Anything not yet checked, or checked and
found questionable, goes under "Open questions" and is **not** acted upon until
confirmed by the author.

---

## 1. Model basis (from the paper — the ground truth)

- Frailty: `Z_jk = alpha_j + eps_jk`, for group `j`, interval `k`.
  - `alpha_j ~ Gamma(mu1/nu, 1/nu)`, iid across `j`.
  - `eps_jk ~ Gamma(mu2/gamma_k, 1/gamma_k)`, iid across `j` for fixed `k`.
  - Identifiability constraint: `mu1 + mu2 = 1` (so `mu2 = 1 - mu1` is not a free
    parameter), imposing `E[Z_jk] = 1` for all `k`.
- Hazard: `h_ijk(t) = Z_jk * exp(beta^T x_ij + phi_k)` for `t` in interval `I_k`.
- Parameter vector `p = [phi_1..phi_L, beta_1..beta_R, mu1, nu, gamma_1..gamma_L]`,
  length `n_p = 2L + R + 2`.
- `e_ijk` (paper Eq. 3): time exposure of unit `ij` inside interval `k` — 0 if the
  event/censoring time is before the interval, full interval width if strictly after,
  partial (`t_ij - a_{k-1}`) if inside it.
- `A_ijk = e_ijk * exp(beta^T x_ij + phi_k)`; `A_j.k = sum_i A_ijk`; `A_j.. = sum_{i,k} A_ijk`.
- Full group log-likelihood `ll_j` (paper Eq. 11 / Eq. 2 summed over groups): three
  additive pieces — (1) linear term in `phi`, `beta`, minus the `mu1/nu` log term with
  `A_j..`; (2) a sum over `k` of `mu2/gamma_k` log terms with `A_j.k`; (3) a sum over
  `k` of `log( sum_{l=0}^{d_j.k} choose(d_j.k, l) * [gamma ratios] )`.
- Posterior frailty (paper Eqs. 4-7, derived in Appendix C by analogy with the
  time-invariant shared-gamma-frailty empirical Bayes estimate):
  - `alpha_hat_j = (mu1_hat/nu_hat + N_j) / (1/nu_hat + H_j,.)`
  - `eps_hat_jk = (mu2_hat/gamma_hat_k + N_j(I_k)) / (1/gamma_hat_k + H_j,.(I_k))`
  - Both then divided by their respective max over `j` (or `j,k`) before summing into
    `Z_hat_jk`, so that the reported `Z_hat` has the right order of magnitude relative
    to `E[Z]=1`.
  - `H_j,.(I_k) = sum_i e_ijk * Y_ij(I_k) * exp(phi_k + beta^T x_ij)`, where `Y_ij(I_k)`
    is the at-risk indicator (paper doesn't give a closed-form for `Y_ij`, only the
    verbal definition: "1 if unit i in group j is at risk of facing the event in
    interval I_k").
- Optimization (paper Algorithm 3 / Section 4.3): Powell's method in `n_p` dimensions,
  implemented as repeated 1-D maximizations (R's `optimize()`, i.e. Brent's method)
  along one parameter at a time, cycling through all `n_p` directions, repeated for
  `n_total_run = n_p + n_extrarun` "runs", stopping early if the log-likelihood change
  between runs drops below `tol_ll`.
- Standard errors (paper Algorithm 5): diagonal-only Hessian via a centered
  finite-difference second derivative at the optimum, one parameter at a time,
  `se_p = 1/sqrt(-H_pp)`.
- Baseline hazard (paper Section 6.2, text only): "obtained by applying the exponential
  operation to `phi_hat`" — text does not mention any rescaling.
- Conditional survival (paper Eq. 9): `S_ij(t_ij | Z_jk) = exp( - exp(beta^T x_ij) *
  sum_k Z_jk * exp(phi_k) * Delta_I_k )` — uses the raw `exp(phi_k)`, not any rescaled
  version.

---

## 2. Code walkthrough vs. paper — confirmed matches

Traced line by line, these reproduce the paper's formulas exactly:

- `time_int_eval()` (`R/Utils.R`) == paper Eq. 3 (`e_ijk`), with `time_axis[k]` playing
  the role of `a_{k-1}` and `time_axis[k+1]` the role of `a_k` (1-indexed shift).
- `ll_AdPaik_centre_1D` / `ll_AdPaik_centre_eval` (`R/AdPaikModel.R`) == paper Eq. 11,
  term by term:
  - `loglik1` == the linear term + `-(mu1/nu) log(1 + nu * A_j..)` piece.
  - `loglik2` == `sum_k -(mu2/gamma_k) log(1 + gamma_k * A_j.k)`.
  - `loglik3`/`loglik4` == the `sum_k log(sum_l ...)` piece. The code computes
    `res6 = Gamma(mu1/nu+l)/Gamma(mu2/gamma_k)` and multiplies by
    `res_gamma3/res_gamma1 = Gamma(mu2/gamma_k+d_j.k-l)/Gamma(mu1/nu)` — at first glance
    this looked like the numerator/denominator of the two gamma ratios got swapped
    across each other, but the product is the same by commutativity:
    `(Gamma(mu1/nu+l)/Gamma(mu2/gamma_k)) * (Gamma(mu2/gamma_k+d_j.k-l)/Gamma(mu1/nu))`
    `== (Gamma(mu2/gamma_k+d_j.k-l)/Gamma(mu2/gamma_k)) * (Gamma(mu1/nu+l)/Gamma(mu1/nu))`,
    which is exactly the paper's Eq. 11 term. Verified algebraically — not a bug.
- The main `AdPaikModel()` optimization `while` loop matches paper Algorithm 3's
  structure: outer loop over runs, inner loop optimizing one parameter at a time via
  `optimize(..., maximum = TRUE)`, `params[index]` updated in place immediately (so
  later directions in the same run see the just-updated earlier parameters — this is
  standard coordinate-ascent / Powell-style behaviour and matches the paper's
  description), convergence via `actual_tol_ll <= tol_ll`, `Status` flag set to `FALSE`
  only if the loop runs out of `n_run` iterations.
- `n_run <- n_params + n_extrarun` == paper's `n_total-run = n_run + n_extra-run` where
  paper's own `n_run` symbol equals `sum(n_p)` (confusingly reused name in the paper
  text itself, not a code issue).
- `frailty_Sd_internal()` (`R/frailty_sd.R`) == paper Section 3.3:
  `var(Z_jk) = mu1*nu + mu2*gamma_k` (full) or `mu2*gamma_k` (time-varying only, `flag_full=FALSE`).
- `post_frailty_internal()` (`R/posterior_frailty.R`) == paper Eqs. 4-6 structurally:
  same numerator/denominator shape for `alpha_hat`/`eps_hat`, same division by the
  respective max, same variance formula (`var = estimate / (par2 + H)`, then divided by
  `max^2`).
- `survivalAdPaik()` (`R/survival.R`) uses `result$OptimalParameters[1:NIntervals]`
  directly (raw `phi`, not the rescaled `bas_hazard()` output) — consistent with paper
  Eq. 9, which uses raw `exp(phi_k)`.
- `extract_dummy_variables()` == paper's stated dummy-coding convention (first level
  alphabetically is the reference, `n-1` dummies for `n` levels).

---

## 3. Author decisions (2026-07-07) and resulting changes

Author reviewed Q1-Q5 below and gave explicit instructions. Summary of decisions and
what was actually changed, kept at the top since it supersedes the "open question"
framing for these five items:

- **Q1 (Y_risk)** — author: "fix this thing, in the way you think is best."
  Fixed in `R/posterior_frailty.R`, `extract_event_data()`: replaced the two-branch
  condition (which excluded the failure interval itself and relied on
  `max(time_to_event)` as a censoring-sentinel proxy) with the single condition
  `time_to_event[j] >= time_axis[k]` — a unit is at risk in interval `k` iff its
  event/censoring time has not occurred before the interval starts. This naturally
  includes the failure interval (partial exposure) and no longer depends on comparing
  against the dataset-wide max. Only affects posterior frailty
  estimates/variances (`post_frailty_est()`, `post_frailty_var()`, and downstream
  plots/CI) — does **not** touch the log-likelihood, the optimized `beta`/`phi`/`mu1`/
  `nu`/`gamma`, or the reported log-likelihood/AIC (those are built from a separate
  `dropout_matrix`/`e_matrix` pair in `AdPaikModel()`, untouched).

- **Q2 (params_se typo)** — author: "fix this typo."
  Fixed in `R/params_se_CI.R`, `params_se()`: corrected `values_plus_h` -> `value_plus_h`
  so the upper-bound clamp actually applies; changed the `if / else if` to two
  independent `if`s so both bounds can clamp at once for very narrow ranges; replaced
  the fixed-step `(ll_plus + ll_minus - 2*ll)/h_dd^2` formula with the general
  unequal-step centered second-derivative formula using the actual (possibly clamped)
  `h_plus`/`h_minus`, which reduces exactly to the original formula whenever neither
  side is clamped. Added a guard: if a parameter sits at (or past) its declared bound
  so that `h_plus <= 0` or `h_minus <= 0`, the second derivative can't be estimated and
  `se[p] <- 1e-4` (same sentinel already used for the pre-existing `Inf`/`-Inf` case).
  This changes standard errors for any parameter that was actually within `h_dd`
  (default `1e-3`) of a declared range bound — plausibly `nu` and/or `mu1` in the
  paper's own worked example (Section 6.5 already discusses `nu` pinning its upper
  bound). Not yet re-run against the published numbers — see §5.

- **Q3 (bas_hazard normalization)** and **Q4 (direction-order scheme)** — author will
  add explanatory text to the paper; short candidate sentences drafted in §3bis below.
  No code change (nothing was wrong; these were documentation gaps).

- **Q5 (undocumented `res6==0` safeguard)** — author: "do also the substitution."
  Read as: the existing `res6` guard (against `Gamma(mu1/nu+l)/Gamma(mu2/gamma_k)`
  underflowing to exactly `0`) only protected *one* of the two paired gamma-ratios in
  the same product; the other, `Gamma(mu2/gamma_k+d_j.k-l)/Gamma(mu1/nu)`, had no
  equivalent guard. Added a mirrored guard (`res7`, same `1e-10` substitution) in both
  `ll_AdPaik_centre_1D` and `ll_AdPaik_centre_eval` (`R/AdPaikModel.R`) — these two
  functions are textually identical copies of the same formula (one is used inside the
  1-D `optimize()` calls, the other to evaluate the log-likelihood at a full parameter
  vector), so the same edit was applied to both.

### Q3bis / Q4bis. Draft sentences for the paper

Not inserted anywhere — for you to place and adapt as needed.

**For Q3 (baseline hazard normalization), candidate addition to Section 6.2** (right
after "The piecewise linear estimated baseline hazard can be obtained by applying the
exponential operation to $\hat\phi$."):

> For visualization purposes, `bas_hazard()` rescales $\exp(\hat\phi_k)$ so that the
> resulting piecewise-constant function integrates to one over the time domain (i.e.
> $\sum_k \Delta_{I_k} \exp(\hat\phi_k) = 1$ after rescaling); the unnormalized values
> $\exp(\hat\phi_k)$ are used instead wherever the baseline hazard enters a model
> quantity, such as the conditional survival function of Eq. 9.

**For Q4 (direction-order scheme), candidate addition to Section 4.3**, after "...and
then repeats the same procedure but using another set of ordered directions.":

> Concretely, the first $\sum n_p$ runs use cyclic permutations of the natural
> parameter order (run $r$ starts from direction $r$ and wraps around through all
> $\sum n_p$ directions); any further runs, up to $n_{\text{extra-run}}$, instead use
> independently generated random permutations of the $\sum n_p$ directions.

---

## 4. Open questions — history (superseded by §3 for Q1-Q5, kept for the record)

### Q1. At-risk indicator `Y_risk` may zero out the failure interval itself
`extract_event_data()` in `R/posterior_frailty.R` (used only for posterior frailty
estimates, i.e. `alpha_hat`, `eps_hat`, `Z_hat` and their variances — **not** used by the
log-likelihood optimization itself, which builds its own `dropout_matrix`/`e_matrix`
directly in `AdPaikModel()`):

```r
if((time_to_event[j] < max(time_to_event)) & (time_to_event[j] > time_axis[k+1]))
  Y_risk[j,k] <- 1
else if(time_to_event[j] == max(time_to_event))
  Y_risk[j,k] <- 1
```

For a unit whose event happens strictly inside interval `k` (i.e. `time_axis[k] <=
time_to_event[j] < time_axis[k+1]`), neither branch fires, so `Y_risk[j,k] = 0` for
that interval — even though `e_ijk[j,k] = time_to_event[j] - time_axis[k] > 0` (the
partial exposure time before the event). Since
`cum_hazard_jk[j,k] <- e_ijk[j,k] * Y_risk[j,k] * exp(...)`, this multiplies that
partial exposure by zero, dropping the unit's contribution to `H_j,.(I_k)` for the
very interval in which they fail.

Standard survival-analysis convention (and my reading of paper Eq. 17 / Algorithm 6
step 7, `H_ijk = e_ijk * Y_ijk * exp(phi_k + x_ij*beta)`) is that a unit *is* at risk
for the partial time leading up to its own event within that interval — that's the
whole reason `e_ijk` captures a partial width instead of just 0/full. If so, `Y_risk`
should be 1 in the failure interval too, and the current code would be
under-counting `H_j,.(I_k)` (and therefore inflating `eps_hat_jk`/`alpha_hat_j`, since
they're `(.../ (par2 + H))`) specifically for group/interval combinations with more
events.

**Question for you:** is `Y_risk = 0` in the failure interval intentional (a specific
modeling choice I'm not aware of), or should it be 1? This only affects posterior
frailty estimates/variances, not the estimated `beta`/`phi`/`mu1`/`nu`/`gamma`
themselves or the reported log-likelihood/AIC.

### Q2. `params_se()` range-clamp looks like it has a typo
`R/params_se_CI.R`, inside the per-parameter loop:

```r
if(value_plus_h > params_range_max[p])
  values_plus_h <- params_range_max[p]     # <- assigns to "values_plus_h" (with s),
                                             #    a variable never read again
else if(value_minus_h < params_range_min[p])
  value_minus_h <- params_range_min[p]
```

The first branch assigns to `values_plus_h` (plural) instead of `value_plus_h`
(singular) — the actual `value_plus_h` used two lines later in
`params_plus[p] <- value_plus_h` is never clamped to the upper range bound. Also,
structurally, this is an `if / else if`, so even if the typo were fixed, a parameter
that is simultaneously within `h_dd` of *both* bounds (possible only for a very narrow
range) would only get one side clamped, and the finite-difference step size would
become asymmetric — but the denominator `h_dd * h_dd` in the second-derivative formula
a few lines below is not adjusted for that asymmetry.

**Question for you:** do you want this looked at more closely (i.e. does it actually
change any reported standard error for the paper's worked example), or is it out of
scope for this revision? I have not checked whether any of the paper's example
parameters actually sit within `h_dd` (default `1e-3`) of a declared bound — that's a
one-line check I can do on request, but I'm not running it unprompted.

### Q3. `bas_hazard()` rescales `exp(phi)`; paper text doesn't mention it
`bas_hazard_internal()` (`R/bas_hazard.R`) divides `exp(phi)` by the total area under
the piecewise-constant curve (`sum_k (time_axis[k+1]-time_axis[k]) * exp(phi_k)`)
before returning it. I checked this numerically against the paper's own printed
example output (Section 6.2, `bas_hazard(result)` -> the 10 numbers starting
`0.248..., 0.157..., 0.661...`): multiplying each by its interval width and summing
gives `1.0000` to 4 decimal places — so **the code and the paper's worked example are
numerically consistent**, but the paper's prose ("piecewise linear estimated baseline
hazard can be obtained by applying the exponential operation to phi_hat") never
mentions this area-normalization step. `survivalAdPaik()` correctly uses the raw
(non-normalized) `exp(phi)` per Eq. 9, so the two functions deliberately use two
different scalings of the same `phi`.

**Question for you:** is this intentional (i.e. `bas_hazard()` is meant to report a
*shape*, normalized like a density, while the survival function needs the actual
unnormalized hazard), and if so, should the paper say so explicitly? This is a
documentation gap, not a code bug as far as I can tell, but I want to confirm before
suggesting any text for the revision.

### Q4. Direction-order scheme in the Powell loop is an interpretation, not a verified match
`AdPaikModel()` builds `RunIndexes` so that the first `n_p` runs are cyclic shifts of
`[1..n_p]` (starting position `i` for run `i`), and any further runs (the
`n_extrarun` of them) use independent random permutations. The paper only says (Section
4.3) that the procedure "repeats ... using another set of ordered directions" without
specifying how that set is generated. The code's scheme is a reasonable reading, but I
have not found anything in the paper that pins it down exactly.

**Question for you:** was this cyclic-shift-then-random scheme a deliberate design
choice, and is it worth spelling out explicitly in the paper (Section 4.3 currently
leaves it implicit)? Not proposing a change, just flagging that I can't verify it
against the text beyond "plausible reading."

### Q5. Undocumented numerical safeguard in the log-likelihood
Inside `ll_AdPaik_centre_1D`/`ll_AdPaik_centre_eval`:
```r
res6 <- res_gamma4/res_gamma2
if(res6 == 0)
  res6 <- 1e-10
```
This silently substitutes `1e-10` whenever the ratio of two gamma functions evaluates
to exactly `0` (floating-point underflow), to avoid `log(0) = -Inf` propagating out of
`loglik3`. Not in the paper. Likely fine as a stability patch, but I have not checked
how often it actually triggers on the `data_dropout` example (i.e. whether it's a
no-op in practice or is silently altering the likelihood surface near the boundaries
discussed in Section 6.5, where `nu` and other parameters push against their range
limits). Flagging, not fixing.

---

## 5. Remaining files — reviewed, nothing flagged

- `R/plot.R` (`plot_bas_hazard`, `plot_post_frailty_est`, `plot_post_frailty_var`,
  `plot_frailty_sd`, `plot_ll_1D`): pure visualization, reads already-computed
  quantities off the `AdPaik`/vector inputs, no independent computation of anything
  that feeds into the optimization or the posterior frailty numbers. Nothing to flag.
- `R/check.R`, `R/check.result.R`: input/output structural validation only (types,
  lengths, ranges). Nothing that touches the log-likelihood, optimization, or
  posterior-frailty formulas. Nothing to flag.
- `R/coefse.R`, `R/extractors.AdPaik.R`, `R/extract_dummy_variables.R`,
  `R/summary_and_print.R`: thin wrappers/extractors over already-computed
  `AdPaikModel()` output, or dummy-coding utilities. Match the paper's stated
  conventions (Section 6.1 dummy coding; AIC formula in Eq. 12 reproduced exactly in
  `extractAIC.AdPaik`/`AdPaikModel`'s own `AIC <- 2*n_params - 2*optimal_loglikelihood`).

This completes a full pass over every file in `R/` (17 files). Everything not listed
as a discrepancy in §3/§4 above is a confirmed match to the paper.

Note (not a flagged issue, just an observation, out of scope of Q1-Q5): `N_i` in
`extract_event_data()` (`R/posterior_frailty.R`, used for the time-independent
`alpha_hat_j`) still uses `time_to_event[index] < max(time_to_event)` as its
"did this individual actually fail" test — the same reliance on the dataset-wide max
as a censoring-sentinel proxy that Q1's `Y_risk` fix removed. It's a different
variable (event count, not at-risk indicator) and wasn't part of what was authorized
for fixing, so left untouched — flagging only so it doesn't get lost.

- Whether `optimize()`'s Brent search, given only an interval and no starting value,
  ever fails to find the true 1-D maximum in a way that matters (paper just cites
  `optimize`/`optim` with `"Brent"`, code uses `optimize()` directly — consistent, not
  flagged as an open question, just noting it's the literal R function the paper
  names).

## 6. Post-fix sanity run — results

Ran the paper's own worked example (`data_dropout`, same `formula`/`time_axis`/
`categories_range_min/max`, `set.seed(1)`, reduced `n_extrarun = 10` instead of the
default 60 for speed) against the fixed code.

**Unchanged vs. the published paper (as expected — these fixes don't touch the
log-likelihood used for optimization):**
- `Loglikelihood = -2175.135`, `AIC = 4398.27` — exact match to the paper.
- `coef()`: `GenderMale = 0.217802`, `CFUP = -1.270657` — exact match to the paper.
- `Status = TRUE`, converged in 31 runs.

**Changed by the fixes (this is the actual, quantified effect of Q1/Q2/Q5):**
- Standard error of `nu` (parameter index 14 = `L+R+2` with `L=10, R=2`) is now
  **24.236**, vs. the paper's published **218.643** (Section 7's discussion of `nu`'s
  independence from the log-likelihood is built on that 218.643 figure — this
  directly changes the numeric evidence behind that paragraph, though the qualitative
  point — `nu` pinned near its upper range bound, poorly identified — still stands:
  24.236 is still very large relative to `nu`'s point estimate near 1).
- Four `gamma_k` (parameter indices 15, 18, 23, 24, i.e. `gamma_1, gamma_4, gamma_9,
  gamma_10`) get the `se = 1e-4` boundary sentinel. This is the *same outcome* as
  before the fix (they hit it via the pre-existing `Inf`-hessian catch previously;
  now they hit it via the explicit "no room for a two-sided step" guard) — not a
  new finding, just re-derived through the corrected code path.
- `mu1`'s SE is 0.0511 (index 13) — did not hit the boundary sentinel, unremarkable.
- Posterior `Z` range: `[0.687, 1.921]`; `alpha` range `[0.355, 1]`; `eps` range
  `[0.332, 1]` (both properly bounded above by 1 since each is divided by its own
  max) — all finite, non-negative. No quantitative comparison to the paper's exact
  per-group numbers was done (would need the same `n_extrarun=60`/full run and the
  paper's own random seed, which isn't published), but the shape/range matches
  Figure 4's visual range (~0.8-1.8) reasonably closely.
- `bas_hazard()` still integrates to exactly `1` over the time domain — confirms the
  Q3 normalization is unaffected by the other fixes (no regression).
- No `NA`/`Inf`/negative values anywhere in SEs or posterior frailty estimates/variances.

**Caveat:** this was a single run with a reduced `n_extrarun` (10 instead of 60) for
speed, not the paper's exact reproduction settings — good enough to confirm no
crashes/regressions and to get an order-of-magnitude read on the `nu` SE change, but
if the exact new SE numbers are going into the paper's revision text, it should be
re-run with `n_extrarun = 60` (the published default) first.

---

## 7. Editor/referee comment: first-order optimality diagnostic (2026-07-07/08)

Full design discussion and rationale in `Rev_1.md`. Summary of what was actually done,
across three commits' worth of work on branch `revision` (created off `work` at
`12543d8`):

- **New post-hoc gradient diagnostic** (`R/gradient_check.R`, wired into
  `AdPaikModel()` via a new `h_grad = 1e-4` argument): a finite-difference gradient of
  the log-likelihood evaluated once at the reported optimum, stored as
  `result$GradientCheck` and shown by default in `print()`/`summary()`. Chosen over
  porting `frailtypack`'s in-loop analytic-gradient threshold, because that only works
  for a gradient-based optimizer — `AdPaikModel()`'s Powell/Brent search is
  derivative-free by construction, so a post-hoc check (exactly what the referee's own
  wording proposed) is the lower-risk option: read-only, cannot perturb the existing
  optimization or previously reported numbers.
  - Step is always symmetric, shrunk to `min(h_grad, distance to nearer declared
    bound)` — an asymmetric clamped step (full step one side, boundary-clamped tiny
    step the other) was tried first and produces numerically unstable, inflated
    gradient values; the symmetric-shrink fixes this cleanly.
  - Parameters whose optimum sits within `10*h_grad` of a declared range bound are
    reported separately (`GradientCheck$BoundaryAdjacent`), excluded from the headline
    `GradientNorm`/`GradientMaxAbs`. These commonly show large, genuine one-sided
    gradients at a boundary solution (seen for `mu1`, `nu`, and several `gamma_k` in
    the worked example below) — a different, already-documented phenomenon (`nu`
    pinning its upper bound is discussed in the paper's Section 7) from an interior
    optimizer stall, which is what the referee is actually asking to rule out.

- **Found and fixed a real numerical fragility while validating the diagnostic**: the
  very first test run produced `Inf`/`NaN` gradient components. Root cause:
  `ll_AdPaik_centre_1D`/`ll_AdPaik_centre_eval` computed the third line of the
  log-likelihood (paper Eq. 11's `sum_k log(sum_l ...)` term) via raw `gamma()` ratios
  (`res_gamma1..4`, `res6`, `res7`). For a `gamma_k` shrinking towards its lower range
  bound, `mu2/gamma_k` grows past ~171, where `gamma()` overflows to `Inf` in double
  precision — a small perturbation (as small as `h_grad = 1e-4`) is then enough to push
  the evaluation over that edge. Rewrote both functions' third-line computation in
  log-space (`lgamma`/`lchoose` + log-sum-exp, `R/AdPaikModel.R`), which is
  mathematically identical (verified by direct comparison against the old formula:
  matches to ~1e-14 away from the overflow region, matches exactly at the previously
  overflowing point, and stays finite under the exact perturbation that used to blow up
  to `Inf`) but has no overflow edge. This also makes the old ad hoc `res6==0`/`res7==0`
  → `1e-10` substitution (Q5 above) unnecessary; it was removed.

- **Found and fixed an unrelated pre-existing bug surfaced while investigating an
  unexpected result**: `AdPaikModel()`'s own roxygen `@examples` block used
  `categories_range_max <- c(-eps, 0, 1 - eps, 1, 10)` (beta category upper bound `0`),
  but `Examples/ReplicationCode.R` — the actual script that reproduces the paper's
  numbers — uses `0.5`. With the upper bound clamped to `0`, `GenderMale`'s true
  coefficient (`0.2178` in the paper) is infeasible and the optimizer instead lands
  exactly on the boundary (`~-6.6e-7`, prints as `0`). Confirmed by running the
  **unmodified, pre-session code** (commit `12543d8`) with both bounds: `max=0` gives
  the wrong, boundary-clamped `GenderMale≈0`/`LL=-2177.212`; `max=0.5` exactly
  reproduces the paper's `GenderMale=0.217802`/`LL=-2175.135`/`NRun=31`. So this was a
  pre-existing documentation bug, unrelated to (and predating) this session's other
  changes — fixed by changing the example's bound to `0.5` to match
  `Examples/ReplicationCode.R`.

- **Final verification** (`data_dropout`, same formula/`time_axis`, `set.seed(1)`,
  `n_extrarun=10`, corrected example bounds, all fixes above applied together):
  `Loglikelihood = -2175.049`, `AIC = 4398.098`, `coef = (GenderMale=0.2181,
  CFUP=-1.2706)`, converged in `31` runs — matches the paper to 2-3 decimal places (the
  small residual difference is the expected effect of the mathematically-equivalent but
  not bit-identical log-space reformulation, propagated through 31 rounds of
  coordinate-wise search). No `NA`/`Inf`/`NaN` anywhere in the new gradient diagnostic.
  Interior gradient: `‖grad‖₂ = 0.0221`, `max|grad| = 0.0162` (18 of 24 parameters) —
  small, consistent with a genuine interior stationary point. Boundary-adjacent:
  `mu1`, `nu`, `gamma_1`, `gamma_4`, `gamma_9`, `gamma_10` (6 of 24), with large
  one-sided gradients as expected for boundary solutions.

### 7bis. Why those 6 parameters, specifically — investigation (2026-07-08)

Prompted by the author asking whether the 6 boundary-adjacent parameters above are an
error. Not a code bug — the diagnostic is correctly reporting genuine boundary
solutions found by the optimizer — but the mechanism is more nuanced than a single
clean story, and one caveat about the test settings matters before treating this as
settled:

- **`mu1 -> 6e-7` (lower bound).** `alpha_j ~ Gamma(mu1/nu, 1/nu)`: as shape
  `mu1/nu -> 0`, this degenerates to a point mass at `0`. The fit is saying the
  group-constant frailty component `alpha_j` vanishes; essentially all frailty
  heterogeneity is time-varying (`eps_jk`) in this run. This specific finding — `mu1`
  itself collapsing, not just `nu` — does not appear to have been flagged before
  anywhere else in this file.
- **`nu -> 0.99999995` (upper bound).** Once `mu1 ~= 0`, the shape `mu1/nu ~= 0`
  regardless of `nu`'s value, so `nu` becomes only weakly identified and likely just
  rides wherever the coordinate sweep last left it. This *is* the same phenomenon the
  paper's Section 7 already discusses for `nu` (previously quantified via its large SE,
  `24.236`/`218.643` depending on the params_se fix — see §3bis above); the gradient
  diagnostic now shows it directly as a large one-sided slope instead of only via SE.
- **4 of 10 `gamma_k -> 0` (`gamma_1, gamma_4, gamma_9, gamma_10`).** Checked directly
  against `data_dropout` rather than assumed:

  | interval | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
  |---|---|---|---|---|---|---|---|---|---|---|
  | total events | 76 | 47 | 233 | 24 | 53 | 110 | 13 | 79 | 6 | 10 |
  | centres (of 16) with 0 events | 0 | 2 | 0 | 5 | 3 | 0 | 10 | 1 | 11 | 9 |
  | `gamma_k` at boundary? | yes | no | no | yes | no | no | no | no | yes | yes |

  Intervals 9/10 are cleanly explained: almost no centre has any event there, so
  there is no information to estimate a between-centre variance component, and
  `gamma_k -> 0` is the expected MLE behaviour for a variance parameter with no
  evidence of heterogeneity. But this is **not** a complete explanation: interval 7 is
  just as sparse (10/16 centres with zero events) yet its `gamma_7` did *not* collapse
  (landed at `0.137`), while interval 1 has substantial data (76 events, every centre
  represented) yet *did* collapse. So sparsity alone does not fully explain the
  pattern — plausibly a mix of genuine data-sparsity-driven shrinkage (intervals 9/10,
  maybe 4) and coordinate-wise-optimization path-dependence in a non-concave
  24-dimensional space (interval 1, possibly 4) — i.e. partially the exact fragility
  the referee's comment is about.

- **Caveat on these specific numbers:** everything above (this run, and the "Final
  verification" bullet before it) used `n_extrarun = 10`, reduced from the paper's own
  default `n_extrarun = 60`, purely for iteration speed while building the diagnostic.
  Which indices land on a boundary is plausibly sensitive to `n_extrarun`/random seed.
  Treat "`mu1`/`nu`/4 specific `gamma_k` are boundary-pinned" as provisional pending a
  re-run at `n_extrarun = 60` (and ideally more than one seed) before it goes into the
  response letter or the paper as a stated finding.

### 7ter. Future work (deferred): optimizer redesign — block/joint updates or adaptive directions

Author asked (2026-07-08) whether to move from one-parameter-at-a-time coordinate
updates to block-wise or fully joint optimization (grouping `gamma_k`, `beta`, `phi`
together, leaving `mu1`/`nu` as-is), and separately whether to compute cross-derivatives
(off-diagonal Hessian) for proper SEs. **Decision: record the analysis below, do not
design or implement yet** — this would be a change to the optimizer/inference itself
(unlike today's diagnostic-only work), with much higher blast radius on previously
reported numbers, and needs to be scoped as its own dedicated piece of work.

**Structural finding, derived directly from Eq. 11/Eq. 2 of the paper (not assumed):**
- The 10 `gamma_k` are **provably separable** from each other: for fixed `phi`,
  `beta`, `mu1`, `nu`, both `loglik2` and `loglik3`'s inner term reference only
  `A_j.k`/`d_j.k` for that same `k` — there is no `gamma_k`/`gamma_k'` cross term
  anywhere in the log-likelihood. So a joint/block update of all `gamma_k` together
  converges to the *same point* as updating them one at a time in sequence; it buys no
  correctness/convergence-quality benefit, only possibly some wall-clock speed from
  batching evaluations. **Recommendation: don't block `gamma_k` for correctness
  reasons** — the pattern found in §7bis (`gamma_1`/`gamma_4`/`gamma_9`/`gamma_10`
  collapsing) is not something block-updating them would change.
- `beta` and `phi` **are** genuinely coupled to each other (and to `mu1`/`nu`) through
  the shared term `-(mu1/nu)*log(1+nu*A_j..)` in `loglik1`, since
  `A_j.. = sum_{i,k} A_ijk` mixes every `phi_k` and every `beta_r` together. This is
  exactly the kind of off-axis correlation where one-at-a-time coordinate ascent can
  zig-zag inefficiently. `beta` (only 2 params) is cheap/low-risk to block-optimize
  jointly; `phi` (`L` params) is the bigger, more consequential candidate.
- **`mu1`/`nu` — pushed back on excluding these.** §7bis's own finding (`mu1 -> 0`
  makes `nu` nearly unidentifiable, since `mu1/nu ~= 0` regardless of `nu`'s value
  once `mu1` collapses) is a textbook off-axis ridge — precisely the pathology
  one-at-a-time Brent search handles worst, and arguably the pair most in need of a
  joint step, not least.

**Three options discussed, not yet chosen between:**
1. **Targeted blocking** — keep the Powell-cycle/Brent scaffolding, replace the 1-D
   update for `beta`, `phi`, and (recommended, contrary to the original suggestion)
   `mu1`+`nu` with small joint box-constrained optimizations (e.g.
   `optim(method="L-BFGS-B")`) using the finite-difference gradient built today
   (`R/gradient_check.R`); leave `gamma_k` as individual 1-D updates (per the
   separability finding above). Moderate scope, keeps most of the existing structure.
2. **True adaptive-direction Powell** (Press 2007, as actually cited by the paper): the
   current code cycles through a *fixed* set of directions (the `n_p` coordinate axes,
   reordered between sweeps) — it is coordinate-descent-with-Brent, not full Powell's
   method as classically defined. Classical Powell additionally replaces one direction
   per sweep with the net direction moved over that whole sweep, letting the direction
   set adapt toward the objective's actual (possibly non-axis-aligned) valleys/ridges
   over iterations — directly targeting exactly the `mu1`/`nu` ridge, while staying
   derivative-free and matching the paper's own citation. Smaller conceptual change
   than switching optimizer family, but a real implementation (needs an acceptance
   test before swapping in the new direction, else the direction set can lose linear
   independence over many sweeps — a known pitfall in naive implementations).
3. **Full multivariate solve** (`optim(method="L-BFGS-B")` over all `n_p` parameters at
   once) — replaces the coordinate-descent scaffolding entirely, handles all
   cross-coupling uniformly without hand-picking blocks, and yields a full numerical
   Hessian (via a follow-up `optimHess()`/finite-difference pass) for proper
   non-diagonal-only SEs as a natural byproduct — answering the separate
   cross-derivative-SE question in the same stroke. Largest rewrite, highest risk to
   previously-published numbers (would need the same rigorous
   old-vs-new-formula-equivalence verification done for the log-space log-likelihood
   fix in §7, but for the optimizer itself, which is harder to verify since it changes
   *which* optimum is found, not just how a fixed quantity is computed).

No design/implementation performed on any of these three — deferred to a future
session at the author's request.

---

## 8. Editor/referee comments: near-flat likelihood regions & cross-derivative SEs (2026-07-08)

Full reasoning, draft response-letter text, and real numbers in `Rev_2_3.md`; code
summary here for the technical record.

- **`R/hessian_check.R`** (new, internal): full (non-diagonal) finite-difference
  Hessian at a given parameter vector, using the same symmetric boundary-aware step
  as `gradient_check()`. Off-diagonal terms via the standard 4-point mixed-partial
  formula. Inverts only the submatrix of non-boundary-adjacent ("interior")
  parameters for the covariance/SE computation (a single near-degenerate boundary
  row/column can corrupt a full-matrix inverse even for otherwise well-behaved
  parameters, unlike the diagonal-only approximation); the `mu1`/`nu` raw 2x2
  sub-block is checked and reported separately regardless of the interior/boundary
  split, since Section 3.3's `mu1`/`nu` algebraic coupling makes it of particular
  interest.
- **`R/AdPaikModel.R`**: new opt-in argument `full_hessian_se = FALSE`. When `TRUE`,
  calls `hessian_check()` once after the existing diagonal-only `params_se()`; adds
  `"HessianCheck"` to the returned object (`NULL` when `FALSE`; appended after
  `"GradientCheck"`, position 25, safe per `check.result()`'s positional check — see
  §7). Kept opt-in rather than default or a fully separate function: cost scales
  `O(n_p^2)` (`~32s` extra at this worked example's `n_p=24`), which is fine here but
  could matter for the larger `n_p` this package's own efficiency claims target;
  default behavior/cost stays exactly as before.
- **`R/summary_and_print.R`**: `summary()`/`print()` show the full-Hessian diagnostic
  when present (`HessianCheck != NULL`), including a corrected regressor SE line in
  `summary()`.
- **Verified against the paper's worked example** (`n_extrarun=60`,
  `Examples/ReplicationCode.R` settings): `full_hessian_se=TRUE` reaches the identical
  optimum as `full_hessian_se=FALSE` (no effect on the fit itself), and the resulting
  interior covariance matrix (18 of 24 parameters) is positive definite. The
  cross-derivative correction is material and parameter-type-specific:
  `SE_full/SE_diag` ratios range `1.00-2.11`; `phi_k` +10-44%, `gamma_k` (non-boundary)
  +0.3-6% (negligible, consistent with their provable pairwise separability, §7ter),
  `beta_CFUP` +45%, **`beta_GenderMale` +111%** (SE `0.0518 -> 0.1093`). Recomputed 95%
  CI for `GenderMale` widens from `(0.117, 0.320)` to `(0.004, 0.432)` — still excludes
  0, materially closer to it than the diagonal-only approximation suggested. The raw
  `mu1`/`nu` 2x2 Hessian sub-block is not negative definite (no valid correlation/SE
  for that pair from any Hessian-based method, not just the diagonal shortcut) —
  independent corroboration of the boundary/ridge finding in §7bis.
- **On an "independence" justification for the diagonal shortcut** (author question):
  checked directly against the paper text. The only independence statements are
  model-level (`alpha_j`/`eps_jk` as independent random variables, Section 3.3;
  conditional independence of units given the frailty). Section 4.4's own stated
  justification for the diagonal-only Hessian is purely computational, not an
  independence argument. Model-level independence of random effects does not imply
  independence of their *estimators*' sampling error — today's own numbers show `beta`
  and `phi` are correlated (shared `A_j..` term) and the `mu1`/`nu` block is not even
  well-behaved, contradicting any such assumption if it was implicitly relied upon.

- **Multi-seed stability check** (`Examples/MultiSeedStabilityCheck.R`, new): the
  editor/referee's near-flat-regions comment explicitly asks "how frequently" this
  occurs and about the role of random initialization — a single seed can't answer
  that. Refit the worked example (`full_hessian_se=TRUE`) from 5 different seeds
  (`1, 2, 42, 123, 2024`), run in parallel. Result: `Loglikelihood` agrees to within
  `3e-5` across all 5, `beta` estimates agree to within `1e-4`, and the **same** set
  of 6 boundary-adjacent parameters (`mu1`, `nu`, `gamma_1`, `gamma_4`, `gamma_9`,
  `gamma_10`) is flagged in every single run, with near-identical `SE_ratio`
  corrections — despite the number of coordinate sweeps to converge varying (23-35).
  For this dataset, the boundary/near-flat behavior is a robust, reproducible feature
  of the likelihood surface given the data, not a random-seed artifact; no evidence
  of the optimizer landing on materially different, competing optima across these 5
  restarts (though this cannot rule out a distant optimum elsewhere in the
  24-dimensional space — a local-stability check, not a global search). Full table
  in `Rev_2_3.md`.
