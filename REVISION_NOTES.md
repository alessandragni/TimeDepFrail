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
