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

## 3. Open questions — need author confirmation before any code change

These are things I found by direct comparison of code vs. paper where either the code
looks inconsistent with the stated formula, or the paper doesn't specify enough to
tell. **Not fixing anything until you weigh in.**

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

## 4. Not yet reviewed

- `R/plot.R`, `R/check.result.R` — plotting/validation-only, lower priority for "is the
  optimization correct" but will look if relevant.
- Whether `optimize()`'s Brent search, given only an interval and no starting value,
  ever fails to find the true 1-D maximum in a way that matters (paper just cites
  `optimize`/`optim` with `"Brent"`, code uses `optimize()` directly — consistent, not
  flagged as an open question, just noting it's the literal R function the paper
  names).
