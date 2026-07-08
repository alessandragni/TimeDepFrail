# Revision 1 — Editor/Referee comment: convergence diagnostics for the Powell/Brent optimizer

## The comment being addressed

**Editor:** The package relies on a coordinate-wise, derivative-free optimization
strategy based on Powell's method combined with Brent's one-dimensional searches. While
this approach is clearly described and computationally convenient, the resulting
log-likelihood is high-dimensional and not guaranteed to be globally concave.
Convergence is currently assessed primarily via changes in the log-likelihood value. The
referee also raises this point and suggests that additional convergence diagnostics,
such as numerical checks of first-order optimality conditions (e.g., finite-difference
gradients evaluated at the reported optimum), would provide stronger reassurance that
the reported solution corresponds to a genuine stationary point rather than a flat
region of the likelihood surface.

**Referee:** ...it would be helpful to report a numerical check of the first-order
optimality conditions (similar to that provided in `frailtypack`). For example, could
the authors evaluate the finite-difference gradient of the log-likelihood at the
reported optimum and report its norm or maximum absolute component? Such diagnostics
would provide reassurance that the solution corresponds to a genuine stationary point
rather than a plateau of the likelihood surface.

## Chosen strategy: post-hoc finite-difference check, not an in-loop analytic-gradient threshold

`frailtypack::frailtyPenal` uses a Marquardt (gradient/Hessian-based) optimizer with an
internal `LIMderiv` threshold on the *analytic* gradient as one of its own stopping
criteria — that's only available because that optimizer already computes derivatives
every iteration. `AdPaikModel()`'s Powell-coordinate/Brent-line-search optimizer is
derivative-free by construction; there is no analytic gradient to threshold without a
full re-derivation of closed-form partials (digamma terms for `mu1`, `nu`, `gamma_k`,
plus the `phi`/`beta` linear terms) and a change to the optimizer's own stopping rule —
which risks shifting the reported optimum and every downstream number in the paper's
worked example, a large risk for a journal-revision deadline.

The referee's own wording already asks for the cheaper, safer alternative: evaluate a
finite-difference gradient **at the reported optimum, after optimization finishes**, and
report its norm / max absolute component. This is read-only and additive — it cannot
perturb the existing optimization or any previously reported number, consistent with the
additive-fix philosophy already used for the Q1/Q2/Q5 items in `REVISION_NOTES.md`.

**Decision: implement the post-hoc finite-difference diagnostic.**

## Design (confirmed with author)

- Diagnostic is computed automatically inside `AdPaikModel()`, stored in the returned
  `"AdPaik"` object, and shown by default in `print.AdPaik()` / `summary.AdPaik()` output
  (same treatment `StandardErrorParameters` already gets) — no extra call needed by the
  user.
- New argument `h_grad = 1e-4` on `AdPaikModel()`, distinct from the existing `h_dd`
  (used only for the Hessian/SE second-derivative approximation, default `1e-3`).

## Implementation

**1. New internal function `gradient_check()` in a new file `R/gradient_check.R`.**

Superseded during implementation (see "What actually changed vs. this plan" below):
the first version mirrored `params_se()`'s per-side clamp (`value_plus_h`/
`value_minus_h` clamped independently, then a general unequal-step formula). Testing
showed this produces numerically unstable, inflated gradients whenever a parameter
sits very close to (but not exactly at) a declared bound, because the two sides end up
with wildly different step sizes. Replaced with a **symmetric shrunk step**:
`h_p <- min(h_grad, value - range_min, range_max - value)`, applied identically on
both sides, plus a **boundary-adjacency flag** so parameters close to a bound are
reported separately from the interior gradient norm/max. Final behaviour:

- `h_p <- min(h_grad, value - params_range_min[p], params_range_max[p] - value)`.
- If `h_p <= 0` (sits at/past a bound): `gradient[p] <- NA`.
- Otherwise centered difference `(ll_plus - ll_minus)/(2*h_p)`, falling back to a
  one-sided difference if perturbing pushes the log-likelihood itself non-finite on
  one side (see the overflow finding below), or `NA` if non-finite on both sides.
- `is_boundary <- dist_to_nearest_bound < 10 * h_grad`; `BoundaryAdjacent` reports
  those indices; `GradientNorm`/`GradientMaxAbs` are computed over the non-boundary,
  finite components only.
- Returns `list(Gradient, BoundaryAdjacent, GradientNorm, GradientMaxAbs)`.

**2. Wire into `AdPaikModel()` (`R/AdPaikModel.R`)**:
- Add `h_grad = 1e-4` to the function signature (near `h_dd`), with a `@param` roxygen
  entry.
- Right after `optimal_params`/`optimal_loglikelihood` are extracted
  (`AdPaikModel.R:344-345`), call:
  `gradient_check_result <- gradient_check(optimal_params, params_range_min, params_range_max, dataset, centre, time_axis, dropout_matrix, e_matrix, h_grad)`
- Append `"GradientCheck" = gradient_check_result` as the **last** element of
  `return_list` (after `"PosteriorFrailtyVariance"`, i.e. new position 24). Important:
  `check.result()` (`R/check.result.R:64-86`) only positionally validates names for
  `i in 1:length(names_list.AdPaik)` (23 names) — appending after that keeps the
  existing check unmodified. `names_list.AdPaik` itself is left untouched.
- Add a `if (verbose) message(...)` line reporting the gradient norm/max-abs, consistent
  with the existing verbose messages for SE/CI/etc.

**3. Update `R/summary_and_print.R`**:
- `summary.AdPaik()`: add `gradient_norm`/`gradient_max_abs` (rounded) to `summary_list`.
- `print.summary.AdPaik()`: add a line near the convergence status, e.g.
  `"First-order optimality check: max|grad| = <x>, ||grad||_2 = <y>"`.
- `print.AdPaik()`: add an analogous block after the Standard Errors / CI printing,
  printing `x$GradientCheck$GradientNorm` and `x$GradientCheck$GradientMaxAbs` (and
  optionally the full gradient vector, matching how `OptimalParameters`/
  `StandardErrorParameters` are printed in full).

**4. Roxygen/docs**:
- Add `@param h_grad` to `AdPaikModel()`'s roxygen block and extend the `@return`
  `@details` list (`AdPaikModel.R:65-108`) with a `GradientCheck` bullet describing its
  three fields.
- Give `gradient_check()` a roxygen block (title/description/params/return,
  `@keywords internal`), following the exact pattern of `params_se()`'s doc block
  (`R/params_se_CI.R:1-37`) — internal, not exported, but still gets a generated `.Rd`.
- Regenerate `NAMESPACE`/`man/*.Rd` via `roxygen2::roxygenise()` (RoxygenNote 7.3.3,
  confirmed installed locally) rather than hand-editing — no new `export()` entry is
  needed since `gradient_check()` stays internal like `params_se()` (not in current
  `NAMESPACE` export list either).

**5. `REVISION_NOTES.md`**: add a short new section documenting this change (mirroring
the existing "Author decisions" write-up style) — what was added, why, and the exact
default (`h_grad = 1e-4`) — so it's ready to lift into the JSS response-to-reviewers
letter.

## Verification

- `Rscript -e 'roxygen2::roxygenise()'` from the package root, confirm `NAMESPACE` is
  unchanged (no new exports) and new `man/gradient_check.Rd` is generated cleanly.
- `R CMD build` + `R CMD check` (or `devtools::check()`) to confirm no new warnings/notes.
- Run the package's own worked example (`data(data_dropout)`, same formula/time_axis as
  in `AdPaikModel.R`'s `@examples`, reduced `n_extrarun` for speed) and confirm:
  - `result$GradientCheck$Gradient` has length `n_params`, all finite (or the documented
    `NA` only at genuinely degenerate boundary parameters, e.g. `gamma_k` hitting the
    `se=1e-4` sentinel case already known from `REVISION_NOTES.md` §6).
  - `GradientNorm`/`GradientMaxAbs` are small relative to typical log-likelihood scale,
    consistent with a converged optimum.
  - `print(result)` / `summary(result)` show the new diagnostic line without breaking
    existing output.
  - `Loglikelihood`, `AIC`, `coef()` are **unchanged** from before this change (the
    diagnostic is read-only/additive, must not perturb the optimization itself).

## What actually changed vs. this plan

Verification surfaced two findings not anticipated in the original plan, both
addressed before considering this done. Full detail in `REVISION_NOTES.md` §7.

**1. Step-size instability (addressed above).** The first cut of `gradient_check()`
reused `params_se()`'s asymmetric clamp-then-blend formula. On the real worked example
this produced `Inf`/huge finite gradients for several parameters — traced to two
distinct causes and fixed:

**2. A genuine `gamma()`-overflow bug in the log-likelihood itself, found by the
diagnostic.** Three `gamma_k` converged very close to their lower range bound; the
existing `ll_AdPaik_centre_1D`/`ll_AdPaik_centre_eval` computed the log-likelihood's
third term via raw `gamma()` ratios, which overflow to `Inf` once `mu2/gamma_k`
exceeds ~171 — exactly the region a small perturbation there pushes into. Per
author's direction, fixed the root cause (not just masked in the diagnostic):
rewrote both functions' third-line computation in log-space (`lgamma`/`lchoose` +
log-sum-exp), verified mathematically identical to the old formula away from the
overflow edge (~1e-14 agreement over 2000 random draws) and at the previously-broken
point, and confirmed it stays finite under the exact perturbation that used to blow
up. This also obsoletes the ad hoc `res6==0`/`res7==0` → `1e-10` substitution from Q5
in `REVISION_NOTES.md` §3 (removed).

**3. A pre-existing, unrelated bug in `AdPaikModel()`'s own `@examples`.** Chasing an
unexpected `GenderMale ≈ 0` result (should be `0.2178` per the paper) led to finding
that the roxygen example's `categories_range_max` used `0` for the beta category,
while `Examples/ReplicationCode.R` (the actual paper-reproducing script) uses `0.5` —
silently clamping `GenderMale` to the boundary for anyone following the documented
example. Confirmed against the **unmodified, pre-session** code that this predates
all of today's other changes. Fixed the example's bound to `0.5`.

**Final verification, all three fixes together** (`data_dropout`, same
formula/`time_axis`, `set.seed(1)`, `n_extrarun=10`): `Loglikelihood = -2175.049`,
`AIC = 4398.098`, `coef = (GenderMale=0.2181, CFUP=-1.2706)`, `NRun = 31` — matches the
paper's published `-2175.135`/`4398.27`/`(0.2178, -1.2707)`/`31` to 2-3 decimal
places. `GradientCheck`: interior `‖grad‖₂ = 0.0221`, `max|grad| = 0.0162` (18/24
params); boundary-adjacent `mu1`, `nu`, `gamma_1`, `gamma_4`, `gamma_9`, `gamma_10`
(6/24) reported separately with their (large, expected) one-sided gradients.
