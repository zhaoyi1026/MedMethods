# medmethods — mediation analysis methods (Python)

A Python port of the [**MedMethods**](https://github.com/zhaoyi1026/MedMethods) R
package: causal mediation analysis for settings that go beyond a single scalar
mediator. Pure `numpy` / `scipy`.

The R package remains the reference implementation. This port is being brought up
to parity method by method, each one cross-checked numerically against R.

## Status

| Function | Method | Status |
|----------|--------|--------|
| `hetermed()` | Heterogeneous (moderated) mediation effects | ✅ **bit-identical** to R |
| `hdmediation()` | High-dimensional exposures **and** mediators | ✅ **bit-identical** to R |
| `pathlasso()` | Pathway Lasso | ✅ **bit-identical** to R |
| `pcma()` | Principal component mediation analysis | ✅ verified vs R |
| `gmed()` | Mediation with a graph (covariance) mediator | ✅ verified vs R |
| `pathlasso2b()` | Two blocks of high-dimensional mediators | ⏳ planned |
| `cfma_concurrent()`, `cfma_historical()` | Causal functional mediation | ⏳ planned |
| `spcma()`, `mcma_pca()`, `mcma_bk()` | Sparse principal component mediation | ⏳ planned |
| `macc()` | Multilevel mediation, structured confounding | ⏳ planned |
| `gma()` | Granger mediation analysis | ⏳ planned |

The planned functions exist and raise `NotImplementedError` with a pointer to the
R equivalent, so `dir(medmethods)` reflects the full method catalogue.

**29 pytest tests pass** (Python 3.7, numpy 1.16, scipy 1.2).

## Install

```bash
git clone https://github.com/zhaoyi1026/MedMethods.git
pip install ./MedMethods/Python
```

Requires Python ≥ 3.7, `numpy` ≥ 1.16 and `scipy` ≥ 1.2 — nothing else. (The code
avoids `numpy.random.default_rng`, which arrived in numpy 1.17, so it runs on old
environments.)

## Quick start

```python
import medmethods as mm

d = mm.examples.hetermed_example()          # n = 600, treatment coded +/-1
fit = mm.hetermed(d["X"], d["M"], d["Y"], d["Z"])
inf = mm.hetermed_inf(d["X"], d["M"], d["Y"], d["Z"], fit)

fit["beta0"], fit["beta1"]        # 0.808, 0.201   (true 0.8, 0.2)
inf["alpha1"]["SE"]               # asymptotic standard errors
fit["ITE"]["NIE"]                 # per-subject natural indirect effect
```

Every implemented method has a matching generator in `medmethods.examples` that
returns the data plus a `truth` dict, so a full run is two lines. See
[USAGE.md](USAGE.md) for each method.

## Verification

Each port was checked against R on identical data. Where a method is called
"bit-identical", every reported quantity agreed to floating-point noise.

**`hetermed`** — bit-identical:

- coefficients (`alpha0`, `alpha1`, `beta0`, `beta1`, `gamma0`, `gamma1`):
  max |Δ| = 1.6e-15
- asymptotic standard errors: max |Δ| = 9.0e-17
- per-subject NIE / NDE and their delta-method SEs: max |Δ| = 1.1e-14

**`hdmediation`** — bit-identical at `lambda` = 0.5 and 2 (`pi` = 0.5, `phi` = 2,
`delta` = 0.5): `alpha`, `beta`, `gamma` and the pathway effects `IE` all within
1.2e-14, and the port reproduces R's non-convergence flag on the same input.

**`pathlasso`** — bit-identical at `lambda` = 1e-3 and 1e-2: `A`, `B`, `C` and
`AB` within 7.3e-14, with the same iteration counts (1112 and 819).

**`pcma`** — verified:

- `pcma_coef` (path coefficients at fixed projections): bit-identical, 3.9e-14.
- `pcma_d1_base` from an identical starting projection: same iteration count
  (76), projection cosines to R = 1.000000, coefficients within 1.2e-7. The
  residual gap is the root-finder tolerance (`scipy.optimize.brentq` vs R's
  `uniroot`) in the unit-norm ridge step.
- The multi-start search draws its starts from numpy's RNG, so individual starts
  differ from R's; pass `phi0_mat` / `psi0_mat` to force identical ones.

**`gmed`** — verified:

- `gmed_coef` (path coefficients at a fixed projection): bit-identical, 7.8e-14.
- `gmed_d1` from an identical start: projection cosine to R = 0.99997, path
  coefficients within 6e-4. The gap is R's `nlme::lme` M-model fit — see below.

### Two deliberate differences from R

**1. `gmed`'s random-effect split.** R fits the mediator model with
`nlme::lme(score ~ X, random = ~1|ID)` where each subject contributes a *single*
score, so the random-intercept groups are singletons: the marginal likelihood
depends only on the sum of the random-intercept and residual variances, making
their split non-identifiable and R's returned value optimiser-arbitrary. This
port uses a deterministic fit — fixed effects are GLS = OLS (identifiable), with
`blup_shrink` controlling the BLUP. The identifiable quantities (`theta`,
`alpha`, `beta`, `gamma`, `IE`) agree with R; `alpha0_rnd` and `tau2` may not.

**2. `hetermed`'s standard errors are corrected.** R's `fit.inf.OLS` computes
`vecTheta.se <- sqrt(diag(cov.vecTheta))` — already standard errors — and then
builds each coefficient table with `SE = sqrt(Theta.se[...])`, taking the square
root a second time. This port reports the correct standard errors. Two
independent confirmations:

- the R value scales as `n^(-1/4)` rather than `n^(-1/2)`: the ratio
  SE(n=600)/SE(n=2400) is 1.42, where a correct SE gives 2.00;
- the *square* of the R value matches the Monte Carlo SD of the estimator over
  300 replicates, while the R value itself is about 7× too large.

R's NIE/NDE tables use the covariance matrix directly and were always correct.
The R package has since been fixed to match.

## Known behaviour worth knowing

- **`gmed` needs `H = None`** (the default), which uses the average covariance.
  Passing an identity `H` collapses the direction search onto a covariate-free
  background direction. There is a test for this.
- **`gmed`'s b-path is attenuated at finite `Ti`.** The projected log-variance is
  an error-prone measurement of the true log-variance, so the outcome-model slope
  is biased toward zero; it improves as `Ti` grows. Use `gmed_boot()` for
  inference on the indirect effect.
- **The penalised methods are tuning-sensitive** and the transition from dense to
  empty solutions is abrupt. On standardized data `pathlasso` wants `lambda` near
  1e-3. Select over a path with `pathlasso_ksc()` (selection-stability / kappa
  criterion) rather than guessing.
- **`hdmediation` needs `n > q`** — it regresses on all exposures directly. Use
  `hdmediation_pca()` when the number of exposures approaches the sample size.
- `lambda` is a Python keyword, so it is spelled `lam` in `pathlasso()` and
  `hdmediation()`.

## Architecture

```
medmethods/
  _core.py        # shared: soft threshold, covariance cubes, projected scores,
                  #         weighted accumulation, generalized eigen solve, DfD
  hetermed.py     # heterogeneous mediation            [done]
  pcma.py         # principal component mediation       [done]
  gmed.py         # graph (covariance) mediator         [done]
  hdmediation.py  # HD exposures + mediators (ADMM)     [done]
  pathlasso.py    # Pathway Lasso (ADMM) + KSC tuning   [done]
  examples.py     # example-data generators + ground truth
  references.py   # per-method citations; medmethods.references()
  pathlasso2b.py cfma.py spcma.py macc.py gma.py        [planned]
```

## Porting roadmap

In rough order of tractability:

1. **`pathlasso2b`** — the ADMM mirrors `pathlasso` and reuses the same
   closed-form pair proximal step; the only missing piece is the per-column lasso
   for the `Lambda` block, which R delegates to `glmnet`. Needs a self-contained
   coordinate-descent lasso matching glmnet's `(1/2n)‖y − b₀ − Xb‖² + λ‖b‖₁`.
2. **`mcma_pca`** (part of `spcma`) — PCA of the mediator plus a per-component
   single-mediator mediation fit and bootstrap. No special machinery.
3. **`cfma`** — deterministic penalised least squares, but needs the Fourier
   basis, the second-derivative penalty, and the concurrent / historical design
   builders, plus the CV and bootstrap wrappers.
4. **`gma`** (single-level) — closed-form linear algebra including the VAR
   companion-matrix asymptotic variance. The two-level model needs
   `nlme::gls`-style generalized least squares and a profile likelihood.
5. **`macc`** — needs a linear-mixed-model fit (REML) to replace `nlme`/`lme4`.
6. **`spcma`** (sparse) — needs a generalized-lasso (fused lasso) path solver
   equivalent to R's `genlasso`; `mcma_bk` additionally leans on the CRAN
   `mediation` package.

## Running the tests

```bash
pip install "./Python[test]"
cd Python && pytest -q
```

## References

`medmethods.references()` prints all of them; `medmethods.references("gmed")`
prints one. The citations mirror the R package's References section.

## License

GPL-3.
