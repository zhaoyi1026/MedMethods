# Using `medmethods`

One section per implemented method: the model, the expected input shapes, a
runnable example, and what to read off the result. Every example uses the
built-in generator, which also returns the ground truth.

```python
import numpy as np
import medmethods as mm
```

---

## `hetermed` — heterogeneous (moderated) mediation effects

Both the a-path and the b-path are moderated by covariates, so every subject has
their own natural indirect and direct effect.

    M_i = Z_i'(alpha0 + X_i alpha1) + e_i
    Y_i = Z_i'(gamma0 + X_i gamma1) + (beta0 + X_i beta1) M_i + u_i

**Inputs** — `X`: `(n,)` treatment coded **+1 / −1**; `M`: `(n,)` mediator;
`Y`: `(n,)` outcome; `Z`: `(n, p)` covariates whose **first column is ones**.

```python
d   = mm.examples.hetermed_example()             # n = 600
fit = mm.hetermed(d["X"], d["M"], d["Y"], d["Z"])
inf = mm.hetermed_inf(d["X"], d["M"], d["Y"], d["Z"], fit)

fit["beta0"], fit["beta1"]         # 0.808, 0.201     (true 0.8, 0.2)
fit["alpha1"]                      # moderated a-path, one entry per column of Z
fit["ITE"]["NIE"]                  # per-subject indirect effect
inf["alpha1"]["pvalue"]            # asymptotic p-values
inf["NIE"]["SE"]                   # delta-method SE per subject
```

Effects for arbitrary parameter values: `mm.hetermed_ite(X, Z, alpha0, alpha1,
beta0, beta1, gamma0, gamma1)`.

Only the OLS estimator is ported; `method="genlasso"` raises
`NotImplementedError`.

---

## `pcma` — principal component mediation analysis

Orthogonal projections of the exposures and of the mediators define scalar
component scores carrying parallel mediation mechanisms.

    M psi = (X phi) alpha + e
    Y     = (X phi) gamma + (M psi) beta + u,        IE = alpha * beta

**Inputs** — `X`: `(n, p)` exposures; `M`: `(n, q)` mediators; `Y`: `(n,)`.

```python
d   = mm.examples.pcma_example()                 # n = 400, p = 5, q = 10
fit = mm.pcma(d["X"], d["M"], d["Y"], nD=2, ninitial=3)

fit["Phi"]        # (p, nD) exposure projections
fit["Psi"]        # (q, nD) mediator projections
fit["alpha"], fit["beta"], fit["gamma"], fit["IE"]    # one entry per component

# projections are identified only up to sign -- compare with |cosine|
cos = lambda a, b: abs(a @ b) / np.sqrt((a @ a) * (b @ b))
cos(fit["Phi"][:, 0], d["truth"]["Phi"][:, 0])        # ~0.94
```

Bootstrap inference, either inline or at fixed projections:

```python
fit = mm.pcma(d["X"], d["M"], d["Y"], nD=1, boot=True, sims=500)
fit["inference"][0]["SE"]

inf = mm.pcma_inf(d["X"], d["M"], d["Y"], fit["Phi"][:, 0], fit["Psi"][:, 0],
                  sims=500)
inf["names"], inf["estimate"], inf["pvalue"]          # alpha, beta, gamma, IE, DE
```

`mm.pcma_coef(X, M, Y, phi, psi)` gives the path coefficients at projections you
supply.

---

## `gmed` — mediation with a graph (covariance) mediator

The exposure shifts the log-variance of a subject-level covariance mediator along
a projection `theta`, and the outcome depends on that log-variance.

    log(theta' Sigma_i theta) = alpha0 + X_i' alpha + b_i
    Y_i = gamma0 + X_i' gamma + beta * log(theta' Sigma_i theta) + e_i

**Inputs** — `X`: `(n,)` or `(n, nX)` exposure and covariates, **no intercept
column**, first column the exposure of interest; `M`: length-`n` list of
`(T_i, p)` matrices; `Y`: `(n,)`.

```python
d   = mm.examples.gmed_example()                 # n = 100, p = 10, T_i = 500
fit = mm.gmed(d["X"], d["M"], d["Y"], nD=1, ninitial=5)

fit["theta"]                       # (p, nD) mediating directions
fit["alpha"], fit["beta"]          # a-path, b-path per direction
fit["IE"], fit["DE"]               # indirect and direct effects

bt = mm.gmed_boot(d["X"], d["M"], d["Y"], theta=fit["theta"][:, 0], sims=500)
bt["estimate"][3], bt["pvalue"][3]                    # index 3 is IE
```

**Leave `H` at its default.** `H=None` uses the average covariance; an identity
`H` sends the search to a covariate-free background direction.

Use `stop_crt="DfD"` with `DfD_thred` to let the number of directions be chosen
by deviation-from-diagonality instead of fixing `nD`. `mm.gmed_refit(X, M, Y,
Theta)` re-estimates the effects at projections you supply.

The b-path and `IE` are attenuated at finite `T_i` — the projected log-variance is
a noisy measurement of the true one. It improves as `T_i` grows.

---

## `hdmediation` — high-dimensional exposures and mediators

    M = X alpha + E,        Y = X gamma + M beta + e

with a pathway-lasso penalty on the exposure→mediator→outcome products
`mu = alpha diag(beta)`, solved by ADMM.

**Inputs** — `X`: `(n, q)` exposures; `M`: `(n, p)` mediators; `Y`: `(n,)`.
Needs `n > q`; use `hdmediation_pca()` otherwise.

```python
d   = mm.examples.hdmediation_example()          # n = 100, q = 20, p = 20
fit = mm.hdmediation(d["X"], d["M"], d["Y"],
                     lam=2.5, pi=0.5, phi=2, delta=0.5)

fit["IE"]                          # (q, p) pathway effects
np.argwhere(np.abs(fit["IE"]) > 1e-3)            # selected exposure-mediator paths
fit["alpha"], fit["beta"], fit["gamma"]
fit["logLik"]                      # dict: logLik, R1, R2, obj
fit["converge"]                    # check this -- see below
```

Tuning: `lam` overall strength, `pi` the l1-vs-group split, `phi` ridge, `delta`
extra l1. The dense→empty transition is abrupt, so scan `lam` and select by BIC
rather than guessing. `fit["out_scaled"]` holds the fit on the standardized scale.

Principal-component variant when the exposures are high-dimensional:

```python
fit = mm.hdmediation_pca(d["X"], d["M"], d["Y"], adaptive=True, var_prop=0.9,
                         lam=1.0, pi=0.5, phi=2, delta=0.5)
fit["n_pc"], fit["rotation"], fit["var_explained"]
```

---

## `pathlasso` — Pathway Lasso

    M = Z A + E,        Y = Z C + M B + e

The penalty acts on the *products* `A_j B_j` rather than on `A` and `B`
separately, so pathways are selected as units.

**Inputs** — `X`: `(n,)` exposure; `M`: `(n, k)` mediators; `Y`: `(n,)`.

```python
d   = mm.examples.pathlasso_example()            # n = 200, k = 50, signal 0..3
fit = mm.pathlasso(d["X"], d["M"], d["Y"], lam=0.001, omega=0, phi=1,
                   max_itr=3000, tol=1e-8)

fit["A"], fit["B"], fit["C"]       # a-path (1, k), b-path (k, 1), direct effect
fit["AB"]                          # (k,) pathway effects -- what to threshold
np.flatnonzero(np.abs(fit["AB"]) > 0.05 * np.max(np.abs(fit["AB"])))
fit["constraint"]                  # ADMM splitting gaps; should be ~0
```

`lam` lives on a small scale for standardized data (near 1e-3 here); at 1e-2 the
b-path collapses to exactly zero. Choose it by selection stability:

```python
out = mm.pathlasso_ksc(d["X"], d["M"], d["Y"],
                       lam_grid=[1e-4, 1e-3, 1e-2], n_rep=5)
out["lambda_est"], out["mean_kappa"]
```

Check `mean_kappa` before trusting `lambda_est`: the criterion splits the sample
in half, so on a small or high-dimensional problem the two halves may select
almost disjoint sets and every kappa sits near zero (this happens on the bundled
example, where each half has only 100 subjects for 50 mediators). Kappa near zero
means "no more agreement than chance", and the returned `lambda_est` is then
arbitrary — use a larger sample or fewer candidate mediators.

`mm.prox_pair(...)` exposes the closed-form solution of the penalised pair
subproblem if you want to build your own solver around it.

---

## Not yet ported

`pathlasso2b`, `cfma_concurrent`, `cfma_historical`, `spcma`, `mcma_pca`,
`mcma_bk`, `macc` and `gma` exist but raise `NotImplementedError` naming the R
function to use instead. See the roadmap in [README.md](README.md).

```python
mm.macc(...)
# NotImplementedError: macc (multilevel mediation under structured confounding)
# is not yet ported to Python. It maximises a hierarchical likelihood over the
# variance components, which R delegates to nlme/lme4. Use the R implementation
# MedMethods::macc() meanwhile; ...
```

## References

```python
mm.references()            # every method
mm.references("pathlasso")  # one
mm.REFERENCES["gmed"][0]["doi"]
```
