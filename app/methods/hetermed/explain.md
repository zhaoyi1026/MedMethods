## The model

Treatment $X$ coded $\pm 1$, mediator $M$, outcome $Y$, and covariates
$Z$ (first column of ones). Both mediation paths are **moderated** by $Z$:

$$M_i = Z_i'(\alpha_0 + X_i \alpha_1) + e_i$$
$$Y_i = Z_i'(\gamma_0 + X_i \gamma_1) + (\beta_0 + X_i \beta_1) M_i + u_i$$

The $\alpha_1$, $\gamma_1$ and $\beta_1$ terms are the interactions with
treatment. Because they are there, the mediation effects are **subject-specific**:

$$\mathrm{NIE}_i = 2(\beta_0 + \beta_1 X_i)(Z_i'\alpha_1)$$
$$\mathrm{NDE}_i = 2 Z_i'\gamma_1 + 2\beta_1\left(Z_i'\alpha_0 + X_i Z_i'\alpha_1\right)$$

so instead of one average indirect effect you get a distribution of them across
subjects -- which is the point: it tells you *for whom* the treatment works
through the mediator.

## Estimation

The $\pm 1$ coding lets both arms be fitted by stacking their regressions; the
main and moderated coefficients are then the average and half-difference of the
two arms' estimates. Standard errors come from the model-based information matrix,
and the per-subject effects get delta-method intervals.

A generalized-lasso variant is also available, which shrinks the moderation
toward a common effect; it uses a bootstrap for inference and is slower.

## A note on the standard errors

The original code applied `sqrt()` to already-computed standard errors when
building its coefficient tables, which inflated every SE, z-value, p-value and
confidence interval by roughly a factor of seven. That has been corrected, so the
SEs shown here scale correctly with sample size and match a Monte Carlo
reference. The per-subject NIE/NDE intervals were always computed correctly.
