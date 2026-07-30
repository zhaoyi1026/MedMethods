## The model

Trial-level data for each subject, with treatment $Z_t$, mediator $M_t$ and outcome $R_t$:

$$M_t = Z_t A + E_{1t}, \qquad R_t = Z_t C + M_t B + E_{2t}$$

The two error terms are **correlated**, with correlation $\delta$:

$$\mathrm{Cov}(E_{1t}, E_{2t}) = \delta\,\sigma_1\sigma_2$$

That correlation is exactly what unmeasured mediator--outcome confounding looks
like: something you did not measure moves both $M$ and $R$.

## Why the multilevel structure matters

In a **single** level, $\delta$ is *not identifiable* -- the data cannot
distinguish a genuine mediator--outcome path from correlated errors. You have to
supply $\delta$, and the conventional choice of 0 (assume no confounding)
biases the estimate of $B$ toward zero.

With **several subjects**, the variability across subjects identifies $\delta$,
so it can be estimated rather than assumed. That is the contribution of the
method, and it is what the built-in example demonstrates: the fit recovers
$\hat\delta \approx 0.58$ against a true 0.5.

## Effects reported

- `A` treatment to mediator, `B` mediator to outcome, `C` direct effect
- `C2` the total effect
- `AB.prod` the indirect effect as the product $A \times B$
- `AB.diff` the indirect effect as the difference $C_2 - C$

Two- and three-level models fit the subject-level coefficients with a mixed
model; the estimator can be a hierarchical likelihood (`HL`), a two-stage
procedure (`TS`), or `HL` followed by `TS`.
