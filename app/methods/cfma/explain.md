## The model

Treatment, mediator and outcome are **all functions of time** -- one curve per
subject. The coefficients are therefore curves too.

The **concurrent** model links them instant by instant:

$$M(t) = Z(t)\,\alpha(t) + \epsilon_1(t)$$
$$Y(t) = Z(t)\,\gamma(t) + M(t)\,\beta(t) + \epsilon_2(t)$$

giving a time-varying indirect effect $\mathrm{IE}(t) = \alpha(t)\beta(t)$ and
direct effect $\mathrm{DE}(t) = \gamma(t)$.

The **historical influence** model instead lets the effect at time $t$
integrate over the recent past, which is appropriate when a treatment acts with a
lag rather than instantaneously.

## Estimation

The coefficient curves are expanded in a **Fourier basis** and fitted by
penalised least squares, with a roughness penalty on the second derivative.
`nbasis` sets how wiggly the curves may be and `lambda` how strongly curvature is
penalised. Cross-validation and bootstrap companions exist in the package for
choosing `lambda` and getting intervals.

## An identifiability note worth knowing

$\gamma(t)$, the direct effect, is the hardest curve to pin down: if $M(t)$
is nearly a deterministic multiple of $Z(t)$ then the two regressors are almost
collinear and $\gamma(t)$ is only weakly identified. The built-in example
deliberately gives the mediator substantial noise so that all three curves are
estimable; it recovers $\alpha(t)$, $\beta(t)$ and $\mathrm{IE}(t)$ with
correlations above 0.99 against the truth.
