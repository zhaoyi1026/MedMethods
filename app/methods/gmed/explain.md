## The model

The mediator is not a number but a **covariance matrix**: for subject $i$ you
observe a $T_i \times p$ data matrix whose covariance $\Sigma_i$ is the
mediator. Think of a connectivity matrix between $p$ regions.

A projection $\theta$ reduces it to a scalar, the projected log-variance, and
the mediation model runs through that scalar:

$$\log(\theta' \Sigma_i \theta) = \alpha_0 + x_i'\alpha + b_i$$
$$Y_i = \gamma_0 + x_i'\gamma + \beta \log(\theta' \Sigma_i \theta) + e_i$$

with indirect effect $\mathrm{IE} = \alpha_x \beta$ and direct effect
$\gamma_x$. $\theta$ is **learned**, not fixed in advance: the method looks
for the direction in which the exposure moves the covariance and the covariance
moves the outcome.

## Two things that matter in practice

**The constraint matrix.** $\theta$ is normalised by $\theta' H \theta = 1$.
The app uses $H$ = the average covariance, which is the right choice; an
identity $H$ sends the search to whichever direction simply has the largest
variance, which is typically covariate-free background.

**The b-path is attenuated at finite $T_i$.** Each subject's covariance is
estimated with error of order $\sqrt{2/T_i}$, so the projected log-variance is
an error-prone regressor and the outcome-model slope $\beta$ is biased toward
zero. It improves as $T_i$ grows -- on the built-in setting $\hat\beta$ is
about 0.47 at $T_i = 150$ and 0.84 at $T_i = 500$. Use the bootstrap for
inference on the indirect effect rather than reading the point estimate as exact.

## Also available

The same estimator ships as `capmediation()` in the companion
[CAPmethods](https://github.com/zhaoyi1026/CAPmethods) package, which collects the
covariate-assisted principal regression family.
