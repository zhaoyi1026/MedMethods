## The model

Both sides are high-dimensional: $q$ exposures and $p$ mediators.

$$M = X \alpha + E, \qquad Y = X \gamma + M \beta + e$$

Every (exposure, mediator) pair is a candidate pathway, so there are $q \times p$
of them -- 400 in the built-in example.

## The penalty

The criterion minimised is

$$\frac{1}{2}\,L(\alpha,\beta,\gamma)
  + \lambda_1 R_1(\alpha,\beta)
  + \lambda_2 R_2(\mu)
  + \lambda_3 R_3(\gamma)$$

with the loss

$$L(\alpha,\beta,\gamma)
  = \mathrm{tr}\!\Big\{(M - \tilde{X}\alpha)^{\top}(M - \tilde{X}\alpha)\Big\}
  + (Y - \tilde{X}\gamma - M\beta)^{\top}(Y - \tilde{X}\gamma - M\beta)$$

and three penalty terms, each doing a different job.

**$R_1$ selects individual mediators.** It generalises the Pathway
Lasso penalty to $q$ exposures:

$$R_1(\alpha,\beta) =
  \sum_{j=1}^{q}\sum_{k=1}^{p}\Big\{|\alpha_{jk}\beta_k| + c_0\big(\alpha_{jk}^2 + \beta_k^2\big)\Big\}
  + c_1\left(\sum_{j=1}^{q}\sum_{k=1}^{p}|\alpha_{jk}| + \sum_{k=1}^{p}|\beta_k|\right)$$

For a given mediator $M_k$, the product term $\sum_j |\alpha_{jk}\beta_k|$ pushes
*all* paths through $M_k$ to zero together, which is what achieves mediator
selection. The quadratic term is what makes it convex: $|ab| + c_0(a^2+b^2)$ is
convex when $c_0 \ge 1/2$, and the implementation fixes $c_0 = 2$. The final term
is an ordinary lasso on the individual path coefficients.

**$R_2$ selects individual exposures.** A group lasso over each
exposure's row of pathway effects $\mu_{jk} = \alpha_{jk}\beta_k$:

$$R_2(\mu) = \sum_{j=1}^{q}\sqrt{p}
  \left(\sum_{k=1}^{p}\mu_{jk}^2\right)^{1/2}$$

so all paths originating from exposure $\tilde{X}_j$ are shrunk to zero together.

**$R_3$ selects direct effects.** A plain lasso,
$R_3(\gamma) = \sum_{j=1}^{q}|\gamma_j|$.

Introducing $\mu_{jk} = \alpha_{jk}\beta_k$ as a separate parameter turns the
problem into a **sparse group lasso**, which has a closed-form solution, and the
constraint $\mu_{jk} = \alpha_{jk}\beta_k$ is enforced by ADMM.

## How the sidebar maps onto that

| Sidebar | Symbol | Role |
|---------|--------|------|
| Penalty | $\lambda$ | overall strength |
| $\ell_1$ vs group split | $\pi$ | $\lambda_1 = \pi\lambda$, $\lambda_2 = (1-\pi)\lambda$; $\lambda_3 = \lambda$ |
| Ridge weight | $c_0$ | convexity of the product term; $\ge 1/2$ required, paper uses 2 |
| Extra $\ell_1$ | $c_1$ | shrinkage on individual $\alpha_{jk}$, $\beta_k$ |

So $\pi = 1$ puts all the weight on mediator selection and $\pi = 0$ all of it on
exposure selection.

## Two practical cautions

**Dimensions.** The direct estimator regresses on all $q$ exposures at once, so
it needs $n > q$. When the exposures approach the sample size, switch on the
principal-component option, which replaces them by their leading components.

**Read it as selection, not estimation.** The product penalty shrinks the effect
*magnitudes* hard. On the built-in example at the default lambda the three true
pathways are selected out of 400 candidates, but their estimated effects are far
below the true 2.25. Use it to find pathways, then re-estimate their sizes
unpenalised. Lower lambda for less shrinkage -- at the cost of far less sparsity
(2 keeps 50 paths, 2.5 keeps 22, 5 keeps none).
