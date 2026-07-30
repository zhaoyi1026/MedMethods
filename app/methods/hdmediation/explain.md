## The model

Both sides are high-dimensional: $q$ exposures and $p$ mediators.

$$M = X \alpha + E, \qquad Y = X \gamma + M \beta + e$$

Every (exposure, mediator) pair is a candidate pathway, so there are $q \times p$
of them -- 400 in the built-in example. The pathway effects form the matrix

$$\mu = \alpha \, \mathrm{diag}(\beta)$$

and the penalty acts on $\mu$, combining an elementwise $\ell_1$ term with a
**sparse-group** term across each exposure's row (so an exposure can be switched
off entirely), plus ridge terms on $\alpha$ and $\beta$. Solved by ADMM on
the splitting $\mu = \alpha\,\mathrm{diag}(\beta)$.

`pi` sets the balance between the elementwise and group terms, `phi` the ridge
weight, and `delta` an extra $\ell_1$ on $\alpha$ and $\beta$.

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
