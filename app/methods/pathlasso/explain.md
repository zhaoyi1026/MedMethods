## The model

One exposure $Z$, many candidate mediators $M$ ($k$ columns), one outcome:

$$M = Z A + E, \qquad Y = Z C + M B + e$$

Mediator $j$ carries indirect effect $A_j B_j$.

## What is different about the penalty

An ordinary lasso would penalise $A$ and $B$ separately. That is awkward,
because a mediator matters only if **both** of its paths are non-zero: a large
$A_j$ with $B_j = 0$ is not mediation. Pathway Lasso penalises the
**product**

$$\sum_j |A_j B_j| \;+\; \text{ridge and } \ell_1 \text{ terms}$$

so whole pathways are selected or dropped as units. The optimisation is done by
ADMM, splitting the product constraint off into its own variable.

## Choosing lambda

`lambda` lives on a **small scale** for standardized data -- around $10^{-3}$
on the built-in example. The transition from a dense solution to an empty one is
abrupt: at $10^{-2}$ the b-path collapses to exactly zero. Adding `omega`
(separate $\ell_1$ shrinkage) can also zero out the whole solution even at
small values.

Rather than guessing, tick *Pick lambda by selection stability*: the sample is
split in half repeatedly and the agreement between the two halves' selected sets
is measured by Cohen's kappa. **Check the kappa values before trusting the
choice** -- if they sit near zero, the two halves are selecting nearly disjoint
sets and the criterion is uninformative at that sample size.

## Reading the results

The built-in example has four signal mediators with products $\pm 4$ among
fifty. At the default lambda the fit recovers exactly those four, with products
of 3.86, -3.94, -4.27 and 3.77, and shrinks the other 46 to zero.
