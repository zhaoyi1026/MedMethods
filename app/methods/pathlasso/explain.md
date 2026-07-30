## The model

One exposure $Z$, many candidate mediators $M$ ($k$ columns), one outcome:

$$M = Z A + E, \qquad Y = Z C + M B + e$$

Mediator $j$ carries indirect effect $A_j B_j$.

## The penalty

An ordinary lasso would penalise $A$ and $B$ separately. That is awkward, because
a mediator matters only if **both** of its paths are non-zero: a large $A_j$ with
$B_j = 0$ is not mediation. Pathway Lasso penalises the *product* directly, so
whole pathways are selected or dropped as units.

The criterion minimised is

$$f(A,B,C) = \frac{1}{2}\,\ell(A,B,C)
  + \lambda\left\{\sum_{j=1}^{K}\Big(|A_j B_j| + \phi\,(A_j^2 + B_j^2)\Big) + |C|\right\}
  + \omega\left\{\sum_{j=1}^{K}\big(|A_j| + |B_j|\big)\right\}$$

with the loss

$$\ell(A,B,C) = \mathrm{tr}\!\Big\{\Omega_1 (M - ZA)^{\top}(M - ZA)\Big\}
  + w_2\,(R - ZC - MB)^{\top}(R - ZC - MB),$$

where $\Omega_1 = \Sigma_1^{-1}$ and $w_2 = \sigma_2^{-2}$. Standardizing the data
to unit scale lets these be replaced by an identity matrix and one, which is what
this implementation does.

Writing the two penalty blocks as $P_1$ and $P_2$, the criterion is
$\frac{1}{2}\ell + \lambda P_1(A,B,C) + \omega P_2(A,B)$. $P_1$ shrinks the
pathway effects $A_jB_j$ and the direct effect $C$; $P_2$ adds separate shrinkage
on the individual $A_j$ and $B_j$, in the spirit of the elastic net.

**The role of $\phi$.** The quadratic term is not incidental — it is what makes
the penalty usable. $|ab|$ alone is *not* convex, and

$$v(a,b) = |ab| + \phi\,(a^2 + b^2)$$

is convex **if and only if $\phi \ge 1/2$**, and strictly convex when
$\phi > 1/2$. That is why the parameter cannot be set below $1/2$ here.

The optimisation is done by ADMM, splitting the product constraint off into its
own variable.

## Tuning parameters

| Sidebar | Symbol | Effect |
|---------|--------|--------|
| Penalty | $\lambda$ | overall strength on the pathway effects and $C$ |
| Ridge weight | $\phi$ | convexity; $\ge 1/2$ required |
| Extra $\ell_1$ on the paths | $\omega$ | separate shrinkage on individual $A_j$, $B_j$ |

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
