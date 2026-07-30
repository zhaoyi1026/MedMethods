## The model

Two mediator blocks **in sequence**. With exposure $X$, upstream block
$M_1$, downstream block $M_2$ and outcome $Y$:

$$M_1 = X \beta + e_1$$
$$M_2 = X \zeta + M_1 \Lambda + e_2$$
$$Y = X \delta + M_1 \theta + M_2 \pi + e$$

Because $M_1$ feeds into $M_2$, the model separates **three families of
pathway**:

| Pathway | Effect |
|---------|--------|
| $X \to M_1 \to Y$ | $\beta_j \theta_j$ |
| $X \to M_2 \to Y$ | $\zeta_k \pi_k$ |
| $X \to M_1 \to M_2 \to Y$ | $\beta_j \Lambda_{jk} \pi_k$ |

The third is what the two-block structure buys you: an effect that is invisible
if you analyse either block on its own.

## Penalties

As in Pathway Lasso the penalty acts on the path **products**, so each family is
selected as units. `kappa` sets the overall strength (applied to all four terms),
`mu` adds separate $\ell_1$ shrinkage on the individual paths, and `nu`
weights the ridge terms. The $M_1 \to M_2$ block $\Lambda$ is fitted with
a lasso.

This method needs **iterations in the thousands**; at a few hundred it reports
non-convergence.

## Reading the results

The built-in example puts signal in the first four mediators of each block. The
results give one table per pathway family, plus a heatmap of the cross-block
effects.
