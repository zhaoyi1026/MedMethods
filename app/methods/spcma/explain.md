## The model

A single exposure $X$, a high-dimensional mediator $M$ ($p$ columns), and a
scalar outcome $Y$. Mediation is assessed not on the raw mediators but on a few
**principal component scores** $\tilde M_j = M \phi_j$:

$$\tilde M_{ij} = \alpha_j X_i + \xi_{ij}, \qquad
  Y_i = \gamma X_i + \sum_j \beta_j \tilde M_{ij} + \eta_i$$

so component $j$ carries indirect effect $\alpha_j \beta_j$.

## What makes the components *sparse*

The loadings $\phi_j$ are estimated with a **fused lasso** penalty, which
penalises differences between neighbouring loadings. The result is *piecewise
constant*: groups of adjacent mediators share a loading, and many are exactly
zero. That is far easier to interpret than a dense principal component, and it
matches the block structure often seen when mediators are ordered (adjacent
voxels, neighbouring regions, spectral bins).

`gamma` controls how much the loadings are additionally shrunk toward zero;
`per.jump` sets where along the fused-lasso path the loadings are read off.

## Reading the results

The built-in example has three true components, each a block of ten consecutive
mediators, with indirect effects 4, -1.5 and 0. The results overlay the
estimated loadings on the truth and report the component indirect effects with
bootstrap intervals. The plain-PCA fit is shown alongside for comparison.
