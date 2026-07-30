## The model

Several exposures and several mediators, reduced to scalar component scores by
two **orthogonal projections**: $\phi$ on the exposures, $\psi$ on the
mediators.

$$M \psi = (X \phi)\,\alpha + e, \qquad
  Y = (X \phi)\,\gamma + (M \psi)\,\beta + u$$

Component $j$ then has a-path $\alpha_j$, b-path $\beta_j$, direct effect
$\gamma_j$ and indirect effect $\alpha_j \beta_j$. The assumption is that
there exist orthogonal directions along which mediation operates in **parallel**,
one mechanism per component.

## Estimation

The projections and effect parameters are estimated jointly by minimising the
negative joint log-likelihood, alternating between the effects, the variances,
and each projection in turn. Each projection update is a ridge-type solve whose
penalty is chosen so the projection has unit norm.

The likelihood surface has local optima, so the search is restarted from several
random projections and the best one is kept -- `ninitial` controls how many.
Higher components are extracted by deflating the data and repeating.

Inference is by bootstrap over subjects at the fitted projections.

## Reading the results

Projections are identified only **up to sign**, so compare them with the truth
using the absolute cosine (as the truth-vs-estimate table does) rather than
element by element. On the built-in example the leading exposure and mediator
projections both recover the truth with cosine around 0.98.
