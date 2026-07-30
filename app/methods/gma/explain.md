## The model

A single time series (or several), with treatment $Z_t$, mediator $M_t$ and
outcome $R_t$:

$$M_t = Z_t A + E_{1t}, \qquad R_t = Z_t C + M_t B + E_{2t}$$

where the error process $(E_{1t}, E_{2t})$ follows a **vector autoregression of
order $p$**, so today's errors depend on the previous $p$ time points. That
is what makes this a *Granger* mediation model: it separates the mediation
structure from the temporal dependence rather than confusing the two.

As in `macc`, a correlation $\delta$ between the two error processes absorbs
unmeasured mediator--outcome confounding.

## Identification

For a **single series** $\delta$ is not identifiable and must be supplied. The
built-in example makes the consequence concrete: the true $B = -1$, and

- at $\delta = 0$ (assume no confounding) the fit gives $\hat B \approx -0.003$;
- at the true $\delta = 0.5$ it gives $\hat B \approx -1.00$.

With **several series** the correlation is identifiable and the two-level model
estimates it.

## VAR lag order

`p` sets the order of the autoregressive error process: how many previous time
points the errors at time $t$ depend on. Choose it the way you would for any
VAR — from the autocorrelation of the residuals, or by comparing fits — and note
that a higher order costs degrees of freedom, since the first $p$ observations
are used to condition on.

Both variance options work at any lag: the asymptotic variance is built from the
VAR companion matrix, and the empirical alternative is available from the same
sidebar.
