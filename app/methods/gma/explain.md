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

`p` sets the autoregressive order. Lags of 2 and above work with asymptotic
variances in this package: the companion-matrix construction used by the
asymptotic variance had mismatched blocks in the original code, which made every
$p \ge 2$ fail, and that has been corrected here.
