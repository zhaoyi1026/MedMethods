# MedMethods Explorer (Shiny app)

A local point-and-click interface to every method in the
[**MedMethods**](https://github.com/zhaoyi1026/MedMethods) R package. Each method
gets its own page: read the model, fit the built-in simulated example, or upload
your own data, then download every result table and plot.

**It runs on your machine.** Uploaded data is read into the R session, used for
the fit, and never written to disk or sent anywhere.

## Install and run

From a clone of the repository:

```r
# 1. the engine (the package is at the repository root)
install.packages(".", repos = NULL, type = "source")
#    or: remotes::install_github("zhaoyi1026/MedMethods")

# 2. the UI packages
install.packages(c("shiny", "bslib", "bsicons", "DT", "plotly",
                   "shinycssloaders", "markdown"))

# 3. go
shiny::runApp("app", launch.browser = TRUE)
```

Or let the launcher do all three:

```bash
Rscript app/run_local.R          # installs what is missing, serves on :7800
Rscript app/run_local.R 8080     # a different port
Rscript app/run_local.R --check  # check prerequisites and exit
```

No compilation is needed — `MedMethods` is pure R.

## Method pages

| Page | Method | Data it wants |
|------|--------|---------------|
| macc | Multilevel mediation under structured confounding | one table of `Z, M, R` (+ an `id` column for the multilevel fit) |
| gma | Granger mediation for time series | one table of `Z, M, R`, one row per time point |
| spcma | Sparse principal component mediation | exposure vector, `n × p` mediators, outcome vector |
| Pathway Lasso | Pathway estimation and selection | exposure vector, `n × k` mediators, outcome vector |
| Multimodal | Two blocks of mediators | exposure, `n × p1` and `n × p2` mediator blocks, outcome |
| HD exposures | High-dimensional exposures **and** mediators | `n × q` exposures, `n × p` mediators, outcome |
| PCMA | Principal component mediation | `n × p` exposures, `n × q` mediators, outcome |
| cfma | Causal functional mediation | three `N × T` curve matrices (`Z`, `M`, `Y`) |
| GMed | Graph (covariance) mediator | a list of `T_i × p` matrices, `n × nX` exposures, outcome |
| HeterMed | Heterogeneous mediation effects | treatment coded ±1, mediator, outcome, `n × p` moderators |

## Using your own data

Each page's sidebar has an **Expected file format** panel listing exactly what
each input needs. Accepted formats:

- **CSV / TSV** — a header row; an obvious leading id column is dropped
  automatically for matrix inputs
- **`.rds`** — any R object of the right shape (this is the practical route for
  GMed's list of per-subject matrices)
- **`.RData` / `.rda`** — the stored object, or a named container from which the
  app picks the object matching the input's name

The fastest way to see the required layout is to switch the data source to
*Built-in example* and use the **Download** buttons — they export the example in
exactly the shape the upload expects, so you can open one and match its structure.

## What each page shows

- **Overview** — the model in maths, what the parameters do, practical cautions,
  and a link to the paper.
- **Run / Demo** — data source (example or upload), a parameter form generated
  from the method's declared parameters, and a data summary.
- **Results** — plots and tables as tabs. Every table has a *Download CSV*
  button, and plots that carry their underlying data offer one too. When the data
  is the built-in example there is also a **Truth vs estimate** table, because
  the generators return the true parameters.

## Notes worth reading before interpreting output

- **`macc` / `gma`:** with a *single* level or a *single* series the error
  correlation `delta` is not identifiable and must be supplied. Leaving it at 0
  assumes no confounding and biases the mediator-to-outcome path toward zero — on
  the `gma` example the true `B = -1` is estimated as `-0.003` at `delta = 0`.
- **`Pathway Lasso` / `HD exposures` / `Multimodal`:** these are tuning-sensitive
  and the jump from a dense solution to an empty one is abrupt. Read them as
  *pathway selection*; the penalties shrink effect magnitudes hard, so re-estimate
  the sizes of the selected pathways unpenalised.
- **`GMed`:** the b-path and indirect effect are attenuated at finite within-
  subject sample size, improving as `T_i` grows. Use the bootstrap for inference.
- **`PCMA` / `spcma` / `GMed`:** projections are identified only up to sign, so
  the truth comparisons use absolute cosine similarity.

## Adding a method

The app is a plugin host. `R/` holds the generic machinery and each method is one
self-contained file that registers itself:

```
app/
  app.R                     entry point; discovers plugins, builds the navbar
  run_local.R               launcher (installs prerequisites, then serves)
  R/registry.R              register_method() / list_methods()
  R/io_helpers.R            upload coercers for the five data shapes used here
  R/ui_helpers.R            parameter spec -> Shiny inputs, table formatting
  R/mod_method.R            ONE generic module that renders any method page
  R/pkg_methods.R           bridge to the installed MedMethods package
  methods/<id>/<id>_method.R   the plugin: example, parse, run, summarize, plots
  methods/<id>/explain.md       the Overview text
  tools/check_plugins.R     headless check of every plugin
```

Drop a new file under `methods/<id>/` that calls `register_method(list(...))` and
it is auto-discovered — no edits to `app.R` or the UI. The spec supplies `id`,
`name`, `full_name`, `data_inputs`, `params`, `example()`, `parse()`, `run()`,
`summarize()`, `plots()` and `explain`.

Unlike the companion CAP app, no private-environment cloning is needed: the
`MedMethods` package exports every method wrapper, so a plugin just calls
`med_fn("gmed")(...)`.

Write inline maths in `explain.md` as `$x$` and display maths as `$$...$$`.
`markdownToHTML()` converts `$x$` into the `\(x\)` delimiters MathJax typesets,
but *strips* a literal `\(`, so the dollar form is the one that works.

## Checking it without a browser

```bash
Rscript app/tools/check_plugins.R app
```

For every plugin this runs `example()` → `run()` (with the declared parameter
defaults) → `summarize()` → `plots()`, and round-trips the example through
`export_example()` → `parse()` so the upload path is exercised with data of
exactly the advertised shape. It prints each method's truth-vs-estimate table.
All ten pass.
