# =============================================================================
# MedMethods Explorer -- Shiny app entry point
# -----------------------------------------------------------------------------
# A local app to demonstrate and run the mediation-analysis methods in the
# MedMethods package. Architecture:
#   R/           generic, method-agnostic machinery (registry + module + I/O)
#   methods/<id> one self-contained plugin per method (registers itself)
#
# To add a method: create methods/<id>/<id>_method.R that calls
# register_method(...). It is auto-discovered and gets its own page.
#
# Run it with:  shiny::runApp("app")   or   Rscript app/run_local.R
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(bsicons)
  library(DT)
  library(plotly)
  library(shinycssloaders)
  library(markdown)
})

# Under shiny::runApp() the working directory is the app directory.
APP_DIR <- getwd()

# ---- load generic machinery (registry FIRST: defines %||%, register_method) --
# local = TRUE keeps definitions in this environment, where APP_DIR lives
# (runApp() evaluates app.R in a sandbox, not the global environment).
source(file.path(APP_DIR, "R", "registry.R"), local = TRUE)
for (f in c("io_helpers.R", "ui_helpers.R", "mod_method.R", "pkg_methods.R")) {
  source(file.path(APP_DIR, "R", f), local = TRUE)
}

# ---- discover & load method plugins -----------------------------------------
method_files <- list.files(file.path(APP_DIR, "methods"),
                           pattern = "_method\\.R$",
                           recursive = TRUE, full.names = TRUE)
for (mf in method_files) {
  tryCatch(source(mf, local = TRUE),
           error = function(e) warning(sprintf("Failed to load %s: %s",
                                               mf, conditionMessage(e))))
}
METHODS <- list_methods()

# ---- the full method catalogue, in the order used by the navbar and home ----
MED_FAMILY <- list(
  list(id = "macc",        name = "macc",
       desc = "Multilevel mediation with structured unmeasured mediator-outcome confounding; the error correlation is estimated, not assumed."),
  list(id = "gma",         name = "gma",
       desc = "Granger mediation analysis for time series with VAR(p) autoregressive errors."),
  list(id = "spcma",       name = "spcma",
       desc = "Sparse principal component mediation analysis for high-dimensional mediators."),
  list(id = "pathlasso",   name = "Pathway Lasso",
       desc = "Pathway estimation and selection with high-dimensional mediators; the penalty acts on the a-b products."),
  list(id = "pathlasso2b", name = "Multimodal",
       desc = "Two blocks of high-dimensional mediators, X -> M1 -> M2 -> Y."),
  list(id = "hdmediation", name = "HDMediation",
       desc = "Mediation when BOTH the exposures and the mediators are high-dimensional; selects mediators, exposures and direct effects."),
  list(id = "pcma",        name = "PCMA",
       desc = "Principal component mediation analysis for multiple exposures and multiple mediators."),
  list(id = "cfma",        name = "cfma",
       desc = "Causal functional mediation analysis: functional treatment, mediator and outcome."),
  list(id = "gmed",        name = "GMed",
       desc = "Mediation with a graph (covariance-matrix) mediator."),
  list(id = "hetermed",    name = "HeterMed",
       desc = "Heterogeneous mediation: moderated a- and b-paths give subject-specific indirect effects.")
)

# =============================================================================
# Home page
# =============================================================================
home_ui <- function() {
  cards <- lapply(MED_FAMILY, function(m) {
    spec <- get_method(m$id)
    impl <- !is.null(spec)
    bslib::card(
      class = if (impl) "h-100 border-primary" else "h-100",
      bslib::card_header(tags$b(m$name), " ",
                         if (impl) status_badge(spec$status %||% "ready")
                         else status_badge("planned")),
      bslib::card_body(
        tags$p(class = "small", m$desc),
        if (impl)
          tags$span(class = "small text-success",
                    bsicons::bs_icon("arrow-right-circle"),
                    " Open from the top navigation.")
        else
          tags$span(class = "small text-muted", "Coming soon.")
      )
    )
  })
  ver <- med_version()
  bslib::page_fluid(
    bslib::card(
      bslib::card_body(
        h2("MedMethods Explorer"),
        tags$p(class = "lead",
               "Run causal mediation analyses that go beyond a single scalar mediator."),
        tags$p("Each method has its own page. Read the model on the Overview tab, then on the Run / Demo tab either fit the built-in simulated example (whose ground truth is known, so the results include a truth-vs-estimate table) or upload your own data. Every result table and plot can be downloaded."),
        div(
          class = "alert alert-light border small",
          bsicons::bs_icon("hdd"), tags$b(" Runs locally. "),
          "Your data stays on this machine: uploads are read into the R session, used for the fit, and never written anywhere or sent off the machine."
        ),
        tags$hr(),
        tags$p(class = "text-muted small",
               sprintf("%d of %d methods wired.%s",
                       length(METHODS), length(MED_FAMILY),
                       if (!is.na(ver)) sprintf("  Engine: MedMethods %s.", ver) else ""))
      )
    ),
    h4("Methods"),
    do.call(bslib::layout_column_wrap, c(list(width = 1/2), cards))
  )
}

# =============================================================================
# Assemble navbar: Home + one page per implemented method
# =============================================================================
.fam_order   <- vapply(MED_FAMILY, function(m) m$id, character(1))
.method_ids  <- names(METHODS)
.ordered_ids <- c(intersect(.fam_order, .method_ids),
                  setdiff(.method_ids, .fam_order))
method_navs <- unname(lapply(.ordered_ids, function(id) {
  spec <- METHODS[[id]]
  bslib::nav_panel(title = spec$name, methodUI(spec$id, spec))
}))

ui <- do.call(
  bslib::page_navbar,
  c(
    list(
      title = tagList(bsicons::bs_icon("diagram-3"), " MedMethods"),
      id = "main_nav",
      theme = bslib::bs_theme(version = 5, bootswatch = "flatly",
                              primary = "#2c7fb8"),
      bslib::nav_panel(title = "Home", icon = bsicons::bs_icon("house"),
                       home_ui())
    ),
    method_navs,
    list(
      bslib::nav_spacer(),
      bslib::nav_item(tags$a(
        href = "https://github.com/zhaoyi1026/MedMethods", target = "_blank",
        rel = "noopener noreferrer", class = "nav-link",
        bsicons::bs_icon("github"), " Source"))
    )
  )
)

# =============================================================================
# Server: wire up each method's module
# =============================================================================
server <- function(input, output, session) {
  for (spec in METHODS) {
    methodServer(spec$id, spec)
  }
}

shinyApp(ui, server)
