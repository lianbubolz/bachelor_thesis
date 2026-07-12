## ============================================================
## scripts/_common.R
## ============================================================

## Locate the project root (the directory containing R/config.R),
## whether the script is run from the root, from scripts/, or via
## Rscript with a relative path.
locate_root <- function() {
  cands <- c(".", "..")
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if (length(fa)) {
    cands <- c(cands,
               file.path(dirname(sub("^--file=", "", fa[1])), ".."))
  }
  for (r in cands) {
    if (file.exists(file.path(r, "R", "config.R"))) {
      return(normalizePath(r))
    }
  }
  stop("cannot locate the project root (looked for R/config.R in: ",
       paste(cands, collapse = ", "), "); run from the project root")
}

## Source every file in R/ (alphabetical; all cross-references are
## inside function bodies, so order does not matter).
source_project <- function(root) {
  src <- file.path(root, "R")
  for (f in sort(list.files(src, pattern = "\\.R$", full.names = TRUE))) {
    source(f)
  }
  invisible(src)
}

## Everything after the raw results are obtained: validation, processed
## summaries, LaTeX tables, weight figure, population checks,
## console diagnostics. 
produce_outputs <- function(cfg, res,
                            run_tables = TRUE,
                            run_figure = TRUE,
                            run_checks = TRUE) {
  ensure_output_dirs(cfg)
  validate_results(res, cfg)
  if (run_tables) {
    write_summaries(res, cfg)
    summ <- write_all_tables(res, cfg)
    print_console_diagnostics(res, summ)
  }
  if (run_figure) {
    tryCatch(make_weight_figure(cfg),
             error = function(e) message("figure skipped: ",
                                         conditionMessage(e)))
  }
  invisible(NULL)
}
