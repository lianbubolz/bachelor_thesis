## ============================================================
## utils.R -- small shared helpers
## ============================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

## Collapse condition messages into one short, CSV-safe string.
trunc_msg <- function(msgs, width = 200L) {
  if (length(msgs) == 0L) return("")
  s <- paste(unique(trimws(msgs)), collapse = " | ")
  s <- gsub("[\r\n]+", " ", s)
  if (nchar(s) > width) s <- paste0(substr(s, 1L, width - 3L), "...")
  s
}

capture_conditions <- function(expr) {
  warns <- character(0)
  msgs  <- character(0)
  err   <- ""
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )
  list(value = value,
       status = if (nzchar(err)) "error" else "ok",
       error = err, warns = warns, msgs = msgs)
}

strict_scalar <- function(z, what) {
  v <- suppressWarnings(as.numeric(unlist(z, use.names = FALSE)))
  if (length(v) != 1L || !is.finite(v)) {
    stop(what, " is not a single finite number (got: ",
         paste(utils::head(v, 3), collapse = ", "), ")")
  }
  v
}
