# Dom Bennett
# 25/09/2015
# Tools for debugging

istore <- function (...) {
  # Store variables specific variables in inspectenv
  if (!exists ("inspectenv", envir=.GlobalEnv)) {
    inspectenv <- new.env (parent=.GlobalEnv)
  }
  args <- as.list (match.call ())
  for (i in 2:length (args)) {
    assign (names (args)[i], args[[i]], envir=inspectenv)
  }
}
