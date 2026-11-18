# Collect function calls without mistaking a column name for a special term.
continuous_formula_calls <- function(expression) {
    if (!is.call(expression)) return(character())
    head <- expression[[1L]]
    name <- if (is.symbol(head)) as.character(head) else if (
        is.call(head) && as.character(head[[1L]]) %in% c("::", ":::")) {
        as.character(head[[3L]])
    } else character()
    c(name, unlist(lapply(as.list(expression)[-1L], continuous_formula_calls), use.names = FALSE))
}
