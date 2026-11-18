continuous_check_static_support <- function(support, horizons = support$horizons) {
    # This is a pooled diagnostic for unchanged assignments, not a check of
    # positivity within every history or of support for a particular regime.
    for (horizon in horizons) {
        switched <- support$switched[support$time < horizon]
        n <- length(switched)
        if (!n) next
        count <- sum(switched)
        if (count == n)
            stop("All ", n, " observed treatment opportunities before horizon ", horizon,
                 " switch treatment; no unchanged assignments were observed at these opportunities.")
        if (count / n >= 0.95)
            warning(count, " of ", n, " treatment opportunities before horizon ", horizon,
                    " switch treatment; support for static interventions may be limited.", call. = FALSE)
    }
    invisible(NULL)
}
