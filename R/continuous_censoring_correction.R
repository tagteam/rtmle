# Per-subject censoring-martingale contribution to the one-step estimator.
continuous_censoring_correction <- function(nuisance, th, grid, q_grid, weights) {
    data <- nuisance$data
    info <- nuisance$prepared_data_object$info
    last_event <- info$last_event[match(th, info$time_horizon)]
    correction <- numeric(nrow(data))

    divide_by_survival <- function(numerator, survival) {
        active <- numerator != 0
        if (any(!is.finite(survival[active]) | survival[active] <= 0))
            stop("Nonpositive event survival in the continuous-time censoring correction.")
        numerator[active] <- numerator[active] / survival[active]
        numerator
    }
    for (k in seq_len(last_event)) {
        # The last recursion step pools subsequent events through the terminal
        # event. Its q and survival models must use the same endpoint.
        pooled <- k == last_event
        suffix <- paste0(if (pooled) "_pooled_" else "_", k)
        entry <- data[[paste0("time_", k - 1L)]]
        next_time <- data[[paste0("time", suffix)]]
        censored <- data[[paste0("event", suffix)]] == "C"
        exit <- pmin(next_time, th)
        eligible <- which(weights[, k] > 0 & entry < exit)
        if (!length(eligible)) next

        q <- q_grid[[k]]
        q_tau <- q[, length(grid)]
        previous_censoring_survival <- rep(1, nrow(data))
        for (j in seq_along(grid)[-1L]) {
            left <- grid[j - 1L]
            right <- grid[j]
            rows <- eligible[entry[eligible] < right & exit[eligible] > left]
            if (!length(rows)) next
            u <- pmin(right, exit[rows])
            from <- pmax(left, entry[rows])
            # q(u) includes events at u. Approach an off-grid time from the
            # right so the approximation also includes a possible atom at u.
            q_difference <- q_tau[rows] - q[rows, j]
            jump_rows <- which(censored[rows] & next_time[rows] == u)
            survival <- continuous_grid_survival(nuisance, k, pooled, rows, u,
                                                  from = from, jump_rows = jump_rows)
            # Each cell covers (from,u]. This conditional dropout probability
            # is the hazard increment for the grid's product-integral model.
            hazard <- 1 - survival$censoring / previous_censoring_survival[rows]
            if (any(!is.finite(hazard)) || any(hazard < -sqrt(.Machine$double.eps)))
                stop("Censoring survival must be nonincreasing over each event interval.")
            previous_censoring_survival[rows] <- survival$censoring
            # A censoring-hazard atom inside a cell must use event survival
            # before the atom, including for subjects who outlast the cell.
            # For the compensator use S(from); for the observed jump use S(u-).
            change <- -divide_by_survival(q_difference * pmax(0, hazard), survival$event)
            change[jump_rows] <- change[jump_rows] +
                divide_by_survival(q_difference[jump_rows], survival$event_at_jump)
            correction[rows] <- correction[rows] + weights[rows, k] * change / survival$censoring
        }
    }
    if (any(!is.finite(correction)))
        stop("Nonfinite continuous-time censoring correction.")
    correction
}
