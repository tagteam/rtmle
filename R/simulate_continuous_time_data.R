#' Simulate continuous-time treatment and event data
#'
#' Simulates binary treatment A and a time-varying covariate L, with
#' treatment opportunities, covariate updates, censoring, an outcome, and a
#' competing event. Waiting times follow exponential proportional-hazard
#' models, except regularly scheduled treatment opportunities. The output
#' can be separated into the outcome, censoring, competing-risk and covariate
#' tables accepted by \code{add_continuous_data}; baseline data are added with
#' \code{add_baseline_data}.
#'
#' @param n Positive integer sample size.
#' @param effects Named list of named coefficient lists. Event-specific
#'   coefficients have suffixes identifying the event index; use
#'   \code{aggregate_parameters = TRUE} to reuse the first event's coefficients.
#' @param static_intervention Treatment value 0 or 1 imposed at baseline and
#'   subsequent treatment opportunities. NULL uses observed treatment assignment.
#'   Specifying an intervention also disables censoring.
#' @param static_intervention_baseline Treatment value at baseline, or NULL to
#'   draw baseline treatment from its model. Overridden by static_intervention.
#' @param baseline_rate_list Named list of finite non-negative rates A, L, C,
#'   Y, and D. Zero disables the corresponding exponential clock.
#' @param max_fup Positive maximum follow-up time.
#' @param visitation_interval Non-negative mean waiting time for scheduled visits.
#' @param visitation_sd Non-negative standard deviation for scheduled visits.
#'   Negative sampled waiting times are set to zero. One time unit is added
#'   between successive events.
#' @param discretize_age Logical; draw age as 50 or 70 rather than uniformly
#'   between 40 and 90.
#' @param no_competing_events Logical; disable competing events.
#' @param uncensored Logical; disable censoring.
#' @param K Positive maximum number of follow-up events, including the final
#'   terminal or administrative event.
#' @param limit_event_A,limit_event_L Maximum number of treatment and covariate
#'   events. Non-negative integers or Inf.
#' @param aggregate_parameters Logical; reuse first-event coefficients.
#' @param static_intervention_stop Reserved; must be NULL.
#' @param time_varying_covariate_cumulative_effect Logical; increment L at
#'   covariate events rather than draw a new binary covariate.
#' @return A list containing one-row-per-subject baseline_data and event-level
#'   timevarying_data with columns id, time, event, A, and L.
#' @examples
#' simulate_continuous_time_data(10)
#' @export
simulate_continuous_time_data <- function(n,
                                          effects = list(
                                              alpha_A_0 = list(intercept = 0, age = 0.002),
                                              alpha_A_1 = list(intercept = 0.3, age = 0.002, L = -0.07, time = 0),
                                              alpha_A_2 = list(intercept = 0.3, age = 0.002, L = -0.07, time = 0),
                                              beta_l_1 = list(A = -0.2, age = 0.015),
                                              beta_l_2 = list(A = -0.2, age = 0.015),
                                              beta_c_1 = list(A = 0, L = 0, age = 0),
                                              beta_c_2 = list(A = 0, L = 0, age = 0),
                                              beta_c_3 = list(A = 0, L = 0, age = 0),
                                              beta_y_1 = list(A = -0.15, L = 0.02, age = 0.025),
                                              beta_y_2 = list(A = -0.15, L = 0.02, age = 0.025),
                                              beta_y_3 = list(A = -0.15, L = 0.02, age = 0.025),
                                              beta_d_1 = list(A = -1.2, L = 0.4, age = 0.015),
                                              beta_d_2 = list(A = -1.2, L = 0.4, age = 0.015),
                                              beta_d_3 = list(A = -1.2, L = 0.4, age = 0.015)
                                          ),
                                          static_intervention = NULL,
                                          static_intervention_baseline = 1,
                                          baseline_rate_list = list(
                                              A = 0.005,
                                              L = 0.001,
                                              C = 0.00005,
                                              Y = 0.0001,
                                              D = 0.00015
                                          ),
                                          max_fup = 900,
                                          visitation_interval = 360,
                                          visitation_sd = 5,
                                          discretize_age = FALSE,
                                          no_competing_events = FALSE,
                                          uncensored = FALSE,
                                          K = 3,
                                          limit_event_A = 1,
                                          aggregate_parameters = FALSE,
                                          static_intervention_stop = NULL,
                                          time_varying_covariate_cumulative_effect = TRUE,
                                          limit_event_L = 1) {
    positive_integer <- function(value, name) {
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value < 1 || value != floor(value) || value > .Machine$integer.max) {
            stop(name, " must be a positive integer.")
        }
    }
    nonnegative_number <- function(value, name, integer = FALSE,
                                   allow_infinite = FALSE) {
        if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
            (!allow_infinite && !is.finite(value)) || value < 0 ||
            (integer && is.finite(value) && value != floor(value))) {
            stop(name, " must be a non-negative ",
                 if (integer) "integer" else "number",
                 if (allow_infinite) " or Inf" else "", ".")
        }
    }
    positive_integer(n, "n")
    positive_integer(K, "K")
    nonnegative_number(limit_event_A, "limit_event_A", TRUE, TRUE)
    nonnegative_number(limit_event_L, "limit_event_L", TRUE, TRUE)
    nonnegative_number(max_fup, "max_fup")
    if (max_fup == 0) stop("max_fup must be positive.")
    nonnegative_number(visitation_interval, "visitation_interval")
    nonnegative_number(visitation_sd, "visitation_sd")
    flags <- list(discretize_age = discretize_age,
                  no_competing_events = no_competing_events,
                  uncensored = uncensored,
                  aggregate_parameters = aggregate_parameters,
                  time_varying_covariate_cumulative_effect =
                      time_varying_covariate_cumulative_effect)
    for (name in names(flags)) {
        value <- flags[[name]]
        if (!is.logical(value) || length(value) != 1L || is.na(value)) {
            stop(name, " must be TRUE or FALSE.")
        }
    }
    interventions <- list(static_intervention = static_intervention,
                          static_intervention_baseline = static_intervention_baseline)
    for (name in names(interventions)) {
        value <- interventions[[name]]
        if (!is.null(value) &&
            (!is.numeric(value) || length(value) != 1L || is.na(value) ||
             !value %in% c(0, 1))) {
            stop(name, " must be 0, 1, or NULL.")
        }
    }
    if (!is.null(static_intervention_stop)) {
        stop("static_intervention_stop is not supported.")
    }
    event_names <- c("A", "L", "C", "Y", "D")
    if (!is.list(baseline_rate_list) ||
        !setequal(names(baseline_rate_list), event_names) ||
        length(baseline_rate_list) != length(event_names) ||
        !all(vapply(baseline_rate_list, function(value) {
            is.numeric(value) && length(value) == 1L &&
                is.finite(value) && value >= 0
        }, logical(1)))) {
        stop("baseline_rate_list must contain finite non-negative rates named A, L, C, Y, and D.")
    }
    if (!is.list(effects) || !length(effects) || is.null(names(effects)) ||
        anyDuplicated(names(effects))) {
        stop("effects must be a named list of coefficient lists.")
    }
    for (name in names(effects)) {
        coefficients <- effects[[name]]
        if (!is.list(coefficients) || is.null(names(coefficients)) ||
            anyDuplicated(names(coefficients)) ||
            !all(vapply(coefficients, function(value) {
                is.numeric(value) && length(value) == 1L && is.finite(value)
            }, logical(1)))) {
            stop("effects$", name, " must be a named list of finite coefficients.")
        }
    }
    require_coefficients <- function(name, variables) {
        if (!all(variables %in% names(effects[[name]]))) {
            stop("effects$", name, " must specify ", paste(variables, collapse = ", "),
                 ". Use aggregate_parameters = TRUE to reuse the first event's effects.")
        }
    }
    if (aggregate_parameters && K > 1L) {
        for (prefix in c("alpha_A_", "beta_l_", "beta_c_", "beta_y_", "beta_d_")) {
            for (j in seq.int(2L, K)) {
                effects[[paste0(prefix, j)]] <- effects[[paste0(prefix, 1L)]]
            }
        }
    }
    for (j in seq_len(K)) {
        require_coefficients(paste0("beta_y_", j), c("A", "L", "age"))
        if (!uncensored && is.null(static_intervention)) {
            require_coefficients(paste0("beta_c_", j), c("A", "age"))
        }
        if (!no_competing_events) {
            require_coefficients(paste0("beta_d_", j), c("A", "L", "age"))
        }
        if (j < K) {
            require_coefficients(paste0("beta_l_", j), c("A", "age"))
            if (is.null(static_intervention)) {
                require_coefficients(paste0("alpha_A_", j),
                                     c("intercept", "L", "time", "age"))
            }
        }
    }
    if (is.null(static_intervention_baseline) && is.null(static_intervention)) {
        require_coefficients("alpha_A_0", c("intercept", "age"))
    }
    exponential_time <- function(n, rate, eta) {
        u <- stats::runif(n)
        if (rate == 0) rep(Inf, n) else -log(u) / (rate * exp(eta))
    }
    id <- time <- event <- age <- L <- A <- L_0 <- A_0 <- entrytime <-
        n_A_events <- n_L_events <- new_A <- NULL
    if (!is.null(static_intervention)) {
        static_intervention_baseline <- static_intervention
    }
    if (!discretize_age) {
        age <- stats::runif(n, 40, 90)
    } else {
        age <- sample(c(50, 70), n, replace = TRUE)
    }
    # Start each subject at time zero with L = 0 and an initial treatment.
    pop <- data.table(
        id = 1:n,
        age = age,
        L = as.numeric(rep(0, n)),
        A = as.numeric(rep(NA, n)),
        time = numeric(n),
        event = rep("0", n)
    )
    pop[, L_0 := 0]

    if (is.null(static_intervention_baseline)) {
        pop[, A_0 := stats::rbinom(n, 1, stats::plogis(effects$alpha_A_0$intercept + effects$alpha_A_0$age * age))]
    } else if (static_intervention_baseline %in% c(0, 1)) {
        pop[, A_0 := static_intervention_baseline]
    } else {
        stop("Intervention must be 0, 1, or NULL")
    }
    pop[, L := L_0]
    pop[, A := A_0]

    people_atrisk <- pop[, data.table::data.table(id, entrytime = time, age, L_0, A_0, A, L)]
    if (!is.null(static_intervention)) {
        uncensored <- TRUE
    }
    # Keep nonterminal records separately from the final record for each person.
    fup_info <- NULL
    has_terminal <- NULL

    j <- 1
    people_atrisk[, n_A_events := 0]
    people_atrisk[, n_L_events := 0]
    while (j <= K && nrow(people_atrisk) > 0) {
        # Draw competing waiting times for the next event. Treatment visits
        # following a treatment event are scheduled; other treatment visits
        # use an exponential clock. The last iteration allows only terminal events.
        if (j < K) {
            if (j == 1) {
                treatment_event <- rep(TRUE, nrow(people_atrisk))
            } else {
                treatment_event <- people_atrisk$event == "A"
            }

            max_event_reached_A <- people_atrisk$n_A_events >= limit_event_A
            a_time <- rep(NA, nrow(people_atrisk))
            a_time[max_event_reached_A] <- Inf
            a_time[treatment_event & !max_event_reached_A] <- pmax(0, visitation_interval +
                stats::rnorm(nrow(people_atrisk[treatment_event & !max_event_reached_A]), 0, visitation_sd))
            a_time[!treatment_event & !max_event_reached_A] <- exponential_time(
                n = nrow(people_atrisk[!treatment_event & !max_event_reached_A]),
                rate = baseline_rate_list$A,
                eta = 0
            )

            max_event_reached_L <- people_atrisk$n_L_events >= limit_event_L
            l_time <- rep(Inf, nrow(people_atrisk))
            l_time[!max_event_reached_L] <- exponential_time(
                n = nrow(people_atrisk[!max_event_reached_L]),
                rate = baseline_rate_list$L,
                eta = effects[[paste0("beta_l_", j)]]$A * people_atrisk[!max_event_reached_L, A] +
                    effects[[paste0("beta_l_", j)]]$age * people_atrisk[!max_event_reached_L, age]
            )
        } else {
            a_time <- rep(max_fup + 1, nrow(people_atrisk))
            l_time <- rep(max_fup + 1, nrow(people_atrisk))
        }
        if (!uncensored) {
            c_time <- exponential_time(
                n = nrow(people_atrisk),
                rate = baseline_rate_list$C,
                eta = effects[[paste0("beta_c_", j)]]$A * people_atrisk$A +
                    effects[[paste0("beta_c_", j)]]$age * people_atrisk$age
            )
        } else {
            c_time <- rep(max_fup + 1, nrow(people_atrisk))
        }
        y_time <- exponential_time(
            n = nrow(people_atrisk),
            rate = baseline_rate_list$Y,
            eta = effects[[paste0("beta_y_", j)]]$A * people_atrisk$A +
                effects[[paste0("beta_y_", j)]]$L * people_atrisk$L +
                effects[[paste0("beta_y_", j)]]$age * people_atrisk$age
        )
        if (!no_competing_events) {
            d_time <- exponential_time(
                n = nrow(people_atrisk),
                rate = baseline_rate_list$D,
                eta = effects[[paste0("beta_d_", j)]]$A * people_atrisk$A +
                    effects[[paste0("beta_d_", j)]]$L * people_atrisk$L +
                    effects[[paste0("beta_d_", j)]]$age * people_atrisk$age
            )
        } else {
            d_time <- rep(max_fup + 1, nrow(people_atrisk))
        }

        waiting_times <- cbind(a_time, l_time, c_time, y_time, d_time)
        next_event <- max.col(-waiting_times, ties.method = "first")
        people_atrisk[, event := factor(next_event, levels = 1:5,
                                        labels = c("A", "L", "C", "Y", "D"))]
        # One time unit separates events; truncate at the follow-up limit.
        people_atrisk[, time := waiting_times[cbind(seq_len(.N), next_event)] + entrytime + 1]
        people_atrisk[time > max_fup, event := "tauend"]
        people_atrisk[time > max_fup, time := max_fup]
        is_terminal <- !(people_atrisk$event %in% c("A", "L"))

        # Remove terminal events before updating treatment and covariate histories.
        has_terminal <- rbind(has_terminal, people_atrisk[is_terminal, data.table::data.table(id,
                                                                                              time = time,
                                                                                              event = event,
                                                                                              age,
                                                                                              L_0,
                                                                                              A_0,
                                                                                              A,
                                                                                              L
                                                                                              )])
        people_atrisk <- people_atrisk[!is_terminal]
        if (!is.null(static_intervention)) {
            people_atrisk[event == "A", new_A := static_intervention]
        } else {
            people_atrisk[event == "A", new_A := stats::rbinom(
                                                            .N, 1,
                                                            stats::plogis(effects[[paste0("alpha_A_", j)]]$intercept +
                                                                  effects[[paste0("alpha_A_", j)]]$L * L +
                                                                  effects[[paste0("alpha_A_", j)]]$time * time +
                                                                  effects[[paste0("alpha_A_", j)]]$age * age)
                                                        )]
        }
        people_atrisk[event == "A", n_A_events := n_A_events + 1]
        if (!time_varying_covariate_cumulative_effect) {
            people_atrisk[event == "L", L := stats::rbinom(.N, 1, 0.5)]
        } else {
            people_atrisk[event == "L", L := L + 1]
        }
        people_atrisk[event == "L", n_L_events := n_L_events + 1]
        people_atrisk[event == "A", A := new_A]

        fup_info <- rbind(fup_info, people_atrisk[, names(pop), with = FALSE], fill = TRUE)
        people_atrisk[, entrytime := time]
        j <- j + 1
    }
    pop <- rbind(has_terminal, fup_info)
    setkey(pop, id, time, event)
    baseline_vars <- c("age", "A_0", "L_0")
    baseline_data <- pop[, c("id", baseline_vars), with = FALSE]
    baseline_data <- baseline_data[!duplicated(baseline_data$id)]
    timevarying_data <- pop[, setdiff(names(pop), baseline_vars), with = FALSE]
    list(baseline_data = baseline_data, timevarying_data = timevarying_data)
}
