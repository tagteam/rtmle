test_that("pooled Cox censoring predictions use the appropriate clock at each cutoff", {
    # The shorter horizon pools at event 2; the longer horizon pools at event 3.
    # Prior event times differ by person, as do the baseline censoring risks.
    w <- data.table::data.table(
        id = 1:8, A_0 = 1, L_0 = c(0, 0, 1, 1, 0, 1, 1, 0),
        time_0 = 0, event_0 = "A",
        time_1 = c(20, 80, 25, 40, 5, 10, 15, 35), event_1 = "A", A_1 = 1,
        time_2 = c(50, 110, 55, 70, 35, 50, 45, 65),
        event_2 = c("A", "A", "A", "A", "C", "Y", "C", "Y"), A_2 = 1,
        time_3 = c(95, 120, 90, 110, NA, NA, NA, NA),
        event_3 = c("Y", "C", "C", "Y", NA, NA, NA, NA), A_3 = 1,
        time_pooled_2 = c(95, 120, 90, 110, 35, 50, 45, 65),
        event_pooled_2 = c("Y", "C", "C", "Y", "C", "Y", "C", "Y"),
        time_pooled_3 = c(95, 120, 90, 110, NA, NA, NA, NA),
        event_pooled_3 = c("Y", "C", "C", "Y", NA, NA, NA, NA))
    w[, `:=`(L_1 = L_0, L_2 = L_0, L_3 = L_0)]
    prepared <- list(
        wide_data = w,
        info = data.table::data.table(time_horizon = c(100, 720),
                                      last_event = c(2L, 3L), is_censored = TRUE),
        time_covariates = c("A", "L"), baseline_covariates = c("A_0", "L_0"),
        data_marginal_censoring = w[, .(id, A_0, L_0, time = time_pooled_2,
                                         event = event_pooled_2)])

    # Compute the Cox survival ratio independently of rtmle's learner.
    manual_survival <- function(train_time, train_event, train_L, new_L, start, end) {
        training <- data.frame(time = train_time, censored = train_event == "C", L = train_L)
        fit <- survival::coxph(survival::Surv(time, censored) ~ L, data = training)
        hazard <- survival::basehaz(fit, centered = FALSE)
        cumulative <- function(t, left_limit = FALSE) {
            index <- findInterval(t, hazard$time)
            if (left_limit) index <- index - as.integer(t %in% hazard$time)
            c(0, hazard$hazard)[index + 1L]
        }
        lp <- stats::predict(fit, newdata = data.frame(L = new_L), type = "lp", reference = "zero")
        as.numeric(exp(-(cumulative(end, TRUE) - cumulative(start)) * exp(lp)))
    }

    options <- list(Markov = "L")
    conditional <- rtmle:::continuous_propensity(prepared, model_options = options,
        censoring_learner = list(name = "coxph", fun = "learn_coxph", marginal_censoring = FALSE))
    marginal <- rtmle:::continuous_propensity(prepared, model_options = options,
        censoring_learner = list(name = "coxph", fun = "learn_coxph", marginal_censoring = TRUE))

    for (k in 2:3) {
        rows <- which(w[[paste0("event_", k - 1L)]] %in% c("A", "L"))
        endpoint <- w[[paste0("time_pooled_", k)]][rows]
        start <- w[[paste0("time_", k - 1L)]][rows]
        event <- w[[paste0("event_pooled_", k)]][rows]
        previous_L <- w[[paste0("L_", k - 1L)]][rows]

        # Conditional Cox is fitted on time since each person's prior event.
        elapsed <- endpoint - start
        conditional_expected <- manual_survival(elapsed, event, previous_L, previous_L,
                                                rep(0, length(rows)), elapsed)
        # Marginal Cox is fitted once on absolute terminal times, then the
        # absolute-time survival ratio is evaluated from start to endpoint.
        marginal_expected <- manual_survival(w$time_pooled_2, w$event_pooled_2, w$L_0,
                                             w$L_0[rows], start, endpoint)
        column <- paste0("survival_censoring_pooled_", k)
        expect_equal(as.numeric(conditional$data[[column]][rows]), conditional_expected, tolerance = 1e-12)
        expect_equal(as.numeric(marginal$data[[column]][rows]), marginal_expected, tolerance = 1e-12)
    }
})
