continuous_example <- function() {
    list(
        baseline_data = data.frame(id = 1:4, A_0 = c(0, 0, 1, 1), L_0 = 0),
        timevarying_data = data.frame(id = 1:4, time = 2,
                                     event = c("tauend", "tauend", "Y", "tauend"),
                                     A = c(0, 0, 1, 1), L = 0)
    )
}

continuous_add_example_events <- function(x, data) {
    events <- data.table::as.data.table(data$timevarying_data)
    point_data <- function(mark, variable) {
        out <- data.table::copy(events[events$event == mark, c("id", "time"), with = FALSE])
        data.table::set(out, j = "value", value = events[[variable]][events$event == mark])
        out
    }
    add_continuous_data(x,
        outcome_data = events[events$event == "Y", c("id", "time"), with = FALSE],
        competing_data = events[events$event == "D", c("id", "time"), with = FALSE],
        censored_data = events[events$event %in% c("C", "tauend"), c("id", "time", "event"), with = FALSE],
        timevar_data = list(A = point_data("A", "A"), L = point_data("L", "L")),
        treatment_events = "A", covariate_events = "L")
}

continuous_add_example_data <- function(x, data) {
    x <- add_baseline_data(x, data$baseline_data)
    continuous_add_example_events(x, data)
}

continuous_analysis <- function(data = continuous_example(), horizons = 2,
                                regimes = c("Always_A", "Never_A"), ...) {
    x <- rtmle_init(c(0, horizons), name_id = "id", name_time = "time",
                    name_outcome = "Y", name_competing = "D",
                    name_censoring = "C", continuous_time = TRUE)
    x <- continuous_add_example_data(x, data)
    x <- prepare_rtmle_data(x, ...)
    x <- regime(x, "Always_A", intervention = 1, treatment_variables = "A")
    x <- regime(x, "Never_A", intervention = 0, treatment_variables = "A")
    x <- target(x, "risk", regimes = regimes)
    model_formula(x, verbose = FALSE)
}
