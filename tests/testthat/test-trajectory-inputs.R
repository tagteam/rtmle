test_that("trajectory inputs report missing identifiers, grids, and treatment histories", {
    x <- rtmle_init(0:2, name_id = "id", name_outcome = "Y")
    expect_error(plot_treatment_trajectories(list()), "class 'rtmle'")
    expect_error(plot_treatment_trajectories(x, NA_character_), "non-empty character")
    expect_error(plot_treatment_trajectories(x), "Cannot determine treatment variables")
    expect_error(plot_treatment_trajectories(x, "A"), "no subject ids")
    x <- add_baseline_data(x, data.frame(id = 1, age = 60))
    expect_error(plot_treatment_trajectories(x, "A"), "No treatment data")
    x$long_data$timevar_data <- list(A = data.frame(id = 1, value = 1))
    expect_error(plot_treatment_trajectories(x, "A"), "must contain either")
    x$time_grid_scale <- c(0, 0)
    expect_error(plot_treatment_trajectories(x, "A"), "increasing finite")
    x$names$id <- ""
    expect_error(plot_treatment_trajectories(x, "A"), "valid subject identifier")
})

test_that("calendar treatment histories use each participant's follow-up start", {
    x <- rtmle_init(c(0, 4), name_id = "id", name_outcome = "Y")
    start <- as.Date(c("2020-01-01", "2020-02-01"))
    x <- add_baseline_data(x, data.frame(id = 1:2, entry = start))
    x <- add_long_data(x, outcome_data = data.frame(id = 1:2, date = start + 4),
                       timevar_data = list(A = data.frame(
                           id = 1:2, date = start + c(1, 2), value = c(TRUE, FALSE))))
    expect_error(plot_treatment_trajectories(x, "A"), "subject-specific start_followup_date")
    missing_baseline <- x
    missing_baseline$data <- NULL
    expect_error(plot_treatment_trajectories(missing_baseline, "A"),
                 "subject-specific start_followup_date")
    x$names$start_followup_date <- "entry"
    p <- plot_treatment_trajectories(x)
    expect_equal(p$layers[[2L]]$data$time, c(1, 2))
    expect_setequal(as.character(p$layers[[2L]]$data$state), c("A=1", "A=0"))
    x$data$baseline_data$entry[2L] <- NA
    expect_error(plot_treatment_trajectories(x), "contain missing values")
})

test_that("wide and long trajectory tables handle missing subjects and zero follow-up", {
    x <- rtmle_init(0:2, name_id = "id", name_outcome = "Y")
    x <- add_baseline_data(x, data.frame(id = 1:3, age = 53:55))
    # Older objects can store a single combined wide table.
    x$data$timevar_data <- data.frame(id = 1:2, A_0 = c(TRUE, FALSE), A_1 = c(NA, TRUE))
    p <- plot_treatment_trajectories(x)
    expect_equal(as.character(p$data$state), c("A=1", "A=Missing", "A=0", "A=1",
                                               "A=Missing", "A=Missing"))
    x$names$treatment_variables <- "A"
    x$long_data$timevar_data <- data.frame(id = 1:2, date = c(0, 2), value = c(0, 1))
    x$long_data$outcome_data <- data.frame(id = 1:3, time = c(0, 1, 2))
    p <- plot_treatment_trajectories(x)
    expect_setequal(as.character(unique(p$data$panel)), c("Continuous time", "Discrete time"))
    expect_false(1L %in% p$data$subject)
    subject_ends <- p$data$end[p$data$subject == 2]
    expect_gt(length(subject_ends), 0L)
    expect_true(all(subject_ends <= 1))
    x$long_data$outcome_data$time <- 0
    expect_error(plot_treatment_trajectories(x), "do not contain any plot-able trajectories")
})

test_that("trajectory scales handle many states and entirely missing histories", {
    x <- rtmle_init(0:1, name_id = "id", name_outcome = "Y")
    x$data$timevar_data <- list(A = data.frame(id = 1:10, A_0 = 1:10))
    p <- plot_treatment_trajectories(x)
    expect_length(unique(ggplot2::ggplot_build(p)$data[[2L]]$colour), 10L)
    # R 3.5 uses rainbow() because grDevices has no hcl.colors().
    legacy_plot <- with_test_dependencies(plot_treatment_trajectories,
        exists = function(x, envir, inherits) {
            if (x == "hcl.colors" && identical(envir, asNamespace("grDevices"))) return(FALSE)
            base::exists(x, envir = envir, inherits = inherits)
        })
    palette <- ggplot2::ggplot_build(legacy_plot(x))$data[[2L]]$colour
    expect_setequal(palette, grDevices::rainbow(10, start = 0, end = .85))
    x$data$timevar_data$A$A_0 <- NA_integer_
    p <- plot_treatment_trajectories(x)
    expect_identical(levels(p$data$state), "A=Missing")
    expect_equal(unique(ggplot2::ggplot_build(p)$data[[2L]]$colour), "#BDBDBD")
    names(x$data$timevar_data$A)[2L] <- "A_5"
    expect_error(plot_treatment_trajectories(x), "do not contain any plot-able trajectories")
})

test_that("trajectory detection handles register histories before regimes are defined", {
    x <- rtmle_init(0:2, name_id = "id", name_outcome = "Y")
    x$long_data$timevar_data <- list(
        exposure = data.frame(id = 1, start_date = 0, end_date = 2),
        assay = data.frame(id = 1, date = 0:2, value = c(20, 50, 80)),
        visit = data.frame(id = 1, date = 1),
        unavailable = NULL
    )
    p <- plot_treatment_trajectories(x)
    expect_setequal(levels(p$data$state), c("exposure=0", "exposure=1"))
    expect_error(plot_treatment_trajectories(x, "unknown"), "No treatment data")
})

test_that("legacy trajectory objects fall back to available dates and grid labels", {
    x <- rtmle_init(0:1, name_id = "id", name_outcome = "Y")
    x$long_data$timevar_data <- list(A = data.frame(id = 1, date = 0, value = TRUE))
    x$names$time <- "visit"
    x$long_data$outcome_data <- data.frame(id = 1, time = .5)
    x$long_data$censored_data <- data.frame(id = 1)
    x$long_data$competing_data <- data.frame(id = integer(), time = numeric())
    x$time_grid_scale <- NULL
    x$time_grid <- NULL
    x$time_grid_labels <- NULL
    x$time_unit <- "Year"
    p <- plot_treatment_trajectories(x)
    expect_equal(p$data$end, .5)
    expect_equal(p$labels$x, "Year")
    expect_equal(p$scales$get_scales("x")$labels, 0:1)
    x$long_data$outcome_data$time <- NA_real_
    p <- plot_treatment_trajectories(x)
    expect_equal(p$data$end, 1)
    x$long_data$timevar_data$A$id <- NULL
    expect_error(plot_treatment_trajectories(x, "A"), "no.*plot-able trajectories")
})

test_that("wide trajectories use prepared data when a historical table lacks its id", {
    x <- rtmle_init(0:1, name_id = "id", name_outcome = "Y")
    x$data$timevar_data <- list(A = data.frame(A_0 = 5))
    x$prepared_data <- data.frame(id = 1, A_0 = 1L)
    p <- plot_treatment_trajectories(x, "A")
    expect_equal(as.character(p$data$state), "A=1")
})
