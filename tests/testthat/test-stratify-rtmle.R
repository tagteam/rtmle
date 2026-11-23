with_baseline_group <- function(x, variable, values) {
    id <- x$names$id
    baseline <- data.table::copy(x$data$baseline_data)
    baseline[[variable]] <- values[match(baseline[[id]], x$prepared_data[[id]])]
    x <- add_baseline_data(x, baseline)
    x$prepared_data <- NULL
    prepare_rtmle_data(x)
}

test_that("stratification checks the variable and interval boundaries", {
    data("rtmle_object", package = "rtmle", envir = environment())
    unprepared <- rtmle_init(time_grid = 0:1, name_id = "id", name_outcome = "Y")
    expect_error(stratify.rtmle(list(), age), "x must be an rtmle object")
    expect_error(stratify.rtmle(unprepared, age), "prepare_rtmle_data")
    expect_error(stratify.rtmle(rtmle_object), "Supply a variable")
    expect_error(stratify.rtmle(rtmle_object, "unknown"), "must be a column")
    expect_error(stratify.rtmle(rtmle_object, c("age", "sex")), "must be a column")
    for (breaks in list(50, c(0, 50, 50), c(50, 0), c(0, NA), c("0", "50"))) {
        expect_error(stratify.rtmle(rtmle_object, age, breaks = breaks),
                     "strictly increasing breaks")
    }
    expect_error(stratify.rtmle(rtmle_object, age, breaks = c(-2, -1)),
                 "No participants fall into a stratum")
})

test_that("stratification requires a baseline covariate", {
    data("rtmle_object", package = "rtmle", envir = environment())
    expect_error(stratify.rtmle(rtmle_object, A_0), "baseline covariate")
    for (variable in c("changeSBP_1", "stroke_1", "id")) {
        expect_error(stratify.rtmle(rtmle_object, variable, breaks = c(-Inf, Inf)),
                     "baseline covariate")
    }
})

test_that("numeric sex indicator can be stratified without breaks", {
    data("rtmle_object", package = "rtmle", envir = environment())
    by_sex <- stratify.rtmle(rtmle_object, sex, time_horizon = 1:2)

    expect_identical(names(by_sex), c("0", "1"))
    expect_null(attr(by_sex, "breaks"))
    original_sex <- rtmle_object$prepared_data$sex
    for (stratum in names(by_sex)) {
        fit <- by_sex[[stratum]]
        rows <- original_sex == as.numeric(stratum)
        expect_identical(fit$prepared_data$sex, original_sex[rows])
        expect_identical(fit$prepared_data$id, rtmle_object$prepared_data$id[rows])
        expect_equal(sort(unique(fit$estimate$Main_analysis$Time_horizon)), 1:2)
    }

    x <- rtmle_object
    x$prepared_data$sex[1:2] <- c(NA_real_, NaN)
    expect_warning(
        by_sex <- stratify.rtmle(x, sex, time_horizon = 1, targets = "Treatment effect"),
        "Excluded 2 participant", fixed = TRUE
    )
    expect_identical(names(by_sex), c("0", "1"))
    expect_equal(sum(attr(by_sex, "counts")), length(original_sex) - 2L)
})

test_that("numeric strata follow observed values when breaks are NULL", {
    data("rtmle_object", package = "rtmle", envir = environment())
    group <- rep_len(c(10, -2, 3.5), nrow(rtmle_object$prepared_data))
    x <- with_baseline_group(rtmle_object, "group", group)
    fits <- stratify(x, group, breaks = NULL, time_horizon = 1, targets = "Treatment effect")

    expect_identical(names(fits), c("-2", "3.5", "10"))
    expect_null(attr(fits, "breaks"))
    for (value in names(fits)) {
        rows <- group == as.numeric(value)
        expect_identical(fits[[value]]$prepared_data$id, x$prepared_data$id[rows])
        expect_identical(fits[[value]]$prepared_data$group, group[rows])
    }
})

test_that("stratification warns below the minimum stratum size and still fits", {
    data("rtmle_object", package = "rtmle", envir = environment())
    labels <- c("small", "boundary", "large")
    sizes <- c(29L, 30L, nrow(rtmle_object$prepared_data) - 59L)
    group <- factor(rep(labels, sizes), levels = labels)
    x <- with_baseline_group(rtmle_object, "group", group)

    expect_error(stratify.rtmle(x, group, breaks = c(0, 1)), "Set breaks to NULL")
    expect_warning(
        fits <- stratify.rtmle(x, group, time_horizon = 1, targets = "Treatment effect"),
        "Small strata.*fewer than 30 participants.*small \\(n = 29\\)"
    )
    expect_identical(unname(attr(fits, "counts")), sizes)
    expect_gt(nrow(fits[["small"]]$estimate$Main_analysis), 0L)

    group[which(group == "large")[1L]] <- "small"
    x <- with_baseline_group(rtmle_object, "group", group)
    expect_silent(stratify.rtmle(x, group, time_horizon = 1,
                                targets = "Treatment effect"))
})

test_that("the number-of-strata warning counts only nonempty groups", {
    data("rtmle_object", package = "rtmle", envir = environment())
    group <- factor(rep_len(letters[1:6], nrow(rtmle_object$prepared_data)),
                    levels = letters[1:7])
    x <- with_baseline_group(rtmle_object, "group", group)
    expect_warning(
        fits <- stratify.rtmle(x, group, time_horizon = 1, targets = "Treatment effect"),
        "6 nonempty strata"
    )
    expect_identical(names(fits), letters[1:6])

    group[group == "f"] <- "e"
    x <- with_baseline_group(rtmle_object, "group", group)
    expect_silent(fits <- stratify.rtmle(x, group, time_horizon = 1,
                                        targets = "Treatment effect"))
    expect_identical(names(fits), letters[1:5])
})

test_that("stratified summaries and plots respect targets and reference regimes", {
    data("rtmle_object", package = "rtmle", envir = environment())
    fits <- stratify.rtmle(rtmle_object, age, breaks = c(-Inf, 50, Inf),
                          time_horizon = 1:2)
    target <- "Treatment effect"
    reference <- list("Treatment effect" = "Never_A")
    estimates <- summary(fits, targets = target, reference = reference, digits = 2)
    expect_identical(unique(as.character(estimates$Target)), target)
    expect_identical(names(estimates)[1L], "age")
    expect_identical(levels(estimates$age), names(fits))
    for (stratum in names(fits)) {
        separate <- summary(fits[[stratum]], targets = target, reference = reference,
                            digits = 2)
        combined <- estimates[estimates$age == stratum, names(separate),
                              with = FALSE]
        expect_equal(combined, separate)
    }

    output <- capture.output(printed <- withVisible(print(fits, targets = target)))
    expect_false(printed$visible)
    expect_identical(printed$value, fits)
    expect_true(any(grepl("^[[:space:]]*age[[:space:]]+Participants", output)))
    expect_true(any(grepl("^[[:space:]]*age[[:space:]]+Target", output)))
    expect_false(any(grepl("Stratum", output, fixed = TRUE)))
    expect_false(any(grepl("Treatment contrast", output, fixed = TRUE)))

    p <- plot(fits, targets = target, reference = reference,
              target_parameter = "Risk_difference", conf_int = FALSE)
    differences <- estimates[estimates$Target_parameter == "Risk_difference"]
    expect_equal(p$data$Estimate, differences$Estimate)
    expect_equal(p$data$Series, rep("Always_A vs Never_A", nrow(differences)))
    expect_false(any(vapply(p$layers, function(layer) {
        inherits(layer$geom, "GeomRibbon")
    }, logical(1))))
    panels <- ggplot2::ggplot_build(p)$layout$layout
    expect_identical(as.character(panels$Stratum), names(fits))
    expect_identical(unique(as.character(panels$Target)), target)

    ratios <- ggplot2::autoplot(fits, targets = target, reference = reference,
                                target_parameter = "Risk_ratio")
    expect_equal(ratios$data$Estimate,
                 estimates$Estimate[estimates$Target_parameter == "Risk_ratio"])
    risks <- ggplot2::autoplot(fits, targets = target)
    baseline <- risks$data[risks$data$Time_horizon == 0]
    expect_equal(baseline$Estimate, rep(0, 2 * length(fits)))
    expect_true(any(vapply(risks$layers, function(layer) {
        inherits(layer$geom, "GeomRibbon")
    }, logical(1))))
    expect_s3_class(ggplot2::ggplotGrob(risks), "gtable")

    layout <- ggplot2::ggplot_build(plot(fits))$layout$layout
    expect_identical(as.character(layout$Stratum), rep(names(fits), each = 2L))
    expect_equal(layout$ROW, c(1L, 1L, 2L, 2L))
    expect_equal(layout$COL, c(1L, 2L, 1L, 2L))
})

test_that("a stratifying variable named Target does not replace the target column", {
    data("rtmle_object", package = "rtmle", envir = environment())
    x <- with_baseline_group(rtmle_object, "Target", factor(rtmle_object$prepared_data$sex))
    fits <- stratify(x, Target, time_horizon = 1, targets = "Treatment effect")

    estimates <- summary(fits)
    expect_identical(names(estimates)[1L], "Target.1")
    expect_identical(levels(estimates$Target.1), names(fits))
    expect_identical(unique(as.character(estimates$Target)), "Treatment effect")
    expect_s3_class(ggplot2::ggplotGrob(plot(fits)), "gtable")
})
