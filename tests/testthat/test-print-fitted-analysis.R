test_that("printing a fitted analysis identifies its data, learner and diagnostics", {
    data("rtmle_object", package = "rtmle", envir = environment())
    output <- paste(capture.output(print(rtmle_object)), collapse = "\n")
    expect_match(output, "Baseline data:[ ]+n=307, p=3")
    expect_match(output, "Censoring variable:[ ]+dropout")
    expect_match(output, "Competing risk:[ ]+death")
    expect_match(output, "Learner:[ ]+learn_glm")
    expect_match(output, "Fitted time horizons:[ ]+1, 2, 3, 4")
    expect_match(output, "constant_outcome_variables", fixed = TRUE)
    expect_match(output, "Results:", fixed = TRUE)
})

test_that("a stratifying variable named Participants has a distinct count column", {
    data("rtmle_object", package = "rtmle", envir = environment())
    baseline <- data.table::copy(rtmle_object$data$baseline_data)
    baseline$Participants <- factor(ifelse(baseline$age < 55, "younger", "older"))
    x <- add_baseline_data(rtmle_object, baseline)
    x$prepared_data <- NULL
    x <- prepare_rtmle_data(x, verbose = FALSE)
    groups <- stratify(x, Participants, time_horizon = 1, targets = "Treatment effect")
    expect_output(print(groups), "Participants[[:space:]]+Participant count")
    expect_equal(unname(attr(groups, "counts")), as.integer(table(baseline$Participants)))
})
