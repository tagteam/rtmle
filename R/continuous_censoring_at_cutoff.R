# Predict survival just after a cutoff on the absolute follow-up clock.
continuous_censoring_at_cutoff <- function(nuisance, k, pooled, rows, cutoff) {
    if (!any(nuisance$prepared_data_object$info$is_censored)) return(rep(1, length(rows)))

    models <- nuisance$censoring_prediction
    specification <- models$specification
    endpoint <- paste0(if (pooled) "time_pooled_" else "time_", k)
    start <- paste0("time_", k - 1L)
    call_learner <- function(args) do.call(specification$learner$fun, args)

    if (specification$marginal || specification$time_varying) {
        prediction <- nuisance$data[rows]
        data.table::set(prediction, j = endpoint, value = rep_len(cutoff, length(rows)))
    }
    if (specification$time_varying) {
        return(continuous_predict_censoring(models$shared, prediction, endpoint, start,
                                            specification, call_learner, left_limit = FALSE))
    }
    if (specification$marginal) {
        model <- models$marginal
        training <- nuisance$prepared_data_object$data_marginal_censoring
    } else {
        model <- models$conditional[[endpoint]]
        training <- model$training
        prediction <- training[match(rows, model$rows)]
        data.table::set(prediction, j = endpoint, value = cutoff - nuisance$data[[start]][rows])
        start <- NULL
    }
    result <- call_learner(c(list(character_formula = model$formula, data = training,
        intervened_data = prediction, time_variable = endpoint, time_ref = start,
        left_limit = FALSE, reuse_fit = model$fit, save_fitted_objects = TRUE), specification$learner$args))
    continuous_nuisance_predictions(result, length(rows), "Censoring", positive = TRUE)
}
