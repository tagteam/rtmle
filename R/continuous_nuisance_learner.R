continuous_nuisance_learner <- function(learner) {
    if (is.function(learner)) learner <- list(name = "custom", fun = learner)
    parse_learners(learner)
}

is_ipcw_glm_learner <- function(learner) {
    identical(learner$fun, "learn_ipcw_glm") || identical(learner$fun, learn_ipcw_glm)
}
