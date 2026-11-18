continuous_nuisance_learner <- function(learner) {
    if (is.function(learner)) learner <- list(name = "custom", fun = learner)
    parse_learners(learner)
}
