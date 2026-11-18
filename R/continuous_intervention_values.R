# Evaluate a regime at one event index using only history available at that
# index. The intervention table follows the same format as discrete regimes.
continuous_intervention_values <- function(data, regime, node) {
    continuous_intervention_policy(data, regime, node)$value
}
