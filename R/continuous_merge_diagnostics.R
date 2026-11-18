continuous_merge_diagnostics <- function(existing, additional) {
    for (name in names(additional)) {
        if (!NROW(additional[[name]])) next
        existing[[name]] <- data.table::rbindlist(
            list(existing[[name]], additional[[name]]), fill = TRUE)
    }
    existing
}
