### clear_rtmle_rows.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Nov 23 2026 (09:34) 
## Version: 
## Last-Updated: Nov 23 2026 (09:35) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
subset_rtmle_rows <- function(x, rows) {
    id_name <- x$names$id
    x$prepared_data <- x$prepared_data[rows, , drop = FALSE]
    if (NROW(x$followup) > 0L) {
        followup_rows <- match(x$prepared_data[[id_name]], x$followup[[id_name]])
        if (anyNA(followup_rows)) stop("Selected participants are missing follow-up data.")
        x$followup <- x$followup[followup_rows, , drop = FALSE]
    }
    for (regime_name in names(x$regimes)) {
        x$regimes[[regime_name]]$intervention_match <- NULL
    }
    x$censoring_use_regime <- NULL
    clear_rtmle_fit(x)
}

######################################################################
### clear_rtmle_rows.R ends here
