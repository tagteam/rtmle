#' The function transforms the simulated event data containing time dependent covariates
#' into the tstart tstop format used by the coxph function.
#'
#' @title Transformation of Event Data
#'
#' @param data Event data: a data frame containing an ID, Time, Delta, L0 and L column.
#'
#' @return Event data in a tstart tstop format
#' @export
#'
#' @examples
#' data <- sim_event_data(10)
#' trans_int_data(data)

trans_int_data <- function(data) {

  k <- ID <- tstart <- tstop <- Time <- NULL

  data <- copy(data)

  # Creating a k variable
  data[, k:= stats::ave(ID, ID, FUN = seq_along)]
  max_k <- max(data$k)

  data_k <- list()
  data_k[[1]] <- data[data$k == 1,]

  data_k[[1]][, tstart := 0]
  data_k[[1]][, tstop := Time]

  for(i in 2:max_k){
    data_k[[i]] <- data[data$k == i,]

    data_k[[i]][, tstart := data_k[[i-1]][data_k[[i-1]]$ID %in% data_k[[i]]$ID,]$tstop]
    data_k[[i]][, tstop := Time]
  }

  res <- do.call(rbind, data_k)

  setkey(res, ID)

  return(res)
}
