exposure_time <- function(data,
                          grid,
                          name,
                          id,
                          threshold = NA,
                          start_name = "date",
                          end_name = "previous_date"){

    gg = copy(grid)
    setnames(gg,c(start_name,end_name),c("end_interval","start_interval"))

    out = discretize_timevarying_exposure(data = data,
                                          grid = gg,
                                          name = name,
                                          point_exposure = FALSE,
                                          threshold = threshold,
                                          id = id)
    return(out)
}
