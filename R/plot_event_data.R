#' @title Plot Event Data
#'
#' @param data Event data: a data frame containing an ID, Time, Delta, L0 and L column.
#'
#' @return A plot of the survival data
#' @export
#'
#' @examples
#' data <- sim_event_data(10)
#' plot_event_data(data)

plot_event_data <- function(data) {
    ID <- Time <- max_time <- Delta <- NULL
    data$ID <- as.factor(data$ID)

  data <- data %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(max_time = max(Time)) %>%
    dplyr::arrange(max_time) %>%
    dplyr::mutate(ID = factor(ID, levels = ID)) %>%
    dplyr::right_join(data, by = "ID") %>%
    dplyr::select(-max_time)

  plotdata <- rbind(data, data.frame("ID" = unique(data$ID), L0 = unique(data$L0),
                                      "Time" = 0, "Delta" = "start", L = 0, A = unique(data$A)))

  diff_events <- length(unique(plotdata$Delta))
  cols <- c("green4", "blue1", "orange1", "red2")
  shapess <- c(20, 20, 20 , 20)

  ggplot2::ggplot(plotdata) +
    ggplot2::geom_line(ggplot2::aes(x = Time, y = ID, group = ID), color = "grey60", size = 0.7) +
    ggplot2::geom_point(ggplot2::aes(x = Time, y = ID, shape = factor(Delta), color = factor(Delta)), size = 2.5, data = data, alpha = 0.8) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::scale_shape_manual(values = shapess[1:diff_events]) +
    ggplot2::scale_color_manual(values = cols[1:diff_events]) +
    ggplot2::labs(
      title = "Event Data",
      x = "Time",
      y = "Patient ID",
      shape = "Event Type",
      color = "Event Type") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      axis.text.y = ggplot2::element_blank(),
      legend.position = "top"
    )
}
