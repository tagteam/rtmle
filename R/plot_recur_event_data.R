#' @title Plot Recurrent Event Data
#'
#' @param data Survival data: a data frame containing an ID, Time, Delta and L0 column. Possibly a L1 and A column.
#'
#' @return A plot of the survival data
#' @export
#'
#' @examples
#' data <- sim_recur_event_data(10)
#' plot_recur_event_data(data)

plot_recur_event_data <- function(data) {

  if("L1" %in% colnames(data)){
    plotdata <- rbind(data, data.frame("ID" = unique(data$ID), L0 = unique(data$L0),
                                       "Time" = 0, "Delta" = "start", L1 = 0, A = unique(data$A))) %>%
      arrange(Time)
  }
  else{
    plotdata <- rbind(data, data.frame("ID" = unique(data$ID), L0 = unique(data$L0),
                                       "Time" = 0, "Delta" = "start"))  %>%
      arrange(Time)
  }

  diff_events <- length(unique(plotdata$Delta))
  cols <- c("darkorange", "blue3", "lightgreen", "darkred", "pink")
  shapess <- c(8, 20, 20 , 8, 20)

  ggplot2::ggplot(plotdata) +
    ggplot2::geom_line(ggplot2::aes(x = Time, y = ID, group = ID), color = "grey60", size = 0.7) +  # Light gray lines for each ID
    ggplot2::geom_point(ggplot2::aes(x = Time, y = ID, shape = factor(Delta), color = factor(Delta)), size = 3, data = data) + # Use shapes and colors for Delta
    ggplot2::theme_minimal(base_size = 15) +  # Increase base font size for readability
    ggplot2::scale_shape_manual(values = shapess[1:diff_events]) +  # Customize shapes for different events
    ggplot2::scale_color_manual(values = cols[1:diff_events]) +  # Customize colors for events
    ggplot2::labs(
      title = "Survival Data",
      x = "Time",
      y = "Patient ID",
      shape = "Event Type",
      color = "Event Type"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),  # Center and bold title
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),  # Add space to y-axis title
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),  # Add space to x-axis title
      legend.position = "top"  # Place legend on top for easy viewing
    )
}
