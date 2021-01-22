# https://github.com/jpwrobinson/funk/blob/master/R/themeSleek_function.R


#' A gg theme Function
#'
#' This function is for publication quality ggplots. It is a copy of Sean Anderson's ggsidekick which is no longer supported, and Cameron Freshwater added argument for top, bottom, middle for multipanel functionality.
#' @param 
#' @keywords ggplot
#' @export
#' @examples
#' theme_sleek()


theme_sleek <- function(base_size = 11, base_family = "", position = "standard") {
  half_line <- base_size/2
  q <- theme_light(base_size = 11, base_family = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "grey30"),
      strip.text.y = element_text(colour = "grey30"),
      axis.text = element_text(colour = "grey30"),
      axis.title = element_text(colour = "grey30"),
      legend.title = element_text(colour = "grey30", size = rel(1.1)),
      panel.border = element_rect(fill = NA, colour = "grey70", size = 1),
      legend.key.size = unit(1, "lines"),
      legend.text = element_text(size = rel(1), colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(colour = "grey30", size = rel(1)),
      plot.subtitle = element_text(colour = "grey30", size = rel(.85))
    )
  if (position == "bottom") {
    q <- q + theme(strip.background = element_blank(),
                   strip.text.x = element_blank(),
                   strip.text.y = element_blank(),
                   axis.text.y = element_text(size = 0.85 * axisSize),
                   axis.text.x = element_text(size = 0.85 * axisSize),
                   axis.title = element_text(size = axisSize)
    )
  }
  if (position == "top") {
    q <- q + theme(strip.text = element_text(size = axisSize),
                   axis.text.y = element_text(size = 0.85 * axisSize),
                   axis.text.x = element_blank(),
                   axis.title.y = element_text(size = axisSize),
                   axis.title.x = element_blank()
    )
  }
  if (position == "topWithX") {
    q <- q + theme(strip.text = element_text(size = axisSize),
                   axis.text.y = element_text(size = 0.85 * axisSize),
                   axis.text.x = element_text(size = 0.85 * axisSize),
                   axis.title.y = element_text(size = axisSize),
                   axis.title.x = element_blank()
    )
  }
  if (position == "mid") {
    q <- q + theme(strip.background = element_blank(),
                   strip.text.x = element_blank(),
                   strip.text.y = element_blank(),
                   axis.text.y = element_text(size = 0.85 * axisSize),
                   axis.text.x = element_blank(),
                   axis.title.y = element_text(size = axisSize),
                   axis.title.x = element_blank()
    )
  }
  return(q)
}