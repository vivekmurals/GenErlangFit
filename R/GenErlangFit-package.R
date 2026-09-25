#' @keywords internal
#' @importFrom stats IQR approx density dgamma ecdf integrate optimize
#' @importFrom stats pgamma rexp rgamma setNames
#' @importFrom utils head tail modifyList
#' @importFrom ggplot2 ggplot aes geom_histogram geom_line geom_step
#' @importFrom ggplot2 geom_vline scale_fill_manual scale_color_manual
#' @importFrom ggplot2 labs theme theme_gray theme_minimal element_blank
#' @importFrom ggplot2 element_text after_stat
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

utils::globalVariables(c("x", "y", "Type", "Statistic", "density"))
