#' Generate simple summary plot
#'
#' Generic plotting function to visualize the prior and posterior distributions
#' of standardized effect sizes given input values. Also plotted are the 95%
#' credible intervals and median point estimate of the standardized effect size.
#'
#' @param t t-statistic computed from the sample
#' @param n1 integer value for the size of the first sample
#' @param n2 integer value for the size of the second sample
#' @param paired logical indicator for whether the samples are paired (true)
#'   or independent (false)
#' @param mu_delta hyperprior on the location of the standardized effect size
#' @param gamma hyperprior on the scale of the standardized effect size
#' @param kappa hyperprior on the degrees of freedom for the t-distribution prior
#'   on the standardized effect size
#' @param x_range numeric vector of length 2 providing the lower and upper bounds
#'   of possible standardized effect sizes to plot
#'
#' @returns the posterior and prior distributions under informed priors as a
#'   ggplot2 object
#'
point_plot <- function(t, n1, n2 = NULL, paired = TRUE, mu_delta, gamma, kappa,
                       x_range = c(-1, 1)) {

  # obtain point estimates of interest
  sum <- ci_t(t, n1, n2, paired, mu_delta, gamma, kappa)

  # create plot
  ggplot() +
    scale_x_continuous(
      limits = c(x_range[1], x_range[2]),
      breaks = seq(x_range[1], x_range[2], sum(abs(x_range)) / 8)
    ) +
    stat_function(
      fun = function(x) 1 / gamma * dt( (x - mu_delta) / gamma, kappa),
      linewidth = 1.15,
      linetype = "dotted",
      aes(color = "Prior")
    ) +
    stat_function(
      geom = "area",
      fun = function(x) posterior_t(x, t, n1, n2, paired, mu_delta, gamma, kappa),
      xlim = c(sum[1], sum[3]),
      alpha = 0.20,
      fill = RColorBrewer::brewer.pal(4, "Dark2")[1]
    ) +
    stat_function(
      fun = function(x) posterior_t(x, t, n1, n2, paired, mu_delta, gamma, kappa),
      linewidth = 1.5,
      aes(color = "Posterior")
    ) +
    annotate(
      "segment",
      x = sum,
      xend = sum,
      y = c(0, 0, 0),
      yend = posterior_t(sum, t, n1, n2, paired, mu_delta, gamma, kappa),
      linewidth = c(1.25, 1, 1.25),
      color = RColorBrewer::brewer.pal(4, "Dark2")[1],
      linetype = c("solid", "dashed", "solid")
    ) +
    annotate(
      "text",
      x = sum + c(-0.07, -1 * sign(sum[2]) * 0.07, 0.07),
      y = posterior_t(sum, t, n1, n2, paired, mu_delta, gamma, kappa) + 0.1,
      label = paste("delta ==", round(sum, 2)),
      parse = TRUE,
      size = 12/.pt
    ) +
    labs(
      y = "Density",
      x = expression(paste("Effect size (", delta, ")")),
      col = ""
    ) +
    scale_color_brewer(
      palette = "Dark2"
    ) +
    theme_classic() +
    theme(
      legend.text = element_text(family = "sans", size = 12),
      legend.position = c(0.9, 0.9),
      axis.line = element_line(linewidth = 1),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}
