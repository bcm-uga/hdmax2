##' Summary plot for ACME
##'
##' This function draw a summary plot of ACME (average causal mediation effect)
##'
##' @param ACME the table of ACME from the wrap_mediation() function
##' @return
##' Summary plot for ACME
##'
##'
##' @export
##' @author Basile Jumentier
##' @examples
##'
##' # see wrap_mediation example
##'
##' @import ggplot2
##'

plot_summary_ACME <- function(ACME) {
  
  # for check problem
  res_med$ACME = ACME
  p <- ggplot(res_med$ACME, aes(est, stats::reorder(feat, est), color = pval <= 0.05, shape = pval <= 0.05)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = CI_2.5, xmax = CI_97.5)) +
    geom_point(size = 0.8) +
    theme_bw() +
    xlab("ACME (Average Causal Mediation Effect)") +
    ylab("CpG") +
    theme(panel.border = element_blank(),
          panel.spacing = unit(0.01, "lines"),
          axis.ticks = element_blank()) +
    scale_color_manual(values = c("black", "red"))
  
  print(p)
}