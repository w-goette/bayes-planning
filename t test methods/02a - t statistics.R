#' Customized t.test function
#'
#' Simple function for calculating either independent (`paired = FALSE`) or
#' dependent (`paired = TRUE`) samples t-test. Compared to the t.test() function
#' in base R, this function assumes equal variances, performs only a two-tailed
#' test, and returns just the t-statistic, p-value, and standardized mean
#' difference.
#'
#' @param x1 numeric vector of values corresponding to observations in sample 1
#' @param x2 numeric vector of values corresponding to observations in sample 2
#' @param paired logical indicator for whether the two samples are dependent or
#'   independent of each other. Default is TRUE.
#' @param mu continuous real numeric value representing a point-estimate null
#'   value to test against. Default is 0 (traditional nil hypothesis test).
#'
#' @returns a numeric vector containing the t-statistic, associated p-
#'   value, and standardized effect size
#'
t_test <- function(x1, x2, paired = TRUE, mu = 0) {

  # extract summary statistics
  nu <- ifelse(
    paired,
    length(x1) - 1,
    length(x1) + length(x2) - 2
  )

  diff <- ifelse(
    paired,
    mean(x2 - x1) - mu,
    mean(x1) - mean(x2) - mu
  )

  sd <- ifelse(
    paired,
    sd(x2 - x1),
    sqrt( ((length(x1) - 1) * var(x1) + (length(x2) - 1) * var(x2)) / nu )
  )

  se <- ifelse(
    paired,
    sd / sqrt( length(x1) ),
    sd * sqrt( 1 / length(x1) + 1 / length(x2) )
  )

  # calculate t-statistic
  t <- diff / se

  # obtain associated p-value
  p <- 2 * pt(-abs(t), nu)

  # calculate standardized effect size
  d <- ifelse(
    paired,
    diff / ( sd / sqrt(2 * (1 - cor(x1, x2))) ),
    diff / sd
  )

  # store results in a vector
  out <- c(t, p, d)

  # return result
  return( out )
}

#' Find smallest effect size detectable with specified power
#'
#' Borrowing from the frequentist equivalence literature, particularly TOSTs,
#' equivalence bounds can be defined in terms of minimally detectable effects
#' under varying power assumptions. This function utilizes traditional power
#' analysis formulas to identify the smallest effect size for which a sample
#' size has some pre-specified power to detect.
#'
#' @param n1 number of observations in sample 1
#' @param n2 number of observations in sample 2
#' @param paired logical indicator for whether the power calculation should be
#'   for a paired (TRUE) or independent (FALSE) samples t-test
#' @param power numeric value between 0 and 1 reflecting the desired power of
#'   the t-test. Default is 0.33
#' @param alpha numeric value between 0 and 1 reflecting the desired Type I
#'   error rate of the expected NHST. Default is 0.05
#' @param two_tailed logical indicator for whether the expected NHST is a test
#'   of a one- or two-tailed hypothesis. Default is TRUE
#'
#' @returns real numeric vector corresponding to the lower and upper bounds of
#'   the standardized effect size centered at zero that could be detected at the
#'   desired power level
#'
min_effect <- function(n1, n2, paired, power = 0.33, alpha = 0.05,
                       two_tailed = TRUE) {

  # calculate needed summary values
  nu <- ifelse(
    paired,
    n1 - 1,
    n1 + n2 - 2
  )

  hm <- 1 / mean( c(1/n1, 1/n2) )

  div <- ifelse(
    paired,
    1,
    2
  )

  cv <- ifelse(
    two_tailed,
    qt(1 - alpha / 2, nu),
    qt(1 - alpha, nu)
  )

  # apply solver to the power function
  out <- uniroot(function(d) {
    # compute power at given d
    c_power <- 1 - ( pt(cv, nu, sqrt(hm / div) * d) -
                       pt(-cv, nu, sqrt(hm / div) * d) )

    # find difference in computed and requested power
    c_power - power
  }, c(1e-6, 10))$root

  # return result
  return( c(-out, out) )
}

#' Define equivalence region for a practice effect
#'
#' prac_rope() is used to convert an a priori estimate of the practice effect
#' (i.e., the expected gain in raw score from prior exposure) into equivalence
#' bounds for the null standardized effect size. The resulting bounds can be
#' treated as a region of practical equivalence (ROPE) for any observed mean
#' difference between paired sample data when the expectation is that any such
#' difference is due to known/predictable practice effects.
#'
#' Put differently, the ROPE computed from this calculation represents an
#' informed prior over the null hypothesis. The bounds for the standardized
#' effect size reflect current uncertainty in the estimated practice effect
#' between samples. In the case that the posterior for a paired samples t-test
#' results in a standardized effect estimate that is entirely or largely within
#' this interval, the conclusion would be that the apparent mean difference
#' between the two testing sessions is consistent with practice effects rather
#' than a method effect. Posterior standardized effect estimates outside of the
#' interval would suggest that there were less than expected practice effects or
#' that there are additional causes of pre- versus post-test differences beyond
#' practice effects.
#'
#' @param prac_effect continuous real numeric value representing the expected
#'   difference in scores due to practice effects
#' @param x1 numeric vector of observations for sample 1
#' @param x2 numeric vector of observations for sample 2
#' @param ci continuous real numeric value between 0 and 1 indicating the
#'   confidence interval size desired for the implied standardized effect size
#'   (e.g., 0.90 will give the 90% confidence interval)
#' @param null_centered logical indicator for whether the practice effect will
#'   be included as a point estimate of the null (TRUE). When FALSE (default),
#'   the resulting bounds are centered around the effect size rather than zero.
#'
#' @returns numeric vector of length 2 consisting of the lower and upper bounds,
#'   respectively, of the standardized effect size
#'
prac_rope <- function(prac_effect, x1, x2, ci = 0.90, null_centered = FALSE) {

  # calculate needed values
  n <- length(x1)
  nu <- n - 1
  hn <- 1 / n
  sd <- sd(x1 - x2) / sqrt(2 * (1 - cor(x1, x2)))
  cv <- c(     (1 - ci) / 2,
           1 - (1 - ci) / 2 )

  # convert practice effect to standardized effect size
  d <- prac_effect / sd

  # calculate standard error of the effect size
  se <- sqrt( (1 / n + d^2 / (2 * n)) * 2 * (1 - cor(x1, x2)) )

  # find quantiles from the non-centered t-distribution
  ncp <- c(qt(cv[1], nu, sqrt(hn) * d),
           qt(cv[2], nu, sqrt(hn) * d))

  # calculate upper and lower bounds of the standardized effect size
  out <- c((d * !(null_centered)) + ncp[1] * se,
           (d * !(null_centered)) + ncp[2] * se)

  # return result
  return( out )
}

#' Scale between standardized practice effect and raw scale practice effect
#'
#' A convenience function for converting between values in a simulation-based
#' population generated on a z-scaled latent variable to the (optionally)
#' discrete raw scale values of a simulated sample. In a simulated test/re-test
#' design, the standardized effect size added to the re-test "scores" may need
#' to represent a practice effect, but the practice effect is typically going
#' to be represented as raw score differences rather than a standardized effect.
#' This function simply takes a desired raw score practice effect and converts
#' this to the appropriate standardized effect scale that can be used to generate
#' the z-scaled latent scores in the population.
#'
#' @param prac_effect continuous real numeric value representing he expected
#'   difference in scores due to practice effects
#' @param sd standard deviation of raw scores in the population
#' @param r12 test-retest reliability in the population
#'
#' @returns continuous real numeric value representing the standardized effect
#'   size associated with the provided raw score practice effect
#'
scale_prac <- function(prac_effect, sd, r12) {

  ( prac_effect / sqrt(2 * sd^2 - 2 * r12 * sd^2) ) * sqrt(2 * (1 - r12))

}
