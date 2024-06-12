#' This script contains lightly edited code from https://osf.io/bsp6z
#' Any credit or recognition of this code should be directed to the above source



#' Integral function for BF10 of d given informed prior on d
#'
#' The integral appears as the numerator of equation (5) of Gronau et al. (2020)
#' and is taken with respect to delta (standardized effect size). The integrand
#' consists of two different student t-distributions. The first is a non-centered
#' t-distribution reflecting the posterior of the obtained t-statistic in the
#' sample. The non-centered parameter (ncp) of the distribution is the product
#' of the squared root of the effect sample size and delta.
#'
#' The second t-distribution reflects the prior on delta given three
#' hyperparameters: mu_delta (location), gamma (scale), and kappa (degrees of
#' freedom). In short, mu_delta is the prior estimate of the standardized effect
#' size, gamma is the scale or dispersion of certainty for the point estimated
#' mu_delta, and kappa reflects overall confidence in the implied distribution
#' under these hyperparameters. As kappa approaches infinity, the t-distribution
#' converges to a standard normal distribution. When kappa is smaller, the tails
#' of the distribution are heavier, permitting the possibility that larger or
#' smaller values of delta to have greater prior probability.
#'
#' In the absence of prior information regarding the standardized effect size,
#' setting mu_delta = 0, gamma = 1, and kappa = 1 results in a generally fair
#' default prior. This arrangment produces a Cauchy distribution (i.e., a t-
#' distribution where df = 1) centered at zero and with variance of one.
#'
#' @param delta continuous real value for the standardized effect size
#' @param t t-statistic computed from the sample
#' @param n effective sample size of the t-test
#' @param nu degrees of freedom for the t-test
#' @param mu_delta hyperprior on the location of the standardized effect size
#' @param gamma hyperprior on the scale of the standardized effect size
#' @param kappa hyperprior on the degrees of freedom for the t-distribution prior
#'   on the standardized effect size
#'
#' @returns point-wise probability density of delta for the values entered in
#'   the remainder of the function
#'
#' @references
#' Gronau, Q. F., Ly, A., & Wagenmakers, E.-J. (2020). Informed Bayesian t-tests.
#' The American Statistician, 74, 137-143. https://doi.org/10.1080/00031305.2018.1562983
#'
integrand_t <- function(delta, t, n, nu, mu_delta, gamma, kappa) {

  suppressWarnings(
    dt(x = t, df = nu, ncp = sqrt(n) * delta) *
      1 / gamma * dt( (delta - mu_delta) / gamma, df = kappa)
  )
}


#' Calculation of the posterior probability distribution at specific points
#'
#' An integral function needed for computing the cumulative density of the
#' posterior distribution for t. Integration is conducted with respect to the
#' standardized effect size.
#'
#' @param delta continuous real value for the standardized effect size
#' @param t t-statistic computed from the sample
#' @param n1 integer value for the size of the first sample
#' @param n2 integer value for the size of the second sample
#' @param paired logical indicator for whether the samples are paired (true)
#'   or independent (false)
#' @param mu_delta hyperprior on the location of the standardized effect size
#' @param gamma hyperprior on the scale of the standardized effect size
#' @param kappa hyperprior on the degrees of freedom for the t-distribution prior
#'   on the standardized effect size
#' @param rel_tol performance-related real value set to the acceptable tolerance
#'   bound on the calculated integral function. High values will produce less
#'   accurate approximations of the true integral while very low values will
#'   produce more accurate approximations at the cost of increased computation
#'   time. The relative cost of compuation to increased precision drops off
#'   rapidly.
#'
#' @returns point-wise posterior probability of delta for the values entered in
#'   the remainder of the function
#'
posterior_t <- function(delta, t, n1, n2 = NULL, paired = TRUE, mu_delta, gamma,
                        kappa, rel_tol = .Machine$double.eps^0.25) {
  # calculate effective sample size and degrees of freedom
  n_eff <- ifelse(paired,
                  n1, n1 * n2 / (n1 + n2))
  nu <- ifelse(paired,
               n1 - 1, n1 + n2 - 2)

  # calulation of numerator
  num <- suppressWarnings(
    dt(x = t, df = nu, ncp = sqrt(n_eff) * delta) *
      1 / gamma * dt( (delta - mu_delta) / gamma, df = kappa)
  )

  # calculation of denominator
  den <- integrate(integrand_t, lower = -Inf, upper = Inf,
                   t = t, n = n_eff, nu = nu,
                   mu_delta = mu_delta, gamma = gamma, kappa = kappa,
                   rel.tol = rel_tol)$value

  # calculation of posterior probability
  out <- num / den
  out[is.na(out)] <- 0

  # return point-wise probability
  return( out )
}

#' Calculate the posterior cumulative probability distribution function
#'
#' Combines the posterior_t() and integrand_t() functions to empirically find
#' the cumulative probability of a particular effect size under the posterior
#' distribution obtained via the remaining arguments.
#'
#' @param x effect size of interest for the cumulative probability
#' @param t t-statistic computed from the sample
#' @param n1 integer value for the size of the first sample
#' @param n2 integer value for the size of the second sample
#' @param paired logical indicator for whether the samples are paired (true)
#'   or independent (false)
#' @param mu_delta hyperprior on the location of the standardized effect size
#' @param gamma hyperprior on the scale of the standardized effect size
#' @param kappa hyperprior on the degrees of freedom for the t-distribution prior
#'   on the standardized effect size
#' @param rel_tol performance-related real value set to the acceptable tolerance
#'   bound on the calculated integral function. High values will produce less
#'   accurate approximations of the true integral while very low values will
#'   produce more accurate approximations at the cost of increased computation
#'   time. The relative cost of compuation to increased precision drops off
#'   rapidly.
#'
#' @returns cumulative posterior probability of delta for the values entered in
#'   the remainder of the function
#'
cdf_t <- function(x, t, n1, n2 = NULL, paired = TRUE, mu_delta, gamma, kappa,
                  rel_tol = .Machine$double.eps^0.25) {

  # perform integration up to x
  out <- integrate(posterior_t, lower = -Inf, upper = x, t = t, n1 = n1, n2 = n2,
                   paired = paired, mu_delta = mu_delta, gamma = gamma,
                   kappa = kappa, rel.tol = rel_tol)$value

  # catch numerical errors
  if (out > 1 & out < 1.001) {
    out <- 1
  }

  # return cumulative probability
  return( out )
}

#' Quantile function of the posterior distribution of delta
#'
#' Quantiles from the posterior distribution are estimated via the Newton-
#' Raphson method. The procedure is iterative and relies on a convergence
#' tolerance (`tol`). This convergence tolerance is separate from the integral
#' tolerance (`rel_tol`). As the estimation of the quantile value continues
#' until the convergence tolerance has been met, reducing `tol` will increase
#' overall computation time but also produces greater overall numerical
#' precision, though the trade-off of time to accuracy rapidly diminshes.
#'
#' An alternative means of improving accuracy and controlling estimation time is
#' by setting a maximum number of iterations. If the convergence tolerance is
#' met before the maximum number of iterations set by the user, then the function
#' returns the result and terminates without needing to use whatever number of
#' remaining iterations. As the convergence tolerance is reduced, the likelihood
#' increases that the function will complete the maximum number of iterations
#' before convergence is reached.
#'
#' Note that no warnings or information messages are pushed to the user in this
#' function. As a result, a user will not definitively know whether this
#' function returns a result that reached the convergence tolerance or hit the
#' maximum number of iterations. The default settings for `tol` and `max_iter`
#' are adequate for most purposes, but if these defaults are changed, then it is
#' encouraged that users experiment with the settings beforehand to ensure that
#' convergence is being reliably reached.
#'
#' @param q probability value corresponding to the desired quantile
#' @param t t-statistic computed from the sample
#' @param n1 integer value for the size of the first sample
#' @param n2 integer value for the size of the second sample
#' @param paired logical indicator for whether the samples are paired (true)
#'   or independent (false)
#' @param mu_delta hyperprior on the location of the standardized effect size
#' @param gamma hyperprior on the scale of the standardized effect size
#' @param kappa hyperprior on the degrees of freedom for the t-distribution prior
#'   on the standardized effect size
#' @param tol real decimal value typically on the order of 1e-3 to 1e-6 that
#'   defines the acceptable difference between estimates obtained on the current
#'   versus previous iteration before the difference is considered negligible
#' @param max_iter maximum number of iterations to perform before stopping the
#'   estimation process
#' @param rel_tol performance-related real value set to the acceptable tolerance
#'   bound on the calculated integral function. High values will produce less
#'   accurate approximations of the true integral while very low values will
#'   produce more accurate approximations at the cost of increased computation
#'   time. The relative cost of compuation to increased precision drops off
#'   rapidly.
#'
#' @returns value of delta corresponding to the desired quantile for the values
#'   entered in the remainder of the function
#'
quantile_t <- function(q, t, n1, n2 = NULL, paired = TRUE, mu_delta, gamma,
                       kappa, tol = 0.0001, max_iter = 100,
                       rel_tol = .Machine$double.eps^0.25) {

  # initial search for reasonable starting value
  delta <- seq(-2.5, 2.5, length.out = 500)
  dens <- posterior_t(delta, t = t, n1 = n1, n2 = n2, paired = paired,
                      mu_delta = mu_delta, gamma = gamma, kappa = kappa)
  x_new <- delta[which.max(dens)]

  # initalize arguments for the while() statement
  x_cur <- Inf
  i <- 1

  # iterate through estimation until tol or max_iter are reached
  while (abs(x_cur - x_new) > tol && i < max_iter) {

    # update current best estimate
    x_cur <- x_new

    # increment precision for a new estimate
    x_new <- x_cur -
      (cdf_t(x_cur, t = t, n1 = n1, n2 = n2, paired = paired,
             mu_delta = mu_delta, gamma = gamma, kappa = kappa,
             rel_tol = rel_tol) - q) /
      posterior_t(x_cur, t = t, n1 = n1, n2 = n2, paired = paired,
                  mu_delta = mu_delta, gamma = gamma, kappa = kappa)

    # increment the iteration indicator
    i <- i + 1
  }

  # return resulting delta estimate
  return( x_new )
}

#' Obtain credible intervals with robust point estimate of effect size
#'
#' Utilizes the quantile_t() function to return point estimates of the lower and
#' upper bounds of the requested credible interval as well as the median value
#' of delta from the posterior distribution. The median point estimate is more
#' appropriate than the posterior mean given that the non-centered t-distributions
#' are skewed with the severity of skew dependent on the non-centrality
#' parameter.
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
#' @param int continuous real value between 0 and 1 corresponding to the desired
#'   credible interval (e.g., `int = 0.95` will return the 95% credible interval)
#' @param type character value reflecting the kind of test/interval wanted. Valid
#'   options include "two-sided" (interval over positive and negative values),
#'   "plus-sided" (interval distributed only over positive effects), and
#'   "min-sided" (interval distributed only over negative effects).
#' @param tol real decimal value typically on the order of 1e-3 to 1e-6 that
#'   defines the acceptable difference between estimates obtained on the current
#'   versus previous iteration before the difference is considered negligible
#' @param max_iter maximum number of iterations to perform before stopping the
#'   estimation process
#' @param rel_tol performance-related real value set to the acceptable tolerance
#'   bound on the calculated integral function. High values will produce less
#'   accurate approximations of the true integral while very low values will
#'   produce more accurate approximations at the cost of increased computation
#'   time. The relative cost of compuation to increased precision drops off
#'   rapidly.
#'
#' @returns value of delta corresponding to the desired quantile for the values
#'   entered in the remainder of the function
#'
ci_t <- function(t, n1, n2 = NULL, paired = TRUE, mu_delta, gamma, kappa,
                 int = 0.95, type = "two-sided", tol = 0.0001, max_iter = 100,
                 rel_tol = .Machine$double.eps^0.25) {

  # convert requested credible interval to associated p-values
  lower  <- (1 - int)/2
  upper  <- int + (1 - int)/2
  median <- .5

  # adjust interval bounds based on one-sided test requests
  if (type == "plus-sided") {

    post_less0 <- cdf_t(x = 0, t = t, n1 = n1, n2 = n2, paired = paired,
                        mu_delta = mu_delta, gamma = gamma, kappa = kappa,
                        rel_tol = rel_tol)

    lower  <- post_less0 + (1 - post_less0) * lower
    upper  <- post_less0 + (1 - post_less0) * upper
    median <- post_less0 + (1 - post_less0) * median

  } else if (type == "min-sided") {

    post_less0 <- cdf_t(x = 0, t = t, n1 = n1, n2 = n2, paired = paired,
                        mu_delta = mu_delta, gamma = gamma, kappa = kappa,
                        rel_tol = rel_tol)

    lower  <- post_less0 * lower
    upper  <- post_less0 * upper
    median <- post_less0 * median

  }

  # find standardized effect size associated with the bounds
  lb <- quantile_t(lower, t = t, n1 = n1, n2 = n2, paired = paired,
                   mu_delta = mu_delta, gamma = gamma, kappa = kappa,
                   rel_tol = rel_tol)
  ub <- quantile_t(upper, t = t, n1 = n1, n2 = n2, paired = paired,
                   mu_delta = mu_delta, gamma = gamma, kappa = kappa,
                   rel_tol = rel_tol)
  mdn <- quantile_t(median, t = t, n1 = n1, n2 = n2, paired = paired,
                    mu_delta = mu_delta, gamma = gamma, kappa = kappa,
                    rel_tol = rel_tol)

  # combine results as a vector
  out <- c(lb, mdn, ub)

  # return result
  return( out )
}

#' Calculate the Bayes factor for an observed effect size
#'
#' The calculation of the Bayes factor for an informed Bayesian t-test is given
#' by equation (5) of Gronau et al. (2020).
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
#' @param rel_tol performance-related real value set to the acceptable tolerance
#'   bound on the calculated integral function. High values will produce less
#'   accurate approximations of the true integral while very low values will
#'   produce more accurate approximations at the cost of increased computation
#'   time. The relative cost of compuation to increased precision drops off
#'   rapidly.
#'
#' @returns BF10 (evidence for the alternative) and BF01 (evidence for the null)
#'   for the standardized effect size observed in a sample given the priors
#'
#' @references
#' Gronau, Q. F., Ly, A., & Wagenmakers, E.-J. (2020). Informed Bayesian t-tests.
#' The American Statistician, 74, 137-143. https://doi.org/10.1080/00031305.2018.1562983
#'
bf_t_point <- function(t, n1, n2 = NULL, paired = TRUE, mu_delta, gamma, kappa,
                       rel_tol = .Machine$double.eps^0.25) {

  # calculate effective sample size and degrees of freedom
  n_eff <- ifelse(paired,
                  n1, n1 * n2 / (n1 + n2))
  nu <- ifelse(paired,
               n1 - 1, n1 + n2 - 2)

  # calculate the numerator of the Bayes factor
  num <- integrate(integrand_t, lower = -Inf, upper = Inf, t = t, n = n_eff,
                   nu = nu, mu_delta = mu_delta, gamma = gamma, kappa = kappa,
                   rel.tol = rel_tol)$value

  # calculate denominator of the Bayes factor
  den <- dt(x = t, df = nu)

  # get result as a vector
  out <- c(num / den, den / num)

  # return result
  return( out )
}

#' Calculate the Bayes factor for an equivalence interval
#'
#' Extension of bf_point() function specific to cases where an equivalence
#' interval is available.
#'
#' @param t t-statistic computed from the sample
#' @param n1 integer value for the size of the first sample
#' @param n2 integer value for the size of the second sample
#' @param paired logical indicator for whether the samples are paired (true)
#'   or independent (false)
#' @param interval numeric vector of length 2 providing the lower and upper
#'   bounds, respectively, of the equivalence interval for the standardized
#'   effect size
#' @param mu_delta hyperprior on the location of the standardized effect size
#' @param gamma hyperprior on the scale of the standardized effect size
#' @param kappa hyperprior on the degrees of freedom for the t-distribution prior
#'   on the standardized effect size
#' @param rel_tol performance-related real value set to the acceptable tolerance
#'   bound on the calculated integral function. High values will produce less
#'   accurate approximations of the true integral while very low values will
#'   produce more accurate approximations at the cost of increased computation
#'   time. The relative cost of compuation to increased precision drops off
#'   rapidly.
#'
#' @returns BF10 (evidence for the alternative) and BF01 (evidence for the null)
#'   for the standardized effect size observed in a sample given the priors
#'
bf_t_inter <- function(t, n1, n2 = NULL, paired = TRUE, interval, mu_delta, gamma,
                       kappa, rel_tol = .Machine$double.eps^0.25) {

  # calculate effective sample size and degrees of freedom
  n_eff <- ifelse(paired,
                  n1, n1 * n2 / (n1 + n2))
  nu <- ifelse(paired,
               n1 - 1, n1 + n2 - 2)

  # calculate the probability terms over the equivalence intervals
  post_dens <- cdf_t(x = interval[2], t = t, n1 = n1, n2 = n2, paired = paired,
                     mu_delta = mu_delta, gamma = gamma, kappa = kappa,
                     rel_tol = rel_tol) -
    cdf_t(x = interval[1], t = t, n1 = n1, n2 = n2, paired = paired,
          mu_delta = mu_delta, gamma = gamma, kappa = kappa, rel_tol = rel_tol)

  prior_dens <- (1 / gamma) * pt((interval[2] - mu_delta / gamma), kappa) -
    (1 / gamma) * pt((interval[1] - mu_delta) / gamma, kappa)

  # confirm numerical stability of estimates
  if ( post_dens < 0 ) {
    post_dens <- 0
  }

  # calculate numerator and denominator of the Bayes factor
  num <- post_dens / prior_dens
  den <- (1 - post_dens) / (1 - prior_dens)

  # get result as a vector
  out <- c(den / num, num / den)

  # return result
  return( out )
}
