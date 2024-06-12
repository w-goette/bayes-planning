#' Conversion to a beta-binomial variable
#'
#' norm2bbin() takes a z-scaled numeric vector or real and converts the value(s)
#' so as to follow a beta-binomially distributed variable with a particular
#' mean and standard deviation that is specified by the user.
#'
#' The function borrows from mixed-continuous copula methods in that the input
#' values are converted to the probabilities according to the standard normal.
#' These probabilities are then passed to the quantile function of the beta-
#' binomial in order to obtain the discrete value of this distribution that is
#' equivalent, with respect to the probability function, as the inputted z-score.
#'
#' The function can be thought of as projecting a continuous, latent variable
#' with a standard normal marginal distribution onto a simple probability-to-
#' quantile mapping function to obtain a beta-binomial marginal distribution.
#' The relative ordering of the scores is maintained, but absolute differences
#' between scores are not since the beta-binomial distribution can be seen as
#' a kind of binning of scores, limiting variability by restricting differences
#' to only integer values.
#'
#' @param q z-scaled quantile to convert to beta-binomial
#' @param mean desired mean of the resulting beta-binomial
#' @param sd desired standard deviation of the resulting beta-binomial
#' @param size number of items/highest raw score possible
#'
#' @returns discrete, lower- and upper-bounded score equivalents of the input
#'   z-scores
#'
norm2bbin <- function(q, mean, sd, size) {
  # convert score to probability relative to input mean and standard deviation
  p <- pnorm(q)

  # obtain beta-binomial shape parameters
  shape <- bb_shape(mean = mean, sd = sd, size = size)

  # convert probability to the beta-binomial scale
  qbbin(p = p, size = size, alpha = shape[1], beta = shape[2])
}

#' Conversion to a negative binomial variable
#'
#' norm2nbin() takes a z-scaled numeric vector or real and converts the value(s)
#' so as to follow a negative binomially distributed variable with a particular
#' mean and standard deviation that is specified by the user.
#'
#' The function borrows from mixed-continuous copula methods in that the input
#' values are converted to the probabilities according to the standard normal.
#' These probabilities are then passed to the quantile function of the negative
#' binomial in order to obtain the discrete value of this distribution that is
#' equivalent, with respect to the probability function, as the inputted z-score.
#'
#' The function can be thought of as projecting a continuous, latent variable
#' with a standard normal marginal distribution onto a simple probability-to-
#' quantile mapping function to obtain a negative binomial marginal distribution.
#' The relative ordering of the scores is maintained, but absolute differences
#' between scores are not since the negative binomial distribution can be seen
#' as a kind of binning of scores, limiting variability by restricting
#' differences to only integer values.
#'
#' WARNING: The negative binomial distribution is only defined when the variance
#'          is larger than the mean. If the variance is similar to the mean, then
#'          the Poisson distribution is likely the more appropriate marginal.
#'
#' @param q z-scaled quantile to convert to negative binomial
#' @param mean desired mean of the resulting negative binomial
#' @param sd standard deviation of the resulting negative binomial
#'
#' @returns discrete, lower-bounded score equivalents of the input z-scores
#'
norm2nbin <- function(q, mean, sd) {
  # convert score to probability relative to input mean and standard deviation
  p <- qnorm(q)

  # obtain negative binomial size parameter
  #Note: defined only when the variance > mean
  size <- mean^2 / (sd^2 - mean)

  #convert probability to the negative binomial scale
  qnbinom(p, size = size, mu = mean)
}

#' Conversion to a Poisson variable
#'
#' norm2pois() takes a z-scaled numeric vector or real and converts the value(s)
#' so as to follow a Poisson distributed variable with a particular mean that is
#' specified by the user. By definition, the mean and variance of a Poisson
#' distribution are the same.
#'
#' The function borrows from mixed-continuous copula methods in that the input
#' values are converted to the probabilities according to the standard normal.
#' These probabilities are then passed to the quantile function of the
#' Poisson in order to obtain the discrete value of this distribution that is
#' equivalent, with respect to the probability function, as the inputted z-score.
#'
#' The function can be thought of as projecting a continuous, latent variable
#' with a standard normal marginal distribution onto a simple probability-to-
#' quantile mapping function to obtain a Poisson marginal distribution. The
#' relative ordering of the scores is maintained, but absolute differences
#' between scores are not since the Poisson distribution can be seen as a kind
#' of binning of scores, limiting variability by restricting differences to only
#' integer values.
#'
#' @param q z-scaled quantile to convert to Poisson
#' @param mean desired mean of the resulting Poisson
#'
#' @returns discrete, lower-bounded score equivalents of the input z-scores
#'
norm2pois <- function(q, mean) {
  # convert score to probability relative to input mean and standard deviation
  p <- qnorm(q)

  #convert probability to the Poisson scale
  qpois(p, lambda = mean)
}



#' Likelihood function of a beta-binomial distribution
#'
#' dbbin() returns the likelihood (or optionally, log-likelihood) of a beta-
#' binomial distribution with specific size and shape parameters. Base R stats
#' does not include the beta-binomial family, requiring its implementation here.
#'
#' @param x integer vector or real for which the likelihood is needed
#' @param size number of trials/items in the binomial component
#' @param alpha first shape parameter for the beta component
#' @param beta second shape parameter for the beta component
#' @param log logical for whether log-likelihood should be returned (false by
#'   default)
#'
#' @returns (log-)probability of x value(s) according to the defined beta-
#'   binomial
#'
dbbin <- function(x, size, alpha, beta, log = FALSE) {
  lpmf <- rep(-Inf, length(x))

  lpmf <- lchoose(size, x) +
    lbeta(x + alpha, size - x + beta) -
    lbeta(alpha, beta)
  if ( log ) {
    return(lpmf)
  } else {
    return( exp(lpmf) )
  }
}

#' Probability mass function of a beta-binomial distribution
#'
#' pbbin() function encodes the probability mass function (pmf) of the beta-
#' binomial distribution. When `lower.tail = TRUE`, the result is the cumulative
#' mass function (cmf). When `lower.tail = FALSE`, the result is the
#' complementary cumulative mass function (ccmf).
#'
#' @param q integer vector or real for which the cmf or ccmf is needed
#' @param size number of trials/items in the binomial component
#' @param alpha first shape parameter for the beta component
#' @param beta second shape parameter for the beta component
#' @param lower.tail logical for whether to include the lower (true) or upper
#'   (false) tail of the distribution
#' @param log.p logical for whether the resulting probability should be log-
#'   scaled or not (false by default)
#'
#' @returns (log-)probability of the q value(s) according to the defined beta-
#'   binomial
#'
pbbin <- function(q, size, alpha, beta, lower.tail = TRUE, log.p = FALSE) {
  q <- trunc(q)
  p <- sapply(q, function(x) {
    sum(
      vapply(
        0:x, dbbin, FUN.VALUE = numeric(1L), size = size, alpha = alpha, beta = beta
      )
    )
  })

  if ( !lower.tail ) {
    p <- 1 - p
  }

  if ( log.p ) {
    return( log(p) )
  } else {
    return(p)
  }
}

#' Quantile function of a beta-binomial distribution
#'
#' qbbin() is an implementation of the quantile function for the general beta-
#' binomial distribution. This implementation is by simple brute force, making
#' it's general efficiency drop as the `size` argument increases as this is
#' associated with a greater number of individual steps to check.
#'
#' @param p vector or real probability(-ies) of the desired quantile
#' @param size number of trials/items in the binomial component
#' @param alpha first shape parameter for the beta component
#' @param beta second shape parameter for the beta component
#' @param lower.tail logical for whether to include the lower (true) or upper
#'   (false) tail of the distribution
#' @param log.p logical for whether the resulting probability should be log-
#'   scaled or not (false by default)
#'
#' @returns value in the specified beta-binomial distribution associated with
#'   the requested quantiles p
#'
qbbin <- function(p, size, alpha, beta, lower.tail = TRUE, log.p = FALSE) {
  if ( log.p ) {
    p <- exp(p)
  }

  q <- sapply(p, function(x) {
    if ( 0 > x | x > 1 ) {
      return( NaN )
    }

    y <- 0

    while ( x > pbbin(y, size = size, alpha = alpha, beta = beta) ) {
      y <- y + 1
    }

    y
  })

  if ( lower.tail ) {
    return(q)
  } else {
    return(size - q)
  }
}

#' Helper function for obtaining beta shape parameters
#'
#' bb_shape() calculates the needed alpha and beta shape parameters that, when
#' used as part of a beta-binomial distribution with given size, will produce
#' a distribution with the requested mean and standard deviation.
#'
#' WARNING: The calculation is just a simple moment matching computation, which
#'          means that there are some use cases where the alpha and beta
#'          parameters could be negative. As the shape parameters of a beta
#'          distribution cannot be negative, this function simply takes the
#'          absolute value of whatever estimates are obtained. While this
#'          prevents any errors from arising, it can result in misleading and
#'          inaccurate results when sensible means, standard deviations, and
#'          sizes are not requested.
#'
#' @param mean desired average value of the resulting distribution
#' @param sd desired standard deviation of the resulting distribution
#' @param size number of trials/items in the resulting distribution
#'
#' @returns vector of length 2 providing the values for alpha and beta,
#'   respectively, needed to produce a beta-binomial distribution with the
#'   requested first and second moments and number of trials/size.
#'
bb_shape <- function(mean, sd, size) {
  var <- sd^2
  t <- (size * var) / mean - size + mean
  a <- (size * mean - var - mean^2) / t
  b <- (size - mean) * (size - (var + mean^2) / mean) / t
  return( abs(c(a, b)) )
}
