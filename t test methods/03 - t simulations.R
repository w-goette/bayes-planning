#' Simulation of independent samples study
#'
#' NOTE: Use of this function requires that the packages `tidyverse` and `parallel`
#'       have been or can be installed. The function will try to install any
#'       missing packages, but it is not guaranteed that this will be successful
#'       as machines not owned by the user may require special privileges to
#'       download R packages.
#'
#' Generic wrapper function for a complete simulation study of an independent
#' samples design. Passing numeric vectors - c(...) - to the function arguments
#' d, n1, n2, mu_delta, gamma, kappa, power, and alpha will result in the implicit
#' production of fully crossed conditions over which a total of `n_iter` simulations
#' will be run.
#'
#' It is important to note that the conditions are fully crossed as this can
#' result in a very rapid increase in computation time. For example, if one were
#' to set `d = c(0.8, 0.2)`, `n1 = c(100, 200)`, `n2 = (100, 200)` and leave the
#' remaining adjustable arguments at their defaults, then this produces 8 unique
#' conditions that will each be run for `n_iter` simulations: 2 x 2 x 2 = 8. If
#' one additional value were requested for each of d, n1, and n2, then the number
#' of conditions is now 3 x 3 x 3 = 27. If one were interested in then also
#' examining the effect of different values of the hyperpriors `mu_delta`, `gamma`,
#' and `kappa`, then even with just two unique values for each of these arguments
#' the function will run 2(d) x 2(n1) x 2(n2) x 2(mu_delta) x 2(gamma) x 2(kappa)
#' = 2^6 = 64 conditions each `n_iter` times. For this reason, it is recommended
#' that this function only be used when there are specific study design or effect
#' conditions that want to be checked (e.g., for sensitivity analyses).
#'
#' In comparison to recommendations for Bayes factor design analyses, the default
#' `n_iter = 1000` is an order of magnitude too small (general recommendation
#' is `n_iter = 10000`). For the aforementioned exponential increase in the
#' number of conditions possible in this general wrapper function, the default
#' `n_iter` is set lower to minimize potential computation times. Once a set of
#' reasonable conditions for a study design are identified, it is recommended
#' that `n_iter` be increased to at least 10000 (though more is also acceptable)
#' before final recommendations regarding sample sizes be made.
#'
#' @param n_iter integer for simulation iterations to run. Final guidance should
#'   be based on at least n_iter = 10000.
#' @param bf minimum Bayes factor desired as evidence for equivalence (or effect)
#'   by the end of the study. Following interpretation guidelines by Raftery:
#'      1 < BF <= 3:   Weak evidence
#'      3 < BF <= 20:  Positive evidence
#'     20 < BF <= 150: Strong evidence
#'          BF  > 150: Very strong evidence
#' @param d numeric value or vector for the true standardized effect size in the
#'   population
#' @param n1 integer value or vector for the number of observations to draw from
#'   the population in sample 1
#' @param n2 integer value or vector for the number of observations to draw from
#'   the population in sample 2
#' @param dist name of the distribution to use for transforming the continuous
#'   latent z-scaled values in the population to observable raw scores. Options
#'   include "bbin" for beta-binomial, "nbin" for negative binomial, and "pois"
#'   for Poisson.
#' @param mn average raw score in the population
#' @param sd standard deviation of raw scores in the population
#' @param size integer of items/trials on the simulated tests. Only relevant when
#'   dist = "bbin"
#' @param mu_delta numeric value or vector giving the location/mean/expectations
#'   of the standardized effect sizes to be used with an informed Bayesian t-test.
#'   Recommended default when no prior information is available is 0.
#' @param gamma numeric value or vector giving the scale/dispersion of the
#'   standardized effect sizes to be used with an informed Bayesian t-test. Common
#'   default when no prior information is available is sqrt(2)/2.
#' @param kappa numeric value or vector giving the degrees of freedom/confidence
#'   in the prior belief for the standardized effect sizes to be used in the
#'   informed Bayesian t-test. Recommended default when no prior information is
#'   available is 1.
#' @param h0 logical indicator for whether Bayes factors should be returned in
#'   terms of evidence for the null (TRUE; BF01 = no difference between groups)
#'   or for the alternative (FALSE; BF10 = non-equivalence of groups).
#' @param power numeric value or vector used in finding the equivalence bounds
#'   for the group means. The equivalence bounds will be defined as the smallest
#'   standardized effect size detectable in the sample with this level of power.
#'   This specification is derived from frequentist equivalence testing,
#'   particularly TOST methods, in which boundaries for "practically equivalent"
#'   are defined in terms of detectable true effects.
#' @param alpha numeric value or vector also used in finding the equivalence
#'   bounds for the group means. As these bounds are based in frequentist
#'   power analysis, the Type I error rate (alpha) needs to be provided.
#' @param n_cores integer value specifying the number of CPU cores to include in
#'   running these simulations. Generally speaking, using more cores will result
#'   in faster results, though there are diminishing returns. Using all available
#'   cores is not recommended as this may cause Windows PCs to crash and will
#'   generally cause a significant increase in computer memory usages that may
#'   limit usability of the machine while the function is running. The number of
#'   cores available on a particular machine can be obtained using the detectCores()
#'   function from the `parallel` package -- parallel::detectCores(). If a large
#'   number of conditions will be examined, a somewhat safe option may be to
#'   specify n_cores = parallel::detectCores() - 1.
#'
#' @returns a dataframe with rows corresponding to each condition examined in
#'   the simulation and columns giving general summaries of the results of the
#'   condition arguments. Most importantly, the final two columns of the dataframe
#'   are the average number of iterations (total iterations given by n_iter) in
#'   which the calculated Bayes factor met the minimum Bayes factor threshold.
#'   These averages can be thought of as analogous to study power in general
#'   simulation designs: i.e., 0.8963 mean that 89.63% of the simulated samples
#'   yielded evidence for the hypothesis of interest that was at least as strong
#'   as the minimum threshold entered in the `bf` argument of the function.
#'
bf_t_indendent <- function(n_iter = 1000, bf = 3, d, n1, n2, dist = "bbin", mn,
                           sd, size, mu_delta = 0, gamma = 1, kappa = 1,
                           h0 = TRUE, power = 0.33, alpha = 0.05, n_cores = 1) {

  # ensure needed libraries are installed and loaded
  dependencies <- c("tidyverse", "parallel")
  sapply(dependencies, function(x) {
    if ( !require(x, character.only = TRUE) ) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })

  # allow for parallel processing
  options(mc.cores = n_cores)

  # extract the fully crossed conditions
  conditions <- expand_grid(n1 = n1, n2 = n2, d = d, power = power, alpha = alpha,
                            mu_delta = mu_delta, gamma = gamma, kappa = kappa) %>%
    as.data.frame()

  # extract number of unique populations implied
  n_pop <- length(unique(conditions$d))

  # add column to conditions indexing the unique populations
  conditions$id_pop <- 0
  for ( i in 1 : n_pop ) {
    conditions[which(conditions$d == unique(conditions$d)[i]), "id_pop"] <- i
  }

  # generate population data for each condition
  pop <- lapply(1:n_pop, function(i) {
    gen_pop_t(d = unique(conditions$d)[i], paired = FALSE)
  })

  # iterate simulation and include summary values at end
  out <- lapply(1:n_iter, function(x) {

    # calculate values of interest (possibly in parallel)
    parallel::mclapply(1:nrow(conditions), function(i) {

      # obtain samples from populations
      sample <- get_sample_t(pop[[ conditions$id_pop[i] ]],
                             n1 = conditions$n1[i],
                             n2 = conditions$n2[i],
                             paired = FALSE,
                             dist = dist, mn = mn, sd = sd, size = size)

      # calculate t-statistic
      freq <- t_test(x1 = sample[[1]], x2 = sample[[2]], paired = FALSE)

      # find equivalence region implied by user inputs
      bounds <- min_effect(n1 = conditions$n1[i],
                           n2 = conditions$n2[i],
                           paired = FALSE,
                           power = conditions$power[i],
                           alpha = conditions$alpha[i])

      # calculate BFs for point and interval estimates
      bf_null <- bf_t_point(t = freq[1],
                            n1 = conditions$n1[i],
                            n2 = conditions$n2[i],
                            paired = FALSE,
                            mu_delta = conditions$mu_delta[i],
                            gamma = conditions$gamma[i],
                            kappa = conditions$kappa[i])

      bf_equiv <- bf_t_inter(t = freq[1],
                             n1 = conditions$n1[i],
                             n2 = conditions$n2[i],
                             paired = FALSE,
                             interval = bounds,
                             mu_delta = conditions$mu_delta[i],
                             gamma = conditions$gamma[i],
                             kappa = conditions$kappa[i])

      # combine results in a named vector
      out <- c(
        "n1" = conditions$n1[i],
        "n2" = conditions$n2[i],
        "true_d" = conditions$d[i],
        "power" = conditions$power[i],
        "alpha" = conditions$alpha[i],
        "mu_delta" = conditions$mu_delta[i],
        "gamma" = conditions$gamma[i],
        "kappa" = conditions$kappa[i],
        "emp_d" = freq[3],
        "raw_diff" = mean(sample[[1]]) - mean(sample[[2]]),
        "bf_null" = ifelse(h0, bf_null[2], bf_null[1]),
        "bf_equiv" = ifelse(h0, bf_equiv[2], bf_equiv[1])
      )

      # return results
      return( out )
    }) %>%
      bind_rows()

  }) %>%
    bind_rows() %>%
    group_by(n1, n2, true_d, power, alpha, mu_delta, gamma, kappa) %>%
    summarize(
      mean_stdfx = mean(emp_d),
      mean_diff  = mean(raw_diff),
      prop_null  = mean(bf_null >= bf),
      prop_equiv = mean(bf_equiv >= bf),
      .groups = "keep"
    ) %>%
    as.data.frame()

  # return results
  return( out )
}

#' Simulation of test/retest samples study
#'
#' NOTE: Use of this function requires that the packages `tidyverse` and `parallel`
#'       have been or can be installed. The function will try to install any
#'       missing packages, but it is not guaranteed that this will be successful
#'       as machines not owned by the user may require special privileges to
#'       download R packages.
#'
#' Generic wrapper function for a complete simulation study of a paired/repeated
#' samples design. Passing numeric vectors - c(...) - to the function arguments
#' d, n1, n2, mu_delta, gamma, kappa, power, and alpha will result in the implicit
#' production of fully crossed conditions over which a total of `n_iter` simulations
#' will be run.
#'
#' It is important to note that the conditions are fully crossed as this can
#' result in a very rapid increase in computation time. For example, if one were
#' to set `d = c(0.8, 0.2)`, `n = c(100, 200)` and leave the remaining adjustable
#' arguments at their defaults, then this produces 4 unique conditions that will
#' each be run for `n_iter` simulations: 2 x 2 = 4. If one additional value were
#' requested for each of d and n, then the number of conditions is now 3 x 3 = 6.
#' If one were interested in then also examining the effect of different values
#' of the hyperpriors `mu_delta`, `gamma`, and `kappa`, then even with just two
#' unique values for each of these arguments the function will run 2(d) x 2(n) x
#' 2(mu_delta) x 2(gamma) x 2(kappa) = 2^5 = 32 conditions each `n_iter` times.
#' For this reason, it is recommended that this function only be used when there
#' are specific study design or effect conditions that want to be checked (e.g.,
#' for sensitivity analyses).
#'
#' In comparison to recommendations for Bayes factor design analyses, the default
#' `n_iter = 1000` is an order of magnitude too small (general recommendation
#' is `n_iter = 10000`). For the aforementioned exponential increase in the
#' number of conditions possible in this general wrapper function, the default
#' `n_iter` is set lower to minimize potential computation times. Once a set of
#' reasonable conditions for a study design are identified, it is recommended
#' that `n_iter` be increased to at least 10000 (though more is also acceptable)
#' before final recommendations regarding sample sizes be made.
#'
#' @param n_iter integer for simulation iterations to run. Final guidance should
#'   be based on at least n_iter = 10000.
#' @param bf minimum Bayes factor desired as evidence for equivalence (or effect)
#'   by the end of the study. Following interpretation guidelines by Raftery:
#'      1 < BF <= 3:   Weak evidence
#'      3 < BF <= 20:  Positive evidence
#'     20 < BF <= 150: Strong evidence
#'          BF  > 150: Very strong evidence
#' @param d numeric value or vector for the true standardized effect size in the
#'   population
#' @param prac_effect numeric value or vector reflecting the expected increase
#'   in raw score attributable to previous exposure/practice with the test. When
#'   this argument is non-zero, the effect is converted to a standardized effect
#'   size and added to d. Thus, if one wants to simulate expected score differences
#'   when the only effect on scores is a practice effect, then d should be zero.
#'   Setting both prac_effect and d to non-zero values/vectors will result in
#'   simulating population effects in which part of mean differences is due to
#'   the effects of practice and another part is due to other group differences.
#' @param n integer value or vector for the number of observations to draw from
#'   the population
#' @param dist name of the distribution to use for transforming the continuous
#'   latent z-scaled values in the population to observable raw scores. Options
#'   include "bbin" for beta-binomial, "nbin" for negative binomial, and "pois"
#'   for Poisson.
#' @param mn average raw score in the population
#' @param sd standard deviation of raw scores in the population
#' @param size integer of items/trials on the simulated tests. Only relevant when
#'   dist = "bbin"
#' @param r12 numeric value or vector test-retest correlation of the measure
#' @param rope_ci numeric value or vector between 0 and 1 reflecting the span of
#'   the confidence interval desired for equivalence testing via a region of
#'   practical equivalence (ROPE)
#' @param nil logical indicator for whether the effect estimates should be for
#'   a nil hypothesis (i.e., null = 0) or relative to the practice effect. If
#'   nil = TRUE, then the expected practice effect is removed from estimates and
#'   the equivalence testing is conducted against the remaining effect being zero.
#'   If nil = FALSE, then equivalence testing is conducted against the whether
#'   the observed effect is different from just the practice effects. Regardless
#'   of usage, the hyperprior mu_delta should remain scaled to a nil effect as
#'   it will be appropriately adjusted internally.
#' @param mu_delta numeric value or vector giving the location/mean/expectations
#'   of the standardized effect sizes to be used with an informed Bayesian t-test.
#'   Recommended default when no prior information is available is 0.
#' @param gamma numeric value or vector giving the scale/dispersion of the
#'   standardized effect sizes to be used with an informed Bayesian t-test. Common
#'   default when no prior information is available is sqrt(2)/2.
#' @param kappa numeric value or vector giving the degrees of freedom/confidence
#'   in the prior belief for the standardized effect sizes to be used in the
#'   informed Bayesian t-test. Recommended default when no prior information is
#'   available is 1.
#' @param h0 logical indicator for whether Bayes factors should be returned in
#'   terms of evidence for the null (TRUE; BF01 = no difference between groups)
#'   or for the alternative (FALSE; BF10 = non-equivalence of groups).
#' @param n_cores integer value specifying the number of CPU cores to include in
#'   running these simulations. Generally speaking, using more cores will result
#'   in faster results, though there are diminishing returns. Using all available
#'   cores is not recommended as this may cause Windows PCs to crash and will
#'   generally cause a significant increase in computer memory usages that may
#'   limit usability of the machine while the function is running. The number of
#'   cores available on a particular machine can be obtained using the detectCores()
#'   function from the `parallel` package -- parallel::detectCores(). If a large
#'   number of conditions will be examined, a somewhat safe option may be to
#'   specify n_cores = parallel::detectCores() - 1.
#'
#' @returns a dataframe with rows corresponding to each condition examined in
#'   the simulation and columns giving general summaries of the results of the
#'   condition arguments. Most importantly, the final two columns of the dataframe
#'   are the average number of iterations (total iterations given by n_iter) in
#'   which the calculated Bayes factor met the minimum Bayes factor threshold.
#'   These averages can be thought of as analogous to study power in general
#'   simulation designs: i.e., 0.8963 mean that 89.63% of the simulated samples
#'   yielded evidence for the hypothesis of interest that was at least as strong
#'   as the minimum threshold entered in the `bf` argument of the function.
#'
bf_t_retest <- function(n_iter = 1000, bf = 3, prac_effect, d, n, dist = "bbin",
                        mn, sd, size, r12, rope_ci = 0.90, nil = FALSE,
                        mu_delta = 0, gamma = 1, kappa = 1, h0 = TRUE, n_cores = 1) {

  # ensure needed libraries are installed and loaded
  dependencies <- c("tidyverse", "parallel")
  sapply(dependencies, function(x) {
    if ( !require(x, character.only = TRUE) ) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })

  # allow for parallel processing
  options(mc.cores = n_cores)

  # extract the fully crossed conditions
  conditions <- expand_grid(n = n, d = d, prac_effect = prac_effect, r12 = r12,
                            rope_ci = rope_ci, mu_delta = mu_delta,
                            gamma = gamma, kappa = kappa) %>%
    as.data.frame() %>%
    mutate(
      prac_scaled = scale_prac(prac_effect, sd, r12),
      total_d = d + prac_scaled
    )

  # add column to conditions indexing the unique populations
  conditions$id_pop <- as.numeric(as.factor(paste0(conditions$total_d, conditions$r12, sep = "_")))

  # get number of unique populations implied by settings
  n_pop <- length(unique(conditions$id_pop))

  # generate population data for each condition
  pop <- lapply(1 : n_pop, function(i) {
    gen_pop_t(d = conditions[which(conditions$id_pop == i)[1], "total_d"],
              paired = TRUE,
              r12 = conditions[which(conditions$id_pop == i)[1], "r12"])
  })

  # iterate simulation and include summary values at end
  out <- lapply(1:n_iter, function(x) {

    # calculate values of interest (possibly in parallel)
    parallel::mclapply(1:nrow(conditions), function(i) {

      # obtain samples from populations
      sample <- get_sample_t(pop[[ conditions$id_pop[i] ]],
                             n1 = conditions$n[i],
                             n2 = conditions$n[i],
                             paired = TRUE,
                             dist = dist, mn = mn, sd = sd, size = size)

      # respond to nil request
      if ( nil ) {
        # calculate t-statistic
        freq <- t_test(x1 = sample[[1]], x2 = sample[[2]], paired = TRUE,
                       mu = conditions$prac_effect[i])

        # find equivalence region implied by user inputs
        bounds <- prac_rope(prac_effect = conditions$prac_effect[i],
                            x1 = sample[[1]], x2 = sample[[2]],
                            ci = conditions$rope_ci[i],
                            null_centered = TRUE)

        # shift mu_delta location
        mu_delta <- conditions$mu_delta[i]

      } else {
        # calculate t-statistic
        freq <- t_test(x1 = sample[[1]], x2 = sample[[2]], paired = TRUE)

        # find equivalence region implied by user inputs
        bounds <- prac_rope(prac_effect = conditions$prac_effect[i],
                            x1 = sample[[1]], x2 = sample[[2]],
                            ci = conditions$rope_ci[i],
                            null_centered = FALSE)

        # shift mu_delta location
        mu_delta <- conditions$mu_delta[i] + scale_prac(conditions$prac_effect[i],
                                                        sd(c(sample[[1]], sample[[2]])),
                                                        cor(sample[[1]], sample[[2]]))
      }

      # calculate BFs for point and interval estimates
      bf_null <- bf_t_point(t = freq[1],
                            n1 = conditions$n[i],
                            paired = TRUE,
                            mu_delta = mu_delta,
                            gamma = conditions$gamma[i],
                            kappa = conditions$kappa[i])

      bf_equiv <- bf_t_inter(t = freq[1],
                             n1 = conditions$n[i],
                             paired = TRUE,
                             interval = bounds,
                             mu_delta = mu_delta,
                             gamma = conditions$gamma[i],
                             kappa = conditions$kappa[i])

      # combine results in a named vector
      out <- c(
        "n" = conditions$n[i],
        "ind_d" = conditions$d[i],
        "prac_d" = conditions$prac_scaled[i],
        "total_d" = conditions$total_d[i],
        "r12" = conditions$r12[i],
        "rope_ci" = conditions$rope_ci[i],
        "mu_delta" = conditions$mu_delta[i],
        "gamma" = conditions$gamma[i],
        "kappa" = conditions$kappa[i],
        "emp_d" = freq[3],
        "raw_diff" = mean(sample[[2]] - sample[[1]]),
        "bf_null" = ifelse(h0, bf_null[2], bf_null[1]),
        "bf_equiv" = ifelse(h0, bf_equiv[2], bf_equiv[1])
      )

      # return results
      return( out )
    }) %>%
      bind_rows()

  }) %>%
    bind_rows() %>%
    group_by(n, ind_d, prac_d, total_d, r12, rope_ci, mu_delta, gamma, kappa) %>%
    summarize(
      mean_stdfx = mean(emp_d),
      mean_diff  = mean(raw_diff),
      prop_null  = mean(bf_null >= bf),
      prop_equiv = mean(bf_equiv >= bf),
      .groups = "keep"
    ) %>%
    as.data.frame()

  # return results
  return( out )
}

#' Generate population data
#'
#' Helper function for simulating population data given user inputs.
#'
#' @param n number of simulated observations to produce that will define the
#'   true size of the population
#' @param d numeric value representing the standardized effect size to apply to
#'   the two groups
#' @param paired logical indicator for whether the population reflects an
#'   independent or paired design
#' @param r12 the test-retest correlation for paired designs
#'
#' @returns an n x 2 matrix of z-scaled values representing the population from
#'   which study samples could be drawn
#'
gen_pop_t <- function(n = 1e8, d, paired = TRUE, r12 = 0.80) {

  if ( paired ) {
    L <- chol(matrix(c(1, r12, r12, 1), 2, 2))
    z <- matrix(rnorm(2 * n), n, 2)
    x <- z %*% L
    x[, 2] <- x[, 2] + d
    return( x )

  } else {
    x1 <- rnorm(n)
    x2 <- rnorm(n) - d
    return( cbind(x1, x2) )
  }
}

#' Obtain random sample from simulated population
#'
#' Given a population, extract a random sample of observations and scale the
#' latent continuous variable to raw scores given the inputs.
#'
#' @param pop an N x 2 matrix of latent z-scaled values where each of N rows
#'   corresponds to one member of the population and the 2 columns define values
#'   for two populations (either independent or dependent)
#' @param n1 integer value reflecting the number of observations to draw from the
#'   first column of the population matrix to define sample 1
#' @param n2 integer value reflecting the number of observations to draw from the
#'   second column of the population matrix to define sample 2
#' @param paired logical indicator for whether the population/sample is a paired
#'   (TRUE) or independent (FALSE) design
#' @param dist name of the distribution to use for transforming the continuous
#'   latent z-scaled values in the population to observable raw scores. Options
#'   include "bbin" for beta-binomial, "nbin" for negative binomial, and "pois"
#'   for Poisson.
#' @param mn average raw score in the population
#' @param sd standard deviation of raw scores in the population
#' @param size integer of items/trials on the simulated tests. Only relevant when
#'   dist = "bbin"
#'
#' @returns a list with element 1 containing the vector of observations in
#'   sample 1 and element 2 the vector of observations in sample 2
#'
get_sample_t <- function(pop, n1, n2, paired = TRUE, dist, mn, sd, size) {

  # extract random sample from the population
  if ( paired ) {
    draw <- sample(nrow(pop), n1)
    x1 <- pop[draw, 1]
    x2 <- pop[draw, 2]

  } else {
    x1 <- pop[sample(nrow(pop), n1), 1]
    x2 <- pop[sample(nrow(pop), n2), 2]
  }

  # apply re-scaling based on requested distribution type
  if ( dist == "bbin" ) {
    x1 <- norm2bbin(x1, mn, sd, size)
    x2 <- norm2bbin(x2, mn, sd, size)

  } else if ( dist == "nbin" ) {
    x1 <- norm2nbin(x1, mn, sd)
    x2 <- norm2nbin(x2, mn, sd)

  } else if ( dist == "pois" ) {
    x1 <- norm2pois(x1, mn)
    x2 <- norm2pois(x2, mn)
  }

  # return result
  return( list(x1, x2) )
}
