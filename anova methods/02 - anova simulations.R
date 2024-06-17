#' Simulation of mixed effects samples study
#'
#' NOTE: Use of this function requires that the packages `tidyverse`, `parallel`,
#'       `BayesFactor`, and `bayestestR` have been or can be installed. The
#'       function will try to install any missing packages, but it is not
#'       guaranteed that this will be successful as machines not owned by the
#'       user may require special privileges to download R packages.
#'
#' Generic wrapper function for a complete simulation study of a mixed effects
#' repeated measures ANOVA. The simulation follows a design in which individuals
#' are tested two times, under two modalities (face-to-face or videoconference),
#' and are randomly assigned to the order in which they will complete testing
#' under these two modalities. The simulation arguments permit specification of
#' a main effect for the modality, a main effect for testing session (i.e., a
#' practice effect), and an interaction effect between modality and testing
#' session (i.e., an effect for the order in which individuals completed their
#' two testing sessions). Passing numeric vectors - c(...) - to the function
#' arguments n, prior, mode_effect, prac_effect, int_effect, ran_sd, and r12
#' will result in the implicit production of fully crossed conditions over which
#' a total of `n_iter` simulations will be run.
#'
#' It is important to note that the conditions are fully crossed as this can
#' result in a very rapid increase in computation time. For this reason, it is
#' recommended that this function only be used when there are specific study or
#' effect conditions that want to be checked (e.g., for sensitivity analyses).
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
#' @param n integer value or vector for the number of observations to draw from
#'   the population. Since the ANOVA design is a 2x2, the total sample size will
#'   be n * 4 (two levels of the between-subject factor "order" x two levels of
#'   the within-subject factor "session")
#' @param mode_effect numeric value or vector of the standardized mean difference
#'   between the two levels of "mode"/"modality" (i.e., face-to-face vs video-
#'   conference)
#' @param prac_effect numeric value or vector reflecting the expected increase
#'   in raw score attributable to previous exposure/practice with the test. When
#'   this argument is non-zero, the effect is converted to a standardized effect
#'   size internally
#' @param int_effect numeric value or vector reflecting the standardized mean
#'   difference incurred from the interaction of modality and session number
#' @param ran_sd positive numeric value or vector for the standard deviation of
#'   the random intercepts in the mixed model. If this is not of interest, then
#'   it should be set to 1. This argument is potentially useful in cases where
#'   there is a question of whether individual (or within-subject) variability
#'   dominates overall performance variability. This scenario may arise, for
#'   example, in cases where participants are highly heterogeneous and have test
#'   performance differences for reasons other than the study design's main effects
#' @param prior character value or vector for the prior scaling on the fixed
#'   effect sizes. Options include "medium", "wide", and "ultrawide." These
#'   effectively set the uncertainty of the priors on the effect sizes and thus
#'   can be used for sensitivity analyses. Holding all other conditions of the
#'   simulation constant, using `prior = c("medium", "wide", "ultrawide")` will
#'   allow one to see how much of an effect the prior's informativeness has on
#'   the resulting proportion of results where the BF threshold is met
#' @param dist name of the distribution to use for transforming the continuous
#'   latent z-scaled values in the population to observable raw scores. Options
#'   include "bbin" for beta-binomial, "nbin" for negative binomial, and "pois"
#'   for Poisson
#' @param mn average raw score in the population
#' @param sd standard deviation of raw scores in the population
#' @param size integer of items/trials on the simulated tests. Only relevant when
#'   dist = "bbin"
#' @param r12 numeric value or vector test-retest correlation of the measure
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
#'   condition arguments. The Bayes factors summarized in these results are
#'   the model-averaged Bayes factors for the inclusion of the main effects of
#'   modality, session number/practice, and the interaction/order of testing.
#'   In effect, these Bayes factors reflect the amount of evidence from the data
#'   and the model for keeping a particular effect in the model. In cases where
#'   the effect of a particular variable is set to zero in the simulation, then
#'   the Bayes factors for the inclusion of that variable's effect should be <1.
#'   The proportion of all results where the Bayes factor for the variable's
#'   inclusion was at least as large as the input minimum threshold can be
#'   safely interpreted as analogous to study power in cases where the variable
#'   was simulated to have an effect. When the variable was simulated to have no
#'   effect, then the proportion of Bayes factors that met the requested minimum
#'   threshold of evidence can be interpreted as analogous to the Type I error
#'   rate (i.e., situations in which evidence was found for an effect that does
#'   not really exist).
#'
bf_mixed_aov <- function(n_iter = 1000, bf = 3, n, mode_effect, prac_effect,
                         int_effect, ran_sd, prior = "medium", dist = "bbin",
                         mn, sd, size, r12, n_cores = 1) {

  # ensure needed libraries are installed and loaded
  dependencies <- c("tidyverse", "parallel", "BayesFactor", "bayestestR")
  sapply(dependencies, function(x) {
    if ( !require(x, character.only = TRUE) ) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })

  # allow for parallel processing
  options(mc.cores = n_cores)

  # extract the fully crossed conditions
  conditions <- expand_grid(n = n, prior = prior, mode_effect = mode_effect,
                            prac_effect = prac_effect, int_effect = int_effect,
                            ran_sd = ran_sd, r12 = r12) %>%
    as.data.frame() %>%
    mutate(
      prac_effect = ( prac_effect / sqrt(2 * sd^2 - 2 * r12 * sd^2) ) * sqrt(2 * (1 - r12))
    )

  # add column to conditions indexing the unique populations
  conditions$id_pop <- as.numeric(as.factor(paste0(conditions$prior,
                                                   conditions$mode_effect,
                                                   conditions$prac_effect,
                                                   conditions$int_effect,
                                                   conditions$ran_sd,
                                                   conditions$r12, sep = "_")))

  # get number of unique populations implied by settings
  n_pop <- length(unique(conditions$id_pop))

  # generate population data for each condition
  pop <- lapply(1 : n_pop, function(i) {
    gen_pop_aov(fix_mode = conditions[which(conditions$id_pop == i)[1], "mode_effect"],
                fix_prac = conditions[which(conditions$id_pop == i)[1], "prac_effect"],
                fix_inter = conditions[which(conditions$id_pop == i)[1], "int_effect"],
                ran_id = conditions[which(conditions$id_pop == i)[1], "ran_sd"],
                r12 = conditions[which(conditions$id_pop == i)[1], "r12"])
  })

  # iterate simulation and include summary values at end
  out <- lapply(1:n_iter, function(x) {

    # calculate values of interest (possibly in parallel)
    parallel::mclapply(1:nrow(conditions), function(i) {

      # obtain samples from populations
      sample <- get_sample_aov(pop[[ conditions$id_pop[i] ]],
                               n = conditions$n[i], dist = dist, mn = mn,
                               sd = sd, size = size)

      # calculate BFs for multiple models
      bf <- generalTestBF(y ~ mode * prac + id, data = sample %>%
                            mutate(id = as.factor(id)), progress = FALSE,
                          whichRandom = "id", neverExclude = "id",
                          whichModels = "all", rscaleFixed = conditions$prior[i])

      # calculate marginal BFs for model terms
      bf <- bayesfactor_inclusion(bf)

      # convert BFs to named vector
      bf_inc <- exp(bf$log_BF)
      names(bf_inc) <- rownames(bf)

      # combine results in a named vector
      out <- c(
        "n" = conditions$n[i],
        "prior" = conditions$prior[i],
        "ff_vc_effect" = conditions$mode_effect[i],
        "t1_t2_effect" = conditions$prac_effect[i],
        "order_effect" = conditions$int_effect[i],
        "person_variability" = conditions$ran_sd[i],
        "test_retest" = conditions$r12[i],
        "bf_modality" = as.numeric(bf_inc["mode"]),
        "bf_session" = as.numeric(bf_inc["prac"]),
        "bf_order" = as.numeric(bf_inc["mode:prac"])
      )

      # return results
      return( out )
    }) %>%
      bind_rows()

  }) %>%
    bind_rows() %>%
    mutate(
      bf_modality = as.numeric(bf_modality),
      bf_session = as.numeric(bf_session),
      bf_order = as.numeric(bf_order)
    ) %>%
    group_by(n, prior, ff_vc_effect, t1_t2_effect, order_effect,
             person_variability, test_retest) %>%
    summarize(
      mean_modality = mean(bf_modality),
      sd_modality = sd(bf_modality),
      prop_modality = mean(bf_modality >= bf),

      mean_session = mean(bf_session),
      sd_session = sd(bf_session),
      prop_session = mean(bf_session >= bf),

      mean_order = mean(bf_order),
      sd_order = sd(bf_order),
      prop_order = mean(bf_session >= bf),
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
#' @param fix_mode numeric value representing the standardized mean difference
#'   between individuals at different levels of the modality factor
#' @param fix_prac numeric value representing the standardized mean difference
#'   between individuals at different levels of the testing session factor
#' @param fix_inter numeric value representing the standardized mean difference
#'   incurred by the interaction of modality and session (i.e., an order effect)
#' @param ran_id positive numeric value reflecting the standard deviation of the
#'   random intercepts given to each participant
#' @param r12 the test-retest correlation for paired designs
#'
#' @returns a long-format dataframe indexing population ids, levels of modality
#'   and testing session, and the outcome variable as a standardized z-scaled
#'   latent value
#'
gen_pop_aov <- function(n = 1e6, fix_mode, fix_prac, fix_inter, ran_id, r12 = 0.80) {

  # create design matrix
  x <- data.frame(
    mode_t1 = rep(c(-0.5, 0.5), each = n / 2),
    mode_t2 = rep(c(0.5, -0.5), each = n / 2),
    prac_t1 = rep(-0.5),
    prac_t2 = rep(0.5)
  )

  # calculate bivariate points
  S <- tcrossprod(c(ran_id, ran_id)) * matrix(c(1, r12, r12, 1), 2, 2)
  L <- chol(S)
  z <- matrix(rnorm(2 * n), n, 2)
  y <- z %*% L

  # compute linear predictions
  x$y_t1 <- y[, 1] + fix_mode * x$mode_t1 + fix_prac * x$prac_t1 + fix_inter * x$mode_t1 * x$prac_t1
  x$y_t2 <- y[, 2] + fix_mode * x$mode_t2 + fix_prac * x$prac_t2 + fix_inter * x$mode_t2 * x$prac_t2

  # add noise
  x$y_t1 <- x$y_t1
  x$y_t2 <- x$y_t2

  # convert data to long-format
  out <- data.frame(
    id = factor(rep(1:n, 2)),
    mode = c(x$mode_t1, x$mode_t2),
    prac = c(x$prac_t1, x$prac_t2),
    order = factor(ifelse(x$mode_t1 == -0.5, "A", "B")),
    y = c(x$y_t1, x$y_t2)
  )

  # return result
  return( out )
}

#' Obtain random sample from simulated population
#'
#' Given a population, extract a random sample of observations and scale the
#' latent continuous variable to raw scores given the inputs.
#'
#' @param pop a long-format dataframe defining population member ids, group
#'   assignments for modality and session, a helper interaction indicator, and
#'   standardized z-scaled latent ability estimates
#' @param n integer value reflecting the number of observations to draw from the
#'   each cell of the implied 2x2 ANOVA design
#' @param dist name of the distribution to use for transforming the continuous
#'   latent z-scaled values in the population to observable raw scores. Options
#'   include "bbin" for beta-binomial, "nbin" for negative binomial, and "pois"
#'   for Poisson.
#' @param mn average raw score in the population
#' @param sd standard deviation of raw scores in the population
#' @param size integer of items/trials on the simulated tests. Only relevant when
#'   dist = "bbin"
#'
#' @returns a dataframe with n * 4 rows selected randomly from the population
#'   dataframe
#'
get_sample_aov <- function(pop, n, dist, mn, sd, size) {

  # subset population into groups by order factor
  g1 <- subset(pop, pop$order == "A")
  g2 <- subset(pop, pop$order == "B")

  # draw random ids from each group
  id1 <- sample(unique(g1$id), n)
  id2 <- sample(unique(g2$id), n)

  # combine random samples from each group
  out <- rbind(g1[which(g1$id %in% id1), ],
               g2[which(g2$id %in% id2), ])

  # apply re-scaling based on requested distribution type
  if ( dist == "bbin" ) {
    out$y <- norm2bbin(out$y, mn, sd, size)
  } else if ( dist == "nbin" ) {
    out$y <- norm2nbin(out$y, mn, sd)
  } else if ( dist == "pois" ) {
    out$y <- norm2pois(out$y, mn)
  }

  # return result
  return( out )
}
