# This script contains formula for computing effect sizes
library(testit)
library(gsubfn)


# Helper function for testing
sample_n <- function() sample.int(140, size = 1) + 10

# Between participants ----
# Compute pooled SD
sd_pooled <- function(n1, n2, sd1, sd2) 
  sqrt(((n1 - 1) * sd1 ^ 2 + (n2 - 1) * sd2 ^ 2) / (n1 + n2 - 2))

# Compute variance of d
bw_var_d <- function(n1, n2, d) 
  ((n1 + n2) / (n1 * n2)) + 
  ((d ^ 2) / (2 * (n1 + n2)))


# Compute d from means and sds
bw_abs_d_from_means_sds <- function(n1, n2, m1, m2, sd1, sd2){
  d <- abs((m1 - m2) / sd_pooled(n1, n2, sd1, sd2))
  vard <- bw_var_d(n1, n2, d)
  
  return(list(d = d, vard = vard, r = NA))
}

# Compute d from t
bw_abs_d_from_t <- function(n1, n2, t) {
  d <- abs(t / sqrt(1 / (1 / n1 + 1 / n2)))
  vard <- bw_var_d(n1, n2, d)
  
  return(list(d = d, vard = vard, r= NA))
}

# Compute d from F
bw_abs_d_from_F <- function(n1, n2, F_val) {
  d <- sqrt(F_val /  (1 / (1 / n1 + 1 / n2)))
  vard <- bw_var_d(n1, n2, d)
  
  return(list(d = d, vard = vard, r = NA))
}


# Test that the formulas work
check_between <- function(n1, n2, tol = 0.00001){
  
  # Draw group means
  mu1 <- rnorm(1)
  mu2 <- rnorm(1)
  
  # Draw group SDs
  sigma1 <- abs(rnorm(1))
  sigma2 <- abs(rnorm(1))
  
  # Random data
  x1 <- rnorm(n1, mean = mu1, sd = sigma1)
  x2 <- rnorm(n2, mean = mu2, sd = sigma2)
  
  # d from means
  d1 <- bw_abs_d_from_means_sds(n1, n2, mean(x1), mean(x2), sd(x1), sd(x2))$d
  
  # t value
  t <- t.test(x1, x2, var.equal = T)$statistic
  
  # d from t
  d2 <- bw_abs_d_from_t(n1, n2, t)$d
  
  # F value
  F_val <- t ^ 2
  
  # d from F
  d3 <- bw_abs_d_from_F(n1, n2, F_val)$d
  
  assert("The three between d formulas don't return the same value",
         (abs(d1 - d2) < tol) & (abs(d1 - d3) < tol))
}

check_between(sample_n(), sample_n())

# Within participant ----

## Extract r for within subject designs ----

extract_r_from_sdd_sds <- function(sd_d, sd1, sd2) {
  cov <- (sd1 ^ 2 + sd2 ^ 2 - sd_d ^ 2) / 2
  r <- cov / (sd1 * sd2)
  
  return(r)
}

extract_r_from_t_means_sds <- function(n, t, m1, m2, sd1, sd2) {
  sed <- abs((m1-m2) / t)
  sd_d <- sed * sqrt(n)
  r <- extract_r_from_sdd_sds(sd_d, sd1, sd2)

  return(r)
}

# extract_r_from_F_means_sds <- function(n, F_val, m1, m2, sd1, sd2) {
#   t <- sqrt(F_val)
#   
#   return(extract_r_from_t_means_sds(n, t, m1, m2, sd1, sd2))
# }


### Functions to test r extraction ----

draw_within_data <- function(n) {
  is_positive_definite <- F
  
  while (is_positive_definite == F) {
    # Draw group means
    mu1 <- rnorm(1)
    mu2 <- rnorm(1)
    
    # Draw group SDs
    sigma1 <- abs(rnorm(1))
    sigma2 <- abs(rnorm(1))
    
    # Draw correlation
    rho <-  rbeta(1, 2, 2)
    
    # Compute covariace
    cov <- rho * sigma1^2 * sigma2^2
    
    # Construct covariance matrix
    Sigmasq <- matrix(c(sigma1 ^ 2, cov, cov, sigma2 ^ 2), 2, 2)
    
    # Test if positive definite
    eS <- eigen(Sigmasq, symmetric = TRUE)
    ev <- eS$values
    
    is_positive_definite <- all(ev >= -1e-06 * abs(ev[1L]))
  }
  
  # Draw data
  x <-
    MASS::mvrnorm(n = n,
                  mu = c(mu1, mu2),
                  Sigma = Sigmasq)
  
  return(x)
}

summarize_within <- function(x) {
  
  # Compute sample means
  m1 <- mean(x[,1])
  m2 <- mean(x[,2])
  
  # Compute sample standard deviations
  sd1 <- sd(x[,1])
  sd2 <- sd(x[,2])
  
  # Compute sample SD of difference
  sd_d <- sd(x[,1] - x[,2])
  
  # Compute r
  r <- cor(x[,1], x[,2])
  
  # Compute t value
  t <- t.test(x[,1] - x[,2], var.equal = T)$statistic
  
  return(list(m1, m2, sd1, sd2, sd_d, r, t))
}


check_extract_r <- function(n, tol = 0.00001){
  
  # Draw data
  x <- draw_within_data(n)
  
  # Summarize data
  list[m1, m2, sd1, sd2, sd_d, r, t] = summarize_within(x)
  
  # Extract r from SD_d and SDs
  r_from_sdd_sds <- extract_r_from_sdd_sds(sd_d, sd1, sd2)
  
  # Extract r from t, means and sds
  r_from_t_means_sds <- extract_r_from_t_means_sds(n, t, mean(x[,1]), 
                                                   mean(x[,2]),
                                                   sd1,
                                                   sd2)
  
  assert("Extracting r formulas don't agree with true r",
         (abs(r - r_from_sdd_sds) < tol) & 
           (abs(r - r_from_t_means_sds) < tol))
}

check_extract_r(sample_n())

## Compute within-subject effect sizes ----
# https://cloud.r-project.org/web/packages/TOSTER/vignettes/SMD_calcs.html#cohens-drm-correlation-corrected
# https://rdrr.io/cran/TOSTER/src/R/cohend_calcs.R#sym-get_ncp_t 
w_abs_d <- function(n, m1, m2, sd1, sd2, r, method, t) {
  
  if (method == "rm_toster_default") {
    # Using Lakens' package. W/o Hedge's g correction
    out <- TOSTER:::d_est_pair(n, m1, m2, sd1, sd2, r, type = "d", denom = "rm")
    # this is the same as: TOSTER:::d_est_pair(n, m1, m2, sd1, sd2, r, type = "d", denom = "rm", smd_ci = "goulet")
    # using the Cousineau and Goulet-Pelletier (2021) 
    d <- abs(out$d)
    
    se_d <- out$d_sigma
    vard = se_d ^ 2
    
  } else if (method == "rm_stackx") {
    # https://stats.stackexchange.com/questions/256053/variance-of-cohens-d-for-within-subjects-designs and https://aaroncaldwell.us/TOSTERpkg/articles/SMD_calcs.html 
    d <- (abs(m1-m2)/sqrt(sd1^2 + sd2^2 - 2*r*sd1*sd2)) * sqrt(2*(1-r))
    
    vard <- ((1/n + d^2/(2*n))*2*(1-r))
    
  } else if (method == "rm_toster_eq_nct") { 
    # noncentral t-distribution. this will approximately match the results of Buchanan et al. (2019) and Ben-Shachar, Lüdecke, and Makowski (2020)
    d <- (abs(m1-m2)/sqrt(sd1^2 + sd2^2 - 2*r*sd1*sd2)) * sqrt(2*(1-r))
    
    df = n-1

    lambda <- 1/n # non-centrality parameter
    
    vard <-((df/(df-2))*((2*(1-r))/n) * (1+d^2*(n/(2*(1-r)))) - (d^2/1^2))

    #   
    # } else if (method == "rm_toster_eq_nct_hedges") { 
    #   
    #   d_orig <- ((m1-m2)/sqrt(sd1^2 + sd2^2 - 2*r*sd1*sd2)) * sqrt(2*(1-r))
    #   
    #   df = n-1
    #   
    #   lambda <- 1/n # non-centrality parameter - but i think this is used to for the CIs
    #   
    #   J < exp ( gamma(df/2) - log(sqrt(df/2)) - gamma ((df-1) /2) )
    #   
    #   d <- d_orig*J
    #   
    #   se_d <-sqrt((df/(df-2))*((2*(1-r))/n) * (1+d^2*(n/(2*(1-r)))) - (d^2/J^2))
    
  } else if (method == "rm_toster_nct") { 
    ## same as rm_toster_eq_nct
    out <- TOSTER:::d_est_pair(n, m1, m2, sd1, sd2, r, type = "d", denom = "rm", smd_ci = "nct")
    
    d <- abs(out$d)
    
    se_d <- out$d_sigma
    vard = se_d ^ 2
    
    # this method uses this variance calc: 
    # df=n-1
    # Calculate the standard error (σ_SMD)
    # sigma_SMD <- sqrt(
    #   (df / (df - 2)) * 
    #     ((2 * (1 - r)) / n) * 
    #     (1 + (d_rm^2 * (n / (2 * (1 - r))))) - (d_rm^2/1^2)
    # )
    
    
  # } else if (method == "av_cheung") {
  #   # Using Mike Cheung's method 
  #   # https://stats.stackexchange.com/questions/256053/variance-of-cohens-d-for-within-subjects-designs
  #   
  #   ## Arrange sds and r into sample covariance matrix
  #   R <- matrix(c(1, r, r, 1), nrow = 2)
  #   Cov <- cor2cov(R, c(sd1, sd2), names=c("x_drug","x_placebo"))
  #   
  #   ## Arrange means
  #   Mean <- c(m1, m2, 13)
  #   
  #   # Lavaan model
  #   model1 <- '# Label the sds with sd1 and sd2
  #          eta_pre =~ sd1*x_drug
  #          eta_post =~ sd2*x_placebo
  #          # Fix the error variances at 0
  #          x_drug ~~ 0*x_drug
  #          x_placebo ~~ 0*x_placebo
  #          # Label the means with m1 and m2 
  #          x_drug ~ m1*1
  #          x_placebo ~ m2*1
  #          # Define d_av
  #          d_av := (m2-m1)/((sd1+sd2)/2)'
  #   
  #   ## Fit the model
  #   fit1 <- cfa(model1, sample.cov=Cov, sample.mean=Mean, 
  #               sample.nobs=n, std.lv=TRUE, 
  #               sample.cov.rescale=FALSE)
  #   
  #   ## Extract d_av and its SE
  #   d_av <- parameterEstimates(fit1)[12, c(5,6)]
  #   
  #   d <- abs(d_av$est)
  #   
  #   se_d <- d_av$se
  #   vard = se_d^2
    
  # } else if (method == "av_algina") { # OLD VERSION: NOT CORRECT: I THINK THIS IS BIRD (2003)
    # https://www.sciencedirect.com/science/article/pii/S0149763417307674?via%3Dihub#fn0015
    # https://journals.sagepub.com/doi/pdf/10.1177/0013164403256358
    # d <- abs(m1-m2)/sqrt((sd1^2+sd2^2)/2)
    # 
    # alpha <- 0.05  # Significance level
    # 
    # # Calculate t-value
    # t_value <- qt(1 - alpha / 2, df = n - 1)
    # 
    # # Calculate the standard error (SE)
    # # SE <- sqrt((2 * (sd1^2 + sd2^2 - 2 * r)) / (n * (sd1^2 + sd2^2))) # incorrectly specified in the Hirst et al paper.
    # SE <- sqrt((2 * (sd1^2 + sd2^2 - 2 * r * sd1 * sd2)) / (n * (sd1^2 + sd2^2))) # correct from the original paper
    # 
    # 
    # # Calculate Cohen's d_av confidence interval
    # CI_low <- d - t_value * SE
    # CI_up <- d + t_value * SE
    # 
    # # Calculate the variance of dav
    # vard <- ((CI_up - CI_low) / (2 * 1.96))^2

    
  } else if (method == "av_algina") { # THIS IS THE CORRECT CALCULATION, derived from Cousineau & Goulet-Pelletier 2018 paper
    d <- abs(m1-m2)/sqrt((sd1^2+sd2^2)/2)
    
    # W <- geometric.mean(c(sd1^2, sd2^2)) / mean(c(sd1^2, sd2^2))
    # rW <- r * W
    # tCI <- conf.limits.nct(d * sqrt(n/(2*(1-rW))), n-1, conf.level = 0.95)
    df = n-1
    tstat <- d * sqrt((n*(sd1^2+sd2^2))/ (2*(sd1^2+sd2^2-2*sd1*sd2*r))) # the noncentrality parameter (e.g., observed t-value) of interest.
    tCI <- conf.limits.nct(tstat, df, conf.level = 0.95)
    
    tCI.low <- tCI$Lower.Limit
    tCI.hig <- tCI$Upper.Limit
    scaling_factor <- sqrt( 2 * (sd1^2 + sd2^2 - 2 * sd1 * sd2 * r) / (n * (sd1^2 + sd2^2)) )

    limits <- c(tCI.low, tCI.hig) / sqrt((n*(sd1^2+sd2^2))/ (2*(sd1^2+sd2^2-2*sd1*sd2*r)))
    # limits <- c(tCI.low, tCI.hig) * scaling_factor # HOW IT'S SPECIFIED IN ALGINA & KESERMAN
    # limits <- c(tCI.low, tCI.hig) * sqrt((2*(1-rW))/n) # SIMPLIFIED IN Cousineau 2021
    
    
    # tCI.low <- tCI$Lower.Limit
    # tCI.hig <- tCI$Upper.Limit
    # limits <- c(tCI.low, tCI.hig) / sqrt(n/(2*(1-rW)))

    # vard <- ((limits[2] - limits[1]) / (2 * 1.96))^2
    vard <- ((limits[2] - limits[1]) / (2 * 1.96))^2
    
  } else if (method == "rm_algina") { # SAME AS ABOVE BUT USES d_rm as d
    d <- (abs(m1-m2)/sqrt(sd1^2 + sd2^2 - 2*r*sd1*sd2)) * sqrt(2*(1-r))
    
    W <- geometric.mean(c(sd1^2, sd2^2)) / mean(c(sd1^2, sd2^2))
    rW <- r * W
    tCI <- conf.limits.nct(d * sqrt(n/(2*(1-rW))), n-1, conf.level = 0.95)
    tCI.low <- tCI$Lower.Limit
    tCI.hig <- tCI$Upper.Limit
    limits <- c(tCI.low, tCI.hig) / sqrt(n/(2*(1-rW)))
    
    vard <- ((limits[2] - limits[1]) / (2 * 1.96))^2
    
    
  } else if (method == "dz") {
    d = abs(t)/sqrt(n)

    vard = (1/n) + d^2/(2*n)
    
    # out_z <- TOSTER:::d_est_pair(n, m1, m2, sd1, sd2, r, type = "d", denom = "z", smd_ci = "nct")
    # 
    # d <- abs(out_z$d)
    # vard =out_z$d_sigma^2
  }
  
  return(list(
    d = d,
    vard = vard
  ))
}


# Compute d_rm from t
w_abs_d_from_t <- function(n, m1, m2, sd1, sd2, t, method, r_correlation) {
  res <- w_abs_d(n, m1, m2, sd1, sd2, r_correlation, method, t)
  return(res)
}

# Compute d_rm from F
w_abs_d_from_F <- function(n, m1, m2, sd1, sd2, F_val, method, r_correlation, t) {
  res <- w_abs_d(n, m1, m2, sd1, sd2, r_correlation, method, t)
  return(res)
}

check_within <- function(n, tol = 0.00001, method){
  
  # Draw within data
  x <- draw_within_data(n)
  
  # Summarize data
  list[m1, m2, sd1, sd2, sd_d, r, t] = summarize_within(x)
  t=NA
  # d from means
  list[d1, vard1] <- w_abs_d(n, m1, m2, sd1, sd2, r, method, t)
  
  # d from t
  list[d2, vard2] <- w_abs_d_from_t(n, m1, m2, sd1, sd2, t, method, r_correlation)
  
  # F value
  F_val <- t ^ 2
  
  # d from F
  list[d3, vard3] <- w_abs_d_from_F(n, m1, m2, sd1, sd2, F_val, method, r_correlation, tstat)
  
  
  assert("The three within d formulas don't return the same value",
         (abs(d1 - d2) < tol) & (abs(d1 - d3) < tol))
}

# check_within(sample_n(), method = "rm")
# 
# check_within(sample_n(), method = "av")

check_within_method_correlation <- function(n_studies){
  
  one_study <- function(id){
    
    # Draw study size
    n <- sample_n()
    
    # Draw within data
    x <- draw_within_data(n)
    t = NA
    # Summarize data
    list[m1, m2, sd1, sd2, sd_d, r, t] = summarize_within(x)
    
    # Compute d_rm (the one originally specified by Yaniv)
    d_rm <- w_abs_d(n, m1, m2, sd1, sd2, r, "rm_toster_default", t)
    
    # Compute d_av (cheung)
    # d_av <- w_abs_d(n, m1, m2, sd1, sd2, r, "av_cheung", t)
    
    return(data.frame(id = id,
                      n = n,
                      r = r,
                      sd1 = sd1,
                      sd2 = sd2,
                      abs_diff = abs(m1 - m2),
                      d_rm = d_rm$d,
                      vard_rm = d_rm$vard,
                      d_av = d_av$d,
                      vard_av = d_av$vard))
  }
  
  # Create data.frame of simulated studies
  studies <- do.call(rbind, lapply(1:n_studies, one_study))
  
  ## Plots for checking the properties of d_rm and d_av against each other
  ggplot(studies, aes(x = d_rm, y = d_av)) + geom_point() +
    geom_abline(linetype = 2) +
    xlim(0,1) + ylim(0,1)
  
  ggplot(studies, aes(x = r, y = d_av-d_rm)) + geom_point()
  
  ggplot(studies, aes(x = sd1 / sd2, y = d_av-d_rm)) + geom_point()
  
  ggplot(studies, aes(x = vard_rm, y = vard_av)) + geom_point()
  
  assert(paste0("Correlation between computed d_rm and d_av is only r=", 
                round(cor(studies$d_rm, studies$d_av), digits = 2)),
         cor(studies$d_rm, studies$d_av) > 0.85)
  
  assert(paste0("Correlation between computed SEs of d_rm and d_av is only r=", 
                round(cor(sqrt(studies$vard_rm), sqrt(studies$vard_av)), digits = 2)),
         cor(sqrt(studies$vard_rm), sqrt(studies$vard_av)) > 0.7)
  return(studies)
}

# check_within_method_correlation(200)

compare_effectsize_calcs <- function(n1, n2, m1, m2, sd1, sd2, r, design, tstat){
  if (design == 'w') { 
    # Compute d_rm (the one originally specified by Yaniv)
    d_rm <- w_abs_d(n1, m1, m2, sd1, sd2, r, "rm_toster_default", tstat)
    
    # Compute d_av (cheung)
    # d_av_cheung <- w_abs_d(n1, m1, m2, sd1, sd2, r, "av_cheung", tstat)
    
    # Compute rm_stackx 
    d_rm_stackx <- w_abs_d(n1, m1, m2, sd1, sd2, r, "rm_stackx", tstat)
    
    # Compute rm_toster_eq_nct 
    d_rm_toster_eq_nct <- w_abs_d(n1, m1, m2, sd1, sd2, r, "rm_toster_eq_nct", tstat)
    
    # Compute rm_toster_nct 
    d_rm_toster_nct <- w_abs_d(n1, m1, m2, sd1, sd2, r, "rm_toster_nct", tstat)
    
    # Compute av_algina
    d_av_algina <- w_abs_d(n1, m1, m2, sd1, sd2, r, "av_algina", tstat)
    
    # Compute rm_algina
    d_rm_algina <- w_abs_d(n1, m1, m2, sd1, sd2, r, "rm_algina", tstat)
    
    
    # Compute dz
    if (!is.na(tstat)) { 
      dz<- w_abs_d(n1, m1, m2, sd1, sd2, r, "dz", tstat)
      d_z = dz$d
      vard_z = dz$vard
    } else {
      d_z = NA
      vard_z = NA
    }
    
  }
  # Compute between
  d_bw <- bw_abs_d_from_means_sds(n1, n2, m1, m2, sd1, sd2)
  
  
  if (design == 'bw') {
    return(list( d_rm_orig = NA,
                 vard_rm_orig = NA,
                 d_rm_stackx = NA,
                 vard_rm_stackx = NA,
                 d_rm_toster_eq_nct = NA,
                 vard_rm_toster_eq_nct = NA,
                 d_rm_toster_nct = NA,
                 vard_rm_toster_nct = NA,
                 # d_av_cheung = NA,
                 # vard_av_cheung = NA, 
                 d_av_algina = NA,
                 vard_av_algina = NA,
                 d_rm_algina = NA,
                 vard_rm_algina = NA,
                 d_bw = d_bw$d,
                 vard_bw = d_bw$vard,
                 d_z = NA,
                 vard_z = NA))
  }else {
    return(list( d_rm_orig = d_rm$d,
                 vard_rm_orig = d_rm$vard,
                 d_rm_stackx = d_rm_stackx$d,
                 vard_rm_stackx = d_rm_stackx$vard,
                 d_rm_toster_eq_nct = d_rm_toster_eq_nct$d,
                 vard_rm_toster_eq_nct = d_rm_toster_eq_nct$vard,
                 d_rm_toster_nct = d_rm_toster_nct$d,
                 vard_rm_toster_nct = d_rm_toster_nct$vard,
                 # d_av_cheung = d_av_cheung$d,
                 # vard_av_cheung = d_av_cheung$vard, 
                 d_av_algina = d_av_algina$d,
                 vard_av_algina = d_av_algina$vard,
                 d_rm_algina = d_rm_algina$d,
                 vard_rm_algina = d_rm_algina$vard,
                 d_bw = d_bw$d,
                 vard_bw = d_bw$vard,
                 d_z = d_z,
                 vard_z = vard_z))
    
  }
}






# Helper functions to unify computation for within and bw
abs_d_from_means_sds <- function(design,
                                 n1,
                                 n2,
                                 m1,
                                 m2,
                                 sd1,
                                 sd2){
  if (design == "bw"){
    return(bw_abs_d_from_means_sds(n1,
                                   n2,
                                   m1,
                                   m2,
                                   sd1,
                                   sd2))
  } else {
    return(list(d=NA, vard=NA)) # This method doesn't work for within
  }
} 


abs_d_from_t <- function(design,
                         n1,
                         n2,
                         m1,
                         m2,
                         sd1,
                         sd2,
                         t,
                         w_method,
                         r_correlation,
                         tstat){
  
  if(design == "bw") {
    return(bw_abs_d_from_t(n1, n2, tstat))
  } else {
    return(w_abs_d_from_t(n1, m1, m2, sd1, sd2, t, w_method, r_correlation))
  }
}

abs_d_from_F <- function(design,
                         n1,
                         n2,
                         m1,
                         m2,
                         sd1,
                         sd2,
                         F_val,
                         w_method,
                         r_correlation,
                         tstat){
  if(design == "bw") {
    return(bw_abs_d_from_F(n1, n2, F_val))
  } else {
    return(w_abs_d_from_F(n1, m1, m2, sd1, sd2, F_val, w_method, r_correlation, tstat))
  }
}