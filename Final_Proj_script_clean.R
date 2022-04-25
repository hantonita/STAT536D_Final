require(rjags)
require(MCMCvis)


## FUNCTIONS
## -----------------------------------------------------------------------------
## Some adapted from STAT 536D lectures by Paul Gustafson

logit <- function(p) { log(p)-log(1-p) }
expit <- function(z) { 1/(1+exp(-z)) }

sim_full_data <- function(n) {
  set.seed(n)

  ### balanced case-control study
  y <- sample(c(rep(0, n/2), rep(1, n/2)))
  
  ### true OR of 2
  gamma.0.TR <- 0.1
  gamma.1.TR <- expit(logit(gamma.0.TR)+log(2))
  
  x <- rbinom(n, size=1, prob=(1-y)*gamma.0.TR + y*gamma.1.TR)
  ### independent surrogate
  xstr <- rbinom(n, size=1, prob=(1-x)*(1-0.85)+x*(0.75))
  
  dta_mnar <- data.frame(x=x, xstr=xstr, y=y)
  
  dta_mnar
}

delete_data <- function(dta_full, prop_miss_X1, prop_miss_X0) {
  dta_mnar <- tibble::rowid_to_column(dta_full, "index")
  
  dta_mnar_exp <- dta_mnar[dta_mnar$x == 1,] 
  dta_mnar_unexp <- dta_mnar[dta_mnar$x == 0,]
  
  # delete prop_miss_X1 of exposed
  n_exp_rm <- ceiling(prop_miss_X1*nrow(dta_mnar_exp))
  # prop_miss_X0 of unexposed
  n_unexp_rm <- ceiling(prop_miss_X0*nrow(dta_mnar_unexp))
  
  if (n_exp_rm > 0) {
    x_rm_exp <- sample(dta_mnar_exp$index, n_exp_rm)
    dta_mnar[x_rm_exp,]$x <- NA
  }
  
  if (n_unexp_rm > 0) {
    x_rm_unexp <- sample(dta_mnar_unexp$index, n_unexp_rm)
    dta_mnar[x_rm_unexp,]$x <- NA
  }
  
  dta_mnar <- dta_mnar[, !(names(dta_mnar) %in% c("index"))]
  dta_mnar
}

perform_simulation <- function(dta_mnar) {
  
  genmod1.string <- "model{

    ### prior distribution
    gamma.0 ~ dunif(0,1)
    gamma.1 ~ dunif(0,1)
    sens ~ dunif(0.5, 1)
    spec ~ dunif(0.5, 1)
  
    trgt <- logit(gamma.1)-logit(gamma.0)
  
    ### statistical model
    for (i in 1:n) {
      x[i] ~ dbern(pr.x[i])
      pr.x[i] <- (1-y[i])*gamma.0 + y[i]*gamma.1
         
      xstr[i] ~ dbern(pr.xstr[i])
      pr.xstr[i] <- (1-x[i])*(1-spec) + x[i]*sens
    }  
  }"
  
  ### generative model, data go in
  mod1 <- jags.model(textConnection(genmod1.string),
                     data=list(x=dta_mnar$x, y=dta_mnar$y, xstr=dta_mnar$xstr,       
                               n=dim(dta_mnar)[1]),                            
                     n.chains=3)
  
  ###  MC output comes out
  mod1.JAGS <- coda.samples(mod1, n.iter=10000, thin=10, 
                            variable.names=c("gamma.0","gamma.1","sens","spec","trgt")) 
  
  res_summary <- MCMCsummary(mod1.JAGS)
  res_summary
  
}

## ----------------------------------------------------------------------------


true_trgt=log(2)
prop_miss_X1s <- seq(0, 1, 0.10)
prop_miss_X0s <- seq(0, 1, 0.10)

trgt_summary <- matrix(nrow = 11, ncol = 11, dimnames = list(prop_miss_X0s, prop_miss_X1s))
trgt_summary <- data.frame(trgt_summary)
bias_summary <- matrix(nrow = 11, ncol = 11, dimnames = list(prop_miss_X0s, prop_miss_X1s))
bias_summary <- data.frame(bias_summary)
sd_summary <- matrix(nrow = 11, ncol = 11, dimnames = list(prop_miss_X0s, prop_miss_X1s))
sd_summary <- data.frame(sd_summary)
CI_summary <- matrix(nrow = 11, ncol = 11, dimnames = list(prop_miss_X0s, prop_miss_X1s))
CI_summary <- data.frame(CI_summary)
m.eff_summary <- matrix(nrow = 11, ncol = 11, dimnames = list(prop_miss_X0s, prop_miss_X1s))
m.eff_summary <- data.frame(m.eff_summary)
sens_summary <- matrix(nrow = 11, ncol = 11, dimnames = list(prop_miss_X0s, prop_miss_X1s))
sens_summary <- data.frame(sens_summary)
spec_summary <- matrix(nrow = 11, ncol = 11, dimnames = list(prop_miss_X0s, prop_miss_X1s))
spec_summary <- data.frame(spec_summary)

dta_full <- sim_full_data(n=1000)

for (pm1 in prop_miss_X1s) {
  
  for (pm0 in prop_miss_X0s) {
    
    cat("Proportion Missing in X=1: ", pm1, "\n")
    cat("Proportion Missing in X=0: ", pm0, "\n")
    
    dta_mnar <- delete_data(dta_full = dta_full, prop_miss_X1 = pm1, prop_miss_X0 = pm0)
    sim_results <- perform_simulation(dta_mnar)
    col <- paste0("X", pm1)
    row <- toString(pm0)
    
    trgt_summary[row, col] <- sim_results['trgt', 'mean']
    bias_summary[row, col] <- sim_results['trgt', 'mean'] - true_trgt
    sd_summary[row, col] <- sim_results['trgt', 'sd']
    CI_summary[row, col][[1]] <- list(c(sim_results['trgt', '2.5%'], sim_results['trgt', '97.5%']))
    m.eff_summary[row, col] <- sim_results['trgt', 'n.eff']
    sens_summary[row, col] <- sim_results['sens', 'mean']
    spec_summary[row, col] <- sim_results['spec', 'mean']
  }
  
}

trgt_summary
bias_summary
sd_summary
CI_summary
m.eff_summary
sens_summary
spec_summary

save(trgt_summary, file = paste0("Data/clean/simple/res_trgt.RData"))
save(bias_summary, file = paste0("Data/clean/simple/res_bias.RData"))
save(sd_summary, file = paste0("Data/clean/simple/res_sd.RData"))
save(CI_summary, file = paste0("Data/clean/simple/res_CI.RData"))
save(m.eff_summary, file = paste0("Data/clean/simple/res_meff.RData"))
save(sens_summary, file = paste0("Data/clean/simple/res_sens.RData"))
save(spec_summary, file = paste0("Data/clean/simple/res_spec.RData"))

