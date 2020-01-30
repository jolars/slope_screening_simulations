library(owl)
library(progress)

source("R/screening-rules.R")

p <- 10000
n <- 200
k <- 10
q <- 0.1*n/p # slope parameter
gamma <- 1 # oscar parameter

lambda_type <- c("OSCAR", "BH", "lasso")
rho <- c(0, 0.05, 0.1, 0.2)

out <- data.frame()

pb <- progress_bar$new(total = length(lambda_type)*length(rho))

for (j in seq_along(rho)) {
  x <- matrix(rnorm(n*p), n)

  if (rho[j] != 0) {
    x <- x + sqrt(rho[j]/(1 - rho[j])) *rnorm(n)
  }

  beta <- double(p)
  beta[1:k] <- runif(k, -2, 2)
  y <- x %*% beta + rnorm(n)

  for (i in seq_along(lambda_type)) {
    pb$tick()

    lambda <- switch(lambda_type[i],
                     OSCAR = gamma*(p:1 - 1) + 1,
                     BH = qnorm(1 - (1:p)*q/(2*p)),
                     lasso = rep(1, p))

    fit <- owl(x,
               y,
               family = family,
               lambda = lambda,
               intercept = FALSE,
               screening = TRUE)

    beta_hat <- coef(fit)
    active_sets <- fit$active_sets
    n_active <- lengths(active_sets)
    n_true <- colSums(beta_hat != 0)
    true_sets <- apply(beta_hat != 0, 2, which)
    n_violations <-
      mapply(function(a, b) length(setdiff(a, b)), true_sets, active_sets)

    tmp <- data.frame(type = lambda_type[i],
                      rho = rho[j],
                      sigma = fit$sigma/max(fit$sigma),
                      active = n_active,
                      true = n_true,
                      violations = n_violations)
    out <- rbind(out, tmp)
  }
}

saveRDS(out, file.path("data", "sim_efficiency_lambda_methods.rda"))
