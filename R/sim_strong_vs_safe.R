source("R/screening-rules.R")

library(owl)

p <- 5000
n <- 200

out <- data.frame()

set.seed(5421)

for (rho in c(0, 0.1, 0.4)) {
  x <- matrix(rnorm(n*p), n)

  if (rho != 0) {
    x <- x + sqrt(rho/(1 - rho))*rnorm(n)
  }

  beta <- double(p)
  k <- p/4
  beta[1:k] <- rnorm(k)
  y <- x %*% beta + rnorm(n)

  y <- y - mean(y)

  x_scale <- apply(x, 2, norm, "2")
  x <- scale(x, scale = x_scale)

  fit <- owl(x,
             y,
             screening = TRUE,
             intercept = FALSE,
             standardize_features = FALSE,
             lambda = "bh",
             diagnostics = TRUE,
             q = 0.1*min(1, n/p))

  beta_hat <- coef(fit)
  sigma <- fit$sigma*nrow(x)
  lambda <- fit$lambda

  n_lambda <- length(lambda)
  n_sigma <- length(sigma)

  intercept_prev <- 0

  active_sets <- matrix(FALSE, p, n_sigma)

  x_colnorms <- apply(x, 2, norm, "2")

  for (method in c("safe", "strong")) {
    for (m in 2:length(sigma)) {
      beta_prev <- beta_hat[, m-1]

      lambda <- fit$lambda*sigma[m]
      lambda_prev <- fit$lambda*sigma[m-1]
      pseudo_gradient_prev <-  x %*% beta_prev - y
      gradient_prev <- t(x) %*% pseudo_gradient_prev

      active_sets[, m] <- activeSet(x,
                                    y,
                                    lambda,
                                    lambda_prev,
                                    beta_prev,
                                    intercept_prev,
                                    gradient_prev,
                                    pseudo_gradient_prev,
                                    x_colnorms,
                                    method = method)
      prev_active <- beta_prev != 0
      active_sets[, m] <- active_sets[, m] | prev_active
    }

    n_violations <- colSums(!active_sets & (beta_hat != 0))
    active_screened <- apply(active_sets, 2, sum)
    active_true <- apply(beta_hat != 0, 2, sum)

    tmp <- data.frame(method = method,
                      sigma = sigma/max(sigma),
                      active_screened = active_screened,
                      active_true = active_true,
                      n_violations = n_violations,
                      rho = rho)
    out <- rbind(out, tmp)
  }
}

saveRDS(out, file.path("data", "sim_strong_vs_safe.rda"))


