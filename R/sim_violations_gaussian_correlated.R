source("R/screening-rules.R")

library(owl)
library(rdatasets)
library(progress)

n <- 100
p <- c(20, 50, 100, 500, 1000)
n_rep <- 100
rho <- 0.5

n_p <- length(p)

set.seed(733)

out <- data.frame()

pb <- progress_bar$new(total = n_p*n_rep,
                       format = "[:bar] :percent eta: :eta")

for (i in 1:n_p) {
  for (j in 1:n_rep) {
    pb$tick()

    # simulate correlated data
    x <- matrix(rnorm(n*p[i]), n)
    x <- x + sqrt(rho/(1 - rho)) * matrix(rnorm(n), n, p[i])

    k <- floor(p[i]/4)

    beta <- double(p[i])
    beta[1:k] <- sample(c(-2, 2), k, replace = TRUE)

    y <- x %*% beta + rnorm(n)

    x_scale <- apply(x, 2, norm, "2")
    x <- scale(x, scale = x_scale)
    y <- y - mean(y)

    fit <- owl(x,
               y,
               family = "gaussian",
               standardize_features = FALSE,
               lambda = "bh",
               q = 0.01,
               lambda_min_ratio = 1e-2,
               tol_dev_change = 0,
               tol_dev_ratio = 1,
               max_variables = p[i]*2, # always full path
               intercept = FALSE,
               screening = TRUE)

    sigma <- fit$sigma

    beta_hat <- coef(fit)

    y_hat <- x %*% beta_hat

    # collect residual sums of squares and r2
    rss <- colSums((y_hat - matrix(y, nrow = n, ncol = ncol(y_hat)))^2)
    r2 <- 1 - rss/sum((y - mean(y))^2)

    n_penalties <- length(sigma)

    active_sets <- matrix(FALSE, p[i], n_penalties)

    x_colnorms <- apply(x, 2, norm, "2")

    for (method in c("safe", "strong")) {
      for (m in 2:length(sigma)) {
        beta_prev <- beta_hat[, m-1]
        intercept_prev <- 0

        lambda <- fit$lambda*sigma[m]
        lambda_prev <- fit$lambda*sigma[m-1]
        pseudo_gradient_prev <- -(y - x %*% beta_prev)
        gradient_prev <- t(x) %*% pseudo_gradient_prev

        active_sets[, m] <- activeSet(x,
                                      y,
                                      lambda*n,
                                      lambda_prev*n,
                                      beta_prev,
                                      intercept_prev,
                                      gradient_prev,
                                      pseudo_gradient_prev,
                                      x_colnorms,
                                      method = method)

        active_sets[, m] <- active_sets[, m] | beta_prev != 0
      }

      n_violations <- colSums(!active_sets & (beta_hat != 0))
      sigma_ratio <- sigma/max(sigma)

      tmp <- data.frame(p = p[i],
                        method = method,
                        j = j,
                        sigma_ratio = signif(sigma_ratio, 3),
                        n_violations = n_violations)
      out <- rbind(out, tmp)
    }
    out
  }
}

saveRDS(out, file.path("data", "sim_violations_gaussian_correlated.rda"))
