source("R/screening-rules.R")

library(owl)
library(rdatasets)

data <- list(
  dorothea = dorothea,
  arcene = arcene,
  golub = golub,
  gisette = gisette
)

# data <- lapply(data, function(d) {
#   list(x = d$x[, 1:10],
#        y = d$y)
# })

out <- data.frame()

for (i in 1:length(data)) {

  nm <- names(data)[i]
  x <- data[[i]]$x
  y <- data[[i]]$y

  x_scale <- apply(x, 2, norm, "2")

  x <- scale(x[, x_scale != 0], scale = x_scale[x_scale != 0])

  n <- nrow(x)
  p <- ncol(x)

  x_colnorms <- apply(x, 2, norm, "2")

  for (family in c("gaussian", "binomial")) {
    y <- data[[i]]$y

    if (family == "gaussian")
      y <- y - mean(y)
    else
      y <- sign(y - 0.5)

    cat(i, nm, family, "\n")

    fit <- owl(x,
               y,
               family = family,
               standardize_features = FALSE,
               lambda = "bh",
               q = 0.1*min(1, n/p),
               lambda_min_ratio = 0.01,
               intercept = TRUE,
               max_passes = 1e5,
               diagnostics = TRUE,
               screening = TRUE)

    n_penalties <- length(fit$sigma)
    n_kkt_checks <- lengths(fit$violations)[1:n_penalties]

    active_sets <- matrix(FALSE, p, n_penalties)

    for (j in 2:n_penalties) {
      beta_prev <- coef(fit)[-1, j-1]
      intercept_prev <- coef(fit)[1, j-1]

      lambda <- fit$lambda*fit$sigma[j]
      lambda_prev <- fit$lambda*fit$sigma[j-1]

      linear_predictor <- x %*% beta_prev + intercept_prev

      if (family == "gaussian") {
        pseudo_gradient_prev <- linear_predictor - y
      } else {
        # binomial
        pseudo_gradient_prev <- -y / (1 + exp(y * linear_predictor));
      }

      gradient_prev <- t(x) %*% pseudo_gradient_prev

      active_sets[, j] <- activeSet(x,
                                    y,
                                    lambda*nrow(x),
                                    lambda_prev*nrow(x),
                                    beta_prev,
                                    intercept_prev,
                                    gradient_prev,
                                    pseudo_gradient_prev,
                                    x_colnorms,
                                    method = "strong")

      active_sets[, j] <- active_sets[, j] | (beta_prev != 0)
    }

    beta <- coef(fit)[-1, ]

    n_violations <- colSums((!active_sets) & (beta != 0))
    n_active <- colSums(active_sets)
    n_true_active <- colSums(beta != 0)
    n_unique <- apply(beta, 2, function(x) {
      length(unique(abs(x[x != 0])))
    })

    tmp <- data.frame(dataset = nm,
                      family = family,
                      n = nrow(x),
                      p = ncol(x),
                      penalty = seq_along(fit$sigma),
                      n_violations = n_violations,
                      n_active = n_active,
                      n_true_active = n_true_active,
                      n_unique = n_unique,
                      n_kkt_checks = n_kkt_checks)
    out <- rbind(out, tmp)
  }
}

saveRDS(out, file.path("data", "sim_efficiency_violations_real_data.rda"))
