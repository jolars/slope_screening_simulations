# simulations for testing performance of screening rule as ratio between n and p vary

library(owl)

n <- 1000

out <- data.frame()

n_it <- 10
n_ratio <- 10
min_ratio <- 0.01
max_ratio <- 10
ratio <- signif(exp(seq(log(min_ratio),
                        log(max_ratio),
                        length.out = n_ratio)), 1)
q <- 0.1

screening <- c(TRUE, FALSE)

set.seed(1044)

iter <- 0

for (j in seq_along(ratio)) {
  p <- n*ratio[j]
  k <- floor(p*q)

  for (i in seq_len(n_it)) {
    x <- matrix(rnorm(n*p), n)
    beta <- double(p)
    beta[1:k] <- seq(2, 1, length.out = k)

    y <- x %*% beta + rnorm(n)

    for (k in 1:length(screening)) {
      iter <- iter + 1

      cat("iter:", iter, "/", n_it*length(ratio)*2, "\t",
          "ratio:", ratio[j], "\t",
          "screening:", screening[k], "\n")

      time <- system.time({
        fit <- owl(x,
                   y,
                   lambda = "bh",
                   q = q,
                   intercept = FALSE,
                   screening = screening[k])
      })

      tmp <- data.frame(n = n,
                        p = p,
                        ratio = ratio[j],
                        screening = screening[k],
                        time = time[3])
      out <- rbind(out, tmp)
    }
  }
}

saveRDS(out, file.path("data", "sim_performance_np_ratio.rda"))
