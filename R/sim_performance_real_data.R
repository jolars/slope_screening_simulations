library(owl)
library(rdatasets)

data <- list(
  cpusmall = rdatasets::cpusmall,
  physician = rdatasets::physician,
  golub = rdatasets::golub,
  zipcode = rdatasets::zipcode
)

# use only first 200 observations of zipcode set
set.seed(1123)
sind <- sample(length(data$zipcode$y), 200)
data$zipcode <- list(x = data$zipcode$x[sind, ],
                     y = data$zipcode$y[sind])

# data <- lapply(data, function(d) {
#   list(x = d$x[, 1:5],
#        y = d$y)
# })

out <- data.frame()

iter <- 0

for (i in 1:length(data)) {

  dataset <- names(data)[i]
  x <- data[[i]]$x
  y <- data[[i]]$y

  family <- switch(dataset,
                   cpusmall = "gaussian",
                   golub = "binomial",
                   abalone = "poisson",
                   physician = "poisson",
                   zipcode = "multinomial")

  n <- nrow(x)
  p <- ncol(x)

  n_lambda <- switch(family, multinomial = p*(length(unique(y)) - 1), p)

  for (screening in c(TRUE, FALSE)) {
    iter <- iter + 1

    cat("iter:", iter, "/", length(data)*2,
        "\tdata:", dataset, family,
        "\tscreening:", screening, "\n")

    time <- system.time({
      fit <- owl(x,
                 y,
                 family = family,
                 lambda = "bh",
                 q = 0.1*min(1, n/p),
                 screening = screening)
    })

    tmp <- data.frame(dataset = dataset,
                      family = family,
                      screening = screening,
                      n = n,
                      p = p,
                      time = time[3])
    out <- rbind(out, tmp)
  }
}

rownames(out) <- NULL

saveRDS(out, file.path("data", "sim_performance_real_data.rda"))
