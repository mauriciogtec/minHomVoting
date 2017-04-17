library(minHomVoting)
data("onedim_voting")
data("twodim_voting")
popl_size <- 100
popl <- vector("list", popl_size)
father <- input_roll_call(onedim_voting)
plot(father)
predict(father)
compute_homology_rank(father)
for (i in 1:popl_size) {
  popl[[i]] <- mutate_roll_call(father)
}
ac <- sapply(popl, accuracy)
hr <-sapply(popl, compute_homology_rank)
hist(ac)
hist(hr)
plot(popl[[which.min(hr)]])
plot(popl[[which.max(hr)]])


init_popl <- 100
n_gen <- 100
survival <- 0.3
father_rate <- 0.3
mutation <- 0.025
popl <- vector("list", init_popl)
for (i in 1:init_popl) {
  popl[[i]] <- father
}
ac <- sapply(popl, accuracy)
for (i in 1:n_gen) {
  popl <- popl[ac == 1]
  hr <- sapply(popl, compute_homology_rank)
  best <- order(hr)[1:ceiling(survival * length(popl))]
  popl <- popl[best]
  popl_size <- length(popl)
  n_missing <- init_popl - popl_size
  for (j in 1:ceiling(father_rate * length(popl))) {
    popl <- c(popl, list(mutate_roll_call(father, mutation)))
  }
  popl_size <- length(popl)
  for (j in 1:n_missing) {
    ind <- sample.int(popl_size, 2, replace = TRUE)
    popl <- c(popl, list(combine_roll_call(
      mutate_roll_call(popl[[ind[1]]], mutation),
      mutate_roll_call(popl[[ind[2]]], mutation)
    )))
  }
  ac <- sapply(popl, accuracy)
  print(n_missing)
  print(hr)
}
hist(ac)
popl <- popl[ac == 1]
hr <- sapply(popl, compute_homology_rank)
hist(hr)
summary(hr)
plot(popl[[which.min(hr)]])



data(onedim_voting)
father <- input_roll_call(twodim_voting)

popl_size <- 25
n_gen <- 500
mutation_size <- 2
mutation_rate <- 4
crossover_size <- 2
crossover_rate <- 0.5
popl <- vector("list", popl_size)
max_attempts <- 10
min_accuracy <- 1

for (i in 1:popl_size) {
  popl[[i]] <- father
}

for (k in 1:n_gen) {
  # mutation
  lambda <- ceiling(mutation_size * popl_size)
  mutation_popl <- list()
  idx <- sample.int(popl_size, lambda, replace = TRUE)
  for (i in 1:lambda) {
    survived <- FALSE
    attempts <- 1
    while (attempts < max_attempts && !survived) {
      x <- mutate_roll_call(popl[[idx[i]]], mutation_rate)
      accuracy <- compute_accuracy(x)
      if (accuracy >= min_accuracy) {
        mutation_popl <- c(mutation_popl, list(x))
        survived <- TRUE
      }
      attempts <- attempts + 1
    }
    if (!survived) {
      mutation_popl <- c(mutation_popl, list(popl[[idx[i]]]))
    }
  }
  # crossover
  full_popl <- c(popl, mutation_popl)
  eta <- ceiling(crossover_size * popl_size)
  crossover_popl <- list()
  idx <- sample.int(popl_size + lambda, 2*eta, replace = TRUE)
  for (i in 1:eta) {
    survived <- FALSE
    attempts <- 0
    while (attempts < max_attempts && !survived) {
      x <- crossover_roll_call(
        full_popl[[idx[2*i - 1]]],
        full_popl[[idx[2*i]]],
        c("intersect", "pool"),
        crossover_rate
      )
      if (compute_accuracy(x) >= min_accuracy) {
        crossover_popl <- c(crossover_popl, list(x))
        survived <- TRUE
      }
    }
    if (!survived) {
      mutation_popl <- c(mutation_popl, list(popl[[idx[i]]]))
    }
  }

  # selection
  full_popl <- c(popl, mutation_popl, crossover_popl)
  hr <- sapply(full_popl, compute_homology_rank)
  best <- order(hr)[1:popl_size]
  popl <- full_popl[best]
  print(hr[best])
  if (min(hr) == 0) break
}
hr <- sapply(popl, compute_homology_rank)
best <- order(hr)[1:popl_size]
summary(hr)
plot(popl[[best[1]]])


# Compare to classical mds
distmat <- matrix(0, nrow = nrow(twodim_voting), ncol = nrow(twodim_voting))
for (i in 1:(nrow(distmat) - 1)) {
  for (j in i:nrow(distmat)) {
    distmat[i, j] <- distmat[j, i] <- sum(twodim_voting[i, ] != twodim_voting[j, ])
  }
}
res <- cmdscale(distmat, k = 2)
plot(res)
text(res, label = row.names(twodim_voting))


# Larger data
library(minHomVoting)
library(metodosMultivariados2017)
data("senado_votacion")
d <- senado_votacion[ ,-(1:3)]
d[is.na(d)] <- 0
d[d == -1] <- 0
father <- input_roll_call(d)



