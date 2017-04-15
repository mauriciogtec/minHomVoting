
data("onedim_voting")
data("twodim_voting")
popl_size <- 100
popl <- vector("list", pop_size)
father <- input_roll_call(twodim_voting)
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
survival <- 0.75
mutation <- 0.025
popl <- vector("list", init_popl)
for (i in 1:init_popl) {
  popl[[i]] <- father
}
ac <- sapply(popl, accuracy)
for (i in 1:n_gen) {
  die <- ac < 1
  popl <- popl[!die]
  hr <- sapply(popl, compute_homology_rank)
  best <- order(hr)[1:ceiling(survival * length(popl))]
  popl <- popl[best]
  popl_size <- length(popl)
  n_missing <- init_popl - popl_size
  for (j in 1:n_missing) {
    ind <- sample.int(popl_size, 2, replace = TRUE)
    popl <- c(popl, list(combine_roll_call(
      popl[[ind[1]]],
      mutate_roll_call(popl[[ind[2]]], mutation)
    )))
  }
  ac <- sapply(popl, accuracy)
}
hist(ac)
popl <- popl[ac == 1]
hr <-sapply(popl, compute_homology_rank)
hist(hr)
plot(popl[[which.min(hr)]])
