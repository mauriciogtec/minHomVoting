library(pscl)
library(minHomVoting)
library(wnominate)
data("onedim_voting")
data("sen90")
# rc <- rollcall(onedim_voting, codes = list(yea = 1, nay = 0))
# rc <- c(rc, list())
t <- sen90$votes
t[(t %in% 1:3)] <- 1
t[(t %in% 4:6)] <- 0
t[!(t %in% 0:1)] <- NA
set.seed(110104)
t <- t[sample.int(nrow(t), 100), sample.int(ncol(t), 50)]
t <- na.omit(t)
info <- sen90$legis.data[match(row.names(t) , row.names(sen90$legis.data)), ]


father <- input_roll_call(t)
plot(father)
predict(father)
compute_homology_rank(father)



init_popl <- 500
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
