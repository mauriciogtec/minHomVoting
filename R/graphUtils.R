#' @title all nodes
#' @description creates the nodes of voters and issues for the roll_call class graph
#' @details the roll_call class is made of a graph which consists of a noder per voter and two
#' nodes for each issue voted representing the yes and no alternatives.
#' @examples
#' data(onedim_voting)
#' gen_nodes(onedim_voting)
#' @export
gen_nodes <- function(d) {
  label = c(
    row.names(d),
    paste0("y", 1:ncol(d)),
    paste0("n", 1:ncol(d))
  )
  list(
    voter = 1:nrow(d),
    y = (nrow(d) + 1):(nrow(d) + ncol(d)),
    n = (nrow(d) + ncol(d) + 1):(nrow(d) + 2*ncol(d)),
    label = label
  )
}



#' @title initial adj matrix
#' @description receives roll call data and creates an adj matrix
#' this matrix creates an adjacency matrix from roll call data that has perfect prediction
#' but usually has the highest possible dimension.
#' @export
gen_adj_matrix <- function(d) {
  nodes <- gen_nodes(d)
  n <- length(nodes$label)
  adj_mat <- matrix(
    0,
    nrow = n,
    ncol =  n,
    dimnames = list(nodes$label, nodes$label)
  )
  for (i in 1:nrow(d)) {
    votes_yes <- as.logical(d[i, ])
    adj_mat[i, nodes$y[votes_yes]] <- 1
    adj_mat[i, nodes$n[!votes_yes]] <- 1
  }
  adj_mat
}


#' @title voting class constructor
#' @description defines the class \code{roll_call}
#' @details returns a list of class {roll_call} with elements
#' @param x a data.frame with each row representing a voter and each column the result of one roll call vote
#' @return
#' \itemize{
#'     \item \code{nodes} names of all the nodes
#'     \item \code{voter} the indexes of voters in \code{nodes}
#'     \item \code{Y} the indexes of 'yes' nodes for each issue in \code{nodes}
#'     \item \code{N} the indexes of 'no' nodes for each in \code{nodes}
#' }
#' @export
input_roll_call <- function(d) {
  nodes <- gen_nodes(d)
  adj_mat <- gen_adj_matrix(d)
  li <- c(nodes, list(votes = d, adj_mat = adj_mat))
  class(li) <- "roll_call"
  li
}

#' @title homology rank (first betti number)
#' @description computes the main fitness function, the rank of the first homology group
#' @param roll_call a list of class roll call, assumed to be generated using input_roll_call
#' @return the first betti number of the graph
#' @export
compute_homology_rank <- function(roll_call) {
  inc <- adj_to_inc(roll_call$adj_mat)
  nullity <- ncol(inc) - Matrix::rankMatrix(inc)
  as.integer(nullity)
}

#' @title plot roll call graph
#' @description  plot roll call graph
#' @param roll_call a list of class roll call, assumed to be generated using input_roll_call
#' @export
plot.roll_call <- function(roll_call,
                           method = c("networkD3", "igraph"),
                           ...) {
  g <- igraph::graph_from_adjacency_matrix(roll_call$adj_mat)
  method <- method[1]
  if (method == "igraph") {
    plot(g, edge.arrow.size = 0)
  }
  if (method == "networkD3") {
    n <- length(roll_call$label)
    group <- integer(n)
    node_colour <- character(n)
    group[roll_call$voter] <- "voter"
    group[roll_call$y] <- "y"
    group[roll_call$n] <- "n"
    colour_scale <- htmlwidgets::JS("d3.scaleOrdinal().domain(['voter', 'y', 'n']).range(['blue', 'green', 'red']);")
    radius_calculation = htmlwidgets::JS(" Math.sqrt(d.nodesize)+20")
    gd3 <- networkD3::igraph_to_networkD3(g, group = group)
    networkD3::forceNetwork(Links = gd3$links, Nodes = gd3$nodes,
                 Source = 'source', Target = 'target',
                 NodeID = 'name', Group = 'group', opacity = 0.7,
                 fontSize = 20, colourScale = colour_scale,
                 bounded = TRUE, opacityNoHover = FALSE)
  }
}

#' @title predict roll call predict
#' @description  predict call graph
#' @param roll_call a list of class roll call, assumed to be generated using input_roll_call
#' @return data frame of predictions
#' @export
predict.roll_call <- function(roll_call) {
  g <- igraph::graph_from_adjacency_matrix(roll_call$adj_mat, "upper")
  dists <- igraph::distances(g, v = roll_call$voter, to = c(roll_call$y, roll_call$n))
  pred <- dists[ ,seq_along(roll_call$y)] / dists[ ,-seq_along(roll_call$y)] <= 1
  pred
}

#' @title accuracy roll call predict
#' @description  accuracy call graph
#' @param roll_call a list of class roll call, assumed to bea generated using input_roll_call
#' @return data frame of accuracy
#' @export
compute_accuracy <- function(roll_call) {
  pred <- predict(roll_call)
  correct <- pred == roll_call$votes
  sum(correct, na.rm = TRUE) / length(correct)
}

#' @title mutate
#' @description  mutate
#' @param roll_call a list of class roll call, assumed to bea generated using input_roll_call
#' @param p probability/frequency of mutation
#' @return mutated roll call
#' @export
# mutate_roll_call <- function(roll_call, p = 0.05) {
#   roll_call_copy <- roll_call
#   change_mat <- matrix(
#     sample(0:1, length(roll_call$adj_mat), c(1-p, p), replace = TRUE),
#     nrow = nrow(roll_call$adj_mat)
#   )
#   change_mat <-  (change_mat + t(change_mat)) %% 2
#   roll_call_copy$adj_mat <- (roll_call$adj_mat + change_mat) %% 2
#   roll_call_copy
# }
mutate_roll_call <- function(roll_call, rate = 1) {
  roll_call_copy <- roll_call
  n <- length(roll_call$label)
  n_mutations <- 1 + rgeom(1, prob = 1 / (1 + rate - 1))
  idx <- sample.int(n * (n - 1) / 2, n_mutations)
  x <- roll_call_copy$adj_mat[upper.tri(roll_call_copy$adj_mat)][idx]
  roll_call_copy$adj_mat[upper.tri(roll_call_copy$adj_mat)][idx] <- (x - 1) %% 2
  roll_call_copy
}

#' @title crossover_roll_call
#' @description  crossover_roll_call
#' @param roll_call a list of class roll call, assumed to bea generated using input_roll_call
#' @param p probability/frequency of mutation
#' @return crossover_roll_call
#' @export
crossover_roll_call <- function(roll_call1, roll_call2, type = c("join", "intersect", "pool"), pool_rate = 0.25) {
  new_roll_call <- roll_call1
  type <-  sample(type, 1)
  if (type == "join") {
    new_roll_call$adj_mat <- as.integer(roll_call1$adj_mat + roll_call2$adj_mat > 1)
  }
  if (type == "intersect") {
    new_roll_call$adj_mat <- roll_call1$adj_mat * roll_call2$adj_mat
  }
  if (type == "pool") {
    s <- rbinom(length(roll_call1$label), 1, pool_rate)
    new_roll_call$adj_mat[ ,as.logical(s)] <- roll_call2$adj_mat[  ,as.logical(s)]
  }
  new_roll_call
}

#' @title optimise graph with min homology
#' @description  optimise
#' @param roll_call a list of class roll call, assumed to bea generated using input_roll_call
#' @param population_size
#' @param n_generations
#' @param survival_rate
#' @param mutation_tate
#' @param max_time
#' @return optimized graph
#' @export
min_hom <- function(roll_call,
                    population_size = 250,
                    n_generations = 50,
                    survival_rate = 0.75,
                    mutation_rate = 0.05,
                    max_time = 300) {

}
init_popl <- 500
n_gen <- 100
survival <- 0.75
mutation <- 0.05