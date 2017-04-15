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
    adj_mat[nodes$y[votes_yes], i] <- 1
    adj_mat[i, nodes$n[!votes_yes]] <- 1
    adj_mat[nodes$n[!votes_yes], i] <- 1
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
plot.roll_call <- function(roll_call) {
  g <- igraph::graph_from_adjacency_matrix(roll_call$adj_mat)
  plot(g, edge.arrow.size = 0)
}

#' @title predict roll call predict
#' @description  predict call graph
#' @param roll_call a list of class roll call, assumed to be generated using input_roll_call
#' @return data frame of predictions
#' @export
predict.roll_call <- function(roll_call) {
  g <- igraph::graph_from_adjacency_matrix(roll_call$adj_mat)
  dists <- igraph::distances(g)
  pred <- data.frame(
    row.names = roll_call$label[roll_call$voter]
  )
  for (j in seq_along(roll_call$y)) {
    new_col <- apply(dists[roll_call$voter, c(roll_call$y[j], roll_call$n[j])], 1, which.min)
    new_col <- as.integer(new_col == 1)
    pred[names(roll_call$votes)[j]] <- new_col
  }
  pred
}

#' @title accuracy roll call predict
#' @description  accuracy call graph
#' @param roll_call a list of class roll call, assumed to bea generated using input_roll_call
#' @return data frame of accuracy
#' @export
accuracy <- function(roll_call) {
  pred <- predict(roll_call)
  correct <- pred == roll_call$votes
  sum(correct) / length(correct)
}

#' @title mutate
#' @description  mutate
#' @param roll_call a list of class roll call, assumed to bea generated using input_roll_call
#' @param p probability/frequency of mutation
#' @return mutated roll call
#' @export
mutate_roll_call <- function(roll_call, p = 0.05) {
  roll_call_copy <- roll_call
  change_mat <- matrix(
    sample(0:1, length(roll_call$adj_mat), c(1-p, p), replace = TRUE),
    nrow = nrow(roll_call$adj_mat)
  )
  change_mat <-  change_mat + t(change_mat)
  roll_call_copy$adj_mat <- (roll_call$adj_mat + change_mat) %% 2
  roll_call_copy
}

#' @title combine_roll_call
#' @description  combine_roll_call
#' @param roll_call a list of class roll call, assumed to bea generated using input_roll_call
#' @param p probability/frequency of mutation
#' @return combine_roll_call
#' @export
combine_roll_call <- function(roll_call1, roll_call2) {
  new_roll_call <- roll_call1
  new_roll_call$adj_mat <- roll_call1$adj_mat * roll_call2$adj_mat
  new_roll_call
}