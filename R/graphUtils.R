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

#' @title initial incidence matrix
#' @description receives roll call data and creates an incidence matrix
#' @details the incidence matrix \code{E} has entry \code{e_{ij}} equal to
#' \itemize{
#'     \code{1} if an edge starts in node \code{i} and finishes in \code{j}
#'     \code{-1} if an edge ends in node \code{i} and finishes in \code{j}
#'     \code{0} if nodes \code{i} and \code{j} are not connected
#' }
#' this matrix creates an incidence matrix from roll call data that has perfect prediction
#' but usually has the highest possible dimension.
#' @examples
#' data(onedim_voting)
#' create_incidence_matrix(onedim_voting)
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
  li <- c(nodes, list(adj_mat = adj_mat))
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
}