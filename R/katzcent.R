#' Find the Katz centrality (Katz Status Index)
#'
#' The Katz centrality for node i is:
#' \deqn{x_{i}=\alpha \sum_{j}A_{ij}x_{j}+\beta}{x(i)=alpha * sum(A(ij)*x(j), j) + beta}
#' where \eqn{A}{A} is the adjacency matrix of the graph G with eigenvalues \eqn{\lambda}{lambda}. The parameter \eqn{\beta}{beta} controls the initial centrality and \eqn{\alpha < \frac{1}{\lambda_{max}}}{alpha < 1/lambda(max)}. 
#' @details 
#' Katz centrality computes the relative influence of a node within a network by measuring the number of the immediate neighbors (first degree nodes) and also all other nodes in the network that connect to the node under consideration through these immediate neighbors. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Katz_Centrality}{Katz Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param alpha The alpha parameter, which must be between 0.0-0.2. The default is 0.1.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @author Algorithm adapted from CentiBin with thanks Dirk Koschutzki. (Junker, Bjorn H. 2006).
#' @references Newman, Mark. Networks: an introduction. Oxford University Press, 2010.
#' @references Junker, Bjorn H., Dirk Koschutzki, and Falk Schreiber. "Exploration of biological network centralities with CentiBiN." BMC bioinformatics 7.1 (2006): 219.
#' @examples
#' g <- barabasi.game(20)
#' katzcent(g)
#' @export

katzcent <- function (graph, vids = V(graph), alpha = 0.1){
  .cs.checkPreconditions(graph, c("loopFree"))
  vids <- .cs.as.igraph.vs(graph, vids)
  alpha <- as.double(alpha)
  if(alpha < 0.0 || alpha > 0.2) stop("Valid alpha is between 0.0-0.2", call. = FALSE)
  adjacencyMatrix <- as.matrix(get.adjacency(graph, names=FALSE))
  aMnrow <- nrow(adjacencyMatrix)
  # check alpha
  maxEigenvalue <- (1 / eigen(adjacencyMatrix)$values[1])
  if (alpha <= 0 || alpha >= maxEigenvalue) {
    stop("Invalid alpha value.", call. = FALSE)
  }
  res <- solve(diag(x = 1, nrow = aMnrow) - (alpha * t(adjacencyMatrix))) %*% matrix(1, nrow = aMnrow, ncol = 1)
  res <- res[,1]
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name
  }
  res[vids]
}