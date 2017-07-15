#' Find the markov centrality score
#'
#' The Markov centrality score uses the concept of a random walk through the graph to calculate the centrality of each vertex.
#' @details
#' The method uses the mean first-passage time from every vertex to every other vertex to produce a score for each vertex. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Markov_Centrality}{Markov Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the markov centrality values are returned.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @author Original code from Bioconductor SANTA package (Cornish AJ, 2014)
#' @references White, S. & Smyth, P. Algorithms for estimating relative importance in networks. Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining, 2003. ACM, 266-275.
#' @references Cornish AJ and Markowetz F (2014). "SANTA: Quantifying the Functional Content of Molecular Networks." PLOS Computational Biology, 10(9), pp. e1003808. http://dx.doi.org/10.1371/journal.pcbi.1003808.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2))
#' markovcent(g)
#' @export

markovcent <- function (graph, vids = V(graph)) {
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  edge.attr = NULL
  average.distances = FALSE
  c <- clusters(graph)
  clustN <- c$no
  clustMem <- c$membership
  vertN <- vcount(graph)
  M <- matrix(rep(Inf, vcount(graph)^2), vertN, vertN)
  adj <- get.adjacency(graph, attr = edge.attr, sparse = FALSE)
  max.dist <- max(adj)
  min.dist <- min(adj[adj > 0])
  for (i in 1:clustN) {
    locs <- which(clustMem == (i))
    M[locs, locs] <- suppressWarnings(.cs.markov.MFPTfct(delete.vertices(graph, 
        (1:vertN)[-locs]), edge.attr = edge.attr, max.dist = max.dist, min.dist = min.dist))
  }
  if (average.distances) {
    lt <- lower.tri(M)
    tmp <- apply(cbind(M[lt], t(M)[lt]), 1, mean)
    M[lt] <- tmp
    M <- t(M)
    M[lt] <- tmp
  }
  diag(M) <- 0
  M[M == Inf] <- 2 * max(M[M != Inf])
  M <- 1/apply(M, 2, mean)
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(M) <- V(graph)$name
  }  
  M[vids]
}

.cs.markov.MFPTfct <- function (g, edge.attr = NULL, max.dist = NULL, min.dist = NULL){
  adj <- get.adjacency(g, attr = edge.attr, sparse = FALSE)
  if (is.null(max.dist))
    max.dist <- max(adj)
  if (is.null(min.dist)) 
    min.dist <- min(adj[adj > 0])
  adj[adj > 0] <- max.dist - adj[adj > 0] + min.dist
  A <- adj/apply(adj, 1, sum)
  if (all(is.nan(A))) 
    A[1, 1] <- 1
  I <- diag(vcount(g))
  U <- matrix(1/nrow(A), nrow(A), nrow(A))
  pi <- U[1, ] %*% solve(I - A + U)
  e <- rep(1, nrow(A))
  Z <- solve(I - A - e %*% pi)
  Zdg <- I * Z
  E <- matrix(1, vcount(g), vcount(g))
  return ( (I - Z + E %*% Zdg) %*% (I * as.vector(1/pi)) )
}