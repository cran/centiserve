#' Find the communicability betweenness centrality
#'
#' The communicability betweenness of a node r is:
#' \deqn{\omega_{r} = \frac{1}{C} \sum_{p}\sum_{q}\frac{G_{prq}}{G_{pq}}, p\neq q,p\neq r, q\neq r}{omega(r) = 1/C * sum(sum(G(prq)/G(pq), q), p), p!=q,p!=r,q!=r}
#' where where \eqn{G_{prq} = (e^{A})_{pq} - (e^{A+E(r)})_{pq}}{G(prq) = e^A(pq) - e^(A+E(r))(pq)} is the number of walks involving node \eqn{r}{r}, \eqn{G_{pq} = (e^{A})_{pq}}{G(pq) = e^A(pq)} is the number of closed walks starting at node p and ending at node \eqn{q}{q}, and \eqn{C = (n-1)^{2}-(n-1)}{C = (n-1)^2 - (n-1)} is a normalization factor equal to the number of terms in the sum.
#' @details 
#' Communicability betweenness measure makes use of the number of walks connecting every pair of nodes as the basis of a betweenness centrality measure. \cr
#' The resulting \eqn{\omega_{r}}{omega(r)} takes values between zero and one. The lower bound cannot be attained for a connected graph, and the upper bound is attained in the star graph. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Communicability_Betweenness_Centrality}{Communicability Betweenness Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param normalized Logical scalar, whether to calculate the normalized score.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @author Algorithm adapted from NetworkX 1.9 (Hagberg, A. 2008).
#' @references Estrada, Ernesto, Desmond J. Higham, and Naomichi Hatano. "Communicability betweenness in complex networks." Physica A: Statistical Mechanics and its Applications 388.5 (2009): 764-774.
#' @references Hagberg, Aric, Pieter Swart, and Daniel S Chult. Exploring network structure, dynamics, and function using NetworkX. No. LA-UR-08-05495; LA-UR-08-5495. Los Alamos National Laboratory (LANL), 2008.
#' @examples
#' g <- graph(c(1,2,2,3,2,6,6,5,3,5,3,4,5,4,4,7), directed=FALSE)
#' communibet(g)
#' @export

communibet <- function (graph, vids = V(graph), normalized = FALSE){
  if(suppressWarnings(requireNamespace("expm", logical.return=T, warn.conflicts=F, quietly=T, verbose=F)!=TRUE)){
    stop("The 'expm' package needed for this function to work. Please install it.", call. = FALSE)
  }
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  n <- vcount(graph)
  A <- as.matrix(get.adjacency(graph, names=FALSE))
  expA = expm::expm(A)
  res <- double()
  for(v in V(graph)[vids]){
    # remove row and col of node v
    comrow <- A[v,]
    comcol <- A[,v]
    A[v,] <- 0
    A[,v] <- 0
    B <- (expA - expm::expm(A)) / expA
    # sum with row/col of node v and diag set to zero
    B[v,] <- 0
    B[,v] <- 0
    B <- B - diag(diag(B))
    res <- append(res, sum(B))
    # put row and col back
    A[v,] <- comrow
    A[,v] <- comcol
  }
  if (as.logical(normalized) && n > 2){
    res <- res * (1.0/((n-1.0)**2-(n-1.0)))
  }	
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}