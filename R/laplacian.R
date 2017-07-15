#' Find the laplacian centrality
#'
#' The Laplacian centrality with respect to v is:
#' \deqn{C_{v}^{L}=(\Delta E)_{v}=d_{G}^{2}(v)+d_{G}(v)+2\sum_{v_{i}\in N(v)}d_{G}(v_{i})}{C(v,L)=Delta E(v)=d(G)^2(v)+d(G)(v)+2 * sum(d(G)(v(i)), v(i) in N(v))}
#' where G is a graph of \eqn{n}{n} vertices, \eqn{N(v)}{N(v)} is the set of neighbors of \eqn{v}{v} in G and \eqn{d_{G}(v_{i})}{d(G)(v(i))} is the degree of \eqn{v_{i}}{v(i)} in G. 
#' @details 
#' Laplacian centrality is a simple centrality measure that can be calculated in linear time. It is defined as the drop in the Laplacian energy (i.e. sum of squares of the eigenvalues in the Laplacian matrix) of the graph when the vertex is removed. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Laplacian_Centrality}{Laplacian Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constatnt, it specifies how to use the direction of the edges if a directed graph is analyzed. For 'out' only the outgoing edges are followed. For 'in' all vertices from which the source vertex is reachable in at most order steps are counted. 'all' ignores the direction of the edges. This argument is ignored for undirected graphs.
#' @param loops Logical; whether the loop edges are also counted.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Qi, Xingqin, et al. "Laplacian centrality: A new centrality measure for weighted networks." Information Sciences 194 (2012): 240-253.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2))
#' laplacian(g)
#' @export

laplacian <- function (graph, vids = V(graph), mode = c("all", "out", "in"), loops = TRUE){
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  res <- integer();
  for (aNode in V(graph)[vids]){
    aNeighbors <- neighborhood(graph, 1, nodes=aNode, mode=mode[1])[[1]][-1]
    deg <- 0
    for (nb in aNeighbors){
      deg <- deg + degree(graph, nb, mode=mode[1], loops=as.logical(loops))
    }
    aDeg <- degree(graph, aNode, mode=mode[1], loops=as.logical(loops))
    res <- append(res, aDeg^2 + aDeg + 2*deg )
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}