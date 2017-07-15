#' Find the cross-clique connectivity (centrality)
#'
#' The cross-clique connectivity \eqn{X(v)}{X(v)} of a node is the number of cliques to which belongs. A node with a high \eqn{X(v)}{X(v)} value is called a highly cross-connected node.
#' @details 
#' Note: Directed graph considered as undirected ones and multiple edges and loops are ignored. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Cross-Clique_Connectivity}{Cross-Clique Connectivity}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Faghani, M., and U. Nguyen. "A Study of XSS Worm Propagation and Detection Mechanisms in Online Social Networks." (2013): 1-1.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2), directed=FALSE)
#' crossclique(g)
#' @export

crossclique <- function(graph, vids = V(graph)){
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  #Note: Directed graph considered as undirected ones and multiple edges and loops are ignored.
  res <- as.vector(table(unlist(cliques(graph))))
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name
  }  
  res[vids]
}