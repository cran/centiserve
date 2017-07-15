#' Find the residual closeness centrality
#'
#' Residual closeness centrality defined as:
#' \deqn{C_{i}=\sum_{j\neq i}\frac{1}{2^{d(i,j)}}}{C_i=sum(1/2^d(i,j), j!=i)}
#' @details 
#' This function calculate closeness of a vertex as Dangalchev defination. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Residual_Closeness_Centrality}{Residual Closeness Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constant, gives whether the shortest paths to or from the given vertices should be calculated for directed graphs. If out then the shortest paths from the vertex, if in then to it will be considered. If all, the default, then the corresponding undirected graph will be used, ie. not directed paths are searched. This argument is ignored for undirected graphs.
#' @param weights Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Dangalchev, Chavdar. "Residual closeness in networks." Physica A: Statistical Mechanics and its Applications 365.2 (2006): 556-564.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2))
#' closeness.residual(g)
#' @export

closeness.residual <- function (graph, vids = V(graph), mode = c("all", "out", "in"), weights = NULL){
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  res <- double();
  sp <- shortest.paths(graph, v=V(graph), mode=mode[1], weights=weights)
  for(v in V(graph)[vids]){
    res <- append(res, sum(1/(2^sp[v,])))
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}