#' Find the variant (Latora) closeness centrality in a disconnected graph
#'
#' Variant (Latora) closeness centrality defined as:
#' \deqn{\sum_{i \neq v}\frac{1}{d(v,i)}}{sum(1/d(v,i), i != v)}
#' @details 
#' This variant (sum of inversed distances to all other nodes instead of the inversed of the sum of distances to all other nodes) applicable to both connected and unconnected graphs. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Closeness_Centrality}{Closeness Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constant, gives whether the shortest paths to or from the given vertices should be calculated for directed graphs. If out then the shortest paths from the vertex, if in then to it will be considered. If all, the default, then the corresponding undirected graph will be used, ie. not directed paths are searched. This argument is ignored for undirected graphs.
#' @param weights Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).
#' @param normalized Logical scalar, whether to calculate the normalized score.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Latora V., Marchiori M., Efficient behavior of small-world networks, Physical Review Letters, V. 87, p. 19, 2001.
#' @references Opsahl, Tore, Filip Agneessens, and John Skvoretz. "Node centrality in weighted networks: Generalizing degree and shortest paths." Social Networks 32.3 (2010): 245-251.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2))
#' closeness.latora(g)
#' @export

closeness.latora <- function (graph, vids = V(graph), mode = c("all", "out", "in"), weights = NULL, normalized = FALSE) {
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  res <- double()
  sp <- shortest.paths(graph, v=V(graph), mode=mode[1], weights=weights)
  for (v in V(graph)[vids]) res <- append(res, sum(1/sp[v, sp[v,]!=0]))
  if(as.logical(normalized)) res <- res/(vcount(graph)*(vcount(graph)-1))
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}