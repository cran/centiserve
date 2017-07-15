#' Find the barycenter centrality score
#'
#' Barycenter scores are calculated as 1 / (total distance from vertex v to all other vertices) in a strongly connected network.
#' @details 
#' There are 2 types of distance centrality scores, Closeness Centrality and Barycenter Centrality. \cr
#' Barycenter Centrality for vertex \eqn{v}{v} defined as:
#' \deqn{1 / (total distance from v to all other vertices)}{1 / (total distance from v to all other vertices)}
#' Closeness scores are calculated using the formula \eqn{1 / (average distance from vertex v to all other vertices)}{1 / (average distance from vertex v to all other vertices)} and Barycenter scores are calculated as \eqn{1 / (total distance from vertex v to all other vertices)}{1 / (total distance from vertex v to all other vertices)}. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Barycenter_Centrality}{Barycenter Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices.
#' @param mode Character constant, gives whether the shortest paths to or from the given vertices should be calculated for directed graphs. If out then the shortest paths from the vertex, if in then to it will be considered. If all, the default, then the corresponding undirected graph will be used, ie. not directed paths are searched. This argument is ignored for undirected graphs.
#' @param weights Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Viswanath, Meghana. Ontology-based automatic text summarization. Diss. University of Georgia, 2009.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2), directed=FALSE)
#' barycenter(g)
#' @import igraph
#' @export

barycenter <- function (graph, vids = V(graph), mode = c("all", "out", "in"), weights = NULL){
  .cs.checkPreconditions(graph, c("connected", "stronglyConnected"))
  vids <- .cs.as.igraph.vs(graph, vids)
  res <- double()
  sp <- shortest.paths(graph, v=V(graph), mode=mode[1], weights=weights)
  for (v in V(graph)[vids]){
    res <- append(res, 1/sum(sp[v, sp[v,]!=Inf]))
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}