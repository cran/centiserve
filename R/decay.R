#' Find the decay centrality of a given vertex
#'
#' Decay centrality of a given vertex \eqn{x}{x} of a graph G is define as:
#' \deqn{\sum_{y \in V(G)}\sigma ^{d(x,y)}}{sum(sigma ^ d(x,y), y in V(G))}
#' where \eqn{d(x,y)}{d(x,y)} denotes the distance between \eqn{x}{x} and \eqn{y}{y} and \eqn{\sigma \in (0,1)}{sigma in (0,1)} is a parameter.
#' @details 
#' Decay centrality is a centrality measure based on the proximity between a choosen vertex and every other vertex weighted by the decay. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Decay_Centrality}{Decay Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constant, gives whether the shortest paths to or from the given vertices should be calculated for directed graphs. If out then the shortest paths from the vertex, if in then to it will be considered. If all, the default, then the corresponding undirected graph will be used, ie. not directed paths are searched. This argument is ignored for undirected graphs.
#' @param weights Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).
#' @param decay A decay parameter which the default is 0.5.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Jana Hurajova, Silvia Gago and Tomas Madaras, Decay Centrality, 15th Conference of Kosice Mathematicians. Herl'ny 2.-5. aprila 2014.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2), directed=FALSE)
#' decay(g)
#' @export

decay <- function (graph, vids = V(graph), mode = c("all", "out", "in"), weights = NULL, decay = 0.5){
  .cs.checkPreconditions(graph, c("stronglyConnected"))
  vids <- .cs.as.igraph.vs(graph, vids)
  decay <- as.double(decay)
  if(decay <= 0 || decay >= 1) stop("The decay parameter must be between 0 and 1", call. = FALSE)
  sp <- shortest.paths(graph, mode=mode[1], weights=weights)
  res <- rowSums(decay ^ sp[,])
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name
  }  
  res[vids]
}