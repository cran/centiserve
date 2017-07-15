#' Find the geodesic k-path centrality
#'
#' Geodesic K-path centrality counts neighbours as those that are on a geodesic path less than "k" away. 
#' @details 
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Geodesic_K-Path_Centrality}{Geodesic K-Path Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constant, gives whether the shortest paths to or from the given vertices should be calculated for directed graphs. If out then the shortest paths from the vertex, if in then to it will be considered. If all, the default, then the corresponding undirected graph will be used, ie. not directed paths are searched. This argument is ignored for undirected graphs.
#' @param weights Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).
#' @param k The k parameter. The default is 3.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Borgatti, Stephen P., and Martin G. Everett. "A graph-theoretic perspective on centrality." Social networks 28.4 (2006): 466-484.
#' @examples
#' g <- barabasi.game(100)
#' geokpath(g)
#' @export

geokpath <- function (graph, vids = V(graph), mode = c("all", "out", "in"), weights = NULL, k = 3){
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  k <- as.integer(k)
  if(k <= 0) stop("The k parameter must be greater than 0.", call. = FALSE)
  res <- integer()
  sp <- shortest.paths(graph, v=V(graph), mode=mode[1], weights=weights)
  for (v in V(graph)[vids]){
    res <- append(res, length(sp[v, sp[v,] <= k]) - 1);
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}