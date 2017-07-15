#' Find the lin centrality in a graph
#'
#' Lin centrality of a node \eqn{x}{x} with a nonempty coreachable set is:
#' \deqn{\frac{\left|\left\{y|d(x,y)<\infty \right\}\right|^2}{\sum_{d(x,y)<\infty}d(x,y)}}{|{y|d(x,y)<infty}|^2/sum(d(x,y), d(x,y)<infty)}
#' where 
#' @details 
#' Lin centrality consider closeness not the inverse of a sum of distances, but rather the inverse of the average distance, which entails a first multiplication by the number of coreachable nodes. This change normalizes closeness across the graph. Now, however, we want nodes with a larger coreachable set to be more important, given that the average distance is the same, so we multiply again by the number of coreachable nodes. Nodes with an empty coreachable set have centrality 1 by definition. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Lin_Centrality}{Lin Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constant, gives whether the shortest paths to or from the given vertices should be calculated for directed graphs. If out then the shortest paths from the vertex, if in then to it will be considered. If all, the default, then the corresponding undirected graph will be used, ie. not directed paths are searched. This argument is ignored for undirected graphs.
#' @param weights Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Lin, Nan. Foundations of social research. New York: McGraw-Hill, 1976.
#' @references Boldi, Paolo, and Sebastiano Vigna. "Axioms for centrality." Internet Mathematics just-accepted (2014): 00-00.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2))
#' lincent(g)
#' @export

lincent <- function (graph, vids = V(graph), mode = c("all", "out", "in"), weights = NULL){
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  res <- double()
  sp <- shortest.paths(graph, v=V(graph), mode=mode[1], weights=weights)
  for (v in V(graph)[vids]){
    linrow <- sp[v, sp[v,]!=Inf]
    res <- append(res, ((length(linrow)-1)^2) / sum(linrow))
  }
  res[res==Inf] <- 1
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}