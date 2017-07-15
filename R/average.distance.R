#' Find the average distance of a node
#'
#' This function return average distance of a node in a strongly connected and loop free graph. 
#' @details 
#' Average distance of node \eqn{u}{u} to the rest of nodes in the net defined as:
#' \deqn{C_{u}=\frac{\sum_{w\in V}dis(u,w)}{n-1}}{C(u)=sum(dis(u,w)/(n-1), w in V)}
#' It is invers of closeness centrality. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Average_Distance}{Average Distance}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constant, gives whether the shortest paths to or from the given vertices should be calculated for directed graphs. If out then the shortest paths from the vertex, if in then to it will be considered. If all, the default, then the corresponding undirected graph will be used, ie. not directed paths are searched. This argument is ignored for undirected graphs.
#' @param weights Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references del Rio, Gabriel, Dirk Koschutzki, and Gerardo Coello. "How to identify essential genes from molecular networks?." BMC systems biology 3.1 (2009): 102.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2), directed=FALSE)
#' averagedis(g)
#' @import igraph
#' @export

averagedis <- function (graph, vids = V(graph), mode = c("all", "out", "in"), weights = NULL) {
  .cs.checkPreconditions(graph, c("stronglyConnected", "loopFree"))
  vids <- .cs.as.igraph.vs(graph, vids)
  res <- double()
  n <- vcount(graph) + 1
  sp <- shortest.paths(graph, v=V(graph), mode=mode[1], weights=weights)
  for (v in V(graph)[vids]) res <- append(res, sum(sp[v,])/n)
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}