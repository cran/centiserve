#' Find the radiality centrality in a graph
#'
#' The radiality is a node centrality index and will give high centralities to vertices that are a short distance to every other vertex in its reachable neighborhood compared to its diameter. \cr
#' @details 
#' Radiality centrality defined as:
#' \deqn{C_{rad}(v)=\frac{\sum_{w\in V}(d+1-d(v,w))}{n-1}}
#' where \eqn{d}{d} is diameter of graph G with \eqn{n}{n} vertices and \eqn{d(v,w)}{d(v,w)} is distance between vertex \eqn{v}{v} and \eqn{w}{w}. \cr
#' The radiality of a node \eqn{v}{v} is calculated by computing the shortest path between the node \eqn{v}{v} and all other nodes in the graph. The value of each path is then subtracted by the value of the diameter +1 (G+1) and the resulting values are summated. Finally, the obtained value is divided for the number of nodes -1 (n-1). The radiality should be always compared to the closeness and to the eccentricity: a node with high eccentricity + high closeness+ high radiality is a consistent indication of a high central position in the graph. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Radiality_Centrality}{Radiality Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constant, gives whether the shortest paths to or from the given vertices should be calculated for directed graphs. If out then the shortest paths from the vertex, if in then to it will be considered. If all, the default, then the corresponding undirected graph will be used, ie. not directed paths are searched. This argument is ignored for undirected graphs.
#' @param weights Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Wolfram Research, Inc., Mathematica, Version 10.0, Champaign, IL (2014). http://reference.wolfram.com/language/ref/RadialityCentrality.html
#' @references Scardoni, G., Laudanna, C., Tosadori, G., Fabbri, F. & Faizaan, M. CentiScaPe: Network centralities for Cytoscape. http://www.cbmc.it/~scardonig/centiscape/CentiScaPefiles/CentralitiesTutorial.pdf
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2))
#' radiality(g)
#' @export

radiality <- function (graph, vids = V(graph), mode = c("all", "out", "in"), weights = NULL){
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  res <- double()
  numVertices <- vcount(graph)
  diam <- diameter(graph)
  sp <- shortest.paths(graph, mode=mode[1], weights=weights)
  for(v in V(graph)[vids]){
    rad <- 0.0
    for(vv in V(graph)){
      if (sp[v, vv] == Inf) {
        rad <- Inf
        break
      }else{
        rad <- rad + (diam + 1.0 - sp[v, vv])
      }
    }
    res <- append(res, rad / (numVertices - 1.0))
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}