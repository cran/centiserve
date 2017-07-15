#' Find the centroid value of graph vertices
#'
#' Centroid value C_{cen}(v) for node v defined as:
#' \deqn{C_{cen}(v) : = min{f(v, w) : w \in V{v}}}{Ccen(v) : = min(f(v, w) : w in V(v))}
#' where \eqn{f(v, w) : = \gamma_{v}(w) - \gamma_{w}(v)}{f(v, w) : = gamma(v) (w) - gamma(w) (v)}, and \eqn{\gamma_{v}(w)}{gamma(v) (w)} is the number of vertex closer to \eqn{v}{v} than to \eqn{w}{w}.
#' @details 
#' Centroid value computed by focusing the calculus on couples of nodes \eqn{(v,w)}{(v,w)} and systematically counting the nodes that are closer (in term of shortest path) to \eqn{v}{v} or to \eqn{w}{w}. The calculus proceeds by comparing the node distance from other nodes with the distance of all other nodes from the others, such that a high centroid value indicates that a node \eqn{v}{v} is much closer to other nodes. Thus, the centroid value provides a centrality index always weighted with the values of all other nodes in the graph. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Centroid_value}{Centroid value}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constant, gives whether the shortest paths to or from the given vertices should be calculated for directed graphs. If out then the shortest paths from the vertex, if in then to it will be considered. If all, the default, then the corresponding undirected graph will be used, ie. not directed paths are searched. This argument is ignored for undirected graphs.
#' @param weights Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @author Algorithm adapted from CentiLib (Grabler, Johannes, 2012).
#' @references Scardoni, Giovanni, Michele Petterlini, and Carlo Laudanna. "Analyzing biological network parameters with CentiScaPe." Bioinformatics 25.21 (2009): 2857-2859.
#' @references Grabler, Johannes, Dirk Koschutzki, and Falk Schreiber. "CentiLib: comprehensive analysis and exploration of network centralities." Bioinformatics 28.8 (2012): 1178-1179.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2), directed=FALSE)
#' centroid(g)
#' @export

centroid <- function (graph, vids = V(graph), mode = c("all", "out", "in"), weights = NULL){
  .cs.checkPreconditions(graph, c("stronglyConnected"))
  vids <- .cs.as.igraph.vs(graph, vids)
  numVertices <- vcount(graph)
  gammaValues <- matrix(0, numVertices, numVertices)
  for(u in V(graph)){
    for(v in V(graph)){
      for(w in V(graph)){
        if(shortest.paths(graph, v=u, to=w, mode=mode[1], weights=weights)[1,1] < shortest.paths(graph, v=v, to=w, mode=mode[1], weights=weights)[1,1]) gammaValues[u, v] <- gammaValues[u, v] + 1
      }
    }
  }
  fValues <- matrix(0, numVertices, numVertices)
  for (u in V(graph)) {
    for (v in V(graph)) {
      fValues[u, v] <- gammaValues[u, v] - gammaValues[v, u]
    }
  }
  res <- double()
  for (w in V(graph)[vids]) {
    result <- 1.7*(10^308)
    for (i in V(graph)) {
      if (fValues[w, i] < result && u != i) {
        result <- fValues[w, i]
      }
    }
    res <- append(res, result)
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}