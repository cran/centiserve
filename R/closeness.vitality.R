#' Find the closeness vitality centrality in a strongly connected graph
#'
#' Closeness vitality of a node is the change in the sum of distances between all node pairs when excluding that node.
#' @details 
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Closeness_Vitality}{Closeness Vitality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character string, defined the types of the paths used for measuring the distance in directed graphs. 'in' measures the paths to a vertex, 'out' measures paths from a vertex, all uses undirected paths. This argument is ignored for undirected graphs.
#' @param weights Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Brandes, U. & Erlebach, T. 2005. Network Analysis: Methodological Foundations, U.S. Government Printing Office.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2,1,4), directed=FALSE)
#' closeness.vitality(g)
#' @export

closeness.vitality <- function (graph, vids=V(graph), mode = c("all", "out", "in"), weights = NULL) {
  .cs.checkPreconditions(graph, c("connected", "stronglyConnected"))
  vids <- .cs.as.igraph.vs(graph, vids)
  res <- double()
  if (!is.null(weights) && any(!is.na(weights))) {
    E(graph)$weight <- as.numeric(weights)
  }else if(is.null(weights)) {
    E(graph)$weight <- 1
  }
  totalWI <- .cS.wienerIndex(graph, mode[1])
  for(v in V(graph)[vids]){
    subgraph <- delete.vertices(graph, v)
    if(!is.connected(subgraph, mode="strong")) stop("Subgraph of graph is not strongly connected.", call. = FALSE)
    res <- append(res, totalWI - .cS.wienerIndex(subgraph, mode[1]))
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}

# calculate wiener index for graph
# @param graph A igraph object
# @param mode
# @param weights
.cS.wienerIndex <- function (graph, mode) {
  return (sum(1/closeness(graph, vids=V(graph), mode=mode)));
}