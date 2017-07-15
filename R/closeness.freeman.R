#' Find the closeness centrality in a strongly connected graph
#'
#' Freeman closeness centrality defined as:
#' \deqn{\frac{1}{\sum_{i\neq v}d(v,i)}}{1/sum( d(v,i), i != v)}
#' @details 
#' Because closeness is infinite if there is no path between two vertex so freeman closeness require a strongly connected graph. In igraph if there is no (directed) path between vertex \eqn{v}{v} and \eqn{i}{i} then the total number of vertices is used in the formula instead of the path length. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Closeness_Centrality}{Closeness Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character string, defined the types of the paths used for measuring the distance in directed graphs. 'in' measures the paths to a vertex, 'out' measures paths from a vertex, all uses undirected paths. This argument is ignored for undirected graphs.
#' @param weights Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).
#' @param normalized Logical scalar, whether to calculate the normalized score.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @author Use igraph package closeness function.
#' @references Freeman, Linton C. "Centrality in social networks conceptual clarification." Social networks 1.3 (1979): 215-239.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2), directed=FALSE)
#' closeness.freeman(g)
#' @export

closeness.freeman <- function (graph, vids = V(graph), mode = c("all", "out", "in"), weights = NULL, normalized = FALSE) {
  .cs.checkPreconditions(graph, c("connected", "stronglyConnected"));
  vids <- .cs.as.igraph.vs(graph, vids)
  closeness(graph, vids=vids, mode=mode[1], weights=weights, normalized=normalized)
}