#' Find the pairwise disconnectivity index 
#'
#' The pairwise disconnectivity index of vertex \eqn{v}{v}, \eqn{Dis(v)}{Dis(v)} defined as:
#' \deqn{Dis(v)=\frac{N_{0}-N_{-v}}{N_{0}}=1-\frac{N_{-v}}{N_{0}}}{Dis(v)=(N(0)-N(-v))/N(0)=1-(N(-v)/N(0))}
#' where \eqn{N_{0}}{N(0)} is the total number of ordered pairs of vertices in a network that are connected by at least one directed path of any length. It is supposed that \eqn{N_{0}}{N(0)} > 0, i.e., there exists at least one edge in the network that links two different vertices. \eqn{N_{-v}}{N(-v)} is the number of ordered pairs that are still connected after removing vertex \eqn{v}{v} from the network, via alternative paths through other vertices.
#' @details 
#' The pairwise disconnectivity defined as index of vertex \eqn{v}{v}, \eqn{Dis(v)}{Dis(v)}, as the fraction of those initially connected pairs of vertices in a network which become disconnected if vertex \eqn{v}{v} is removed from the network. The pairwise disconnectivity index quantifies how crucial an individual element is for sustaining the communication ability between connected pairs of vertices in a network that is displayed as a directed graph. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Pairwise_Disconnectivity_Index}{Pairwise Disconnectivity Index}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Potapov, Anatolij P., Bjorn Goemann, and Edgar Wingender. "The pairwise disconnectivity index as a new metric for the topological analysis of regulatory networks." BMC bioinformatics 9.1 (2008): 227.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2))
#' pairwisedis(g)
#' @export

pairwisedis <- function (graph, vids = V(graph)){
  .cs.checkPreconditions(graph, c("directed"))
  vids <- .cs.as.igraph.vs(graph, vids)
  res <- double()
  sp <- shortest.paths(graph, v=V(graph), mode="out", weights=NA)
  allPaths <- length(sp[sp!=Inf]) - vcount(graph)
  for(v in V(graph)[vids]){
    g <- delete.vertices(graph, v)
    sp <- shortest.paths(g, v=V(g), mode="out", weights=NA)
    paths <- length(sp[sp!=Inf]) - vcount(g)
    res <- append(res, (allPaths-paths)/allPaths)
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}