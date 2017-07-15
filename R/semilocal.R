#' Find the semi local centrality (or local centrality)
#'
#' The local centrality \eqn{CL(v)}{CL(v)} of node \eqn{v}{v} is defined as:
#' \deqn{C_{L}(v)=\sum_{u\in \Gamma (v)}Q(u)}{C_L(v)=sum(Q(u), u in Gamma(v))}
#' where
#' \deqn{Q(u)=\sum_{w\in \Gamma (u)}N(w)}{Q(u)=sum(N(w), w in Gamma(u))}
#' and \eqn{\Gamma (u)}{Gamma (u)} is the set of the nearest neighbors of node \eqn{u}{u} and \eqn{N(w)}{N(w)} is the number of the nearest and the next nearest neighbors of node \eqn{w}{w}.
#' @details 
#' The local centrality is proposed aiming at identifying the influencers in undirected network, it can be applied to directed network as well with a modified definition of \eqn{N(w)}{N(w)}. Of course, for directed network, \eqn{N(w)}{N(w)} should be the number of the nearest and next nearest upstream nodes of node \eqn{w}{w}. \cr
#' Local centrality measure is likely to be more effective to identify influential nodes than degree centrality measure as it utilizes more information, while it has much lower computational complexity than the betweenness and closeness centralities. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Semi_Local_Centrality}{Semi_Local Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constatnt, it specifies how to use the direction of the edges if a directed graph is analyzed. For 'out' only the outgoing edges are followed. For 'in' all vertices from which the source vertex is reachable in at most order steps are counted. 'all' ignores the direction of the edges. This argument is ignored for undirected graphs.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Chen, Duanbing, et al. "Identifying influential nodes in complex networks." Physica a: Statistical mechanics and its applications 391.4 (2012): 1777-1787.
#' @examples
#' g <- barabasi.game(10)
#' semilocal(g)
#' @export

semilocal <- function (graph, vids = V(graph), mode = c("all", "out", "in")){
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  res <- integer()
  for (v in V(graph)[vids]){
    vNeighbors <- neighborhood(graph, 1, nodes=v, mode[1])[[1]][-1]
    sl <- 0
    for (vv in vNeighbors){
      vvNeighbors <- neighborhood(graph, 1, nodes=vv, mode[1])[[1]][-1]
      for (vvv in vvNeighbors){
        sl <- sl + length(neighborhood(graph, 2, nodes=vvv, mode[1])[[1]][-1])
      }
    }
    res <- append(res, sl)
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}