#' Find the lobby index (centrality)
#'
#' The l-index or lobby index of a node \eqn{x}{x} is the largest integer \eqn{k}{k} such that \eqn{x}{x} has at least \eqn{k}{k} neighbors with a degree of at least \eqn{k}{k}. 
#' @details 
#' For more detail at \href{http://www.centiserver.org/?q1=centrality&q2=Lobby_Index}{Lobby Index}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constatnt, it specifies how to use the direction of the edges if a directed graph is analyzed. For 'out' only the outgoing edges are followed. For 'in' all vertices from which the source vertex is reachable in at most order steps are counted. 'all' ignores the direction of the edges. This argument is ignored for undirected graphs.
#' @param loops Logical; whether the loop edges are also counted.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Korn, A., A. Schubert, and A. Telcs. "Lobby index in networks." Physica A: Statistical Mechanics and its Applications 388.11 (2009): 2221-2226.
#' @examples
#' g <- random.graph.game(20, 3/10)
#' lobby(g)
#' @export

lobby <- function (graph, vids = V(graph), mode = c("all", "out", "in"), loops = TRUE){
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  d <- degree(graph, mode=mode[1], loops=as.logical(loops))
  res <- numeric()
  nd <- numeric()
  for (v in V(graph)[vids]){
    n <- neighborhood(graph, 1, nodes=v, mode=mode[1])
    for (vv in n){
      nd <- c(nd, d[vv])
    }
    res <- c(res, .cs.getHIndex(nd))
    nd <- numeric()
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res 
}

.cs.getHIndex <- function (m){
  mlen <- length(m)
  s <- numeric(mlen)
  for (i in 1:mlen){
    imin <- min(mlen, m[i])
    s[imin] <- s[imin] + 1
  }
  sum <- 0
  for (i in length(s):1){
    sum <- sum + s[i]
    if (sum >= i) return (i)
  }
  return (0)
}