#' Find the variant (Latora) closeness centrality in a disconnected graph
#'
#' The diffusion degree of a node is defined as the cumulative contribution score of the node itself and its neighbors.
#' @details 
#' Diffusion degree \eqn{C_{DD}}{C(DD)} of node \eqn{v}{v} defined as:
#' \deqn{C_{DD}(v)=\lambda _{v} * C_{D}(v)+\sum_{i\in neighbors(v)}\lambda _{i} * C_{D}(i)}{C_DD(v)=lambda(v) * C_D(v)+sum(lambda(i) * C_D(i), i in neighbors(v))}
#' where \eqn{C_{D}}{C(DD)} is degree of of vertex and \eqn{\lambda}{lambda} is propagation probability of vertex. \cr
#' In a diffusion process, a node \eqn{v}{v} with propagation probability \eqn{\lambda _{v}}{lambda(v)}, can activate its neighbor \eqn{u}{u} with probability \eqn{\lambda _{v}}{lambda(v)}. \cr
#' When the diffusion process propagates to the next level, active neighbors of \eqn{v}{v} will try to activate their inactive neighbors. Thus the cumulative contribution in the diffusion process by neighbors of \eqn{v}{v} will be maximized when all of its neighbors will be activated in the previous step. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Diffusion_Degree}{Diffusion Degree}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constatnt, it specifies how to use the direction of the edges if a directed graph is analyzed. For 'out' only the outgoing edges are followed. For 'in' all vertices from which the source vertex is reachable in at most order steps are counted. 'all' ignores the direction of the edges. This argument is ignored for undirected graphs.
#' @param lambda Possibly a numeric vector giving propagation probability of vertices. The default is 1 for all vertices.
#' @param loops Logical; whether the loop edges are also counted.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Pal, Sankar K., Suman Kundu, and C. A. Murthy. "Centrality Measures, Upper Bound, and Influence Maximization in Large Scale Directed Social Networks." Fundamenta Informaticae 130.3 (2014): 317-342.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2))
#' diffusion.degree(g)
#' @export

diffusion.degree <- function (graph, vids = V(graph), mode = c("all", "out", "in"),
                              loops = TRUE, lambda = 1){
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  d <- degree(graph, mode=mode[1], loops=as.logical(loops)) * as.numeric(lambda)
  res <- numeric()
  for (v in V(graph)[vids]){
    n <- neighborhood(graph, 1, nodes=v, mode=mode[1])
    sm <- 0
    for (vv in n){
      sm <- sum(sm + d[vv])
    }
    res <- append(res, sm)
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }
  res
}