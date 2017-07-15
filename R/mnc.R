#' Find the maximum neighborhood component (MNC)
#'
#' Maximum Neighborhood Component defined as:
#' \deqn{MNC(v)=\left|V(MC(v))\right|}{MNC(v)=|V(MC(v))|}
#' where where MC(v) is a maximum connected component of the \eqn{G[N(v)]}{G[N(v)]} and \eqn{G[N(v)]}{G[N(v)]} is the induced subgraph of G by \eqn{N(v)}{N(v)} and \eqn{N(v)}{N(v)} is neighborhoods of node \eqn{v}{v}.
#' @details 
#' The neighborhood of a node \eqn{v}{v}, nodes adjacent to \eqn{v}{v}, induce a subnetwork \eqn{N(v)}{N(v)}. The score of node \eqn{v}{v}, \eqn{MNC(v)}{MNC(v)}, is defined to be the size of the maximum connected component of \eqn{N(v)}{N(v)}. The neighborhood \eqn{N(v)}{N(v)} is the set of nodes adjacent to \eqn{v}{v} and does not contain node \eqn{v}{v}. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=MNC_Maximum_Neighborhood_Component}{MNC-Maximum Neighborhood Component}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constatnt, it specifies how to use the direction of the edges if a directed graph is analyzed. For 'out' only the outgoing edges are followed. For 'in' all vertices from which the source vertex is reachable in at most order steps are counted. 'all' ignores the direction of the edges. This argument is ignored for undirected graphs.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Lin, Chung-Yen, et al. "Hubba: hub objects analyzer-a framework of interactome hubs identification for network biology." Nucleic acids research 36.suppl 2 (2008): W438-W443.
#' @examples
#' g <- random.graph.game(20, 3/10)
#' mnc(g)
#' @export

mnc <- function (graph, vids = V(graph), mode = c("all", "out", "in")){
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  res <- integer();
  for(v in V(graph)[vids]){
    g <- induced.subgraph(graph, neighborhood(graph, 1, nodes=v, mode=mode[1])[[1]][-1])
    csize <- clusters(g,  mode="strong")$csize
    if(length(csize) == 0){
      res <- append(res, 0)
    }else{
      res <- append(res, max(csize))
    }
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}