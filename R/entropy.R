#' Find the entropy centrality in a graph
#'
#' Entropy centrality measures centrality of nodes depending on their contribution to the entropy of the graph.
#' @details 
#' The centrality entropy measures \eqn{H_{ce}}{H(ce)} of a graph G, defined as:
#' \deqn{H_{ce}(G)=-\sum_{i=1}^{n}\gamma(v_{i})\times log_{2}\gamma(v_{i})}{H(ce)G=-sum(gamma(v(i)) X log2(gamma(v(i))), i=1,n)}
#' where  \eqn{\gamma(v_{i})=\frac{paths(v_{i})}{paths(v_{1}, v_{2}, ..., v_{M})}}{gamma(v(i))=paths(v(i))/paths(v(1), v(2), ..., v(M))} where \eqn{paths(v_{i})}{paths(v(i))} is the number of geodesic paths from node \eqn{v_{i}}{v(i)} to all the other nodes in the graph and \eqn{paths(v_{1}, v_{2}, ..., v_{M})}{paths(v(1), v(2), ..., v(M))} is the total number of geodesic paths M that exists across all the nodes in the graph. \cr
#' The centrality entropy provides information on the degree of centrality for a node in the graph. Those nodes that will split the graph in two or that will reduce substantially the number of paths available to reach other nodes when removed, will have a higher impact in decreasing the total centrality entropy of a graph. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Entropy_Centrality}{Entropy Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constant, gives whether the shortest paths to or from the given vertices should be calculated for directed graphs. If out then the shortest paths from the vertex, if in then to it will be considered. If all, the default, then the corresponding undirected graph will be used, ie. not directed paths are searched. This argument is ignored for undirected graphs.
#' @param weights Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Ortiz-Arroyo, Daniel, and DM Akbar Hussain. "An information theory approach to identify sets of key players." Intelligence and Security Informatics. Springer Berlin Heidelberg, 2008. 15-26.
#' @examples
#' g <- erdos.renyi.game(10, 1/10)
#' entropy(g)
#' @export

entropy <- function (graph, vids = V(graph), mode = c("all", "out", "in"), weights = NULL){
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  if (!is.null(weights) && any(!is.na(weights))) {
    E(graph)$weight <- as.numeric(weights)
  }else if(is.null(weights)) {
    E(graph)$weight <- 1
  }
  res <- double();
  for(v in V(graph)[vids]){
    g <- graph - V(graph)[v]
    sp <- shortest.paths(g, v=V(g), mode=mode[1])
    paths <- (length(sp[sp!=Inf]) - vcount(g))/2
    H <- 0.0
    if(paths > 0){
      for(vv in V(g)){
        Y <- (length(sp[vv, sp[vv,]!=Inf])-1)/paths
        if(Y > 0) H <- H + (Y * log2(Y))
      }   
    }
    res <- append(res, -H)
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}