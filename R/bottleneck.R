#' Find the BottleNeck centrality score
#'
#' BottleNeck Centrality for vertex v defined as:
#' \deqn{BN(v) = \sum_{s\in v} P_{s}(v)}{BN(v) = sum(P(s)(v), s in v)}
#' Let \eqn{T_{s}}{T(s)} be a shortest path tree rooted at node \eqn{s}{s}.
#' \eqn{P_{s}(v) = 1}{P(s)(v) = 1} if more than \eqn{|V(T{s})|/4}{|V(T(s))|/4} paths from node \eqn{s}{s} to other nodes in \eqn{T_{s}}{T(s)} meet at the vertex \eqn{v}{v}, otherwise \eqn{P_{s}(v) = 0}{P(s)(v) = 0}. 
#' @details 
#' For each node \eqn{v}{v} in the graph, construct a tree \eqn{T_{v}}{T(v)} of shortest paths from that node to all other nodes in the graph. For a node \eqn{v}{v}, \eqn{n_{v}}{n(v)} is the number of nodes that are directly or indirectly connected to node \eqn{v}{v} (i.e. the tree \eqn{T_{v}}{T(v)} contains \eqn{n_{v}}{n(v)} nodes). So extract all nodes \eqn{w}{w} on the above defined tree \eqn{T_{v}}{T(v)} of shortest paths from node \eqn{v}{v}, such that more than \eqn{n_{v}/4}{n(v)/4} paths from \eqn{v}{v} to other nodes in the tree meet at node \eqn{w}{w}. Nodes \eqn{w}{w} extracted in this way represent 'bottle necks' of the shortest path tree \eqn{T_{v}}{T(v)} rooted at node \eqn{v}{v}, since at least \eqn{n_{v}/4}{n(v)/4} paths of the \eqn{n_{v}-node}{n(v)-node} tree \eqn{T_{v}}{T(v)} 'meet' at \eqn{w}{w}. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=BottleNeck}{BottleNeck}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constant, gives whether the shortest paths to or from the given vertices should be calculated for directed graphs. If out then the shortest paths from the vertex, if in then to it will be considered. If all, the default, then the corresponding undirected graph will be used, ie. not directed paths are searched. This argument is ignored for undirected graphs.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Przulj, N., Dennis A. Wigle, and Igor Jurisica. "Functional topology in a network of protein interactions." Bioinformatics 20.3 (2004): 340-348.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2))
#' bottleneck(g)
#' @export

bottleneck <- function (graph, vids = V(graph), mode = c("all", "out", "in")){
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  res <- matrix(0, vcount(graph), 1)
  for(v in V(graph)){
    s <- table(unlist(get.all.shortest.paths(graph, from=v, to=V(graph), mode=mode[1], weights=NA)$res))
    for(i in which(s > (length(s)/4), arr.ind = T)){
      if(i != v) res[i,1] <- res[i,1] + 1
    }
  }
  res <- res[,1]
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name
  }  
  res[vids]
}