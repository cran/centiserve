#' Find current-flow closeness centrality
#'
#' Current-flow closeness centrality is defined by:
#' \deqn{C_{cc}(s)=\frac{n-1}{\sum_{s\neq t}p_{st}(s)-p_{st}(t)} for all: s \in V}{C_cc(s)=(n-1)/sum(p_st(s), s!=t) - p_st(t) for all: s in V}
#' where \eqn{(n-1)}{(n-1)} is a normalizing factor, \eqn{p_{st}(s)}{p_st(s)} is the absolute electrical potential of vertex s based on the electrical current supply from vertex \eqn{s}{s} to vertex \eqn{t}{t}, and \eqn{p_{st}(s) - p_{st}(t)}{p_st(s) - p_st(t)} corresponds to the effective resistance typically measured as voltage, which can be interpreted as an alternative measure of distance between \eqn{s}{s} and \eqn{t}{t}.
#' @details 
#' The closeness index based on shortest paths can also be transformed to a measure based on electrical current. For the electrical current model set, Brandes et al. developed an alternative measure of the distance between two vertices \eqn{s}{s} and \eqn{t}{t}, which is defined as the difference of their electrical potentials. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Current-Flow_Closeness_Centrality}{Current-Flow Closeness Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param weights Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @author Note: This implementation is based on Daniel Fleischer's implementation for yFiles and JMP which convert to java implementation for CentiLib by Johannes Graessler and Dirk Koschuetzki.
#' @references Brandes, Ulrik, and Daniel Fleischer. Centrality measures based on current flow. Springer Berlin Heidelberg, 2005.
#' @references Grabler, Johannes, Dirk Koschutzki, and Falk Schreiber. "CentiLib: comprehensive analysis and exploration of network centralities." Bioinformatics 28.8 (2012): 1178-1179.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2), directed=FALSE)
#' closeness.currentflow(g)
#' @export

closeness.currentflow <- function (graph, vids = V(graph), weights = NULL) {
  .cs.checkPreconditions(graph, c("undirected", "connected", "loopFree"))
  vids <- .cs.as.igraph.vs(graph, vids)
  if (!is.null(weights) && any(!is.na(weights))) {
    E(graph)$weight <- as.numeric(weights)
  }else if(is.null(weights)) {
    E(graph)$weight <- 1
  }
  numVertices <- vcount(graph)
  resultVector <- double(numVertices)
  solveTemporaryVector <- matrix(nrow=numVertices - 1, ncol=1)
  kroneckerVector <- matrix(nrow=numVertices - 1, ncol=1)
  # L = D - A
  L <- diag(degree(graph, mode="all", loops=FALSE), nrow=numVertices, ncol=numVertices) - get.adjacency(graph, attr="weight")
  L <- L[-1,-1];
  for (i in 1:(numVertices - 1)) {
    solveTemporaryVector[,1] <- 0
    kroneckerVector[,1] <- 0
    kroneckerVector[i,1] <- 1
    solveTemporaryVector <- solve(L, kroneckerVector)
    resultVector[1] <- resultVector[1] + solveTemporaryVector[i, 1]
    resultVector[i + 1] <-  resultVector[i + 1] + solveTemporaryVector[i, 1]
    for (j in 1:(numVertices - 1)) {
      resultVector[i + 1] <- resultVector[i + 1] + (solveTemporaryVector[i, 1] - 2 * solveTemporaryVector[j, 1])
      resultVector[j + 1] <- resultVector[j + 1] + solveTemporaryVector[i, 1]
    }
  }
  resultVector <- (numVertices - 1) / resultVector
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(resultVector) <- V(graph)$name
  }  
  resultVector[vids]
}