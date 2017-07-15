#' Find the topological coefficient of a node in a undirected graph
#'
#' The topological coefficient is a relative measure for the extent to which a node shares neighbors with other nodes.
#' @details 
#' Topological coefficient \eqn{T_{n}}{T(n)} of a node \eqn{n}{n} with \eqn{k_{n}}{k(n)} neighbors defined as:
#' \deqn{T_{n}=\frac{avg(J(n,m))}{k_{n}}}{T(n)=avg(J(n,m))/k(n)}
#' where \eqn{J(n,m)}{J(n.m)} is defined for all nodes \eqn{m}{m} that share at least one neighbor with \eqn{n}{n}. The value \eqn{J(n,m)}{J(n,m)} is the number of neighbors shared between the nodes \eqn{n}{n} and \eqn{m}{m}, plus one if there is a direct link between \eqn{n}{n} and \eqn{m}{m}. \cr
#' Nodes that have one or no neighbors are assigned a topological coefficient of zero. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Topological_Coefficient}{Topological Coefficient}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Assenov, Yassen, et al. "Computing topological parameters of biological networks." Bioinformatics 24.2 (2008): 282-284.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2), directed=FALSE)
#' topocoefficient(g)
#' @export

topocoefficient <- function (graph, vids = V(graph)){
  .cs.checkPreconditions(graph, c("undirected"))
  vids <- .cs.as.igraph.vs(graph, vids)
  tcs <- double()
  for (aNode in V(graph)[vids]){
    aNeighbors <- neighborhood(graph, 1, aNode)[[1]][-1]
    comNeNodes <- vector()
    tc <- 0
    for (nb in aNeighbors){
      currentComNeNodes <- neighborhood(graph, 1, nb)[[1]][-1]
      for (n in currentComNeNodes){
        if(n != aNode){
          tc <- tc + 1
          if(!(n %in% comNeNodes)){
            comNeNodes <- append(comNeNodes, n)
            if ( n %in% aNeighbors ){
              tc <- tc + 1
            }
          }
        }
      }
    }
    tcs <- c(tcs, tc / (length(comNeNodes) * length(aNeighbors)) )
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(tcs) <- V(graph)$name[vids]
  }  
  tcs
}