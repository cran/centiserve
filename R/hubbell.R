#' Find the Hubbell centrality or the Hubbell Index
#'
#' Hubbell centrality defined as:
#' \deqn{C_{h} = E + WC_{h}}{C(h) = E + WC(h)}
#' where \eqn{E}{E} is some exogeneous input and \eqn{W}{w} is a weight matrix derived from the adjancancy matrix \eqn{A}{A}.
#' @details 
#' This centrality value is defined by means of a weighted and loop allowed network. The weighted adjacency matrix \eqn{W}{w} of a network G is asymmetric and contains real-valued weights for each edge. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Hubbell_Index}{Hubbell Index}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param weights Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).
#' @param weightfactor The weight factorLogical which must be greater than 0. The defualt is 0.5.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @author Algorithm adapted from CentiLib (Grabler, Johannes, 2012).
#' @references Hubbell, Charles H. "An input-output approach to clique identification." Sociometry (1965): 377-399.
#' @references Grabler, Johannes, Dirk Koschutzki, and Falk Schreiber. "CentiLib: comprehensive analysis and exploration of network centralities." Bioinformatics 28.8 (2012): 1178-1179.
#' @examples
#' g <- barabasi.game(100)
#' hubbell(g)
#' @export

# Original code from CentiLib

hubbell <- function (graph, vids = V(graph), weights = NULL, weightfactor = 0.5){
  .cs.checkPreconditions(graph, c("connected", "loopFree"))
  vids <- .cs.as.igraph.vs(graph, vids)
  if (!is.null(weights) && any(!is.na(weights))) {
    E(graph)$weight <- as.numeric(weights)
  }else if(is.null(weights)) {
    E(graph)$weight <- 1
  }
  vlen <- vcount(graph)
  weightFactor <- as.double(weightfactor)
  if(is.na(weightFactor)) stop("Hubbell index centrality is not solvable for this graph. Weight factor must be >0", call. = FALSE)
  if(weightFactor <= 0) stop("Hubbell index centrality is not solvable for this graph. Weight factor must be >0", call. = FALSE)
  weightedAdjacencyMatrix <- as.matrix(get.adjacency(graph, names=FALSE, attr="weight")) * weightFactor
  ev <- eigen(weightedAdjacencyMatrix)$values
  if(length(ev[ev > 1.0]) > 0){
    stop("Hubbell index centrality is not solvable for this graph.", call. = FALSE)
  }
  # res = ((I-W)^-1)*vectorE
  res <- solve(diag(x = 1, nrow = vlen) - weightedAdjacencyMatrix) %*% matrix(1, nrow = vlen, ncol = 1)       
  res <- res[,1]
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name
  }  
  res[vids]
}