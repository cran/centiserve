#' Find the SALSA as 'hub' or 'authority' score
#'
#' The Stochastic Approach for Link-Structure Analysis (SALSA) is combination of HITS and PageRank which creates a neighborhood graph using authority and hub pages and links and create a bipartite graph of the authority and hub pages in the neighborhood graph.
#' @details 
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=SALSA}{SALSA}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param score Character constant, gives which score should be calculated and must be one of 'hub' or 'authority'. The default is 'hub'.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Lempel, Ronny, and Shlomo Moran. "SALSA: the stochastic approach for link-structure analysis." ACM Transactions on Information Systems (TOIS) 19.2 (2001): 131-160.
#' @examples
#' g <- barabasi.game(10)
#' salsa(g)
#' @export

salsa <- function (graph, vids = V(graph), score = c("hub", "authority")){
  .cs.checkPreconditions(as.undirected(graph), c("stronglyConnected"))
  vids <- .cs.as.igraph.vs(graph, vids)
  if(missing(score)){
    score <- "hub"
  }else{
    if(length(score) > 1 || (as.character(score)!="hub" && as.character(score)!="authority")){
      stop("The score should be one of 'hub' or 'authority'.", call. = FALSE)
    }  
  }
  L <- as.matrix(get.adjacency(graph, names=FALSE))
  Lr <- prop.table(L, 1)
  Lr[is.nan(Lr)] <- 0
  Lc <- prop.table(L, 2)
  Lc[is.nan(Lc)] <- 0
  rm(L)
  #H <- Lr %*% t(Lc)
  #A <- t(Lc) %*% Lr
  if(as.character(score=="hub")){
    M <- Lr %*% t(Lc)
  }else{
    M <- t(Lc) %*% Lr
  }
  rm(Lr, Lc)
  res <- eigen(M)$values
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name
  }  
  res[vids]
}