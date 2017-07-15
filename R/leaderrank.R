#' Find the LeaderRank in a directed graph
#'
#' This function find the LeaderRank in a directed graph
#' @details 
#' Given a network consisting of N nodes and M directed links, a ground node connected with every node by a bidirectional link is added. Then, the network becomes strongly connected and consists of N+1 nodes and M+2N links (a bidirectional link is counted as two links with inverse directions). LeaderRank directly applies the standard random walk process to determine the score of every node. Accordingly, if the score of node \eqn{i}{i} at time step \eqn{t}{t} is \eqn{si(t)}{si(i)}, the dynamics can be described by an iterative process as:
#' \deqn{s_{i}(t+1)=\sum_{j=1}^{N+1}\frac{a_{ji}}{k_{j}^{out}}s_{j}(t)}{s_i(t+1)=sum(a(ji)/k(j, out) * s_j(t), j=1, N+1)}
#' where \eqn{a_{ji}}{a(ji)} is the element of the corresponding (N + 1)-dimensional adjacency matrix, which equals 1 if there is a directed link from \eqn{j}{j} to \eqn{i}{i} and 0 otherwise, and \eqn{k_{j}^{out}}{k(j)^out} is the out-degree of node \eqn{j}{j}. The process starts with the initialization where all node scores are 1 and will soon converge to a unique steady state denoted as \eqn{s_{i}^{\infty}, (i = 1, 2, ..., N, N+1)}{s(i)^ infty, (i = 1, 2, ..., N, N+1)}. LeaderRank ranks all nodes according to \eqn{s_{i}^{\infty}}{s(i)^infty}, and the nodes with larger final scores are considered to be more influential in spreading. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=LeaderRank}{LeaderRank}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Lu, Linyuan, et al. "Leaders in social networks, the delicious case." PloS one 6.6 (2011): e21202.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2))
#' leaderrank(g)
#' @import Matrix
#' @export

leaderrank <- function (graph, vids = V(graph)){
  #if(suppressWarnings(library("Matrix", logical.return=T, warn.conflicts=F, quietly=T, verbose=F)!=TRUE)){
  #  stop("The 'Matrix' package needed for this function to work. Please install it.", call. = FALSE)
  #}
  .cs.checkPreconditions(graph, c("directed"))
  vids <- .cs.as.igraph.vs(graph, vids)
  N <- vcount(graph)
  EG1 <- matrix(0, nrow=N, ncol=2)
  EG2 <- matrix(0, nrow=N, ncol=2)
  for (i in 1:N){
    EG1[i,1] = N+1
    EG1[i,2] = i
  }
  EG2[,1] = EG1[,2]
  EG2[,2] = EG1[,1]
  EG1 = rbind(get.edgelist(graph, names=FALSE), EG1, EG2)
  EG1 = Matrix::sparseMatrix(EG1[,1], EG1[,2], x=1)
  EG2 <- matrix(0, nrow=N+1, ncol=2)
  for (j in 1:(N+1)){
    EG2[j,1] = j
    EG2[j,2] = 1/sum(EG1[j,])
  }
  EG2 = Matrix::sparseMatrix(EG2[,1], EG2[,1], x=EG2[,2])
  EG1 = EG2 %*% EG1
  EG2 <- matrix(1, nrow=N+1, ncol=1)
  EG2[N+1,1] <- 0
  error <- 10000
  error_threshold <- 0.00002
  EG1 <- t(EG1)
  while ( error > error_threshold){
    M <- EG2
    EG2 <- EG1 %*% EG2
    error = sum(abs(EG2 - M)/M)/(N+1)
  }
  EG2 <- EG2 + (EG2[N+1]/N)
  EG2 <- EG2[1:N]
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(EG2) <- V(graph)$name
  }  
  EG2[vids]
}