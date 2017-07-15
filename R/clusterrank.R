#' Find the ClusterRank ranks in a graph
#'
#' Mathematically, the ClusterRank score \eqn{s_{i}}{s_i} of node \eqn{i}{i} is defined as:
#' \deqn{s_{i} = f(c_{i})\sum_{j\in \tau _{i}}(k_{out}^{j}+1)}{s_i = f(c_i) sum(k_out,j+1, j in tau _i)}
#' where the term f(c_{i}) accounts for the effect of i's local clustering and the term '+1' results from the contribution of \eqn{j}{j} itself. \cr
#' Here \eqn{f(c_{i}) = 10^{-c_{i}}}{f(c_i) = 10^(-c_i)} 
#' @details 
#' ClusterRank is a local ranking algorithm which takes into account not only the number of neighbors and the neighbors' influences, but also the clustering coefficient. \cr
#' ClusterRank can also be applied to undirected networks where the superiority of ClusterRank is significant compared with degree centrality and k-core decomposition. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=ClusterRank}{ClusterRank}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param directed Logical scalar, whether to directed graph is analyzed. This argument is ignored for undirected graphs.
#' @param loops Logical; whether the loop edges are also counted.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Chen, Duan-Bing, et al. "Identifying influential nodes in large-scale directed networks: the role of clustering." PloS one 8.10 (2013): e77455.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2,2,5,5,3,4,1,4,3,1,6,6,3,3,6,2,6,5,6))
#' clusterrank(g)
#' @export

clusterrank <- function (graph, vids = V(graph), directed = TRUE, loops = TRUE){
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  res <- double()
  if(directed && is.directed(graph)){
    cc <- .cs.transitivity.directed(graph)
    mode <- "out"
  }else{
    cc <- transitivity(graph, type="local")
    mode <- "all"
  }
  for (v in V(graph)[vids]){
    if(is.nan(cc[v])){
      res <- append(res, NaN)
    }else{
      vNeighbors <- neighborhood(graph, 1, v, mode=mode)[[1]][-1]
      cr <- 0
      for (vv in vNeighbors){
        cr <- cr + degree(graph, v=vv, mode=mode, loops=as.logical(loops)) + 1
      }
      res <- append(res, cr * cc[v])
    }
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}

.cs.transitivity.directed <- function(graph){
  res <- double()
  for(v in V(graph)){
    nei <- neighborhood(graph, 1, nodes=v, mode="out")[[1]][-1]
    if(length(nei) < 2){
      res <- append(res, NaN)
    }else{
      g <- induced.subgraph(graph, nei)
      res <- append(res, ecount(g)/(vcount(g)*(vcount(g)-1)))
    }
  }
  res
}