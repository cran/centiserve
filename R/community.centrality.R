#' Find the community-based node centrality
#'
#' This function returns community-based node centrality measures. 
#' @details 
#' The "commweight" type weights each community that a node belongs to by how similar that community is to each of the other communities to which the node also belongs. 
#' For node \eqn{i}{i} the community centrality is:
#' \deqn{C_{c}(i)=\sum_{i \in j}^{N}(1 - \frac{1}{m}\sum_{i \in j\cap k}^{m}S(j,k))}{C_c(i)=sum(1 - 1/m*sum(S(j,k), i in j\cap k, m), i in j, N)}
#' where the main sum is over the N communities to which node \eqn{i}{i} belongs, and \eqn{S(j,k)}{S(j,k)} refers to the similarity between community \eqn{j}{j} and \eqn{k}{k}, calculated as the Jaccard coefficient for the number of shared nodes between each community pair, and this is averaged over the \eqn{m}{m} communities paired with community \eqn{j}{j} and in which node \eqn{i}{i} jointly belongs. \cr
#' The "commconn" type weights each community that a node belongs to by how many connections the community forms outside of itself relative to how many connections the community has within itself (the inverse of modularity), so that nodes that belong to more highly connecting communitites will receive a higher community centrality score. For node i the community centrality is:
#' \deqn{C_{c}(i)=\sum_{i \in j}^{N}e_{ij} \frac{\check{e}_{B(j)}}{\check{e}_{W(j)}}}
#' where \eqn{e_{ij}}{e(ij)} is the number of edges node \eqn{i}{i} has in community \eqn{j}{j}, \eqn{\check{e}_{B(j)}=\frac{e_{B(j)}}{n_{j}\bar{d}}}{e'(B(j))=e(B(j))/n(j)\bar{d}} is the number of edges community \eqn{j}{j} makes outside of itself normalised by the number of nodes in community \eqn{j}{j} multiplied by the average degree in the network, and \eqn{\check{e}_{W(j)}=\frac{e_{W(j)}}{n(n-1)/2}}{e'(W(j))=e(W(j))/(n(n-1)/2)} is the number of edges within community \eqn{j}{j} normalised by the total number possible. \cr
#' For more detail see 'linkcomm' package and \href{http://www.centiserver.org/?q1=centrality&q2=Community_Centrality}{Community Centrality}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param type A character string naming the community centrality measure. Can be one of "commweight" or "commconn"
#' @param normalise Logical, whether to normalise community connectedness for "commconn". Defaults to TRUE. Will be ignored for "commweight".
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @author Code obtained from 'linkcomm' package.
#' @references Kalinka, Alex T., and Pavel Tomancak. "linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type." Bioinformatics 27.14 (2011).
#' @examples
#' g <- random.graph.game(20, 3/10)
#' communitycent(g)
#' @export

communitycent <- function (graph, vids = V(graph), type = c("commweight","commconn"), normalise = TRUE){
  if(suppressWarnings(requireNamespace("linkcomm", logical.return=T, warn.conflicts=F, quietly=T, verbose=F)!=TRUE)){
    stop("The 'linkcomm' package needed for this function to work. Please install it.", call. = FALSE)
  }
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  if(missing(type)){
    type <- "commweight"
  }else{
    if(length(type) > 1 || (as.character(type)!="commweight" && as.character(type)!="commconn")){
      stop("The type should be one of 'commweight' or 'commconn'.", call. = FALSE)
    }  
  }
  
  x <- linkcomm::getLinkCommunities(get.edgelist(graph), directed=is.directed(graph), plot=FALSE, verbose=FALSE)
  nodes <- names(x$numclusters)
  clusterids <- unique(x$nodeclusters[x$nodeclusters[, 1] %in% nodes, 2])

  if (type == "commconn") {
    cc <- linkcomm::getCommunityConnectedness(x, clusterids = clusterids, normalise = normalise)
    res <- rep(0, length(nodes))
    for (i in 1:length(nodes)) {
      for (j in 1:length(cc)) {
        ee <- x$edgelist[x$clusters[[clusterids[j]]], 
                         ]
        nn <- length(unique(c(which(ee[, 1] == nodes[i]), 
                              which(ee[, 2] == nodes[i]))))
        res[i] <- sum(res[i], nn * cc[j])
      }
    }
  }else{
    #type = "commweight"
    res <- NULL  
    for (i in 1:length(nodes)) {
      cids <- as.integer(unique(x$nodeclusters[x$nodeclusters[, 1] %in% nodes[i], 2]))
      if(length(cids) > 1) {
        cR <- 1 - linkcomm::getClusterRelatedness(x, clusterids = cids, cluster = FALSE, verbose = FALSE)
        summed <- 0
        for (j in 1:length(cids)) {
          inds <- .cs.getUpperTriIndices(length(cids), which = j)
          summed <- summed + (1 - mean(cR[inds]))
        }
        res[i] <- 1 + summed
      }else if(length(cids) == 1) {
        res[i] <- 1
      }else{
        res[i] <- 0
      }
    }
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name
  }  
  res[vids]
}

# escape .getUpperTriIndices function!
.cs.getUpperTriIndices <- function (numedg, which = 1){
  rows <- NULL
  cols <- NULL
  k <- numedg - 1
  for (i in 1:(numedg - 1)) {
    rows <- append(rows, rep(i, k))
    k <- k - 1
    cols <- append(cols, (i + 1):numedg)
  }
  ret <- union(which(rows == which), which(cols == which))
  return(ret)
}