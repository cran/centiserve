#' Find the density of maximum neighborhood component (DMNC) in a graph
#'
#' The score of node \eqn{v}{v}, \eqn{DMNC(v)}{DMNC(v)}, is defined to be \eqn{\frac{E}{N^{\epsilon}}}{E/(N^epsilon)}:
#' \deqn{\frac{\left|E(MNC(v))\right|}{\left|V(MNC(v))\right|^{\epsilon}}}{|E(MNC(v))|/(|V(MNC(v))|^epsilon)}
#' where for some \eqn{1 \leq \epsilon \leq 2}{1 <= epsilon <= 2}.
#' @details 
#' See Maximum Neighborhood Component (MNC) \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=DMNC-Density_of_Maximum_Neighborhood_Component}{DMNC-Density of Maximum Neighborhood Component}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param mode Character constatnt, it specifies how to use the direction of the edges if a directed graph is analyzed. For 'out' only the outgoing edges are followed. For 'in' all vertices from which the source vertex is reachable in at most order steps are counted. 'all' ignores the direction of the edges. This argument is ignored for undirected graphs.
#' @param epsilon \eqn{\epsilon} parameter which default is 1.67.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Lin, Chung-Yen, et al. "Hubba: hub objects analyzer-a framework of interactome hubs identification for network biology." Nucleic acids research 36.suppl 2 (2008): W438-W443.
#' @examples
#' g <- random.graph.game(20, 3/10)
#' dmnc(g)
#' @export

dmnc <- function (graph, vids = V(graph), mode = c("all", "out", "in"), epsilon = 1.67){
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  epsilon <- as.double(epsilon)
  if(epsilon < 1 || epsilon > 2) stop("The epsilon parameter must be between 1 and 2", call. = FALSE)
  res <- double()
  for(v in V(graph)[vids]){
    g <- induced.subgraph(graph, neighborhood(graph, 1, nodes=v, mode=mode[1])[[1]][-1])
    c <- clusters(g,  mode="strong")
    ec <- ifelse(length(c$csize)==0, 0, ecount(induced.subgraph(graph, which(c$membership%in%which(c$csize==max(c$csize))))))
    if(ec == 0){
      res <- append(res, 0)
    }else{
      res <- append(res, ec/(max(c$csize) ^ epsilon))
    }
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }  
  res
}