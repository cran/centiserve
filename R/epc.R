#' Find the edge percolated component (EPC) in a graph
#'
#' For a node \eqn{v}{v} in G, \eqn{EPC(v)}{EPC(v)} is defined as:
#' \deqn{EPC(v)=\frac{1}{\left|v\right|}\sum_{k=1}^{1000}\sum_{t\in e}\delta_{vt}^{k}}{EPC(v)=1/|v|*sum(sum(delta(vt, k), t\in e), k=1, 1000)}
#' Given a threshold \eqn{(0 \leq the threshold \leq 1)}{(0 <= the threshold <= 1)}, we create 1000 reduced network by asigning a random number between 0 and 1 to every edge and remove edges if their associated random numbers are less than the threshold. \cr
#' Let the \eqn{G_{k}}{G(k)} be the reduced network generated at the \eqn{k_{th}}{k(th)} time reduced process. If nodes \eqn{u}{u} and \eqn{v}{v} are connected in \eqn{G_{k}}{G(k)}, set \eqn{\delta_{vt}^{k}}{delta(vt)^k} to 1; otherwise \eqn{\delta_{vt}^{k}=0}{delta(vt)^k=0}.
#' @details 
#' For an interaction network G, assign a removing probability p to every edge. Let G'be a realization of the random edge removing from G. If nodes \eqn{v}{v} and \eqn{w}{w} are connected in G', set \eqn{d_{vw}}{d(vw)} be 1, otherwise set \eqn{d_{vw}}{d(vw)} be 0. The percolated connectivity of \eqn{v}{v} and \eqn{w}{w}, \eqn{c_{vw}}{c(vw)}, is defined to be the average of \eqn{d_{vw}}{d(vw)} over realizations. The size of percolated component containing node \eqn{v}{v}, \eqn{s_{v}}{s(v)}, is defined to be the sum of \eqn{c_{vw}}{c(vw)} over nodes \eqn{w}{w}. The score of node \eqn{v}{v}, \eqn{EPC(v)}{EPC(v)}, is defined to be \eqn{s_{v}}{s(v)}. \cr
#' More detail at \href{http://www.centiserver.org/?q1=centrality&q2=EPC-Edge_Percolated_Component}{EPC-Edge Percolated Component}
#' @param graph The input graph as igraph object
#' @param vids Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices. 
#' @param threshold The threshold parameter, for filter graph and create reduced one, which must be between 0 and 1. The default is 0.5.
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Mahdi Jalili \email{m_jalili@@farabi.tums.ac.ir}
#' @references Lin, Chung-Yen, et al. "Hubba: hub objects analyzer-a framework of interactome hubs identification for network biology." Nucleic acids research 36.suppl 2 (2008): W438-W443.
#' @references Chen, Shu-Hwa, et al. "cyto-Hubba: A Cytoscape plug-in for hub object analysis in network biology." 20th International Conference on Genome Informatics. 2009.
#' @examples
#' g <- graph(c(1,2,2,3,3,4,4,2))
#' epc(g)
#' @export

epc <- function (graph, vids = V(graph), threshold = 0.5){
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  threshold <- as.double(threshold)
  if(threshold < 0 || threshold > 1) stop("The threshold parameter must be between 0 and 1", call. = FALSE)
  epcT <- integer(vcount(graph))
  for(i in 1:1000){
    g <- graph
    E(g)$weight <- runif(ecount(g), 0, 1)
    g <- g - edge(E(g)[E(g)$weight < threshold])
    membership <- clusters(g)$membership
    epc <- integer()
    for(v in V(graph)){
      epcv <- 0
      for(t in V(graph)){
        if(membership[v] == membership[t]) epcv <- epcv + 1
      }
      epc <- append(epc, epcv)	
    }
    epcT <- epcT + epc
  }
  epcT <- epcT / vcount(graph)
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(epcT) <- V(graph)$name
  }  
  epcT[vids]
}