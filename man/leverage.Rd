% Generated by roxygen2 (4.0.2): do not edit by hand
\name{leverage}
\alias{leverage}
\title{Find the leverage centrality}
\usage{
leverage(graph, vids = V(graph), mode = c("all", "out", "in"),
  loops = TRUE)
}
\arguments{
\item{graph}{The input graph as igraph object}

\item{vids}{Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices.}

\item{mode}{Character constatnt, it specifies how to use the direction of the edges if a directed graph is analyzed. For 'out' only the outgoing edges are followed. For 'in' all vertices from which the source vertex is reachable in at most order steps are counted. 'all' ignores the direction of the edges. This argument is ignored for undirected graphs.}

\item{loops}{Logical; whether the loop edges are also counted.}
}
\value{
A numeric vector contaning the centrality scores for the selected vertices.
}
\description{
Leverage centrality considers the degree of a node relative to its neighbors and operates under the principle that a node in a network is central if its immediate neighbors rely on that node for information.
}
\details{
Leverage centrality of vertex \eqn{i}{i} defined as:
\deqn{l_{i}=\frac{1}{k_{i}}\sum_{N_{i}}\frac{k_{i}-k_{j}}{k_{i}+k_{j}}}
where \eqn{k_{i}}{k(i)} is degree of a given node \eqn{i}{i}, \eqn{k_{j}}{k(j)} is degree of each of its neighbors and \eqn{N_{i}}{N(i)} is all neighbors. \cr
A node with negative leverage centrality is influenced by its neighbors, as the neighbors connect and interact with far more nodes. A node with positive leverage centrality, on the other hand, influences its neighbors since the neighbors tend to have far fewer connections. \cr
More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Leverage_Centrality}{Leverage Centrality}
}
\examples{
g <- graph(c(1,2,2,3,3,4,4,2))
leverage(g)
}
\author{
Mahdi Jalili \email{m_jalili@farabi.tums.ac.ir}
}
\references{
Joyce, Karen E., et al. "A new measure of centrality for brain networks." PLoS One 5.8 (2010): e12200.
}

