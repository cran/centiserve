% Generated by roxygen2 (4.0.2): do not edit by hand
\name{closeness.residual}
\alias{closeness.residual}
\title{Find the residual closeness centrality}
\usage{
closeness.residual(graph, vids = V(graph), mode = c("all", "out", "in"),
  weights = NULL)
}
\arguments{
\item{graph}{The input graph as igraph object}

\item{vids}{Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices.}

\item{mode}{Character constant, gives whether the shortest paths to or from the given vertices should be calculated for directed graphs. If out then the shortest paths from the vertex, if in then to it will be considered. If all, the default, then the corresponding undirected graph will be used, ie. not directed paths are searched. This argument is ignored for undirected graphs.}

\item{weights}{Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).}
}
\value{
A numeric vector contaning the centrality scores for the selected vertices.
}
\description{
Residual closeness centrality defined as:
\deqn{C_{i}=\sum_{j\neq i}\frac{1}{2^{d(i,j)}}}{C_i=sum(1/2^d(i,j), j!=i)}
}
\details{
This function calculate closeness of a vertex as Dangalchev defination. \cr
More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Residual_Closeness_Centrality}{Residual Closeness Centrality}
}
\examples{
g <- graph(c(1,2,2,3,3,4,4,2))
closeness.residual(g)
}
\author{
Mahdi Jalili \email{m_jalili@farabi.tums.ac.ir}
}
\references{
Dangalchev, Chavdar. "Residual closeness in networks." Physica A: Statistical Mechanics and its Applications 365.2 (2006): 556-564.
}

