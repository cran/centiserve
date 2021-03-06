% Generated by roxygen2 (4.0.2): do not edit by hand
\name{geokpath}
\alias{geokpath}
\title{Find the geodesic k-path centrality}
\usage{
geokpath(graph, vids = V(graph), mode = c("all", "out", "in"),
  weights = NULL, k = 3)
}
\arguments{
\item{graph}{The input graph as igraph object}

\item{vids}{Vertex sequence, the vertices for which the centrality values are returned. Default is all vertices.}

\item{mode}{Character constant, gives whether the shortest paths to or from the given vertices should be calculated for directed graphs. If out then the shortest paths from the vertex, if in then to it will be considered. If all, the default, then the corresponding undirected graph will be used, ie. not directed paths are searched. This argument is ignored for undirected graphs.}

\item{weights}{Possibly a numeric vector giving edge weights. If this is NULL, the default, and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute).}

\item{k}{The k parameter. The default is 3.}
}
\value{
A numeric vector contaning the centrality scores for the selected vertices.
}
\description{
Geodesic K-path centrality counts neighbours as those that are on a geodesic path less than "k" away.
}
\details{
More detail at \href{http://www.centiserver.org/?q1=centrality&q2=Geodesic_K-Path_Centrality}{Geodesic K-Path Centrality}
}
\examples{
g <- barabasi.game(100)
geokpath(g)
}
\author{
Mahdi Jalili \email{m_jalili@farabi.tums.ac.ir}
}
\references{
Borgatti, Stephen P., and Martin G. Everett. "A graph-theoretic perspective on centrality." Social networks 28.4 (2006): 466-484.
}

