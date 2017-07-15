.cs.checkPreconditions <- function (g, conditions=c()) {
  if (!is.igraph(g)) {
    stop("Not a graph object", call. = FALSE)
  }
  error <- character();
  for(condition in conditions){
    switch( condition,
            connected={
              if( !is.connected(g, mode="weak") )
                error <- append(error, "Graph is not connected.");
            },
            stronglyConnected={
              if( !is.connected(g, mode="strong") )
                error <- append(error, "Graph is not strongly connected.");
            },
            noMultiple={
              if( any(is.multiple(g)) )
                error <- append(error, "Graph have multiple edges.");
            },
            loopFree={
              if( any(is.loop(g)) )
                error <- append(error, "Graph is not loop free.");
            },
            undirected={
              if( is.directed(g) )
                error <- append(error, "Graph is not undirected.");
            },
            directed={
              if( !is.directed(g) )
                error <- append(error, "Graph is not directed.");
            },
            weighted={
              if(is.null(E(g)$weighte))
                error <- append(error, "Graph is not weighted.");
            },
{}
    )
  }
  if(length(error)>0) stop(paste(error, collapse=" "), call. = FALSE);
}
# Blow functions is one of internal igaph package functions (as.igraph.vs). This is for escape "Unexported object imported by a ':::' call" NOTE.
.cs.as.igraph.vs <- function (graph, v, na.ok = FALSE) 
{
  if (is.character(v) && "name" %in% list.vertex.attributes(graph)) {
    v <- as.numeric(match(v, V(graph)$name))
    if (!na.ok && any(is.na(v))) {
      stop("Invalid vertex names")
    }
    v
  }
  else {
    if (is.logical(v)) {
      res <- as.vector(V(graph))[v]
    }
    else if (is.numeric(v) && any(v < 0)) {
      res <- as.vector(V(graph))[v]
    }
    else {
      res <- as.numeric(v)
    }
    if (!na.ok && any(is.na(res))) {
      stop("Invalid vertex name(s)")
    }
    res
  }
}
