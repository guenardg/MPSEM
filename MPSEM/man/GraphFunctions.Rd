\encoding{utf8}
\name{GraphFunctions}
\alias{ graph-functions }
\alias{ pop.graph }
\alias{ add.vertex }
\alias{ rm.vertex }
\alias{ add.edge }
\alias{ rm.edge }
\alias{ collapse.vertex }
\alias{ Phylo2DirectedGraph }
\title{ Graph creation and manipulation functions }
\description{
  A set of primitive functions for creating and munipulating graphs.
}
\usage{
  pop.graph(n, vertex=list(), label=NULL)
  add.vertex(x, n, vertex=list(), label=NULL)
  add.edge(x, from, to, edge=list(), label=NULL)
  rm.edge(x, id)
  rm.vertex(x, id)
  collapse.vertex(x, id)
  Phylo2DirectedGraph(tp)
}
\arguments{
  \item{ x }{ A \code{graph-class} object. }
  \item{ n }{ Number of vertex to populate a new graph
    (\code{\link{pop.graph}}) or to add to an existing graph
    (\code{\link{add.vertex}}). }
  \item{ vertex }{ List of vertex properties. }
  \item{ edge }{ List of edge properties. }
  \item{ label }{ Labels to be given to edges or vertices. }
  \item{ from }{ Origin of the edge to be added (vertex labels or
    indices). }
  \item{ to }{ Destination of the edge to be added (vertex labels or
    indices). }
  \item{ id }{ Label or index of vertex or edge to be removed. }
  \item{ tp }{ Phylogenetic tree object of class \sQuote{phylo}, as
    defined in \code{\link{ape-package}}. }
}
\details{
  A new graph can be populated with \code{n} vertices using function
  \code{\link{pop.graph}} and vertices can be added later with function
  \code{\link{add.vertex}}. The graphs so created contain no edges; the
  latter are added using function \code{\link{add.edge}}. Vertices and
  edges are removed using functions \code{\link{rm.vertex}} and
  \code{\link{rm.edge}}, respectively.
  
  Function \code{\link{collapse.vertex}} allows one to remove a vertex
  while reestablishing the connections between the vertices located
  above and below that vertex using a new set of edges.
  
  Function \code{\link{Phylo2DirectedGraph}} uses the graph functions to
  convert a rooted phylogenetic tree of class \sQuote{phylo} (see
  \code{\link{ape-package}}) to a directed graph object of
  \code{\link{graph-class}}. It recycles tip labels and creates default
  node labels, if they were absent from the \code{phylo} object, and
  uses them as vertex labels. The resulting acyclic graph (i.e. a
  mono-phylogeny) can then be edited to represent cases that do not have
  a tree topology (poly-phylogenies).
}
\value{
  A \code{\link{graph-class}} object. Objects returned by
  \code{\link{Phylo2DirectedGraph}} have a \code{\link{numeric}} edge
  property called \sQuote{distance} featuring branch lengths and a
  \code{link{logical}} vertex property called \sQuote{species}
  specifying whether a vertex is a tree tip or an internal node.
}
\author{ Guillaume Guénard, Département de sciences biologiques,
  Université de Montréal, Montréal, QC, Canada. }
\references{
  Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic
  eigenvector maps (PEM): a framework to model and predict species
  traits. Meth. Ecol. Evol. In press.
  
  Makarenkov, V., Legendre, L. & Desdevise, Y. 2004. Modelling
  phylogenetic relationships using reticulated networks. Zool. Scr. 33:
  89-96
  
  Blanchet, F. G., Legendre, P. & Borcard, D. 2008. Modelling
  directional spatial processes in ecological data. Ecol. Model. 215:
  325-336
}
\seealso{
  \code{\link{graph-class}}.
}
\examples{
## Populate a graph with 7 vertices labeled A-G having properties x and y:
gr <- pop.graph(n=7,
               vertex=list(x=rnorm(7,0,1),y=rnorm(7,0,1)),
               label=c("A","B","C","D","E","F","G"))
gr
## Adding 3 vertices H, I, and J with property x (y is absent) and a new
## property z (type character), which is unknown for A-G:
gr <- add.vertex(x=gr,
                 n=3,
                 label=c("H","I","J"),
                 vertex=list(x=rnorm(3,0,1),z=c("A","B","C")))
gr
gr$vertex
### Adding 10 edges, labeled E1-E10 and with properties a and b, to the graph:
gr <- add.edge(x=gr,
               from=c("A","B","B","C","C","D","D","E","E","F"),
               to=c("A","C","D","E","F","F","G","H","I","J"),
               edge=list(a=rnorm(10,0,1),b=rnorm(10,0,1)),
               label=paste("E",1:10,sep=""))
gr
gr$edge
### Removing edges 2, 4, and 7 from the graph:
print(rm.edge(gr,id=c(2,4,7)))
### Removing vertices 1, 3, 7, and 10 from the graph:
print(rm.vertex(gr,id=c(1,3,7,10)))
# Notice that the edges that had one of the removed vertex as their
# origin or destination are also removed:
print.default(rm.vertex(gr,id=c(1,3,7,10)))
### Vertex collapsing.
x <- pop.graph(n=9,label=c("A","B","C","D","E","F","G","H","I"))
x <- add.edge(x,from=c("A","A","B","B","C","C","D","D","E","E"),
                to=c("B","C","D","E","E","I","F","G","G","H"),
              label=paste("E",1:10,sep=""),edge=list(length=c(1,2,3,2,1,3,2,2,1,3)))
print.default(x)
for(i in c("A","B","C","D","E","F","G","H","I"))
  print(collapse.vertex(x,id=i))
#
if(require(ape)) {
  tree1 <- read.tree(text=paste("(((A:0.15,B:0.2)N4:0.15,C:0.35)N2:0.25,((D:0.25,E:0.1)N5:0.3,",
                                "(F:0.15,G:0.2)N6:0.3)N3:0.1)N1;",sep=""))
  x <- Phylo2DirectedGraph(tree1)
  print(x)
}
}
