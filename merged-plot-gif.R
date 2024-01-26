#########################################################################
# Network Analysis of Search Trajectory Networks (STN)
# Authors: Gabriela Ochoa, Katherine Malan, Christian Blum
# Date: May 2021
# Visualisation of merged STN network of several algorithms
# Input:  File name with merged STN graph object (RData file)
# Output: Network plots (pdf) saved in current folder
#########################################################################

library(igraph)  # assume it is already installed for producing STNs
library(animation)

# ---------- Processing inputs from command line ----------
# If running in R GUI, turn off the below statement
#args = commandArgs(trailingOnly=TRUE)   # Take command line arguments

if (length(args) < 1) { #  Test if there are two arguments if not, return an error
  stop("One rgument is required: the input folder with stn objects. \\
        A 2nd argument can be given, a size factor to scale nodes and edges", call.=FALSE)
}
infile <- args[1]

if (!file.exists(infile) ){
  stop("Input file does not exist", call.=FALSE)
}

if (length(args) > 1) {
  size_factor <- as.numeric(args[2])
} else {
  size_factor <- 1
}

if (is.na(size_factor)) {
  stop("2nd argument is not a number", call.=FALSE)
}

# Functions for visualization, and Default Colors -------------------------------- 

# Node Colors
best_ncol <- "red"     # Best solution found
end_ncol <- "gray30"   # End of trajectories for each run.
start_ncol <- "gold"   # Start of trajectories
shared_col <- "gray70" #  Visited by more than one algorithms

# Algorithm colors - used for single algorithms - Algorithms will be colored
# in alphabetical order of their name: blue, purple, green, cyan
alg_col <- c("#377eb8", "#633974", "#229954", "#77D9EB")
# alg_col <- c("#229954", "#77D9EB")

#################################################################
# triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color, col=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip, plot=mytriangle)

# Plot Networks 
# N:        Graph object
# tit:      title for type of plot
# ewidthf:  factor for vector with 
# asize:    arrow size for plots
# nsizef:   factor to multiply for node sizes
# ecurv:    curvature of the edges (0 = non, 1 = max)
# mylayout: graph layout as a parameter, as different situations require diff. layouts
# bleg:     add a legend to the plot

plotNet <- function(N, tit, nsizef, ewidthf, asize, ecurv, mylay, bleg=F) 
{
  nsize <- nsizef * V(N)$size
  ewidth <- ewidthf * E(N)$width
  print(tit)
  par(mar = c(1, 1, 1, 1))
  plot(N, layout=mylay, vertex.label="", vertex.size=nsize,
       edge.width=ewidth, edge.arrow.size=asize,
       edge.curved=ecurv, loop.size=0.3)
  if (bleg == T) {
    legend("topleft", legend.txt, pch=legend.shape, col=legend.col, 
           cex=0.9, pt.cex=1.4, bty="n")
  }
}

# Decorate merged STN - combining algorithms

# sfac: multiplicative factor for size of nodes
# type: type of graph drawn by PR(0), F1(1)
stn_decorate <- function(N, type)  {
  # Decoration of nodes- best evaluation read from RData file
  # Color of Nodes 
  V(N)$color <- shared_col # default color of nodes
  for (i in 1:num_alg) {   # Assign color by algorithm order in vector - only color nodes visited by a single alg.
    V(N)[V(N)$Alg == algn[i]]$color <- alg_col[i]
  }
  
  # Take the IDS of the type of nodes for decoration
  start_nodes <- grepl("start", V(N)$Type, fixed=TRUE)
  end_nodes <- grepl("end", V(N)$Type, fixed=TRUE)
  if (type == 0)
  {
    best_nodes <- grepl("bestpr", V(N)$Type, fixed=TRUE)
  }
  else if (type == 1)
  {
    best_nodes <- grepl("bestf1", V(N)$Type, fixed=TRUE)
  }
  
  V(N)[start_nodes]$color = start_ncol  # Color of start nodes
  V(N)[end_nodes]$color = end_ncol  # Color of end of runs nodes
  V(N)[best_nodes]$color = best_ncol    # Color of best nodes
  
  # Frame colors are the same as node colors, white around best to highlight it
  V(N)$frame.color <- V(N)$color
  print(V(N)$color)
  V(N)[V(N)$color == shared_col]$frame.color <- "gray40"
  V(N)[best_nodes]$frame.color <- "white"
  
  # Shape of nodes
  V(N)$shape <- "circle"               # circle is the default shape
  V(N)[start_nodes]$shape <- "square"  # Square for start nodes
  V(N)[end_nodes]$shape <- "triangle"  # Triangle for start nodes
  
  # Size of Nodes Proportional to  incoming degree, (remove loops and multiple edges)
  sim_N <- simplify(N, remove.multiple=T, remove.loops=T)
  V(N)$size <- strength(sim_N, mode="in") + 1   # nodes with strength 0 have at least size 0.8 
  V(N)[end_nodes]$size <- V(N)[end_nodes]$size + 0.3 # Increase a a bit size of end nodes
  V(N)[best_nodes]$size <- V(N)[best_nodes]$size + 0.6   # Increease a bit more the size of  best nodes
  
  # Color of edges
  E(N)$color <- shared_col # default color of edges
  for (i in 1:num_alg) {
    E(N)[E(N)$Alg==algn[i]]$color <- alg_col[i]
  }
  # width of edges proportional to weight - times visited (limit width of self-loop edges)
  E(N)$width <- E(N)$weight
  E(N)[which_loop(N)]$width <- ifelse(E(N)[which_loop(N)]$weight > 3, 2.5, E(N)[which_loop(N)]$weight)
  return(N)
}


#------------------------------------------------------------------------
# Creates a sub-network with nodes with and below a given fitness level 

# type: type of graph drawn by PR(0), F1(1)
subTop <- function(N, type, threshold)
{
  if (type == 0)
  {
    Top <- induced.subgraph(N, V(N)$PR>=threshold)
  }
  else if (type == 1)
  {
    Top <- induced.subgraph(N, V(N)$F1>=threshold)
  }
  return (Top)
}

load(infile, verbose=F)

# Plot legend
legend.col <- c("gold", "gray30", "red", alg_col[1:length(algn)], shared_col)
legend.shape <- c(15,17,16,16,16,16,16,16) # square and circles
legend.txt <- c("Start", "End", "Best", algn, "Shared") # needs to read names of algorithms

stnm <- stn_decorate(stnm, 0)
lfr <-layout.fruchterman.reingold(stnm)

nf <- size_factor
ef <- size_factor*.7  # Edges width factor is 70% of nodes factor

# plotNet(stnm, tit="", nsizef=nf, ewidthf=ef, asize=.16, ecurv=.3, mylay=lfr)

ani.options(interval=1)
saveGIF({
  total_time <- max(E(stnm)$Time)
  dt <- 1
  for(t in seq(3, total_time, dt)){
    net <- delete_edges(stnm, which(E(stnm)$Time > t))
    
    iso <- which(degree(net)==0)
    net <- delete.vertices(net, iso)
    if (length(iso) != 0) {
      l <- lfr[-as.numeric(iso),]
    }
    else {
      l <- lfr
    }
    
    plotNet(net, tit="", nsizef=nf, ewidthf=ef, asize=.16, ecurv=.3, mylay=l)
  }
}, interval = 0.1, movie.name = "demo.gif", ani.width = 1000, ani.height = 1000)
