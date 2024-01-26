#########################################################################
# Network Analysis of Search Trajectory Networks (STN)
# Authors: Gabriela Ochoa, Katherine Malan, Christian Blum
# Date: May 2021
# Visualisation of merged STN network of several algorithms
# Input:  File name with merged STN graph object (RData file)
# Output: Network plots (pdf) saved in current folder
#########################################################################

library(igraph)  # assume it is already installed for producing STNs

# ---------- Processing inputs from command line ----------
# If running in R GUI, turn off the below statement
args = commandArgs(trailingOnly=TRUE)   # Take command line arguments

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

plotNet <- function(N, tit, nsizef, ewidthf, asize, ecurv, mylay, bleg=T) 
{
  nsize <- nsizef * V(N)$size
  ewidth <- ewidthf * E(N)$width
  print(tit)
  par(mar = c(1, 1, 1, 1))
  # plot(N, layout=mylay, vertex.label=round(V(N)$PR*6,0),
  #      vertex.label.cex=1.0, vertex.label.dist=0.9, vertex.label.color="black",
  #      vertex.size=nsize,
  #      edge.width=ewidth, edge.arrow.size=asize,
  #      edge.curved=ecurv, loop.size=0.3)
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
  
  V(N)[start_nodes]$color = start_ncol # Color of start nodes
  V(N)[end_nodes]$color = end_ncol     # Color of end of runs nodes
  V(N)[best_nodes]$color = best_ncol   # Color of best nodes
  
  # Frame colors are the same as node colors, white around best to highlight it
  V(N)$frame.color <- V(N)$color
  V(N)[V(N)$color == shared_col]$frame.color <- "gray40"
  V(N)[best_nodes]$frame.color <- "white"
  
  # Shape of nodes
  V(N)$shape <- "circle"               # circle is the default shape
  V(N)[start_nodes]$shape <- "square"  # Square for start nodes
  V(N)[end_nodes]$shape <- "triangle"  # Triangle for start nodes
  
  # Size of Nodes Proportional to  incoming degree, (remove loops and multiple edges)
  sim_N <- simplify(N, remove.multiple=T, remove.loops=T)
  V(N)$size <- strength(sim_N, mode="in") + 1          # nodes with strength 0 have at least size 0.8 
  V(N)$size <- ifelse(V(N)$size > 3, 2.5, V(N)$size)

  V(N)[end_nodes]$size <- V(N)[end_nodes]$size + 0.3   # Increase a a bit size of end nodes
  V(N)[best_nodes]$size <- V(N)[best_nodes]$size + 0.6 # Increase a bit more the size of  best nodes
  
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

# Plot the merged Network
load(infile, verbose=F)

# Plot legend
legend.col <- c("gold", "gray30", "red", alg_col[1:length(algn)], shared_col)
legend.shape <- c(15,17,16,16,16,16,16,16) # square and circles
legend.txt <- c("Start", "End", "Best", algn, "Shared") # needs to read names of algorithms

stnm_pr <- stn_decorate(stnm, 0)
stnm_f1 <- stn_decorate(stnm, 1)
lkk <-layout.kamada.kawai(stnm_f1)
lfr <-layout.fruchterman.reingold(stnm_f1)

# Produce a sub-graph of the merged STN by pruning by fitness top 25%
zoom_pr <- subTop(stnm_pr, 0, as.numeric(quantile(V(stnm_pr)$PR)[4])) # prune by fitness top 25%, quantile [3]
#zoom_pr <- delete.vertices(zoom_pr, degree(zoom_pr)==0) # Remove isolated nodes
lzfr_pr <- layout.fruchterman.reingold(zoom_pr)
lzkk_pr <- layout.kamada.kawai(zoom_pr)

zoom_f1 <- subTop(stnm_f1, 1, as.numeric(quantile(V(stnm_f1)$F1)[4])) # prune by fitness top 25%, quantile [3]
#zoom_f1 <- delete.vertices(zoom_f1, degree(zoom_f1)==0) # Remove isolated nodes
lzfr_f1 <- layout.fruchterman.reingold(zoom_f1)
lzkk_f1 <- layout.kamada.kawai(zoom_f1)


ofname <- gsub('.{6}$', '', infile) # removes (last 6characters) .RData from file to use as name
ofname <- paste0(ofname,'-plot.pdf')

pdf(ofname, width=7, height=6) 
print(ofname)

nf <- size_factor
ef <- size_factor*.7  # Edges width factor is 70% of nodes factor

plotNet(stnm_pr, tit="", nsizef=nf, ewidthf=ef, asize=.16, ecurv=.3, mylay=lfr)
plotNet(stnm_f1, tit="FR layout (F1)", nsizef=nf, ewidthf=ef, asize=.16, ecurv=.3, mylay=lfr)
# Slightly smaller nodes and edges for the KK layout as it spreads the components and makes the nodes closer to each other
plotNet(stnm_pr, tit="KK Layout (PR)", nsizef=nf*.8, ewidthf=ef*.8, asize=.12, ecurv=.3, mylay=lkk)
plotNet(stnm_f1, tit="KK Layout (F1)", nsizef=nf*.8, ewidthf=ef*.8, asize=.12, ecurv=.3, mylay=lkk)

# Plots of zoomed  -- Increased summed size and edges as the zoomed network has less nodes 
nf <- size_factor*1.5  
ef <- size_factor

plotNet(zoom_pr, tit="Zoomed (top 25%) FR (PR)", nsizef=nf, ewidthf=ef, asize=.35, ecurv=.3, mylay=lzfr_pr, bleg=T)
plotNet(zoom_f1, tit="Zoomed (top 25%) FR (F1)", nsizef=nf, ewidthf=ef, asize=.35, ecurv=.3, mylay=lzfr_f1, bleg=T)
plotNet(zoom_pr, tit="Zoomed (top 25%) KK (PR)", nsizef=nf*.8, ewidthf=ef*.8, asize=.2, ecurv=.3, mylay=lzkk_pr, bleg=T)
plotNet(zoom_f1, tit="Zoomed (top 25%) KK (F1)", nsizef=nf*.8, ewidthf=ef*.8, asize=.2, ecurv=.3, mylay=lzkk_f1, bleg=T)

dev.off()

print("Merged STN number of nodes:")
print(vcount(stnm))