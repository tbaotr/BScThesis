#########################################################################
# Network Analysis Search Trajectory Networks
# Author: Gabriela Ochoa
# Date: May 2021
# STN Visualization, Includes decoration
# Plots for single algorithm, as opposed to merging several algorithms
# Input:  Folder with STN graph objects for a single algorithm
# Output: Network plots (pdf) saved in folder
#########################################################################

# ---------- Processing inputs from command line ----------
# If running in R GUI, turn off the below statement
args = commandArgs(trailingOnly=TRUE) # Take command line arguments

if (length(args) < 1) { #  Test if there are two arguments if not, return an error
  stop("One rgument is required: the input folder with stn objects. \\
        A 2nd argument can be given, a size numeric factor for nodes and edges", call.=FALSE)
}

infolder <- args[1]
if (!dir.exists(infolder) ){
  stop("Input folder does not exist", call.=FALSE)
}

# Create out-folder folder to save STN objects -- rule append "-plot" to input folder
outfolder <- paste0(infolder,"-plot")
if (!dir.exists(outfolder) ){
  dir.create(outfolder)
}
cat("Output folder: ", outfolder, "\n")

if (length(args) > 1) {
  size_factor <- as.numeric(args[2])
} else {
  size_factor <- 1
}

if (is.na(size_factor)) {
  stop("2nd argument is not a number", call.=FALSE)
}

if (!require("igraph", character.only = TRUE)) {
  install.packages("igraph", dependencies = TRUE)
  library("igraph", character.only = TRUE)
}

# Functions for visualization, and Default Colors -------------------------------- 

# Node Colors
best_ncol <- "red"   # Best solution found
std_ncol <- "gray70" # Local Optima Gray
end_ncol <- "gray30" # End of trajectories for each run
start_ncol <- "gold" # Start of trajectories

# Edge Colors
# STNs  model has 3 types of perturbation edges: 3 Types: (i)improvement, (e)equal, (w)worsening
# alpha is for transparency: (as an opacity, 0 means fully transparent,  max (255) opaque)
impru_ecol <- "gray50"
equal_ecol <- rgb(0,0,250, max=255, alpha=180)  # transparent blue for worsening edges
worse_ecol <- rgb(0,250,0, max=255, alpha=180)  # transparent green worsening edges

#################################################################
# Triangle vertex shape: because igraph does not have a native
# triangle shape, a function is provided to have a triangle shape
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


# Legend
legend.txt <- c("Start", "End", "Medium", "Best", "Improve", "Equal", "Worse")
legend.col <- c(start_ncol,end_ncol,std_ncol,best_ncol,impru_ecol,equal_ecol,worse_ecol)
legend.shape <- c(15,17,21,21,NA,NA,NA) # Circles for nodes and NA (no shape) for edges
legend.lty <- c(NA,NA,NA,NA,1,1,1)      # Line style, NA for nodes, solid line for edges

# Plot Networks 
# N:       Graph object
# tit:     type of 
# ewidthf: factor for vector with
# asize:   arrow size for plots
# nsizef:  factor to multiply for node sizes
# ecurv:   curvature of the edges (0 = non, 1 = max)
# mylay:   graph layout as a parameter, as different situations require diff. layouts

plotNet <- function(N, tit, nsizef, ewidthf, asize, ecurv, mylay) 
{
  maxns <- max(V(N)$size)
  if (maxns > 100) { # for very large nodes use sqrt for scaling size
    nsize <- nsizef * sqrt(V(N)$size) + 1
  } else {
    if (maxns  > 10) { # for large nodes, use half of size
      nsize <- nsizef * 0.5*V(N)$size + 1
    } else {
      nsize <- nsizef * V(N)$size
    }  
  }
  ewidth <- ewidthf * E(N)$width
  title <- paste(tit, 'Nodes:',vcount(N), 'Edges:',ecount(N), 'Comp:',components(N)$no)
  
  par(mfrow=c(1,1))
  
  print(title)
  plot(N, layout=mylay, vertex.label="", vertex.size=nsize, main=title,
       edge.width=ewidth, edge.arrow.size=asize, edge.curved=ecurv,
       loop.size=0.3)
  
  legend(-1.2,1, legend.txt, pch=legend.shape, col=legend.col, 
         pt.bg=legend.col, lty=legend.lty, cex=0.7, pt.cex=1.35, bty="n")
}

plot2Net <- function(N1, N2, tit, nsizef, ewidthf, asize, ecurv, mylay) 
{
  ## FIRST GRAPH (PR)
  maxns_1 <- max(V(N1)$size)
  if (maxns_1 > 100) { # for very large nodes use sqrt for scaling size
    nsize_1 <- nsizef * sqrt(V(N1)$size) + 1
  } else {
    if (maxns_1  > 10) { # for large nodes, use half of size
      nsize_1 <- nsizef * 0.5*V(N1)$size + 1
    } else {
      nsize_1 <- nsizef * V(N1)$size
    }  
  }
  ewidth_1 <- ewidthf * E(N1)$width
  title_1 <- 'Recall'
  
  # SECOND GRAPH (SR)
  maxns_2 <- max(V(N2)$size)
  if (maxns_2 > 100) { # for very large nodes use sqrt for scaling size
    nsize_2 <- nsizef * sqrt(V(N2)$size) + 1
  } else {
    if (maxns_2  > 10) { # for large nodes, use half of size
      nsize_2 <- nsizef * 0.5*V(N2)$size + 1
    } else {
      nsize_2 <- nsizef * V(N2)$size
    }  
  }
  ewidth_2 <- ewidthf * E(N2)$width
  title_2 <- 'Precision'
  
  par(mfrow=c(1,2), mar=c(0, 1, 1, 1), xpd=NA)
  
  print(title_1)
  plot(N1, layout=mylay,
       vertex.label=round(V(N1)$PR*18, 0), 
       vertex.label.cex=1.0, vertex.label.dist=1.0, vertex.label.color="black",
       vertex.size=nsize_1, main=title_1,
       edge.width=ewidth_1, edge.arrow.size=asize, edge.curved=ecurv,
       loop.size=0.3)
  
  print(title_2)
  plot(N2, layout=mylay,
       vertex.label=round(V(N2)$SR, 2), 
       vertex.label.cex=1.0, vertex.label.dist=1.0, vertex.label.color="black",
       vertex.size=nsize_2, main=title_2,
       edge.width=ewidth_2, edge.arrow.size=asize, edge.curved=ecurv,
       loop.size=0.3)
  
  legend(-3.4,-1.15, legend.txt, pch=legend.shape, col=legend.col, horiz=TRUE,
         pt.bg=legend.col, lty=legend.lty, cex=1.2, pt.cex=1.5, bty="n")
}

# Decorate nodes and edges an STN for visualizing a single algorithm STN
# N: Graph object
# type: type of graph drawn by PR(0), SR(1), F1(2)

stn_decorate <- function(N, type) {
  el <- as_edgelist(N)
  names <- V(N)$name
  
  ## get the objective values at each endpoint of an edge
  if (type == 0) {
    prs <- V(N)$PR
    f1 <- prs[match(el[,1],names)]
    f2 <- prs[match(el[,2],names)]
  }
  else if (type == 1) {
    srs <- V(N)$SR
    f1 <- srs[match(el[,1],names)]
    f2 <- srs[match(el[,2],names)]
  }
  else if (type == 2) {
    f1s <- V(N)$F1
    f1 <- f1s[match(el[,1],names)]
    f2 <- f1s[match(el[,2],names)]
  }
  E(N)[which(f2>f1)]$Type = "improving" # improving edges
  E(N)[which(f2<f1)]$Type = "worsening" # worsening edges
  E(N)[which(f2==f1)]$Type = "equal"    # equal fitness edges
  
  # Coloring nodes and edges. Also give size to nodes
  E(N)$color[E(N)$Type=="improving"] = impru_ecol
  E(N)$color[E(N)$Type=="equal"] = equal_ecol
  E(N)$color[E(N)$Type=="worsening"] = worse_ecol
  
  # width of edges proportional to weight - times visited (limit width of self-loop edges)
  E(N)$width <- E(N)$weight
  E(N)[which_loop(N)]$width <- ifelse(E(N)[which_loop(N)]$weight > 3, 2.5, E(N)[which_loop(N)]$weight)
  
  start_nodes <- grepl("start", V(N)$Type, fixed=TRUE)
  end_nodes <- grepl("end", V(N)$Type, fixed=TRUE)
  if (type == 0 || type == 1)
  {
    best_nodes <- grepl("bestpr", V(N)$Type, fixed=TRUE)
  }
  else if (type == 2)
  {
    best_nodes <- grepl("bestf1", V(N)$Type, fixed=TRUE)
  }
  
  # Color of Nodes
  V(N)$color <- std_ncol               # default color of nodes
  V(N)[start_nodes]$color = start_ncol # Color of start nodes
  V(N)[end_nodes]$color = end_ncol     # Color of end of runs nodes
  V(N)[best_nodes]$color = best_ncol   # Color of best nodes
  
  # Shape of nodes
  V(N)$shape <- "circle"              # circle is the default shape
  V(N)[start_nodes]$shape = "square"  # Square for start nodes
  V(N)[end_nodes]$shape = "triangle"  # Triangle for start nodes
  
  # Frame colors are the same as node colors. White  frame for best nodes to highlight them
  V(N)$frame.color <- V(N)$color
  V(N)[best_nodes]$frame.color <- "white"
  
  # Size of Nodes Proportional to incoming degree, (remove loops and multiple edges)
  sim_N <- simplify(N, remove.multiple=T, remove.loops=T)
  V(N)$size <- strength(sim_N, mode="in") + 1 # nodes with strength 0 have at least size 1 
  V(N)[best_nodes]$size = V(N)[best_nodes]$size + 0.5 # Increase a bit size of best node
  
  return(N)
}

stn_plot <- function(inst) {
  print(inst)
  fname <- paste0(infolder,"/",inst)
  load(fname, verbose=F)
  STN_pr <- stn_decorate(STN, 0)
  STN_sr <- stn_decorate(STN, 1)
  STN_f1 <- stn_decorate(STN, 2)
  lkk <- layout.kamada.kawai(STN_f1)
  lfr <- layout.fruchterman.reingold(STN_f1)
  
  fname <- gsub('.{5}$','',inst) # removes (last 5 characters, RData) from file to use as name
  fname <- paste0(outfolder,"/",fname,"pdf")
  pdf(fname, width=10, height=6)
  print(fname)
  plot2Net(STN_pr, STN_sr, tit="FR Layout.", nsizef=size_factor, ewidthf=size_factor*.5, asize=.3, ecurv=.3, mylay=lfr)
  plot2Net(STN_pr, STN_sr, tit="KK Layout.", nsizef=size_factor*.8, ewidthf=size_factor*.4, asize=.2, ecurv=.3, mylay=lkk)
  dev.off()
  return(vcount(STN))
}

# ---- Process all datasets in the given inpath folder ----------------
dataf <- list.files(infolder)
print(infolder)
print(dataf)
nsizes <- lapply(dataf,stn_plot) # plot all STNs in folder, return node sizes
print("Numer of nodes in STNs ploted:")
print(as.numeric(nsizes))