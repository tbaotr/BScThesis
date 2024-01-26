#########################################################################
# Network Analysis of Search Trajectory Networks (STN)
# Authors: Gabriela Ochoa, Katherine Malan, Christian Blum
# Date: May 2021
# Computing metrics of STN networks given on a folder 
# Input:  Folder name with files containing STN graph object (RData file)s
# Output: CSV file with metrics saved in current folder
#########################################################################

library(igraph)  # assume it is already installed for producing STNs

# ---------- Processing inputs from command line ----------
# If running in R GUI, turn off the below statement
args = commandArgs(trailingOnly=TRUE)   # Take command line arguments

if (length(args) < 1) { #  Test if there are two arguments if not, return an error
  stop("One rgument is required: the input folder with STN objects.", call.=FALSE)
}

infolder <- args[1]
if (!dir.exists(infolder) ){
  stop("Input folder does not exist", call.=FALSE)
}

#--------------------------------------------------------------------------
# Create dataframe with metrics
# instance:   Name of the file 
# nodes:      Total number of nodes
# edges:      Total number of edges
# nbest:      Number of best nodes (nodes with equal or lower than given best evaluation), zero if none exist
# nend:       Number of nodes at the end of trajectories (excluding the best nodes)
# components: Number of connected components
# The following metrics only apply if the number of best > 0, otherwise they are NA
# strength: Normalized strength (incoming weighted degree) of best nodes - normalized with the number of runs
# plength:  verge of the shortest path length from start nodes to the best node, NA if non best exist
# npaths:   Number of shortest paths to best optima

col_types = c("character", "integer", "integer", "integer", "integer", "integer", 
              "integer", "numeric", "integer", "integer", "numeric", "integer", 
              "integer")

col_names = c("instance", "nodes", "edges", "nbest_pr", "nbest_f1", "nend", 
              "components", "strength_pr", "plength_pr", "npaths_pr", "strength_f1", 
              "plength_f1", "npaths_f1")

metrics <- read.table(text="", colClasses=col_types, col.names=col_names)


# ---- Process all datasets in the given inpath folder ----------------
instances <- list.files(infolder)

i <- 1 # index to store in dataframe
for (inst in instances) {
  print(inst)
  fname <- paste0(infolder,"/",inst)
  load(fname, verbose=F)
  iname <- gsub('.{6}$','',inst) # removes (last 5 characters, .RData) from file to use as name
  metrics[i,"instance"] <- iname
  metrics[i,"nodes"] <- vcount(STN)
  metrics[i,"edges"] <- ecount(STN)
  bestpr_ids <- which(grepl("bestpr", V(STN)$Type, fixed=TRUE)) # ids of best nodes
  metrics[i,"nbest_pr"] <- length(bestpr_ids)
  bestf1_ids <- which(grepl("bestf1", V(STN)$Type, fixed=TRUE)) # ids of best nodes
  metrics[i,"nbest_f1"] <- length(bestf1_ids)
  start_ids <- which(grepl("start", V(STN)$Type, fixed=TRUE))   # ids of start nodes
  end_ids <- which(grepl("end", V(STN)$Type, fixed=TRUE))       # ids of end nodes
  metrics[i,"nend"] <- length(end_ids)
  metrics[i,"components"] <- components(STN)$no
  if (length(bestpr_ids) > 0) { # if there are nodes with best-known evaluation
    bestpr_str <- sum(strength(STN, vids=bestpr_ids, mode="in"))  # incoming strength of best
    metrics[i,"strength_pr"] <- round(bestpr_str/nruns,4)
    dg <- distances(STN, v=start_ids, to=bestpr_ids, mode="out", weights=NULL)
    d <- dg[is.finite(dg)] # Remove Inf values from distance matrix d
    metrics[i,"plength_pr"] <- round(mean(d),4) # average length of shortest path to best
    metrics[i,"npaths_pr"] <- length(d)         # Number of shortest paths to best
  } else {
    metrics[i,"plength_pr"] <- NA   # average length of shortest path to best
    metrics[i,"npaths_pr"] <- 0     # Number of shortest paths to best
  }
  if (length(bestf1_ids) > 0) { # if there are nodes with best-known evaluation
    bestf1_str <- sum(strength(STN, vids=bestf1_ids, mode="in"))  # incoming strength of best
    metrics[i,"strength_f1"] <- round(bestf1_str/nruns,4)
    dg <- distances(STN, v=start_ids, to=bestf1_ids, mode="out", weights=NULL)
    d <- dg[is.finite(dg)] # Remove Inf values from distance matrix d
    metrics[i,"plength_f1"] <- round(mean(d),4) # average length of shortest path to best
    metrics[i,"npaths_f1"] <- length(d)         # Number of shortest paths to best
  } else {
    metrics[i,"plength_f1"] <- NA # average length of shortest path to best
    metrics[i,"npaths_f1"] <- 0   # Number of shortest paths to best
  }
  i <- i + 1
}

# Save metrics as .csv file
# Create out-folder folder to save STN objects -- rule append "-plot" to input folder
ofname <- paste0(infolder,"-metrics.csv")
write.csv(metrics, file = ofname)