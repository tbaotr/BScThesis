#########################################################################
# Network Analysis of Search Trajectory Networks (STN)
# Authors: Gabriela Ochoa, Katherine Malan, Christian Blum
# Date: May 20201
# Computing metrics of a merged STN network given as input
# Input:  File .RData containing Merged STN graph object 
# Output: CSV file with metrics saved in current folder
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
  size_factor <- as.integer(args[2])
} else {
  size_factor <- 1
}

if (is.na(size_factor)) {
  stop("2nd argument is not a number", call.=FALSE)
}

#--------------------------------------------------------------------------
# Create dataframe with metrics
# instance:   name of the file 
# nodes:      number of nodes
# edges:      total number of edges
# nshared:    number of nodes visited by more than one algorithm
# nbest:      number of best nodes (nodes with equal or lower than given best evaluation)
# nend:       number of nodes at the end of trajectories (excluding the best nodes)
# components: Number of connected components
# strength:   normalized strength (incoming weighted degree) of best nodes - normalized with the number of runs.  

col_types =  c("character", "integer", "integer", "integer", 
               "integer", "integer", "integer", "integer", 
               "numeric", "numeric")

col_names =  c("instance", "nodes", "edges", "nshared", 
               "nbest_pr", "nbest_f1", "nend", "components", 
               "strength_pr", "strength_f1")

metrics  <- read.table(text="", colClasses=col_types, col.names=col_names)

# ---- Process all datasets in the given input folder ----------------

i <- 1 # index to store in dataframe. Single column dataset

cat ("Input file: ",infile,"\n")
load(infile, verbose=F)

iname <- gsub('.{6}$','',infile) # removes  (last 6 characters, .RData) from file to use as name
metrics[i,"instance"] <- iname
metrics[i,"nodes"] <- vcount(stnm)
metrics[i,"edges"] <- ecount(stnm)

# Take the IDS of the nodes for metric computation
end_nodes <- which(grepl("end", V(stnm)$Type, fixed=TRUE))
bestpr_nodes <- which(grepl("bestpr", V(stnm)$Type, fixed=TRUE))
bestf1_nodes <- which(grepl("bestf1", V(stnm)$Type, fixed=TRUE))

metrics[i,"nshared"] <- length(which(V(stnm)$Shared == TRUE))
metrics[i,"nbest_pr"] <- length(bestpr_nodes)
metrics[i,"nbest_f1"] <- length(bestf1_nodes)
metrics[i,"nend"] <- length(end_nodes)
metrics[i,"components"] <- components(stnm)$no
# Strength metric
bestpr_str <- sum(strength(stnm, vids=bestpr_nodes, mode="in")) # incoming strength of best
bestf1_str <- sum(strength(stnm, vids=bestf1_nodes, mode="in")) # incoming strength of best
metrics[i,"strength_pr"] <- round(bestpr_str/(num_alg*nruns),4) # normalized by total number of runs
metrics[i,"strength_f1"] <- round(bestf1_str/(num_alg*nruns),4) # normalized by total number of runs

# Save metrics as .csv file
ofname <- paste0(iname,"-metrics.csv")   # file name and -metrics
cat ("Output file: ", ofname, "\n" )
write.csv(metrics, file = ofname)