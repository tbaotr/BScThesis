#########################################################################
# Network Analysis of Search Trajectory Networks (STN)
# Authors: Gabriela Ochoa, Katherine Malan, Christian Blum
# Date: May 2021
# STN construction for single algorithms. 
# Input:  Folder containing text file trace of runs, number of runs
# Output: STN graph objects - saved in output folder 
#########################################################################

# ---------- Processing inputs from command line ----------
# If running in R GUI, turn off the below statement
args = commandArgs(trailingOnly=TRUE) # Take command line arguments

#  Test if there are two arguments if not, return an error
if (length(args) < 1) {
  stop("The first argument is required, arguments 2 is optional: \ 
       1) Name of the input folder 
       2) The number of runs from the data filesto be used. This should be a number between 1 up to total number of runs within in the raw data files. 
       If no argument is given, the largest run number in the input files is used."
       , call.=FALSE)
}

infolder <- args[1]
if (!dir.exists(infolder) ){
  stop("Error: Input folder does not exist", call.=FALSE)
}

# Default values of parameters if not given in command line.
best <- 1.
nruns <- NA   # Not given om command line, taken from data
if (length(args) > 1){
  nruns <- as.integer(args[2]) 
  if (is.na(nruns)) {
    stop("Error: 2th argument is not a number", call.=FALSE)
  }
}

# Create out-folder folder to save STN objects  -- rule append "-stn" to input folder
outfolder <- paste0(infolder,"-stn")
if (!dir.exists(outfolder) ){
  dir.create(outfolder)
}
cat("Output folder: ", outfolder, "\n")


## Packages required
# igraph: tools handling graph objects
# plyr:   tools for Splitting, Applying and Combining Data

packages = c("igraph", "plyr")

## If a package is installed, it will be loaded. If any are not, 
## the missing package(s) will be installed from CRAN and then loaded.

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

#-----------------------------------------------------------------------------------------------
# Function for Creating the STN of a given instance/algorithm 
# Read data from text input file and construct the STN network model
# Saves the STN object in a file within the out-folder
# Returns number of nodes of the STN

stn_create <- function(instance) {
  fname <- paste0(infolder,"/",instance)
  print(fname)
  file_ext <- substr(fname, nchar(fname)-3, nchar(fname))
  mysep <- ifelse(file_ext == ".csv",",","")
  trace_all <- read.table(fname, header=T, sep=mysep,
                          colClasses=c("integer","integer","numeric","numeric","character","numeric","numeric","character"),
                          stringsAsFactors=F)
  trace_all <- trace_all[trace_all$Run <= nruns,]
  lnodes <- vector("list", nruns)
  ledges <- vector("list", nruns)
  
  # Count the number of successful run
  succ_num = 0
  for (i in (1:nruns)) {
    trace <- trace_all[which(trace_all$Run==i),]
    if (max(trace$PR1) == best | max(trace$PR2) == best) {
      succ_num <- succ_num + 1
    }
  }
  
  # Store the nodes at the end of runs. Need to detect change of Run
  k = 1
  end_ids <- vector()
  start_ids <- vector()
  n <- nrow(trace_all)-1
  start_ids[1] <- trace_all$Sol1[1]
  for (j in (1:n)) {
    if (trace_all$Run[j] != trace_all$Run[j+1]) { # when the run counter changes
      end_ids[k] <- trace_all$Sol2[j]             # keep the name of the end solution
      start_ids[k+1] <- trace_all$Sol1[j+1]       # keep the name of the next start solution
      k = k+1
    }
  }
  end_ids[k] <- trace_all$Sol2[j+1] # last end ID -s considered as the end of the trajectory
  
  end_ids <- unique(end_ids) # only unique nodes required
  start_ids <- unique(start_ids)
  
  for (i in (1:nruns)) { # combine all runs in a single network
    trace <- trace_all[which(trace_all$Run==i),c(-1)] # take by run and remove first column run number
    colnames(trace) <- c("time","pr1","sr1","node1","pr2","sr2","node2") # set simpler names to column
    lnodes[[i]] <- rbind(setNames(trace[,c("node1","pr1","sr1")], c("Node","PR","SR")),
                         setNames(trace[,c("node2","pr2","sr2")], c("Node","PR","SR")))
    ledges[[i]] <- trace[,c("node1","node2","time")]
  }
  
  # combine the list of nodes into one dataframe and
  # group by (Node,PR,SR) to identify unique nodes and count them
  nodes <- ddply((do.call("rbind",lnodes)), .(Node,PR,SR), nrow)
  colnames(nodes) <- c("Node","PR","SR","Count")
  nodesu <- nodes[!duplicated(nodes$Node),] # eliminate duplicates from dataframe, in case node ID us duplicated
  # combine the list of edges into one dataframe and
  # group by (node1,node2) to identify unique edges and count them
  edges <- ddply(do.call("rbind", ledges), .(node1,node2,time), nrow)
  colnames(edges) <- c("Start","End","Time","weight")
  
  STN <- graph_from_data_frame(d=edges, directed=T, vertices=nodesu) # Create graph
  STN <- simplify(STN, remove.multiple=F, remove.loops=F)            # Remove self loops
  
  # Calculate F1 score
  V(STN)$F1 <- 2./(1./V(STN)$PR + 1./V(STN)$SR)
  
  # Best Ids with F1 are also Best Ids with PR, but it is not true in the contrary.
  best_pr_ids <- which(V(STN)$PR==best)
  best_f1_ids <- which(V(STN)$F1==best)
  
  # Four types of nodes, useful for visualization: Start, End, Best and Standard.
  V(STN)$Type <- "medium"  # Default type
  V(STN)[end_ids]$Type <- "end"
  V(STN)[start_ids]$Type <- "start"
  V(STN)[best_pr_ids]$Type <- paste(V(STN)[best_pr_ids]$Type, "bestpr", sep="_")
  V(STN)[best_f1_ids]$Type <- paste(V(STN)[best_f1_ids]$Type, "bestf1", sep="_")
  
  fname <- gsub('.{4}$','',instance) # removes (last 4 characters, .ext) from file to use as name
  fname <- paste0(outfolder,"/",fname,"_stn.RData")
  save(STN,nruns,succ_num,best,file=fname) # Store STN, whether it is a minimization problem and the best-known given
  return(vcount(STN))
}

#--------------------------------------------------------------------------------
# Extracts the required data fro the input file
# Input:  String with name of file
# Output: Data frame with trace data, 2) name of output file .Rdata

get_data <- function(instance) {
  file_ext <- substr(instance, nchar(instance)-3, nchar(instance))
  mysep <- ifelse(file_ext == ".csv",",","")
  trd <- read.table(paste0(infolder,"/",instance), header=T, sep=mysep,
                    colClasses=c("integer","integer","numeric","numeric","character","numeric","numeric","character"),
                    stringsAsFactors=F)
  return (trd)
}


# ---- Process all data files in the given input folder ----------------
data_files <- list.files(infolder)

# This is only executed if the nruns parameters are not given
if (is.na(nruns))  {
  dfs <- lapply(data_files,get_data) # Store data in a list, so is an argument to create function 
  if (is.na(nruns)) {
    l <- lapply(dfs, function(x) {x[c("Run")]}) # Extract Runs
    nruns <- max(unlist(l, recursive=T))
    cat("Number of runs in data:", nruns, "\n")
  }
  remove(dfs)
}


nsizes <- lapply(data_files,stn_create) # Applies stn_create function to all files
print("Number of nodes in the STNs created:")
print(as.numeric(nsizes))