#########################################################################
# Network Analysis of Search Trajectory Networks (STN)
# Authors: Gabriela Ochoa, Katherine Malan, Christian Blum
# Date: May 2021
# Construction of merged STN network of several algorithms
# Input:  Folder containing RData files of algorithm STNs to merge 
#         It is assumed that the folder contains one file per algorithm to merge  
# Output: Merged STN graph object - saved in current directory as an .RData file
#############################################################################

# ---------- Processing input from command line ----------
# If running in R GUI, turn off the below statement
args = commandArgs(trailingOnly=TRUE) # Take command line arguments

if (length(args) < 1) { # Test if there is one argument if not, return an error
   stop("One argument is required: the input folder with STNs files to merge.", call.=FALSE)
}

infolder <- args[1]
if (!dir.exists(infolder) ){
  stop("Input folder does not exist", call.=FALSE)
}

## Packages required
# igraph: tools handling graph objects
# tidyr:  tools to help to create tidy data  (we used function "replace_na")

packages = c("igraph", "tidyr")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      suppressWarnings(library(x, character.only = TRUE))
    }
  }
)

# ---- Read input files given folder ----------------
data_files <- list.files(infolder)
cat("Input folder: ",infolder,"\n")
num_alg <- length(data_files) # Number of algorithms

if (num_alg < 2){
  stop("Number of algorithms to merge can only be larger than 1", call.=FALSE)
}  

alg <- vector(mode="list", length=num_alg)       # to keep the STNS of the algorithms 
algn <- vector(mode="character", length=num_alg) # to keep names of the algorithms 

i <- 1
for (f in data_files) {
  alg_name <- strsplit(f,"_")[[1]][1] # Assume that the name of the algorithm is the first string before the "_"
  print(alg_name)
  fname <- paste0(infolder,"/",f)
  load(fname, verbose=F)
  # Add Algorithm as a property of nodes and edges
  alg_name <- paste(alg_name,"(",succ_num,"/",nruns,")")
  V(STN)$Alg <- alg_name
  E(STN)$Alg <- alg_name
  algn[i] <- alg_name
  alg[[i]] <- STN # keep STN in a list to then proceed with unioning the graphs 
  i <- i + 1
}

################################################################################
# We iteratively add each graph to total graph in turn

stnm <- alg[[1]]
for (i in (2:num_alg)) {
  # Union the separate algorithms graphs
  stnm <- stnm %u% alg[[i]]
  
  ## VERTEX PROCESSING ##
  
  # This handles the three main attributes of nodes: PR, SR and Count
  V(stnm)$PR <- apply(cbind(V(stnm)$PR_1, V(stnm)$PR_2), 1, max, na.rm=TRUE)
  V(stnm)$SR <- apply(cbind(V(stnm)$SR_1, V(stnm)$SR_2), 1, max, na.rm=TRUE)
  V(stnm)$F1 <- apply(cbind(V(stnm)$F1_1, V(stnm)$F1_2), 1, max, na.rm=TRUE)
  V(stnm)$Count <- rowSums(cbind(V(stnm)$Count_1, V(stnm)$Count_2), na.rm=TRUE) 
  
  # This concatenate the strings as a given node can have more than one Type and 
  # more than one Alg (coming from each algorithm STNs).
  V(stnm)$Type_1 <- replace_na(V(stnm)$Type_1,"")
  V(stnm)$Type_2 <- replace_na(V(stnm)$Type_2,"")
  V(stnm)$Alg_1 <- replace_na(V(stnm)$Alg_1,"")
  V(stnm)$Alg_2 <- replace_na(V(stnm)$Alg_2,"")
  dft <- data.frame(V(stnm)$Type_1, V(stnm)$Type_2)
  dfa <- data.frame(V(stnm)$Alg_1, V(stnm)$Alg_2)
  dft <- unite(dft,"Type", remove=T, sep="_")
  dfa <- unite(dfa,"Alg", remove=T, sep="")
  V(stnm)$Type <- as.vector(dft$Type) # Contains a concatenation of types of nodes associated with each node
  V(stnm)$Alg <- as.vector(dfa$Alg)   # Contains a concatenation of Algorithms names that visited each node
  
  # Remove vertex attributes that are no longer needed 
  old_vattr = c("PR_1", "PR_2", "SR_1", "SR_2", "F1_1", "F1_2", "Count_1", "Count_2", 
                "Type_1", "Type_2", "Alg_1", "Alg_2")
  for (i in old_vattr) {
    stnm<-delete_vertex_attr(stnm, name=i)
  }
  
  
  ## EDGE PROCESSING ##
  
  # Weights are summed, aggregating the visits of the combined algorithms
  E(stnm)$weight <- rowSums(cbind(E(stnm)$weight_1, E(stnm)$weight_2), na.rm=TRUE)
  E(stnm)$Time <- apply(cbind(E(stnm)$Time_1, E(stnm)$Time_2), 1, min, na.rm=TRUE)
  
  # This concatenate the strings as a given edge can have more than one more than 
  # one Alg (coming from each algorithm STNs).
  E(stnm)$Alg_1 <- replace_na(E(stnm)$Alg_1,"")
  E(stnm)$Alg_2 <- replace_na(E(stnm)$Alg_2,"")
  dfa <- data.frame(E(stnm)$Alg_1, E(stnm)$Alg_2)
  dfa <- unite(dfa,"Alg", remove = T, sep = "")
  E(stnm)$Alg <- as.vector(dfa$Alg)
  
  # Remove edge attributes that are no longer needed 
  old_eattr = c("weight_1", "weight_2", "Alg_1", "Alg_2", "Time_1", "Time_2")
  for (i in old_eattr) {
    stnm <- delete_edge_attr(stnm, name=i)
  }
}

# Detecting shared nodes, that were visited by more than one algorithm
# Keep an attribute for shared nodes, this is later useful for visualisation and metrics
V(stnm)$Shared <- TRUE # Assume shared node
for (i in 1:num_alg) { #  Detect nodes visited by a single Algorithm
  V(stnm)[V(stnm)$Alg == algn[i]]$Shared <- FALSE
}

# Create output file and save the relevant objects
ofname <- paste0(infolder,"-merged.RData")
cat("Output file: ", ofname, "\n")
save(stnm,nruns,num_alg,algn,best,file=ofname)