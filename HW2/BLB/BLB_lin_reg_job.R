
mini <- FALSE
verbose <- FALSE

#============================== Setup for running on Gauss... ==============================#

args <- commandArgs(TRUE)

cat("Command-line arguments:\n")
print(args)

####
# sim_start ==> Lowest possible dataset number
###

###################
sim_start <- 1000
###################

if (length(args)==0){
  sim_num <- sim_start + 1
  set.seed(121231)
} else {
  # SLURM can use either 0- or 1-indexing...
  # Lets use 1-indexing here...
  sim_num <- sim_start + as.numeric(args[1])
  sim_seed <- (762*(sim_num-1) + 121231)
}

cat(paste("\nAnalyzing dataset number ",sim_num,"...\n\n",sep=""))

# Find r and s indices:
s_index<-(sim_num-sim_start-1) %/% 50 + 1
r_index<-(sim_num-sim_start-1) %% 50 + 1

#============================== Run the simulation study ==============================#

# Load packages:
library(BH)
library(bigmemory.sri)
library(bigmemory)
library(biganalytics)

# I/O specifications:
datapath <- "/home/pdbaines/data"
outpath <- "output/"

# mini or full?
if (mini){
	rootfilename <- "blb_lin_reg_mini"
} else {
	rootfilename <- "blb_lin_reg_data"
}

# Filenames:
descriptorfilename <- paste0(rootfilename,".desc")

# Set up I/O stuff:
descriptorfile <- paste(datapath,descriptorfilename,sep="/")

# Attach big.matrix :
if (verbose){
  cat("Attaching big.matrix...\n")
}
mydata <- attach.big.matrix(dget(descriptorfile),backingpath=datapath)

# Remaining BLB specs:
n<-nrow(mydata)
d<-ncol(mydata)-1
gamma<-0.7
b<-n^gamma

# Extract the subset:
set.seed(s_index*17)
samp<-mydata[sample(1:n,b),]

# Reset simulation seed:
set.seed(s_index*r_index)

# Bootstrap dataset:
mywts<-as.numeric(rmultinom(1, n, rep(1,b)/b))

# Fit lm:
myfit<-lm(samp[,d+1]~-1+samp[,1:d],weights=mywts)

# Output file:
outfile = paste0("output/","coef_",sprintf("%02d",s_index),"_",sprintf("%02d",r_index),".txt")

# Save estimates to file:
coef<-myfit$coef
names(coef)<-NULL
write(coef,ncolumns=d,file=outfile)
