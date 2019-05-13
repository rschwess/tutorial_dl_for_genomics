# R helper to visualize Sequence based NN kernels as Heatmaps and PWMs

source("./helper/functions_for_motif_plotting.R")

args <- commandArgs(trailingOnly = TRUE)

weights.dir <- args[1]  # get weights txt file
out.dir <- args[2]  # store here
plot.prefix <- args[3]  # prefix
plot.width <- as.numeric(args[4])
plot.height <- as.numeric(args[5])

# read in all weights in dir -----------------------------------------------------------------
weights <- list.files(weights.dir, pattern="^filter")

print(weights)

sum.threshold <- 0.01
pass.count <- 0
skip.count <- 0

# Read in every weight txt and make Moitf Plots out of it ------------------------------------
for(weight in weights){
  
  # read
  txt <- read.table(paste0(weights.dir, "/", weight), colClasses = rep("numeric", 4))
  
  #sum up and threshold weights
  su <- sum(txt[txt>=0])
  
  if(su >= sum.threshold){
    pass.count <- pass.count + 1
    # get basename
    basename <- gsub(".txt", ".png", weight)
    # transform to ICM like representation
    # transpose for plot format
    m <- t(txt)
    # name rows
    rownames(m) <- c("A", "C", "G", "T")
    # convert to ICM
    m <- weightsToIcm(m, base = 1000)
    # make plot
    q <- plotWeights(m)
    # save
    print(paste("Saving Plot",basename))
    ggsave(q, filename = paste0(out.dir, "/", plot.prefix, "_", basename), width = plot.width, height = plot.height)  
  }else{
    skip.count <- skip.count + 1
  }
  
}

print(paste('Passing:', pass.count, ' Skipped:', skip.count))





