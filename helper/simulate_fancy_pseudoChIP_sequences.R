# 20170112 Make simulated sequences containing COnensus sequences & PFM matches in random sequeneces

library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TFBSTools)
library(JASPAR2018)

source("/home/ron/fusessh/scripts/plot_themes.R")
source("/home/ron/fusessh/scripts/utility/R_visualization/functions_for_motif_plotting.R")
source("/home/ron/fusessh/scripts/utility/R_genomics_helpers/helper_functions_for_genomics.R")

setwd("/home/ron/fusessh/teaching/tutorial_dl_for_genomics/data/")

set.seed(1234)

# get GATA and CTCF motifs icms plots and pwms --------------------------------
ctcf.list <- getMatrixSet(JASPAR2016, opts=list(species=9606, name="CTCF"))  # get human CTCF PWM
ctcf.icm <- toICM(ctcf.list[[1]]) # information content matrix
ctcf.pfm <- ctcf.list[[1]]@profileMatrix # pfm matrix
ctcf.df <- as.data.frame(ctcf.icm@profileMatrix) #get matrix as dataframe
ctcf.p <- plotICM(ctcf.df)
gata.list <- getMatrixSet(JASPAR2016, opts=list(species=9606, name="GATA2"))  # get human gata PWM
gata.icm <- reverseComplement(toICM(gata.list[[1]])) # information content matrix
gata.pfm <- reverseComplement(gata.list[[1]]@profileMatrix) # pfm matrix
gata.df <- as.data.frame(gata.icm@profileMatrix) #get matrix as dataframe
gata.p <- plotICM(gata.df)
max.list <- getMatrixSet(JASPAR2016, opts=list(species=9606, name="MAX"))  # get human gata PWM
max.icm <- reverseComplement(toICM(max.list[[1]])) # information content matrix
max.pfm <- reverseComplement(max.list[[1]]@profileMatrix) # pfm matrix
max.df <- as.data.frame(max.icm@profileMatrix) #get matrix as dataframe
max.p <- plotICM(max.df)
# even more 
esr.list <- getMatrixSet(JASPAR2016, opts=list(species=9606, name="ESR1"))  # get human gata PWM
esr.icm <- reverseComplement(toICM(esr.list[[1]])) # information content matrix
esr.pfm <- reverseComplement(esr.list[[1]]@profileMatrix) # pfm matrix
esr.df <- as.data.frame(esr.icm@profileMatrix) #get matrix as dataframe
esr.p <- plotICM(esr.df)
sp1.list <- getMatrixSet(JASPAR2016, opts=list(species=9606, name="SP1"))  # get human gata PWM
sp1.icm <- reverseComplement(toICM(sp1.list[[1]])) # information content matrix
sp1.pfm <- reverseComplement(sp1.list[[1]]@profileMatrix) # pfm matrix
sp1.df <- as.data.frame(sp1.icm@profileMatrix) #get matrix as dataframe
sp1.p <- plotICM(sp1.df)
nrf1.list <- getMatrixSet(JASPAR2016, opts=list(species=9606, name="NRF1"))  # get human gata PWM
nrf1.icm <- reverseComplement(toICM(nrf1.list[[1]])) # information content matrix
nrf1.pfm <- reverseComplement(nrf1.list[[1]]@profileMatrix) # pfm matrix
nrf1.df <- as.data.frame(nrf1.icm@profileMatrix) #get matrix as dataframe
nrf1.p <- plotICM(nrf1.df)
sox9.list <- getMatrixSet(JASPAR2016, opts=list(species=9606, name="SOX9"))  # get human gata PWM
sox9.icm <- reverseComplement(toICM(sox9.list[[1]])) # information content matrix
sox9.pfm <- reverseComplement(sox9.list[[1]]@profileMatrix) # pfm matrix
sox9.df <- as.data.frame(sox9.icm@profileMatrix) #get matrix as dataframe
sox9.p <- plotICM(sox9.df)
runx2.list <- getMatrixSet(JASPAR2016, opts=list(species=9606, name="RUNX2"))  # get human gata PWM
runx2.icm <- reverseComplement(toICM(runx2.list[[1]])) # information content matrix
runx2.pfm <- reverseComplement(runx2.list[[1]]@profileMatrix) # pfm matrix
runx2.df <- as.data.frame(runx2.icm@profileMatrix) #get matrix as dataframe
runx2.p <- plotICM(runx2.df)


# save motif plots for showcase


gc.add <- "GCGCGC"
at.add <- "ATATAT"

# generate random sequences ---------------------------------------------------
l <- 300
n <- 20000

# generate set of n times 200 bp sequences
# classes: 
# 0: GC rich nothing else; sp1, nrf1, 
# 1:  gata + max + gc + sox9 + runx2
# 2: ctcf + max + at + sox9 + runx2
# 3: gata + nrf1 + max + sp1
# initialize dataframes with class tagged
c0.set <- data.frame(class=rep(0, n))
c1.pwm.set <- data.frame(class=rep(1, n))
c2.pwm.set <- data.frame(class=rep(2, n))
c3.pwm.set <- data.frame(class=rep(3, n))

# create random sequences with pfms and gc and at matches sampled in
c0.set$seq <- sapply(c0.set$class, function(x) makeRandomSequence(l, add.sequence = as.list(
  rep(gc.add, sample(c(5:20), 1)), add.reverse = T)))

# with pfm samples
c1.pwm.set$seq <- sapply(c1.pwm.set$class, function(x) makeRandomSequence(l, add.sequence = rep(gc.add, sample(c(3:10),1)), add.pfm.match = c(
  rep(list(gata.pfm), sample(c(1:3),1)), 
  rep(list(max.pfm), sample(c(1:3),1)),
  rep(list(sp1.pfm), sample(c(1:3),1)),
  rep(list(nrf1.pfm), sample(c(1:3),1))
), add.reverse = T))
c2.pwm.set$seq <- sapply(c2.pwm.set$class, function(x) makeRandomSequence(l, add.sequence = rep(at.add, sample(c(3:10),1)), add.pfm.match = c(
  rep(list(max.pfm), sample(c(1:3),1)), 
  rep(list(ctcf.pfm), sample(c(1:3),1)),
  rep(list(sox9.pfm), sample(c(1:3),1)),
  rep(list(runx2.pfm), sample(c(1:3),1))
), add.reverse = T))
c3.pwm.set$seq <- sapply(c3.pwm.set$class, function(x) makeRandomSequence(l, add.pfm.match = c(
  rep(list(gata.pfm), sample(c(1:3),1)), 
  rep(list(nrf1.pfm), sample(c(1:3),1)),
  rep(list(max.pfm), sample(c(1:3),1)),
  rep(list(sp1.pfm), sample(c(1:3),1))
), add.reverse = T))


# sample training, test and validation sets -----------------------------------
# assemble sets
# easy.set <- rbind(c0.set, c1.easy.set, c2.easy.set, c3.easy.set)
pwm.set <- rbind(c0.set, c1.pwm.set, c2.pwm.set, c3.pwm.set)
  
# randomly resample
# easy.set <- easy.set[sample(c(1:nrow(easy.set))),]
pwm.set <- pwm.set[sample(c(1:nrow(pwm.set))),]

# write entire data as is
# sample training test and validation
# easy.test.set <- easy.set[c(1:1000),]
# easy.valid.set <- easy.set[c(1001:2000),]
# easy.train.set <- easy.set[-c(1:2000),]
pwm.test.set <- pwm.set[c(1:2000),]
pwm.valid.set <- pwm.set[c(2001:4000),]
pwm.train.set <- pwm.set[-c(1:4000),]

# write data as txt
# write.table(easy.test.set, file = "easy_seq_200bp_test_set.txt", sep="\t", col.names=F, row.names = F, quote=F)
# write.table(easy.valid.set, file = "easy_seq_200bp_valid_set.txt", sep="\t", col.names=F, row.names = F, quote=F)
# write.table(easy.train.set, file = "easy_seq_200bp_train_set.txt", sep="\t", col.names=F, row.names = F, quote=F)
write.table(pwm.test.set, file = "pwm_fancy_seq_200bp_test_set.txt", sep="\t", col.names=F, row.names = F, quote=F)
write.table(pwm.valid.set, file = "pwm_fancy_seq_200bp_valid_set.txt", sep="\t", col.names=F, row.names = F, quote=F)
write.table(pwm.train.set, file = "pwm_fancy_seq_200bp_train_set.txt", sep="\t", col.names=F, row.names = F, quote=F)
