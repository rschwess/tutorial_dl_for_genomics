# R functions for plotting Motifs, ICM matrices with ggplot
# Author: Ron Schwessinger
# Date: 13.12.2016

# Plotting ICM matrices using ggplot2 to combine those plots with other ggplots and cowplot etc.
# The functions implement the plotting of the ICM/Motif by having the letters as polygons
# Pretty much what you see if you download an svg from JASPAR
# Usage:
# plotICM(data.frame) to plot the ICM matrix (4 rows A, C, G, T times the number of positions)
#
# Note:
# all alpha like usable so happy about any comments
# *circle function adapted from Joran from stackoverflow
# * don't tell Hadley ;)

# required
require(ggplot2, quietly = TRUE)
require(RColorBrewer, quietly = TRUE)

# THEME -----------------------------------------------------------------------
motif_theme <- theme(
  panel.grid = element_blank(),
  text = element_text(size = 14),
  axis.title.x = element_blank(),
  axis.line = element_line(color="black", size = 0.7),
  axis.line.x = element_line(color="black", size = 0.7),
  axis.line.y = element_line(color="black", size = 0.7),
  plot.margin = unit(c(0.7,0.7,0.7,0.7), "lines"),
  panel.border=element_blank(),
  strip.background = element_blank()
)

# theme for aligned matching sequence plotting
seq_theme <- theme(
  panel.grid = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_text(angle=0, hjust=.6),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  plot.margin = unit(c(0,-2, 0, -2), "lines"),
  panel.border=element_blank(),
  strip.background = element_blank(),
  legend.position = "none"
)

letter.colours <- brewer.pal(5, "Set1")[c(3,2,5,1,4)]

# FUNCTIONS -------------------------------------------------------------------
circleFun <- function(center = c(0,0), diameter = 1, npoints = 100){
  # circle plotting helper function modified from "joran" (stackoverflow)
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  xx[xx < .5] <- xx[xx < .5] + (0.5 - xx[xx < .5]) * .2 #squash a bit
  xx[xx > .5] <- xx[xx > .5] - (xx[xx > .5] - .5) * .2
  return(data.frame(x = xx, y = yy))
}

plotBitBase <- function(letter, scale=1, xintercept=0, yintercept=0){
  # Function to plot a Base letter in Bit content format using ggplot2 and polygons
  # Input:
  #   letter: letter [A,C,G,T]
  #   scale: REAL scale factor to multiply with the coordinates (defaults to 1). 
  #           Use to scale the letter size to the bit content (will only scale in y direction to keep column width)
  #   xintercept: REAL adjust the x-postion (center) of the base letter (default=0) use to place the letter in your plot
  #   yintercept: REAL adjust the y-postion (default=0,) use to stack letters
  # Returns: list of geom_polygon(s) that can be added to any ggplot2 bject/plot
  
  # check letter
  if(!letter %in% c("A", "C", "G", "T")){
    warning(paste0(letter, " is not a valid base letter!"))
    return(NA_real_)
  }
  
  # A
  if(letter == "A"){
    a.pol1 <- data.frame(
      x=c(.1, .25, .35, .65, .75, .9, .6, .4),
      y=c(0, 0, .35, .35, 0, 0, 1, 1)
    )
    a.pol2 <- data.frame(
      x=c(.4, .6, .5),
      y=c(.5, .5, .85)
    )
    # scale
    a.pol1$x <- a.pol1$x + xintercept
    a.pol1$y <- a.pol1$y * scale + yintercept
    a.pol2$x <- a.pol2$x + xintercept
    a.pol2$y <- a.pol2$y * scale + yintercept
    
    l <- list(
      geom_polygon(data=a.pol1, aes(x=x, y=), fill=letter.colours[1]),
      geom_polygon(data=a.pol2, aes(x=x, y=), fill="white")
    )
    
    # C
  }else if(letter == "C"){
    c.pol1 <- circleFun(c(.5,.5), 1, npoints = 100)
    c.pol2 <- circleFun(c(.5,.5), .75, npoints =100)
    c.pol3 <- data.frame(x=c(.75, 1, 1, .75), y=c(.3, .3, .7, .7))
    # scale
    c.pol1$x <- c.pol1$x + xintercept
    c.pol2$x <- c.pol2$x + xintercept
    c.pol3$x <- c.pol3$x + xintercept
    c.pol1$y <- c.pol1$y * scale + yintercept
    c.pol2$y <- c.pol2$y * scale + yintercept
    c.pol3$y <- c.pol3$y * scale + yintercept
    # report list
    l <- list(
      geom_polygon(data=c.pol1, aes(x=x, y=y), fill=letter.colours[2]),
      geom_polygon(data=c.pol2, aes(x=x, y=y), fill="white"),
      geom_polygon(data=c.pol3, aes(x=x, y=y), fill="white")
    )
    # G
  }else if(letter == "G"){
    g.pol1 <- circleFun(c(.5,.5), 1, npoints = 100)
    g.pol2 <- circleFun(c(.5,.5), .75, npoints =100)
    g.pol3 <- data.frame(x=c(.75, 1, 1, .75), y=c(.3, .3, .7, .7))
    g.pol4 <- data.frame(x=c(.75, .9, .9, .5, .5, .75), y=c(0, 0, .45, .45, .35, .35))
    # scale
    g.pol1$x <- g.pol1$x + xintercept
    g.pol2$x <- g.pol2$x + xintercept
    g.pol3$x <- g.pol3$x + xintercept
    g.pol4$x <- g.pol4$x + xintercept
    g.pol1$y <- g.pol1$y * scale + yintercept
    g.pol2$y <- g.pol2$y * scale + yintercept
    g.pol3$y <- g.pol3$y * scale + yintercept
    g.pol4$y <- g.pol4$y * scale + yintercept
    #report list
    l <- list(
      geom_polygon(data=g.pol1, aes(x, y), fill=letter.colours[3]),
      geom_polygon(data=g.pol2, aes(x=x, y=y), fill="white"),
      geom_polygon(data=g.pol3, aes(x=x, y=y), fill="white"),
      geom_polygon(data=g.pol4, aes(x=x, y=y), fill=letter.colours[3])
    )
    # T
  }else if(letter == "T"){
    t.pol <- data.frame(
      x=c(.425, .575, .575, .9, .9, .1, .1, .425),
      y=c(0, 0, .85, .85, 1, 1, .85, .85)
    )
    # scale
    t.pol$x <- t.pol$x + xintercept
    t.pol$y <- t.pol$y * scale + yintercept
    # report list
    l <- list(geom_polygon(data=t.pol, aes(x=x, y=y), fill=letter.colours[4]))
  }
  
  return(l)
  
}

plotICM <- function(icm){
  # Plot function to generate a polygon based base letter representation of the ICM matrix/motif
  # Input:
  #   icm: 4 * X dataframe (rows: A, C, G, T, columns: as many positions as there are in the motif)
  # Returns: ggplot2 object with the motif plot
  
  # space for icm data frame checking
  if(nrow(icm) != 4){
    warnings("Must be a 4 row (A,C,G,T) data frame")
    return(NA_character_)
  }
  
  # convert data frame to data frame listing the bit content of letters stagged on top of each other with the largst on top
  stagged.icm <- apply(icm, 2, function(x){
    temp <- data.frame(from=rep(0,4), to=rep(0,4)) # init from_to table
    row.names(temp) <- c("A", "C", "G", "T")
    ordered <- order(x) # order entries
    cum <- 0 # init cumulative value
    # run over bases and count up cum and set from to values for plot
    for(i in c(1:4)){
      temp[ordered[i], "from"] <- cum
      cum <- cum + x[ordered[i]]
      temp[ordered[i], "to"] <- cum
    }
    temp <- unlist(temp)
    return(temp)
  })
  # Convert those to a from - to plot valued, long dataframe for ggplot2
  stagged.icm <- t(stagged.icm)
  temp.df <- data.frame(
    pos=rep(c(1:ncol(icm)), times=4),
    base=rep(c("A", "C", "G", "T"), each=ncol(icm)),
    from=c(stagged.icm[,"from1"], stagged.icm[,"from2"], stagged.icm[,"from3"], stagged.icm[,"from4"]),
    to=c(stagged.icm[,"to1"], stagged.icm[,"to2"], stagged.icm[,"to3"], stagged.icm[,"to4"])
  ) # convert to dataframe
  
  # lay plot base
  p <- ggplot(data.frame(x=factor(c(0,ncol(icm)), levels=c(0,ncol(icm))), y=c(0, max(temp.df$to))), aes(x=x, y=y)) + 
    labs(x="pos", y="bits") + xlim(.5, nrow(temp.df)/4+.5) + ylim(0,2) + theme_bw() + motif_theme 
  
  # add base letters using the position as xintercept, the bit content as scale 
  # and the summed up bit content of the lower letters as yintercept
  for(i in c(1:nrow(temp.df))){
    p <- p + plotBitBase(temp.df$base[i], xintercept = temp.df$pos[i]-1+.5, yintercept = temp.df$from[i], scale= temp.df$to[i] - temp.df$from[i])  
  }
  return(p)
}

plotAlignedSeq <- function(seq, id="", with.variant=FALSE, variant.base="N", variant.pos=0){
  # Function to plot a matching/aligned sequence (with variant) to combine/plot under a sequence motif
  # Input:
  #   seq: character string - sequence input
  #   id: CHAR string to label the sequence/variant
  #   with.variant: TRUE/FALSE - select if to plot with or without a vairant base
  #   variant.base: variant base [A,C,G,T,X] use "X" for deletions
  #   variant.pos: INT position of the variant base
  # Returns: ggplot2 object of plotted/aligned sequence (e.g. use with cowplot::plot_grid for aligned combination of plots)
  
  # make sequence df
  temp.df <- data.frame(
    pos=c(1:nchar(seq)),
    base=factor(unlist(strsplit(seq, "")), levels=c("A", "C", "G", "T")),
    height=rep(1, nchar(seq))
  )
  
  if(with.variant){
    # make  SNP df
    temp.snp.df <- data.frame(
      pos = variant.pos,
      base = variant.base,
      height= .85
    )
  }
  # make plot
  q <- ggplot(temp.df, aes(y=height, x=pos, label=base, fill=base)) + 
    geom_label(size=6) + 
    ylab(id) + 
    seq_theme + 
    xlim(.5, nchar(seq)+.5) +
    scale_fill_manual(values = letter.colours)  
  
  # add variant if specified
  if(with.variant){
    q <- q + geom_label(data=temp.snp.df, aes(y=height, x=pos, label=base), size=5.5) + 
      ylim(.75, 1.1)
  }else{
    q <- q + ylim(.9, 1.1) #set no variant ranges
  }
  return(q)
}

# HELPER FUNCTIONS ------------------------------------------------------------
reverseBase <- function(b){
  r <- b
  # helper function to reverser base
  if(b == "A"){
    r <- "T"
  }else if(b == "C"){
    r <- "G" 
  }else if(b == "G"){
    r <- "C" 
  }else if(b == "T"){
    r <- "A"
  }
  return(r)
}

# SANDBOX FUNCTIONS -----------------------------------------------------------
plotMotifWithVarMatches <- function(chr, pos, id, ref.base, alt.base, extend, pwm, icm, bsgenome, min.score){
  # Wrapper funtion: make a cowplo t combined version for single SNPs of interest
  # Requires Biostrings, the respective bsgenome, TFBSTools, ggplot2, cowplot, RColorBrewer
  # Input:
  #   chr: CHAR chromosome
  #   pos: INT  1-based position in genome
  #   id: CHAR rsID or similar to label the variant
  #   ref.base: CHAR single base reference
  #   alt.base: CHAR single base alternative base
  #   extend: INT number of base pairs around variant to extend the sequence and scan pwm against
  #   pwm: TFBSTools PWM object
  #   icm: TFBSTool ICM object
  #   bsgenome: biostrings BSGenome to retrieve the reference sequence from
  #   min.score: REAL/PERCETNAGE minimum relative score a PWM match has to achieve to be reported
  # Returns: cowplot combining the motif plot and aligned sequnces with variant plots
  
  require(Biostrings)
  require(TFBSTools)
  require(ggplot2)
  require(RColorBrewer)
  require(cowplot)
  
  # placeholder to check input
  
  # init empty vector of potential var/alt positions 
  position.var <- c()
  
  # 0) make the icm plot
  icm.plot <- plotICM(icm.matrix) # plot icm as part of later plots
  
  # 1) get longer surrounding sequence
  ref.seq <- as.character(getSeq(bsgenome, as.character(chr), start=as.numeric(pos)-extend, end=as.numeric(pos)+extend))
  # 2) make alternative sequence with alt.base
  alt.seq <- ref.seq
  temp.base <- substring(alt.seq, extend+1, extend+1)  # check if middle equals ref base
  if(temp.base != ref.base){
    print(paste0("Warning ... ", id, "'s major ref base not in reference genome"))
    if(temp.base == alt.base){
      print("Reference base matches the alt.base ... substituting ...")
      substring(alt.seq, extend+1, extend+1) <- as.character(ref.base)
      alt.base <- ref.base
      ref.base <- temp.base
    }else{
      return(NA_character_)
    }
  }else{
    substring(alt.seq, extend+1, extend+1) <- as.character(alt.base)  # substitute for alternative base
  }
  
  # 3) match sequence against motif
  motif.match <- searchSeq(pwm, ref.seq, seqname="seq1", min.score=min.score, strand="*")
  # for each match
  num.matches <- length(motif.match@views) # get number of matches
  seq.plot.list <- list() # init list
  seq.plot.list[[1]] <- icm.plot
  
  # 4) make sequence aligned plot parts for every matching sequence
  if(num.matches == 0){
    # what to do when no matches found
    warnings("No PWM matches to sequence above threshold!")
    return(NA_character_)
  }else{
    for(i in c(1:num.matches)){
      # selected  scoring matches top score down
      j <- order(motif.match@score, decreasing=TRUE)[i]
      
      match.offset <- min(motif.match@views@ranges[j])  # get offset
      match.strand <- motif.match@strand[j]  # get strand
      temp.var <- alt.base
      if(match.strand == "-"){
        match.seq <- as.character(reverseComplement(motif.match@views[[j]]))
        temp.var <- reverseBase(temp.var)
        match.offset <- max(motif.match@views@ranges[j])  # get offset
        match.offset <- match.offset - extend
      }else{
        match.seq <- as.character(motif.match@views[[j]])  # get sequence part according to strand info 
        match.offset <- min(motif.match@views@ranges[j])  # get offset
        match.offset <- extend + 1 - match.offset + 1
      }
      # check if variant is within the matching sequence
      variant.in.match <- TRUE
      if((match.offset > extend + 1) | (match.offset + ncol(factor.icm) < extend + 1)){
        variant.in.match <- FALSE
      }
      match.seq.plot <- plotAlignedSeq(
        seq=match.seq, 
        id=paste0(id, ".",i," ", match.strand), 
        with.variant = variant.in.match, 
        variant.base = temp.var, 
        variant.pos = match.offset)
      seq.plot.list[[i+1]] <- match.seq.plot
      
      position.var <- c(position.var, match.offset)  # save potential var position
    }
  }
  
  # make the combined plot
  combined.plot <- do.call(
    plot_grid, c(
      seq.plot.list, list(
        nrow=length(seq.plot.list),
        align="v",
        rel_heights=c(2.5, rep(1, length(seq.plot.list)-1))
        )
      )
    )
  
  # report plot and relative variant positions
  return(list(plot=combined.plot, position.var=position.var))
  
}

plotWeights <- function(w){
  # Plot function to generate a polygon based base letter representation of a kernel weight matrix of shape 4(bases)*X
  # Input:
  #   w: 4 * X dataframe (rows: A, C, G, T, columns: as many positions as there are in the kernel)
  # Returns: ggplot2 object with the motif plot

  # space for icm data frame checking
  if(nrow(w) != 4){
    warnings("Must be a 4 row (A,C,G,T) data frame")
    return(NA_character_)
  }

  # convert data frame to data frame listing the weight content of letters stagged on top or below one another
  stagged.w <- apply(w, 2, function(x){
    temp <- data.frame(from=rep(0,4), to=rep(0,4)) # init from_to table
    row.names(temp) <- c("A", "C", "G", "T")
    ordered <- order(x) # order entries
    pos.cum <- 0 # init cumulative value
    neg.cum <- 0 # init cumulative value
    # run over bases and count up cum and set from to values for plot
    for(i in c(1:4)){
      if(x[ordered[i]] >= 0){
        temp[ordered[i], "from"] <- pos.cum
        pos.cum <- pos.cum + x[ordered[i]]
        temp[ordered[i], "to"] <- pos.cum  
      }else if(x[ordered[i]] < 0){
        temp[ordered[i], "to"] <- neg.cum
        neg.cum <- neg.cum + x[ordered[i]]
        temp[ordered[i], "from"] <- neg.cum  
      }
    }
    temp <- unlist(temp)
    return(temp)
  })
  # Convert those to a from - to plot valued, long dataframe for ggplot2
  stagged.w <- t(stagged.w)
  temp.df <- data.frame(
    pos=rep(c(1:ncol(w)), times=4),
    base=rep(c("A", "C", "G", "T"), each=ncol(w)),
    from=c(stagged.w[,"from1"], stagged.w[,"from2"], stagged.w[,"from3"], stagged.w[,"from4"]),
    to=c(stagged.w[,"to1"], stagged.w[,"to2"], stagged.w[,"to3"], stagged.w[,"to4"])
  ) # convert to dataframe

  # lay plot base
  p <- ggplot(data.frame(x=factor(c(0,ncol(w)), levels=c(0,ncol(w))), y=c(0, max(temp.df$to))), aes(x=x, y=y)) +
    geom_hline(yintercept=0) +
    labs(x="pos", y="weight") + xlim(.5, nrow(temp.df)/4+.5) + #ylim(-1,1) + 
    theme_bw() + motif_theme

  # add base letters using the position as xintercept, the bit content as scale
  # and the summed up bit content of the lower letters as yintercept
  for(i in c(1:nrow(temp.df))){
    p <- p + plotBitBase(temp.df$base[i], xintercept = temp.df$pos[i]-1+.5, yintercept = temp.df$from[i], scale= temp.df$to[i] - temp.df$from[i])
  }
  return(p)
}
