library(gridExtra)
library(lattice)
library(dplyr)
createworld <- function(xl = 50,
                        yl = 50,
                        tn = 200,
                        ra_mean = 2,
                        ra_sd = 1,
                        allgroups = 1,
                        groupprob = 1,
                        plot_rank_abundance = FALSE,
                        capacity = 100,
                        distfakt = 0.3,
                        divcapacity = 20,
                        usediv = T) {
  
  # define geographic midepoints of the taxa
  tx <- sample(1:xl,tn,replace = T)
  ty <- sample(1:yl,tn,replace = T)
  
  # maximal abundance of taxa
  tam <- round((rlnorm(tn,ra_mean,ra_sd)))
  tam[which(tam == 0)] <- 1
  
  # assign groups to taxa
  mygroups <- sample(allgroups,size = tn,replace = T,prob = groupprob)
  groupnames1 <- allgroups
  
  # plot the rank-abundance (maxima) of taxa, color-coded for groups
  if(plot_rank_abundance == T) plot(sort(tam, decreasing = T), pch = 19, col = mygroups)
  
  # distribution of taxa:
  
  # array to save the world
  td <- array(0,dim = c(xl,yl,tn))
  
  
  # fill in the distribution into the world - run across the grid filling genera to cells, with occurrence limit (capacity)
  if(usediv == T){ for(xn in 1:xl) {
    for(yn in 1:yl){
      # for every genus a number, depending on max abundance and distance from geogr. midpoint of the genus
      distfunct <- function(abu, n1) {
        dist <- 1+sqrt((abs(tx[n1]-xn))^2 + (abs(ty[n1]-yn))^2)
        return(1/(1+exp(1)^(distfakt*dist-log(abu))))
      }
      probfactors <- sapply(1:tn, function(f) tam[f]*distfunct(abu = tam[f], n1 = f))
      # create ra distribution
      radist <- sort(round(rlnorm(divcapacity,ra_mean,ra_sd)),decreasing = T)
      radist[which(radist==0)] <- 1
      td[xn,yn,sample(1:tn,size = divcapacity,replace = F, prob = probfactors)] <- radist
      
    }
    if((xn/10) %in% 1:1000) print(paste("processing x   ",xn))
    if((xn/10000) == 1) print(paste("Stopped counting at x   ",xn))
  }
  }
  
  # fill in the distribution into the world - run across the grid filling genera to cells, with occurrence limit (capacity)
  if(usediv == F){ for(xn in 1:xl) {
    for(yn in 1:yl){
      # for every genus a number, depending on max abundance and distance from geogr. midpoint of the genus
      distfunct <- function(abu, n1) {
        dist <- 1+sqrt((abs(tx[n1]-xn))^2 + (abs(ty[n1]-yn))^2)
        return(1/(1+exp(1)^(distfakt*dist-log(abu))))
      }
      probfactors <- sapply(1:tn, function(f) tam[f]*distfunct(abu = tam[f], n1 = f))
      taxtab <- table(sample(1:tn,size = capacity, replace = T, prob = probfactors)) # sample the genera with the probfactors
      td[xn,yn,as.numeric(names(taxtab))] <- c(taxtab)
    }
    if((xn/10) %in% 1:1000) print(paste("processing x   ",xn))
    if((xn/10000) == 1) print(paste("Stopped counting at x   ",xn))
  }
  }
  
  
  return(list(td,mygroups,tam))
}


plotworld <- function(world, metric = "occurrences", logarithmic = F) {
  
  thisworld <- world[[1]]
  if(logarithmic == T) thisworld <- log10(world[[1]]+0.1)
  
  if (metric == "occurrences") {
    out <-   levelplot(apply(thisworld[,,],c(1,2),function(x) sum(x)),col.regions = rev(heat.colors(100)), 
                       main = "Occurrences per cell", xlab = "longitude", ylab = "latitude")
    
  }
  
  if (metric == "diversity") {
    out <-   levelplot(apply(thisworld[,,],c(1,2),function(x) length(x[which(x>0)])),col.regions = rev(heat.colors(100)), 
                       main = "Diversity in cells", xlab = "longitude", ylab = "latitude")
  }
  
  if (metric == "evenness") {
    Hmax <- log(apply(thisworld[,,],c(1,2),function(x) length(x[which(x>0)])))
    totaloccs <- apply(world[[1]][,,],c(1,2),function(x) sum(x))
    P <- world[[1]]
    for(pi in 1:dim(world[[1]])[3]) P[,,pi] <- world[[1]][,,pi]/apply(world[[1]][,,],c(1,2),function(x) sum(x))
    H_cont <- -P*log(P)
    H_cont[which(is.na(H_cont))] <- 0
    H <- apply(H_cont,c(1,2),sum)
    Eve <- H/Hmax
    out <-   levelplot(Eve,col.regions = rev(heat.colors(100)), 
                       main = "Evenness in cells", xlab = "longitude", ylab = "latitude")
  }
  
  out
}

plotworldtaxon <- function(world,taxonid, logarithmic = F) {
  thisworld <- world[[1]]
  if(logarithmic == T) thisworld <- log10(world[[1]]+0.1)
  
  
  out <-   levelplot(apply(thisworld[,,taxonid],c(1,2),function(x) sum(x)),col.regions = rev(heat.colors(100)), 
                     main = "Occurrences per cell", xlab = "longitude", ylab = "latitude")
  
  
  out
}


sampleworld <- function(world,
                        mtimebins = 10,      # repeat the sampling process in different "time bins"
                        ncells = rep(100,10),   # number of cells sampled for each time bin. Needs to be
                        # a vector with length mtimebins
                        x_range = NULL,  # optional: restrict the longitude sampled
                        # needs to be a list of vectors of length mtimbins.
                        # Each list element needs to specify the start and the end of the range.
                        # For example: x_range = rep(list(c(25,35)[1])), 10) restricts the longitude
                        # in all 10 time bins to between 25 and 35.
                        y_range = NULL) {
  
  if(is.null(x_range)) x_range <- rep(list(c(1,dim(world[[1]])[1])), mtimebins)
  if(is.null(y_range)) y_range <- rep(list(c(1,dim(world[[1]])[2])), mtimebins)
  
  
  # names of the time bins (just use numbers)
  alltimebins <- as.character(1:mtimebins)
  
  
  # matrix for the data
  occmat <- as.data.frame(matrix(NA,nrow=0,ncol =6))
  
  # fill the matrix
  for(m in 1:mtimebins){
    lonlat <- apply(expand.grid(x_range[[m]][1]:x_range[[m]][2],y_range[[m]][1]:y_range[[m]][2]),1,function(x) paste(x[1],x[2],sep = "_"))
    lonlats <- sample(lonlat,size = ncells[m], replace = F)
    lat1 <- as.numeric(sapply(lonlats,function(x) strsplit(x,split = "_")[[1]][2]))
    lon1 <- as.numeric(sapply(lonlats,function(x) strsplit(x,split = "_")[[1]][1]))
    timebinout <- matrix(NA,nrow=0,ncol =6)
    for(n in 1:ncells[m]){
      celltax <- which(world[[1]][lon1[n],lat1[n],] != 0) # taxon name
      celltaxn <- world[[1]][lon1[n],lat1[n],which(world[[1]][lon1[n],lat1[n],] != 0)] # frequency of these taxa
      groupscell <- world[[2]][celltax]
      out <- cbind(
        unlist(c(mapply(function(x,y) rep(x,each = y),x = groupscell, y = celltaxn))), # group
        unlist(c(mapply(function(x,y) rep(x,each = y),x = celltax, y = celltaxn))), # taxon
        rep(m,sum(celltaxn)), # timebin
        rep(lat1[n],sum(celltaxn)), # latitude
        rep(lon1[n],sum(celltaxn)), #longitude
        rep(lonlats[n],sum(celltaxn)) #cell name
        
      )
      timebinout <- rbind(timebinout,out) 
    }
    
    occmat <- rbind(occmat,timebinout)  
    print(c("time = ", m))
  }
  colnames(occmat) <- c("group", "taxon", "timebin", "latitude", "longitude", "cell")
  # head(occmat)
  
  occmat <- (occmat %>% mutate_if(is.factor,function(x) as.character(x)))
  occmat[c("group", "taxon", "timebin", "latitude", "longitude", "cell"),] <- 
    (occmat[c("group", "taxon", "timebin", "latitude", "longitude", "cell"),] %>% 
       mutate_if(is.character,function(x) as.numeric(x)))
  
  return(occmat)
}


plotsampleworld <- function(data, xl,yl, metric = "occurrences", timebin1 = 1) {
  timesub <- subset(data, timebin == timebin1)
  if (metric == "occurrences") {
    out <- matrix(0,nrow = xl, ncol = yl)
    mf <- sapply(unique(timesub$cell), function(x) {
      c(as.numeric(unique(timesub$longitude[which(timesub$cell == x)])),
        as.numeric(unique(timesub$latitude[which(timesub$cell == x)])),
        nrow(timesub[which(timesub$cell == x),]))
      
    }
    )
    out[cbind(c(mf[1,]),c(mf[2,]))] <- c(mf[3,])
    
    out2 <- levelplot(out,col.regions = rev(heat.colors(100)), main = "Occurrences per cell", xlab = "longitude", ylab = "latitude")
    
  }
  if (metric == "diversity") {
    out <- matrix(0,nrow = xl, ncol = yl)
    mf <- sapply(unique(timesub$cell), function(x) {
      c(as.numeric(unique(timesub$longitude[which(timesub$cell == x)])),
        as.numeric(unique(timesub$latitude[which(timesub$cell == x)])),
        length(unique(timesub$taxon[which(timesub$cell == x)])))
      
    }
    )
    out[cbind(c(mf[1,]),c(mf[2,]))] <- c(mf[3,])
    
    out2 <- levelplot(out,col.regions = rev(heat.colors(100)), main = "Diversity in cells", xlab = "longitude", ylab = "latitude")
    
  }
  out2
}

#
# Define the SCOR function
#
################
####
#
getSCOR <- function(dataframe, timebinnames, groupnames = "none", useallcells = T, celldownsample = F, subcells = NA,
                    taxondownsample = F, ntaxa = NA,
                    samplingrep = 1, n_plus_one = T) { 
  #
  ### define the SCOR function to use in the larger function
  #
  simpleSCOR <- function(df1) {
    # number of occupied cells in the current time bin
    if (useallcells == F) n = length(unique(df1$cell))
    if (useallcells == T) n = timebincells
    # number of cells, for each genus
    gen_locs = apply(tapply(df1$cell,list(df1$taxon,df1$cell), length),
                     1, function(x) length(which(!(is.na(x)))))
    if (n_plus_one == T) p = gen_locs/(n+1) # proportional occupancy for each genus, with n+1 to avoid Inf for omnipresent genera
    if (n_plus_one == F) p = gen_locs/(n) # proportional occupancy for each genus, as in Hannisdal 2012
    
    lambda = -log(1-p) # lambda
    # output: a vector of SCOR, variance of SCOR, standard error of SCOR, number of cells and number of gener
    c(SCOR = sum(lambda), SCORvar = sum(p/((1 - p)*n)), SCORse = sqrt(sum(p/((1 - p)*n)))/sqrt(n),
      Ncells = n, Ngenera = length(unique(df1$taxon)))
  }
  #
  #
  # save the output
  out <- array(NA,dim = c(length(timebinnames), 5, length(groupnames), samplingrep), 
               dimnames = list(timebinnames,c("SCOR","SCORvar","SCORse","Ncells","Ngenera"),groupnames, 1:samplingrep))
  #
  ### Lop for samplingrep sampling repetitions (for cell subsampling)  
  for(s in 1:samplingrep) {
    ### A loop with one repetition for each of the time bins
    for(t in 1:length(timebinnames)){
      #
      #
      # create a subset for the current timebin, and exclude groups that have not been mentioned in groupnames
      tempsub <- subset(dataframe,timebin == timebinnames[t] & group %in% groupnames)
      #
      ### In case you want uniform cell number per time bin: downsample to minimum cell number (ncells):
      if (celldownsample == T)  tempsub <- subset(tempsub,cell %in% base::sample(unique(tempsub$cell), size = subcells, replace = F))
      #
      ### In case you want uniform taxon number per time bin: downsample to minimum taxon number (ntaxon):
      if (taxondownsample == T)  tempsub <- subset(tempsub,taxon %in% base::sample(unique(tempsub$taxon), size = ntaxa, replace = F))
      #
      #
      # Are there groups you want to use?
      if (groupnames[1] != "none"){
        #
        # all cells in the timebin
        timebincells <- length(unique(tempsub$cell))
        # group loop
        for(g in 1:length(groupnames)) {
          groupsub <- subset(tempsub,group == groupnames[g])
          # group scor
          out[t,,g,s] <- simpleSCOR(groupsub)
        } # close the group loop
        #
        # if there are no groups:
      } else {
        timebincells <- length(unique(tempsub$cell))
        out[t,,,s] <- simpleSCOR(tempsub)
      }
    }# cloese the time loop
  } # close the cell subsampling loop
  out
}


plotSCOR <- function(scor) {
  layout(matrix(c(1,2,3,4,5,6),nrow=3), width=c(4,1,4,1,4,1))
  par(mar=c(5,4,2.5,0)) #No margin on the right side
  
  ### SCOR
  plot(0,0,xlim = c(1,length(scor[,1,1,1])), ylim = c(0,max(c(scor[,1,,]))), xlab = "timebin index", ylab = "SCOR",
       main = "SCOR")
  colors2 <- rainbow(ncol(scor[,,1,1]), alpha = 0.8)
  # abline(v = mean(as.numeric(alltimebins)), lty = 2)
  for(i in 1:dim(scor)[3]) points(scor[,1,i,1], type = "o", col = colors2[i], lwd = 3)
  
  ###
  # Genera
  plot(0,0,xlim = c(1,length(scor[,1,1,1])), ylim = c(0,max(c(scor[,5,,]))), xlab = "timebin index", ylab = "genera",
       main = "genera")
  colors2 <- rainbow(ncol(scor[,,1,1]), alpha = 0.8)
  # abline(v = mean(as.numeric(alltimebins)), lty = 2)
  for(i in 1:dim(scor)[3]) points(scor[,5,i,1], type = "o", col = colors2[i], lwd = 3)
  
  ###
  # Cells
  plot(0,0,xlim = c(1,length(scor[,1,1,1])), ylim = c(0,max(c(scor[,4,,]))), xlab = "timebin index", ylab = "cells",
       main = "cells")
  points(scor[,4,1,1], type = "o", col = "black", lwd = 3)
  
  
  # legend
  par(mar=c(5,0.1,4,2)) #No margin on the left side
  plot(c(0,0),type="n", axes=F, xlab="", ylab="") # empty plot for the legend
  if(dim(scor)[3] > 1) legend("center", colnames(scor[,1,,]), col = colors2, lwd = 3, cex = 0.8) 
  #
  plot(c(0,0),type="n", axes=F, xlab="", ylab="") # empty plot for the legend
  if(dim(scor)[3] > 1) legend("center", colnames(scor[,1,,]), col = colors2, lwd = 3, cex = 0.8) 
}


repeatsampleSCOR <- function(
  #
  world,
  repetitions = 10,
  # for the resampling
  mtimebins = 20,
  ncells = rep(100,20),
  x_range = rep(list(c(1,xl)), 20),
  y_range = rep(list(c(1,yl)), 20),
  # for SCOR
  timebinnames, groupnames = "none", useallcells = T, celldownsample = F, subcells = NA,
  taxondownsample = F, ntaxa = NA,
  samplingrep = 1, n_plus_one = T
  
) {
  scor_out <- array(NA,dim = c(length(timebinnames), 5, length(groupnames), samplingrep, repetitions), 
                    dimnames = list(timebinnames,c("SCOR","SCORvar","SCORse","Ncells","Ngenera"),groupnames, 1:samplingrep,
                                    1:repetitions))  
  for (nrep in 1:repetitions){
    sworld <- sampleworld(world = world, mtimebins = mtimebins, ncells = ncells, x_range = x_range, y_range = y_range)
    
    scor_out[,,,,nrep] <-  getSCOR(sworld, timebinnames = timebinnames, groupnames = groupnames, useallcells = useallcells, 
                                   celldownsample = celldownsample, subcells = subcells, 
                                   taxondownsample = taxondownsample, ntaxa = ntaxa, 
                                   samplingrep = samplingrep, n_plus_one = n_plus_one)
    print(paste("this is repetition   ", nrep))
  }
  return(scor_out)
}

error_polygon <- function(ep,en,tstart,tend,tmid,color) {
  
  polygon( c(tstart, tmid, tend, tend, rev(tmid), tstart),
           c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
           border = NA, col = color)
  
}

plotrepeatSCOR <- function(scor, endext = 0.2, metric = "se", nmet = 2) {
  layout(matrix(c(1,2,3,4,5,6),nrow=3), width=c(4,1,4,1,4,1))
  par(mar=c(5,4,2,0)) #No margin on the right side
  
  ### SCOR
  
  sc1 <-  apply(scor,c(1,2,3,4),mean)
  if(metric == "se") scsd <- apply(scor,c(1,2,3,4),function(x) sd(x) /sqrt(length(test3[1,1,1,1,])))
  if(metric == "sd") scsd <- apply(scor,c(1,2,3,4),function(x) sd(x))
  
  colors2 <- rainbow(length(sc1[1,1,,1]), alpha = 0.8)
  colors3 <- rainbow(length(sc1[1,1,,1]), alpha = 0.2)
  
  plot(0,0,xlim = c(1,length(sc1[,1,1,1])), ylim = c(0,max(c(scor[,1,,,]))), xlab = "timebin index", 
       ylab = "SCOR")
  for(i in 1:length(scor[1,1,,1,1])) error_polygon(ep = sc1[,1,i,1]+nmet*scsd[,1,i,1], en = sc1[,1,i,1]-nmet*scsd[,1,i,1],
                                                   tstart = 1-endext,tend=endext+length(sc1[,1,i,1]),
                                                   tmid=1:length(sc1[,1,i,1]),  color = colors3[i])
  for(i in 1:length(scor[1,1,,1,1])) points(sc1[,1,i,1], type = "o", col = colors2[i], lwd = 3)
  
  ###
  # Genera
  plot(0,0,xlim = c(1,length(sc1[,5,1,1])), ylim = c(0,max(c(scor[,5,,,]))), xlab = "timebin index", 
       ylab = "genera")
  
  for(i in 1:length(scor[1,1,,1,1])) error_polygon(ep = sc1[,5,i,1]+nmet*scsd[,5,i,1], en = sc1[,5,i,1]-nmet*scsd[,5,i,1],
                                                   tstart = 1-endext,tend=endext+length(sc1[,5,i,1]),
                                                   tmid=1:length(sc1[,5,i,1]),  color = colors3[i])
  for(i in 1:length(scor[1,1,,1,1])) points(sc1[,5,i,1], type = "o", col = colors2[i], lwd = 3)
  
  ###
  # Cells
  plot(0,0,xlim = c(1,length(sc1[,4,1,1])), ylim = c(0,max(c(scor[,4,,,]))), xlab = "timebin index", 
       ylab = "cells")
  
  error_polygon(ep = sc1[,4,1,1]+nmet*scsd[,4,1,1], en = sc1[,4,1,1]-nmet*scsd[,4,1,1],
                tstart = 1-endext,tend=endext+length(sc1[,4,1,1]),
                tmid=1:length(sc1[,4,1,1]),  color = colors3[i])
  points(sc1[,4,1,1], type = "o", col = "black", lwd = 3)
  
  
  # legend
  par(mar=c(5,0.1,4,2)) #No margin on the left side
  plot(c(0,0),type="n", axes=F, xlab="", ylab="") # empty plot for the legend
  if(dim(scor)[3] > 1) legend("center", names(scor[1,1,,1,1]), col = colors2, lwd = 3, cex = 0.8) 
  #
  plot(c(0,0),type="n", axes=F, xlab="", ylab="") # empty plot for the legend
  if(dim(scor)[3] > 1) legend("center", names(scor[1,1,,1,1]), col = colors2, lwd = 3, cex = 0.8) 
}
