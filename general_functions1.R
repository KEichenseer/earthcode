
# Not in
'%!in%' <- function(x,y)!('%in%'(x,y))

# logit
logit <- function(lo) log(lo/(1-lo)) 
logit2 <- function(lo) {
  res <- lo # a vector to save updated results
  closest_to_boundaries <- min(min(res[which(res>0)]), 1-max(res[which(res<1)])) # the smallest distance from 0 or 1
  res[which(res == 0)] <- 0.5*closest_to_boundaries # replace 0
  res[which(res == 1)] <- 1 - 0.5*closest_to_boundaries # replace 1
  
  log(res/(1-res)) # calculate the logit with the new data
}

# transparency function
#note: always pass alpha on the 0-255 scale
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

# custom function for the mode

modus <- function(x, na.rm = TRUE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  
  tab <- tabulate(match(x, ux)); ux[tab == max(tab)]
  
  return(ux[tab == max(tab)])
  
}

log_mod <- function(lo) {
  res <- lo # a vector to save updated results
  closest_to_zero <- min(res[which(res>0)]) # the smallest distance from 0 or 1
  res[which(res == 0)] <- 0.5*closest_to_zero # replace 0
  log(res)
}


gend <- function(x) {
  saveacf <- acf(x,plot = F)
  savegendif <- NULL
  for (i in 1:(length(x)-1)) savegendif[i] <- x[i+1] - saveacf$acf[2]*x[i]
  savegendif
}


minmax <- function(x) c(min(x,na.rm = T),max(x, na.rm = T))


### Three-point-average

threepav <- function(x) {
  xl <- length(x)
  save3x <- NULL
  for(i in 2:(xl-1)) {
    save3x[i-1] <- mean(x[c(i-1,i,i+1)],na.rm = T)
  }
  save3x
}


### Section correlations

seccor <- function(a,b) {
  output <- matrix(NA,nrow = 7, ncol = 2)
  output[1,] <- unlist(cor.test(a[1:45],b[1:45], method = "spearman")[c(3,4)])
  output[2,] <- unlist(cor.test(a[1:21],b[1:21], method = "spearman")[c(3,4)])
  output[3,] <- unlist(cor.test(a[22:45],b[22:45], method = "spearman")[c(3,4)])
  output[4,] <- unlist(cor.test(a[1:17],b[1:17], method = "spearman")[c(3,4)])
  output[5,] <- unlist(cor.test(a[18:21],b[18:21], method = "spearman")[c(3,4)])
  output[6,] <- unlist(cor.test(a[22:29],b[22:29], method = "spearman")[c(3,4)])
  output[7,] <- unlist(cor.test(a[30:45],b[30:45], method = "spearman")[c(3,4)])
  output
  
}

gendseccor <- function(a,b) {
  output <- matrix(NA,nrow = 7, ncol = 2)
  output[1,] <- unlist(cor.test(gend(a[1:45]),gend(b[1:45]), method = "spearman")[c(3,4)])
  output[2,] <- unlist(cor.test(gend(a[1:21]),gend(b[1:21]), method = "spearman")[c(3,4)])
  output[3,] <- unlist(cor.test(gend(a[22:45]),gend(b[22:45]), method = "spearman")[c(3,4)])
  output[4,] <- unlist(cor.test(gend(a[1:17]),gend(b[1:17]), method = "spearman")[c(3,4)])
  output[5,] <- unlist(cor.test(gend(a[18:21]),gend(b[18:21]), method = "spearman")[c(3,4)])
  output[6,] <- unlist(cor.test(gend(a[22:29]),gend(b[22:29]), method = "spearman")[c(3,4)])
  output[7,] <- unlist(cor.test(gend(a[30:45]),gend(b[30:45]), method = "spearman")[c(3,4)])
  output
  
}



stage_seccor <- seccor <- function(a,b) {
  output <- matrix(NA,nrow = 10, ncol = 2)
  output[1,] <- unlist(cor.test(a[1:85],b[1:85], method = "spearman")[c(3,4)])
  output[2,] <- unlist(cor.test(a[1:38],b[1:38], method = "spearman")[c(3,4)])
  output[3,] <- unlist(cor.test(a[39:85],b[39:85], method = "spearman")[c(3,4)])
  output[4,] <- unlist(cor.test(a[1:29],b[1:29], method = "spearman")[c(3,4)])
  output[5,] <- unlist(cor.test(a[30:38],b[30:38], method = "spearman")[c(3,4)])
  output[6,] <- unlist(cor.test(a[39:52],b[39:52], method = "spearman")[c(3,4)])
  output[7,] <- unlist(cor.test(a[53:68],b[53:68], method = "spearman")[c(3,4)])
  output[8,] <- unlist(cor.test(a[68:85],b[68:85], method = "spearman")[c(3,4)])
  output[9,] <- unlist(cor.test(a[39:67],b[39:67], method = "spearman")[c(3,4)])
  output[10,] <- unlist(cor.test(a[53:85],b[53:85], method = "spearman")[c(3,4)])
  output
  
}

stage_gendseccor <- function(a,b) {
  output <- matrix(NA,nrow = 10, ncol = 2)
  output[1,] <- unlist(cor.test(gend(a[1:85]),gend(b[1:85]), method = "spearman")[c(3,4)])
  output[2,] <- unlist(cor.test(gend(a[1:38]),gend(b[1:38]), method = "spearman")[c(3,4)])
  output[3,] <- unlist(cor.test(gend(a[39:85]),gend(b[39:85]), method = "spearman")[c(3,4)])
  output[4,] <- unlist(cor.test(gend(a[1:29]),gend(b[1:29]), method = "spearman")[c(3,4)])
  output[5,] <- unlist(cor.test(gend(a[30:38]),gend(b[30:38]), method = "spearman")[c(3,4)])
  output[6,] <- unlist(cor.test(gend(a[39:52]),gend(b[39:52]), method = "spearman")[c(3,4)])
  output[7,] <- unlist(cor.test(gend(a[53:67]),gend(b[53:67]), method = "spearman")[c(3,4)])
  output[8,] <- unlist(cor.test(gend(a[68:85]),gend(b[68:85]), method = "spearman")[c(3,4)])
  output[9,] <- unlist(cor.test(gend(a[39:67]),gend(b[39:67]), method = "spearman")[c(3,4)])
  output[10,] <- unlist(cor.test(gend(a[53:85]),gend(b[53:85]), method = "spearman")[c(3,4)])
  output
  
}
# proportional standard error

pse <- function(p,n) sqrt((p*(1-p))/n)

# a function to evenlyt extend the range of two points
fairextend <- function(m,p) {
  output <- NULL
  output[1] <- m[1] - (m[2]-m[1])*p
  output[2] <- m[2] + (m[2]-m[1])*p
  output
}
fairextend(c(2,10),0.1)

# a function to de-mean and stitch together two time series between t and (t-1)
stitch <- function(a,t) {
  a1 <- a[1:(t-1)]
  a2 <- a[t:length(a)]
  
  b1 <- a1 - mean(a1)
  b2 <- a2 - mean(a2)
  
  c(b1,b2)
}


# a function to linearly detrend a time series
ldt <- function(x,t) {
  lmx <- lm(x~t)
  model_formula <- function(t) lmx$coefficients[1] + lmx$coefficients[2]*t
  x - model_formula(t)
}


# functoin to draw polygons of dataframes, proportions from 0 to 1
polygon( c(tstart, tmid, tend, tend, rev(tmid), tstart),
         c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
         border = NA, col = color)
z <- cbind(c(0.1,0.2,0.1,0.2,0.5),c(0.1,0.1,0.1,0.0,0.2),c(0.5,0,0.2,0.1,0.2),c(0.3,0.7,0.6,0.7,0.1))
z
col <- c("blue","red","yellow","green")
#apply(z[,1:i],1,sum)
t <- 1:5
tstart <- 0
tend <- 6
polyplot01 <- function(z,t,tstart,tend,col) {
  nc <- ncol(z)
  nr <- nrow(z)
  plot(t,rep(0,nr),type = "n",ylim = c(0,1), xlim = c(tstart,tend), xaxs = "i", yaxs = "i")
  for(i in 1:nc) {
    lpos <- apply(as.matrix(z[,1:i]),1,sum) - z[,i]
    upos <- apply(as.matrix(z[,1:i]),1,sum)
    polygon(c(tstart,t,tend,tend,rev(t),tstart), 
            c(lpos[1],lpos,lpos[nr],upos[nr],rev(upos),upos[1]),
            border = NA, col = col[i])
  }
}

#polyplot01(z,t,tstart,tend,col)


# function  to estimate missing b from a, based on log-log relationship

logestimator <- function(a,b) {
  llma <- lm(log(a) ~ log(b))
  model_formulaa <- function(b) llma$coefficients[1] + llma$coefficients[2]*log(b)
  
  llmb <- lm(log(b) ~ log(a))
  model_formulab <- function(a) llmb$coefficients[1] + llmb$coefficients[2]*log(a)
  
  a_estimated <-  model_formulaa(b[which(is.na(a))])
  a_out <- a
  a_out[which(is.na(a))] <- exp(1)^a_estimated
  
  b_estimated <-  model_formulab(a[which(is.na(b))])
  b_out <- b
  b_out[which(is.na(b))] <- exp(1)^b_estimated
  
  cbind(a_out,b_out)
  
  
}
rawestimator <- function(a,b) {
  llma <- lm(a ~ b)
  model_formulaa <- function(b) llma$coefficients[1] + llma$coefficients[2]*b
  
  
  llmb <- lm(b ~ a)
  model_formulab <- function(a) llmb$coefficients[1] + llmb$coefficients[2]*a
  
  
  a_estimated <-  model_formulaa(b[which(is.na(a))])
  a_out <- a
  a_out[which(is.na(a))] <- a_estimated
  
  b_estimated <-  model_formulab(a[which(is.na(b))])
  b_out <- b
  b_out[which(is.na(b))] <- b_estimated
  
  cbind(a_out,b_out)
  
}

error_polygon <- function(ep,en,tstart,tend,tmid,color) {
  
  polygon( c(tstart, tmid, tend, tend, rev(tmid), tstart),
           c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
           border = NA, col = color)
  
}

error_polygon2 <- function(ep,en,tstart,tend,tmid,color) {
  nacheck <- is.na(ep)
  segm <- rle(nacheck)
  starts <- c(1,cumsum(segm$lengths[1:(length(segm$lengths)-1)])+1)
  ends <- cumsum(segm$lengths)
  for(i in 1:length(segm$values)) {
    if(!(is.na(segm$values[i]))) {
      ind <- starts[i]:ends[i]
      tmidi <- tmid[ind]
      tstarti <- tstart[ind]
      tendi <- tend[ind]
      
      polygon( c(tstarti[1], tmidi, tendi[length(tendi)], tendi[length(tendi)], rev(tmidi), tstarti[1]),
               c((ep[ind])[1],ep[ind], (ep[ind])[length(ep[ind])], (en[ind])[length(en[ind])], rev(en[ind]), (en[ind])[1]),
               border = NA, col = color)
      
    }
  }
}

vo <- c(2,4,8,16,32)
to <- c(1,2,3,4,6)
tt <- c(0,1,4,4,4.5,8)
l <- 5
#### Intermediate values at certain points
inter <- function(vo,to,tt) {
  vn <- rep(NA,length(tt))
  for (l in 1:length(tt)) {
    if (tt[l] <= min(to)) vn[l] <- vo[which(to == min(to))]
    if (tt[l] > min(to) & tt[l] < max(to)) {
      if(tt[l] %in% to) {vn[l] <- vo[which(to == tt[l])]
      } else {
        lower <- which(to < tt[l])
        higher <- which(to > tt[l])
        l1 <- to[lower[length(lower)]]
        h1 <- to[higher[1]]
        dur <- h1  - l1
        pos <- (h1-tt[l])/dur
        vn[l] <- (pos*vo[lower[length(lower)]]+((1-pos)*vo[higher[1]]))
      }
    }
    if (tt[l] >= max(to)) vn[l] <- vo[which(to == max(to))]
    
    
  }
  vn 
}
library(scales)
linp_t <- function(t,a,b,tmin,tmax) {
  # t = time at which to interpolate
  # v1t = value_1_time
  # v2t = value_2_time
  # a = value1
  # b = value2
  
  out <- (a+(b-a)*(my_rescale(c(t,tmin,tmax),c(0,1))[1:length(t)]))
  
  out
}
a <- 2.5
b <- 22.5
tmin <- 0
tmax = 100
t <- bicarb$time[3:46]
#plot(t,linp_t(t,a,b,0,180),type = "o")

#summary(lm(linp_t(t,a,b,0,180)~t))


a <- c(1,3,5,7,3,7,3)
t <- c(1,2,3,4,5,7,10)
tn <- c(0,1,3,3.5,3.6,4.5,8,1,10.2,4,7,8)
#plot(t,a,type = "o", xlim = c(0,10), ylim = c(-2, 10))
#v <- inter(a,t,tn)
#points(tn,v,type = "p", col = "orange")

polygon_between_curves <- function(x1,x2,tstart,tend,tmid,color1 = rgb(1,0,0,1/8), color2 = rgb(0,0,1,1/8)) {
  above<-x1>x2
  # Points always intersect when above=TRUE, then FALSE or reverse
  intersect.points<-which(diff(above)!=0)
  # Find the slopes for each line segment.
  x1.slopes<-(x1[intersect.points+1]-x1[intersect.points])/(tmid[intersect.points+1]-tmid[intersect.points])
  x2.slopes<-(x2[intersect.points+1]-x2[intersect.points])/(tmid[intersect.points+1]-tmid[intersect.points])
  # Find the intersection for each segment.
  x.points<- tmid[intersect.points] + ((x2[intersect.points] - x1[intersect.points]) / (x1.slopes-x2.slopes))
  
  
  y.points<- x1[intersect.points] + ((x1.slopes)*(x.points-tmid[intersect.points]))
  
  
  xy <- cbind(x.points,y.points)
  x1.1 <- cbind(tmid,x1)
  x2.1 <- cbind(tmid,x2)
  
  x1.1_big <- x1.1
  x1.1_big[which(x1.1[,2] < x2.1[,2]),2] <- x2.1[which(x1.1[,2] < x2.1[,2]),2]
  
  x1.1_raw <- rbind(x1.1,xy)
  x1.1_raw <- x1.1_raw[order(x1.1_raw[,1]),]
  
  x1.1_big <- rbind(x1.1_big,xy)
  x1.1_big <- x1.1_big[order(x1.1_big[,1]),]
  polygon( c(tstart, x1.1_big[,1], tend, tend, rev(x1.1_big[,1]), tstart),
           c(x1.1_big[1,2],x1.1_big[,2], x1.1_big[nrow(x1.1_big),1], x1.1_raw[nrow(x1.1_raw),1], rev(x1.1_raw[,2]), x1.1_raw[1,2]),
           border = NA, col = color1)
  
  x2.1_big <- x2.1
  x2.1_big[which(x2.1[,2] < x1.1[,2]),2] <- x1.1[which(x2.1[,2] < x1.1[,2]),2]
  
  x2.1_raw <- rbind(x2.1,xy)
  x2.1_raw <- x2.1_raw[order(x2.1_raw[,1]),]
  
  x2.1_big <- rbind(x2.1_big,xy)
  x2.1_big <- x2.1_big[order(x2.1_big[,1]),]
  polygon( c(tstart, x2.1_big[,1], tend, tend, rev(x2.1_big[,1]), tstart),
           c(x2.1_big[1,2],x2.1_big[,2], x2.1_big[nrow(x2.1_big),1], x2.1_raw[nrow(x2.1_raw),1], rev(x2.1_raw[,2]), x2.1_raw[1,2]),
           border = NA, col = color2)
  
}


myloess <- function(var,t,span) {
  
  dat <- data.frame(x =  t, y = var)
  
  
  predict(loess(y~x, dat, span = span),method = "loess()")
}

x <- c(0.1,0.4,-1.224)
y <- c(22,80)
l <- 3
n <- 1
find_labs <- function(x,n,l) {
  labs <- round(seq(min(x), max(x), (max(x)-min(x))/(l-1)),n)
  labs[which(labs > (min(x)-(max(x) - min(x))*0.02) & labs < (max(x)+(max(x) - min(x))*0.02))]
}
find_labs(x,1,7)

### Custom additive multiple regression


mlr <- function(a1,b1,c1 = NA, d1= NA, y1,aic = "both", xnames = c("a", "b", "c", "d"), yname = "y") {
  if (!(is.na(d1[1]))) {
    dat <- as.data.frame(cbind(a1, b1, c1, d1, y1))
    colnames(dat) = c(xnames, yname)
    model <- lm(formula = y ~ a+ b + c + d, data = dat)
    print(summary(model))
    print(step(model, direction=aic))
    print(summary(step(model, direction=aic)))
    
    
  }
  if (!(is.na(c1[1])) & is.na(d1[1])) {
    dat <- as.data.frame(cbind(a1, b1, c1, y1))
    colnames(dat) = c(xnames[1:3], yname)
    model <- lm(formula = y ~ a+ b + c, data = dat)
    print(summary(model))
    print(step(model, direction=aic))
    print(summary(step(model, direction=aic)))
    
  }
  if (!(is.na(b1[1])) & is.na(c1[1]) & is.na(d1[1])) {
    dat <- as.data.frame(cbind(a1, b1, y1))
    colnames(dat) = c(xnames[1:2], yname)
    model <- lm(formula = y ~ a + b, data = dat)
    print(summary(model))
    print(step(model, direction=aic))
    print(summary(step(model, direction=aic)))
  }
}


### get the formula of a linear model



modelf <- function(x,y) {
  model <- lm(y~x)
  c(model[[1]][1], model[[1]][2])
}

getval <- function(t)  ab[1] + ab[2]*t




range_of_line <- function(vec,overdraw = 0.1) {
  c(min(vec) - (max(vec) - min(vec))*overdraw,
    getval(min(vec) - (max(vec) - min(vec))*overdraw),
    max(vec) + (max(vec) - min(vec))*overdraw,
    getval(max(vec) + (max(vec) - min(vec))*overdraw))
}


customsegment <- function(coord,col = "black", lwd  = 1, lty = 1) {
  segments(coord[1],coord[2],coord[3],coord[4],col = col, lwd = lwd, lty = lty)
}

my_rescale <- function(val,arange) {
  amin <- min(val,na.rm = T)
  amax <- max(val,na.rm = T)
  span <- arange[2]-arange[1]
  (val-amin)*span/(amax-amin)+arange[1]
}
  
