### Very general plotting script on period level

#### These need to be specified:
# first_period
# last_period
# age_at
# ylimit
# xlable
# ylable
# lablesize
# title
# xlabfactor
# boxfactor
#####
ytop <- max(ylimit)
ybase <- min(ylimit)
labelsize <- 1
xlabsize <- labelsize
  ylabel <- ""

  y1 <- ylimit[1] # lower boundary of plot area
  y2 <- ylimit[1]-boxfactor*(ylimit[2]-ylimit[1]) # position of x axis


all_period <- c("Cb","O", "S", "D", "C", "P", 
            "T", "J", "K", "Pg", "N")

all_period_start <- c(541,485.4,443.8,419.2,358.9,298.9,251.902,201.3,145,66,23.03)
all_period_end <- c(485.4,443.8,419.2,358.9,298.9,251.902,201.3,145,66,23.03,0)
all_periods <- c("Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian","Triassic","Jurassic","Cretaceous", "Paleogene", "Neogene")
periodindex <- match(last_period, all_periods):match(first_period, all_periods)
firstindex <- periodindex[1]
lastindex <- periodindex[length(periodindex)]
ylabpush = -0.14009*(all_period_start[firstindex]-all_period_end[lastindex])
  
periodbase <- -all_period_start[periodindex]
periodtop <- -all_period_end[periodindex]
periodmid <- (periodbase + periodtop) / 2

period <- all_period[periodindex]


plot(0,0, xlim = c(min(periodbase), 0), ylim = ylimit, xaxs="i",
     xlab="", xaxt="n", ylab = NA, yaxs="i", type="n", axes = FALSE, main = title)

for (i in 1:length(periodbase)) {
  polygon(c(periodbase[i],periodtop[i],periodtop[i],periodbase[i]),xpd=TRUE,
          c(y1,y1,y2,y2), col= "white", border= "black")
}




abline(h=y1)

text(periodmid, mean(c(y1,y2)), 
     labels = period, srt = 0, adj = c(0.5,0.5), xpd = TRUE, cex=lablesize)

text(ylabpush+periodbase[1], mean(ylimit),
     labels = ylable, srt = 90, adj = c(0.5,0.5), xpd = TRUE, cex=lablesize)


axis_at <- seq(0, -all_period_start[firstindex], by = -age_at)

lablesize <- 0.9

xlaxtadjust <- 2.71828^(-lablesize)
axis(side=1, at = axis_at, labels = -axis_at, tick = TRUE, line = NA,
     pos = y2, outer = FALSE, font = NA, lty = "solid",
     lwd = 1, lwd.ticks = 1, col = NULL, col.ticks = NULL,
     hadj = NA, padj = -xlaxtadjust, cex.axis=lablesize)

text(mean(c(min(periodbase), 0)),y1 -xlabfactor*(y1-y2),labels = xlable, xpd = TRUE)

abline(h = ylimit, v = 0)

