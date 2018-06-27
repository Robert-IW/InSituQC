
library(bfast)
library(tidyverse)
library(forecast)
library(stlplus)
library(kohonen)
library(grid)
library(gridExtra)

col.par = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#CBD588", "#5F7FC7",
            "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
            "#8A7C64", "#599861")


# this function looks for long NA segments that cannot be interpolated satisfactorily
# the longest set of observations
strip.na <- function(x){
  ind.na <- is.na(x)
  con.seq <- rle(ind.na)
  ind <- which(con.seq$value==T & con.seq$length>365)
  if (length(ind)==0){
    print(paste("No long NA sequence in ",stat.name))
    return(seq(1,length(x),1))                                            # if no NA sequences longer than a year
  }
  if (length(con.seq$lengths)!=ind[length(ind)]){       # if last index does not include > 365 NA
    ind <- c(ind, length(con.seq$lengths))
  }
  cnt=1
  ind.st <- matrix(nrow=(length(ind)-1),ncol=2)
  for (z in 1:(length(ind)-1)){
    ind.st[z,1] <- sum(con.seq$length[cnt:(ind[z]-1)])
    ind.st[z,2] <- sum(con.seq$length[(ind[z]+1):(ind[z+1]-1)])
    cnt <- ind[z]+1
  }
  ind.vec <- as.vector(c(ind.st[,1],ind.st[nrow(ind.st),ncol(ind.st)]))
  seg <- which(ind.vec==max(ind.vec))[1]
  if (ind.vec[seg] < 365*3){
    print(paste("Not enough data in ",stat.name))
    next
  }
  if (seg>1){
    ind.extr <- seq(sum(con.seq$lengths[1:ind[seg-1]]),sum(con.seq$lengths[1:(ind[seg]-1)]),1)
  } else {
    ind.extr <- seq(1,sum(con.seq$lengths[1:(ind[seg]-1)]),1)
  }
  print(paste("Data extracted from ", stat.name))
  return(ind.extr)
}


do.insitu <- T
do.sst <- F

data.dir1 <- "/media/robert/KELP-HDD-Portable/SST/GHRSST/L4/"
data.dir3 <- "/ts/ESdata/"
dir.matrix <- matrix(c("AVHRR-OI","MUR","MW_IR_v4.0","OSTIA","CMC"))

if (do.insitu){
  load("~/R/Data/SACTN_daily_v4.2.Rdata")
  load("~/R/Data/site_list_v4.2.RData")
  df.data <- SACTN_daily_v4.2
  df.data.stor <- df.data
} else if (do.sst){
  load(data.file1)
  load(data.file2)
  df.data.stor <- data.in
  df.data <- data.in
  site_list <- site.loc
  rm(data.in, site.loc)
}

df.data <- df.data %>%
  mutate(month=format(date, "%m")) %>%
  group_by(index, month)

bfast.list <- list()
station.nona.list <- list()

for (j in 1:nrow(site_list)){
  
  bfast.mon.list <- list()
  stat.nona <- list()
  
  stat.name <- as.character(site_list$index[site_list$index[j]])
  
  data.in <- df.data.stor[which(df.data.stor$index == site_list$index[[j]]),]
  
  # record data with trend for CP detection
  df.data.trend <- data.frame("date"=data.in$date, 
                              "temp.org"=data.in$temp,
                              "doy"=as.numeric(format(data.in$date, "%j")),
                              "dom"=format(data.in$date, "%d"),
                              "month"=as.numeric(format(data.in$date, "%m")),
                              "year"=format(data.in$date, "%Y"),
                              stringsAsFactors = F)
  
  df.data.trend <- df.data.trend %>% 
    filter(!is.na(doy))
  
  if (length(unique(df.data.trend$year)) < 5){
    print(paste("Time series too short for ", stat.name))
    next
  }
  
  # remove long NA periods or extract useful lengths
  x <- data.in$temp
  y <- strip.na(x)
  data.in <- data.in[y,]
  
  # fill in missing values using double STL
  doy <- format(data.in$date[1], "%j")
  yr <- format(data.in$date[1], "%Y")
  tot.yr <- length(unique(format(data.in$date, "%Y")))
  
  data.ts <- ts(data.in$temp, frequency=365,
                start = c(as.numeric(yr),as.numeric(doy)))
  
  if (length(data.ts) < 365){
    next
  }
  
  na.ind <- which(is.na(data.ts))
  
  if (length(na.ind)/length(data.ts)>0.5){
    print(paste("Too many NA ", stat.name))
    next
  } else if (length(na.ind) > 0){
    data.nona <- na.interp(data.ts,lambda=NULL)
    if (tot.yr < 4){
      t.w <- (tot.yr-1)*365
      s.w <- t.w
    } else {
      t.w <- 5*365
      s.w <- 3*365
    }
    stl.nona <- stlplus(data.nona, s.window = s.w,
                        s.degree = 1, t.window = t.w, n.p = 365)
    trend.ts <- ts(stl.nona$data$trend, frequency=365,
                   start = c(as.numeric(yr),as.numeric(doy)))
    data.nona[na.ind] <- stl.nona$data$seasonal[na.ind]+stl.nona$data$trend[na.ind]
  } else {
    data.nona <- data.ts
    if (tot.yr < 4){
      t.w <- (tot.yr-1)*365
      s.w <- t.w
    } else {
      t.w <- 5*365
      s.w <- 3*365
    }
    stl.nona <- stlplus(data.nona, s.window = s.w,
                        s.degree = 1, t.window = t.w, n.p = 365)
    trend.ts <- ts(stl.nona$data$trend, frequency=365,
                   start = c(as.numeric(yr),as.numeric(doy)))
    
  }
  
  df.data <- data.frame("temp"=data.nona,
                        "temp.org"=data.ts,
                        "date"=data.in$date,
                        "doy"=as.numeric(format(data.in$date, "%j")),
                        "dom"=format(data.in$date, "%d"),
                        "month"=as.numeric(format(data.in$date, "%m")),
                        "year"=format(data.in$date, "%Y"),
                        stringsAsFactors = F)  
  
  for (k in 1:12){
    sub.df <- df.data %>% 
      filter(month==k)
    
    freq <- as.numeric(max(df.data$dom))
    st.yr <- as.numeric(df.data$year[1])
    st.dom <- as.numeric(df.data$dom[1])
    
    # set minimum segment length to 3 months
    bfast.ts <- ts(sub.df$temp, frequency = freq, start=c(st.yr, st.dom))
    if (length(bfast.ts) <= 2*90){
      h.value <- 0.49
    } else {
      h.value <- (3*freq)/length(bfast.ts)
    }
    bfast.res <- bfast(bfast.ts, season = "none", max.iter = 10, h=h.value, level=0.01)
    
    if (bfast.res$output[[1]]$Vt.bp==0){
      stor <- list("series"=bfast.res$Yt, "slopes"=NA, "bps"=NA,
                   "image"=plot(bfast.res, type="trend", main=paste(stat.name, month.name[k])))
    } else {
      stor <- list("series"=bfast.res$Yt, "slopes"=coef(bfast.res$output[[length(bfast.res$output)]]$bp.Vt)[,2],
                   "bps"=bfast.res$output[[2]]$bp.Vt$breakpoints,
                   "image"=plot(bfast.res, type="trend", main=paste(stat.name, month.name[k])))
    }
    
    bfast.mon.list[[k]] <- stor
    stat.nona[[k]] <- sub.df
    rm(stor, bfast.res, bfast.ts)
  }
  
  bfast.list[[j]] <- bfast.mon.list
  station.nona.list[[j]] <- stat.nona
  names(bfast.list)[j] <- stat.name
  names(station.nona.list)[j] <- stat.name
}


# Run SOM on slope coefficients -------------------------------------------

# remove all null elements in list
bfast.list[sapply(bfast.list, is.null)] <- NULL

# for each month of the year do the following
df.som <- list()
for (k in 1:12){
  
  df.month <- list()
  # for each station
  for (m in 1:length(bfast.list)){
    
    stat.name <- names(bfast.list[m])
    
    a <- bfast.list[[m]][[k]]
    
    # if no break points detected get slope of full time series
    if (is.na(a$bps)){
      slopes <- lm(a$series~time(a$series))$coefficients[2]
    } else {
      slopes <- a$slopes
    }
    
    bps <- a$bps
    freq <- frequency(a$series)
    st.dom.yr <- start(a$series)
    en.dom.yr <- end(a$series)
    
    # create a data frame for the data
    yrs.all <- rep(1970:2017, each=freq)
    dom.all <- rep(1:freq, (2017-1970+1))
    date.all.lu <- data.frame("years"=yrs.all, "dom"=dom.all)
    rm(yrs.all, dom.all)
    
    # create date look-up table
    yrs <- c(rep(st.dom.yr[1], freq-st.dom.yr[2]+1),
             rep((st.dom.yr[1]+1):(en.dom.yr[1]-1),each=freq),
             rep(en.dom.yr[1], en.dom.yr[2]))
    doms <- c(seq(st.dom.yr[2],freq,1),
              rep(1:freq, en.dom.yr[1]-st.dom.yr[1]-1),
              seq(1,en.dom.yr[2],1))
    date.lu <- data.frame("years"=yrs, "dom"=doms)
    rm(yrs, doms)
    
    if (is.na(bps)){
      data.vec <- rep(NA, nrow(date.all.lu))
      ind.st <- which(date.all.lu$years==st.dom.yr[1] & date.all.lu$dom==st.dom.yr[2])
      ind.en <- which(date.all.lu$years==en.dom.yr[1] & date.all.lu$dom==en.dom.yr[2])
      data.vec[ind.st:ind.en] <- slopes
    } else {
      # find start month nearest to breakpoints
      bps.new <- numeric()
      for (i in 1:length(bps)){
        dom.old <- date.lu[bps[i],2]
        if (dom.old >= 15){
          dom.new <- which(date.lu$years==(date.lu[bps[i],1]+1) & date.lu$dom==1)
        } else {
          dom.new <- which(date.lu$years==(date.lu[bps[i],1]) & date.lu$dom==1)
        }
        bps.new[i] <- dom.new
      }
      
      # create a vector with coefficients
      data.vec <- rep(NA, nrow(date.all.lu))
      ind.st <- which(date.all.lu$years==st.dom.yr[1] & date.all.lu$dom==st.dom.yr[2])
      for (j in 1:length(bps.new)){
        ind.en <- which(date.all.lu$years==date.lu$years[bps.new[j]] & date.all.lu$dom==date.lu$dom[bps.new[j]])
        data.vec[ind.st:ind.en] <- slopes[j]
        ind.st <- ind.en+1
      }
      ind.en <- which(date.all.lu$years==en.dom.yr[1] & date.all.lu$dom==en.dom.yr[2])
      data.vec[ind.st:ind.en] <- slopes[j+1]
    }
    b <- as.data.frame(data.vec)
    colnames(b) <- stat.name
    df.month[[m]] <- b
  }
  df.som[[k]] <- do.call(cbind, df.month)
}



# SOM on the Bfast Trends -------------------------------------------------

x <- t(do.call(rbind, df.som[1]))

# remove all na columns
x <- x[,-which(colSums(is.na(x))==nrow(x))]

# remove data with less than 5 years (approx 3x 31)
limit <- (ncol(x)-(5*30))/ncol(x)
ind.maxNA <- which(rowSums(is.na(x))/ncol(x) > limit)
x <- x[-ind.maxNA,]

res.som <- som(x, grid = somgrid(4, 4, "hexagonal"),
               rlen=2000, alpha=c(0.05, 0.005), maxNA.fraction=limit)

plot(res.som, main = "Cluster January")
plot(res.som, type="changes")
plot(res.som, type="dist.neighbours")

grp.som <- res.som$unit.classif

plot.som <- cbind(data.frame("index"=rownames(x),"grp"=grp.som),x)
plot.som <- plot.som %>% gather(variable, value, -(grp:index)) %>% arrange(index)
#plot.som <- plot.som %>% gather(variable, value, -(grp:month)) %>% arrange(index)
plot.som$variable <- as.numeric(substr(plot.som$variable,start = 1,stop = 4))

name.grp <- sort(unique(grp.som))
n.grp <- length(name.grp)

plot.list <- list()
cnt = 0
for (i in 1:n.grp){
  cnt = cnt+1
  df.cent <- data.frame(x=plot.som$variable,y=as.vector(res.som$codes[[1]][name.grp[i],]))
  
  plot.list[[cnt]] <- ggplot(plot.som %>% dplyr::filter(grp==name.grp[i]), aes(colour=index)) +
    geom_line(aes(x=variable, y=value, group=index), lwd=1) +
    #geom_line(data=df.cent, aes(x=x, y=y), color="black", lwd=1.5) +
    geom_hline(yintercept = 0, colour="gray30", lwd=.6) +
    ylim(c(-3,3)) + ylab("Temperature (°C)") + xlab("Year") +
    scale_color_manual(name="SOM Cluster Stations",
                       values=col.par[1:24],
                       guide=guide_legend(override.aes=list(alpha=.5, size=2))) +
    #scale_x_continuous(breaks=min(plot.som$variable):max(plot.som$variable)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=9),
          panel.grid.major.x=element_line(colour="gray90",size=0.75))
}
ggsave(file = paste0("~/Desktop/",month.abb[m],".png"),
       arrangeGrob(grobs = plot.list, ncol = 4),  width = 40, height=16)
