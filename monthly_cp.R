library(cowplot)
library(tidyverse)
library(ggpubr)
library(forecast)
library(stlplus)
library(fitdistrplus)
library(gridExtra)
library(RColorBrewer)

#------------------------------------------------------------------------
# Need to be detrended first to remove the potential bias in 
# month on month differences e.g. a positive trend will produce a higher
# positive difference, this will also remove the potential change point
# difference from the probability density

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

# normalize data columns
normalit<-function(m){
  (m - min(m, na.rm=T))/(max(m, na.rm=T)-min(m, na.rm=T))
}

# normalize temperature data for plotting with sd means to range of sd (0 to 2.357977)
normalit2<-function(m){
  (5.4556*(m - min(m, na.rm=T)))/(max(m, na.rm=T)-min(m, na.rm=T))
}

load("~/R/Data/SACTN_daily_v4.2.Rdata")
load("~/R/Data/site_list_v4.2.RData")

do.plot=F
do.mean=T
do.sd=F

# Kmeans Cluster ----------------------------------------------------------
# create month column for SACTN data
df.data <- SACTN_daily_v4.2

df.data <- df.data %>%
  mutate(month=format(date, "%m")) %>%
  group_by(index, month) %>% 
  summarize(mean=mean(temp, na.rm=T),
            med=median(temp, na.rm=T),
            sd=sd(temp, na.rm=T),
            max=max(temp, na.rm=T),
            min=min(temp, na.rm=T),
            range=max(temp, na.rm=T)-min(temp, na.rm=T))


# group by index and month and create summary statistics
df.data.kmeans <- df.data %>% 
  gather(variable, value, -(index:month)) %>%
  unite(temp, month, variable) %>%
  spread(temp, value)

# cluster the data using all monthly stats
ind.na <- apply(df.data.kmeans, 1, function(x) any(is.na(x)))  # find rows containing NA
df.data.kmeans <- df.data.kmeans[!ind.na,]
site_list <- site_list[!ind.na,]

df.data.kmeans.nor <- df.data.kmeans[,-1] %>% mutate_all(funs(normalit))

results.5 <- kmeans(df.data.kmeans.nor, 5, nstart=25)$cluster
p1 <- ggplot(data = site_list, aes(x = lon, y = lat, colour = as.factor(results.5))) +
  borders() + xlab("Longitude") + ylab("Latitude")
  geom_point() +
  labs(colour = "Cluster") +
  coord_equal(xlim = c(15, 35), ylim = c(-37, -27)) +
  ggtitle(paste0("Number of Clusters = ", 5)) +
  theme(panel.grid.major = element_line(colour = "grey"))

ggsave(p1, filename = "~/Desktop/StationCP/ClusterGroups.png", width = 16, height=8)
rm(p1, p2)

if (do.plot){
  # visualize the clusters
  clust_i <- function(h) {
    ggplot(data = site_list, 
           aes(x = lon, y = lat, 
               colour = as.factor(kmeans(df.data.kmeans.nor, h)$cluster))) +
      borders() +
      geom_point() +
      labs(colour = "cluster") +
      coord_equal(xlim = c(15, 35), ylim = c(-37, -27)) +
      ggtitle(paste0("clust = ", h))
  }
  
  clust_6 <- clust_i(6)
  clust_5 <- clust_i(5)
  clust_4 <- clust_i(4)
  clust_3 <- clust_i(3)
  
  # choose cluster of 5 as this includes sufficient time series
  ggarrange(clust_6, clust_5, clust_4, clust_3, common.legend = T)
}

# Get group monthly change stats ------------------------------------------
# for each group subset the data by station, obtain the monthly change stats
cluster.list <- list()
clustertrend.list <- list()
station.list <- list()
stationtrend.list <- list()
statmon.list <- list()
statmontrend.list <- list()

for (i in 1:5){
  
  print(paste("Working on Cluster ",i))
  
  A <- seq(1970,2018,1)
  B <- c("median", "sd")
  
  cols <- as.character(paste(rep(A, each = length(B)), B, sep = "_"))
  rm(A,B)
  
  df.all <- data.frame(matrix(ncol = length(cols)+2))
  colnames(df.all) <- c("index","month",cols)
  df.all.trend <- df.all
  
  station.list.temp <- list()
  stationtrend.list.temp <- list()
  
  statmon.list.temp <- list()
  statmontrend.list.temp <- list()
  ind.grp <- which(results.5==i)
  
  # run through each station and exclude those with > 365 consecutive NA
  for (j in 1:length(ind.grp)){
    
    stat.name <- as.character(site_list$index[ind.grp[j]])
    #print(paste("Working on Station ",stat.name))
    data.in <- SACTN_daily_v4.2[which(SACTN_daily_v4.2$index == site_list$index[ind.grp[j]]),]
    
    # record data with trend for CP detection
    df.data.trend <- data.frame("date"=data.in$date, 
                                "temp.org"=data.in$temp,
                                "doy"=as.numeric(format(data.in$date, "%j")),
                                "dom"=format(data.in$date, "%d"),
                                "month"=as.numeric(format(data.in$date, "%m")),
                                "year"=format(data.in$date, "%Y"),
                                stringsAsFactors = F)
    
    if (length(unique(df.data.trend$year)) < 5){
      print(paste("Time series too short for ", stat.name))
      next
    }
    
    # if the data frame ends before a full month of Jan, delete the incomplete month
    end.day <- df.data.trend$doy[nrow(df.data.trend)]
    if (end.day < 31){
      df.data.trend <- df.data.trend[-(nrow(df.data.trend)-((end.day-1):0)),]
    }
    rm(end.day)
    
    # if the data frame begins after a full month of Dec, delete the incomplete month
    begin.day <- df.data.trend$doy[1]
    if (begin.day > 334){
      # find first occurence of doy=1
      ind.begin <- which(df.data.trend$doy==1)[1]
      df.data.trend <- df.data.trend[-(1:(ind.begin-1)),]
    }
    rm(begin.day, ind.begin)
    
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
    }
    
    plot(x, typ="l", main=paste("Original ", stat.name))
    plot.ts(data.nona, col="red",main=paste("Extracted ", stat.name), ylab="Temperature", xlab="Years")
    lines(data.ts, col="black")

    # remove the trend from the data for monthly mean change
    data.nona <- data.nona-trend.ts
    
    # create data frame with no missing values
    df.data <- data.frame("temp"=data.nona,
                          "date"=data.in$date,
                          "doy"=as.numeric(format(data.in$date, "%j")),
                          "dom"=format(data.in$date, "%d"),
                          "month"=as.numeric(format(data.in$date, "%m")),
                          "year"=format(data.in$date, "%Y"),
                          stringsAsFactors = F)
    
    # if the data frame ends before a full month of Jan, delete the incomplete month
    end.day <- df.data$doy[nrow(df.data)]
    if (end.day < 31){
      df.data <- df.data[-(nrow(df.data)-((end.day-1):0)),]
    }
    rm(end.day)
    
    # if the data frame begins after a full month of Dec, delete the incomplete month
    begin.day <- df.data$doy[1]
    if (begin.day > 334){
      # find first occurence of doy=1
      ind.begin <- which(df.data$doy==1)[1]
      df.data <- df.data[-(1:(ind.begin-1)),]
    }
    rm(begin.day, ind.begin)
    
    # create year/month summaries
    df.data.yrmon <- df.data %>%
      group_by(year, month) %>% 
      summarize(mean=mean(temp, na.rm=T),
                med=median(temp, na.rm=T),
                sd=sd(temp, na.rm=T),
                max=max(temp, na.rm=T),
                min=min(temp, na.rm=T))
    
    st.mon <- as.numeric(sprintf("%02d",df.data.yrmon$month[1]))
    st.year <- as.numeric(df.data.yrmon$year[1])
    end.year <- max(as.numeric(df.data.yrmon$year))
    
    df.data.yrmon$mean <- ts(df.data.yrmon$mean, frequency=12, start=c(st.year, st.mon))
    df.data.yrmon$med <- ts(df.data.yrmon$med, frequency=12, start=c(st.year, st.mon))
    df.data.yrmon$sd <- ts(df.data.yrmon$sd, frequency=12, start=c(st.year, st.mon))
    
    # create monthly summaries for training data
    df.data.mon <- df.data %>%
      dplyr::select(dom, month, year, temp) %>% 
      mutate(year.dom=paste0(year,".",dom)) %>%
      mutate(temp=as.numeric(temp)) %>% 
      dplyr::select(month, year.dom, temp) %>% 
      spread(year.dom, temp) %>% 
      dplyr::select(-month)
    
    df.data.mon$mean <- apply(df.data.mon,1,mean,na.rm=T)
    df.data.mon$sd <- apply(df.data.mon[,!names(df.data.mon) %in% "mean"],1,sd,na.rm=T)
    
    # create a data frame for storage of month-year stats
    A <- as.character(seq(st.year,end.year,1))
    B <- c("median", "sd")
    
    cols <- as.character(paste(rep(A, each = length(B)), B, sep = "_"))
    rm(st.mon, st.year, end.year)
    
    # for training
    df.temp <- data.frame(matrix(ncol = length(cols), nrow = 12))
    colnames(df.temp) <- cols
    
    # create year/month summaries for test data (trend)
    df.data.yrmon.trend <- df.data.trend %>%
      group_by(year, month) %>% 
      summarize(mean=mean(temp.org, na.rm=T),
                med=median(temp.org, na.rm=T),
                sd=sd(temp.org, na.rm=T),
                max=max(temp.org, na.rm=T),
                min=min(temp.org, na.rm=T))
    
    st.mon <- as.numeric(sprintf("%02d",df.data.yrmon.trend$month[1]))
    st.year <- as.numeric(df.data.yrmon.trend$year[1])
    end.year <- max(as.numeric(df.data.yrmon.trend$year))
    
    df.data.yrmon.trend$mean <- ts(df.data.yrmon.trend$mean, frequency=12, start=c(st.year, st.mon))
    df.data.yrmon.trend$med <- ts(df.data.yrmon.trend$med, frequency=12, start=c(st.year, st.mon))
    df.data.yrmon.trend$sd <- ts(df.data.yrmon.trend$sd, frequency=12, start=c(st.year, st.mon))

    # create monthly summaries for testing data (with trend)
    df.data.mon.trend <- df.data.trend %>%
      dplyr::select(dom, month, year, temp.org) %>% 
      mutate(year.dom=paste0(year,".",dom)) %>%
      mutate(temp.org=as.numeric(temp.org)) %>% 
      dplyr::select(month, year.dom, temp.org) %>% 
      spread(year.dom, temp.org) %>% 
      dplyr::select(-month)
    
    df.data.mon.trend$mean <- apply(df.data.mon.trend,1,mean,na.rm=T)
    df.data.mon.trend$sd <- apply(df.data.mon.trend[,!names(df.data.mon.trend) %in% "mean"],1,sd,na.rm=T)
    
    # create a data frame for storage of month-year stats
    A.trend <- as.character(seq(st.year,end.year,1))
    B <- c("median", "sd")
    
    cols.trend <- as.character(paste(rep(A.trend, each = length(B)), B, sep = "_"))
    rm(st.mon, st.year, end.year)
    
    # for testing
    df.test <- data.frame(matrix(ncol = length(cols.trend), nrow = 12))
    colnames(df.test) <- cols.trend
    
    sub.df <- df.data.mon[,1:(ncol(df.data.mon)-2)]
    col.names <- substr(x = colnames(sub.df), start=1, stop=4)
    statmon.list.temp[[j]] <- sub.df
    rm(sub.df)
    
    sub.df <- df.data.mon.trend[,1:(ncol(df.data.mon.trend)-2)]
    col.names.trend <- substr(x = colnames(sub.df), start=1, stop=4)
    statmontrend.list.temp[[j]] <- sub.df
    rm(sub.df)
    
    # for each year get monthly mean and sd
    for (k in 1:12){
      
      mon.name <- month.name[k]
      
      mon <- as.matrix(df.data.mon[k,1:(ncol(df.data.mon)-2)])
      ind <- which(is.na(mon))

      if (length(ind)>0){
        mon <- mon[-ind]
        col.names.mon <- col.names[-ind]
      } else {
        mon <- as.numeric(mon)
        col.names.mon <- col.names
      }
      
      ind.yr <- match(unique(col.names.mon), col.names)
      ind.2 <- which(df.data.yrmon$month==k)
      med <- df.data.yrmon$med[ind.2]
      mean <- df.data.yrmon$mean[ind.2]
      sd <- df.data.yrmon$sd[ind.2]
      stat.stor <- c(rbind(med,sd))
      
      # get index of month start and end year of series from all years
      ind.st <- which(A==col.names.mon[1])
      ind.end <- which(A==tail(col.names.mon,n=1))
      df.temp[k,(ind.st*2-1):(ind.end*2)] <- as.vector(stat.stor)
      
      #-------------------------------------------------------------------- trended
      mon.trend <- as.matrix(df.data.mon.trend[k,1:(ncol(df.data.mon.trend)-2)])
      
      # check if the first month in NA (change this to more than 15 of 31 NA)
      if (all(is.na(mon.trend[1,1:31]))){
        mon.trend <- mon.trend[1,-(1:31)]
        col.names.t <- col.names.trend[-(1:31)]
      } else {
        mon.trend <- as.numeric(mon.trend)
        col.names.t <- col.names.trend
      }
      
      if (all(is.na(mon.trend[length(mon.trend)-(30:0)]))){
        mon.trend <- mon.trend[-(length(mon.trend)-(30:0))]
        col.names.t <- col.names.t[-(length(col.names.t)-(30:0))]
      }

      # check if last 31 days ar NA
      ind.yr <- match(unique(col.names.t), col.names.trend)
      
      ind.2 <- which(df.data.yrmon.trend$month==k)
      med <- df.data.yrmon.trend$med[ind.2]
      mean <- df.data.yrmon.trend$mean[ind.2]
      sd <- df.data.yrmon.trend$sd[ind.2]
      
      stat.stor.trend <- c(rbind(med,sd))
      
      # get index of month start and end year of series from all years
      ind.st <- which(A.trend==col.names.t[1])
      ind.end <- which(A.trend==tail(col.names.t,n=1))
      df.test[k,(ind.st*2-1):(ind.end*2)] <- as.vector(stat.stor.trend)
      
      if (do.plot){
        # for plotting
        med.seg <- data.frame(x1=ind.yr-1, y1=med,
                              x2=c(ind.yr[2:length(ind.yr)],length(mon)), y2=med)
        mean.seg <- data.frame(x1=ind.yr-1, y1=mean,
                               x2=c(ind.yr[2:length(ind.yr)],length(mon)), y2=mean)
        sd.seg <- data.frame(x1=ind.yr-1, y1=sd,
                             x2=c(ind.yr[2:length(ind.yr)],length(mon)), y2=sd)
        
        #par(mar = c(5,5,2,5))                                          # to plot a second y axis
        #par(mfrow=c(2,1), mar=c(4,4,4,1))
        layout(matrix(c(1,1,2), nrow = 3, ncol = 1, byrow = TRUE))
        plot(mon, typ="l", main=paste(mon.name, stat.name), xaxt="n", ylab="Temperature", xlab="")
        #axis(1, at=ind.yr, labels=col.names[ind.yr])
        abline(v=ind.yr, col="gray30")
        segments(med.seg$x1, med.seg$y1, med.seg$x2, med.seg$y2, col="red", lwd=2)
        segments(mean.seg$x1, mean.seg$y1, mean.seg$x2, mean.seg$y2, col="green", lwd=2)
        legend("bottomleft", c("Median","Mean"),
               lty=c(1,1), lwd=c(2,2),col=c("red","green"),
               text.width=c(0,1,1.5),
               inset=c(0,.4,0), xpd=T, horiz=T, bty='n')
        
        plot(mon, ylim=c(0,2), type="n", ylab="Standard Deviation", xaxt="n", xlab="Years")
        segments(sd.seg$x1, sd.seg$y1, sd.seg$x2, sd.seg$y2, col="blue", lwd=2)
        axis(1, at=ind.yr, labels=col.names[ind.yr])
        abline(v=ind.yr, col="gray30")
      }
    }
    # condense months for each station
    df.temp <- cbind("index"=stat.name, "month"=1:12, df.temp)
    df.all <- bind_rows(df.all, df.temp)
    station.list.temp[[j]] <- df.temp
    
    df.test <- cbind("index"=stat.name, "month"=1:12, df.test)
    df.all.trend <- bind_rows(df.all.trend, df.test)
    stationtrend.list.temp[[j]] <- df.test
    rm(df.temp)
  }
  # condense stations into single data.frame
  cluster.list[[i]] <- df.all
  clustertrend.list[[i]] <- df.all.trend
  
  station.list[[i]] <- station.list.temp
  stationtrend.list[[i]] <- stationtrend.list.temp
  
  statmon.list[[i]] <- statmon.list.temp
  statmontrend.list[[i]] <- statmontrend.list.temp
  rm(df.all, station.list.temp, statmon.list.temp)
}


# Summarize the Cluster Data ----------------------------------------------

# Summaries of the data by month for indivual years and all year are obtained.
# Only the latter are used for the change point detection
# Changes are detected for the non-detrended data
for (c in 1:length(cluster.list)){
  cluster.ind <- c
  
  # data for the population parameters
  cluster.1 <- cluster.list[[c]][-1,]
  
  if (do.mean){
    temp1 <- cluster.1 %>%
      dplyr::select(month, matches("median"))
  } else if (do.sd) {
    temp1 <- cluster.1 %>%
      dplyr::select(month, matches("sd"))
  }
  
  temp2 <- temp1
  temp2[,2:ncol(temp2)] <- NA
  temp3 <- temp2
  
  # normalize and get first difference
  for (i in 1:nrow(temp1)){
    temp2[i,2:ncol(temp2)] <- normalit(temp1[i,2:ncol(temp1)])
    temp3[i,3:ncol(temp3)] <- diff(as.numeric(temp2[i,2:ncol(temp2)]), lag=1)
  }
  
  temp4 <- temp3 %>% 
    group_by(month) %>% 
    summarise_all(funs(mean(., na.rm=T), sd(.,na.rm=T), n=sum(!is.na(.))))
  
  # plot out the monthly mean and confidence/error bars for each month
  # create a more user friendly df for ggplot
  temp5 <- gather(temp4, month)
  t1 <- nrow(temp5)/3
  t2 <- t1+1
  t3 <- t1+t2
  temp6 <- cbind(rep(1970:2018, each=12), rep(1:12, length.out=t1),
                 temp5[1:t1,2],temp5[t2:(t3-1),2],temp5[t3:nrow(temp5),2])
  colnames(temp6) <- c("year","month","mean", "sd", "count")
  temp6$month <- month.abb[temp6$month]
  temp6$month <- factor(temp6$month, levels=c("Jan","Feb","Mar","Apr",
                                              "May","Jun","Jul","Aug",
                                              "Sep","Oct","Nov","Dec"))
  
  # the monthly mean and sd for the cluster
  temp7 <- data.frame(month=1:12, mean=rep(NA, 12), sd=rep(NA, 12))
  for (i in 1:12){
    temp7$mean[i] <- mean(as.matrix(subset(temp3, temp3$month==i)[,-1]), na.rm=T)
    temp7$sd[i] <- sd(as.matrix(subset(temp3, temp3$month==i)[,-1]), na.rm=T)
  }
  
  ggplot(temp6, aes(mean, group=month, color=as.factor(month))) +
    geom_density(alpha=.2) +
    labs(x = 'Data', y = 'Density') +
    ggtitle(paste("Density Plots of Monthly Year-on-Year Change - Cluster ",c))
  
  ggplot(data.frame(x = c(-1, 1)), aes(x)) + 
    mapply(function(mean, sd, col) {
      stat_function(fun = dnorm, args = list(mean = mean, sd = sd), col = col)
    }, 
    mean = temp7$mean, 
    sd = temp7$sd,
    col = brewer.pal(n=12, "Paired")
    )
  
  temp8 <- data.frame(mean=mean(as.matrix(temp3[,2:ncol(temp3)]), na.rm=T))
  temp8$sd <- sd(as.matrix(temp3[,2:ncol(temp3)]), na.rm=T)
  
  if (do.mean){
    plot.tit1 <- paste("Cluster ",c, " Group Mean")
    filename1 <- paste0("~/Desktop/StationCP/cluster_",c,"_Series_Mean.png")
  } else if (do.sd){
    plot.tit1 <- paste("Cluster ",c, " Group SD")
    filename1 <- paste0("~/Desktop/StationCP/cluster_",c,"_Series_SD.png")
  }
  
  plot.tit2 <- paste("Cluster ",c, " Number of Series per Year")
  
  p1 <- ggplot(data=temp6, aes(y=mean, x=year, group=month)) +
    xlab("") +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size=.2) +
    geom_hline(aes(yintercept = 0), color="red", size=.5) +
    #geom_vline(data=xint, aes(xintercept = year, color=factor(grp)), size=.3, show.legend = F) +
    facet_grid(month ~ ., margins = F, scales = "fixed") +
    ggtitle(plot.tit1) +
    #scale_fill_discrete(guide=FALSE) +
    #scale_x_date(date_labels = "%Y", date_breaks = "2 year") +
    scale_color_manual(values=group.colors) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          panel.grid.major.x=element_line(colour="gray30",size=0.75),
          panel.background = element_rect(fill = 'gray90', colour = 'red'))
  
  p2 <- ggplot(data=temp6, aes(y=count, x=year, group=year)) +
    geom_bar(stat = "summary", fun.y = "mean", width=.3, fill="gray30") +
    geom_hline(yintercept = c(10,20), col="gray50") +
    ggtitle(plot.tit2) +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          panel.background = element_rect(fill = 'gray90', colour = 'red'))
  
  #plot_grid(p1, p2, ncol=1, nrow=2, align="v", rel_heights = c(1,.2), axis="lr")

  #grid.arrange(grobs=plot.list, ncol=1)
  ggsave(file = filename1,
         plot_grid(p1, p2, ncol=1, nrow=2, align="v", rel_heights = c(1,.2), axis="lr"),
         width = 12, height=16)
  rm(p1,p2)
  
  test <- gather(temp3, year, value, -month)
  
  ggplot(data=test, aes(value, group=month)) +
    #geom_histogram(binwidth = 0.005) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(test$value,na.rm=T), sd = sd(test$value, na.rm=T)),
                  col="red") +
    facet_grid(month ~ ., margins = FALSE, scales = "free") +
    ggtitle(paste("Cluster ",c," Temperature Difference Frequency by Month"))
  
  # Note this is the mean of the Median or Std Dev Difference
  ggplot(data=temp6, aes(mean)) +
    geom_histogram(binwidth = 0.02) +
    ggtitle(paste("Cluster ",c," Temperature Difference Frequency"))
  
  
  # Detect Changes for Individual Series in Cluster
  # For each series in cluster obtain a monthly series,
  # normalize the data,
  # obtain the difference and 
  # determine the probability of the difference given the cluster monthly probability density
  
  # NOTE: dnorm determines the height of a given value assuming normal distribution with given mean and sd
  #       pnorm determines the probability that a random number will be less than a given number (tail=F givens p larger)
  test <- as.vector(as.matrix(temp3[,2:ncol(temp3)]))
  test <- test[!is.na(test)]
  
  descdist(test, discrete = F)  # the results indicate normal or logistic
  # fit.norm <- fitdist(test, "norm")
  # fit.logist <- fitdist(test, "logist")
  # fit.norm$bic
  # fit.logist$bic
  
  # get the 'population' mean and std if year to year change
  sample.mu <- temp8$mean
  sample.sd <- temp8$sd
  
  # get the z-score for each month difference for the cluster
  # data for the evaluation parameters
  cluster.1b <- clustertrend.list[[c]][-1,]
  
  if (do.mean){
    temp1b <- cluster.1b %>%
      dplyr::select(month, matches("median"))
  } else if (do.sd){
    temp1b <- cluster.1b %>%
      dplyr::select(month, matches("sd"))
  }
  
  temp2b <- temp1b
  temp2b[,2:ncol(temp2b)] <- NA
  temp3b <- temp2b
  
  # normalize and get first difference
  for (i in 1:nrow(temp1b)){
    temp2b[i,2:ncol(temp2b)] <- normalit(temp1b[i,2:ncol(temp1b)])
    temp3b[i,3:ncol(temp3b)] <- diff(as.numeric(temp2b[i,2:ncol(temp2b)]), lag=1)
  }
  
  temp9 <- (temp3b[,2:ncol(temp3b)]-sample.mu)/sample.sd
  
  # determine potential change points for each series
  stat.names <- clustertrend.list[[c]][,1]
  stat.names <- stat.names[-1]
  temp9 <- cbind(stat.names, temp9)
  stat.mon <- temp3b[,1]
  stat.unique <- as.vector(unique(site_list$index[which(results.5==c)]))
  
  # remove those names where the statmon.list is empty (rejected stations)
  ind.tmp <- sapply(statmontrend.list[[c]], function(x) is_empty(x))
  ind.s <- which(ind.tmp==F)
  
  for (s in ind.s){
    
    print(s)
    stat <- stat.unique[s]
    
    data.ts <- SACTN_daily_v4.2 %>% 
      dplyr::select(date,index,temp) %>% 
      filter(index==stat)

    doy <- format(data.ts$date[1], "%j")
    yr <- format(data.ts$date[1], "%Y")
    data.plot <- ts(data.ts$temp, frequency=365,
                  start = c(as.numeric(yr),as.numeric(doy)))
    
    stat.daily <- statmontrend.list[[c]][[s]]       # this should be the original data
    col.dates <- as.numeric(substr(colnames(stat.daily),start = 1,stop = 4))
    colnames(stat.daily) <- col.dates
    stat.daily <- as.data.frame(t(stat.daily), stringsAsFactors = F)
    colnames(stat.daily) <- month.abb[1:12]
    stat.daily <- suppressWarnings(cbind("year"=as.numeric(rownames(stat.daily)), stat.daily))
    
    stat.zscore <- temp9[temp9$stat.names==stat,][,-1]

    #colnames(stat.zscore) <- col.dates
    #stat.zscore <- stat.zscore %>% 
    #  select_if(~!all(is.na(.)))
    
    # remove all empty columns before and after first and last data
    st <- which(colSums(is.na(stat.zscore))<nrow(stat.zscore))[1]
    en <- tail(which(colSums(is.na(stat.zscore))<nrow(stat.zscore)),n=1)
    #stat.zscore <- stat.zscore[,st:en]
    stat.zscore <- stat.zscore[,st:en]
    col.dates <- unique(as.numeric(substr(colnames(stat.zscore),start = 1,stop = 4)))
    stat.zscore <- cbind(NA, stat.zscore)
    colnames(stat.zscore) <- c(col.dates[1]-1,col.dates)
    rm(st, en)
    #colnames(stat.zscore) <- unique(col.dates)
    stat.zscore <- as.data.frame(t(stat.zscore), stringsAsFactors = F)
    colnames(stat.zscore) <- month.abb[1:12]
    stat.zscore <- cbind("year"=as.numeric(rownames(stat.zscore)), stat.zscore)
    
    stat.mean <- stationtrend.list[[c]][[s]]
    if (do.mean){
      stat.mean <- stat.mean %>% 
        dplyr::select(matches("median"))
    } else if (do.sd){
      stat.mean <- stat.mean %>% 
        dplyr::select(matches("sd"))
    }

    colnames(stat.mean) <- unique(as.numeric(substr(colnames(stat.mean),start = 1,stop = 4)))
    stat.mean <- as.data.frame(t(stat.mean), stringsAsFactors = F)
    colnames(stat.mean) <- month.abb[1:12]
    stat.mean <- cbind("year"=as.numeric(rownames(stat.mean)),stat.mean)
    
    new <- as.data.frame(as.numeric(stat.daily[,1]))
    colnames(new) <- "year"
    
    stat.mean <- merge(new, stat.mean)
    stat.zscore <- merge(new, stat.zscore)
    
    ind.tick <- which(!duplicated(stat.daily$year))
    ind.label <- as.character(stat.daily$year[ind.tick])
    
    ind.x <- ind.tick
    ind.xend <- c(ind.tick[-1],nrow(stat.mean))
    
    yr.fun <- function (x) {
      out <- tryCatch(list(unique(stat.zscore$year[x])),
                      error = function(e) NULL)
      return(out)
    }
    
    ind.fun <- function(x){
      out <- tryCatch(sapply(x, function(y) which(stat.zscore==y)[1]),
                      error = function(e) NULL)
      return(out)
    }
    
    plot.list=list()
    plot.add=list()
    
    tick.date <- which(format(data.ts$date,"%d")=="01" & format(data.ts$date, "%m")=="01")
    plot.add[[1]] <- 
      ggplot(data.ts, aes(date, temp)) +
      geom_line() +
      ggtitle(paste(stat, " -  Original Daily Temperture")) +
      xlab("Date") +
      ylab("Temperature (°C)") +
      scale_x_date(date_breaks="1 year",date_labels=format("%Y")) +
      geom_vline(xintercept=as.numeric(data.ts$date[tick.date]), color="gray30", size=.8) +
      theme(plot.title = element_text(size = 20, face = "bold", hjust=.5),
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    for (m in 1:12){
      
      stat.data <- as.data.frame(cbind(stat.daily[,1], stat.daily[,m+1]), stringsAsFactors = F)
      colnames(stat.data) <- c("year","month")
      ind.y <- stat.mean[,m+1][ind.tick]
      ind.yend <- ind.y
      df.seg <- data.frame(year=as.factor(ind.label),x=ind.x, xend=ind.xend, y=ind.y, yend=ind.yend)
      
      ind.0.1 <- which(abs(stat.zscore[,m+1]) > 1.6449 & abs(stat.zscore[,m+1]) <= 1.96)
      if (length(ind.0.1)==0){
        ind.yr.0.1 <- NA
      } else {
        yr.0.1 <- yr.fun(ind.0.1)
        ind.yr.0.1 <- ind.fun(yr.0.1)
      }
      
      ind.0.05 <- which(abs(stat.zscore[,m+1]) > 1.96 & abs(stat.zscore[,m+1]) <= 2.5758)
      if (length(ind.0.05)==0){
        ind.yr.0.05 <- NA
      } else {
        yr.0.05 <- yr.fun(ind.0.05)
        ind.yr.0.05 <- ind.fun(yr.0.05)
      }
      
      ind.0.01 <- which(abs(stat.zscore[,m+1]) > 2.5758)
      if (length(ind.0.01)==0){
        ind.yr.0.01 <- NA
      } else {
        yr.0.01 <- yr.fun(ind.0.01)
        ind.yr.0.01 <- ind.fun(yr.0.01)
      }
      
      if (do.mean){
        ylab.plot <- "Temperature (°C)"
        ylim.plot <- c(10,28)
      } else if (do.sd){
        ylab.plot <- "Std Dev"
        ylim.plot <- c(0,2.5)
      }
      
      plot.list[[m]] <- 
        ggplot(data=stat.data, aes(x=(1:nrow(stat.data)), y=month, group=year)) +
        {if (do.mean) geom_line(size=.8)} +
        {if (do.sd) geom_line(y=normalit2(stat.data$month),size=.8)} +
        ylab(ylab.plot) +
        ylim(ylim.plot) +
        ggtitle(month.name[m]) +
        geom_vline(xintercept = ind.tick, color="white", size=.5) +
        geom_segment(data=df.seg,  aes(x=x, xend=xend, y=y, yend=yend), col="blue", size=1.5, alpha=.7) +
        scale_x_continuous(name="Years", breaks=ind.tick, labels=ind.label) +
        {if (!is.na(ind.yr.0.01)) geom_vline(xintercept = ind.yr.0.01, color="red", size=1)} +
        {if (!is.na(ind.yr.0.05)) geom_vline(xintercept = ind.yr.0.05, color="darkorange", size=1)} +
        {if (!is.na(ind.yr.0.1)) geom_vline(xintercept = ind.yr.0.1, color="yellow", size=1)} +
        geom_vline(data=data.frame(x=-Inf),aes(xintercept=x, color="red"), alpha=1, show.legend=TRUE)+
        geom_vline(data=data.frame(x=-Inf),aes(xintercept=x, color="green"), alpha=0, show.legend=TRUE)+
        geom_vline(data=data.frame(x=-Inf),aes(xintercept=x, color="blue"), alpha=0, show.legend=TRUE) +
        scale_colour_manual(name = 'Probability', 
                            values =c(blue="red",green="darkorange",red="yellow"),
                            labels = c('p < 0.01','p < 0.05', 'p < 0.1')) +
        theme(legend.position="top", legend.box = "horizontal",
              axis.text.x = element_text(angle = 45, hjust = 1, size=10),
              panel.background = element_rect(fill = 'gray90', colour = 'red'))
    }
    
    plot.grob <- c(plot.add, plot.list)
    
    filename1 <- gsub(" ","",stat)
    filename1 <- gsub("/","_",filename1)
    
    if (do.mean){
      filename1 <- paste0("~/Desktop/StationCP/cluster_",c,"_",filename1,"_meanCP.png")
    } else if (do.sd){
      filename1 <- paste0("~/Desktop/StationCP/cluster_",c,"_",filename1,"_sdCP.png")
    }
    
    #grid.arrange(grobs=plot.list, ncol=1)
    ggsave(file = filename1,
           arrangeGrob(grobs = plot.grob, ncol = 2, layout_matrix = cbind(1:7,c(1,8:13))),
           width = 12, height=16)
    rm(plot.list)
  }
}
