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

load("~/R/Data/SACTN_daily_v4.2.Rdata")
load("~/R/Data/site_list_v4.2.RData")

do.plot=F

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

# normalize data columns
normalit<-function(m){
  (m - min(m, na.rm=T))/(max(m, na.rm=T)-min(m, na.rm=T))
}
df.data.kmeans.nor <- df.data.kmeans[,-1] %>% mutate_all(funs(normalit))

results.5 <- kmeans(df.data.kmeans.nor, 5, nstart=25)$cluster
ggplot(data = site_list, aes(x = lon, y = lat, colour = as.factor(results.5))) +
  borders() +
  geom_point() +
  labs(colour = "cluster") +
  coord_equal(xlim = c(15, 35), ylim = c(-37, -27)) +
  ggtitle(paste0("clust = ", 5))

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
station.list <- list()
statmon.list <- list()

for (i in 1:5){
  
  A <- seq(1970,2018,1)
  B <- c("median", "sd")
  cols <- as.character(paste(rep(A, each = length(B)), B, sep = "_"))
  df.all <- data.frame(matrix(ncol = length(cols)+2))
  colnames(df.all) <- c("index","month",cols)
  
  station.list.temp <- list()
  statmon.list.temp <- list()
  ind.grp <- which(results.5==i)
  
  # run through each station and exclude those with > 365 consecutive NA
  
  for (j in 1:length(ind.grp)){
    
    stat.name <- as.character(site_list$index[j])
    data.in <- SACTN_daily_v4.2[which(SACTN_daily_v4.2$index == site_list$index[j]),]
    
    # remove long NA periods or extract useful lengths
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
      j=1
      ind.st <- matrix(nrow=(length(ind)-1),ncol=2)
      for (i in 1:(length(ind)-1)){
        ind.st[i,1] <- sum(con.seq$length[j:(ind[i]-1)])
        ind.st[i,2] <- sum(con.seq$length[(ind[i]+1):(ind[i+1]-1)])
        j <- ind[i]+1
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
    
    x <- data.in$temp
    plot(x, typ="l", main=paste("Original ", stat.name))
    y <- strip.na(x)
    plot(x[y],typ="l", main=paste("Extracted ", stat.name))
    
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
    
    if (length(na.ind) > 0){
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
      
      if (do.plot){
        plot.ts(data.nona, col="red", main=stat.name, ylab="Temperature", xlab="Years")
        lines(data.ts, col="black")
      }
    } else {
      data.nona <- data.ts
    }
    
    # remove the trend from the data
    data.nona <- data.nona-trend.ts
    
    # create data frame with no missing values
    df.data <- data.frame("temp"=data.nona,
                          "temp.org"=as.numeric(data.ts),
                          "date"=data.in$date,
                          "doy"=as.numeric(format(data.in$date, "%j")),
                          "dom"=format(data.in$date, "%d"),
                          "month"=as.numeric(format(data.in$date, "%m")),
                          "year"=format(data.in$date, "%Y"),
                          stringsAsFactors = F)
    
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
    
    # create monthly summaries
    df.data.mon <- df.data %>%
      dplyr::select(dom, month, year, temp) %>% 
      mutate(year.dom=paste0(year,".",dom)) %>%
      mutate(temp=as.numeric(temp)) %>% 
      dplyr::select(month, year.dom, temp) %>% 
      spread(year.dom, temp) %>% 
      dplyr::select(-month)
    
    df.data.mon$mean <- apply(df.data.mon,1,mean,na.rm=T)
    df.data.mon$sd <- apply(df.data.mon[,!names(df.data.mon) %in% "mean"],1,sd,na.rm=T)
    
    sub.df <- df.data.mon[,1:(ncol(df.data.mon)-2)]
    col.names <- substr(x = colnames(sub.df), start=1, stop=4)
    statmon.list.temp[[j]] <- sub.df
    rm(sub.df)
    
    # create a data frame for storage of month-year stats
    A <- as.character(seq(st.year,end.year,1))
    B <- c("median", "sd")
    
    cols <- as.character(paste(rep(A, each = length(B)), B, sep = "_"))
    
    df.temp <- data.frame(matrix(ncol = length(cols), nrow = 12))
    colnames(df.temp) <- cols
    
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
    rm(df.temp)
  }
  # condense stations into single data.frame
  cluster.list[[i]] <- df.all
  station.list[[i]] <- station.list.temp
  statmon.list[[i]] <- statmon.list.temp
  rm(df.all, station.list.temp, statmon.list.temp)
}


# Summarize the Cluster Data ----------------------------------------------
# Summaries of the data by month for indivual years and all year are obtained. Only the latter are
# used for the change point detection
for (c in 1:length(cluster.list)){
  cluster.ind <- c
  cluster.1 <- cluster.list[[c]][-1,]
  
  temp1 <- cluster.1 %>%
    dplyr::select(month, matches("median"))
  
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
  
  ggplot(data=temp6, aes(y=mean, x=year, group=month)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size=.2) +
    geom_hline(aes(yintercept = 0), color="red", size=.5) +
    #geom_vline(data=xint, aes(xintercept = year, color=factor(grp)), size=.3, show.legend = F) +
    facet_grid(month ~ ., margins = FALSE, scales = "free") +
    ggtitle(paste("Cluster ",c)) +
    #scale_fill_discrete(guide=FALSE) +
    #scale_x_date(date_labels = "%Y", date_breaks = "2 year") +
    scale_color_manual(values=group.colors) +
    theme(axis.text.x = element_text(angle = 45))
  
  test <- gather(temp3, year, value, -month)
  ggplot(data=test, aes(value, group=month)) +
    #geom_histogram(binwidth = 0.005) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(test$value,na.rm=T), sd = sd(test$value, na.rm=T)),
                  col="red") +
    facet_grid(month ~ ., margins = FALSE, scales = "free") +
    ggtitle(paste("Cluster ",c," Temperature Difference Frequency by Month"))
  
  ggplot(data=temp6, aes(mean)) +
    geom_histogram(binwidth = 0.02) +
    ggtitle(paste("Cluster ",c," Temperature Difference Frequency"))
  
  
  # Detect Changes for Individual Series in Cluster -------------------------
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
  
  # get the z-score for each month difference for the cluster
  sample.mu <- temp8$mean
  sample.sd <- temp8$sd
  temp9 <- (temp3[,2:ncol(temp3)]-sample.mu)/sample.sd
  
  # determine potential change points for each series
  
  stat.names <- cluster.list[[1]][,1]
  stat.names <- stat.names[-1]
  stat.mon <- temp3[,1]
  stat.unique <- unique(stat.names)
  
  for (s in 1:length(stat.unique)){
    stat <- stat.unique[s]
    
    stat.daily <- statmon.list[[1]][[s]]
    
    # create data frame with all objects
    start.yr <- as.numeric(substr(colnames(stat.daily[1]),start = 1,stop = 4))
    end.yr <- as.numeric(substr(colnames(stat.daily[ncol(stat.daily)]),start = 1,stop = 4))
    col.dates <-  as.character(rep(start.yr:end.yr, each=31))
    
    colnames(stat.daily) <- col.dates
    stat.daily <- as.data.frame(t(stat.daily), stringsAsFactors = F)
    colnames(stat.daily) <- month.abb[1:12]
    stat.daily <- suppressWarnings(cbind("year"=as.numeric(rownames(stat.daily)), stat.daily))
    
    stat.ind <- which(stat.names==stat)
    stat.zscore <- temp9[stat.ind,]
    stat.zscore <- stat.zscore %>% 
      select_if(~!all(is.na(.)))
    stat.zscore <- cbind(NA, stat.zscore)
    colnames(stat.zscore) <- unique(col.dates)
    stat.zscore <- as.data.frame(t(stat.zscore), stringsAsFactors = F)
    colnames(stat.zscore) <- month.abb[1:12]
    stat.zscore <- cbind("year"=as.numeric(rownames(stat.zscore)), stat.zscore)
    
    stat.mean <- station.list[[1]][[s]]
    stat.mean <- stat.mean %>% 
      dplyr::select(matches("median"))
    colnames(stat.mean) <- unique(col.dates)
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
      
      plot.list[[m]] <- 
        ggplot(data=stat.data, aes(x=(1:nrow(stat.data)), y=month, group=year)) +
        geom_line(size=.8) +
        ylab("Normalized Temperature") +
        ylim(c(-3,7)) +
        ggtitle(paste(stat, "    ", month.name[m], " Detrended Daily Temperture")) +
        geom_vline(xintercept = ind.tick, color="gray30", size=.5) +
        geom_segment(data=df.seg,  aes(x=x, xend=xend, y=y, yend=yend), col="blue", size=3, alpha=.5) +
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
        theme(legend.position="top", legend.box = "horizontal")
    }
    
    filename1 <- gsub(" ","",stat)
    filename1 <- gsub("/","_",filename1)
    filename1 <- paste0("~/Desktop/StationCP/cluster_",c,"_",filename1,".png")
    
    #grid.arrange(grobs=plot.list, ncol=1)
    ggsave(file = filename1,
           arrangeGrob(grobs = plot.list, ncol = 2, layout_matrix = cbind(1:6,7:12)),
           width = 12, height=16)
  }
}
