
library(tidyverse)
library(stlplus)
library(gridExtra)
        


getOutlier <- function(df, mon, mean, sd){
  lo <- mean-(4*sd)
  hi <- mean+(4*sd)
  idx <- which(df$month == mon & df$grad_res > lo & df$grad_res < hi)
  return(idx)
}

getSpike <- function(df, mon, mean, sd){
  lo <- mean-(4*sd)
  hi <- mean+(4*sd)
  idx <- which(df$month == mon & df$grad < lo & df$grad[-1] > hi |
                 df$month == mon & df$grad > hi & df$grad[-1] < lo)
  return(idx)
}

getResolut <- function(df, temp){
  
}
  

load("~/R/Data/SACTN_daily_v4.2.Rdata")
load("~/R/Data/site_list_v4.2.RData")

for (y in 1:nrow(site_list)){
  # sub set by individual series
  site_name <- site_list$index[y]
  data_day <- subset(SACTN_daily_v4.2, index==site_list$index[y])
  data_day$year <- format(data_day$date,"%Y")
  data_day$month <- format(data_day$date,"%m")
  data_day <- data_day %>%
    mutate(grad = c(NA,diff(temp,1)))
  
  # get doy from first date for ts
  doy <- format(data_day$date[1], "%j")
  yr <- format(data_day$date[1], "%Y")
  ts_day <- ts(data_day$temp, frequency=365, start = c(as.numeric(yr), as.numeric(doy)))
  
  tot.day <- nrow(data_day)
  if (tot.day < 365 ){
    next
  }
  # # apply STL to ts
  # tot.yr <- length(unique(data_day$year))
  # if (tot.yr <= 1){
  #   next
  # } else if (tot.yr ==2 & nrow(data_day) < 365) {
  #   next
  # } else if (tot.yr < 4){
  #   t.w <- (tot.yr-1)*365
  #   s.w <- t.w
  # } else {
  #   t.w <- 5*365
  #   s.w <- 3*365
  # }
  # decmp <- stlplus(ts_day, s.window = s.w,
  #                  s.degree = 1, t.window = t.w, n.p = 365)[[1]]
  # data_day$res <- decmp$remainder
  # data_day$raw <- decmp$raw
  # data_day$seas <- decmp$seasonal
  # data_day$trend <- decmp$trend
  # 
  # grad_resid <- diff(data_day$res, lag=1)
  # data_day$grad_res <- NA
  # data_day$grad_res[2:nrow(data_day)] <- grad_resid
  
  stat_monGrad <- data_day %>%
    group_by(month) %>%
    summarise(mean = mean(grad, na.rm=T), sd = sd(grad, na.rm=T))

  sdres_ls <- split(stat_monGrad, seq(nrow(stat_monGrad)))
  data_day$spike <- FALSE
  spike_ls <- lapply(sdres_ls, function(x) getSpike(data_day, x[[1]], x[[2]], x[[3]]))
  
  for (i in 1:12){
    data_day$spike[spike_ls[[i]]] <- TRUE
  }
  
  # get index for NA and replace TRUE with NA
  idx <- is.na(data_day$grad)
  data_day$spike[idx] <- NA
  
  sample_spike <- which(data_day$spike==T)
  
  p1 <- ggplot() + 
    geom_line(data=data_day, aes(date, temp), colour="blue", size=1) +
    geom_point(data=data_day, aes(date, temp), colour="red", size=1.5) +
    geom_vline(xintercept=as.numeric(data_day$date[data_day$spike]), size=.5, colour="darkgreen") +
    scale_x_date(date_labels = "%Y", date_breaks = "2 year") + xlab("") + 
    ylab("Temperature") +
    ggtitle(site_name)
  
  if (length(sample_spike) >= 1){
    p2 <- ggplot() + 
      geom_line(data=data_day[(sample_spike[1]-10):(sample_spike[1]+10),], aes(date, temp), colour="blue", size=1) +
      geom_point(data=data_day[(sample_spike[1]-10):(sample_spike[1]+10),], aes(date, temp), colour="red", size=1.5) +
      geom_vline(xintercept=as.numeric(data_day$date[data_day$spike]), size=.5, colour="darkgreen") +
      scale_x_date(date_labels = "%Y", date_breaks = "2 year") + xlab("") + 
      ylab("Temperature") +
      ggtitle("Temperature")
    
    plot(arrangeGrob(grobs=list(p1,p2), ncol=4, layout_matrix=cbind(1,1,1,2)))
  } else {
    print(p1)
  }
}
