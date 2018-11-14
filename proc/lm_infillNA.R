# Fill in NA intelligently #------------------------------------------------------------------------
# fill in missing values using lm with stations that have r>0.8
# find the longest non-NA segement
# each matching station is used individually to fill in missing data to
#   maximize filled in values
library(tidyverse)
library(lme4)
library(forecast)
library(corrplot)

for (j in 1:length(ind.grp)){
  stat.name <- as.character(site_list$index[ind.grp[j]])
  #stat.name <- rlang::sym(stat.name)
  
  grp.stat <- as.vector(site_list$index[ind.grp])   # get the station names in the group
  grp.stats <- df.data.stor %>%                     # spread each station as a row
    filter(index %in% grp.stat) %>%
    spread(index, temp)
  rm(grp.stat)
  
  # trim the df for the selected station (first and last non-NA)
  query.stat <- rlang::sym(stat.name)
  query.ind <- which(colnames(grp.stats) == query.stat) # first column in date
  grp.stats <- ts.trim(grp.stats, query.ind)            # returns trimmed df
  
  # strip down to the longest section of non-NA
  ind.long <- strip.na(grp.stats[,query.ind, drop=T])     # returns indices to longest non-NA
  grp.sect <- grp.stats[ind.long,]
  
  grp.cor <- cor(grp.sect[,-1], use = "pairwise.complete.obs")  # trimmed to longest section
  corrplot(grp.cor)
  
  # get the strongest correlation stations r > 0.8 (1 is always self)
  stat.corr <- sort(grp.cor[,(query.ind-1)], decreasing = T)
  stat.refs <- stat.corr[2:length(which(stat.corr > .8))]       # returns site and r value
  
  # create time series of original data
  doy <- format(grp.stats$date[1], "%j")
  yr <- format(grp.stats$date[1], "%Y")
  
  # create time series of station in question
  y.in <- as.ts(grp.stats[,query.ind], frequency=365)[,1]
  y.na <- which(is.na(y.in))                                    # record positions of NA in station ts
  x.in <- as.ts(grp.stats[names(stat.refs)], frequency=365)
  #names(x.in) <- stat.refs
  
  # Predict using full time series with linear mixed model ------------------
  # Use the full time series to determine monthly coefficients and adjust predicted
  #   values by the difference in medians of the observed - predicted values
  y.med <- grp.stats %>% 
    dplyr::select(date, stat.name) %>% 
    mutate(y.in = .[[2]], month = format(date, "%m"), year = format(date, "%Y")) %>% 
    group_by(month, year) %>% 
    summarise(month.med = median(y.in, na.rm=T), na.values = sum(is.na(y.in))) %>% 
    arrange(year, month) %>% 
    mutate(date = as.character(paste(month, year, sep = "-"))) %>% 
    ungroup()
  
  y.med$month.med.update <- mapply(function(x, y) ifelse(x >= 20, NA, y),
                                   y.med$na.values, y.med$month.med)
  
  # interpolate missing medians using seasonal information
  y.med.ts <- ts(y.med$month.med.update, frequency = 12)
  y.med.interp <- na.interp(y.med.ts)
  y.med$y.med.interp <- y.med.interp
  
  filename1 <- paste0("~/Desktop/StationInfill/",
                      gsub(gsub(stat.name, pattern=" ", replacement=""), pattern = "/", replacement = "_"),
                      "_MedMonthInterp.png")
  par(mar=c(1,1,1,1))
  png(filename1, width=18, height=12, units="cm", res=300)
  
  ind.st <- which(y.med$month == "01")[1]
  
  plot(as.numeric(y.med.interp), type = "l", col = "red", lwd = 2,
       main = paste("Month Median for", stat.name), xaxt = "n",
       ylab = "Temperature (°C)", xlab = "Time (months)")
  axis(1,at=seq(ind.st,length(y.med.interp),12),labels=y.med$date[seq(ind.st,length(y.med.interp),12)])
  abline(v = seq(ind.st,length(y.med.interp),12))
  legend("topleft", legend = c("Original", "Interpolated", "NA > 20"), 
         inset=c(0,-.1), bty = "n", xpd=TRUE, mar(c(3,3,3,3)),
         cex = .8, pch=19, col=c("gray30","red","gray70"), horiz = T)
  lines(y.med$month.med, col="black", lwd=2)
  points(y.med$month.med, col = "gray70", pch = 19, cex = 1)
  points(as.numeric(y.med.interp), col = "red", pch = 19, cex = 1)
  points(y.med$month.med.update, col = "gray30", pch = 19, cex = 1)
  # abline(v = 1:nrow(y.med), col="gray80")
  # abline(v = which(y.med$month=="01"), col = "darkorange")
  # abline(v = which(y.med$month=="07"), col = "blue")
  
  dev.off()
  
  # create empty data frame to store interpolated data
  interp.stor <- data.frame(numeric(nrow(x.in)), stringsAsFactors = F)
  
  for (i in 1:ncol(x.in)){
    ref.stat <- as.numeric(x.in[,i])
    ref.stat <- data.frame("ref.stat"=ref.stat, "month"=format(grp.stats$date, "%m"))
    # normalize using 90th percentile
    #ref.stat.norm <- normalit.quan(ref.stat)
    #y.in.norm <- normalit.quan(y.in)
    
    # 'naive' lm impute
    df.in.train <- ref.stat %>% 
      mutate(y.in=y.in)
    
    #names(df.in.train) <- c("y.in", "ref.stat", "month")
    #fit.reg <- lm(y.in ~ ., data = df.in.train)
    fit.reg <- lmer(y.in ~ ref.stat + (1|month), data = df.in.train)
    
    df.in.pred <- df.in.train %>% 
      dplyr::select(ref.stat, month)
    #names(df.in.pred) <- ref.name
    pred.reg <- predict(fit.reg, newdata = df.in.pred, allow.new.levels = T)
    
    # create a monthly median from the observation data andl linearly interpolate missing months
    #   with too few data
    cnt = 1
    cnt.add <- 29       # cnt + cnt.add + 1 is the start of the next interation
    pred.adj <- vector(mode = "numeric", length = length(y.in))
    while ((cnt + cnt.add) < (length(y.in) + cnt.add)){   # if == to length(y.in) then already done
      if (cnt + cnt.add > length(y.in)){      # if longer than y.in shorten cnt.add to end at length y.in
        cnt.add <- length(y.in) - cnt
      } 
      # if there is not enough data for median in observations use interpolated median
      if (sum(is.na(y.in[cnt:(cnt+cnt.add)])) > 20){
        interp.dates <- grp.stats$date[cnt:(cnt+cnt.add)]
        interp.monyear <- c(format(interp.dates, "%m")[15],format(interp.dates, "%Y")[15])
        interp.med <- y.med %>% 
          dplyr::filter(., month==interp.monyear[1] & year==interp.monyear[2]) %>% 
          select(y.med.interp)
        diff <- as.numeric(interp.med) - median(pred.reg[cnt:(cnt+cnt.add)], na.rm = T)
        pred.adj[cnt:(cnt+cnt.add)] <- pred.reg[cnt:(cnt+29)] + diff
        print(paste("Too few original data for median, using interpolated for",
                    interp.monyear[1],interp.monyear[2]))
        cnt = cnt + cnt.add + 1
      } else {
        diff <- median(y.in[cnt:(cnt+cnt.add)], na.rm = T) - median(pred.reg[cnt:(cnt+cnt.add)], na.rm = T)
        pred.adj[cnt:(cnt+cnt.add)] <- pred.reg[cnt:(cnt+cnt.add)] + diff
        lines(pred.adj, col="blue")
        cnt = cnt + cnt.add + 1
      }
    }
    
    pred.adj <- unlist(pred.adj)
    interp.stor <- cbind(interp.stor, pred.adj)
    
    err1 <- sqrt(mean((y.in - ref.stat[,1])^2, na.rm = T))
    err2 <- sqrt(mean((y.in - pred.reg)^2, na.rm = T))
    err3 <- sqrt(mean((y.in - pred.adj)^2, na.rm = T))
    
    ref.name <- c(dimnames(x.in)[[2]][i],"month")
    
    filename1 <- paste0("~/Desktop/StationInfill/",
                        gsub(gsub(stat.name, pattern=" ", replacement=""), pattern = "/", replacement = "_"),
                        "_AdjLMM", i, ".png")
    par(mar=c(1,1,1,1))
    png(filename1, width=18, height=12, units="cm", res=300)
    
    ind.st <- which(format(grp.stats$date, "%m") == "01" & format(grp.stats$date, "%d") == "01")
    text.pos <- length(y.in)/5.5
    
    plot(y.in, type="l", lwd=2, main=paste(stat.name,"with ref",ref.name[1]),
         ylab = "Temperature (°C)", xlab="Time", xaxt = "n")
    axis(1,at = ind.st, labels = format(grp.stats$date[ind.st], "%Y"))
    abline(v = ind.st)
    legend("topleft", legend = c("Original", paste("Reference (",round(err1,2),")"),
                                 paste("LMM (",round(err2,2),")"),
                                 paste("Adjusted LMM (",round(err3,2),")")), 
           inset=c(0,-.1), bty = "n", xpd=NA, cex = .8, lty=1, lwd=3, xjust = 0,
           col=c("black","blue","chartreuse4","darkorange"),
           x.intersp = 0.5, horiz = T,  text.width=c(1,text.pos,text.pos,text.pos))
    lines(ref.stat[,1], col="blue")
    lines(pred.reg, col = "chartreuse4")
    lines(pred.adj, col = "darkorange")
    
    dev.off()
  }
  
  # replace first column with observed data
  interp.stor[,1] <- y.in
  # rename arbitarily for unique names
  names(interp.stor) <- as.character(1:ncol(interp.stor))
  y.in.interp <- interp.stor %>% 
    mutate(new = do.call(coalesce, .)) %>% 
    select(new)
  
  filename1 <- paste0("~/Desktop/StationInfill/",
                      gsub(gsub(stat.name, pattern=" ", replacement=""), pattern = "/", replacement = "_"),
                      "_Interpolated", i, ".png")
  
  png(filename1, width=18, height=18, units="cm", res=300)
  #par(mar=c(4,4,3,3), mfrow=c(2,1))
  par(mfrow=c(2,1))
  par(mar=c(2,4,3,3))
  na.n1 <- sum(is.na(y.in))
  plot(y.in, type="l", lwd=2, 
       main=paste(stat.name,"Original (NA =",na.n1,")"),
       ylab = "Temperature (°C)", xlab="", xaxt = "n")
  abline(v = ind.st)
  par(mar=c(4,4,2,3))
  na.n2 <- sum(is.na(y.in.interp))
  plot(y.in.interp, type="l", lwd=2, 
       main=paste(stat.name,"with Interpolated Values (NA = ",na.n2,")"),
       ylab = "Temperature (°C)", xlab="Time", xaxt = "n")
  abline(v = ind.st)
  axis(1,at = ind.st, labels = format(grp.stats$date[ind.st], "%Y"))
  
  dev.off()
}
