# Fill in NA intelligently #------------------------------------------------------------------------
# fill in missing values using lm with stations that have r>0.8
# find the longest non-NA segement
# each matching station is used individually to fill in missing data to
#   maximize filled in values
library(tidyverse)
library(lme4)
library(forecast)

# normalize data columns using quantiles as min / max (does not work due to apparent seasonal
#   $ temperature dependency)
normalit.quan<-function(m){
  minmax <- quantile(m, c(.1,.9), na.rm=T)
  min <- minmax[1]
  max <- minmax[2]
  (m - min(m, na.rm=T))/(max(m, na.rm=T)-min(m, na.rm=T))
}

stat.name <- as.character(site_list$index[ind.grp[2]])
stat.name <- rlang::sym(stat.name)

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
  arrange(year, month)

y.med$month.med.update <- mapply(function(x, y) ifelse(x >= 20, NA, y),
                          y.med$na.values, y.med$month.med)

# interpolate missing medians using seasonal information
y.med.ts <- ts(y.med$month.med.update, frequency = 12)
y.med.interp <- na.interp(y.med.ts)

plot(as.numeric(y.med.interp), type = "l", col = "red", lwd = 2,
     main = paste("Month Median for", stat.name),
     ylab = "Temperature (°C)", xlab = "Time (months)")
legend("topleft", legend = c("Original", "Interpolated", "NA > 20"), 
       inset=c(0,-.1), bty = "n", xpd=TRUE, mar(c(3,3,3,3)),
       cex = 1, pch=19, col=c("gray30","red","gray70"), horiz = T)
lines(y.med$month.med, col="black", lwd=2)
points(y.med$month.med, col = "gray70", pch = 19, cex = 1.5)
points(as.numeric(y.med.interp), col = "red", pch = 19, cex = 1.5)
points(y.med$month.med.update, col = "gray30", pch = 19, cex = 1.5)
abline(v = 1:nrow(y.med), col="gray80")
abline(v = which(y.med$month=="01"), col = "darkorange")
abline(v = which(y.med$month=="07"), col = "blue")

for (i in 1:ncol(x.in)){
  ref.stat <- as.numeric(x.in[,i])
  ref.stat <- data.frame("ref.stat"=ref.stat, "month"=format(grp.stats$date, "%m"))
  # normalize using 90th percentile
  #ref.stat.norm <- normalit.quan(ref.stat)
  #y.in.norm <- normalit.quan(y.in)
  
  ref.name <- c(dimnames(x.in)[[2]][i],"month")
  
  plot(y.in, type="l", lwd=2, main=paste(stat.name,"with ref",ref.name[1]),
       ylab = "Temperature (°C)")
  legend("topleft", legend = c("Original", "Reference", "Naive LM"), 
         inset=c(0,-.1), bty = "n", xpd=TRUE, mar(c(3,3,3,3)),
         cex = 1, lty=1, lwd=3, col=c("black","blue","red"), horiz = T)
  lines(ref.stat[,1], col="blue")
  
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
  pred.adj <- vector(mode = "numeric", length = length(y.in))
  while ((cnt + 30) < length(y.in)){
    # if there is not enough data for median in observations use interpolated median
    
    diff <- median(y.in[cnt:(cnt+29)], na.rm = T) - median(pred.reg[cnt:(cnt+29)], na.rm = T)
    pred.adj[cnt:(cnt+29)] <- pred.reg[cnt:(cnt+29)] - diff
    cnt = cnt + 30
  }
 
  lines(pred.reg, col = "purple")
  lines(pred.adj, col = "darkorange")

  
  
# Train 60 predict month 1 year  ahead ------------------------------------
#   this fails with too few records in the 60 day interval

  cnt = 1
  interp.forward <- vector(mode = "numeric", length = length(y.in))
  while (cnt + 379 < length(y.in)){
    y.in.sub <- y.in[cnt:(cnt+59)]
    x.in.sub <- ref.stat[cnt:(cnt+59),]
    
    if (sum(is.na(x.in.sub$ref.stat))/60 > 0.6 | sum(is.na(y.in.sub))/60 > 0.6){
      print("Too few")
      cnt <- cnt + 30
      next
    }
    
    df.in.train <- data.frame(cbind(y.in.sub, x.in.sub$ref.stat))
    names(df.in.train) <- c("y.in.sub", "ref.stat.sub")
    df.in.pred <- data.frame(ref.stat$ref.stat[(cnt+350):(cnt+379)])
    names(df.in.pred) <- "ref.stat.sub"
    
    fit.reg <- lm(y.in.sub ~ ref.stat.sub, data = df.in.train)
    pred.reg <- predict.lm(fit.reg, newdata = df.in.pred)
    interp.forward[(cnt+350):(cnt+379)] <- pred.reg
    
    cnt <- cnt + 30
  }
  

# Similar to above --------------------------------------------------------
# Predicts moving forward and the predicts moving backwards and looks for
# Largest errors between the two predicted graphs
# Again runs into problems with the number of observations for training
  
  # for each station in x.in which will be in order
  # set window size for training and forward predicting
  size.window <- 60     
  size.interp <- 30
  size.gap <- 350  
  
  # function to interpolate using sliding LM forward/backward, using 4 month window
  #   to predict 30 day ahead
  interp.forback <- function(y.in, ref.stat){
    # create empty vector storing predictions
    interp.forward <- vector(mode = "numeric", length = length(y.in))
    interp.backward <- interp.forward
    
    cnt.st <- 1
    cnt.en <- cnt.st + size.window -1
    
    while (cnt.en < length(y.in)){
      if ((cnt.en + size.interp) < (length(y.in)-(2*size.interp))){
        y.in.sub <- y.in[cnt.st:cnt.en]
        x.in.sub <- ref.stat[cnt.st:cnt.en,]
        df.in.train <- data.frame(cbind(y.in.sub, x.in.sub))
        df.in.pred <- data.frame(ref.stat[cnt.en:(cnt.en + size.interp -1),])
        names(df.in.pred) <- names(df.in.train)[-1]
        
        # if there is insufficient data use previous lm
        if (sum(is.na(df.in.train$y.in.sub)) / nrow(df.in.train) < 0.66){
          fit.reg <- lmer(y.in.sub ~ ref.stat + (1|month), data = df.in.train)
        } else {
          print("Reverting to previous LM")
        }
        
        # if the first interation, use model prediction on first 120 days
        if (cnt.st == 1){
          pred.first <- predict(fit.reg,
                                   newdata=data.frame(df.in.train[,-1]),
                                   allow.new.levels = T)
          interp.forward[cnt.st:cnt.en] <- pred.first
          rm(pred.first)
          cnt.st <- cnt.st + size.interp
          cnt.en <- cnt.st + size.window -1
          print("First interation complete")
          next
        }
        pred.reg <- predict(fit.reg, newdata = df.in.pred,
                            allow.new.levels = T)
        interp.forward[(cnt.en+1):(cnt.en + size.interp)] <- pred.reg
        
        cnt.st <- cnt.st + size.interp
        cnt.en <- cnt.st + size.window -1
        
        # to deal with the last segment variable length
      } else {
        # if 31 to 60 after 120 training data then extend predict to end
        if (length(y.in) - cnt.en - size.interp < 2*size.interp){
          size.interp.tmp <- length(y.in) - cnt.en
          
          y.in.sub <- y.in[cnt.st:cnt.en]
          x.in.sub <- ref.stat[cnt.st:cnt.en,]
          df.in.train <- data.frame(cbind(y.in.sub, x.in.sub))
          df.in.pred <- data.frame(ref.stat[(cnt.en+1):(cnt.en + size.interp.tmp -1),])
          #names(df.in.pred) <- names(df.in.train)[-1]
          
          fit.reg <- lmer(y.in.sub ~ ref.stat + (1|month), data = df.in.train)
          pred.reg <- predict(fit.reg, newdata = df.in.pred,
                                 allow.new.levels = T)
          interp.forward[(cnt.en+1):(cnt.en + size.interp.tmp -1)] <- pred.reg
          print("Reached the end of the series")
          break
          
          # is less than 120 then self predict
        } else if (length(y.in) - cnt.en <= 120){
          y.in.sub <- y.in[cnt.st:length(y.in)]
          x.in.sub <- ref.stat[cnt.st:length(y.in)]
          df.in.train <- data.frame(cbind(y.in.sub, x.in.sub))
          pred.last <- predict(fit.reg, newdata=data.frame(df.in.train[,-1]),
                               allow.new.levels = T)
          interp.forward[cnt.en:length(y.in)] <- pred.last
          rm(pred.last)
          print("Reached the end of the series")
          break
        }
      }
      # reverse and run backward
      interp.forward.org <- interp.forward
      y.in.org <- y.in
      y.in <- rev(y.in.org)
      ref.stat.org <- ref.stat
      ref.stat <- rev(ref.stat)
      
      # reset and plot
      interp.backward <- rev(interp.forward)
      y.in <- y.in.org
      ref.stat <- ref.stat.org
      interp.forward <- interp.forward.org
      
      plot(y.in, type="l")
      lines(interp.forward, col="darkorange")
      lines(interp.backward, col="chartreuse4")
    }
  }
}

plot(y.in, type="l")
lines(interp.forward, col="darkorange")
lines(interp.backward, col="chartreuse4")


# Unknown -----------------------------------------------------------------


# # find which time series are best for filling in missing data
# if (ncol(x.in)>1){
#   pair.comb <- combn(1:ncol(x.in),2,simplify=FALSE)
#   # get the highest number of complete observations where y is NA
#   best.comb <- sapply(pair.comb, function(x) length(na.omit(x.in[y.na, x])))
#   if (length(best.comb)==1 & best.comb == 0){        # if there are no overlapping data in regressors
#     # evaluate the regrssors individually
#     best.comb <- sapply(1:ncol(x.in), function(x) length(na.omit(x.in[y.na, x])))
#     best.comb <- which(best.comb==max(best.comb))
#     x.in <- x.in[,best.comb]
#   } else {
#     best.comb <- which(best.comb==max(best.comb))
#     if (length(best.comb)>1){
#       best.comb <- best.comb[1]
#     }
#     x.in <- x.in[,pair.comb[[best.comb]]]
#   }
# }

# get the detail of the regression station for plot
corr.stats <- stat.refs[names(stat.refs) %in% colnames(x.in)]
subtitle <- sapply(1:length(corr.stats),
                   function(x) paste0(names(corr.stats[x]),
                                      " (r=",round(as.vector(corr.stats[x]),2),")"))
subtitle <- str_c(subtitle, collapse=" ")


df.in.train <- data.frame(cbind(y.in, x.in))
fit.reg <- lm(y.in ~ ., data = df.in.train)
df.in.pred <- data.frame(df.in.train[,-1])
names(df.in.pred) <- names(df.in.train)[-1]
pred.reg <- predict(fit.reg, newdata = df.in.pred)
pred.reg <- ts(pred.reg, frequency=365,
               start = c(as.numeric(yr),as.numeric(doy)))

# fit.reg <- auto.arima(y.in, xreg=cbind(x.in), seasonal=F,
#                       allowdrift = F)
# pred.reg <- ts(predict(fit.reg,
#                        newxreg = cbind(x.in))[[1]], frequency=365,
#                start = c(as.numeric(yr),as.numeric(doy)))
