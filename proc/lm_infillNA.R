# Fill in NA intelligently #------------------------------------------------------------------------
# fill in missing values using lm with stations that have r>0.8
# find the longest non-NA segement
# each matching station is used individually to fill in missing data to
#   maximize filled in values

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

# for each station in x.in which will be in order
# set window size for training and forward predicting
size.window <- 120
size.interp <- 30

for (i in 1:ncol(x.in)){
  ref.stat <- x.in[,i]
  ref.name <- dimnames(x.in)[[2]][i]
  
  plot(y.in, type="l", lwd=2, main=paste(stat.name,"with ref",ref.name),
       ylab = "Temperature (°C)")
  lines(ref.stat, col="blue")

}


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
