library(ggplot2)
library(gridExtra)
library(dplyr)
library(changepoint)
library(changepoint.np)
library(nortest)
library(qcc)

load("~/Data/SACTNdaily_v4.1.Rdata")
load("~/Data/SACTNhourly_v4.1.Rdata")

dataMB_day_DEA <- subset(SACTNdaily_v4.1, site=="Mossel Bay" & src=="DEA")
dataMB_day_SAWS <- subset(SACTNdaily_v4.1, site=="Mossel Bay" & src=="SAWS")
dataMB_day_DEA$SAWS <- dataMB_day_SAWS$temp[match(dataMB_day_DEA$date,dataMB_day_SAWS$date)]    # SAWS data matching DEA

dataMB_day_DEA$diff <- dataMB_day_DEA$SAWS-dataMB_day_DEA$temp

dataMB_day_DEA$year <- format(dataMB_day_DEA$date,"%Y")
dataMB_day_DEA$month <- format(dataMB_day_DEA$date,"%b")

dataMB_day_SAWS$year <- format(dataMB_day_SAWS$date,"%Y")
dataMB_day_SAWS$month <- format(dataMB_day_SAWS$date,"%b")

dataMB_day_DEA <- dataMB_day_DEA %>% mutate(grad=temp-lag(temp))

graddetrend_DEA <- dataMB_day_DEA$seastrend[2:nrow(dataMB_day_DEA)] - dataMB_day_DEA$seastrend[1:nrow(dataMB_day_DEA)-1]
dataMB_day_DEA$graddetrend <- NA
dataMB_day_DEA$graddetrend[2:nrow(dataMB_day_DEA)] <- graddetrend_DEA

dataMB_day_SAWS <- dataMB_day_SAWS %>% mutate(grad=temp-lag(temp))

# get doy from first date
doy <- format(dataMB_day_DEA$date[1], "%j")
yr <- format(dataMB_day_DEA$date[1], "%Y")
ts_DEA <- ts(dataMB_day_DEA$temp, frequency=365, start = c(as.numeric(yr),
                                                           as.numeric(doy)))
decmp_DEA <- stlplus(ts_DEA, s.window="periodic")[[1]]
dataMB_day_DEA$res <- decmp_DEA$remainder
dataMB_day_DEA$raw <- decmp_DEA$raw
dataMB_day_DEA$seas <- decmp_DEA$seasonal
dataMB_day_DEA$trend <- decmp_DEA$trend
dataMB_day_DEA$noseas <- decmp_DEA$remainder + decmp_DEA$trend
dataMB_day_DEA$notrend <- decmp_DEA$raw - decmp_DEA$trend
dataMB_day_DEA$NA_idx <- is.na(dataMB_day_DEA$temp)

# to fill in NA with seasonal

ggplot(data=dataMB_day_DEA) +
  geom_line(aes(x=date, y=raw, colour="raw"), size=1) +
  geom_line(aes(x=date, y=notrend, colour="notrend"), size=1) +
  geom_line(aes(x=date, y=noseas, colour="noseas"), size=.5) +
  geom_line(aes(x=date, y=trend, colour="trend"), size=1) +
  geom_line(aes(x=date, y=seas, colour="seas")) +
  geom_line(aes(x=date, y=res, colour="res"), size=.5) +
  scale_colour_manual("",
                      breaks=c("raw","notrend","noseas","trend","seas","res"),
                      values=c("gray","blue","darkorange","red","green","black"))

ggplot(aes(x = seastrend, fill = month), data = dataMB_day_DEA) +
  geom_histogram(binwidth=.5) +
  facet_grid(month ~ ., margins = FALSE, scales = "free") +
  ggtitle("Mossel Bay DEA Temp.") +
  scale_fill_discrete(guide=FALSE)

ggplot(aes(x = graddetrend, fill = month), data = dataMB_day_DEA) +
  geom_histogram(binwidth=.1) +
  xlim(-3.5,3.5) +
  facet_grid(month ~ ., margins = FALSE, scales = "free") +
  ggtitle("Mossel Bay DEA Temp.") +
  scale_fill_discrete(guide=FALSE)

grad_resid_DEA <- dataMB_day_DEA$res[2:nrow(dataMB_day_DEA)] - dataMB_day_DEA$res[1:nrow(dataMB_day_DEA)-1]
dataMB_day_DEA$grad_res <- NA
dataMB_day_DEA$grad_res[2:nrow(dataMB_day_DEA)] <- grad_resid_DEA

doy <- format(dataMB_day_SAWS$date[1], "%j")
yr <- format(dataMB_day_SAWS$date[1], "%Y")
ts_SAWS <- ts(dataMB_day_SAWS$temp, frequency=365, start = c(as.numeric(yr),as.numeric(doy)))
decmp_SAWS <- stlplus(ts_SAWS, s.window="periodic")[[1]]
dataMB_day_SAWS$res <- decmp_SAWS$remainder

grad_resid_SAWS <- dataMB_day_SAWS$res[2:nrow(dataMB_day_SAWS)] - dataMB_day_SAWS$res[1:nrow(dataMB_day_SAWS)-1]
dataMB_day_SAWS$grad_res <- NA
dataMB_day_SAWS$grad_res[2:nrow(dataMB_day_SAWS)] <- grad_resid_SAWS

write.csv(dataMB_day_DEA, "~/R/Docs/dataMB_day_DEA.csv")
write.csv(dataMB_day_SAWS, "~/R/Docs/dataMB_day_SAWS.csv")

tempUwind <- dataMB_wind$Uwnd[match(dataMB_day_DEA$date, as.Date(dataMB_wind$datetime))]
tempVwind <- dataMB_wind$Vwnd[match(dataMB_day_DEA$date, as.Date(dataMB_wind$datetime))]
dataMB_day_DEA$Uwnd <- tempUwind
dataMB_day_DEA$Vwnd <- tempVwind
rm(tempUwind, tempVwind)

#________________________________________________________________________________________
ggplot() + 
  geom_line(data=dataMB_day_DEA, aes(date, temp, color='DEA'), size=1.5) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") + xlab("") + ylab("Tempearture") +
  geom_line(data=dataMB_day_SAWS, aes(date, temp, color='SAWS'), size=1.5, alpha=0.5) +
  geom_line(data=dataMB_day_DEA, aes(date, diff, color='Difference'), size=1.5) +
  scale_colour_manual("", 
                    breaks = c("DEA", "SAWS", "Difference"),
                    values = c("DEA"="blue", "SAWS"="red", 
                               "Difference"="darkgreen"))

# ________________________________________________________________________________________________________________
qplot(dataMB_day_DEA$temp, geom = "histogram")
qplot(dataMB_day_SAWS$temp, geom = "histogram")

# ________________________________________________________________________________________________________________
seas <- list(c("12","01","02"),c("03","04","05"),c("06","07","08"),c("09","10","11"))

dataMB_day_DEAsub <- lapply(seas, function(x) subset(SACTNdaily_v4.1, site=="Mossel Bay" & src=="DEA" &
                                                    substr(date, 6,7) %in% x))
dataMB_day_SAWSsub <- lapply(seas, function(x) subset(SACTNdaily_v4.1, site=="Mossel Bay" & src=="SAWS" &
                                                     substr(date, 6,7) %in% x))

p1 <- lapply(dataMB_day_DEAsub, function(x) qplot(x$temp, geom = "histogram", xlab="Tenperature",
                                               main="Mossel Bay DEA", xlim=c(10,25), asp=1))
p2 <- lapply(dataMB_day_SAWSsub, function(x) qplot(x$temp, geom = "histogram", xlab = "Temperature",
                                                main="Mossel Bay SAWS",xlim=c(10,25),asp=1))
do.call(grid.arrange,c(p1, ncol=4))
do.call(grid.arrange,c(p2, ncol=4))

# how many numbers end on x.5 and y.0
revtrunc <- function(x) { x - floor(x) }
check <- revtrunc(dataMB_day_SAWS$temp)
qplot(check, geom="bar", binwidth=0.5, xlim=c(-0.25,0.75))

seas_bias <- lapply(dataMB_day_SAWSsub, function(x) revtrunc(x$temp))
p3 <- lapply(seas_bias, function(x) qplot(x, geom="histogram", main="SAWS decimal bias",
                                          binwidth=0.5, xlim=c(-0.25,0.75)))
do.call(grid.arrange, c(p3, ncol=4))
#_______________________________________________________________________________________-
table1 <- data.frame(mean=numeric(),median=numeric(), std=numeric(),stringsAsFactors = F)
t1 <- c(mean(dataMB_day_DEA$temp, na.rm = T), 
                  median(dataMB_day_DEA$temp,na.rm = T), 
                  sd(dataMB_day_DEA$temp,na.rm = T))
t2 <- c(mean(dataMB_day_SAWS$temp, na.rm = T), 
        median(dataMB_day_SAWS$temp,na.rm = T), 
        sd(dataMB_day_SAWS$temp,na.rm = T))
table1[nrow(table1)+1,] <- t1
table1[nrow(table1)+1,] <- t2
row.names(table1) <- c("Mossel Bay-DEA","Mossel Bay-SAWS")

#________________________________________________________________________________________-
stat_DEA <- dataMB_day_DEA %>% group_by(month,year) %>% summarise(mean = mean(temp, na.rm=T),
                                                                      sd = sd(temp, na.rm=T))
stat_SAWS <- dataMB_day_DEA %>% group_by(month,year) %>% summarise(mean = mean(SAWS, na.rm=T),
                                                                  sd = sd(SAWS, na.rm=T))

stat_Uwnd <- dataMB_day_DEA %>% group_by(month,year) %>% summarise(mean = mean(Uwnd, na.rm=T),
                                                                   sd = sd(Uwnd, na.rm=T))
stat_Vwnd <- dataMB_day_DEA %>% group_by(month,year) %>% summarise(mean = mean(Vwnd, na.rm=T),
                                                                   sd = sd(Vwnd, na.rm=T))

ggplot(data = stat_DEA, aes(y=mean, x=year, group=month)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=.2) +
  facet_grid(month ~ ., margins = TRUE, scales = "free") +
  ggtitle("Mossel Bay DEA Temp.") +
  scale_fill_discrete(guide=FALSE)

ggplot(data = stat_SAWS, aes(y=mean, x=year, group=month)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=.2) +
  facet_grid(month ~ ., margins = TRUE, scales = "free") +
  ggtitle("Mossel Bay SAWS Temp.") +
  scale_fill_discrete(guide=FALSE)

ggplot(data = stat_Uwnd, aes(y=mean, x=year, group=month)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=.2) +
  facet_grid(month ~ ., margins = TRUE, scales = "free") +
  ggtitle("Mossel Bay Uwind.") +
  scale_fill_discrete(guide=FALSE)

#######################################################################################
# Determine gradients on seasonal
ggplot(aes(x = grad, fill = factor(month)), data = dataMB_day_DEA) +
  geom_histogram(binwidth=.05) +
  xlim(-4,4) +
  geom_vline(aes(xintercept = 0)) +
  facet_grid(month ~ ., margins = FALSE, scales = "free") +
  ggtitle("Mossel Bay DEA Temp. Grad.") +
  scale_fill_discrete(guide=FALSE)

ggplot(aes(x = grad, fill = factor(month)), data = dataMB_day_SAWS) +
  geom_histogram(binwidth=.05) +
  xlim(-4,4) +
  geom_vline(aes(xintercept = 0)) +
  facet_grid(month ~ ., margins = FALSE, scales = "free") +
  ggtitle("Mossel Bay SAWS Temp. Grad.") +
  scale_fill_discrete(guide=FALSE)

# Determine gradients on residuals
ggplot(aes(x = grad_res, fill = factor(month)), data = dataMB_day_DEA) +
  geom_histogram(binwidth=.05) +
  geom_vline(aes(xintercept = 0)) +
  facet_grid(month ~ ., margins = FALSE, scales = "free") +
  ggtitle("Mossel Bay DEA Res. Grad.") +
  scale_fill_discrete(guide=FALSE)

ggplot(aes(x = grad_res, fill = factor(month)), data = dataMB_day_SAWS) +
  geom_histogram(binwidth=.05) +
  geom_vline(aes(xintercept = 0)) +
  facet_grid(month ~ ., margins = FALSE, scales = "free") +
  ggtitle("Mossel Bay SAWS Res. Grad.") +
  scale_fill_discrete(guide=FALSE)
#######################################################################################
library(nortest)
# fails if month contains lest than 7 values

# function to test for normality
normtest <- function(x){
  if (length(x) < 20){
    return(NA)
  } else if (length(x) >= 20) {
    a <- ad.test(x)
    if (a[[2]] >= 0.5){
      return(TRUE)
    } else if (a[[2]] < 0.5) {
      return(FALSE)
    }
  }
}

norm_DEA <- dataMB_day_DEA %>% group_by(month,year) %>% summarise(test = normtest(temp))

########################################################################################
# Detection of outliers

percres_ls <- split(stat_DEA_percres, seq(nrow(stat_DEA_percres)))

df <- dataMB_day_DEA
mon <- percres_ls[[1]][[1]]
hi <- percres_ls[[1]][[2]]
lo <- percres_ls[[1]][[3]]

df$outlier <- TRUE
getOutlier <- function(df, mon, hi, lo){
  idx <- which(df$month == mon & df$grad_res > lo & df$grad_res < hi)
}

outlier_ls <- lapply(percres_ls, function(x) getOutlier(df, x[[1]], x[[2]], x[[3]]))

for (i in 1:12){
  df$outlier[outlier_ls[[i]]] <- FALSE
}

ggplot(aes(x = seastrend, fill = month), data = dataMB_day_DEA) +
  geom_histogram(binwidth=.5) +
  facet_grid(month ~ ., margins = FALSE, scales = "free") +
  ggtitle("Mossel Bay DEA Temp.") +
  scale_fill_discrete(guide=FALSE)

#######################################################################################
# Change point detection
data_cp_res <- decmp_DEA$remainder[-which(is.na(decmp_DEA$raw))]
data_cp_raw <- decmp_DEA$raw[-which(is.na(decmp_DEA$raw))]
data_cp_noseas <- (decmp_DEA$raw - decmp_DEA$seasonal)[-which(is.na(decmp_DEA$raw))]

plot(data_cp_res, type="l")
plot(data_cp_noseas, type="l")

ts_cp_res <- ts(data_cp_res, frequency=365, start = c(1999, 1))
ts_cp_noseas <- ts(data_cp_noseas, frequency=365, start = c(1999, 1))

hist(data_cp_noseas, breaks=seq(10,28,0.2), main="Mossel Bay-DEA Non-seasonal")
ad.test(data_cp_noseas)

cp_method <- "BinSeg"
cp_mean01 <- cpt.mean(ts_cp_noseas, method=cp_method)
cp_var01 <- cpt.var(ts_cp_noseas, method=cp_method)
cp_meanvar01 <- cpt.meanvar(ts_cp_noseas, method=cp_method)

plot(cp_mean01)
plot(cp_var01)
plot(cp_meanvar01)

cp_method <- "SegNeigh"
cp_mean03 <- cpt.mean(ts_cp_noseas, penalty="None", method=cp_method, Q=5)
cp_var03 <- cpt.var(ts_cp_noseas, penalty="None", method=cp_method, Q=5)
cp_meanvar03 <- cpt.meanvar(ts_cp_noseas, penalty="None", method=cp_method, Q=5)

plot(cp_mean03, main="Mean CP-No Pen")
plot(cp_var03, main="Var CP-No Pen")
plot(cp_meanvar03, main="Mean + Var CP-No Pen")

cp_method <- "SegNeigh"
cp_mean04 <- cpt.mean(ts_cp_noseas, penalty="Manual", pen.value='1.5*log(n)',method=cp_method, Q=5)
cp_var04 <- cpt.var(ts_cp_noseas, penalty="Manual", pen.value='1.5*log(n)', method=cp_method, Q=5)
cp_meanvar04 <- cpt.meanvar(ts_cp_noseas, penalty="Manual", pen.value='1.5*log(n)', method=cp_method, Q=5)

plot(cp_mean04, main="Mean CP-1.5logN")
plot(cp_var04, main="Var CP-1.5logN")
plot(cp_meanvar04, main="Mean + Var CP-1.5logN")

cp_method <- "SegNeigh"
cp_mean05 <- cpt.mean(ts_cp_noseas, penalty = 'None', method=cp_method, Q=5, test.stat='CUSUM')
cp_var05 <- cpt.var(ts_cp_noseas, penalty='None',method=cp_method, Q=5, test.stat='CSS')

plot(cp_mean05, main="Mean -CUSUM")
plot(cp_var05, main="Var -CSS")

cp_method <- "PELT"
cp_mean02 <- cpt.mean(ts_cp_noseas, method=cp_method)
cp_var02 <- cpt.var(ts_cp_noseas, method=cp_method)
cp_meanvar02 <- cpt.meanvar(ts_cp_noseas, method=cp_method)

plot(cp_mean02)
plot(cp_var02)
plot(cp_meanvar02)

cp_np_noseas <- cpt.np(ts_cp_noseas, penalty = "CROPS", pen.value = c(100,301), method = "PELT",
       test.stat = "empirical_distribution", class = FALSE, minseglen = 365)

cp_np_noseas$cpt.out
cp_np_noseas$changepoints

cp_np_res <- cpt.np(ts_cp_res, penalty = "CROPS", pen.value = c(100,301), method = "PELT",
                       test.stat = "empirical_distribution", class = FALSE, minseglen = 365)

cp_np_res$cpt.out
cp_np_res$changepoints

plot(data_cp_noseas, type="l")
abline(v=cp_np_res$changepoints[[1]],col="blue", lwd=4)

plot(data_cp_res, type="l")
abline(v=cp_np_res$changepoints[[1]],col="red", lwd=4)

cp_np@param.est[2]

plot(cp_np)

require(ecp)
data_ediv <- dataMB_day_DEA$res[!dataMB_day_DEA$NA_idx]
mat_ediv <- as.matrix(data_ediv)
cp_ediv <- e.divisive(X=mat_ediv, min.size=182, k=NULL, alpha=2, sig.lvl = 0.9)

ggplot(data=dataMB_day_DEA[!dataMB_day_DEA$NA_idx,]) +
  geom_line(aes(y=res,x=seq(1,length(temp),1))) +
  #geom_line(aes(y=raw,x=seq(1,length(temp),1)), colour="blue") +
  geom_vline(xintercept=cp_ediv$estimates, colour="red", linetype="dashed")

#********************************************************************************
# Apply WMA to temperature time series
data <- as.ts(dataMB_day_DEA$temp)

window_make <- function(len){
  window_point <- c(seq((len+1)/2,2,-1),1,seq(2,(len+1)/2,1))
  return(window_point)
}
# start with 7(1), 15(2) and 31(3) day windows (also 61(4) and 91(5))
window_size <- 91                         # total window length in days
window_point <- window_make(window_size)  # time points before and after t
#weights <- 2/((window_point)+1)          # exponential weights
weights <- 0 + (max(window_point)-(window_point-1))/
  max(window_point)                      # linear weights
#weights <- window_point/window_point     # equal weights

data_flt5_lin <- filter(data, weights/sum(weights), method="convolution", sides=2)

plot(data[1:1000], type="l")
lines(data_flt1[1:1000], col="red")
lines(data_flt2[1:1000], col="darkgreen")
lines(data_flt3[1:1000], col="blue")

plot(data[800:1000], type="l", lwd=3)
lines(data_flt2_lin[800:1000], col="red", lwd=3)         # linear weights
lines(data_flt2_exp[800:1000], col="darkgreen", lwd=3)   # exponential weights
lines(data_flt2_sim[800:1000], col="blue", lwd=3)        # equal weights

plot(data[800:1000], type="l", lwd=3)
lines(data_flt3_lin[800:1000], col="red", lwd=3)         # 31 days
lines(data_flt4_lin[800:1000], col="darkgreen", lwd=3)   # 61 days
lines(data_flt5_lin[800:1000], col="blue", lwd=3)        # 91 days

plot(data, type="l", lwd=3)
lines(data_flt3_lin, col="red", lwd=3)
lines(data_flt4_lin, col="green", lwd=3)   # 61 days
lines(data_flt5_lin, col="blue", lwd=3)  

plot(data, type="l", lwd=3)
lines(data_flt4_lin, col="green", lwd=3)
