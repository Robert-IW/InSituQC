library(tidyverse)
library(forecast)
require(fitdistrplus)
require(changepoint)
library(stlplus)
library(stringr)
library(zoo)


load("~/R/Data/SACTNdaily_v4.1.Rdata")

my.site <- c("Cape Agulhas","Betty's Bay","Mossel Bay","Mossel Bay","Stilbaai","Gordons Bay",
             "De Hoop","Muizenberg","Richards Bay","Roman Rock")
my.src <- c("SAWS","DAFF","DEA","SAWS","SAWS","SAWS","DAFF","SAWS","KZNSB","UWC")

i=10

# Data Preparation --------------------------------------------------------

data.in <- SACTN_daily_v4.2 %>% 
  dplyr::select(date,site,src,temp) %>% 
  filter(site==my.site[i],src==my.src[i])
rm(i)

doy <- format(data.in$date[1], "%j")
yr <- format(data.in$date[1], "%Y")
tot.yr <- length(unique(format(data.in$date, "%Y")))

data.ts <- ts(data.in$temp, frequency=365,
              start = c(as.numeric(yr),as.numeric(doy)))
na.ind <- which(is.na(data.ts))

data.nona <- na.interp(data.ts,lambda=NULL)

s.win.value <- length(unique(df.data$year))
stl.nona <- stlplus(data.nona, s.window=s.win.value,
                    s.degree = 1, t.window=1460, n.p = 365)

data.nona[na.ind] <- stl.nona$data$seasonal[na.ind]+stl.nona$data$trend[na.ind]

par(mfrow=c(2,1), mar=c(4,4,4,1))
plot(data.nona, type="l", col="red", lwd=2,
     ylab="Temperature", xlab="Years")
lines(data.ts, col="black", lwd=2)
legend("bottom", c("Data Interp", "Original"),
       lty=c(1,1), lwd=c(2,2), col=c("red", 'black'),
       inset=c(0,.4,0), xpd=T, horiz=T, bty='n')

df.data <- data.frame("temp"=data.nona,
                      "temp.org"=as.numeric(data.ts),
                      "date"=data.in$date,
                      "doy"=as.numeric(format(data.in$date, "%j")),
                      "dom"=format(data.in$date, "%d"),
                      "month"=as.numeric(format(data.in$date, "%m")),
                      "year"=format(data.in$date, "%Y"),
                      stringsAsFactors = F)

# find which months have more than 10% NA data
df.na.mon <- df.data %>%
  dplyr::select(year,month,temp.org) %>% 
  group_by(year,month) %>% 
  summarize(num.days=n(),num.NA=sum(is.na(temp.org))) %>% 
  mutate(perc.NA=as.integer(num.NA/num.days*100)) %>% 
  mutate(ignore=perc.NA>50)


# Daily Data Smoothing ---------------------------------------------------

#data.smooth <- ma(data.nona, order=31)   # pre smoothing the data makes no difference to trend
data.trend.ma <- ma(data.nona, order=365)
y <- data.nona
x <- seq(1,length(y),1)
data.trend.loess <- loess(y~x, span=.3, family="symmetric")
data.trend.loess <- ts(data.trend.loess$fitted, frequency=365,
                       start = c(as.numeric(yr),as.numeric(doy)))

data.trend <- data.trend.loess
df.data <- df.data %>% 
  mutate(trend=data.trend) %>% 
  mutate(notrend=temp-trend)

# Get seasonal by LOESS vs avg DOY ---------------------------------------------------

# get avg DOY after removing trend
df.data.spread <- df.data %>%
  dplyr::select(doy, year, notrend) %>% 
  spread(year, notrend)

df.data.spread$med <- apply(df.data.spread[,2:tot.yr+1],1,
                            function(x) median(x, na.rm=T))
df.data.spread$sd <- apply(df.data.spread[,2:tot.yr+1],1,
                            function(x) sd(x, na.rm=T))

# smooth the DOY seasonal data
y <- rep(df.data.spread$med,3)
x <- seq(1,length(y),1)
data.doy.loess <- loess(y~x, span=.15, family="symmetric")

plot(y, type="n", xlab="Days")
polygon(c(seq(1,1098,1),rev(seq(1,1098,1))),
        c(data.doy.loess$fitted+df.data.spread$sd,
          rev(data.doy.loess$fitted-df.data.spread$sd)),
        col="gray90", border=NA)
lines(y, lwd=2, col="gray30")
lines(data.doy.loess$fitted, col="red", lwd=2)
legend("bottomleft", c("DOY SD","DOY Median","LOESS Fit"),
       lty=c(1,1), lwd=c(2,2),col=c("gray90","gray30","red"),
       text.width=c(0,100,100),
       inset=c(0,.4,0), xpd=T, horiz=T, bty='n')

df.data.spread <- df.data.spread %>% 
  mutate(med.smoo=data.doy.loess$fitted[367:732])

# asign the seasonal data to day of the year
df.data$seas.avg <- NA
for (i in 1:366){
  idx <- which(df.data$doy==i)
  df.data$seas.avg[idx] <- df.data.spread$med.smoo[i]
}
df.data$seas.avg <- ts(df.data$seas.avg, frequency=365,
                       start = c(as.numeric(yr),as.numeric(doy)))
par(mfrow=c(1,1))
plot(df.data$notrend, type="l", col="gray", lwd=2,
     ylab="Anomnaly", xlab="Years")
lines(df.data$seas.avg, col="red", lwd=2)
legend("bottomleft", c("Original","DOY LOESS Fit"),
       lty=c(1,1), lwd=c(2,2),col=c("gray","red"),
       text.width=c(0,1),
       inset=c(0,1,0), xpd=T, horiz=T, bty='n')

# fit LOESS smoother to original temperature before STL decomp
y <- df.data$temp
x <- seq(1,length(y),1)
data.seas.loess <- loess(y~x, span=.01, family="symmetric")
data.seas.loess <- ts(data.seas.loess$fitted, frequency=365,
                      start = c(as.numeric(yr), as.numeric(doy)))

# fit 7 day MA smoother to original temperature before STL decomp
data.seas.ma <- ma(df.data$temp, order = 31)
data.seas.ma <- ts(data.seas.ma, frequency=365,
                      start = c(as.numeric(yr), as.numeric(doy)))

par(mfrow=c(2,1), mar=c(4,4,4,1))
plot(df.data$temp, type="l", col="gray", lwd=2,
     ylab="Temperature", xlab="Years")
lines(data.seas.loess, col="red", lwd=2)
lines(data.seas.ma, col="blue", lwd=3)
legend("bottomleft", c("Original","LOESS Smooth","MA Smooth"),
       lty=c(1,1), lwd=c(2,2),col=c("gray","red","blue"),
       text.width=c(0,1,1.5),
       inset=c(0,.4,0), xpd=T, horiz=T, bty='n')

s.win.value <- length(unique(df.data$year))
# fit.loess <- stlplus(data.seas.loess, s.window=s.win.value,
#                      s.degree = 1, t.window=1460, n.p = 365)
fit.loess <- stlplus(data.seas.ma, s.window=s.win.value,
                     s.degree = 1, t.window=1460, n.p = 365)

df.data$seas.stl <- ts(fit.loess$data$seasonal, frequency=365,
                       start = c(as.numeric(yr),as.numeric(doy)))

df.data$trend.stl <- ts(fit.loess$data$trend, frequency=365,
                        start = c(as.numeric(yr), as.numeric(doy)))

plot(data.nona, type="l", col="gray", lwd=2, 
     ylab="Temperature", xlab="Years")
lines(data.trend.ma, col="red", lwd=2)            # based on nona
lines(ma(data.trend.ma, order=365), col="chartreuse4", lwd=2)
lines(data.trend.loess, col="blue", lwd=2)        # based on nona
lines(df.data$trend.stl, col="darkorange", lwd=2) # based on MA smoothed
legend("bottomleft",
       c("Original","MA","Smooth MA", "LOESS","STL"),
       lty=c(1,1,1,1), lwd=c(2,2,2,2),
       col=c("gray","red","chartreuse4","blue","darkorange"),
       inset=c(0,.4,0), xpd=T, horiz=T, bty='n',
       text.width=c(0,.5,.4,.6,.6))

plot(df.data$notrend, type="l", col="blue", lwd=2,
     ylab="Anomaly", xlab="Years")                              # LOESS trend
lines(df.data$temp-df.data$trend.stl, col="darkorange", lwd=2)  # STL trend
legend("bottomleft", c("Detrended LOESS","Detrended STL"),
       lty=c(1,1), lwd=c(2,2),col=c("blue","darkorange"),
       text.width=c(0,1.5),
       inset=c(0,.4,0), xpd=T, horiz=T, bty='n')

plot(df.data$temp-df.data$trend.stl, type="l", col="gray", lwd=2,
     ylab="Anomaly", xlab="Years")
lines(df.data$seas.avg, col="blue", lwd=2)
lines(df.data$seas.stl, col="darkorange", lwd=2)
legend("bottomleft", c("STL Detrended","DOY Seaonal LOESS","MA Seasonal STL"),
       lty=c(1,1), lwd=c(2,2),col=c("gray","blue","darkorange"),
       text.width=c(0,1.5,1.8),
       inset=c(0,.4,0), xpd=T, horiz=T, bty='n')

# Get residuals and test AR -----------------------------------------------
# the trend is not removed instead the mean is subtracted before seasonal
# adjustment as the impacts on the change in mean
df.data$res.mean <- ts(df.data$temp - median(df.data$temp) - df.data$seas.stl,
                  frequency=365, start = c(as.numeric(yr),
                                           as.numeric(doy)))

# the residuals show no auto correlation or do they?
# is this important for change point detection?
styear <- 365-as.numeric(doy)
tickseq <- c(0+styear,365+styear,365*2+styear,365*3+styear,
             365*4+styear,365*5+styear)
layout(matrix(c(1,1,2,3),2,2,byrow = T))
plot(df.data$res.mean, type="l", main="Residuals (Deseasonalized by STL)",
     ylab="Anomaly", xlab="Years")
abline(h=0, col="red")
acf(as.numeric(df.data$res.mean), lag.max = 1825,
    main="ACF of Residuals", xlab="Years", xaxt="n")
axis(side=1, at=tickseq, labels=unique(format(data.in$date, "%Y"))[2:7])
pacf(as.numeric(df.data$res.mean), main="Partial ACF of Residuals",
     xlab="Days")

# if they do, try first difference to remove remaining seasonality

# but the square of the residuals show seasonal autocorrelation
#pacf(df.data$res.mean^2, lag.max = 1800)

# to detect change in variance we remove the trend
df.data$res.var <- ts(df.data$temp - df.data$trend - df.data$seas.stl,
                      frequency=365, start = c(as.numeric(yr),
                                               as.numeric(doy)))
# square the residuals to make them positive
df.data$res.var.sq <- df.data$res.var^2

# fit a loess smoother to doy squared residuals
df.res.spread <- df.data %>%
  dplyr::select(doy, year, res.var.sq) %>%
  spread(year, res.var.sq)
# df.res.spread <- df.data %>%
#   dplyr::select(doy, year, res.var) %>%
#   spread(year, res.var)

df.res.spread$res.var.sd <- apply(df.res.spread[,2:tot.yr+1],
                                  1, function(x) sd(x, na.rm=T))
y <- rep(df.res.spread$res.var.sd,3)
x <- seq(1,length(y),1)
data.sd.loess <- loess(y~x, span=.15, family="symmetric")

df.data.spread$res.var.sq.seas <- data.sd.loess$fitted[367:732]

df.data$res.var.sq.seas <- NA
for (i in 1:366){
  idx <- which(df.data$doy==i)
  df.data$res.var.sq.seas[idx] <- df.data.spread$res.var.sq.seas[i]
}
df.data$res.var.sq.seas <- ts(df.data$res.var.sq.seas, frequency=365, start = c(as.numeric(yr),
                                                                as.numeric(doy)))
# and divide residuals by the seasonal squared residuals to
df.data$unit.var <- df.data$res.var.sq / df.data$res.var.sq.seas
# df.data$unit.var <- df.data$res.var / df.data$res.var.sq.seas

par(mfrow=c(2,1), mar=c(4,4,4,1))
plot(df.res.spread$res.var.sd, type="l", lwd=2, col="gray",
     ylab="Variance",xlab="Days")
lines(data.sd.loess$fitted[367:732], col="red", lwd=2)
legend("bottomleft", c("DOY Variance","LOESS Smooth"),
       lty=c(1,1), lwd=c(2,2),col=c("gray","red"),
       text.width=c(0,30),
       inset=c(0,.4,0), xpd=T, horiz=T, bty='n')

plot(df.data$unit.var, type="l", lwd=2, col="gray30",
     main="Deseasonalized Variance (Var / Seas. Var)",
     ylab="Variance", xlab="Years")

pacf(df.data$unit.var, lag.max = 1800)

# plot(log(df.data$unit.var), type="l", lwd=2)
# hist(log(df.data$unit.var))
plot(df.data$unit.var, type="l", lwd=2)

par(mfrow=c(1,2))
hist(df.data$unit.var, breaks = seq(0,15,.5))
hist(log(df.data$unit.var), breaks = seq(-20,10,.5))

descdist(log(as.numeric(df.data$unit.var)))
descdist(as.numeric(df.data$unit.var))

# Daily Change Point Analysis ---------------------------------------------------------
dates <- df.data$date
my.method <- "SegNeigh"
my.penalty <- "Asymptotic"
my.pen.value <- 0.01
my.Q <- 5
my.minseglen <- 365

# the res.mean is normal so test.stat = "Normal"
#y.mean <- as.numeric(df.data$res.mean)
y.mean <- as.numeric(df.data$temp-df.data$seas.stl-df.data$trend.stl)

results.1 <- cpt.mean(y.mean, method=my.method,
                      penalty=my.penalty,
                      pen.value = my.pen.value,
                      Q=my.Q,
                      #minseglen = my.minseglen,
                      test.stat = "Normal")
cpts(results.1)
param.est(results.1)

results.2 <- cpt.var(y.mean, method=my.method,
                     penalty="Asymptotic",
                     pen.value=my.pen.value,
                     Q=my.Q,
                     #minseglen = my.minseglen,
                     test.stat = "Normal")
cpts(results.2)
param.est(results.2)

# y.var is heavy tailed so test.stat = "CSS"
# may be worth transforming by Box-Cox tranformation to Normal
my.method <- "PELT"
my.penalty <- "CROPS"
my.pen.value <- c(0,1)

y.var <- as.numeric(df.data$unit.var)
# y.var <- as.numeric(log(df.data$unit.var))

results.3 <- cpt.mean(y.var, method=my.method,
                      penalty = my.penalty,
                      pen.value = my.pen.value,
                      Q = my.Q,
                      minseglen = my.minseglen,
                      test.stat = "Normal")
cpts(results.3)
param.est(results.3)

results.4 <- cpt.var(y.var, method=my.method,
                     penalty=my.penalty,
                     pen.value=my.pen.value,
                     Q=my.Q,
                     minseglen = my.minseglen,
                     test.stat = "Exponential")
cpts(results.4)
param.est(results.4)

par(mfrow=c(2,2))

plot(results.1, cpt.col="red", xlab="Index", cpt.width=4)
plot(results.2, cpt.col="blue", xlab="Index", cpt.width=4)
plot(results.3, cpt.col="red", xlab="Index", cpt.width=4)
plot(results.4, cpt.col="blue", xlab="Index", cpt.width=4)

par(mfrow=c(2,1))
# plot the cp for the mean and variance based on non-detrended ts
plot(as.numeric(df.data$temp)~dates, type="l", col="black")
abline(v=dates[cpts(results.1)], col='red', lwd=3)
abline(v=dates[cpts(results.2)], col='blue', lwd=2)

# plot the cp for the variance based on var deseasonalized series
plot(as.numeric(df.data$temp)~dates, type="l", col="black")
abline(v=dates[cpts(results.3)], col='red', lwd=3)
abline(v=dates[cpts(results.4)], col='blue', lwd=2)

par()


# Monthly Average Data ----------------------------------------------------
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

par(mfrow=c(1,1))
plot(data.nona, type="l", col="gray", lwd=3, ylab="Temperature", xlab="Years")
lines(df.data.yrmon$mean, col="red", lwd=3)
lines(df.data.yrmon$med, col="chartreuse4", lwd=3)
abline(v=st.year:end.year, col="gray30", lty=2)
legend("bottomleft", c("Mean Monthly", "Median Monthly"),
       lty=c(1,1), lwd=c(2,2), col=c("red", 'chartreuse4'),
       inset=c(0,1,0), xpd=T, horiz=T, bty='n')

par(mfrow=c(2,1), mar=c(4,4,4,1))
plot.ts(df.data.yrmon$med, lwd=2, ylab="Temperature", xlab="Years")
abline(v=st.year:end.year, col="gray30", lty=2)
legend("bottomleft", c("Median Monthly"),
       lty=c(1,1), lwd=c(2,2), col="black",
       inset=c(0,.6,0), xpd=T, horiz=T, bty='n')
plot.ts(df.data.yrmon$sd, lwd=2, ylab="Temperature", xlab="Years")
abline(v=st.year:end.year, col="gray30", lty=2)
legend("bottomleft", c("Std. Dev. Monthly"),
       lty=c(1,1), lwd=c(2,2), col="black",
       inset=c(0,.6,0), xpd=T, horiz=T, bty='n')


# Daily Data by Month -----------------------------------------------------

# create monthly summaries
df.data.mon <- df.data %>%
  dplyr::select(dom, month, year, temp) %>% 
  mutate(year.dom=paste0(year,".",dom)) %>%
  mutate(temp=as.numeric(temp)) %>% 
  dplyr::select(month, year.dom, temp) %>% 
  spread(year.dom, temp) %>% 
  dplyr::select(-month)

# from the daily time series by month, compute CPT for mean and variance of standardized data

# create a series of normalized means and variances based on cpts

# get the average mean and variance across all months to indicate more consistent change

df.data.mon$mean <- apply(df.data.mon,1,mean,na.rm=T)
df.data.mon$sd <- apply(df.data.mon[,!names(df.data.mon) %in% "mean"],1,sd,na.rm=T)

s.win.value <- length(unique(df.data$year))

# STL decomp of median temperature
fit.loess <- stlplus(df.data.yrmon$med, s.window=s.win.value, t.window = 60,
                     s.degree = 1, n.p = 12)
df.data.yrmon$seas.stl <- ts(fit.loess$data$seasonal, frequency=12,
                       start = c(st.year,st.mon))

df.data.yrmon$trend.stl <- ts(fit.loess$data$trend, frequency=12,
                        start = c(st.year, st.mon))
rm(fit.loess)

# STL decomp of SD temperature
fit.loess <- stlplus(df.data.yrmon$sd, s.window="periodic", t.window = 60,
                     s.degree = 1, n.p = 12)
df.data.yrmon$seas.sd <- ts(fit.loess$data$seasonal, frequency=12,
                             start = c(st.year,st.mon))

df.data.yrmon$trend.sd <- ts(fit.loess$data$trend, frequency=12,
                              start = c(st.year, st.mon))


plot(data.nona, type="l", col="gray", lwd=2, 
     ylab="Temperature", xlab="Years")
lines(df.data.yrmon$trend.stl, col="darkorange", lwd=2) # based on MA smoothed
legend("bottomleft",
       c("Original","STL Trend"),
       lty=c(1,1,1,1), lwd=c(2,2,2,2),
       col=c("gray","darkorange"),
       inset=c(0,.4,0), xpd=T, horiz=T, bty='n',
       text.width=c(0,.5,.4,.6,.6))

plot(data.nona, type="l", col="gray", lwd=2, 
     ylab="Temperature", xlab="Years")
lines(df.data.yrmon$trend.stl+df.data.yrmon$seas.stl, col="darkorange", lwd=2) # based on MA smoothed
legend("bottomleft",
       c("Original","STL Season"),
       lty=c(1,1,1,1), lwd=c(2,2,2,2),
       col=c("gray","darkorange"),
       inset=c(0,.4,0), xpd=T, horiz=T, bty='n',
       text.width=c(0,.5,.4,.6,.6))

par(mfrow=c(1,1))
plot(df.data.yrmon$sd, type="l", main="Monthly SD with Trend",
     ylab="Std. Dev", xlab="Years")
lines(df.data.yrmon$seas.sd+df.data.yrmon$trend.sd,
      col="darkorange", lwd=2)

# Monthly Change Point Analysis -------------------------------------------

# get residuals after removing seasonal and not trend
par(mfrow=c(2,1), mar=c(4,4,4,1))

y.mean <- df.data.yrmon$med-mean(df.data.yrmon$trend.stl)-df.data.yrmon$seas.stl
plot(y.mean, type="l", ylab="Anomaly", main="Deseasonalized Monthly Temperature")
abline(v=st.year:end.year, col="gray30", lty=2)
abline(h=0, col="gray30", lty=2)

y.sd <- df.data.yrmon$sd / (df.data.yrmon$seas.sd+mean(df.data.yrmon$trend.sd))
plot(y.sd, type="l", ylab="Ratio SD:E(SD)", main="Monthly SD to E(SD) Ratio")
abline(v=st.year:end.year, col="gray30", lty=2)
abline(h=1, col="gray30", lty=2)

# Remove months with > 50% NA
y.mean[df.na.mon$ignore] <- NA

# Change point analysis on SD ratios
dates.mon <- paste(df.data.yrmon$year,
               sprintf("%02d",df.data.yrmon$month), sep="-")
dates.mon <- as.yearmon(dates.mon,"%Y-%m")

dates.day <- as.yearmon(df.data$date, "%Y-%m-%d")

my.method <- "PELT"
my.penalty <- "Hannan-Quinn"
my.pen.value <- 0.01
my.Q <- 10

results.1 <- cpt.mean(y.mean, method=my.method,
                      penalty=my.penalty,
                      pen.value = my.pen.value,
                      Q=my.Q,
                      minseglen = my.minseglen,
                      test.stat = "CUSUM")
cpts(results.1)
param.est(results.1)

results.2 <- cpt.var(y.mean, method=my.method,
                     penalty=my.penalty,
                     pen.value=my.pen.value,
                     Q=my.Q,
                     minseglen = my.minseglen,
                     test.stat = "CSS")
cpts(results.2)
param.est(results.2)

results.3 <- cpt.mean(y.sd, method=my.method,
                      penalty=my.penalty,
                      pen.value = my.pen.value,
                      Q=my.Q,
                      minseglen = my.minseglen,
                      test.stat = "Normal")
cpts(results.3)
param.est(results.3)

results.4 <- cpt.var(y.sd, method=my.method,
                     penalty=my.penalty,
                     pen.value=my.pen.value,
                     Q=my.Q,
                     minseglen = my.minseglen,
                     test.stat = "Normal")
cpts(results.4)
param.est(results.4)

par(mfrow=c(2,2))

plot(results.1, cpt.col="red", xlab="Years", cpt.width=4,
     main="Mean CPT for Deseasonalized Temperature")
abline(v=st.year:end.year, col="gray30", lty=2)
abline(h=0, col="gray30", lty=2)
plot(results.2, cpt.col="blue", xlab="Years", cpt.width=4,
     main="Variance CPT for Deseasonalized Temperature")
abline(v=st.year:end.year, col="gray30", lty=2)
abline(h=0, col="gray30", lty=2)
plot(results.3, cpt.col="red", xlab="Years", cpt.width=4,
     main="Mean CPT for SD:E(SD)")
abline(v=st.year:end.year, col="gray30", lty=2)
abline(h=1, col="gray30", lty=2)
plot(results.4, cpt.col="blue", xlab="Years", cpt.width=4,
     main="Variance CPT for SD:E(SD)")
abline(v=st.year:end.year, col="gray30", lty=2)
abline(h=1, col="gray30", lty=2)

par(mfrow=c(2,1))
# plot the cp for the mean and variance based on non-detrended ts
plot(df.data$temp, type="l", col="black")
abline(v=dates.mon[cpts(results.1)], col='red', lwd=3)
abline(v=dates.mon[cpts(results.2)], col='blue', lwd=2)
abline(v=st.year:end.year, col="gray30", lty=2)

# plot the cp for the variance based on var deseasonalized series
plot(df.data$temp, type="l", col="black")
abline(v=dates.mon[cpts(results.3)], col='red', lwd=3)
abline(v=dates.mon[cpts(results.4)], col='blue', lwd=2)
abline(v=st.year:end.year, col="gray30", lty=2)

# CPM Change Point --------------------------------------------------------

library(cpm)

# ECP Change Point --------------------------------------------------------

# non parametric change point detection, assumes observations are iid
require(ecp)
mat_ediv <- as.matrix(y.mean)

cp_ediv <- e.divisive(X=mat_ediv,
                      min.size=365,
                      k=NULL,
                      alpha=1,
                      sig.lvl = 0.05,
                      R=199)

ggplot() +
  geom_line(aes(y=mat_ediv[,1],x=seq(1,nrow(mat_ediv),1))) +
  #geom_line(aes(y=raw,x=seq(1,length(temp),1)), colour="blue") +
  geom_vline(xintercept=cp_ediv$estimates, colour="red", linetype="dashed")



df.data$week.mean <- NA
df.data$week.sd <- NA
df.data$seas.week.mean <- NA
df.data$seas.week.sd <- NA

i <- 4

while (i < nrow(df.data)-3){
  df.data$week.mean[i] <- mean(df.data$temp[(i-3):(i+3)], na.rm=T)
  df.data$week.sd[i] <- sd(df.data$temp[(i-3):(i+3)], na.rm=T)
  
  df.data$seas.week.mean[i] <- mean(df.data$seas.me[(i-3):(i+3)], na.rm=T)
  df.data$seas.week.sd[i] <- mean(df.data$seas.sd[(i-3):(i+3)], na.rm=T)
  i <- i+1
}

df.data$cv <- df.data$week.sd / df.data$week.mean

plot(df.data$week.mean, type="l")
lines(df.data$seas.week.mean, col="red")
plot(df.data$week.sd, type="l")
lines(df.data$seas.week.sd, col="red")
