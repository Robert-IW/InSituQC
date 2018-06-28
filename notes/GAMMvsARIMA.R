require(tidyr)
require(ggplot2)
require(mgcv)
require(forecast)

# GAMM, a mixed model with the GAM for the seasonal and trend, and the AR
# process for the residuals

# fill in the NA data from the 'forecast' package for exponential smoothing
doy <- format(dataMB_day_DEA$date[1], "%j")
yr <- format(dataMB_day_DEA$date[1], "%Y")
ts_DEA <- ts(dataMB_day_DEA$temp, frequency=365, start = c(as.numeric(yr),
                                                           as.numeric(doy)))
fit <- na.interp(ts_DEA)

# Fitting a GAMM to seasonal data
gam_DEA <- dataMB_day_DEA[c("date", "temp", "month", "year")]
gam_DEA$fit <- as.vector(fit)
gam_DEA$doy <- as.numeric(format(gam_DEA$date, "%j"))
gam_DEA$datenum <- seq_along(gam_DEA$date)

# This first model does not take into account the autocorrelation
m <- gamm(temp ~ s(doy, bs="cc", k=366) + s(datenum, bs="cr"),
          data=gam_DEA, na.action=na.omit)

summary(m$gam)

layout(matrix(1:2, ncol=2))
plot(m$gam, scale=0)

acf(resid(m$lme), lag.max = 30, main="ACF")
pacf(resid(m$lme), lag.max = 30, main="pACF")
layout(1)
#-----------------------------------------------------------------
# fit a AR(p={1,2,3})
ctrl <- list(niterEM=0, msVerbose=T, optimMethod="L-BFGS-B")

p1 <-
  gamm(
    fit ~ s(doy, bs = "cc", k = 365) + s(datenum, bs = "cr"),   # k is the no. of unique values
    correlation = corARMA(form = ~ 1 | doy, p=1),control = ctrl,
    data = gam_DEA
  )
p2 <-
  gamm(
    fit ~ s(doy, bs = "cc", k = 365) + s(datenum, bs = "cr"),   # k is the no. of unique values
    correlation = corARMA(form = ~ 1 | doy, p=2),control = ctrl,
    data = gam_DEA, na.action = na.omit
  ) # this one fails
p3 <-
  gamm(
    fit ~ s(doy, bs = "cc", k = 365) + s(datenum, bs = "cr"),   # k is the no. of unique values
    correlation = corARMA(form = ~ 1 | doy, p=3), control = ctrl,
    data = gam_DEA
  )

p3_reml <-
  gamm(
    fit ~ s(doy, bs = "cc", k = 365) + s(datenum, bs = "cr"),   # k is the no. of unique values
    correlation = corARMA(form = ~ 1 | year, p=3), method = "REML",
    data = gam_DEA
  )


anova(m$lme,p1$lme,p3$lme)

plot(p1$gam, scale=0)

res <- resid(p3$lme, type = "normalized")
acf(res, lag.max=30, main = "ACF AR(3) residuals")
pacf(res, lag.max=30, main = "pACF AR(3) residuals")
layout(1)

# Fitting a seasonal  ARIMA
a <- p3$gam$y
b <- p3$lme$fitted
df <- data.frame(a=a, b=b, x=seq_along(a))

ggplot(df) + 
  geom_line(aes(x=x, y=a), col="azure4") + 
  geom_line(aes(x=x, y=b.fixed), col="blue", lwd=.8) +
  geom_line(aes(x=x, y=b.g), col="red", lwd=1) +
  geom_line(aes(x=x, y=b.g.0), col="chartreuse3", lwd=1)

# are the resulting residuals normally distributed?
res_gam <- as.numeric(p3$gam$residuals)

#******************************************************************************
lowpass_spline <- smooth.spline(fit, spar=0.05)

plot(fit, col="gray")
lines(lowpass_spline, col="blue", lwd=2)

# smooth the data for the GAMM using 61 day window and linear decay weights
window_size <- 61                         # total window length in days
window_point <- window_make(window_size)  # time points before and after t
#weights <- 2/((window_point)+1)          # exponential weights
weights <- 0 + (max(window_point)-(window_point-1))/
  max(window_point)                      # linear weights
#weights <- window_point/window_point     # equal weights

data_flt4_lin <- filter(ts_DEA, weights/sum(weights), method="convolution", sides=2)
lines(data_flt4_lin, col="red", lwd=2)

# NOTE: The low pass spline is almost identical to the 61 WMA
# apply the GAMM to the low pass spline
fit_spline <- lowpass_spline$y
gam_DEA$fit_sp <- as.vector(fit_spline)
# fit a AR(p={1,2,3})
ctrl <- list(niterEM=0, msVerbose=T, optimMethod="L-BFGS-B")

p0.sp <-
  gamm(
    fit_sp ~ s(doy, bs = "cc", k = 365) + s(datenum, bs = "cr"),
    method = "REML", data = gam_DEA
  )

p1.sp <-
  gamm(
    fit_sp ~ s(doy, bs = "cc", k = 365) + s(datenum, bs = "cr"),   # k is the no. of unique values
    correlation = corARMA(form=~1|year,p=1), method = "REML",
    data = gam_DEA
  )
p2.sp <-
  gamm(
    fit_sp ~ s(doy, bs = "cc", k = 365) + s(datenum, bs = "cr"),   # k is the no. of unique values
    correlation = corARMA(form=~1|year,p=2), method = "REML",
    data = gam_DEA, na.action = na.omit
  )
p3.sp <-
  gamm(
    fit_sp ~ s(doy, bs = "cc", k = 365) + s(datenum, bs = "cr"),   # k is the no. of unique values
    correlation = corARMA(form=~1|year,p=3), method = "REML",
    data = gam_DEA
  )

anova(p1$lme,p2$lme,p3$lme)
anova(p1.sp$lme,p2.sp$lme,p3.sp$lme)

intervals(p1$lme, which = "var-cov")$corStruct

layout(matrix(1:4, ncol = 4))
pacf(resid(p0.sp$lme, type = "normalized"), lag.max = 48, main = "pACF of AR(0)")
pacf(resid(p1.sp$lme, type = "normalized"), lag.max = 48, main = "pACF of AR(1)")
pacf(resid(p2.sp$lme, type = "normalized"), lag.max = 48, main = "pACF of AR(2)")
pacf(resid(p3.sp$lme, type = "normalized"), lag.max = 48, main = "pACF of AR(3)")

arma_res <- auto.arima(resid(p0.sp$lme, type = "normalized"),
                       stationary = TRUE, seasonal = FALSE)

arma_res$coef

# Fitting a seasonal  ARIMA
a <- p1.sp$gam$y
b <- p0.sp$lme$fitted
c <- p1.sp$lme$fitted
d <- p2.sp$lme$fitted
e <- p3.sp$lme$fitted
f <- p3$lme$fitted

df <- data.frame(a=a, b=b, c=c, d=d, e=e, f=f, x=seq_along(a))

ggplot(df) + 
  geom_line(aes(x=x, y=a), col="azure4") + 
  geom_line(aes(x=x, y=b.g.0), col="chartreuse3", lwd=1) +
  geom_line(aes(x=x, y=c.g.0), col="blue", lwd=1) +
  geom_line(aes(x=x, y=d.g.0), col="red", lwd=1) +
  geom_line(aes(x=x, y=e.g.0), col="gold1", lwd=1)

ggplot(df) + 
  geom_line(aes(x=x, y=a), col="azure4") + 
  geom_line(aes(x=x, y=b.fixed), col="blue", lwd=.8) +
  geom_line(aes(x=x, y=b.g), col="red", lwd=1) +
  geom_line(aes(x=x, y=b.g.0), col="chartreuse3", lwd=1) 