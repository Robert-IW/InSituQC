library(mgcv)

test.data <- df.data$temp
#test.data <- df.data$res.mean
#plot(df.data$res.mean, type="l")

plot(test.data, type="l")
# the data are non-stationary with a clear change in mean

# acf includes the intervening t between Tt and Tt-x
# shows both serial and seasonal autocorrelation
acf(as.numeric(test.data), lag.max = 365*10)

# pacf looks at direct covariance with Tt and Tt-x
# below shows significant correlations up to 15 days
pacf(as.numeric(test.data), lag.max = 30) 

# take the first difference to make the data stationary (no trend)
# however this will remove the change in mean
data.diff <- diff(df.data$res.mean, lag=1)
plot(data.diff, type="l")

acf(data.diff)
pacf(data.diff, lag.max = 1800)
pacf(data.diff^2, lag.max = 1800)


# Test GAMM model ---------------------------------------------------------
ctrl <- list(niterEM=0, msVerbose=T, optimMethod="L-BFGS-B")
# Fitting a GAMM to seasonal data
gam_res <- data.in[c("date", "temp")]
gam_res$year <- as.numeric(format(gam_res$date, "%Y"))
gam_res$fit <- as.vector(data.nona)
gam_res$doy <- as.numeric(format(gam_res$date, "%j"))
gam_res$datenum <- seq_along(gam_res$date)

p3_reml <-
  gamm(
    fit ~ s(doy, bs = "cc", k = 365) + s(datenum, bs = "cr"),   # k is the no. of unique values
    correlation = corARMA(form = ~ 1, p=5, q=5), method = "REML",
    data = gam_res
  )

a <- p3_reml$gam$y
b <- p3_reml$lme$fitted
df <- data.frame(a=a, b=b, x=seq_along(a))

ggplot(df) + 
  geom_line(aes(x=x, y=a), col="azure4") + 
  geom_line(aes(x=x, y=b.fixed), col="blue", lwd=.8) +
  geom_line(aes(x=x, y=b.g), col="red", lwd=1) +
  geom_line(aes(x=x, y=b.g.0), col="chartreuse3", lwd=1)

res_gam <- as.numeric(p3_reml$gam$residuals)
plot(res_gam, type="l")

acf(res_gam)
pacf(res_gam)
