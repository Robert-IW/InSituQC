library(data.table)
library(ggplot2)
library(gridExtra)
library(plyr)
library(dplyr)
library(qcc)

theme_update(plot.title = element_text(hjust = 0.5))

# function to convert month to season
mon2seas <- function(month){
  if (month %in% c(12,1,2)){
    return("summer")
  } else if (month %in% c(3,4,5)){
    return("autumn")
  } else if (month %in% c(6,7,8)){
    return("winter")
  } else if (month %in% c(9,10,11)){
    return("spring")
  }
}

# load Mossel Bay SST data -----------------------------------------------------------------

dataMB_avhrr1 <- read.csv(file="/media/robert/KELP-HDD-Portable/SST/GHRSST/L4/AVHRR-OI/ts/AVHRR-OI_L4_Mossel Bay_station_SST_timeseries_5nearest.csv",
                          header = T, stringsAsFactors = F)

dataMB_avhrr2 <- read.csv(file="/media/robert/KELP-HDD-Portable/SST/GHRSST/L4/AVHRR-OI/ts/ESdata/AVHRR-OI_Mossel Bay_SST_timeseries_5nearest.csv",
                          header = T, stringsAsFactors = F)

#dataMB_avhrr3 <- read.csv(file="/media/robert/KelpPhysics1/Backup Satellite Data/SST/GHRSST/L4/AVHRR-OI/ts/AVHRR-OI_L4_Mossel Bay_station_SST_timeseries_5nearest.csv",
#                          header = T, stringsAsFactors = F)

dataMB_SST <- dataMB_avhrr2[c("date","nearest1")]
dataMB_SST$datetime <- as.POSIXct(as.character(dataMB_SST$date), format='%Y%m%d')
dataMB_SST <- dataMB_SST %>% mutate(month = as.numeric(format(datetime, "%m")))

dataMB_SST <- dataMB_SST %>% mutate(seas = sapply(month,mon2seas))
dataMB_SST$seas <- factor(dataMB_SST$seas, levels=c("summer","autumn","winter","spring"))
dataMB_SST$year <- format(dataMB_SST$datetime, format="%Y")

# plot monthly histograms
ggplot(aes(x = nearest1, fill = factor(month)), data = dataMB_SST) +
  geom_histogram(binwidth=.2) +
  xlim(12,25) +
  facet_grid(month ~ ., margins = FALSE, scales = "free") +
  ggtitle("Mossel Bay AVHRR SST.") +
  scale_fill_discrete(guide=FALSE)

# get monthly mean and sd
stat_SST <- dataMB_SST %>% group_by(month,year) %>% summarise(mean = mean(nearest1, na.rm=T),
                                                                 sd = sd(nearest1, na.rm=T))
stat_SST$year <- as.Date(stat_SST$year, format="%Y")

ggplot(data = stat_SST, aes(y=mean, x=year, group=month)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=.2) +
  facet_grid(month ~ ., margins = FALSE, scales = "free") +
  ggtitle("Mossel Bay AVHRR SST.") +
  scale_fill_discrete(guide=FALSE) +
  scale_x_date(date_labels = "%Y", date_breaks = "2 year") +
  theme(axis.text.x = element_text(angle = 45))

########################################################################################3
# load Mossel Bay Wind Vectors - combine each 'local' year file ----------------------------
files <- list.files("/media/robert/KELP-HDD-Portable/SurfaceWind/CCMP/ts/",
           pattern = glob2rx("Mossel*local*"),full.names = T)
dataMB_wind <- adply(files, 1, read.csv, .id = NULL)
# combine date and time
dataMB_wind$datetime <- as.POSIXct(paste(dataMB_wind$date, sprintf("%04.0f", dataMB_wind$time)), format='%Y%m%d %H%M')
# relable time and add month
dataMB_wind <- dataMB_wind %>% mutate(time = factor(time, levels = c(200, 800, 1400, 2000),
                                                 labels = c("02:00","08:00","14:00","20:00")),
                                      month = as.numeric(format(as.Date(datetime), "%m")))

detach(package:plyr)
stat_wndSpd <- dataMB_wind %>% group_by(time) %>% summarise(mean = mean(wndSpd, na.rm=T),
                                                   sd = sd(wndSpd, na.rm=T))

stat_Uwnd <- dataMB_wind %>% group_by(time) %>% summarise(mean = mean(Uwnd, na.rm=T),
                                                            sd = sd(Uwnd, na.rm=T))

stat_Vwnd <- dataMB_wind %>% group_by(time) %>% summarise(mean = mean(Vwnd, na.rm=T),
                                                          sd = sd(Vwnd, na.rm=T))

dataMB_wind <- dataMB_wind %>% mutate(seas = sapply(month,mon2seas))
dataMB_wind$seas <- factor(dataMB_wind$seas, levels=c("summer","autumn","winter","spring"))
dataMB_wind$year <- format(dataMB_wind$datetime, format="%Y")

save(dataMB_wind, file="~/R/Docs/dataMB_wind.Rdata")

ggplot(aes(y = wndSpd, x = seas, fill = time), data = dataMB_wind) +  geom_boxplot() + 
  ggtitle("Mossel Bay Wind Speed")
ggplot(aes(y = Uwnd, x = seas, fill = time), data = dataMB_wind) + geom_boxplot() +
  ggtitle("Mossel Bay U Speed")
ggplot(aes(y = Vwnd, x = seas, fill = time), data = dataMB_wind) + geom_boxplot() +
  ggtitle("Mossel Bay V Speed")

ggplot(aes(x = wndSpd, fill = seas), data = dataMB_wind) +  geom_histogram(binwidth = .5) +
  facet_grid(seas ~ ., margins = FALSE, scales = "free") +
  ggtitle("Mossel Bay Wind Speed")
ggplot(aes(x = Uwnd, fill = seas), data = dataMB_wind) +  geom_histogram(binwidth=.5) + 
  facet_grid(seas ~ ., margins = FALSE, scales = "free") +
  ggtitle("Mossel Bay Uwnd") +
  geom_vline(aes(xintercept = 0), color="red", linetype="dashed", size=1)
ggplot(aes(x = Vwnd, fill = seas), data = dataMB_wind) +  geom_histogram(binwidth=.5) + 
  facet_grid(seas ~ ., margins = FALSE, scales = "free") +
  ggtitle("Mossel Bay Vwnd") + geom_vline(aes(xintercept = 0), color="red", linetype="dashed", size=1)

# summary seasonal statistics for each 6 hour period                                  
stat_wndSpd_seas <- dataMB_wind %>% group_by(time,seas) %>% summarise(mean = mean(wndSpd, na.rm=T),
                                                            sd = sd(wndSpd, na.rm=T))

stat_Uwnd_seas <- dataMB_wind %>% group_by(time,seas) %>% summarise(mean = mean(Uwnd, na.rm=T),
                                                          sd = sd(Uwnd, na.rm=T))

stat_Vwnd_seas <- dataMB_wind %>% group_by(time,seas) %>% summarise(mean = mean(Vwnd, na.rm=T),
                                                          sd = sd(Vwnd, na.rm=T))

# plot monthly mean and sd with time
stat_Uwind <- dataMB_wind %>% group_by(month,year) %>% summarise(mean = mean(Uwnd, na.rm=T),
                                                                 median = median(Uwnd, na.rm=T),
                                                                  sd = sd(Uwnd, na.rm=T))
stat_Uwind$year <- as.Date(stat_Uwind$year, format="%Y")

stat_Vwind <- dataMB_wind %>% group_by(month,year) %>% summarise(mean = mean(Vwnd, na.rm=T),
                                                                 median = median(Vwnd, na.rm=T),
                                                                 sd = sd(Vwnd, na.rm=T))
stat_Vwind$year <- as.Date(stat_Vwind$year, format="%Y")

stat_Wind <- dataMB_wind %>% group_by(month,year) %>% summarise(mean = mean(wndSpd, na.rm=T),
                                                                 median = median(wndSpd, na.rm=T),
                                                                 sd = sd(wndSpd, na.rm=T))
stat_Wind$year <- as.Date(stat_Wind$year, format="%Y")

save(stat_Uwind, file="~/R/Docs/stat_Uwind.Rdata")
save(stat_Vwind, file="~/R/Docs/stat_Vwind.Rdata")
save(stat_Wind, file="~/R/Docs/stat_Wind.Rdata")

ggplot(data = stat_Uwind, aes(y=mean, x=year, group=month)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=.2) +
  geom_hline(aes(yintercept = 0), color="red", size=.5) +
  facet_grid(month ~ ., margins = FALSE, scales = "free") +
  ggtitle("Mossel Bay Uwind.") +
  scale_fill_discrete(guide=FALSE) +
  scale_x_date(date_labels = "%Y", date_breaks = "2 year") +
  theme(axis.text.x = element_text(angle = 45))

ggplot(data = stat_Vwind, aes(y=mean, x=year, group=month)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=.2) +
  geom_hline(aes(yintercept = 0), color="red", size=.5) +
  facet_grid(month ~ ., margins = FALSE, scales = "free") +
  ggtitle("Mossel Bay Vwind.") +
  scale_fill_discrete(guide=FALSE) +
  scale_x_date(date_labels = "%Y", date_breaks = "2 year") +
  theme(axis.text.x = element_text(angle = 45))

#######################################################################################
library()
# create daily mean of sub-daily wind
df1 <- dataMB_wind %>% 
  group_by(date) %>% 
  summarise(Uwind = median(Uwnd))
df2 <- dataMB_wind %>% 
  group_by(date) %>% 
  summarise(Vwind = median(Vwnd))
# create daily mean of sub-daily wind
df3 <- dataMB_wind %>% 
  group_by(date) %>% 
  summarise(Wspeed = median(wndSpd))

# create daily cumulative sum of absolut sub-daily wind
df3 <- dataMB_wind %>% 
  group_by(date) %>% 
  summarise(Uabs = mean(abs(Uwnd)))
df4 <- dataMB_wind %>% 
  group_by(date) %>% 
  summarise(Vabs = mean(abs(Vwnd)))

# merge the data sets on date and convert date
dataMB_cumWind <- Reduce(function(x, y) merge(x, y, all=TRUE), list(df1, df2, df3, df4))
dataMB_cumWind$datetime <- as.POSIXct(as.character(dataMB_cumWind$date), format='%Y%m%d')
dataMB_cumWind <- dataMB_cumWind[,-1]

dataMB_rawWind <- Reduce(function(x, y) merge(x, y, all=TRUE), list(df1, df2, df3))
dataMB_rawWind$datetime <- as.POSIXct(as.character(dataMB_rawWind$date), format='%Y%m%d')
dataMB_match <- merge(dataMB_rawWind, dataMB_SST)
dataMB_DEA_match <- merge(dataMB_rawWind, dataMB_day_DEA, by="date")


p1 <- ggplot(dataMB_match, aes(x=Uwnd, y=nearest1)) +
  geom_point(shape=1)
p2 <- ggplot(dataMB_match, aes(x=Vwnd, y=nearest1)) +
  geom_point(shape=1)
p3 <- ggplot(dataMB_match, aes(x=Uabs, y=nearest1)) +
  geom_point(shape=1)
p4 <- ggplot(dataMB_match, aes(x=Vabs, y=nearest1)) +
  geom_point(shape=1)

grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)

p1 <- ggplot(dataMB_match[1:1000,], (aes(x=datetime, y=Uwnd))) +
  geom_line() +
  geom_hline(yintercept = 0, colour="red")
p2 <- ggplot(dataMB_match[1:1000,], (aes(x=datetime, y=Vwnd))) +
  geom_line() +
  geom_hline(yintercept = 0, colour="red")
p3 <- ggplot(dataMB_match[1:1000,], (aes(x=datetime, y=Uabs))) +
  geom_line()
p4 <- ggplot(dataMB_match[1:1000,], (aes(x=datetime, y=Vabs))) +
  geom_line()

grid.arrange(p1,p3,nrow=2)
grid.arrange(p2,p4,nrow=2)

stat_SSTmean <- dataMB_SST %>% group_by(month) %>% summarise(avgSST = mean(nearest1, na.rm=T))

stat_match <- dataMB_match %>% group_by(month,year) %>% 
  summarise(mean = mean(nearest1, na.rm=T),
            median = median(nearest1, na.rm=T),
            sd = sd(nearest1, na.rm=T))

stat_match$year <- as.Date(stat_match$year, format="%Y")
save(stat_match, file="~/R/Docs/stat_match.Rdata")

ggplot(data = stat_match, aes(y=meanSST, x=year, group=month)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=meanSST-sdSST, ymax=meanSST+sdSST),width=.2) +
  facet_grid(month ~ ., margins = FALSE, scales = "free") +
  ggtitle("Mossel Bay SST.") +
  scale_fill_discrete(guide=FALSE) +
  scale_x_date(date_labels = "%Y", date_breaks = "2 year") +
  theme(axis.text.x = element_text(angle = 45))

ggplot(dataMB_match) +
  geom_line(aes(x=datetime,y=Ucum))

###############################################################################################################
# Detect anomaly by CUSUM and/or EWMA
require(qcc)
dataSST <- data.frame("date"=dataMB_SST$datetime,"sst"=dataMB_SST$nearest1,stringsAsFactors = F)
dataWind <- data.frame("date"=dataMB_wind$datetime, "uwind"=dataMB_wind$Uwnd,"vwind"=dataMB_wind$Vwnd,
                       "dir"=dataMB_wind$wndDir,"spd"=dataMB_wind$wndSpd,stringsAsFactors = F)

# create ts object
doy <- format(dataSST$date[1], "%j")
yr <- format(dataSST$date[1], "%Y")
ts_sst <- ts(dataSST$sst, frequency=365, start = c(as.numeric(yr),as.numeric(doy)))

doy <- format(dataWind$date[1], "%j")
yr <- format(dataWind$date[1], "%Y")
ts_uwind <- ts(dataWind$uwind, frequency=365*4, start = c(as.numeric(yr),as.numeric(doy)))
ts_vwind <- ts(dataWind$vwind, frequency=365*4, start = c(as.numeric(yr),as.numeric(doy)))
ts_spd <- ts(dataWind$spd, frequency=365*4, start = c(as.numeric(yr),as.numeric(doy)))
#ts_dir <- ts(dataWind$dir, frequency=365*4, start = c(as.numeric(yr),as.numeric(doy)))

# check for normality
require(nortest)          
ad.test(dataSST$sst)      # non_normal
ad.test(dataWind$uwind)   # non-normal
ad.test(dataWind$vwind)   # non-normal
ad.test(dataWind$spd)     # non-normal

# fit distribution
require(fitdistrplus)
descdist(dataWind$uwind)

fit.norm <- fitdist(dataSST$sst, "norm")        # not normality, bimodal
plot(fit.norm)
fit.norm <- fitdist(dataWind$uwind, "norm")     # assume normality
plot(fit.norm)
fit.norm <- fitdist(dataWind$vwind, "norm")     # assume normality
plot(fit.norm)
fit.norm <- fitdist(dataWind$spd, "norm")       # not normal
plot(fit.norm)

descdist(dataWind$spd)                          # possible gamma
#fit.gam <- fitdist(dataWind$spd, "gamma")       # not normal
#plot(fit.gam)
fit.wei <- fitdist(dataWind$spd, "weibull")       # not normal
plot(fit.wei)
#fit.wei <- fitdist(dataWind$spd, "lnorm")       # not normal
#plot(fit.wei)

# transform the speed data from weibull to normal as weibull residuals not normal
require(Johnson)
#spd_john <- RE.Johnson(dataWind$spd)
spd_john <- RE.Johnson(decmp_spd$remainder)
hist(spd_john$transformed, breaks = 25, col = rgb(0.9, 0.1, 0.1, 0.5))
fit.norm <- fitdist(spd_john$transformed, "norm")
plot(fit.norm)
qqnorm(spd_john$transformed)
qqline(spd_john$transformed, col = "red")

ts_spd <- ts(spd_john$transformed, frequency=365*4, start = c(as.numeric(yr),as.numeric(doy)))

# decompose series
require(stlplus)
decmp_sst <- stlplus(ts_sst, s.window="periodic")[[1]]
decmp_uwind <- stlplus(ts_uwind, s.window="periodic")[[1]]
decmp_vwind <- stlplus(ts_vwind, s.window="periodic")[[1]]
decmp_spd <- stlplus(ts_spd, s.window="periodic")[[1]]

fit.norm <- fitdist(decmp_sst$remainder, "norm")       # normal
plot(fit.norm)
plot(stlplus(ts_sst, s.window="periodic"),main="SST Decomp")

fit.norm <- fitdist(decmp_uwind$remainder, "norm")     # normal
plot(fit.norm)
plot(stlplus(ts_uwind, s.window="periodic"),main="U Wind Decomp")

fit.norm <- fitdist(decmp_vwind$remainder, "norm")     # normal
plot(fit.norm)
plot(stlplus(ts_vwind, s.window="periodic"),main="V Wind Decomp")

fit.norm <- fitdist(decmp_spd$remainder, "norm")       # normal
plot(fit.norm)
plot(stlplus(ts_spd, s.window="periodic"),main="Wind Speed (Trans.) Decomp")

# 

