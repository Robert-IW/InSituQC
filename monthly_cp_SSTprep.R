# This code will extract all the individual station data for a given SST product and process
# it to a single data frame with columns 'index'= {Station Name/Source}, 'date'={1991-01-01},
# 'temp'={20.0000}

# The source files are in ;KELP-HDD-Portable/SST/GHRSTT/L4' and individual products {AVHRR-OI,
# CMC, G1SST, K10, MUR, MW_IR_v1.0, MW_IR_v4.0, ODYSSEA_SAF, OSTIA} and sub-directories '/ts/ESdata/'
# and file nameing convention 'MW-IR-v4_Betty's Bay_SST_timeseries_5nearest.csv

# Each csv file has the headings {date={19810101}, station={Betty's Bay}}, lon, lat, nearest1}

# OBJECTIVE: For each product, open each station csv and extract names as 'index', 'date' and 'temp'
# OUTPUT: A single data frame to be used by 'monthly_cp.R'

library(tidyverse)
library(dplyr)

data.dir1 <- "/media/robert/KELP-HDD-Portable/SST/GHRSST/L4/"
#data.dir2 <- "AVHRR-OI/"
#data.dir2 <- "MUR/"
#data.dir2 <- "MW_IR_v4.0/"
#data.dir2 <- "OSTIA/"
data.dir2 <- "CMC/"
data.dir3 <- "ts/ESdata/"

data.dir <- paste0(data.dir1, data.dir2, data.dir3)

setwd(data.dir)

dir.list <- as.list(list.files(data.dir))

extr.data <- function(x){
  a <- read.csv(x, stringsAsFactors = F)
  b <- a %>%
    dplyr::select(index=station, date, temp=nearest1) %>%
    mutate(date=as.Date(strptime(date, "%Y%m%d"), format="%Y%m%d"))
}
 
extr.site <-  function(x){
  a <- read.csv(x, stringsAsFactors = F)
  b <- a %>%
    dplyr::select(index=station, lon, lat) %>% 
    filter(row_number()==1)
}
  
temp.list <- lapply(dir.list, function(x) extr.data(x))
data.in <- do.call(rbind, temp.list)

site.list <- lapply(dir.list, function(x) extr.site(x))
site.loc <- do.call(rbind, site.list)

#write.csv(data.in, file = "~/R/Data/StationSat.csv")
save(data.in, file = "~/R/Data/StationSat.Rdata")
save(site.loc, file = "~/R/Data/StationSat_SiteList.Rdata")
