# 2_netCDF2csv.R

library(ncdf4)
library(data.table)
library(plyr)
library(tidyverse)
library(reshape2)
library(lubridate)
library(stringr)
library(raster)
library(mgcv)
library(doMC); doMC::registerDoMC(cores = 4)


# Set bbox ----------------------------------------------------------------

bbox <- data.frame(SA = c(-39, -25, 12, 37), # All of South Africa
                   AC = c(-45, -20, 6.25, 45), # Agulhas Current
                   # BC = c(-45, -5, -60, -25), # Brazil Current (works with gshhs)
                   BC = c(-45, -5, 300, 335), # Brazil Current (works with SST data)
                   EAC = c(-42.5, -15, 145, 160), # East Australian Current
                   GS = c(20, 50, 270, 320), # Gulf Stream (or 270-360 and 320-360 for AVISO+ data)
                   # KC = c(5, 60, 120, 180), # Kuroshio Current (wide)
                   KC = c(20, 45, 120, 175), # Kuroshio Current (narrow)
                   row.names = c("latmin", "latmax", "lonmin", "lonmax"))

# Read OISST --------------------------------------------------------------
#          1         2         3
# 123456789012345678901234567890
# avhrr-only-v2.19810901.nc

# function to extract the dims and data from OISST netCDFs
read_nc <- function(ncFile, region = region, csvDir = csvDir) {
  coords <- bbox[, region]
  nc <- nc_open(ncFile)
  pathLen <- nchar(OISST.dir) + 1 # to account for the "/" that needs to be inserted
  fNameStem <-
    substr(ncFile, pathLen + 1, pathLen + 13)
  fDate <- substr(ncFile, pathLen + 15, pathLen + 22)
  LatIdx <- which(nc$dim$lat$vals > coords[1] & nc$dim$lat$vals < coords[2])
  LonIdx <- which(nc$dim$lon$vals > coords[3] & nc$dim$lon$vals < coords[4])
  sst <- ncvar_get(nc,
                   varid = "sst",
                   start = c(LonIdx[1], LatIdx[1], 1, 1),
                   count = c(length(LonIdx), length(LatIdx), 1, 1)) %>%
    round(3)
  dimnames(sst) <- list(lon = nc$dim$lon$vals[LonIdx],
                        lat = nc$dim$lat$vals[LatIdx])
  nc_close(nc)
  sst <-
    as.data.frame(melt(sst, value.name = "temp"), row.names = NULL) %>%
    mutate(t = ymd(fDate)) %>%
    na.omit()
  fwrite(sst,
         file = paste(csvDir, "/", region, "-", fNameStem, ".", strtDate, "-", endDate, ".csv", sep = ""),
         append = TRUE)
  rm(sst)
}

# setup OISST
  # NB: This directory will be different on other computers
  # It is setup to run off of a AJ's hard drive
OISST.dir <- "/media/rws/Benguela/OceanData/OISSTv2/daily/netCDF"
csv.dir <- "data"

# the list of files
ncList <- list.files(path = OISST.dir, pattern = "*.nc", full.names = TRUE, include.dirs = TRUE)
# ncList <- ncList[1:10] # Tester...
strtDate <- str_sub(ncList[1], start = 66, end = 73)
endDate <- str_sub(ncList[length(ncList)], start = 66, end = 73)

# apply the function
system.time(llply(ncList, read_nc, region = "SA", csvDir = csv.dir, .parallel = TRUE))
# system.time(llply(ncList, read_nc, region = "AC", csvDir = OISST.csv.dir, .parallel = TRUE))
# system.time(llply(ncList, read_nc, region = "BC", csvDir = OISST.csv.dir, .parallel = TRUE))
# system.time(llply(ncList, read_nc, region = "EAC", csvDir = OISST.csv.dir, .parallel = TRUE))
# system.time(llply(ncList, read_nc, region = "KC", csvDir = OISST.csv.dir, .parallel = TRUE))
# system.time(llply(ncList, read_nc, region = "GS", csvDir = OISST.csv.dir, .parallel = TRUE))


# Apply shapefile ---------------------------------------------------------

# OISST
OISST <- fread("data/SA-avhrr-only-v2.19810901-20171231.csv")
colnames(OISST) <- c("lon", "lat", "temp", "date")

load("data/200nm_mask.RData")
pf <- fortify(p)
# ggplot(pf, aes(x = long, y = lat)) +
#   borders() +
#   geom_path(aes(group = group)) +
#   coord_equal(expand = c(0,0), xlim = c(12, 37), ylim = c(-25, - 39))

bound <- list(list(lon = pf$long, lat = pf$lat))

inbound <- with(OISST, inSide(bound, lon, lat))

OISST_inbound <- OISST[inbound,]
ggplot(filter(OISST_inbound, date == "1997-07-18"), aes(x = lon, y = lat)) +
  borders() +
  geom_raster(aes(fill = temp)) +
  coord_equal(expand = c(0,0), xlim = c(12, 37), ylim = c(-25, - 39))


# Read CMC SST ------------------------------------------------------------

# function to extract the dims and data from CMC netCDFs
read_nc <- function(ncDir = ncDir, region = region, csvDir = csvDir) {
  ncList <- list.files(path = paste0(ncDir, "/", region), pattern = "*.nc", full.names = TRUE, include.dirs = TRUE)
  ncFirst <- head(list.files(path = paste0(ncDir, "/", region), pattern = "*.nc", full.names = FALSE), 1)
  ncLast <- tail(list.files(path = paste0(ncDir, "/", region), pattern = "*.nc", full.names = FALSE), 1)
  strtDate <- str_sub(ncFirst, start = 1, end = 8)
  endDate <- str_sub(ncLast, start = 1, end = 8)
  ncFun <- function(ncFile = ncFile, region = region, csvDir = csvDir) {
    nc <- nc_open(ncFile)
    pathLen <- nchar(paste0(ncDir, "/", region)) + 1
    fNameStem <-
      substr(ncFile, pathLen + 37, pathLen + 56)
    fDate <- substr(ncFile, pathLen + 1, pathLen + 8)
    sst <- ncvar_get(nc, varid = "analysed_sst") %>%
      round(4)
    dimnames(sst) <- list(lon = nc$dim$lon$vals,
                          lat = nc$dim$lat$vals)
    nc_close(nc)
    sst <-
      as.data.frame(melt(sst, value.name = "temp"), row.names = NULL) %>%
      mutate(t = ymd(fDate)) %>%
      na.omit()
    fwrite(sst,
           file = paste0(csvDir, "/", region, "-", fNameStem, "-", strtDate, "-", endDate, ".csv"),
           append = TRUE, col.names = FALSE)
    rm(sst)
  }
  llply(ncList, ncFun, region = region, csvDir = csvDir, .parallel = TRUE)
}

ncDir <- "/Volumes/Challenger_Deep/OceanData/MW_OI-REMSS-L4-GLOB-v4.0/netCDF"
MW_OI.v4.0.csv.dir <- "/Volumes/Challenger_Deep/spatial/processed/MW_OI-REMSS-L4-GLOB-v4.0/WBC/daily"
# ncDir <- "/Users/ajsmit/Dropbox/repos/temp/CMC0.2deg-CMC-L4-GLOB-v2.0"
# CMC.out.dir <- "/Users/ajsmit/temp/CMC0.2deg-CMC-L4-GLOB-v2.0"

system.time(read_nc(ncDir, "AC", MW_OI.v4.0.csv.dir))
# repeat GS extraction...

ncDir <- "/Volumes/Challenger_Deep/OceanData/CMC0.2deg-CMC-L4-GLOB-v2.0/netCDF"
CMC.csv.dir <- "/Volumes/Challenger_Deep/spatial/processed/CMC0.2deg-CMC-L4-GLOB-v2.0/WBC/daily"

system.time(read_nc(ncDir, "AC", CMC.csv.dir))



# Read MW_OI-REMSS-L4-GLOB-v4.0 -------------------------------------------


# function to extract the dims and data from netCDFs
read_nc <- function(ncDir = ncDir, region = region, csvDir = csvDir) {
  ncList <- list.files(path = paste0(ncDir, "/", region), pattern = "*.nc", full.names = TRUE, include.dirs = TRUE)
  ncFirst <- head(list.files(path = paste0(ncDir, "/", region), pattern = "*.nc", full.names = FALSE), 1)
  ncLast <- tail(list.files(path = paste0(ncDir, "/", region), pattern = "*.nc", full.names = FALSE), 1)
  strtDate <- str_sub(ncFirst, start = 1, end = 8)
  endDate <- str_sub(ncLast, start = 1, end = 8)
  ncFun <- function(ncFile = ncFile, region = region, csvDir = csvDir) {
    nc <- nc_open(ncFile)
    pathLen <- nchar(paste0(ncDir, "/", region)) + 1
    fNameStem <- "MW_OI-REMSS-L4-GLOB-v4.0"
    fDate <- substr(ncFile, pathLen + 1, pathLen + 8)
    sst <- ncvar_get(nc, varid = "analysed_sst") %>%
      round(4)
    dimnames(sst) <- list(lon = nc$dim$lon$vals,
                          lat = nc$dim$lat$vals)
    nc_close(nc)
    sst <-
      as.data.frame(melt(sst, value.name = "temp"), row.names = NULL) %>%
      mutate(t = ymd(fDate)) %>%
      na.omit()
    fwrite(sst,
           file = paste0(csvDir, "/", region, "-", fNameStem, "-", strtDate, "-", endDate, ".csv"),
           append = TRUE, col.names = FALSE)
    rm(sst)
  }
  llply(ncList, ncFun, region = region, csvDir = csvDir, .parallel = TRUE)
}

ncDir <- "/Volumes/Benguela/OceanData/MW_OI-REMSS-L4-GLOB-v4.0/netCDF"
csvDir <- "/Volumes/Benguela/spatial/processed/MW_OI-REMSS-L4-GLOB-v4.0/WBC/daily"

system.time(read_nc(ncDir, "AC", csvDir))
