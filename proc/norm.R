
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(data.table)
library(fasttime)
library(sp)


# Load and prep the data --------------------------------------------------

## In situ data
load("data/SACTN_daily_v4.2_sans_QC.RData")
SACTN <- SACTN_daily_v4.2 %>% 
  mutate(month = month(date))

load("setup/site_list_v4.2_sans_QC.RData")
site_coords <- site_list %>%
  dplyr::select(lon, lat, index)

## OISST
# OISST <- fread("data/SA-avhrr-only-v2.19810901-20171231.csv")
# colnames(OISST) <- c("lon", "lat", "temp", "date")
# OISST <- OISST %>% 
#   # mutate(month = lubridate::month(date)) #%>% 
#   mutate(month = format(as.Date(fastPOSIXct(date)), "%m"))
load("data/OISSt_inbound.RData")

## Mask
load("data/200nm_mask.RData")


# Normalisation function --------------------------------------------------

normalit<-function(m){
  (m - min(m, na.rm=T))/(max(m, na.rm=T)-min(m, na.rm=T))
}


# Insitu clustering function ----------------------------------------------

# df <- SACTN
is_k_means <- function(df){
  df.data <- data.table(df)[,.(mean = mean(temp, na.rm=T),
                               # med = median(temp, na.rm=T),
                               sd = sd(temp, na.rm=T),
                               max = max(temp, na.rm=T),
                               min = min(temp, na.rm=T)),
                            by = .(index, month)]

  # group by index and month and create summary statistics
  df.data.kmeans <- df.data %>% 
    gather(variable, value, -(index:month)) %>%
    unite(temp, month, variable) %>%
    spread(temp, value) %>% 
    na.omit()
  
  # Calculate k-means
  df.data.kmeans.nor <- df.data.kmeans[,-1] %>%
    mutate_all(funs(normalit)) %>% 
    mutate(cluster = kmeans(., 6, nstart=100)$cluster,
           index = df.data.kmeans$index)%>% 
    left_join(site_coords, by = "index")
  
  return(df.data.kmeans.nor)
}


# SST clustering function -------------------------------------------------

# df <- OISST_inbound
sst_k_means <- function(df){
  df.data <- data.table(df)[,.(mean = mean(temp, na.rm=T),
                               # med = median(temp, na.rm=T),
                               sd = sd(temp, na.rm=T),
                               max = max(temp, na.rm=T),
                               min = min(temp, na.rm=T)),
                            by = .(lon, lat, month)]
  
  # group by index and month and create summary statistics
  df.data.kmeans <- df.data %>% 
    gather(variable, value, -(lon:month)) %>%
    unite(temp, month, variable) %>%
    spread(temp, value) %>% 
    na.omit()
  
  # Calculate k-means
  df.data.kmeans.nor <- df.data.kmeans[,-c(1:2)] %>%
    mutate_all(funs(normalit)) %>% 
    mutate(cluster = kmeans(., 6, nstart=100)$cluster,
           lon = df.data.kmeans$lon,
           lat = df.data.kmeans$lat)
  
  return(df.data.kmeans.nor)
}


# k-means -----------------------------------------------------------------

SACTN_cluster <- is_k_means(SACTN)

OISST_cluster <- sst_k_means(OISST_inbound)

ggplot(data = OISST_cluster, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = as.factor(cluster)), show.legend = F) +
  borders(fill = "grey70") +
  geom_point(data = SACTN_cluster, size = 2.57, shape = 16, colour = "white") +
  geom_point(data = SACTN_cluster, size = 2.57, shape = 16, #colour = "white",
             aes(colour = as.factor(cluster)), show.legend = F) +
  geom_polygon(aes(group = cluster)) +
  coord_equal(expand = 0, xlim = c(12, 37), ylim = c(-25, -38)) +
  scale_colour_viridis_d(option = "C") +
  scale_fill_viridis_d(option = "D")
  

