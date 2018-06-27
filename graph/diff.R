
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(FNN)


# Load data ---------------------------------------------------------------

## In situ data
load("data/SACTN_daily_v4.2_sans_QC.RData")
load("data/SACTN_monthly_v4.2.RData")
load("data/SACTN_yearly_v4.2.RData")
load("data/SACTN_clim_v4.2.RData")

## Site list
load("setup/site_list_v4.2_sans_QC.RData")
site_coords <- site_list %>%
  dplyr::select(lon, lat, index)

## OISST
load("data/OISSt_inbound.RData")

## Bathymetry
load("data/bathy_mask.Rdata")
bathy_mask_sub <- bathy_mask %>% 
  filter(depth <= 0)
rm(bathy_mask)

# Create means ------------------------------------------------------------

breaks <- seq(8,28,2)

pal_1 <- c("#49AEB5", "#41B0AC", "#3EB1A3", "#41B198", "#48B28C", "#52B181",
           "#5EB175", "#6BAF69", "#78AE5E", "#85AC54", "#92A94B", "#9FA643",
           "#ACA23E", "#B99E3C", "#C5993C", "#D09440", "#DA8F45", "#E38A4D", "#EB8456")

OISST_cut <- OISST_inbound %>% 
  dplyr::filter(month %in% c("01", "02", "03", "07", "08", "09")) %>% 
  dplyr::mutate(season = if_else(month %in% c("01", "02", "03"), "Summer", "Winter")) %>% 
  dplyr::group_by(lon, lat, season) %>% 
  dplyr::summarise(mean_temp = mean(temp, na.rm = T)) %>%
  dplyr::mutate(temp_cut = cut(mean_temp, breaks = breaks)) %>% 
  dplyr::ungroup()
# OISST_cut <- as.data.frame(OISST_cut)

SACTN_cut <- SACTN_daily_v4.2 %>% 
  dplyr::mutate(month = month(date)) %>% 
  dplyr::filter(month %in% c(1, 2, 3, 7, 8, 9)) %>% 
  dplyr::mutate(season = if_else(month %in% c(1, 2, 3), "Summer", "Winter")) %>% 
  dplyr::group_by(index, season) %>% 
  dplyr::summarise(mean_temp = mean(temp, na.rm = T)) %>%
  dplyr::mutate(temp_cut = cut(mean_temp, breaks = breaks)) %>% 
  left_join(site_coords, by = "index") %>% 
  dplyr::ungroup() %>% 
  dplyr::select(lon, lat, everything())
# SACTN_cut <- as.data.frame(SACTN_cut)

# Calculate differences ---------------------------------------------------

# NB: Lon/lat must be the first and second columns of the dataframe in that order
df_is <- SACTN_cut
df_sst <- OISST_cut
product_compare <- function(df_is, df_sst){
  
  # Spread the in situ data by season
  df_is_wide <- df_is %>% 
    dplyr::select(-temp_cut) %>% 
    tidyr::spread(key = season, value = mean_temp)
  
  # Spread the sst data by season
  df_sst_wide <- df_sst %>% 
    dplyr::select(-temp_cut) %>% 
    tidyr::spread(key = season, value = mean_temp)
  
  # Find the nearest three matching pixels
  match_index <- knnx.index(data =  as.matrix(df_sst_wide[,1:2]),
                            query = as.matrix(df_is_wide[,1:2]), k = 3)
  
  # Grab the winter temperatures
  df_sst_winter <- data.frame(sub_1 = df_sst_wide$Winter[match_index[,1]],
                           sub_2 = df_sst_wide$Winter[match_index[,2]],
                           sub_3 = df_sst_wide$Winter[match_index[,3]]) %>% 
    dplyr::mutate(sub = rowMeans(.))
  
  # Grab the summer temperatures
  df_sst_summer <- data.frame(sub_1 = df_sst_wide$Summer[match_index[,1]],
                              sub_2 = df_sst_wide$Summer[match_index[,2]],
                              sub_3 = df_sst_wide$Summer[match_index[,3]]) %>% 
    dplyr::mutate(sub = rowMeans(.))
  
  # Join in the mean value for the nearest three pixels
  df_is_res <- df_is_wide %>% 
    dplyr::mutate(Winter_nearest = df_sst_winter$sub,
                  Summer_nearest = df_sst_summer$sub,
                  Winter_dif = Winter - Winter_nearest,
                  Summer_dif = Summer - Summer_nearest) %>% 
    dplyr::select(-c(Summer:Summer_nearest)) %>% 
    gather(key = season, value = dif, 
           -lon, -lat, -index) %>% 
    dplyr::mutate(cut = cut(dif, breaks = seq(-8, 1, 1)))
  
  # Return the final product
  return(df_is_res)
}

SACTN_dif <- product_compare(SACTN_cut, OISST_cut)
SACTN_dif$season[SACTN_dif$season == "Winter_dif"] <- "Winter"
SACTN_dif$season[SACTN_dif$season == "Summer_dif"] <- "Summer"

# Plot --------------------------------------------------------------------

# Comparing Winter and SUmmer SACTN and OISST means
dif_1 <- ggplot(OISST_cut, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = temp_cut)) +
  geom_contour(data = bathy_mask_sub, aes(z = round(depth, -2)), 
               colour = "grey20", alpha = 0.4) +
  borders(fill = "grey70") +
  geom_point(data = SACTN_cut, size = 3.87, 
             shape = 16, colour = "white") +
  geom_point(data = SACTN_cut, size = 2.36, 
             aes(colour = temp_cut), show.legend = F) +
  facet_wrap(~season) +
  coord_equal(expand = 0, xlim = c(12, 37), ylim = c(-25, -38)) +
  # scale_fill_viridis_d() +
  # scale_colour_viridis_d() +
  scale_fill_brewer("Temp. (C)", palette = "Spectral", direction = -1, breaks = breaks) +
  scale_colour_brewer(palette = "Spectral", direction = -1) +
  # scale_colour_gradientn(colours = rainbow(13)) +
  # scale_colour_manual(values = pal_1) +
  # scale_fill_manual("Temp. (C)", values = pal_1) +
  labs(x = "", y = "") +
  theme(legend.direction = "horizontal",
        legend.position = "bottom")
dif_1

# Comparing the difference between coastal pixels
dif_2 <- ggplot(SACTN_dif, aes(x = lon, y = lat)) +
  geom_raster(data = OISST_cut, aes(fill = temp_cut)) +
  geom_contour(data = bathy_mask_sub, aes(z = round(depth, -2)), 
               colour = "grey20", alpha = 0.4) +
  borders(fill = "grey70") +
  geom_point(size = 3.87, shape = 16, colour = "white") +
  geom_point(size = 2.36, aes(colour = dif), show.legend = F) +
  facet_wrap(~season) +
  coord_equal(expand = 0, xlim = c(12, 37), ylim = c(-25, -38)) +
  # scale_fill_viridis_d() +
  # scale_colour_viridis_d() +
  scale_fill_brewer("Temp. (C)", palette = "Spectral", direction = -1, breaks = breaks) +
  # scale_colour_brewer(palette = "Spectral", direction = -1) +
  scale_colour_gradient2(low = "blue", high = "red") +
  # scale_colour_manual(values = pal_1) +
  # scale_fill_manual("Temp. (C)", values = pal_1) +
  labs(x = "", y = "") +
  theme(legend.direction = "horizontal",
        legend.position = "bottom")
dif_2

ggarrange(dif_1, dif_2, nrow = 2, ncol = 1, common.legend = TRUE)

