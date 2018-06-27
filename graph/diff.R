
# Load libraries ----------------------------------------------------------

library(tidyverse)
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


# Create means ------------------------------------------------------------

breaks <- seq(10,28,2)

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
  dplyr::ungroup()
# SACTN_cut <- as.data.frame(SACTN_cut)

# Calculate differences ---------------------------------------------------

match_index <- knnx.index(data =  as.matrix(OISST_cut[,1:2]),
                          query = as.matrix(SACTN_cut[,5:6]), k = 3)
OISST_sub <- data.frame(sub_1 = OISST_cut$mean_temp[match_index[,1]],
                        sub_2 = OISST_cut$mean_temp[match_index[,2]],
                        sub_3 = OISST_cut$mean_temp[match_index[,3]])
OISST_sub <- OISST_sub %>% 
  mutate(sub = sum(sub_1))


# Plot --------------------------------------------------------------------

# Comparing Winter and SUmmer SACTN and OISST means
ggplot(OISST_cut, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = temp_cut)) +
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

