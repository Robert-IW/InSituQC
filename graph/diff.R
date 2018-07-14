
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(lubridate)
library(FNN)
library(akima)


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
bathy_mask_shelf <- bathy_mask_sub %>% 
  filter(depth <= -150, depth >= -250)
bathy_mask_slope <- bathy_mask_sub %>% 
  filter(depth <= -500, depth >= -1000) %>% 
  mutate(depth = round(depth, -2))
bathy_mask_plain <- bathy_mask_sub %>% 
  filter(depth <= -1000, depth >= -6000) %>% 
  mutate(depth = round(depth, -3))

## More bathymetry
bathy <- as.matrix(read.table("data/gebco_sa.asc", sep = "", skip = 6)) # Load in bathymetry plot
info <- read.table("data/gebco_sa.asc", sep = "", nrows = 6) # Load in plot size parameters

xrng <- seq(from = info[3,2], by = info[5,2], length.out = info[1,2]) # Create x range integers
yrng <- seq(from = info[4,2], by = info[5,2], length.out = info[2,2]) # Create y range integers
xy <- expand.grid(xrng, yrng)
bathy[bathy == info[6,2]] <- NA # Renames error values with NA
bathy[bathy > 0] <- NA # Removes any values above the surface of the water
bathy <-  t(bathy[nrow(bathy):1,]) # Adds a column giving row names
bathy <- as.vector(reshape2::melt(bathy, value.name = "z", na.rm = FALSE)$z)
bathy <- as.data.frame(cbind(xy, bathy))
colnames(bathy) <- c("lon", "lat", "depth"); rm(xy)


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

# Join the different dif outputs for plotting
SACTN_combi <- SACTN_cut %>% 
  left_join(SACTN_dif, join = c("lon", "lat", "index", "season")) %>% 
  na.omit() %>% # Remove Hout Bay
  filter(index != "Dyer Island/DEA")


# Interpolate SACTN -------------------------------------------------------

load("data/sa_coast.Rdata")
# This warning may be ignored
res_winter <- data.frame(interpp(x = filter(SACTN_combi, season == "Winter")$lon, 
                                 y = filter(SACTN_combi, season == "Winter")$lat, 
                                 z = filter(SACTN_combi, season == "Winter")$dif,
                                 xo = sa_coast$X, yo = sa_coast$Y, linear = TRUE, 
                                 extrap = FALSE, dupl = "mean"))
res_winter$index <- 1:nrow(res_winter)
res_winter$x <- NULL
res_winter$y <- NULL
res_winter_approx <- data.frame(approx(x = res_winter$z, method = "linear", 
                            xout = res_winter$index[is.na(res_winter$z)]))
colnames(res_winter_approx) <- c("index", "z")
res_winter_approx <- res_winter_approx %>% 
  dplyr::select(index, z)

res_winter_full <- rbind(res_winter, res_winter_approx) %>% 
  na.omit() %>% 
  arrange(index) %>% 
  mutate(season = "Winter")

# Now the same for summer
res_summer <- data.frame(interpp(x = filter(SACTN_combi, season == "Summer")$lon, 
                                 y = filter(SACTN_combi, season == "Summer")$lat, 
                                 z = filter(SACTN_combi, season == "Summer")$dif,
                                 xo = sa_coast$X, yo = sa_coast$Y, linear = TRUE, 
                                 extrap = FALSE, dupl = "mean"))
res_summer$index <- 1:nrow(res_summer)
res_summer$x <- NULL
res_summer$y <- NULL
res_summer_approx <- data.frame(approx(x = res_summer$z, method = "linear", 
                                       xout = res_summer$index[is.na(res_summer$z)]))
colnames(res_summer_approx) <- c("index", "z")
res_summer_approx <- res_summer_approx %>% 
  dplyr::select(index, z)

res_summer_full <- rbind(res_summer, res_summer_approx) %>% 
  na.omit() %>% 
  arrange(index) %>% 
  mutate(season = "Summer")

res_season_full <- rbind(res_winter_full, res_summer_full) %>% 
  spread(key = season, value = z)

# Plot --------------------------------------------------------------------

# Comparing Winter and SUmmer SACTN and OISST means
dif_1 <- ggplot(OISST_cut, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = temp_cut), show.legend = T) +
  # stat_contour(data = bathy, aes(x = lon, y = lat, z = depth, alpha = ..level..), 
  #              col = "grey70", size = 0.2, binwidth = 500, show.legend = F) +
  # geom_contour(data = bathy_mask_shelf, aes(z = depth), 
  #              colour = "grey20", alpha = 0.4, bins = 1) +
  borders(fill = "grey70") +
  # geom_point(data = SACTN_dif, size = 3.87, 
  #            shape = 16, aes(colour = dif)) +
  # geom_point(data = SACTN_combi, size = 2.36, shape = 21, 
  #            aes(fill = temp_cut, colour = dif), stroke = 1) +
  geom_point(data = SACTN_combi, size = 3, shape = 21, 
             aes(fill = temp_cut), colour = "white", stroke = 0.5, show.legend = F) +
  # geom_point(data = SACTN_combi, aes(colour = dif), 
  #            shape = "+", size = 4) +
  facet_wrap(~season, ncol = 2) +
  coord_equal(expand = 0, xlim = c(12, 37), ylim = c(-25, -38)) +
  # scale_fill_viridis_d() +
  # scale_colour_viridis_d() +
  scale_fill_brewer("Temp. (°C)", palette = "Spectral", direction = -1) +
  # scale_colour_brewer(palette = "Spectral", direction = -1) +
  # scale_colour_gradient2(low = "black", high = "white") +
  labs(x = "", y = "") +
  theme(legend.direction = "horizontal",
        legend.position = "bottom")
# dif_1
ggsave(plot = dif_1, filename = "graph/summer_winter_diff.png", height = 4, width = 10)
ggsave(plot = dif_1, filename = "graph/summer_winter_diff.pdf", height = 4, width = 10)

# Line graph showing difference along coast
dif_2 <- ggplot(data = res_season_full, aes(x = index)) +
  geom_ribbon(aes(ymax = Winter, ymin = Summer))
# dif_2

# ggarrange(dif_1, dif_2, ncol = 1, nrow = 2, heights = c(5,1), align = "h", axis = "l")
gridExtra::grid.arrange(dif_1, dif_2, layout_matrix = cbind(c(1,1,1,NA), c(1,1,1,2),
                                                            c(1,1,1,2), c(1,1,1,NA)))
cowplot::plot_grid(dif_1, dif_2, labels = NULL, ncol = 1, 
                   align = 'v', axis = 'l')
       