
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(data.table)
library(lubridate)
library(stlplus)


# Load data ---------------------------------------------------------------

load("data/SACTN_daily_v4.2_sans_QC.RData")
SACTN <- as.data.frame(SACTN_daily_v4.2)
rm(SACTN_daily_v4.2)

load("setup/site_list_v4.2_sans_QC.RData")


# Flag 3SD ----------------------------------------------------------------

# Calculare the 3SD range per site per month
monthly_means <- SACTN %>%
  mutate(month = month(date, abbr = T)) %>%
  group_by(index, month) %>%
  summarise(temp_mean = mean(temp, na.rm = T),
            sd = sd(temp, na.rm = T)) %>%
  ungroup() %>%
  mutate(bottom_limit = temp_mean-sd*3,
         top_limit = temp_mean+sd*3) %>%
  dplyr::select(index, month, bottom_limit, top_limit)

# Apply the flag
SACTN_flags <- SACTN %>% 
  mutate(month = month(date, abbr = T)) %>%
  left_join(monthly_means, by = c("index", "month")) %>% 
  mutate(flag_1 = ifelse(temp > top_limit | temp < bottom_limit, 1, 0)) %>%
  # filter(flag_1 == 1) # 4278
  dplyr::select(-top_limit, -bottom_limit, -month)


# Flag spikes -------------------------------------------------------------

# grads <- SACTN_flags %>% 
#   group_by(index) %>% 
#   mutate(grad = temp - lag(temp),
#          month = month(date),
#          year = year(date)) %>% 
#   # na.omit() %>% 
#   group_by(index, month, year) %>% 
#   summarise(mean = mean(grad, na.rm=T),
#             sd = sd(grad, na.rm=T))

# Tester...
# df <- SACTN %>%
# filter(index == "Mossel Bay/SAWS")
ts_resid <- function(df){
  df <- df %>% 
    mutate(month = month(date))
  
  doy <- format(df$date[1], "%j")
  yr <- format(df$date[1], "%Y")
  
  df_ts <- ts(df$temp, frequency=365, 
             start = c(as.numeric(yr), as.numeric(doy)))

  df_decmp <- stlplus(df_ts, s.window="periodic")[[1]]

  df$dif <- df_decmp$remainder - lag(df_decmp$remainder)

  month_dif <- df %>% 
    mutate(month = month(date)) %>% 
    group_by(month) %>% 
    summarise(sd_dif = sd(dif, na.rm = T)) %>% 
    mutate(top_limit = sd_dif + sd_dif*5,
           bot_limit = sd_dif - sd_dif*5)
  
  df1 <- df %>% 
    mutate(month = month(date)) %>% 
    left_join(month_dif, by = "month") %>% 
    mutate(flag_2 = ifelse(dif > top_limit | dif < bot_limit, 1, 0)) %>% 
    select(-sd_dif, -top_limit, -bot_limit)

# ggplot(df1, aes(x = date, y = temp)) +
#   geom_line() +
#   geom_vline(data = filter(df1, flag_2 == 1),
#              aes(xintercept = date), colour = "red")

return(df1)  
}

