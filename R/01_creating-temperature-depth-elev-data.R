## this script wrangles air and sea surface temperature data, as well as elevation and bathymetry data
library(ncdf4)
library(tidyverse)
library(raster)
library(evobiR)
library(sf)
library(abind)
#library(devtools)
#install_github("adamlilith/enmSdm")
select <- dplyr::select
mean <- base::mean
rasterOptions(maxmemory = 1000000000000)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######       MAKING HIGH AND LOW TEMPERATURE RASTERS     #######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
###### Terrestrial seasonal temperature high and low #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## max and min terrestrial temps:
filename <- paste("large-files/air-temperatures/Complete_TMAX_Daily_LatLong1_1950.nc", sep = "")
ncfile <- nc_open(filename)

## create variables for things needed to use data
lat <- ncvar_get(ncfile, "latitude")
long <- ncvar_get(ncfile, "longitude")

## close the file
nc_close(ncfile)

## create arrays to hold mean temperatures from each 10 year dataset
mean_max <- mean_min <- mean <- array(dim = c(360, 180, 3650))
mean_max[,,] <- NaN
mean_min[,,] <- NaN
mean[,,] <- NaN

rep = 1950
while (rep < 2020) {
  ## open max and minimum files 
  filename_max <- paste("large-files/air-temperatures/Complete_TMAX_Daily_LatLong1_", rep, ".nc", sep = "")
  ncfile_max <- nc_open(filename_max)
  filename_min <- paste("large-files/air-temperatures/Complete_TMIN_Daily_LatLong1_", rep, ".nc", sep = "")
  ncfile_min <- nc_open(filename_min)
  filename_avg <- paste("large-files/air-temperatures/Complete_TAVG_Daily_LatLong1_", rep, ".nc", sep = "")
  ncfile_avg <- nc_open(filename_avg)
  
  ## create variables for data
  date <- ncvar_get(ncfile_max, "date_number")
  arr.anom_max <- ncvar_get(ncfile_max, "temperature")
  arr.clim_max <- ncvar_get(ncfile_max, "climatology")
  arr.anom_min <- ncvar_get(ncfile_min, "temperature")
  arr.clim_min <- ncvar_get(ncfile_min, "climatology")
  arr.anom_tavg <- ncvar_get(ncfile_avg, "temperature")
  arr.clim_tavg <- ncvar_get(ncfile_avg, "climatology")
  
  ## close files
  nc_close(ncfile_max)
  nc_close(ncfile_min)
  nc_close(ncfile_avg)
  
  ## figure out which years in this time frame are leap years
  leap_years <- seq(from = rep, to = (rep+9), by = 1) %% 4 == 0
  
  ## figure out which index represents the end of each year 
  x <- c(1)
  i = 1
  while(i < length(leap_years)) {
    if(leap_years[i] == FALSE) {
      x <- append(x, x[i] + 365)
    }
    else {
      x <- append(x, x[i] + 366)
    }
    i = i+1
  }
  
  leap_years <- data.frame(leap_year = leap_years, index = x)
  
  ## create arrays to store temperatures in 
  arr.temps_max <- array(dim = c(nrow(arr.clim_max), ncol(arr.clim_max), 3650))
  arr.temps_min <- array(dim = c(nrow(arr.clim_min), ncol(arr.clim_min), 3650))
  arr.temps_tavg <- array(dim = c(nrow(arr.clim_tavg), ncol(arr.clim_tavg), 3650))
  
  ## loop through each element (unique pairs of row x column)
  row =  1
  while(row < nrow(arr.anom_max) + 1) {
    col = 1
    while(col < ncol(arr.anom_max) + 1) {
      print(paste("On rep number: ", rep, ", row: ", row, ', col:', col, sep = "", " :-)"))
      ## retrieve climatology and anomaly in the cell 
      anom_max <- arr.anom_max[row, col, ]
      anom_min <- arr.anom_min[row, col, ]
      anom_tavg <- arr.anom_tavg[row, col, ]
      
      year <- 1
      this_index <- leap_years
      
      ## for leap years, remove element 60 within that year in anomaly data 
      while(year < nrow(leap_years) + 1) {
        if (this_index$leap_year[year] == TRUE) {
          index <- this_index$index[year]
          anom_max <- anom_max[-index+60]
          anom_min <- anom_min[-index+60]
          anom_tavg <- anom_tavg[-index+60]
          
          this_index$index <- this_index$index - 1
        }
        year = year + 1
      }
      ## extend climatology and add anomaly to climatology to get temperature 
      clim_max <- rep(arr.clim_max[row, col, ], times = 10)
      temps_max <- clim_max + anom_max
      clim_min <- rep(arr.clim_min[row, col, ], times = 10)
      temps_min <- clim_min + anom_min
      clim_tavg <- rep(arr.clim_tavg[row, col, ], times = 10)
      temps_tavg <- clim_tavg + anom_tavg
      
      ## store temperatures for that cell 
      arr.temps_max[row, col, ] <- temps_max
      arr.temps_min[row, col, ] <- temps_min
      arr.temps_tavg[row, col, ] <- temps_tavg
      
      ## calculate average temp in each cell for each day over the period and store 
      x = 1
      iteration = (rep-1950)/10
      while (x < 366) {
        mean_max[row, col, x+365*iteration] <- mean(arr.temps_max[row, col, seq(x, 3650, 365)], na.rm = TRUE)
        mean_min[row, col, x+365*iteration] <- mean(arr.temps_min[row, col, seq(x, 3650, 365)], na.rm = TRUE)
        mean[row, col, x+365*iteration] <- mean(arr.temps_tavg[row, col, seq(x, 3650, 365)], na.rm = TRUE)
        x = x +1
      }
      col = col + 1
    }
    row = row + 1
  }
  rep = rep + 10
}

# saveRDS(mean_min, "large-files/air-temperatures/terrestrial_mean-min.rds")
# saveRDS(mean_max, "large-files/air-temperatures/terrestrial_mean-max.rds")
# saveRDS(mean, "large-files/air-temperatures/terrestrial_mean.rds")
mean_min <- readRDS("large-files/air-temperatures/terrestrial_mean-min.rds")
mean_max <- readRDS("large-files/air-temperatures/terrestrial_mean-max.rds")
mean <- readRDS("large-files/air-temperatures/terrestrial_mean.rds")

## create temperature set matrices and lists to store data in  
## add columns for temps species with cold and hot dormancy
mean_max_new <- mean_min_new <- mean_new <- array(dim = c(360, 180, 365))
final_max <- final_min <- final_avg <- list()

element = 1
row = 1
while(row < nrow(mean_max) + 1) {
  col = 1
  while (col < ncol(mean_max) +1) {
    x = 1
    while (x < 366) {
      ## get mean min and max temp for each cell over all 70 years:
      mean_max_new[row, col, x] <- mean(mean_max[row, col, seq(x, 3650, 365)], 
                                        na.rm = TRUE)
      mean_min_new[row, col, x] <- mean(mean_min[row, col, seq(x, 3650, 365)], 
                                        na.rm = TRUE)
      mean_new[row, col, x] <- mean(mean[row, col, seq(x, 3650, 365)], 
                                        na.rm = TRUE)
      x = x+1
    }
    ## get maximum mean daily temperature and minimum mean daily temperature for cell:
    max <- max(mean_max_new[row, col,], na.rm=TRUE)
    min <- min(mean_min_new[row, col,], na.rm=TRUE)
    avg <- mean(mean_new[row, col,], na.rm=TRUE)
    
    ## get mean temp over the 30 days before and after to min and max temps:
    max_index2 <- first(which(mean_max_new[row, col,] == max))
    min_index2 <- first(which(mean_min_new[row, col,] == min))
    
    acc_max <- append(mean_max_new[row, col,], mean_max_new[row, col,1:183]) 
    acc_min <- append(mean_min_new[row, col,], mean_min_new[row, col,1:183]) 
    max_index1 <- ifelse(max_index2 - 30 < 0, max_index2-30+183, max_index2-30)
    min_index1 <- ifelse(min_index2 - 30 < 0, min_index2-30+183, min_index2-30)
    
    if (!is.infinite(max) & !is.infinite(min)) {
      acc_max <- mean(acc_max[max_index1:max_index2], na.rm = TRUE)
      acc_min <- mean(acc_min[min_index1:min_index2], na.rm = TRUE)
    } 
    else {
      acc_max = acc_min = NA
    }
  
    ## create df for temps and add column for month and month index
    daily_highs <- data.frame(day = c(1:365), temp = mean_max_new[row, col,], 
                              month = factor(c(rep("Jan", 31),
                                               rep("Feb", 28),
                                               rep("Mar", 31),
                                               rep("Apr", 30),
                                               rep("May", 31),
                                               rep("Jun", 30),
                                               rep("Jul", 31),
                                               rep("Aug", 31),
                                               rep("Sep", 30),
                                               rep("Oct", 31),
                                               rep("Nov", 30),
                                               rep("Dec", 31)), levels = c("Jan", "Feb", "Mar",
                                                                           "Apr", "May", "Jun",
                                                                           "Jul", "Aug", "Sep", 
                                                                           "Oct", "Nov", "Dec")), 
                              month_index = c(rep(1, 31),
                                              rep(2, 28),
                                              rep(3, 31),
                                              rep(4, 30),
                                              rep(5, 31),
                                              rep(6, 30),
                                              rep(7, 31),
                                              rep(8, 31),
                                              rep(9, 30),
                                              rep(10, 31),
                                              rep(11, 30),
                                              rep(12, 31)))
    daily_lows <- daily_highs %>%
      mutate(temp = mean_min_new[row, col,])
    
    ## get most extreme temperature high and low temp from each month, and create 'circular' data frame
    monthly_highs <- daily_highs %>%
      group_by(month) %>%
      do(mutate(., temp = max(.$temp, na.rm = TRUE))) %>%
      select(-day) %>%
      unique(.) %>%
      rbind(., .[1:6,])
    
    monthly_lows <- daily_lows %>%
      group_by(month) %>%
      do(mutate(., temp = min(.$temp, na.rm = TRUE))) %>%
      select(-day) %>%
      unique(.) %>%
      rbind(., .[1:6,])
    
    ## find 6 months with the highest mean extreme daily high temp
    sw_high <- SlidingWindow(monthly_highs$temp, window = 6, FUN = sum, step = 1)
    
    ## find 6 months with the lowest mean extreme daily low temp
    sw_low <- SlidingWindow(monthly_lows$temp, window = 6, FUN = sum, step = 1)
    
    ## get index of months to block out
    hot_start <- which(sw_high == max(sw_high))[1]
    cold_start <- which(sw_low == min(sw_low))[1]
    
    hot_end <- hot_start + 5
    cold_end <- cold_start + 5
    
    ## get month indecies to block out
    if(hot_end > 12) {
      hot_end = hot_end - 12
      hot_indecies <- c(hot_start:12, 1:hot_end)
    }
    else {
      hot_indecies <- c(hot_start:hot_end)
    }
    
    if(cold_end > 12) {
      cold_end = cold_end - 12
      cold_indecies <- c(cold_start:12, 1:cold_end)
    }
    else {
      cold_indecies <- c(cold_start:cold_end)
    }
    
    ## block out temperatures in 6 hottest consecutive months and find max
    h_d_max_high <- daily_highs %>%
      filter(!month_index %in% hot_indecies) %>%
      select(temp) %>%
      max()
    
    h_d_min_low <- daily_lows %>%
      filter(!month_index %in% hot_indecies) %>%
      select(temp) %>%
      min()
    
    ## block out temperatures in 6 coldest consecutive months and find min
    c_d_min_low <- daily_lows %>%
      filter(!month_index %in% cold_indecies) %>%
      select(temp) %>%
      min()
    
    c_d_max_high <- daily_highs %>%
      filter(!month_index %in% cold_indecies) %>%
      select(temp) %>%
      max()
    
    ## calculate mean temp before this extreme:
    max_index2 <- first(which(mean_max_new[row, col,] == h_d_max_high))
    min_index2 <- first(which(mean_min_new[row, col,] == c_d_min_low))
    
    acc_max_dormancy <- append(mean_max_new[row, col,], mean_max_new[row, col,1:26]) 
    acc_min_dormancy <- append(mean_min_new[row, col,], mean_min_new[row, col,1:26]) 
    max_index1 <- ifelse(max_index2 - 3 < 0, max_index2-3+26, max_index2-3)
    min_index1 <- ifelse(min_index2 - 3 < 0, min_index2-3+26, min_index2-3)
    
    if (!is.infinite(max) & !is.infinite(min)) {
      acc_max_dormancy <- mean(acc_max_dormancy[max_index1:max_index2], na.rm = TRUE)
      acc_min_dormancy <- mean(acc_min_dormancy[min_index1:min_index2], na.rm = TRUE)
    } 
    else {
      acc_max_dormancy = acc_min_dormancy = NA
    }
    
    
    final_max[[element]] <- c(lat[col], long[row], max, h_d_max_high, c_d_max_high, 
                              acc_max, acc_max_dormancy)
    final_min[[element]] <- c(lat[col], long[row],  min, c_d_min_low, h_d_min_low, 
                              acc_min,acc_min_dormancy)
    final_avg[[element]] <- c(lat[col], long[row],  avg)
    
    print(paste("done column ", col, " of row ", row, " :-)", sep = ''))
    element = element + 1
    col = col + 1
  }
  row = row + 1
}

final_max <- data.frame(do.call(rbind, final_max), stringsAsFactors = FALSE) 
colnames(final_max) <- c("latitude", "longitude", "seasonal_high_temp",
                         "hot_dormancy_6mo", "cold_dormancy_6mo", "hot_acc_temp",
                         "hot_acc_temp_dormancy")
final_max <- filter(final_max, !is.infinite(seasonal_high_temp)) %>%
  select(longitude, latitude, seasonal_high_temp, hot_dormancy_6mo, cold_dormancy_6mo, 
         hot_acc_temp, hot_acc_temp_dormancy)

final_min <- data.frame(do.call(rbind, final_min), stringsAsFactors = FALSE) 
colnames(final_min) <- c("latitude", "longitude", "seasonal_low_temp",
                         "cold_dormancy_6mo", 'hot_dormancy_6mo', "cold_acc_temp", 
                         "cold_acc_temp_dormancy")
final_min <- filter(final_min, !is.infinite(seasonal_low_temp)) %>%
  select(longitude, latitude, seasonal_low_temp, cold_dormancy_6mo, hot_dormancy_6mo,
         cold_acc_temp, cold_acc_temp_dormancy)

final_avg <- data.frame(do.call(rbind, final_avg), stringsAsFactors = FALSE) 
colnames(final_avg) <- c("latitude", "longitude", "mean_temp")
final_avg <- filter(final_avg, !is.infinite(mean_temp)) %>%
  select(longitude, latitude, mean_temp)

## save data:
write.csv(final_max, "data-processed/intermediate-files/terrestrial_seasonal-max-temps_6mo-dormancy.csv", row.names = FALSE)
write.csv(final_min, "data-processed/intermediate-files/terrestrial_seasonal-min-temps_6mo-dormancy.csv", row.names = FALSE)
write.csv(final_avg, "data-processed/intermediate-files/terrestrial_avg-temp.csv", row.names = FALSE)

## is heat more variable at high latitudes?
first <- 1
second <- 60
sd <- array(dim = c(360,180,1))
while (second < 361) {
  rep = 1950
  while (rep < 2020) {
    ## open max files 
    filename_max <- paste("large-files/air-temperatures/Complete_TMAX_Daily_LatLong1_", rep, ".nc", sep = "")
    ncfile_max <- nc_open(filename_max)
    
    ## create variables for data
    arr.anom_max <- ncvar_get(ncfile_max, "temperature")[first:second, 1:180, ]
    
    ## close files
    nc_close(ncfile_max)
    
    if(rep == "1950") {
      temps <- arr.anom_max
    }
    else {
      temps <- abind(temps, arr.anom_max)
    }
    
    print(paste("On rep number:", rep))
    rep = rep + 10
  }
  
  ## calculate and store sd:
  sd[first:second, 1:180,] <- apply(temps[1:60,1:180,], c(1,2), sd, na.rm = T)
  
  first = first + 60
  second = second + 60
}

library(reshape2)
df <- melt(sd)[,-3]
colnames(df) <- c("x", "y", "sd")

r = rasterFromXYZ(df)
plot(r)

var_heat_map <- ggplot(df, aes(x = x, y = y, fill = sd)) + geom_raster() +
  coord_fixed() +
  theme_minimal() + 
  labs(x = "Latitude", y = "Longitude", 
       fill = "Standard deviation\nof hottest daily\nair temperature (°C)") +
  scale_fill_gradient2(mid = "white", high = "#b45346", low = 'steelblue', midpoint = 0,
                       na.value = "transparent") +
  scale_x_continuous(labels = c(-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150),
                     breaks = c(30, 60, 90, 120, 150,180,210,240,270,300,330)) + 
  scale_y_continuous(labels = c(-90, -60, -30, 0, 30, 60, 90),
                     breaks = c(0, 30, 60, 90,120,150,180)) +
  theme(legend.title = element_text(size = 9))

ggsave(var_heat_map, path = "figures/extended-data/", filename = "sd_heat_map.png", width = 9, height = 7)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######   Marine seasonal temperature high and low   #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in mean daily SST from 1982-2020
## make raster layers that keep track of of highest and lowest values and of acclimation temps
filenames <- c('sst.day.mean.1982.nc', 'sst.day.mean.1983.nc','sst.day.mean.1984.nc','sst.day.mean.1985.nc','sst.day.mean.1986.nc',
               'sst.day.mean.1987.nc', 'sst.day.mean.1988.nc','sst.day.mean.1989.nc','sst.day.mean.1990.nc','sst.day.mean.1991.nc',
               'sst.day.mean.1992.nc', 'sst.day.mean.1993.nc','sst.day.mean.1994.nc','sst.day.mean.1995.nc','sst.day.mean.1996.nc',
               'sst.day.mean.1997.nc', 'sst.day.mean.1998.nc','sst.day.mean.1999.nc','sst.day.mean.2000.nc','sst.day.mean.2001.nc',
               'sst.day.mean.2002.nc', 'sst.day.mean.2003.nc','sst.day.mean.2004.nc','sst.day.mean.2005.nc','sst.day.mean.2006.nc',
               'sst.day.mean.2007.nc', 'sst.day.mean.2008.nc','sst.day.mean.2009.nc','sst.day.mean.2010.nc','sst.day.mean.2011.nc',
               'sst.day.mean.2012.nc', 'sst.day.mean.2013.nc','sst.day.mean.2014.nc','sst.day.mean.2015.nc','sst.day.mean.2016.nc',
               'sst.day.mean.2017.nc', 'sst.day.mean.2018.nc','sst.day.mean.2019.nc','sst.day.mean.2020.nc')
file = 1
while(file < length(filenames)) {
  print(paste("On year ", file+1981, " :-)", sep = ""))
  
  ## read in:
  cur <- nc_open(paste("large-files/sea-surface-temperatures/", filenames[[file]], sep = ""))
  sst <- ncvar_get(cur, "sst")
  lat <- ncvar_get(cur, "lat")
  lon <- ncvar_get(cur, "lon")
  nc_close(cur)
  
  max_cur <- min_cur <- sst[,,1]
  max_cur[1:1440,1:720] <- apply(sst[1:1440,1:720,], c(1,2), max, na.rm = T)
  min_cur[1:1440,1:720] <- apply(sst[1:1440,1:720,], c(1,2), min, na.rm = T)

  ## if this raster has hotter/colder daily temperature, store and calculate acclimation temp
  if(file == 1) {
    max_sst = max_cur
    min_sst = min_cur
    acc_hot = max_cur
    acc_cold = min_cur
  }
  else {
    max_sst[max_sst < max_cur & !is.infinite(max_cur)] = max_cur[max_sst < max_cur & !is.infinite(max_cur)]
    min_sst[min_sst > min_cur & !is.infinite(min_cur)] = min_cur[min_sst > min_cur & !is.infinite(min_cur)] 
  }
  
  ## if the temps get updated, calculate the new acclimation temperature
  ## calculate acclimatization temp for each cell:
  x = 1
  while(x < 1441) {
    y=1
    while(y < 721) {
      if(max_sst[x,y] == max_cur[x,y] & !is.infinite(max_cur[x,y]))  {
        ## calc acclimation temp
        max_index2 <- which(sst[x,y,] == max_sst[x,y])
        max_index1 <- ifelse(max_index2 - 7 <= 0, max_index2-7+365, max_index2-7)
        
        ## compute the max daily air temperature 7 days before hottest day
        acc_max <- append(sst[x, y,], sst[x, y, 1:183]) 
        
        ## compute the max daily sea surface temperature 7 days before hottest day
        if (max_index1 > max_index2) {
          acc_max <- max(acc_max[append(c(max_index1:365),1:(max_index2-1))], na.rm = TRUE)
        }
        else {
          acc_max <- max(acc_max[max_index1:(max_index2-1)], na.rm = TRUE)
        }
        acc_hot[x,y] = acc_max
      }
      else if(min_sst[x,y] == min_cur[x,y] & !is.infinite(min_cur[x,y])) {
        ## calc acclimation temp
        min_index2 <- which(sst[x,y,] == min_sst[x,y])
        min_index1 <- ifelse(min_index2 - 49 <= 0, min_index2-49+365, min_index2-49)
        
        ## compute the max daily air temperature 7 days before hottest day
        acc_min <- append(sst[x, y,], sst[x, y, 1:183])
        
        ## compute the min daily sea surface temperature 49 days before coldest day
        if (min_index1 > min_index2) {
          acc_min <- min(acc_min[append(c(min_index1:365),c(1:min_index2-1))], na.rm = TRUE)
        }
        else {
          acc_min <- min(acc_min[min_index1:(min_index2-1)], na.rm = TRUE)
        }
        acc_cold[x,y] = acc_min
      }
      print(paste("On cell ", x, " ", y, " :-)", sep = ""))
      y = y + 1
    }
    x = x + 1
  }

  file = file + 1
}

## save: 
# saveRDS(max_sst, "data-processed/intermediate-files/max_sst.rds")
# saveRDS(min_sst, "data-processed/intermediate-files/min_sst.rds")
# saveRDS(acc_hot, "data-processed/intermediate-files/acc_hot.rds")
# saveRDS(acc_cold, "data-processed/intermediate-files/acc_cold.rds")

max_sst <- readRDS("data-processed/intermediate-files/max_sst.rds")
min_sst <- readRDS("data-processed/intermediate-files/min_sst.rds")
acc_hot <- readRDS("data-processed/intermediate-files/acc_hot.rds")
acc_cold <- readRDS("data-processed/intermediate-files/acc_cold.rds")

## get rid of infinities
max_sst[is.infinite(max_sst)] <- NA
min_sst[is.infinite(min_sst)] <- NA

## rasterize and aggregate  to 1 degree by 1 degree resolution:
r <- raster(nrow=720, ncol=1440, xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))

df <- melt(max_sst)[]
colnames(df) <- c("x", "y", "max_sst")
df$x <- rep(lon, n = 720)
df$y <- rep(lat, each = 1440)
df$x[which(df$x > 180)] <- df$x[which(df$x > 180)] - 360
rast = rasterFromXYZ(df)

sst_max_new <- aggregate(rast, fact = 4, fun = base::max, na.rm = TRUE)

df <- melt(min_sst)[]
colnames(df) <- c("x", "y", "min_sst")
df$x <- rep(lon, n = 720)
df$y <- rep(lat, each = 1440)
df$x[which(df$x > 180)] <- df$x[which(df$x > 180)] - 360
rast = rasterFromXYZ(df)

sst_min_new <- aggregate(rast, fact = 4, fun = base::min, na.rm = TRUE)

df <- melt(acc_hot)[]
colnames(df) <- c("x", "y", "acc_hot")
df$x <- rep(lon, n = 720)
df$y <- rep(lat, each = 1440)
df$x[which(df$x > 180)] <- df$x[which(df$x > 180)] - 360
rast = rasterFromXYZ(df)

hot_acc <- aggregate(rast, fact = 4, fun = base::max, na.rm = TRUE)

df <- melt(acc_cold)[]
colnames(df) <- c("x", "y", "acc_cold")
df$x <- rep(lon, n = 720)
df$y <- rep(lat, each = 1440)
df$x[which(df$x > 180)] <- df$x[which(df$x > 180)] - 360
rast = rasterFromXYZ(df)

cold_acc <- aggregate(rast, fact = 4, fun = base::min, na.rm = TRUE)

## turn into data frames and add empty columns for hot and cold dormancy
hot_acc <- data.frame(rasterToPoints(hot_acc)) 
colnames(hot_acc) <- c("longitude", "latitude", "hot_acc_temp")

final_max <- data.frame(rasterToPoints(sst_max_new)) 
colnames(final_max) <- c("longitude", "latitude", "seasonal_high_temp")

final_max <- left_join(final_max, hot_acc) %>%
  filter(!is.infinite(seasonal_high_temp)) %>%
  mutate(hot_dormancy_6mo = NA, cold_dormancy_6mo = NA, hot_acc_temp_dormancy = NA) %>%
  select(longitude, latitude, seasonal_high_temp, hot_dormancy_6mo, cold_dormancy_6mo, 
         hot_acc_temp, hot_acc_temp_dormancy)

cold_acc <- data.frame(rasterToPoints(cold_acc)) 
colnames(cold_acc) <- c("longitude", "latitude", "cold_acc_temp")

final_min <- data.frame(rasterToPoints(sst_min_new)) 
colnames(final_min) <- c("longitude", "latitude", "seasonal_low_temp")

final_min <- left_join(final_min, cold_acc) %>%
  filter(!is.infinite(seasonal_low_temp)) %>%
  mutate(cold_dormancy_6mo = NA, hot_dormancy_6mo = NA, cold_acc_temp_dormancy = NA) %>%
  select(longitude, latitude, seasonal_low_temp, cold_dormancy_6mo, hot_dormancy_6mo, cold_acc_temp, 
         cold_acc_temp_dormancy)

## save these datasets as seasonal high and low temps
write.csv(final_max, "data-processed/intermediate-files/marine_seasonal-max-temps_6mo.csv", row.names = FALSE)
write.csv(final_min, "data-processed/intermediate-files/marine_seasonal-min-temps_6mo.csv", row.names = FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######                   RESAMPLING ELEVATION                 #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
###### Resample minimum and maximum elevation in each grid cell   #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in EarthEnv data of minimum elevations at 10km resolution
minelev <- raster("data-raw/elevation-data/elevation_10KMmi_GMTEDmi.tif")
minelev[is.infinite(minelev)] <- NA
##plot(minelev, asp = 1)
## read in max
maxelev <- raster("data-raw/elevation-data/elevation_10KMma_GMTEDma.tiff")
maxelev[is.infinite(maxelev)] <- NA
##plot(maxelev, asp = 1)

## aggregate to 1deg x 1deg resolution, selecting minimum value of grid cells
minelev <- aggregate(minelev, fact = 12, fun = min, na.rm = TRUE)
##plot(minelev, asp = 1)
maxelev <- aggregate(maxelev, fact = 12, fun = max, na.rm = TRUE)
##plot(maxelev, asp = 1)

## crop raster to include only terrestrial areas using the terrestrial mask
t_mask <- raster("data-processed/masks/raster_terr_mask.grd")
minelev <- extend(minelev, t_mask) ## extend to same extent as temperature data 
minelev <- mask(minelev, t_mask) ## masks temps outside of the terrestrial mask 
maxelev <- extend(maxelev, t_mask) ## extend to same extent as temperature data 
maxelev <- mask(maxelev, t_mask) ## masks temps outside of the terrestrial mask 

## ## stack and write out to file:
elevs <- stack(minelev, maxelev)
names(elevs) <- c("elev_min", "elev_max")
writeRaster(elevs, "data-processed/masks/elev-and-bath-layers/elevs_minmax.grd", 
            overwrite = TRUE, format = "raster")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######                RESAMPLING BATHYMETRY               #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in bathymetry data and layer informing where land vs ocean occurs:
bath <- raster("large-files/bathymetric-data/gebco_2020_netcdf/GEBCO_2020.nc")
tid <- brick("large-files/bathymetric-data/gebco_2020_tid_netcdf/GEBCO_2020_TID.nc")
#plot(bath, asp = 1)
#plot(tid, asp = 1)

## set cells on land (0 in mask) to have elevation/depth of 0:
stack <- stack(bath, tid)

fun <- function(stack) {stack[[1]][stack[[2]] == 0] <- 0; return(stack[[1]]) }
bath <- raster::calc(stack, fun)
#plot(bath, asp = 1)

## save really long step:
writeRaster(bath, "large-files/bathymetric-data/bath_intermediate.grd", 
            overwrite = TRUE, format = "raster")
bath <- raster("large-files/bathymetric-data/bath_intermediate.grd")

tid <- aggregate(tid, fact = 240, fun = max, na.rm = TRUE) ## by assigning max value to aggregate, if aggregate area has any bit of ocean in it (value > 0) then the cell will not be considered land
tid[!tid == 0] <- 1 ## create land mask where land = 0, ocean = 1
#plot(tid, asp = 1)
writeRaster(tid, "data-processed/intermediate-files/TID_agg.grd", 
            overwrite = TRUE, format = "raster")

## aggregate data to same resolution as temperature data:
depth_max <- aggregate(bath, fact = 240, fun = max, na.rm = TRUE) ## assign cell a value of the maximum depth within the aggregate cells 
depth_max[tid == 0] <- NA
writeRaster(depth_max, "data-processed/intermediate-files/bath_agg_max.grd", 
            overwrite = TRUE, format = "raster")

depth_min <- aggregate(bath, fact = 240, fun = min, na.rm = TRUE) ## assign cell a value of the minimum depth within the aggregate cells 
depth_min[tid == 0] <- NA
writeRaster(depth_min, "data-processed/intermediate-files/bath_agg_min.grd", 
            overwrite = TRUE, format = "raster")

## stack and write out to file:
depths <- stack(depth_min, depth_max)
names(depths) <- c("depth_min", "depth_max")
writeRaster(depths, "data-processed/masks/elev-and-bath-layers/depths_minmax.grd", 
            overwrite = TRUE, format = "raster")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######     Create general depth restriction layers        #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## create boolean raster informing which areas are:
##      - have max that is less than -200m elev
##      - have min that is greater than 0 
##      - are defined as being on land in tid
depth <- depth_min
depth[depth_max < -200] <- NA
depth[tid == 0] <- NA
depth[depth_min > 0] <- NA
#plot(depth, asp = 1)

## looks good, save it and a mask version of it!
depth_mask <- depth
depth_mask[!is.na(depth_mask)] <- 1 ## change all areas where elev > -200m to 1
#plot(depth_mask, asp = 1)
writeRaster(depth, "data-processed/masks/elev-and-bath-layers/raster_200mdepth.grd", 
            overwrite = TRUE, format = "raster") ## with depth values 
writeRaster(depth_mask, "data-processed/masks/elev-and-bath-layers/raster_200mdepth_mask.grd", 
            overwrite = TRUE, format = "raster") ## mask 


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######   Rasterize marine & terrestrial seasonal temperatures #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
            res = 1)

## read in seasonal high and low temp data:
terr_seasonal_high <- read.csv("data-processed/intermediate-files/terrestrial_seasonal-max-temps_6mo-dormancy.csv")
terr_seasonal_low <- read.csv("data-processed/intermediate-files/terrestrial_seasonal-min-temps_6mo-dormancy.csv")
terr_mean <- read.csv("data-processed/intermediate-files/terrestrial_avg-temp.csv")

## rasterize:
raster_terr_high <- rasterize(terr_seasonal_high[, 1:2], r, terr_seasonal_high[,3:7], 
                              fun=base::mean) 
raster_terr_high[is.infinite(raster_terr_high)] <- NA
names(raster_terr_high) <- c("seasonal_high_temp", "hot_dormancy_6mo", "cold_dormancy_6mo", 
                             "hot_acc_temp", "hot_acc_temp_dormancy")
##plot(raster_terr_high[[1]], asp = 1)

raster_terr_low <- rasterize(terr_seasonal_low[, 1:2], r, terr_seasonal_low[,3:7], fun=base::mean)
raster_terr_low[is.infinite(raster_terr_low)] <- NA
names(raster_terr_low) <-  c("seasonal_low_temp", "cold_dormancy_6mo", "hot_dormancy_6mo", 
                             "cold_acc_temp", "cold_acc_temp_dormancy")
##plot(raster_terr_low[[1]], asp = 1)

raster_terr_mean <- rasterize(terr_mean[, 1:2], r, terr_mean[,3], fun=base::mean)
raster_terr_mean[is.infinite(raster_terr_mean)] <- NA
names(raster_terr_mean) <-  c("mean_temp")
raster_terr_mean <- mask(raster_terr_mean, raster_terr_low[[1]])
##plot(raster_terr_mean, asp = 1)

## restrict by zoogeographical realms:
cmec <- read_sf("data-raw/polygons/CMEC regions & realms/newRealms.shp")

# simplify the object to make it faster to use
cmec <- cmec %>% 
  st_simplify(dTolerance = 0.01) %>% 
  group_by(Realm) %>% 
  summarise()

## rasterize using Adam's function:
rast <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
               crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
               res = 1, vals = 1)

source("R/shp2rast.R") 
for (i in 1:nrow(cmec)) {
  cur <- cmec[i,]
  shp <- as_Spatial(cur)
  
  if (i==1) {
    zoo_rast <- shp2rast(shp = shp, rast = rast) 
    zoo_rast[!is.na(zoo_rast)] = i
  }
  else {
    temp <- shp2rast(shp = shp, rast = rast)
    temp[!is.na(temp)] = i
    zoo_rast <- addLayer(zoo_rast, temp)
  }
}

mosaic <- zoo_rast[[1]]
for (i in 1:nlayers(zoo_rast)-1) {
  mosaic <- mosaic(mosaic, zoo_rast[[i+1]], fun = base::min)
}
plot(mosaic)

## get rid of cells in realm that there is no temperature data for:
#plot(!is.na(mosaic) & is.na(raster_terr_high[[1]])) # get rid of these 

mosaic <- mask(raster_terr_high[[1]], mosaic)

raster_terr_high_zoo <- mask(raster_terr_high, mosaic)
raster_terr_low_zoo <- mask(raster_terr_low, mosaic)
raster_terr_mean_zoo <- mask(raster_terr_mean, mosaic)

## write out:
writeRaster(raster_terr_low_zoo, "data-processed/temperature-data/terrestrial/raster_terr_low.grd", 
            overwrite = TRUE, format = "raster")
writeRaster(raster_terr_high_zoo, "data-processed/temperature-data/terrestrial/raster_terr_high.grd", 
            overwrite = TRUE, format = "raster")
writeRaster(raster_terr_mean_zoo, "data-processed/temperature-data/terrestrial/raster_terr_mean.grd", 
            overwrite = TRUE, format = "raster")


## write out mask layer for use in restricting realized ranges:
raster_terr_mask <- raster_terr_high_zoo[[1]]
raster_terr_mask[!is.na(raster_terr_mask)] = 1
##plot(raster_terr_mask, asp = 1)
writeRaster(raster_terr_mask, "data-processed/masks/raster_terr_mask.grd", 
            overwrite = TRUE, format = "raster")

## write out new version of temp data, excluding cells outside of realm data:
terr_seasonal_high_zoo <- as.data.frame(rasterToPoints(raster_terr_high_zoo))
terr_seasonal_low_zoo <- as.data.frame(rasterToPoints(raster_terr_low_zoo))
names(terr_seasonal_high_zoo)[1:2] <- names(terr_seasonal_low_zoo)[1:2] <- c("longitude", "latitude")
write.csv(terr_seasonal_high_zoo, 
          "data-processed/intermediate-files/terrestrial_seasonal-max-temps_6mo-dormancy_zoo.csv", 
          row.names = FALSE)
write.csv(terr_seasonal_low_zoo, 
          "data-processed/intermediate-files/terrestrial_seasonal-min-temps_6mo-dormancy_zoo.csv", 
          row.names = FALSE)

## read in seasonal high, low and mean temp data:
marine_seasonal_high <- read.csv("data-processed/intermediate-files/marine_seasonal-max-temps_6mo.csv") 
marine_seasonal_low <- read.csv("data-processed/intermediate-files/marine_seasonal-min-temps_6mo.csv")

## rasterize:
raster_marine_high <- rasterize(marine_seasonal_high[, 1:2], r, marine_seasonal_high[,3:7], fun=base::mean)
raster_marine_high[is.infinite(raster_marine_high)] <- NA
names(raster_marine_high) <-  names(raster_terr_high) 
##plot(raster_marine_high[[1]], asp = 1)

raster_marine_low <- rasterize(marine_seasonal_low[, 1:2], r, marine_seasonal_low[,3:7], fun=base::mean)
raster_marine_low[is.infinite(raster_marine_low)] <- NA
names(raster_marine_low) <- names(raster_terr_low) 
##plot(raster_marine_low[[1]], asp = 1)

## write out mask layer for use in restricting realized ranges:
raster_marine_mask <- raster_marine_low
raster_marine_mask[!is.na(raster_marine_mask)] = 1
##plot(raster_marine_mask[[1]], asp = 1)
writeRaster(raster_marine_mask[[1]], "data-processed/masks/raster_marine_mask.grd", 
            overwrite = TRUE, format = "raster")

## write out:
writeRaster(raster_marine_low, "data-processed/temperature-data/marine/raster_marine_low.grd", 
          overwrite = TRUE, format = "raster")
writeRaster(raster_marine_high, "data-processed/temperature-data/marine/raster_marine_high.grd", 
         overwrite = TRUE, format = "raster")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######   Create intertidal seasonal temperature rasters       #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## restrict intertidal species to continental shelf around land masses:
bath <- raster("data-processed/masks/elev-and-bath-layers/raster_200mdepth_mask.grd")

## clump and turn into polygons 
cs_clumps <- bath %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(cs_clumps)

## buffer land by 1 cell:
t_poly <- raster_terr_high[[1]]
t_poly[!is.na(t_poly)] <- 1
t_poly <- t_poly %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
smooth <- smoothr::smooth(t_poly, method = "ksmooth", smoothness = 10) 
st_crs(smooth) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
## plot(st_geometry(smooth), axes = TRUE) 

buffer_land <- st_buffer(smooth, dist = 1)
plot(buffer_land)

## filter out clumps that don't intersect with the terrestrial mask:
intersects <- st_intersects(cs_clumps, buffer_land, sparse = FALSE)[,]
intersects <- filter(cs_clumps, intersects == TRUE)
plot(intersects)

## get temps within these areas:
raster_marine_low <- stack("data-processed/temperature-data/marine/raster_marine_low.grd")
raster_marine_high <- stack("data-processed/temperature-data/marine/raster_marine_high.grd")
raster_terr_low <- stack("data-processed/temperature-data/terrestrial/raster_terr_low.grd")
raster_terr_high <- stack("data-processed/temperature-data/terrestrial/raster_terr_high.grd")
raster_terr_mean <- stack("data-processed/temperature-data/terrestrial/raster_terr_mean.grd")

# merge marine and terr temps, giving priority to SST in cells with both 
raster_intertidal_low <- merge(raster_marine_low[[1:5]], raster_terr_low[[1:5]]) 
names(raster_intertidal_low) <- names(raster_marine_low)[1:5]
raster_intertidal_high <- merge(raster_marine_high[[1:5]], raster_terr_high[[1:5]])
names(raster_intertidal_high) <- names(raster_marine_high)[1:5]
raster_intertidal_mean <- merge(raster_marine_mean[[1]], raster_terr_mean[[1]])
names(raster_intertidal_mean) <- names(raster_marine_mean)[1]

# mask by intertidal coninental shelf areas:
raster_intertidal_low <- mask(raster_intertidal_low, intersects)
raster_intertidal_high <- mask(raster_intertidal_high, intersects)
plot(raster_intertidal_low)
plot(raster_intertidal_high)

## write out:
writeRaster(raster_intertidal_low, "data-processed/temperature-data/intertidal/raster_intertidal_low.grd", 
            overwrite = TRUE, format = "raster")
writeRaster(raster_intertidal_high, "data-processed/temperature-data/intertidal/raster_intertidal_high.grd", 
            overwrite = TRUE, format = "raster")

## write out mask layer for use in restricting realized ranges:
raster_intertidal_mask <- raster_intertidal_low
raster_intertidal_mask[!is.na(raster_intertidal_mask)] = 1
##plot(raster_intertidal_mask[[1]], asp = 1)
writeRaster(raster_intertidal_mask[[1]], "data-processed/masks/raster_intertidal_mask.grd", 
            overwrite = TRUE, format = "raster")
