load(file=".RData")

setwd("~/PM")

install.packages(pkgs= "sf", type = "source", configure.args = c("--with-gdal-config=/bin/gdal-config", "--with-proj-lib=/usr/local/lib64", "--with-proj-include=/usr/local/include", "--with-geos-config=/bin/geos-config"))
install.packages("Rcpp")
install.packages(pkgs= "terra", type = "source", configure.args = c("--with-gdal-config=/bin/gdal-config"))
install.packages("raster")
# install.packages(pkgs= "gstat", type = "source", configure.args = c("--with-gdal-config=/bin/gdal-config", "--with-proj-lib=/usr/local/lib64", "--with-proj-include=/usr/local/include", "--with-geos-config=/bin/geos-config"))
# install.packages(pkgs= "sftime", type = "source", configure.args = c("--with-gdal-config=/bin/gdal-config", "--with-proj-lib=/usr/local/lib64", "--with-proj-include=/usr/local/include", "--with-geos-config=/bin/geos-config"))
# install.packages(pkgs= "stars", type = "source", configure.args = c("--with-gdal-config=/bin/gdal-config", "--with-proj-lib=/usr/local/lib64", "--with-proj-include=/usr/local/include", "--with-geos-config=/bin/geos-config"))
library(sf)
library(geojsonsf)
library(readr)
library(dplyr)
library(stars)
library(ggplot2)
library(qmap)
library(raster)
library(perm)
library(sp)
library(spacetime)
library(gstat)
library(lattice)


## Data Preparation
data <- geojson_sf("data/ms_Finedust_PM25.geojson")

data$timestamp <- parse_datetime(data$timestamp)

library(ggplot2)

tmp_muenster <- tempfile(fileext = ".geojson")
download.file(
  "https://opendata.stadt-muenster.de/sites/default/files/stadtbezirke-muenster.geojson",
  tmp_muenster
)
münster_sf <- read_sf(tmp_muenster)
münster_borders <- st_as_sf(münster_sf$geometry, crs=4326)
ms_bbox <- st_bbox(münster_borders)

data_oct_feb = filter(data, 
  (data$timestamp > as.POSIXct("2024-10-01 00:00:00", tz="UTC") 
  & data$timestamp < as.POSIXct("2025-02-01 00:00:00", tz="UTC")))

data_nov_feb = filter(data, 
  (data$timestamp > as.POSIXct("2024-11-01 00:00:00", tz="UTC") 
  & data$timestamp < as.POSIXct("2025-02-01 00:00:00", tz="UTC")))


# get stationary data and crop by münster bbox
data_stationary <- read.csv("data/stationary/stationary_pm25.csv",
                            colClasses=c("numeric","numeric","character","character","character","numeric","character")) |>
  st_as_sf(coords=c("x_col_name"="lon","y_col_name"="lat"), crs=4326)
data_stationary <- data_stationary[
    st_coordinates(data_stationary)[, 1] >= ms_bbox$xmin & st_coordinates(data_stationary)[, 1] <= ms_bbox$xmax &
    st_coordinates(data_stationary)[, 2] >= ms_bbox$ymin & st_coordinates(data_stationary)[, 2] <= ms_bbox$ymax, 
]
data_stationary$timestamp <- as.POSIXct(data_stationary$createdAt, tryFormats = c("%Y-%m-%dT%H:%M:%OS"))
data_stationary <- data_stationary[data_stationary$value <= 100,]
# maybe remove outliers over 100

borders_plot <- ggplot() +
  geom_sf(data = münster_borders, fill = NA) + theme_minimal() + 
  geom_sf(data = data_stationary, aes(color = value), size = 1) +
  scale_color_continuous(limits = c(0, 75),   # Limit color scale
                         low = "green", high = "red") +
  coord_sf(xlim = c(ms_bbox["xmin"], ms_bbox["xmax"]),   
           ylim = c(ms_bbox["ymin"], ms_bbox["ymax"]), expand = FALSE) 
print(borders_plot)

borders_plot <- ggplot() +
  geom_sf(data = münster_borders, fill = NA) + theme_minimal() + 
  geom_sf(data = data_oct_feb, aes(color = value), size = 1) +
  scale_color_continuous(limits = c(0, 75),   # Limit color scale
                        low = "green", high = "red") +
  coord_sf(xlim = c(ms_bbox["xmin"], ms_bbox["xmax"]),   
           ylim = c(ms_bbox["ymin"], ms_bbox["ymax"]), expand = FALSE) 
print(borders_plot)

# final_plot <- ggplot() +
#   geom_sf(data = münster_sf, fill = NA, color = "black") +
#   geom_sf(data = data$value, aes(fill = kde_value), alpha = 1) +
#   scale_fill_viridis_c(option = "D", name = "Density") +
#   theme_minimal()
#
#
# ggplot(data = data) +
#   geom_sf(aes(color = value), size = 1) +          # Use geom_sf for sf objects
#   labs(title = "Geographic Plot of Measurement Points",
#        x = "Longitude",
#        y = "Latitude") +
#   theme_minimal()      
#
ggplot(data = data_nov_feb, aes(x = timestamp, y = value, color=box_id)) +
  geom_point() +
  labs(title = "Measurements over Time", x = "Timestamp", y = "Measurement") +
  theme_minimal() +
  scale_color_discrete(name ="box_id") + 
  ylim(0,100) +
  #scale_y_log10() + 
  theme(legend.position="none")

data_stationary_plot = data_stationary[seq(1, nrow(data_stationary), 100), ]
ggplot(data = data_stationary_plot, aes(x = timestamp, y = value, color=boxName)) +
  geom_point() +
  labs(title = "Measurements over Time", x = "Timestamp", y = "Measurement") +
  theme_minimal() +
  scale_color_discrete(name ="box_id") + 
  ylim(0,100)
  #scale_y_log10()+
  #theme(legend.position="none")
remove(data_stationary_plot)

################################################################################
# COMPARE MOBILE TO STATIONARY DATA AND NORMALIZE TO SOME INDEX
#  this is because the sensors are not necessarily mounted in accordance to 
#  regulation, so using measured values could cause some confusion.

# how would we do that?
# when comparing the histograms of both datasets, the mobile sensor data seems to be
# slightly skewed towards higher values. so normalization could be useful...
hist(data_stationary$value, breaks=500, xlim=c(0,100))
hist(data_nov_feb$value, breaks=500, xlim=c(0,100))
summary(data_stationary$value)
summary(data_nov_feb[data_nov_feb$value <= 100,]$value)

# use quantile mapping to account for the bias of the bike sensors being mounted
#  at a lower elevation and being operated near roads
# install.packages("qmap")
library(qmap)

qm_fit <- fitQmapQUANT(data_stationary$value, data_nov_feb$value)
data_nov_feb$pm_idx <- doQmap(data_nov_feb$value, qm_fit)

# we call this the aligned PM index now because it doesnt reflect measurement values anymore
summary(data_nov_feb$pm_idx)


## PLANNING
# what to consider as input data
#   we have the bike sensor data. that certainly
#   compile into 4 seasonal bins -- not possible
#   normalize to 125x125m grid
#   if that works, also consider adding NDVI and road network

################################################################################
# CREATE GRID FOR DATA ACCUMULATION
#install.packages("sf")
#install.packages("raster")
library(sf)
library(raster)

#city_borders <- st_read("path_to_your_shapefile.shp")

#combine stationary and mobile datasets
data_nov_feb$stationary <- rep(FALSE, nrow(data_nov_feb))
data_stationary$stationary <- rep(TRUE, nrow(data_stationary))

combined_data <- data.frame(
  boxName = append(data_stationary$boxName, data_nov_feb$box_id),
  boxId = append(data_stationary$boxId, data_nov_feb$box_id),
  value = append(data_stationary$value, data_nov_feb$pm_idx),
  geometry = append(data_stationary$geometry, data_nov_feb$geometry),
  timestamp = append(data_stationary$timestamp, data_nov_feb$timestamp),
  stationary = append(data_stationary$stationary, data_nov_feb$stationary)
)
combined_data$unit=rep(data_stationary$unit[1], nrow(combined_data))
combined_data <- st_as_sf(combined_data)

# Create a grid of n x n meters
grid_size <- 500  # in meters
# Get the bounding box of the shapefile

# use epsg3857
ms_bbox <- st_transform(ms_bbox, 3857)
münster_borders <- st_transform(münster_borders, 3857)
münster_borders <- st_union(münster_borders)
combined_data <- st_transform(combined_data, 3857)

raster_template <- raster(
  xmn = ms_bbox["xmin"], xmx = ms_bbox["xmax"],
  ymn = ms_bbox["ymin"], ymx = ms_bbox["ymax"],
  res = c(grid_size, grid_size), crs = 3857
)
# Convert the raster to a polygon grid
grid <- rasterToPolygons(raster_template)
# Transform `grid` to an sf object
grid_sf <- st_as_sf(grid)
# Clip the grid by the city borders
#grid_clipped <- st_intersection(grid_sf, münster_borders)
grid_clipped <- grid_sf

#segment data into days. and calculate mean grid cells for each
combined_data$date <- trunc(combined_data$timestamp, "day") |> as.character()
dates <- unique(combined_data$date)
dates <- dates[order(as.Date(unique(combined_data$date), format="%Y-%m-%d"))]

# just using the entire day mean for stationary sensors will result in wildly different
# values to the cell means produced by mobile sensors, since there will at best
# be a minute of those within a 500x500m cell. There are two possible approaches now
# 1. increase spatial resolution
# 2. omit stationary data for periods where no mobile measurements were made.
# going with the latter for now as vgm calculation already takes very long on 
# a 1-day-resolution.

#add daytime column to combined data
combined_data$daytime <- format(combined_data$timestamp, "%H:%M:%S")
combined_data$daytime <- as.POSIXct(format(combined_data$timestamp, "%H:%M:%S"),
                                    format="%H:%M:%S") #date is irrelevant here

# view a histogram of where most measurements happen
hist(combined_data$daytime, breaks=100,
     xlab = deparse1(substitute(combined_data$daytime)),
     plot = TRUE, freq = FALSE,
     start.on.monday = TRUE, format, right = TRUE)
hist(combined_data[combined_data$stationary==TRUE,]$daytime, breaks=100,
     xlab = deparse1(substitute(combined_data$daytime)),
     plot = TRUE, freq = FALSE,
     start.on.monday = TRUE, format, right = TRUE)
hist(combined_data[combined_data$stationary==FALSE,]$daytime, breaks=100,
     xlab = deparse1(substitute(combined_data$daytime)),
     plot = TRUE, freq = FALSE,
     start.on.monday = TRUE, format, right = TRUE)

#filter for data-points in peak bicycle hours.
# we take 07:00-09:30 & 15:00-17:00
combined_data_peak <-
  combined_data[(combined_data$daytime >= as.POSIXct("07:00", f="%H:%M") & combined_data$daytime <= as.POSIXct("09:30", f="%H:%M")) 
                | combined_data$daytime >= as.POSIXct("15:00", f="%H:%M") & combined_data$daytime <= as.POSIXct("17:00", f="%H:%M"),]

cell_means <- list()
cell_mean_centroids <- list()

# for (day in dates){
#   day_data <- combined_data_peak[combined_data_peak$date == day,]
#   #day_data <- combined_data[combined_data$date == day,]
#   measurement_grid <- st_join(grid_clipped, day_data, join = st_contains)
# 
#   stationary_exists <- (length(day_data$value[day_data$stationary == TRUE]) > 0)
#   mobile_exists <- (length(day_data$value[day_data$stationary == FALSE]) > 0)
# 
#   #visualise histograms for each
#   par(mfrow = c(3, 1))
# 
#   # hist(day_data$value,
#   #      breaks=100, xlim=c(0,100),
#   #      main = paste("all ", day),
#   #      xlab = ("measured and adjusted µg/m³")
#   # )
#   #
#   # if(stationary_exists){
#   #   hist(day_data$value[day_data$stationary == TRUE],
#   #        breaks=100, xlim=c(0,100),
#   #        main = paste("stationary ", day),
#   #        xlab = ("measured µg/m³")
#   #        )
#   # }
#   #
#   # if(mobile_exists){
#   #   hist(day_data$value[day_data$stationary == FALSE],
#   #        breaks=100, xlim=c(0,100),
#   #        main = paste("bike ", day),
#   #        xlab = ("adjusted µg/m³")
#   #   )
#   # }
# 
#   # Calculate daily mean by grid cell
#   measurement_grid <- measurement_grid %>%
#     filter(!is.na(value))
# 
#   grid_with_values <- measurement_grid %>%
#     group_by(geometry) %>%
#     summarize(mean_value = mean(value, na.rm = TRUE)) %>%
#     ungroup()
#   grid_with_values$date <- day
#   grid_with_values$timeSlice <- as.POSIXct(grid_with_values$date)
# 
#   cell_means[[length(cell_means)+1]] <- grid_with_values
# 
#   centroids_with_means <- st_centroid(grid_with_values)
# 
#   cell_mean_centroids[[length(cell_mean_centroids)+1]] <- centroids_with_means
# }

for (day in dates){
  #split into two groups for peak hours
  day_data <- combined_data_peak[combined_data_peak$date == day,]
  #day_data <- combined_data[combined_data$date == day,]

  morning_data <- day_data[(day_data$daytime >= as.POSIXct("07:00", f="%H:%M")
                        & day_data$daytime <= as.POSIXct("09:30", f="%H:%M")),]
  afternoon_data <- day_data[(day_data$daytime >= as.POSIXct("15:00", f="%H:%M")
                        & day_data$daytime <= as.POSIXct("17:00", f="%H:%M")),]

  # Calculate peak mean by grid cell
  m_grid_morn <- st_join(grid_clipped, morning_data, join = st_contains) |>
    filter(!is.na(value))
  m_grid_aftn <- st_join(grid_clipped, afternoon_data, join = st_contains) |>
    filter(!is.na(value))

  grid_with_values_morn <- m_grid_morn |>
    group_by(geometry) |>
    summarize(mean_value = mean(value, na.rm = TRUE)) |>
    ungroup()
  grid_with_values_aftn <- m_grid_aftn |>
    group_by(geometry) |>
    summarize(mean_value = mean(value, na.rm = TRUE)) |>
    ungroup()

  grid_with_values_morn$date <- day
  grid_with_values_morn$timeSlice <- as.POSIXct(
    paste(grid_with_values_morn$date, "08:15"), f="%Y-%m-%d %H:%M" ) #midpoint
  grid_with_values_aftn$date <- day
  grid_with_values_aftn$timeSlice <- as.POSIXct(
    paste(grid_with_values_aftn$date, "16:00"), f="%Y-%m-%d %H:%M") #midpoint

  cell_means[[length(cell_means)+1]] <- grid_with_values_morn
  cell_means[[length(cell_means)+1]] <- grid_with_values_aftn

  centr_with_means_morn <- st_centroid(grid_with_values_morn)
  centr_with_means_aftn <- st_centroid(grid_with_values_aftn)

  cell_mean_centroids[[length(cell_mean_centroids)+1]] <- centr_with_means_morn
  cell_mean_centroids[[length(cell_mean_centroids)+1]] <- centr_with_means_aftn
}

par(mfrow = c(1, 1))

################################################################################
# EXAMINE FOR WEEKDAY-WEEKEND DIFFERENCES
#  is the difference significant? does a periodicity exist?
#  if so, the data can be divided in ways to let these 

# add a weekday column
combined_data$weekday <- weekdays(as.Date(combined_data$date), abbreviate=TRUE)

par(mfrow = c(2, 1))
hist(combined_data[!(combined_data$weekday == "Sat" | combined_data$weekday == "Sun"),]$value, 
     breaks=500, xlim=c(0,10),
     main = paste("Weekdays"),
     xlab = ("adjusted µg/m³")
)
hist(combined_data[combined_data$weekday == "Sat" | combined_data$weekday == "Sun",]$value, 
     breaks=500, xlim=c(0,10),
     main = paste("Weekends"),
     xlab = ("adjusted µg/m³")
)
# #exploring the difference
summary(combined_data[!(combined_data$weekday == "Sat" | combined_data$weekday == "Sun"),]$value)
summary(combined_data[(combined_data$weekday == "Sat" | combined_data$weekday == "Sun"),]$value)

# run a permutation test to see if the difference is significant
#install.packages("perm")
library(perm)
# TODO: reevaluate this test. i think n may be poorly chosen
set.seed(123)

sampled_weekday <- sample(combined_data[!(combined_data$weekday == "Sat" | combined_data$weekday == "Sun"),]$value, 30000)
sampled_weekend <- sample(combined_data[(combined_data$weekday == "Sat" | combined_data$weekday == "Sun"),]$value, 10000)
summary(sampled_weekday)
summary(sampled_weekend)
sampled_combined <- c(sampled_weekday, sampled_weekend)
observed_diff <- mean(sampled_weekend) - mean(sampled_weekday)

# input data for permutation test

permTS(
  #sampled_weekend, sampled_weekday,
  combined_data[!(combined_data$weekday == "Sat" | combined_data$weekday == "Sun"),]$value,
  combined_data[(combined_data$weekday == "Sat" | combined_data$weekday == "Sun"),]$value,
  method = "exact.mc",
  mcsamples = 100000
  )

n_weekday <- length(sampled_weekday)
n_weekend <- length(sampled_weekend)
n <- 40000


n_permutations <- 10000
# Initialize a vector to store the permuted statistics
permuted_diff <- numeric(n_permutations)

pb <- txtProgressBar(min=0, max=n_permutations, initial = 0)

# Permutation test loop
# this takes a little while
for (i in 1:n_permutations) {
  setTxtProgressBar(pb,i)
  # Shuffle the combined dataset
  shuffled_data <- sample(sampled_combined, n)
  
  # Split into two new groups
  permuted_weekday <- shuffled_data[1:n_weekday]
  permuted_weekend <- shuffled_data[(n_weekday + 1):n]
  
  # Calculate the difference in means
  permuted_diff[i] <- mean(permuted_weekend, na.rm=TRUE) - mean(permuted_weekday, na.rm=TRUE)
}
close(pb)
weekend_weekday_p_value <- mean(abs(permuted_diff) >= abs(observed_diff))
# Display results
cat("Observed difference:", observed_diff, "\n")
cat("P-value:", weekend_weekday_p_value, "\n")

# p-value close to 0
# -> statistically significant difference between weekends and weekdays.
# the observed difference between weekend and weekday means is around 1.01
# looking back at the weekend-weekday histogram comparison, weekdays have a lot
# more collected data in the lower range. between 0 and 1 than weekends.
# This is counter-intuitive because an assumption would be that higher traffic
# on weekdays results in higher PM pollution.
# several factors could have played into this. Weather conditions, holidays,
# heating in winter, etc. 
# unfortunately it would pass beyond the scope of this seminar project to 
# investigate further.

#clean up
remove(observed_diff, n_weekday, n_weekend, shuffled_data, n, n_permutations, 
       permuted_diff, permuted_weekday, permuted_weekend, i, sampled_weekday,
       sampled_weekend, sampled_combined, pb)

################################################################################
# PREPARE THE SPATIOTEMPORAL INTERPOLATION  
#  I am working by the example of spatiotemporal interpolation for irregular
#  data as provided in the lecture. This part of the lecture was rather short.
library(sp)
library(spacetime)
library(gstat)
library(lattice)
# needs more than just a stars object

# we create an STIDF object. it is a class for unstructured spatiotemporal data.
# add POSIXct column for time slices we consider (like days or 6h intervals)

# create something usable from cell-means-centroid
c_means_centr <- do.call(rbind.data.frame, cell_mean_centroids)
c_means_centr <- st_transform(c_means_centr, 3857)

#subset the data to test out the kriging parameters without having to wait long.
c_means_centr <- c_means_centr[c_means_centr$timeSlice >= "2024-11-01" & c_means_centr$timeSlice <= "2025-01-01",]

#combined_data$timeSlice <- as.POSIXct(combined_data$date)
#combined_data_sp <- as(combined_data, "Spatial")

stidf = STIDF(
  as(c_means_centr$geometry,"Spatial"), 
  c_means_centr$timeSlice, 
  data.frame(z=c_means_centr$mean_value)
  )
stplot(stidf, main="plot of some PM data samples")

library(sftime)
sft = st_as_sftime(stidf) 

#create the 500m res spatial grid

par(mfrow = c(1, 1))
#the originl grid the data was accumulated on
plot(st_geometry(grid_sf), col = NA, border = 'blue')
plot(st_geometry(münster_borders), add = TRUE, col = NA, border = 'red')

grd <- SpatialGrid(
  points2grid(as(st_centroid(grid_sf), "Spatial"))
  ) |> as("SpatialPixels")

plot(grd, add = TRUE, col="green")
#fits!

#grd <- st_set_crs(grd, st_crs(c_means_centr))

#because we currently have regular time slices, we can work with this
# but definitely change it up once you use irregular ones
tgrd <- c_means_centr$timeSlice
stf <- STF(grd, tgrd)

library(stars)
st <- st_as_stars(stf) |> st_set_crs(st_crs(sft))

timelag <- difftime(unique(c_means_centr$timeSlice)[1:length(unique(c_means_centr$timeSlice))],unique(c_means_centr$timeSlice)[1]) /60 /60

timelag <- difftime(unique(c_means_centr$timeSlice)[1:24],unique(c_means_centr$timeSlice)[1]) /60 /60
# Empirical variogram
# it takes too long to calculate one for the entire dataset at once
# I'm hoping to be able to sufficiently fit a variogram to an empirical vgm
# that was computed from subsets of the data
emp_variogram <- variogramST(z ~ 1, data = stidf, tunit = "hours",
                             tlags=timelag,
                             progress = TRUE,
                             na.omit = TRUE,
                             assumeRegular = FALSE
                             )

saveRDS(emp_variogram, "vgm_full_22tl.RDATA")
plot(emp_variogram, wireframe=F)

ggplot(data = c_means_centr, aes(x = timeSlice, y = mean_value)) +
  geom_point() +
  labs(title = "Measurements over Time", x = "Timestamp", y = "Measurement") +
  theme_minimal() +
  ylim(0,100)
#scale_y_log10()+
#theme(legend.position="none")

####
# FITTING A VGM
emp_variogram <- readRDS("./vgm_novdec_20tl.RDATA")
plot(emp_variogram)
plot(emp_variogram, wireframe=T)

#due to the irregularity of the data, a sum-metric or simplified sum-metric
# model would be smart.
# how do i select the right parameters? inform yourself by Gräler et al. 2016

# define a variogram model
# oriented by this guide: https://r-video-tutorial.blogspot.com/2015/08/spatio-temporal-kriging-in-r.html,
# and Gräler et al. 2016

# lower and upper bounds
pars.l <- c(sill.s = 10, range.s = 500, nugget.s = 0,
            sill.t = 20, range.t = 1, nugget.t = 0,
            sill.st = 0, range.st = 10, nugget.st = 0, 
            anis = 1000)
pars.u <- c(sill.s = 200, range.s = 5000, nugget.s = 100,
            sill.t = 100, range.t = 60, nugget.t = 100,
            sill.st = 100, range.st = 5000, nugget.st = 100,
            anis = 3000)

#metric model
metric <- vgmST("metric", joint = vgm(50,"Mat", 2000, 0), stAni=200) 
vgm_metric <- fit.StVariogram(emp_variogram, metric, method="L-BFGS-B",lower=pars.l)
plot(emp_variogram)
plot(emp_variogram, vgm_metric)
attr(metric_Vgm, "MSE")

#this is still largely guesstimation by eyeballing the sample variogram plot
#the following was taken off of 
sumMetricModel <- vgmST("sumMetric",
                        space = vgm(10, "Sph", 1000, 1),
                        time = vgm(24, "Exp", 30, 1),
                        joint = vgm(20, "Sph", 300, 2.5),
                        stAni = 1000, #m/h
)

vgm_sumMetric <- fit.StVariogram(emp_variogram, sumMetricModel, fit.method = 7,
                                 stAni = 1000, method = "L-BFGS-B",
                                 control = list(parscale = c(1, 100, 1, 1, 0.5, 1,
                                                             1, 100, 1, 100),
                                                maxit=10000),
                                 lower = pars.l,
                                 upper = pars.u,
                                 tunit="hours")
plot(emp_variogram)
plot(emp_variogram, vgm_sumMetric)
attr(vgm_sumMetric, "MSE")
#cant seem to get below MSE of 210

SimplesumMetric <- vgmST("simpleSumMetric",
                         space = vgm(5,"Sph", 500, 0),
                         time = vgm(500,"Sph", 500, 0), 
                         joint = vgm(1,"Sph", 500, 0), 
                         nugget=1, stAni=500)


vgm_simpleSumMetric <- fit.StVariogram(emp_variogram, SimplesumMetric,method = "L-BFGS-B",
                                       lower = c(sill.s = 10, range.s = 500, nugget.s = 1,
                                                 sill.t = 20, range.t = 20, nugget.t = 2,
                                                 sill.st = 0, range.st = 500, nugget.st = 1,
                                                 anis = 4000),)
plot(emp_variogram, vgm_simpleSumMetric)
attr(vgm_simpleSumMetric, "MSE")

#simpleSumMetric delivers lowest error

###########################################
# ST-KRIGING
#create the 500m res spatial grid

par(mfrow = c(1, 1))
#the originl grid the data was accumulated on
plot(st_geometry(grid_sf), col = NA, border = 'blue')
plot(st_geometry(münster_borders), add = TRUE, col = NA, border = 'red')

grd <- SpatialGrid(
  points2grid(as(st_centroid(grid_sf), "Spatial"))
) |> as("SpatialPixels")

plot(grd, add = TRUE, col="green")
#fits!

#grd <- st_set_crs(grd, st_crs(c_means_centr))

# this step is repeated until all time slices are interpolated
# can take several hours for each month
tgrd <- unique(c_means_centr$timeSlice)[21:60]
tgrd <- unique(c_means_centr$timeSlice)[61:90]
tgrd <- unique(c_means_centr$timeSlice)[91:122]
tgrd <- unique(c_means_centr$timeSlice)[123:152]
tgrd <- unique(c_means_centr$timeSlice)[152:184]
stf <- STF(grd, tgrd)

library(stars)
st <- st_as_stars(stf) |> st_set_crs(st_crs(sft))

krig_sft_0115_0131 <- krigeST(z~1, sft, st, modelList = vgm_simpleSumMetric, nmax=10, computeVar = T)
#stplot(locKrig_sft)
#mean_krig <- st_apply(locKrig_sft )
#plot(krig_sft_1201_1215)
#plot(st_geometry(münster_borders), add = TRUE, col = NA, border = 'red')
write_stars(krig_sft_0115_0131, "interpolation_0101_0115.tif")



# load all computed tiffs
# a bunch of them have been computed on the cluster
nov1_krig <- read_stars('interpolation_1101-1110.tiff')
nov2_krig <- read_stars('interpolation_1111-1130.tiff')
dec1_krig <- read_stars('interpolation_1201-1215.tiff')
dec2_krig <- read_stars('interpolation_1216_1231.tiff')
jan1_krig <- read_stars('interpolation_0101_0115.tiff')
jan2_krig <- read_stars('interpolation_0116_0131.tiff')

krig <- c(nov1_krig, nov2_krig, dec1_krig, dec2_krig, jan1_krig, jan2_krig, along=3)
plot(krig)
krig_mean <- st_apply(krig, MARGIN = c("x", "y"), FUN = mean, na.rm = TRUE)
plot(krig_mean)
write_stars(krig_mean, "interpolated_means.tif")

road_network <- read_sf("data/bicycleinfrastructure_MS.geojson") |> st_transform(3857)

road_raster_template <- raster(
  xmn = ms_bbox["xmin"], xmx = ms_bbox["xmax"],
  ymn = ms_bbox["ymin"], ymx = ms_bbox["ymax"],
  res = c(20, 20), crs = 3857
)
road_network <- road_network[st_geometry_type(road_network) %in% c("LINESTRING", "MULTILINESTRING"), ]
#road_raster <- rasterize(road_network, road_raster_template, fun = sum)

#cut road-raster with interpolated means
road_raster_updated <- st_as_stars(road_raster)
resampled_krig <- st_warp(krig_mean, dest = road_raster_updated, method = "bilinear", use_gdal = TRUE)
mask <- is.na(road_raster_updated[[1]])
resampled_krig_cut <- resampled_krig
resampled_krig_cut[[1]][mask] <- NA
plot(resampled_krig_cut)

# Save the raster to a file
write_stars(resampled_krig_cut, "road_raster.tif", format = "GTiff")
plot(st_transform(resampled_krig_cut, 4326))
write_stars(st_transform(resampled_krig_cut, 4326), "road_raster_4326.tif", format = "GTiff")
plot(road_raster)
