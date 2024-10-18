#!/usr/bin/env Rscript

# parse args
args = commandArgs(trailingOnly=TRUE)
f  <- args[1]
utm_zone <- args[2]

# libraries
library(plyr)
library(dplyr)
library(terra)
library(sf)
library(lidR)
library(future)
plan(multisession, workers=8)

#---------------------------------------------------------------
# functions for later
# rounds to even number
round_to_even <- function(n) {
  rounded <- round(n)
  if (rounded %% 2 == 1) {  # if odd, make it even
    if (n - floor(n) < 0.5) {
      return(rounded - 1)  # round down
    } else {
      return(rounded + 1)  # or round up
    }
  } else {
    return(rounded)
  }
}


# function to generate custom tree ID based on UTM zone and coordinates
generate_custom_id <- function(utm, geometry) {
  coords <- st_coordinates(geometry)
  x <- coords[, "X"]
  y <- coords[, "Y"]
  # Apply rounding to the nearest even integer
  x_even <- round_to_even(x)
  y_even <- round_to_even(y)
  return(paste(utm, x_even, y_even, sep = "_"))
}


# ------------------------------------------------------------
# paths that stay the same
chm_dir = 'chm'
if (!dir.exists(chm_dir)){
  dir.create(chm_dir)
}

ttops_dir = 'ttops'
if (!dir.exists(ttops_dir)){
  dir.create(ttops_dir)
}

crowns_dir = 'crowns'
if (!dir.exists(crowns_dir)){
  dir.create(crowns_dir)
}


tryCatch({
  # get basename
  base_name <- tools::file_path_sans_ext(basename(f))
  log_file <- paste0(base_name, '_log.txt')

  print('Readling laz...')
  las <- readLAS(
    f,
    filter='-drop_withheld',
    select='xyzirc'
    )
  
  print('filtering laz...')
  las <- filter_duplicates(las)
  
  # send las_check to logfile
  print('checkinging laz...')
  sink(log_file)
  las_check(las)
  sink()

  # normalize and drop unwanted by height
  normed <- normalize_height(las, tin())
  normed <- filter_poi(normed, Z >= 3 & Z <= 120)

  # make chm
  chm <- rasterize_canopy(
    normed,
    res=0.5,
    pitfree(
      thresholds=c(0, 10, 20),
      max_edge=c(0, 1.5))
    )
  
  #  smooth
  w <- matrix(1, 3, 3)
  smoothed <- terra::focal(chm, w, fun=mean, na.rm=TRUE)
  
  smooth_chm_path <- paste0(site_dir, '/chm', '/', base_name, '.tif')
  writeRaster(smoothed, smooth_chm_path, overwrite=TRUE)
  
  # find ttops with height dependent window size
  wf <- function(x) {
    y <- abs(x/8)
    y[x <= 32] <- 4
    y[x > 80] <- 10
    return(y)
  }
  ttops <- locate_trees(smoothed, lmf(wf))

  #  replace treeID with unique treeID
  ttops$uniqueID <- mapply(generate_custom_id, utm = utm_zone, geometry = ttops$geometry)
  
  # segment
  segs = segment_trees(
    normed,
    dalponte2016(smoothed, ttops)
  )
  
  crowns <- crown_metrics(
  segs,
  func=.stdmetrics,
  geom='concave'
  )
  
  crowns <- merge(
    x=crowns,
    y=st_drop_geometry(ttops)[, c("treeID", "uniqueID")],
    by="treeID",
    all.x=TRUE
    )
  
  # filter out trees shorter than 5m
  crowns <- filter(crowns, zq95 >= 5)
  
  # write vectors
  st_write(
    ttops,
    file.path(
      ttops_dir,
      paste0(base_name, '_ttops.gpkg')
    )
  )
  
  # write crowns
  st_write(
    crowns,
    file.path(
      crowns_dir,
      paste0(base_name, '_crowns.gpkg')
    )
  )

}, error = function(e) {
  sink(log_file)
  print(paste('Tile', f, 'is bad! It was all like:'))
  print('- - - - - - - - - - - - - - - - - - - - - - -')
  print('')
  print(e$message)
  print('')
  print('- - - - - - - - - - - - - - - - - - - - - - -')
  
})

