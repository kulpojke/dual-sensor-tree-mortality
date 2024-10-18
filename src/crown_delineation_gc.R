#!/usr/bin/env Rscript

# parse args
args = commandArgs(trailingOnly=TRUE)
laz <- args[1]
tif <- args[2]
utm_zone <- args[3]

package.install(c('plyr', 'dplyr', 'terra', 'sf', 'liddR', 'future'))

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

round_to_even <- function(n) {
  #' rounds to even number
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


generate_custom_id <- function(utm, geometry) {
  #' generates custom tree ID based on UTM zone and coordinate
  coords <- st_coordinates(geometry)
  x <- coords[, 'X']
  y <- coords[, 'Y']
  # Apply rounding to the nearest even integer
  x_even <- round_to_even(x)
  y_even <- round_to_even(y)
  return(paste(utm, x_even, y_even, sep = '_'))
}


# ------------------------------------------------------------
# paths that stay the same
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
  stage <- 'making paths'
  base_name <- tools::file_path_sans_ext(basename(laz))
  log_file <- paste0(base_name, '_log.txt')

# read las, drop -Z points bc this is HAG
  print('Readling laz...')
  stage <- 'readLAS'
  las <- readLAS(
    laz,
    filter='-drop_withheld -drop_z_below 0',
    select='xyzirc'
    )
  
  print('filtering duplicate points...')
  stage <- 'filter_duplicates'
  las <- filter_duplicates(las)
  
  # send las_check to logfile
  print('checkinging laz...')
  stage <- 'las_check'
  sink(log_file)
  las_check(las)
  sink()

  # laz is already normalized, drop unwanted by Z (Z=HAG)
  stage <- 'filter_poi'
  normed <- filter_poi(las, Z >= 3 & Z <= 120)

  # read chm
  stage <- 'rast(tif)'
  chm <- rast(tif)
  
  # find ttops with height dependent window size
  wf <- function(x) {
    y <- abs(x/8)
    y[x <= 32] <- 4
    y[x > 80] <- 10
    return(y)
  }
  stage <- 'locate_trees(chm, lmf(wf))'
  ttops <- locate_trees(chm, lmf(wf))

  #  replace treeID with unique treeID
  stage <- 'generate_custom_id'
  ttops$uniqueID <- mapply(
    generate_custom_id,
    utm = utm_zone,
    geometry = ttops$geometry
    )
  
  # segment
  stage <- 'segment_tree'
  segs = segment_trees(
    normed,
    dalponte2016(chm, ttops)
  )
  
  stage <- 'crown_metric'
  crowns <- crown_metrics(
    segs,
    func=.stdmetrics,
    geom='concave'
  )
  
  stage <- 'merge'
  crowns <- merge(
    x=crowns,
    y=st_drop_geometry(ttops)[, c('treeID', 'uniqueID')],
    by='treeID',
    all.x=TRUE
    )
  
  # filter out trees shorter than 5m
  stage <- 'filte'
  crowns <- filter(crowns, zq95 >= 5)
  
  # write vectors
  stage <- 'st_write(ttops)'
  st_write(
    ttops,
    file.path(
      ttops_dir,
      paste0(base_name, '_ttops.gpkg')
    )
  )
  
  # write crowns
  stage <- 'st_write(crowns)'
  st_write(
    crowns,
    file.path(
      crowns_dir,
      paste0(base_name, '_crowns.gpkg')
    )
  )

}, error = function(e) {
  sink(log_file)
  print(paste('Tile', laz, 'is bad! During', stage, 'it was all like:'))
  print('- - - - - - - - - - - - - - - - - - - - - - -')
  print('')
  print(e$message)
  print('')
  print('- - - - - - - - - - - - - - - - - - - - - - -')
  
})

