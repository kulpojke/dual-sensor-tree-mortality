#!/usr/bin/env Rscript

# parse args
args = commandArgs(trailingOnly=TRUE)
site_dir  <- args[1] 
laz_dir   <- args[2]
cold_storage <- args[3]
utm_zone <- args[4]

print(paste('cold storage:', cold_storage))
print(paste('utm zone:', utm_zone))

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

# to move file from one drive to another
move_file_across_drives <- function(src, dst) {
  # copy the file to the target location
  copy_success <- file.copy(src, dst)
  
  # if the copy was successful, delete the original file
  if (copy_success) {
    remove_success <- file.remove(src)
    
    if (remove_success) {
      cat("File moved successfully.\n")
      return(TRUE)
    } else {
      cat("File copied but original could not be deleted.\n")
      return(FALSE)
    }
  } else {
    cat("File could not be copied.\n")
    return(FALSE)
  }
}

# ------------------------------------------------------------
# paths that stay the same
chm_dir = paste0(site_dir, '/', 'chm')
if (!dir.exists(chm_dir)){
  dir.create(chm_dir)
}

ttops_dir = paste0(site_dir, '/', 'ttops')
if (!dir.exists(ttops_dir)){
  dir.create(ttops_dir)
}

crowns_dir = paste0(site_dir, '/', 'crowns')
if (!dir.exists(crowns_dir)){
  dir.create(crowns_dir)
}

segs_dir = paste0(site_dir, '/', 'segmented_pc')
if (!dir.exists(segs_dir)){
  dir.create(segs_dir)
}

# list tiles
pattern <- '\\.la[sz]$'
tiles <- list.files(path=laz_dir, pattern=pattern, full.names=TRUE)

# run on each tile
for (f in tiles) {
  # get basename
  base_name <- tools::file_path_sans_ext(basename(f))

  print('Reading points')
  las <- readLAS(
    f,
    filter='-drop_withheld',
    select='xyzirc'
    )
  
  las <- filter_duplicates(las)
  las_check(las)

  # normalize and drop unwanted by height
  print('Normalizing')
  normed <- normalize_height(las, tin())
  normed <- filter_poi(normed, Z >= 3 & Z <= 120)

  # make chm
  print('creating chm')
  chm <- rasterize_canopy(
    normed,
    res=0.5,
    pitfree(
      thresholds=c(0, 10, 20),
      max_edge=c(0, 1.5))
    )
  
  #  smooth
  print('smooting chm...')
  w <- matrix(1, 3, 3)
  smoothed <- terra::focal(chm, w, fun=mean, na.rm=TRUE)
  
  smooth_chm_path <- paste0(site_dir, '/chm', '/', base_name, '.tif')
  writeRaster(smoothed, smooth_chm_path, overwrite=TRUE)
  
  # find ttops with height dependent window size
  print( 'finding ttops...')
  
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
  print('segmenting...')
  segs = segment_trees(
    normed,
    dalponte2016(smoothed, ttops)
  )
  
  print('crown metrics...')
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
  
  #write vectors
  print('writing ttops')
  st_write(
    ttops,
    file.path(
      site_dir,
      'ttops',
      paste0(base_name, '.gpkg')
    )
  )
  
  print('writing crowns')
  st_write(
    crowns,
    file.path(
      site_dir,
      'crowns',
      paste0(base_name, '.gpkg')
    )
  )

  print('writing segmented point cloud to cold storage')
  print(file.path(
    cold_storage,
    'retiled_lidar',
    paste0(base_name, '.laz')
  ))
  
  writeLAS(
    segs,
    file.path(
      cold_storage,
      'retiled_lidar',
      paste0(base_name, '.laz')
      )
  )

  print('moving lidar tile to cold storage')
  
  dst <- file.path(
    cold_storage,
    'segmented_lidar',
    paste0(basename(f))
  )
  
  move_file_across_drives(f, dst)
}
