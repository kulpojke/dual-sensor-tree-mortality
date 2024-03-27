#!/usr/bin/env Rscript

# parse args
args = commandArgs(trailingOnly=TRUE)
site_dir  <- args[1] 
laz_dir   <- args[2] 

# libraries
library(plyr)
library(dplyr)
library(terra)
library(sf)
library(lidR)
library(future)
plan(multisession, workers=8)

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
  basename <- tools::file_path_sans_ext(basename(f))

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
  
  smooth_chm_path <- paste0(site_dir, '/chm', '/', basename, '.tif')
  writeRaster(smoothed, smooth_chm_path, overwrite=TRUE)
  
  # find ttops with height dependent window size
  print( 'finding ttops...')
  
  f <- function(x) {
    y <- abs(x/8)
    y[x <= 32] <- 4
    y[x > 80] <- 10
    return(y)
  }
  ttops <- locate_trees(smoothed, lmf(f))
  
  # segment
  print('segmenting...')
  segs = segment_trees(
    normed,
    dalponte2016(smoothed, ttops),
    attribute='IDdalponte'
  )
  
  print('crown metrics...')
  crowns <- crown_metrics(
    segs,
    func=.stdmetrics,
    attribute='IDdalponte',
    geom='concave'
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
      paste0(basename, '.gpkg')
    )
  )
  
  print('writing crowns')
  st_write(
    crowns,
    file.path(
      site_dir,
      'crowns',
      paste0(basename, '.gpkg')
    )
  )

   print('writing segmented point cloud')
  writeLAS(
    segs,
    file.path(
      segs_dir,
      paste0(basename, '.laz')
      )
  )
  
  
}
