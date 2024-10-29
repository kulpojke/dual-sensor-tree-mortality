#!/usr/bin/env Rscript

# libraries
library(fs)
library(plyr)
library(dplyr)
library(terra)
library(sf)
library(lidR)
library(future)
plan(multisession, workers=8)

# parse args
args = commandArgs(trailingOnly=TRUE)
laz_dir <- path(args[1])
chm_dir <- path(args[2])
ttops_dir <- path(args[3])
crowns_dir <- path(args[4])
segs_dir <- path(args[5])
log_dir <- path(args[6])
utm_zone <- path(args[7])

# make dirs if they DNE
dirs <- list(
  laz_dir,
  chm_dir,
  ttops_dir,
  crowns_dir,
  segs_dir,
  log_dir
)

for (d in dirs) {
  if (!dir.exists(d)){
    dir.create(d)
  }
}

# Initialize a global log variable
LOG <- character()

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


# function to append a message to the log
add_to_log <- function(message) {
  LOG <<- c(LOG, message)
}


# function to write contents of LOG to a file and clear it
write_to_log <- function(log_dir) {
  #  path to log file
  log_path <- file.path(log_dir, "R_log.txt")
  
  # open connection in append mode
  con <- file(log_path, open = "a")
  
  # write to file and close connection
  writeLines(LOG, con = con)
  close(con)
  
  # Clear the LOG
  LOG <<- character()
}


main <- function(complete) {
  # run on each tile
  for (f in dir_ls(laz_dir, glob = "*.laz")) {
    # get basename
    base_name <- tools::file_path_sans_ext(basename(f))

    # if file has already benn processed, skip it
    if (base_name %in% complete) {
      add_to_log('-------------------------------------------')
      add_to_log(paste(base_name, 'has already been processed. Skipping.'))
      add_to_log(format(Sys.time(), '%Y-%m-%d %H:%M:%S'))
      add_to_log('-------------------------------------------')
      write_to_log(log_dir)
      next

      } else {

      # make ttops, crowns, etc...
      tryCatch({
        
        add_to_log('-------------------------------------------')
        add_to_log(base_name)
        add_to_log(format(Sys.time(), '%Y-%m-%d %H:%M:%S'))
        add_to_log('-------------------------------------------')
        print(LOG)
        write_to_log(log_dir)

        print('Reading points')
        las <- readLAS(
          f,
          filter='-drop_withheld',
          select='xyzirc'
          )
        
        #filter duplicate points and check integrity
        las <- filter_duplicates(las)
        print('las chaeck..')
        las_check_output <- capture.output(las_check(las))
        add_to_log('blah blah')
        print('adding to log')
        if (is.character(las_check_output)) {
          add_to_log(las_check_output)
        } else {
          add_to_log("Error: las_check did not return a character vector.")
        }
        write_to_log(log_dir)

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
        
        chm_path <- file.path(chm_dir, paste0(base_name, '.tif'))
        writeRaster(smoothed, chm_path, overwrite=TRUE)
        
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
        ttops$uniqueID <- mapply(
          generate_custom_id,
          utm = utm_zone,
          geometry = ttops$geometry
        )
        
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
            ttops_dir,
            paste0(base_name, '.gpkg')
          )
        )
        
        print('writing crowns')
        st_write(
          crowns,
          file.path(
            crowns_dir,
            paste0(base_name, '.gpkg')
          )
        )

        print('writing segmented point cloud to segs_dir')
        
        writeLAS(
          segs,
          file.path(
            segs_dir,
            paste0(base_name, '.laz')
            )
        )

        add_to_log(paste(
          f,
          'succesfully processed at: ',
          format(Sys.time(), '%Y-%m-%d %H:%M:%S')
          )
        )
        
        write_to_log(log_dir)
      

      }, error = function(e) {
        add_to_log('- - - - - - - - - - - - - - - - - - - - - - -')
        add_to_log(paste('problem with tile', f))
        add_to_log(format(Sys.time(), '%Y-%m-%d %H:%M:%S'))
        add_to_log(e$message)
        add_to_log('- - - - - - - - - - - - - - - - - - - - - - -')
        write_to_log(log_dir)
      })
    } 
  }
}

# ----------------------------------------------------------

# counter (start at one since this is R)
i <- 1

# list and count completed tiles
complete <- path_ext_remove(path_file(dir_ls(crowns_dir, glob = "*.gpkg")))
n_done <- length(complete)
n_laz <- length(path_ext_remove(path_file(dir_ls(laz_dir, glob = "*.laz"))))

while (i < 10) {
  if (n_laz > n_done) {
    main(complete)
    complete <- path_ext_remove(path_file(dir_ls(crowns_dir, glob = "*.gpkg")))
    n_done <- length(complete)
    n_laz <- length(path_ext_remove(path_file(dir_ls(laz_dir, glob = "*.laz"))))
    i <- 1
  } else {
    # check every minute in case list is still growing
    complete <- path_ext_remove(path_file(dir_ls(crowns_dir, glob = "*.gpkg")))
    n_done <- length(complete)
    n_laz <- length(path_ext_remove(path_file(dir_ls(laz_dir, glob = "*.laz"))))
    print(paste('Waiting', i, 'of 10 minutes for list to grow.'))
    i <- i + 1
    Sys.sleep(60)
  }
}

