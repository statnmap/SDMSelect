#' A function to crop white margins of a PNG image
#'
#' @param x path to the PNG image
#' @param pixel number of white pixels lines to keep
#' @export

Rm_WhiteMargins <- function(x, pixel = 2)
{
  # Cut the output image to remove dirty white margins from corrplot
  img <- png::readPNG(x)

  img.test.row <- apply(img, 3, function(layer) {
    apply(layer, 1, function(i) {(sum(i != 1) > 0)})
  }) %>%
    apply(., 1, function(i) {(sum(i) > 0)})
  rowMin <- max(min(which(img.test.row[1:round(length(img.test.row) / 2)])) - (1 + pixel), 1)
  rowMax <- min(max(c(1:length(img.test.row))[
    round(length(img.test.row) / 2) : length(img.test.row)][
      which(img.test.row[(length(img.test.row) / 2) : length(img.test.row)])]) + 1 + pixel,
    length(img.test.row))

  img.test.col <- apply(img, 3, function(layer) {
    apply(layer, 2, function(i) {(sum(i != 1) > 0)})
  }) %>%
    apply(., 1, function(i) {(sum(i) > 0)})
  colMin <- max(min(which(img.test.col[1:round(length(img.test.col) / 2)])) - (1 + pixel), 1)
  colMax <- min(max(c(1:length(img.test.col))[
    round(length(img.test.col) / 2) : length(img.test.col)][
      which(img.test.col[(length(img.test.col) / 2) : length(img.test.col)])]) + 1 + pixel,
    length(img.test.col))

  # Remove rows and cols with white pixels from the original image
  img <- img[rowMin:rowMax, colMin:colMax,]
  png::writePNG(img, target = paste0(gsub(".png$", "", x), "_crop.png"))
  rm(img)
}

#' Find the most precise data of a vector
#'
#' @param x vector of values
#' @return precision as a power of 10: 10^n
#' @export

prec_data <- function(x) {

  # Get the most precise number
  x.dec <- strsplit(as.character(x), split = ".", fixed = TRUE)

  if (max(lengths(x.dec)) == 1) {
    # if there are no number after comma
    n.plus <- max(nchar(as.character(x)))
    x2 <- x / 10^n.plus

    x.dec <- strsplit(as.character(x2), split = ".", fixed = TRUE)

    n.minus <- max(unlist(lapply(x.dec,
                                 function(y) {
                                   nchar(as.character(as.numeric(y[2])))
                                 })), na.rm = TRUE)

    n <- n.plus - n.minus
  } else {

    n <- -1 * max(unlist(lapply(x.dec,
                                function(y) {
                                  nchar(as.character(as.numeric(y[2])))
                                })), na.rm = TRUE)
  }
  return(n)
}

#' Find area extent
#' Verify if this function is not in my spatial script
#' @param x SpatialPointsDataFrame
#' @param saveWD directory where to save Extent
#'
#' @importFrom raster xmin xmax ymin ymax
#'
#' @export

CreateExtent <- function(x, saveWD) {
  # Create / Load extent of the studied area based on dataset
  ZoneExt <- raster::extent(x)*1.2 # Can be drawn by hand : drawExtent()
  save(ZoneExt, file = paste0(saveWD, "/AreaExtent.RData"))
  load(paste0(saveWD, "/AreaExtent.RData"))
  # Transform as a polygon
  mpol1 <- rbind(c(xmin(ZoneExt), ymin(ZoneExt)), c(xmin(ZoneExt), ymax(ZoneExt)),
                 c(xmax(ZoneExt), ymax(ZoneExt)), c(xmax(ZoneExt), ymin(ZoneExt)))
  mpol1
  ZoneEtude_pol <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(mpol1)), "Study")),
                                       proj4string = sp::CRS(sp::proj4string(x)))
}

#' Find an equilibrated layout for a figure according to number of panels to include
#'
#' @param x number of panels of the future plot
#'
#' @details As screens are larger than height, nrow is set not to exceed ncol.
#' Rules are : ncol >= nrow, ncol <= nrow + 2, ntot >= x
#' 100*100 is the maximum window. It is surely too big for a screen, but who
#' knows what you do with this function... Be sure your window is big enough !
#'
#'
#' @return
#' nrow
#' @export

Fig_split <- function(x) {

  dim.xy <- expand.grid(nrow = 1:min(x,10), ncol = 1:min(x,10)) %>%
    dplyr::mutate(ntot = nrow * ncol) %>%
    dplyr::mutate(diff.fig = ncol - nrow) %>%
    dplyr::mutate(diff.x = ntot - x) %>%
    dplyr::mutate(diff.sum = diff.x) %>%
    dplyr::filter(ncol >= nrow, ncol <= nrow + 2, ntot >= x) %>%
    dplyr::arrange(diff.sum, ntot)

  return(as.list(dim.xy[1, 1:2]))
}


#' @title Test if projections are equals
#'
#' @description Some projection are the same but written differently
#' or with less useless arguments.
#' This function test if two projections are equal while accounting for
#' multiple writings, in particular when EPSG can not be retrieved
#' @param x proj4string
#' @param y proj4string

is_ProjEqual <- function(x, y)
{
  res <- FALSE
  # Test for supported epsg
  x.EPSG <- as.numeric(rgdal::showEPSG(x))
  y.EPSG <- as.numeric(rgdal::showEPSG(y))
  epsg <- FALSE
  if (!is.na(x.EPSG) & !is.na(y.EPSG)) {
    if (x.EPSG == y.EPSG) {
      epsg <- TRUE
    }
  }

  if (x == y | epsg) {
    res <- TRUE
  } else {
    # Test for other cases
    x.split.all <- x.split <- strsplit(x, " ")[[1]]
    y.split.all <- y.split <- strsplit(y, " ")[[1]]
    in.lat <- 0
    if (length(grep("lat_1", x.split.all)) != 0
        & length(grep("lat_2", x.split.all)) != 0) {
      in.lat <- 1
      x.split <- x.split.all[-c(grep("lat_1", x.split.all),
                                grep("lat_2", x.split.all))]
    }
    if (length(grep("lat_1", y.split.all)) != 0
        & length(grep("lat_2", y.split.all)) != 0) {
      in.lat <- in.lat + 1
      y.split <- y.split.all[-c(grep("lat_1", y.split.all),
                                grep("lat_2", y.split.all))]
    }

    # Test: ordered differently or included in
    l.xy.all <- which(x.split.all %in% y.split.all)
    l.yx.all <- which(y.split.all %in% x.split.all)
    if (in.lat != 0) {
      l.xy <- which(x.split %in% y.split)
      l.yx <- which(y.split %in% x.split)
    }
    if (length(l.xy.all) == length(x.split.all) |
        length(l.yx.all) == length(y.split.all))
    {
      res <- TRUE
    } else if (in.lat != 0) {
      if (length(l.xy) == length(x.split) |
          length(l.yx) == length(y.split)) {
        # One projection is included in the other,
        # this is just a problem of name
        res <- TRUE
        # Verify if there are "lat_" to be tested
        if (in.lat == 1) {
          # projections are different because only one has reference points
          res <- FALSE
        }
        if (in.lat == 2) {
          # Must test if there is no inversion of lat_1 and lat_2
          # For lambert 93, test with :
          # x <- "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
          # y <- "+proj=lcc +lat_2=49 +lat_1=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs"
          message("Two projections having lat_1 and lat_2 switched are currently considered different.")
          res <- FALSE
        }
      }
    }
  }
  return(res)
}
