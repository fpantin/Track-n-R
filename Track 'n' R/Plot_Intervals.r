################################################################################
#                                                                              #
#         Plot growth, surface fitting or division at regular intervals        #
#                          from raw PointTracker data                          #
#                                                                              #
#                             Florent Pantin, 2013                             #
#                                                                              #
################################################################################


## Load all useful functions and libraries ##

source("D:/Track_n_R.r")


## Load parameters ##

para <- read.table("D:/Parameters.txt", sep = "/")

wd <- gsub("\\\\", "/", as.character(para[1,]))
csvFile <- gsub("\\\\", "/", as.character(para[2,]))
ifelse (paste(wd, "Data", sep = "") == substr(csvFile, 1, nchar(paste(wd, "Data", sep = ""))),
        csvFile <- gsub(paste(wd, "Data/", sep = ""), "", csvFile, fixed = T),
        stop("The csv file should be located in the 'Data' folder of the project directory.", call. = F))

if (para[3,] == "NULL") { meanAngle <- NULL } else { meanAngle <- as.numeric(as.character(para[3,])) }
alignToPetioleLaminaBoundary <- as.logical(para[4,]) # Default: T
if (para[5,] == "NULL") { cellAtPetioleLaminaBoundary <- NULL } else { cellAtPetioleLaminaBoundary <- as.integer(as.character(para[5,])) }
if (para[6,] == "NULL") { leafShape <- NULL } else { leafShape <- as.logical(para[6,]) }
before <- as.logical(para[7,]) # Default: F
k.growth <- as.integer(as.character(para[8,])) # Default: 2
black <- as.logical(para[9,]) # Default: T
xlim <- NULL # this parameter is not allowed to be defined by the user
if (para[10,] == "NULL") { ylim <- NULL } else { ylim <- c(as.numeric(as.character(para[10,])), as.numeric(as.character(para[11,]))) }
if (para[12,] == "NULL") { scaleBar <- NULL } else { scaleBar <- as.numeric(as.character(para[12,])) }
tick <- as.logical(para[13,]) # Default: F
if (para[14,] == "NULL") { Percent <- NULL } else { Percent <- as.logical(para[14,]) }
round.zlim <- as.logical(para[15,]) # Default: T
if (para[16,] == "NULL") { zlim <- NULL } else { zlim <- c(as.numeric(as.character(para[16,])), as.numeric(as.character(para[17,]))) }
fix.min <- as.logical(para[18,]) # Default: T
fix.max <- as.logical(para[19,]) # Default: T
colorPaletteOfGFtbox <- as.logical(para[20,]) # Default: T
growthScale <- as.logical(para[21,]) # Default: T
drawTicks <- as.logical(para[22,]) # Default: T
if (para[23,] == "NULL") { txt.legend <- NULL } else { txt.legend <- as.character(para[23,]) }
if (para[24,] == "NULL") { anisotropy <- NULL } else { anisotropy <- as.logical(para[24,]) }
aniso.threshold <- as.numeric(as.character(para[25,])) # Default: 0.05
if (para[26,] == "NULL") { aniso.lwd <- NULL } else { aniso.lwd <- as.numeric(as.character(para[26,])) }
aniso.lwd.constant <- as.logical(para[27,]) # Default: F
fit.last <- as.logical(para[28,]) # Default: F

div.parameter <- as.character(para[29,]) # Default: "div&comp"
round.zlimDiv <- as.logical(para[30,]) # Default: T
if (para[31,] == "NULL") { zlimDiv <- NULL } else { zlimDiv <- c(as.numeric(as.character(para[31,])), as.numeric(as.character(para[32,]))) }
fix.minDiv <- as.logical(para[33,]) # Default: T
fix.maxDiv <- as.logical(para[34,]) # Default: T
colorPaletteOfGFtboxDiv <- as.logical(para[35,]) # Default: T
growthScaleDiv <- as.logical(para[36,]) # Default: T
drawTicksDiv <- as.logical(para[37,]) # Default: T
if (para[38,] == "NULL") { txt.legendDiv <- NULL } else { txt.legendDiv <- as.character(para[38,]) }

n.pred <- as.numeric(as.character(para[39,])) # Default: 100
polynomDegree <- as.numeric(as.character(para[40,])) # Default: 3
rangeParameter <- as.numeric(as.character(para[41,])) # Default: 0.001
if (para[42,] == "NULL") { contour.levels <- NULL } else { contour.levels <- seq(from = as.numeric(as.character(para[42,])), to = as.numeric(as.character(para[43,])), by = as.numeric(as.character(para[44,]))) }
leafShapeKriging <- as.logical(para[45,]) # Default: F
if (para[46,] == "NULL") { Reslim <- NULL } else { Reslim <- c(as.numeric(as.character(para[46,])), as.numeric(as.character(para[47,]))) }

gro <- as.logical(para[48,]) # Default: T
kri <- as.logical(para[49,]) # Default: T
div <- as.logical(para[50,]) # Default: T
interval <- as.numeric(as.character(para[51,])) # Default: 24
Image0 <- as.character(para[52,]) # Default: first image
initial.stage <- as.logical(para[53,]) # Default: T
ini <- as.logical(para[54,]) # Default: T
show.cell.number <- as.logical(para[55,]) # Default: F


## Set working directory ##

setwd(wd)


## Import raw data ##

c(InfoVertices, Vertices, Cells, Divisions) := import.rawData(csvFile)


## Decide whether or not plotting leaf outlines ##

if (is.null(leafShape))
  {
  check_leafShape <- T
  for (img in as.character(InfoVertices$Images))
    {
    if (!file.exists(paste("Leaf Shape/", substr(img, 1, nchar(img) - 3), "txt", sep = ""))) { check_leafShape <- F }
    }
  ifelse (check_leafShape, leafShape <- T, leafShape <- F)
  }


## Compute 'meanAngle' if not provided ##

if (is.null(meanAngle))
  {
  ifelse (leafShape, meanAngle <- get.meanAngle(InfoVertices), meanAngle <- 0)
  }


##  Find dimensions of the last image if required ##

if (fit.last)
  {
  time0 <- InfoVertices$Time[InfoVertices$Images == Image0]
  timeFinal <- InfoVertices$Time[nrow(InfoVertices)]
  wanted.intervals <- seq(time0, timeFinal, interval)
  if (wanted.intervals[length(wanted.intervals)] + interval < timeFinal + interval/3)
    {
    wanted.intervals <- c(wanted.intervals, wanted.intervals[length(wanted.intervals)] + interval)
    }
  actualTime1.last <- InfoVertices$Time[which.min((InfoVertices$Time - wanted.intervals[length(wanted.intervals)-1])^2)]
  actualTime2.last <- InfoVertices$Time[which.min((InfoVertices$Time - wanted.intervals[length(wanted.intervals)])^2)]
  Image1.last <- as.character(InfoVertices$Image[InfoVertices$Time == actualTime1.last])
  Image2.last <- as.character(InfoVertices$Image[InfoVertices$Time == actualTime2.last])
  
  c(Growth, Shapes) := recompute.GrowthAndShapes(Cells, Divisions, Vertices, InfoVertices, Image1.last, Image2.last, meanAngle)

  arg <- list(Image1 = Image1.last, Image2 = Image2.last, before = before, meanAngle = meanAngle, leafShape = leafShape,
              alignToPetioleLaminaBoundary = alignToPetioleLaminaBoundary,
              Shapes = Shapes, Cells = Cells, Vertices = Vertices, InfoVertices = InfoVertices)
  if (!is.null(cellAtPetioleLaminaBoundary)) { arg$cellAtPetioleLaminaBoundary <- cellAtPetioleLaminaBoundary }
  c(VertX, VertY, ShapesX, ShapesY, LeafShape, LeafSpline) := do.call(scale.Objects, arg)

  if(leafShape)
    {
    xlim <- c(-1.01*max(abs(LeafSpline$x)), 1.01*max(abs(LeafSpline$x)))
    ylim <- c(min(LeafSpline$y) - 0.01*abs(min(LeafSpline$y)), 1.01*max(LeafSpline$y))
    }
  else
    {
    xlim <- c(min(VertX, na.rm = T), max(VertX, na.rm = T))
    ylim <- c(min(VertY, na.rm = T) - 0.01*abs(min(VertY, na.rm = T)), 1.01*max(VertY, na.rm = T))
    }
  }
  
  
## Select the initial image from which divisions start to be processed ##

ifelse(ini, ini <- which(as.character(InfoVertices$Image) == Image0), ini <- 1)


## Arguments ##

  # Compulsory arguments #

  arg <- list(Cells = Cells, Divisions = Divisions, Vertices = Vertices, InfoVertices = InfoVertices, meanAngle = meanAngle,
              alignToPetioleLaminaBoundary = alignToPetioleLaminaBoundary,
              interval = interval, Image0 = Image0, initial.stage = initial.stage,
              leafShapeKriging = F,
              before = before, k.growth = k.growth, PNG = T, black = black, tick = tick,
              round.zlim = round.zlim, fix.min = fix.min, fix.max = fix.max,
              colorPaletteOfGFtbox = colorPaletteOfGFtbox, growthScale = growthScale, drawTicks = drawTicks,
              aniso.threshold = aniso.threshold, aniso.lwd.constant = aniso.lwd.constant,
              n.pred = n.pred, polynomDegree = polynomDegree, rangeParameter = rangeParameter,
              cellularScale = F, plotResidual = F, exportCellValues = F,
              div.parameter =  div.parameter, round.zlimDiv = round.zlimDiv, fix.minDiv = fix.minDiv, fix.maxDiv = fix.maxDiv,
              colorPaletteOfGFtboxDiv = colorPaletteOfGFtboxDiv, growthScaleDiv = growthScaleDiv, drawTicksDiv = drawTicksDiv,
              wd = wd, File = csvFile, ini = ini, show.cell.number = show.cell.number) # edit 2016/01/15

            
  # Optional arguments #

  if (!is.null(leafShape)) { arg$leafShape <- leafShape }
  if (!is.null(cellAtPetioleLaminaBoundary)) { arg$cellAtPetioleLaminaBoundary <- cellAtPetioleLaminaBoundary }
  if (!is.null(xlim)) { arg$xlim <- xlim }
  if (!is.null(ylim)) { arg$ylim <- ylim }
  if (!is.null(scaleBar)) { arg$scaleBar <- scaleBar }
  if (!is.null(Percent)) { arg$Percent <- Percent }
  if (!is.null(zlim)) { arg$zlim <- zlim }
  if (!is.null(txt.legend)) { arg$txt.legend <- txt.legend }
  if (!is.null(anisotropy)) { arg$anisotropy <- anisotropy }
  if (!is.null(aniso.lwd)) { arg$aniso.lwd <- aniso.lwd }
  if (!is.null(contour.levels)) { arg$contour.levels <- contour.levels }
  if (!is.null(zlimDiv)) { arg$zlimDiv <- zlimDiv }
  if (!is.null(txt.legendDiv)) { arg$txt.legendDiv <- txt.legendDiv }


## Loop on each plot type ##

plot.types <- vector()
if (gro) { plot.types <- c(plot.types, "growth") }
if (kri) { plot.types <- c(plot.types, "kriging") }
if (div) { plot.types <- c(plot.types, "division") }

my_img <- 0

for (plot.type in plot.types)
  {
  # Set the current plot type within the arguments #
  arg$plot.type <- plot.type

  # Plot and save the growth map #
  img_count <- do.call(plot.intervals, arg)
  
  # Counts the plotted images #
  my_img <- my_img + img_count

  # Additional maps for kriging #
  if (plot.type == "kriging")
    {
    arg$cellularScale <- T
    img_count <- do.call(plot.intervals, arg)
    my_img <- my_img + img_count

    arg$plotResidual <- T
    arg$zlim <- NULL
    if (!is.null(Reslim)) { arg$zlim <- Reslim }
    img_count <- do.call(plot.intervals, arg)
    my_img <- my_img + img_count

    if (leafShapeKriging)
       {
       arg$cellularScale <- F
       arg$plotResidual <- F
       arg$leafShapeKriging <- T
       arg$zlim <- NULL
       if (!is.null(zlim)) { arg$zlim <- zlim }
       img_count <- do.call(plot.intervals, arg)
       my_img <- my_img + img_count
       }
    }
  }


## Save the number of images which have been plotted ##

write.table(my_img, "D:/count_img.txt", row.names = F, col.names = F)
