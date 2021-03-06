################################################################################
#                                                                              #
#                Replot growth from processed PointTracker data                #
#                                                                              #
#                             Florent Pantin, 2013                             #
#                                                                              #
################################################################################


## Load all useful functions and libraries ##

source("D:/Track_n_R.r")


## Load parameters ##

para <- read.table("D:/Parameters.txt", sep = "/")

wd <- gsub("\\\\", "/", as.character(para[1,]))
xlsFile <- gsub("\\\\", "/", as.character(para[2,]))
ifelse (paste(wd, "Growth", sep = "") == substr(xlsFile, 1, nchar(paste(wd, "Growth", sep = ""))),
        xlsFile <- gsub(paste(wd, "Growth/", sep = ""), "", xlsFile, fixed = T),
        stop("The xls file should be located in the 'Growth' folder of the project directory.", call. = F))

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


## Set working directory ##

setwd(wd)


## Find used images ##

UsedImages <- read.table(paste("Growth", xlsFile, sep = "/"), sep = ",", nrows = 1, skip = 3)
UsedImages <- UsedImages[!is.na(UsedImages) & UsedImages != ""][-1]


## Decide whether or not plotting leaf outlines ##

if (is.null(leafShape))
  {
  check_leafShape <- T
  if (is.null(meanAngle)) # if 'meanAngle' is not provided, need to check all images as they will be used by 'get.meanAngle()'
    {
    dataImport <- read.table(paste("Growth", xlsFile, sep = "/"), sep = ",", header = F, fill = T, blank.lines.skip = F, col.names = 1:1000)
    Info <- dataImport[which(dataImport[1] == "Data") + 2:5, 1:length(which(!is.na(dataImport[(which(dataImport[1] == "Data") + 5), ])))]
    InfoVertices <- data.frame(Images = NULL, ShiftX = NULL, ShiftY = NULL, AngleShift = NULL, Time = NULL, ScalingX = NULL, ScalingY = NULL)
    for (imag in seq(2, ncol(Info), 2)) { InfoVertices <- rbind(InfoVertices, data.frame(Images = as.character(Info[1, imag]), ShiftX = as.numeric(as.character(Info[2, imag])), ShiftY = as.numeric(as.character(Info[2, imag+1])), AngleShift = as.numeric(as.character(Info[3, imag])), Time = as.numeric(as.character(Info[3, imag+1])), ScalingX = as.numeric(as.character(Info[4, imag])), ScalingY = as.numeric(as.character(Info[4, imag+1])))) }
    InfoVertices$Images <- as.factor(InfoVertices$Images)
    for (img in as.character(InfoVertices$Images))
      {
      if (!file.exists(paste("Leaf Shape/", substr(img, 1, nchar(img) - 3), "txt", sep = ""))) { check_leafShape <- F }
      }
    }  
  else # if 'meanAngle' is provided, only the desired images are checked because 'get.meanAngle()' won't be called
    {
    for (my_img in 1:nrow(UsedImages))
      {
      img <- as.character(UsedImages[my_img, 1])
      if (!file.exists(paste("Leaf Shape/", substr(img, 1, nchar(img) - 3), "txt", sep = ""))) { check_leafShape <- F }
      }
    }
  ifelse (check_leafShape, leafShape <- T, leafShape <- F)
  }


##  Find dimensions of the last image if required ##

if (fit.last)
  {
  Image1.last <- UsedImages[length(UsedImages)-1]
  Image2.last <- UsedImages[length(UsedImages)]

  arg <- list(File = xlsFile, Image1 = Image1.last, Image2 = Image2.last, leafShape = leafShape)
  if (!is.null(leafShape)) { arg$leafShape <- leafShape }
  c(Image.dim, InfoVertices, meanAngle, Shapes, Vertices, Cells, Divisions, GrowthComputationParameters, Growth, GrowthWall) := do.call(import.processedData, arg)

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

  
## Loop on each image couple ##

for (my_img in 1:(length(UsedImages)-1))
  {

  # Name of the images #
  
  Image1 <- UsedImages[my_img]
  Image2 <- UsedImages[my_img+1]


  # Import processed data #

  arg <- list(File = xlsFile, Image1 = Image1, Image2 = Image2, leafShape = leafShape)
  if (!is.null(meanAngle)) { arg$meanAngle <- meanAngle }

  c(Image.dim, InfoVertices, meanAngle, Shapes, Vertices, Cells, Divisions, GrowthComputationParameters, Growth, GrowthWall) := do.call(import.processedData, arg)


  # Compulsory arguments #

  arg <- list(InfoVertices = InfoVertices, Vertices = Vertices, Cells = Cells, Divisions = Divisions,
              meanAngle = meanAngle, Growth = Growth, Shapes = Shapes,
              Image1 = Image1, Image2 = Image2,
              alignToPetioleLaminaBoundary = alignToPetioleLaminaBoundary,
              before = before, k.growth = k.growth, PNG = T, black = black, tick = tick,
              round.zlim = round.zlim, fix.min = fix.min, fix.max = fix.max,
              colorPaletteOfGFtbox = colorPaletteOfGFtbox, growthScale = growthScale, drawTicks = drawTicks,
              aniso.threshold = aniso.threshold, aniso.lwd.constant = aniso.lwd.constant)


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


  # Plot and save the growth map #

  if (!file.exists("Graphical Outputs")) { dir.create("Graphical Outputs") }

  if      (k.growth == 2) { k.g <- "karea" }
  else if (k.growth == 3) { k.g <- "kmaj" }
  else if (k.growth == 4) { k.g <- "kmin" }
  else if (k.growth == 5) { k.g <- "theta" }
  else if (k.growth == 6) { k.g <- "phi" }
  else if (k.growth == 7) { k.g <- "anisotropy" }
  else if (k.growth == 8) { k.g <- "kpertoml" }
  else if (k.growth == 9) { k.g <- "kml" }
  else if (k.growth == 10) { k.g <- "cellarea" }
  else if (k.growth == 11) { k.g <- "sectorarea" }

  actualTime1 <- InfoVertices$Time[InfoVertices$Images == Image1]
  actualTime2 <- InfoVertices$Time[InfoVertices$Images == Image2]

  filename <- paste("Graphical Outputs/Growth from PT__", k.g, "__",
                    substr(Image1, 1, nchar(Image1) - 4), " - ", substr(Image2, 1, nchar(Image2) - 4), " (",
                    round(actualTime1, 1), " - ", round(actualTime2, 1), " h).png", sep = "")
  if (nchar(filename) > 259 - nchar(wd)) { filename <- paste(substr(filename, 1, 259 - nchar(wd) - 4), ".png", sep = "") }
  if (file.exists(filename)) { file.remove(filename) }
  size.in.pix <- 2000
  png(filename, width = size.in.pix, height = size.in.pix, res = 0.15*size.in.pix)
  do.call(plot.growth, arg)
  dev.off()
  }


## Save the number of images which have been plotted ##

write.table(my_img, "D:/count_img.txt", row.names = F, col.names = F)
