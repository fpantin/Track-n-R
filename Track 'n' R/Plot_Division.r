################################################################################
#                                                                              #
#         Plot cell division and competence from raw PointTracker data         #
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
div.parameter <- as.character(para[8,]) # Default: "div&comp"
black <- as.logical(para[9,]) # Default: T
xlim <- NULL # this parameter is not allowed to be defined by the user
if (para[10,] == "NULL") { ylim <- NULL } else { ylim <- c(as.numeric(as.character(para[10,])), as.numeric(as.character(para[11,]))) }
if (para[12,] == "NULL") { scaleBar <- NULL } else { scaleBar <- as.numeric(as.character(para[12,])) }
tick <- as.logical(para[13,]) # Default: F
round.zlim <- as.logical(para[14,]) # Default: T
if (para[15,] == "NULL") { zlim <- NULL } else { zlim <- c(as.numeric(as.character(para[15,])), as.numeric(as.character(para[16,]))) }
fix.min <- as.logical(para[17,]) # Default: T
fix.max <- as.logical(para[18,]) # Default: T
colorPaletteOfGFtbox <- as.logical(para[19,]) # Default: T
growthScale <- as.logical(para[20,]) # Default: T
drawTicks <- as.logical(para[21,]) # Default: T
if (para[22,] == "NULL") { txt.legend <- NULL } else { txt.legend <- as.character(para[22,]) }
fit.last <- as.logical(para[23,]) # Default: F
ini <- as.logical(para[24,]) # Default: F
show.cell.number <- as.logical(para[25,]) # Default: F

## Set working directory ##

setwd(wd)


## Import raw data ##

c(InfoVertices, Vertices, Cells, Divisions) := import.rawData(csvFile)


## Decide whether or not plotting leaf outlines ##

ProcessedImages <- read.table("D:/ImageList.txt")
if (nrow(ProcessedImages) == 1)
  {
  ProcessedImages[2, 1] <- ProcessedImages[1, 1]
  #if (div.parameter %in% c("Area_At_Division", "Symmetry", "Cell_Cycle_Duration")) { growthScale <- F } # ignore colour scale because all cells will be empty
  if (div.parameter %in% c("Area_At_Division", "Symmetry", "Cell_Cycle_Duration", "Number_Of_Divisions")) { growthScale <- F } # ignore colour scale because all cells will be empty # edit 2015/12/15
  }
if (is.null(leafShape))
  {
  check_leafShape <- T
  if (is.null(meanAngle)) # if 'meanAngle' is not provided, need to check all images as they will be used by 'get.meanAngle()'
    {
    for (img in as.character(InfoVertices$Images))
      {
      if (!file.exists(paste("Leaf Shape/", substr(img, 1, nchar(img) - 3), "txt", sep = ""))) { check_leafShape <- F }
      }
    }  
  else # if 'meanAngle' is provided, only the desired images are checked because 'get.meanAngle()' won't be called
    {
    for (my_img in 1:nrow(ProcessedImages))
      {
      img <- as.character(ProcessedImages[my_img, 1])
      if (!file.exists(paste("Leaf Shape/", substr(img, 1, nchar(img) - 3), "txt", sep = ""))) { check_leafShape <- F }
      }
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
  Image1.last <- as.character(ProcessedImages[nrow(ProcessedImages)-1, 1])
  Image2.last <- as.character(ProcessedImages[nrow(ProcessedImages), 1])
  
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

ifelse(ini, ini <- which(as.character(InfoVertices$Image) == ProcessedImages[1, 1]), ini <- 1)


## Loop on each image couple ##

for (my_img in 1:(nrow(ProcessedImages)-1))
  {

  # Name of the images #
  
  Image1 <- as.character(ProcessedImages[my_img, 1])
  Image2 <- as.character(ProcessedImages[my_img+1, 1])


  # Compute divisions #
  
  c(Growth, Shapes) := recompute.GrowthAndShapes(Cells, Divisions, Vertices, InfoVertices, Image1, Image2, meanAngle)
  arg <- list(InfoVertices = InfoVertices, Vertices = Vertices, Cells = Cells, Divisions = Divisions,
              meanAngle = meanAngle, Growth = Growth, Shapes = Shapes, alignToPetioleLaminaBoundary = alignToPetioleLaminaBoundary, 
              Image1 = Image1, Image2 = Image2, leafShape = leafShape, before = before, File = csvFile, ini = ini)
  if (!is.null(cellAtPetioleLaminaBoundary)) { arg$cellAtPetioleLaminaBoundary <- cellAtPetioleLaminaBoundary }
  Div <- do.call(process.DivWithInterval, arg)
  

  # Compulsory arguments #

  arg <- list(InfoVertices = InfoVertices, Vertices = Vertices, Cells = Cells, Divisions = Divisions,
              meanAngle = meanAngle, Shapes = Shapes, Div = Div,
              alignToPetioleLaminaBoundary = alignToPetioleLaminaBoundary,
              Image1 = Image1, Image2 = Image2,
              before = before, PNG = T, black = black, tick = tick,
              div.parameter = div.parameter, round.zlim = round.zlim, fix.min = fix.min, fix.max = fix.max,
              #colorPaletteOfGFtbox = colorPaletteOfGFtbox, growthScale = growthScale, drawTicks = drawTicks)
              colorPaletteOfGFtbox = colorPaletteOfGFtbox, growthScale = growthScale, drawTicks = drawTicks, ini = ini, show.cell.number = show.cell.number) # edits 2014/09/11 and 2016/01/15


  # Optional arguments #

  if (!is.null(leafShape)) { arg$leafShape <- leafShape }
  if (!is.null(cellAtPetioleLaminaBoundary)) { arg$cellAtPetioleLaminaBoundary <- cellAtPetioleLaminaBoundary }
  if (!is.null(xlim)) { arg$xlim <- xlim }
  if (!is.null(ylim)) { arg$ylim <- ylim }
  if (!is.null(scaleBar)) { arg$scaleBar <- scaleBar }
  if (!is.null(zlim)) { arg$zlim <- zlim }
  if (!is.null(txt.legend)) { arg$txt.legend <- txt.legend }


  # Plot and save the growth map #

  if (!file.exists("Graphical Outputs")) { dir.create("Graphical Outputs") }

  if      (div.parameter == "div&comp") { d.p <- "division and competence" }
  else if (div.parameter == "div") { d.p <- "division only" }
  else if (div.parameter == "comp") { d.p <- "competence only" }
  else if (div.parameter == "progeny") { d.p <- "progeny" } # edit 2014/09/11
  else if (div.parameter == "CellArea") { d.p <- "area of all cells" }
  else if (div.parameter %in% c("AreaTime1", "AreaTime2")) { d.p <- "area of all sectors" }
  else if (div.parameter == "Area_At_Division") { d.p <- "cell area at division" }
  else if (div.parameter == "Symmetry") { d.p <- "symmetry of division" }
  else if (div.parameter == "Cell_Cycle_Duration") { d.p <- "duration of cell cycle" }
  else if (div.parameter == "Number_Of_Divisions") { d.p <- "number of divisions" } # edit 2015/12/15

  actualTime1 <- InfoVertices$Time[InfoVertices$Images == Image1]
  actualTime2 <- InfoVertices$Time[InfoVertices$Images == Image2]

  filename <- paste("Graphical Outputs/Division__", d.p, "__",
                    substr(Image1, 1, nchar(Image1) - 4), " - ", substr(Image2, 1, nchar(Image2) - 4), " (",
                    round(actualTime1, 1), " - ", round(actualTime2, 1), " h).png", sep = "")
  if (nchar(filename) > 259 - nchar(wd)) { filename <- paste(substr(filename, 1, 259 - nchar(wd) - 4), ".png", sep = "") }
  if (file.exists(filename)) { file.remove(filename) }
  size.in.pix <- 2000
  png(filename, width = size.in.pix, height = size.in.pix, res = 0.15*size.in.pix)
  do.call(plot.division, arg)
  dev.off()
  }


## Save the number of images which have been plotted ##

write.table(my_img, "D:/count_img.txt", row.names = F, col.names = F)
