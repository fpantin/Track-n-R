################################################################################
#                                                                              #
#           Creates maps of cell lineages from raw PointTracker data           #
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
if (para[4,] == "NULL") { leafShape <- NULL } else { leafShape <- as.logical(para[4,]) }
if (para[5,] == "NULL") { target.cells <- NULL } else { target.cells <- as.numeric(as.character(para[5,])) }
all.cells <- as.logical(para[6,])
all.lineages <- as.logical(para[7,])
combinedPDF <- as.logical(para[8,])
reorient <- as.logical(para[9,])
n.col <- as.integer(as.character(para[10,]))
main.plot <- as.character(para[11,])
display.vertices <- as.logical(para[12,])
create.csv <- as.logical(para[13,])
optimal.interval <- as.numeric(as.character(para[14,]))
ini <- as.logical(para[15,]) # Default: F
makeAVI <- as.logical(para[16,]) # Default: F
compr <- as.character(para[17,]) # Default: JPEG
fps <- as.integer(as.character(para[18,])) # Default: 10


## Set working directory ##

setwd(wd)


## Import raw data ##

c(InfoVertices, Vertices, Cells, Divisions) := import.rawData(csvFile)


## Check if the cells should be as in PointTracker or reoriented according to Track 'n' R ##

ProcessedImages <- read.table("D:/ImageList.txt")
if (reorient)
  {
  # Check leaf outlines
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

  # Compute 'meanAngle' if not provided
  if (is.null(meanAngle))
    {
    ifelse (leafShape, meanAngle <- get.meanAngle(InfoVertices), meanAngle <- 0)
    }
  }

if (!reorient)
  {
  leafShape <- F
  meanAngle <- 0
  }
  

## Select the initial image from which divisions start to be processed ##

ifelse(ini, ini <- which(as.character(InfoVertices$Image) == ProcessedImages[1, 1]), ini <- 1)


## Compute the number of lines in the montage ##

n.row <- ceiling(nrow(ProcessedImages)/n.col)


## Create directory for graphics ##

if (!file.exists("Lineages")) { dir.create("Lineages") }


## Set temporary output directory for cell coordinates ##

outputDir <- "D:/CellCoord_temp"


## Initialize ImageJ ##

  # Load the R library
  #library("rJava") # already loaded by sourcing 'PointTracker_to_R.r'

  # Initialize Java with a heap size of 1950 Mo
  .jinit(parameters = "-Xmx1950m")

  # Add the ImageJ *.jar to the classpath
  .jaddClassPath("C:/Program Files/ImageJ/ij.jar")

  # Create a new object
  IJ <- .jnew("ij/IJ")


## Find target cells ##

if (is.null(target.cells))
  {
  if (all.cells)
    {
    target.cells <- as.numeric(substr(as.character(Cells$Cell), 6, nchar(as.character(Cells$Cell)))[!is.nan(rowMeans(Cells[, 2:ncol(Cells)], na.rm = T))])
    }
  else if (all.lineages)
    {
    cellsAtFirstImage <- find.Cells(InfoVertices$Images[1], Cells, Divisions, Vertices, InfoVertices)
    target.cells <- as.numeric(substr(cellsAtFirstImage, 6, nchar(cellsAtFirstImage)))
    }
  }


## Prepare PDF (if combined PDF for all cells) ##

if (combinedPDF)
  {
  if (all.cells) { outputame <- "Lineages/All cells" }
  else if (all.lineages) { outputName <- "Lineages/All lineages" }
  else { outputName <- paste("Lineages/Single cell__Cell", target.cells) }
  #pdf(paste(outputName, ".pdf", sep = ""), width = 10, height = 10*n.row/n.col)
  pdf(paste(outputName, ".pdf", sep = ""), width = 10, height = 10*n.row/(n.col+2)) # edit 2014/09/11
  #par(mar = c(1,0,0,0), mfrow = c(n.row, n.col))
  par(mar = c(1,0,0,0), mfrow = c(n.row, n.col+2)) # edit 2014/09/11
  }


## Loop on each target cell ##

for (target.cell in target.cells)
  {
  # Prepare PDF (if one PDF per target cell)
  if (!combinedPDF)
    {
    if (all.cells) { outputName <- paste("Lineages/All cells__Cell", target.cell) }
    else if (all.lineages) { outputName <- paste("Lineages/All lineages__Lineage ", which(target.cell == target.cells) , " (Cell ", target.cell, ")", sep = "") }
    else { outputName <- paste("Lineages/Single cell__Cell", target.cell) }
    #pdf(paste(outputName, ".pdf", sep = ""), width = 10, height = 10*n.row/n.col)
    pdf(paste(outputName, ".pdf", sep = ""), width = 10, height = 10*n.row/(n.col+2)) # edit 2014/09/11
    #par(mar = c(1,0,0,0), mfrow = c(n.row, n.col))
    par(mar = c(1,0,0,0), mfrow = c(n.row, n.col+2)) # edit 2014/09/11
    }

  # Plot cells
  draw.Lineage(target.cell, ProcessedImages, InfoVertices, Divisions, Cells, Vertices, leafShape, meanAngle, outputDir,
               n.col, n.row, main.plot, display.vertices, makeAVI, compr, fps)

  # Close PDF (if one PDF per target cell)
  if (!combinedPDF) { dev.off() }

  # Call the ImageJ macro
  IJ$runMacroFile("D:/Analyze_Lineages_called.txt");

  # Rename the montage and the movie (if any)
  if (!combinedPDF | !(all.cells | all.lineages))
    {
    file.rename("Lineages/Montage.png", paste(outputName, ".png", sep = ""))
    if (makeAVI) { file.rename("Lineages/Montage.avi", paste(outputName, ".avi", sep = "")) }
    }
  else
    {
    if (all.cells)
      {
      file.rename("Lineages/Montage.png", paste(outputName, "__Cell ", target.cell, ".png", sep = ""))
      if (makeAVI) { file.rename("Lineages/Montage.avi", paste(outputName, "__Cell ", target.cell, ".avi", sep = "")) }
      }
    else if (all.lineages)
      {
      file.rename("Lineages/Montage.png", paste(outputName, "__Lineage ", which(target.cell == target.cells) , " (Cell ", target.cell, ").png", sep = ""))
      if (makeAVI) { file.rename("Lineages/Montage.avi", paste(outputName, "__Lineage ", which(target.cell == target.cells) , " (Cell ", target.cell, ").avi", sep = "")) }
      }
    }
  }


## Close PDF (if combined PDF for all cells) ##

if (combinedPDF) { dev.off() }


## Create the csv file with information about division (if required) ##

if (create.csv)
  {
  arg <- list(InfoVertices = InfoVertices, Vertices = Vertices, Cells = Cells, Divisions = Divisions,
              meanAngle = meanAngle, alignToPetioleLaminaBoundary = F,
              optimal.interval = optimal.interval, leafShape = leafShape, CSV = T, File = csvFile, ini = ini)
  do.call(process.AllDiv, arg)
  }


## Delete output directory for cell coordinates ##

unlink(outputDir, recursive = TRUE)


## Save the number of images which have been plotted ##

write.table(length(target.cells), "D:/count_img.txt", row.names = F, col.names = F)