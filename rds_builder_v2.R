#open libraries
library(flowWorkspace)
library(dplyr)
library(ggplot2)
library(cytofkit)
library(data.table)

# Functions ---------------------------------------------------------------

getStats <- function(x, ...)UseMethod("getStats")

getStats.GatingSetList <- function(x, ...){
  getStats.GatingSet(x, ...)
}

#' @export
#' @rdname getStats
getStats.GatingSet <- function(x, ...){
  res <-  lapply(x, function(gh){
    getStats(gh, ...)
    
  })
  rbindlist(res, idcol = "name")
}

getStats.GatingHierarchy <- function(x, nodes = NULL, type = "count", inverse.transform = FALSE, stats.fun.arg = list(), ...){
  gh <- x
  if(is.null(nodes))
    nodes <- getNodes(gh, ...)
  res <- sapply(nodes, function(node){
    if(is.character(type))
    {
      type <- match.arg(type, c("count", "percent"))
      if(type == "count")
      {
        res <- getTotal(gh, node)
        names(res) <- "count"
      }else if(type == "percent")
      {
        res <- getProp(gh, node)
        names(res) <- "percent"
      }else
        stop("unsupported stats type: ", type)
    }else{
      fr <- getData(gh, y = node)
      if(inverse.transform)
      {
        trans <- getTransformations(gh, inverse = TRUE)
        if(length(trans)==0)
          stop("No inverse transformation is found from the GatingSet!")
        trans <- transformList(names(trans), trans)
        fr <- transform(fr, trans)
      }
      thisCall <- quote(type(fr))
      thisCall <- as.call(c(as.list(thisCall), stats.fun.arg))
      
      res <- eval(thisCall)
    }
    
    as.data.table(t(res))
  }, simplify = FALSE)
  rbindlist(res, idcol = "Path")
  
}
pop.MFI <- function(fr){
  pd <- pData(parameters(fr))
  pd <- data.table(pd)
  #pd <- pd[!is.na(desc), ] Mike Jiang said to delete that bit... problem is that it has highlighted the problem, it can't find names for the parameters
  chnls  <- pd[, name]
  markers <- pd[, desc]
  
  res <- colMedians(exprs(fr)[, chnls, drop = FALSE])
  names(res) <- markers
  res
}

multiply <- function(x){
  Reduce("*",x)
  
}

pop.quantiles <- function(fr){ 
  chnls <- colnames(fr) 
  res <- matrixStats::colQuantiles(exprs(fr), probs = 0.50) 
  names(res) <- chnls 
  res 
} 

pop.medians <- function(fr){ 
  chnls <- paste("MedianFI", colnames(fr), sep = "_")
  res <- matrixStats::colMedians(exprs(fr)) 
  names(res) <- chnls 
  res 
} 

pop.means <- function(fr){ 
  chnls <- colnames(fr) 
  res <- matrixStats::colMeans2(exprs(fr)) 
  names(res) <- chnls 
  res 
} 

geometric.means <- function(fr){
  chnls <- paste("gMFI", colnames(fr), sep = "_")
  exprs(fr)[exprs(fr) <= 0] <- NA #excludes all negative values as is necessary for the next step but this is causing discrepancies with FlowJo
  res <- apply(exprs(fr), 2, function(x)exp(mean(log(x),na.rm = TRUE)))
  names(res) <- chnls
  res
}



# Body --------------------------------------------------------------------


#Open a workspace
setwd("H:/Programming/R/FlowLink/TSK04_Trans_002/")
sFile <- file.path("4-Data Analysis (Excel, Batch Report, Prism)/180228 RC3-5/Panel 2.wsp")
ws <- NULL
ws <- openWorkspace(wsFile)


#Parse full stain and FMO data for each donor into it's own RDS file
##User says which FCS file is which
###For now this is hardcoded below:
fileList <- c("Specimen_001_F1_F01_001.fcs", "Specimen_001_F2_F02_002.fcs", "Specimen_001_F3_F03_003.fcs", "Specimen_001_G1_G01_005.fcs", "Specimen_001_G2_G02_006.fcs", "Specimen_001_G3_G03_007.fcs", "Specimen_001_H1_H01_009.fcs", "Specimen_001_H2_H02_010.fcs", "Specimen_001_H3_H03_011.fcs")
donorList <- c("RCC_001_DTC", "RCC_001_DTC", "RCC_001_DTC", "RCC_002_DTC", "RCC_002_DTC", "RCC_002_DTC", "RCC_003_DTC", "RCC_003_DTC", "RCC_003_DTC") #need a proper donor naming system and to have that in the 
panelList <- c("panel2", "panel2", "panel2", "panel2", "panel2", "panel2", "panel2", "panel2", "panel2")
stainList <- c("Full stain", "CD25 FMO", "CD38 FMO", "Full stain", "CD25 FMO", "CD38 FMO", "Full stain", "CD25 FMO", "CD38 FMO")
FCSfiles <- data.frame(fileList,donorList,panelList,stainList)
FCSfiles <- data.frame(lapply(FCSfiles, as.character), stringsAsFactors=FALSE)
subsetList <- list()
allSamples <- "All Samples"

for (i in 1:nrow(getFJWSubsetIndices(ws, key=NULL, value = NULL, 1))){
  temp_gs <- parseWorkspace(ws, name=allSamples, path="2-Raw Data/180228 RC3-5/TSK04_trans_002 RC3-5 Panel 2/Stain Plate", subset=i)
  subsetForList <- getPopStats(temp_gs)
  subsetForList$Population <- getPopStats(temp_gs)$Population
  subsetForList$Path <- getPopStats(temp_gs, path="full")$Population
  subsetForList <- merge(subsetForList,getStats(temp_gs, getNodes(temp_gs), type = pop.medians, inverse.transform = TRUE),by=c("name","Path"))
  subsetForList <- merge(subsetForList,getStats(temp_gs, getNodes(temp_gs), type = geometric.means, inverse.transform = TRUE),by=c("name","Path")) 
  subsetForList <- as.data.frame(lapply(subsetForList, function(x) {
    gsub("\\.fcs..*", ".fcs", x)}), stringsAsFactors = FALSE)
  subsetForList <- merge(subsetForList, FCSfiles, by.x = "name", by.y = "fileList")
  subsetNameInList <- subsetForList$donorList[1]
  if ("RCC_001_DTC" %in% names(subsetList)){
    subsetList[[subsetNameInList]] <- rbind(subsetList[[subsetNameInList]],subsetForList)
  }else{
    subsetList[[subsetNameInList]] <- subsetForList
  }
}
lapply(subsetList,function(x){
  filepath_RDS = paste0(x$donorList[1],"_", x$panelList[1], ".rds")
  saveRDS(x,file=filepath_RDS)
})


### need to figure out how to automate the receptor density calculations.. 
### Simplest way would be to input the data for the standard curve in the upload R shiny
### the next step would be for R to generate the standard curve and then an extra column of receptor density from it


