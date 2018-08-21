
# Choices (ASSIGN INPUT FROM GUI) -----------------------------------------------------------------
##Overlay expression? y/n
expressionToggle <- as.character("y")

##Which parameter?
expressionParameter <- as.character("FMOSubtractSet.Adjusted_CD25_PE_gMFI")

##Stack populations based on parent? y/n
stackToggle <- as.character("y")

##Group bars based on indication? y/n
indicationToggle <- as.character("n")


# Libraries ---------------------------------------------------------------
#library(flowWorkspace)
library(plyr)
library(dplyr)
library(ggplot2)
#library(ggthemes)
#library(extrafont)
library(scales)
#library(hmisc)
#library(cytofkit)
#library(scales)
library(wesanderson)
library(RColorBrewer)
library(data.table)


# Functions ---------------------------------------------------------------
squish_trans <- function(from, to, factor) {
  
  trans <- function(x) {
    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squished", trans, inv))
}

statSetting <- function(x, y, col){
  z <- ddply(x, y, .fun = function(xx){
    c(mean= mean(xx[,col],na.rm=TRUE),
      N = length(xx[,col]), sd = sd(xx[,col],na.rm=TRUE),se = sd(xx[,col],na.rm=TRUE) / sqrt(length(xx[,col]))) })
  return(z)
}

# Generate fullTable -----------------------------------------------
#Load files
setwd("H:/Programming/R/FlowLink/TSK04_Trans_002/") ##Change to Linux server location
fileList <- c("RCC_001_DTC_panel2.rds","RCC_002_DTC_panel2.rds","RCC_003_DTC_panel2.rds") ##ASSIGN INPUT FROM GUI
dat_list = lapply(fileList, function (x) data.table(readRDS(x)))
fullTable = rbindlist(dat_list, fill = TRUE)

#remove expression that causes trouble -> âˆ’
fullTable <- as.data.frame(lapply(fullTable, function(x) {
  gsub("âˆ’", "_", x)}), stringsAsFactors = FALSE)

#Make additional columns numeric
numericColumns <- names(fullTable)[!names(fullTable) %in%  c("name", "Path", "Population", "Parent", "Donor", "sample", "donorList", "panelList", "stainList")]
for(numericColumn in numericColumns) {
  fullTable[[numericColumn]] <- as.numeric(fullTable[[numericColumn]])
}

#Add extra statistic
fullTable <- mutate(fullTable,frequency_of_Parent = (Count/ParentCount)*100)



# Subset fulltable and subtract FMOs --------------------------------------------------------
##Get subset of data
plotSet <- subset(fullTable, Population == "B cells" | Population == "CD14-CD3-" |Population == "Myeloid DC" | Population == "Plasmacytoid DC" | Population == "Dendritic cells" |Population == "Lin-" | Population == "CD56 high" | Population == "CD56 low" | Population == "NK cells" | Population == "Monocytes" | Population == "Classical Monocytes" | Population == "Intermediate Monocytes" | Population == "Non-classical monocytes" | Population == "T cells")
plotSet <- mutate(plotSet,ID = paste(donorList, Population, sep = '_'))
plotSet <- mutate(plotSet, Indication = strsplit(ID,'_')[[1]][1])


##subtract FMOs
FMOSubtractSet <- plotSet[c("ID","name","Path","Population","Parent","Count","ParentCount", "donorList", "panelList", "stainList", "frequency_of_Parent", "gMFI_Comp.PE.A")]
FMOSubtractSet <- subset(FMOSubtractSet, stainList == "Full stain"|stainList == "CD25 FMO")
FMOSubtractSet <- reshape(FMOSubtractSet, idvar = "ID", timevar = "stainList", direction = "wide")

pathSplits<-data.frame(do.call('rbind', strsplit(as.character(FMOSubtractSet$'Path.Full stain'),'/',fixed=TRUE)))
for(i in 1:ncol(pathSplits)){
  if(length(unique(pathSplits[,i]))>1){
    #create the extra column in FMOSubtractSet that has the right FMO value to subtract etc
    FMOSubtractSet$bello <- pathSplits[,i]
    dfFMO <- merge(data.frame(FMOSubtractSet$ID, FMOSubtractSet$`Population.Full stain`,FMOSubtractSet$`gMFI_Comp.PE.A.CD25 FMO`),data.frame(unique(pathSplits[,i])), by.x = "FMOSubtractSet..Population.Full.stain.", by.y = "unique.pathSplits...i..")
    FMOSubtractSet$bello <- paste0(FMOSubtractSet$`donorList.Full stain`,"_",FMOSubtractSet$bello) 
    FMOSubtractSet <- merge(dfFMO, FMOSubtractSet, by.x = "FMOSubtractSet.ID", by.y = "bello")
    FMOSubtractSet <- mutate(FMOSubtractSet, Adjusted_CD25_PE_gMFI = ifelse(`gMFI_Comp.PE.A.Full stain` - `FMOSubtractSet..gMFI_Comp.PE.A.CD25.FMO.`>0,`gMFI_Comp.PE.A.Full stain` - `FMOSubtractSet..gMFI_Comp.PE.A.CD25.FMO.`,0))
    plotSet <- subset(plotSet, stainList == "Full stain")
    plotSet <- merge(plotSet,data.frame(FMOSubtractSet$ID,FMOSubtractSet$Adjusted_CD25_PE_gMFI),by.x = "ID",by.y = "FMOSubtractSet.ID")
    break
  }
}

plotSet <- subset(plotSet, Population == "B cells" | Population == "Myeloid DC" | Population == "Plasmacytoid DC" | Population == "Lin-" | Population == "CD56 high" | Population == "CD56 low" | Population == "Classical Monocytes" | Population == "Intermediate Monocytes" | Population == "Non-classical monocytes" | Population == "T cells") ##ASSIGN INPUT FROM GUI
#plotSet <- subset(plotSet, donorList == "RCC_001_DTC") ##ASSIGN INPUT FROM GUI

##Dataframe per donor and dynamically figure out x-axis labels
tableList <- list()
for (i in 1:length(unique(FMOSubtractSet$`donorList.Full stain`))){
  plotTableName <- unique(FMOSubtractSet$`donorList.Full stain`)[i]
  tableListItem <- subset(plotSet, donorList == unique(FMOSubtractSet$`donorList.Full stain`)[i])
  change <- tableListItem[duplicated(tableListItem$Parent),]$Parent
  tableListItem <- within(tableListItem, { xLabels = ifelse(Parent %in% change, Parent, Population) })
  tableList[[plotTableName]] <- tableListItem
}

# INDIVIDUAL PLOTS --------------------------------------------
##plot each donor
myPlots <- list()
columnNames <- names(tableList[[1]])
for (i in 1:length(tableList)){
  local({
  i <- i
  tempSet <- as.data.frame(tableList[i])
  colnames(tempSet) <- columnNames
  stackToggle <- assign(stackToggle, ifelse(stackToggle == "y", "xLabels", "Population"))
  if(expressionToggle=="y"){
    plot <- ggplot(tempSet, aes(x=get(stackToggle), y=Count, fill=get(expressionParameter))) + geom_bar(stat="identity") +  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ labs(fill=paste0(expressionParameter), title=paste0("proportions and CD25 receptor density of immune \n cell populations in ", names(tableList)[i])) + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank()) + scale_fill_gradient(low = "white", high = "darkred", limits=c(0,as.numeric(paste0(max(as.numeric(lapply(tableList, function(x){max(x[[expressionParameter]])})))))))
  } else if(expressionToggle=="n"){
    plot <- ggplot(tempSet, aes(x=get(stackToggle), y=Count, fill=Population)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title= paste0("Proportions of immune cell populations in ", names(tableList)[i])) + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank()) +scale_fill_brewer(palette="Set3")
  }
  myPlots[[i]] <<- plot
  print(myPlots[i])
  })
}


# POOL PLOTS --------------------------------------------------------------
plotSet <- rbindlist(tableList)

##get statistics for the pooled donors
y_variable <- as.character("Count")
statSet <- statSetting(plotSet, "Population", y_variable)
statSet[[y_variable]]= with(statSet,mean)
statSet <- merge(statSet,unique.data.frame(data.frame(plotSet$Population,plotSet$xLabels)), by.x = "Population", by.y = "plotSet.Population")
names(statSet)[7] <- "xLabels"


stackToggle <- assign(stackToggle, ifelse(stackToggle == "y", "xLabels", "Population"))
if(indicationToggle=="y"){
  poolPlot <- ggplot(plotSet, aes(x=get(stackToggle), y=get(y_variable), fill=Population)) + geom_bar(data=statSet,stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +scale_fill_brewer(palette="Set3") + geom_dotplot(binaxis='y', binwidth = 1, stackdir='center') + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) + labs(y = paste0(y_variable), title= paste0("Proportions of immune cell populations in ")) + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank()) +facet_wrap(~Indication)
} else if(indicationToggle=="n"){
  poolPlot <- ggplot(plotSet, aes(x=get(stackToggle), y=get(y_variable), fill=Population)) + geom_bar(data=statSet,stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +scale_fill_brewer(palette="Set3") + geom_dotplot(binaxis='y', binwidth = 1, stackdir='center') + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) + labs(y = paste0(y_variable), title= paste0("Proportions of immune cell populations in ")) + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank())
}
print(poolPlot)




# Test code ---------------------------------------------------------------

#FACS plots
#Fancy FACS plots - these don't work from the RDS files, only the GatingSets
#for more see http://bioconductor.org/packages/release/bioc/vignettes/flowWorkspace/inst/doc/plotGate.html
plotGate(gs_full, "aTreg cells (Fr. II)", xbin = 0, main = "Graph title")

#t-SNE
## install
##source("https://bioconductor.org/biocLite.R")
##biocLite("cytofkit")

## run new analysis
cytofkit_GUI()

## open old analysis (or recently finished analysis)
cytofkitShinyAPP()
