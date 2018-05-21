#open libraries
library(flowWorkspace)
library("plyr")
library("dplyr")
library(ggplot2)
library(ggthemes)
library(extrafont)
library(scales)
library(hmisc)
library(cytofkit)

#Open a workspace
setwd("H:/Programming/R/FlowLink/TSK04_Trans_002/")
wsFile <- file.path("4-Data Analysis (Excel, Batch Report, Prism)/180228 RC3-5/Panel 1.wsp")
ws <- openWorkspace(wsFile)

#parsing the data
#name specifies the sample group. perhaps this should be user controlled 
#e.g. "select the full stain group" then "select the CD25 FMO group" etc...
gs_full<- parseWorkspace(ws, name=7, path="2-Raw Data/180228 RC3-5/TSK04_trans_002 RC3-5 Panel 1/Stain Plate");
gs_CD25FMO<- parseWorkspace(ws, name=2, path="2-Raw Data/180228 RC3-5/TSK04_trans_002 RC3-5 Panel 1/Stain Plate");
gs_CD38FMO<- parseWorkspace(ws, name=3, path="2-Raw Data/180228 RC3-5/TSK04_trans_002 RC3-5 Panel 1/Stain Plate");

#Test section - need a for loop to get ds for each subset
subsetList <- list()
for (j in 1:rapply(getSampleGroups(ws),function(x)length(unique(x)))[[1]]){ #what a convoluted way of getting the number 7!!
  subsetList <- list()
  for (i in 1:nrow(getFJWSubsetIndices(ws, key=NULL, value = NULL, j))){
    subsetNameInList <- paste0("gs_",i)
    subsetForList <- getPopStats(parseWorkspace(ws, name=j, path="2-Raw Data/180228 RC3-5/TSK04_trans_002 RC3-5 Panel 1/Stain Plate", subset=i))
    subsetList[[subsetNameInList]] <- subsetForList
  }
  assign(paste0("group",j),subsetList)
  subset <- NULL
}

#Test section - need a for loop to get ds for each subset
subsetList <- list()
for (j in 1:rapply(getSampleGroups(ws),function(x)length(unique(x)))[[1]]){ #what a convoluted way of getting the number 7!!
  for (i in 1:nrow(getFJWSubsetIndices(ws, key=NULL, value = NULL, j))){
    subsetNameInList <- paste0("gs_",i)
    subsetForList <- getPopStats(parseWorkspace(ws, name=j, path="2-Raw Data/180228 RC3-5/TSK04_trans_002 RC3-5 Panel 1/Stain Plate", subset=i))
    subsetList[[subsetNameInList]] <- subsetForList
  }
}

#send the gs classes to the server as individual .RDS files for each patient and run the rest of the script from there
#Going to have to think about how to split and join samples so that they sit on the server sensibly i.e. one file for wseach patient that contains full stain and FMOs

#making data objects
ds_full<-getPopStats(gs_full)
ds_CD25FMO<-getPopStats(gs_CD25FMO)
ds_CD38FMO<-getPopStats(gs_CD38FMO)

#Making select data tables
tab_full <- subset(ds_full, Population == "aTreg cells (Fr. II)/CD25+")
tab_CD25FMO <- subset(ds_CD25FMO, Population == "aTreg cells (Fr. II)/CD25+")
tab_CD38FMO <- subset(ds_CD38FMO, Population == "aTreg cells (Fr. II)/CD25+")

#Ensuring control and test groups are distinguishable and merging the tables
tab_full$group <- as.character(as.numeric(tab_full$group))
tab_CD25FMO$group <- as.character(as.numeric(tab_CD25FMO$group))
tab_CD38FMO$group <- as.character(as.numeric(tab_CD38FMO$group))
tab_full$group <- "Full Stain"
tab_CD25FMO$group <- "CD25 FMO"
tab_CD38FMO$group <- "CD38 FMO"
tab_merge <- rbind(tab_full, tab_CD25FMO, tab_CD38FMO)


#Getting extra statistics
tab_merge <- mutate(tab_merge,frequency_of_Parent = (Count/ParentCount)*100)

#summary table
cdata <- ddply(tab_merge, c("group"), summarise,
               N    = length(frequency_of_Parent),
               mean = mean(frequency_of_Parent),
               frequency_of_Parent = mean(frequency_of_Parent),
               sd   = sd(frequency_of_Parent),
               se   = sd / sqrt(N)
)

#plot the bar graph
plot1 <- ggplot(tab_merge, aes(x=group, y=frequency_of_Parent, )) + geom_dotplot(binaxis='y', binwidth = 1, stackdir='center')+ stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) + stat_summary(fun.y=mean, geom="point", color="red") + geom_bar(data=cdata,stat="identity", color="black", alpha=0.2) + theme_classic()+labs(title="Proportion of aTreg that are CD25+")
print(plot1)

#plot the Stacked columns
##plot2 <- ggplot(tab_merge, aes(x=group, y=frequency_of_Parent, )) + geom_dotplot(binaxis='y', binwidth = 1, stackdir='center')+ stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) + stat_summary(fun.y=mean, geom="point", color="red") + geom_bar(data=cdata,stat="identity", color="black", alpha=0.2) + theme_classic()+labs(title="Proportion of aTreg that are CD25+")
print(plot1)


#Fancy FACS plots 
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
