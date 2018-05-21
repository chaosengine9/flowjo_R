#to open a workspace
setwd("H:/Programming/R/FlowLink/TSK04_Trans_002/")
library(flowWorkspace)
wsFile <- file.path("4-Data Analysis (Excel, Batch Report, Prism)/180228 RC3-5/Panel 1.wsp")
ws <- openWorkspace(wsFile)

#different ways of looking at the workspace
ws
getSamples(ws)
getSampleGroups(ws)
sn <- "Specimen_001_A1_A01_001.fcs"
getKeywords(ws, sn)[1:5]

#parsing the data - going to have to make this detectable
gs_full<- parseWorkspace(ws, name=7, path="2-Raw Data/180228 RC3-5/TSK04_trans_002 RC3-5 Panel 1/Stain Plate");
gs_CD25FMO<- parseWorkspace(ws, name=2, path="2-Raw Data/180228 RC3-5/TSK04_trans_002 RC3-5 Panel 1/Stain Plate");
gs_CD38FMO<- parseWorkspace(ws, name=3, path="2-Raw Data/180228 RC3-5/TSK04_trans_002 RC3-5 Panel 1/Stain Plate");
class(gs)

#different ways of looking at the gating set - using gs[2] rather than whole gs as the gating hierarchy for the first file is different and messes it up
gs
sampleNames(gs)
plot(gs_full)
getNodes(gs, path = 1)
getNodes(gs, path = "auto")

#handling the seperate populations
nodelist <- getNodes(gs, path = "auto")
nodelist
#is.vector(nodelist) - good technique for checking things are what you think they are

node <- nodelist[8]
g <- getGate(gs, node)
g

ds_full<-getPopStats(gs_full)
ds_CD25FMO<-getPopStats(gs_CD25FMO)
ds_CD38FMO<-getPopStats(gs_CD38FMO)

#This would be where you would let users select populations to look at but for PoC..
tab_full <- subset(ds_full, Population == "aTreg cells (Fr. II)/CD25+")
tab_CD25FMO <- subset(ds_CD25FMO, Population == "aTreg cells (Fr. II)/CD25+")
tab_CD38FMO <- subset(ds_CD38FMO, Population == "aTreg cells (Fr. II)/CD25+")

#wrangling the dataframe
tab_full$group <- as.character(as.numeric(tab_full$group))
tab_CD25FMO$group <- as.character(as.numeric(tab_CD25FMO$group))
tab_CD38FMO$group <- as.character(as.numeric(tab_CD38FMO$group))
tab_full$group <- "Full Stain"
tab_CD25FMO$group <- "CD25 FMO"
tab_CD38FMO$group <- "CD38 FMO"
tab_merge <- rbind(tab_full, tab_CD25FMO, tab_CD38FMO)
library("plyr")
library("dplyr")
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
library(ggplot2)
library(hmisc)
plot1 <- ggplot(tab_merge, aes(x=group, y=frequency_of_Parent, )) + geom_dotplot(binaxis='y', binwidth = 1, stackdir='center')+ stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) + stat_summary(fun.y=mean, geom="point", color="red") + geom_bar(data=cdata,stat="identity", color="black", alpha=0.2) + theme_classic()+labs(title="Proportion of aTreg that are CD25+")
print(plot1)

plot1 <- ggplot(tab_merge[!excludevars], aes(x=group, y=frequency_of_Parent, )) + geom_point(stat="identity") +geom_errorbar(stat="identity") + geom_bar(stat="identity") + theme_classic() +labs(title="Frequency bar chart")
print(plot1)

ds_subset1_cols <- ds_fullset[, c("name", "Population", "Count", "ParentCount")]
ds_subset1 <- subset(ds_subset1_cols, Population == "aTreg cells (Fr. II)/CD25+") #two-step process can be halved! 
library(ggplot2)
plot1 <- ggplot(ds_subset1, aes(x=name, y=Count, )) + geom_bar(stat="identity") + labs(title="Frequency bar chart")
print(plot1)

#plots - for more see http://bioconductor.org/packages/release/bioc/vignettes/flowWorkspace/inst/doc/plotGate.html
plotGate(gs_full, "aTreg cells (Fr. II)", xbin = 0, main = "Graph title")

plotGate(gs[[2]], gpar = list(nrow = 5))
gh <- gs[[4]]
plot(gh)
getProp(gh,node)
getPopStats(gh)

#data wrangling
fs <- getData(gs)
class(fs)

aTreg <- getData(gs, node)


nrow(aTreg[[7]])

getData(gs);
table(getIndices(gh,node))
