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

#parsing the data
gs <- parseWorkspace(ws, name=1, path="2-Raw Data/180228 RC3-5/TSK04_trans_002 RC3-5 Panel 1/Stain Plate");
class(gs)

#different ways of looking at the gating set - using gs[2] rather than whole gs as the gating hierarchy for the first file is different and messes it up
gs
sampleNames(gs)
plot(gs)
getNodes(gs, path = 1)
getNodes(gs, path = "auto")

#handling the seperate populations
nodelist <- getNodes(gs, path = "auto")
nodelist
#is.vector(nodelist) - good technique for checking things are what you think they are

node <- nodelist[8]
g <- getGate(gs, node)
g

ds_fullset<-getPopStats(gs)
ds_subset1_cols <- ds_fullset[, c("name", "Population", "Count", "ParentCount")]
ds_subset1 <- subset(ds_subset1_cols, Population == "aTreg cells (Fr. II)/CD25+") #two-step process can be halved! 
library(ggplot2)
plot1 <- ggplot(ds_subset1, aes(x=name, y=ParentCount, )) + geom_bar(stat="identity") + labs(title="Frequency bar chart")
print(plot1)

#plots - for more see http://bioconductor.org/packages/release/bioc/vignettes/flowWorkspace/inst/doc/plotGate.html
plotGate(gs[c(1,4,7)], "aTreg cells (Fr. II)", xbin = 0, main = "Graph title")

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
