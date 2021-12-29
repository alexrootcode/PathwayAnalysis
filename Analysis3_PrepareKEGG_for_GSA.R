# Script to convert the KEGG flat file to use with the GSA package

rm(list = ls())
setwd("/Users/arr47/Documents/Dev/KEGG_Dev")

k1 <- read.table("KEGG_PathwayHierarchy_Flatfile.txt", sep="\t", header=TRUE)
head(k1)

# summarize the tree
# hierarchy 1
h1 <- data.frame(table(k1$HighestPathway))
h1[order(h1$Freq), ]

# hierarchy 2
h2 <- data.frame(table(k1$HigherPathway))
h2[order(h2$Freq), ]

# hierarchy 3
h3 <- data.frame(table(k1$Pathway))
h3[order(h3$Freq, decreasing = TRUE), ][1:10, ]

# how many genesets are there?
length(unique(k1$Pathway))
ps <- unique(k1$Pathway) # get the pathways
gs <- vector("list", length = length(ps))
i <- 1

paste(k1[k1$Pathway %in% ps[i], ]$KO, sep="")

for (i in 1:length(ps)){
  gs[[i]] <- k1[k1$Pathway %in% ps[i], ]$KO
  #names(gs[[i]]) <- ps[i]
}
gs[[452]]
unique(gs[452])
gs.names <- ps


# create a dummy matrix where clearly Photosynthesis is over-expressed
# and clearly nitrogen metabolism is under-expressed
head(k1)
kphoto <- k1[grep("00195 Photosynthesis", k1$Pathway), ]$KO
knitro <- k1[grep("Nitrogen", k1$Pathway), ]$KO

length(unique(k1$KO))

# Create a toy data set with 3 cases and 3 controls
toydf <- data.frame(KO = unique(k1$KO), case1 = rnorm(n = length(unique(k1$KO)), mean = 0, sd = 0.01),
                    case2 = rnorm(n = length(unique(k1$KO)), mean = 0, sd = 0.001),
                  case3 = rnorm(n = length(unique(k1$KO)), mean = 0, sd = 0.001),
                    ctl1 = rnorm(n = length(unique(k1$KO)), mean = 0, sd = 0.001),
                    ctl2 = rnorm(n = length(unique(k1$KO)), mean = 0, sd = 0.001),
                    ctl3 = rnorm(n = length(unique(k1$KO)), mean = 0, sd = 0.001) )
head(toydf)

# For the KOs in the photo and nitro pathways, artificially bump them out or down
toydf[toydf$KO %in% kphoto, ]$case1 <- toydf[toydf$KO %in% kphoto, ]$case1 + 5
toydf[toydf$KO %in% kphoto, ]$case2 <- toydf[toydf$KO %in% kphoto, ]$case2 + 5
toydf[toydf$KO %in% kphoto, ]$case3 <- toydf[toydf$KO %in% kphoto, ]$case3 + 5

toydf[toydf$KO %in% knitro, ]$case1 <- toydf[toydf$KO %in% knitro, ]$case1 - 5
toydf[toydf$KO %in% knitro, ]$case2 <- toydf[toydf$KO %in% knitro, ]$case2 - 5
toydf[toydf$KO %in% knitro, ]$case3 <- toydf[toydf$KO %in% knitro, ]$case3 - 5

# ok, so now there is some signal, although not much
library(GSA)

head(toydf)
toym <- as.matrix(toydf[ ,c(2:7)])
rownames(toym) <- unique(k1$KO)
head(toym)

# Runs the test
GSA.obj<-GSA(x = toym, 
             y = c(2,2,2,1,1,1), 
             genenames=unique(k1$KO), 
             genesets=gs,  
             resp.type="Two class unpaired", nperms=100)
GSA.obj

GSA.listsets(GSA.obj, geneset.names=gs.names, FDRcut=.5)






