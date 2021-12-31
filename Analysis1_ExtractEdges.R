# Attempt to get a network from the kgml!

# library(KEGGREST)
# k1 <- keggGet(dbentries="ko01100"
#           , option = c("kgml"))
# setwd("/Users/arr47/Documents/Dev")
# save(k1, file = "KEGGtest.txt")
# writeLines(k1, con="KEGGtest.txt")

setwd("/Users/arr47/Documents/Dev")
fileName <- "KEGGtest.txt"
conn <- file(fileName,open="r")
linn <-readLines(conn)
length(linn)
linn[1]

c1 <- linn[grep(pattern = "cpd:", x = linn)]
head(c1)
tail(c1)

k1 <- linn[grep(pattern = "ko:", x = linn)]
k1[1:10]

library(stringr)
l1 <- str_extract_all(string = k1, pattern = "ko:K[0-9]+")
l2 <- str_extract_all(string = k1, pattern = "rn:R[0-9]+")

length(l1)
length(l2)
kom <- data.frame(KO = rep(NA, count), RN = rep(NA, count))
count = 1
for (i in 1:length(l1)){
  nk <- length(l1[[i]])
  nr <- length(l2[[i]])
  for (j in 1:nk){
    for (k in 1:nr){
      print(c( l1[[i]][j], l2[[i]][k]))
      if(length( l2[[i]][k] ) > 0){
        kom[count, "KO"] <- l1[[i]][j]
        kom[count, "RN"] <- l2[[i]][k]
        count = count + 1
      }
    }
  }
}
count
head(kom) # these are all the reactions so the question is what are the substrates and products for these
kom <- kom[!is.na(kom$RN), ]
nrow(kom)

r1 <- linn[grep(pattern = "<reaction id=", x = linn)]
s1 <- linn[grep(pattern = "<substrate id=", x = linn)]
p1 <- linn[grep(pattern = "<product id=", x = linn)]

length(linn)
linn[40008]
# make it 20,000
m1 <- matrix(data=NA, nrow = length(linn), ncol=7)
for (i in 1:length(linn)){
  #print(linn[i])
  if(grepl(pattern = "<reaction id=", x = linn[i])){
    m1[i,1] <- linn[i]
    if(!grepl(pattern = "</reaction", x = linn[i+1])){
      m1[i,2] <- linn[i+1]
    }
    
    if(!grepl(pattern = "</reaction", x = linn[i+2])){
      m1[i,3] <- linn[i+2]
    }
    
    if(!grepl(pattern = "</reaction", x = linn[i+3])){
      m1[i,4] <- linn[i+3]
    }
    
    if(!grepl(pattern = "</reaction", x = linn[i+4])){
      m1[i,5] <- linn[i+4]
    }
    
    if(!grepl(pattern = "</reaction", x = linn[i+5])){
      m1[i,6] <- linn[i+5]
    }
    
    if(!grepl(pattern = "</reaction", x = linn[i+6])){
      m1[i,7] <- linn[i+6]
    }
  }
}
head(m1)
d1 <- as.data.frame(m1)
d1 <- d1[!is.na(d1$V1), ]
head(d1)
summary(d1)



head(d1)

# try to parse the substrate-product pairs
d1[is.na(d1$V4), c("V4", "V5", "V6", "V7")] <- "nada"
d1[is.na(d1$V5), c("V5", "V6", "V7")] <- "nada"
d1[is.na(d1$V6), c("V6" ,"V7")] <- "nada"
setwd("/Users/arr47/Documents/Dev")
write.table(d1, "IntermediateReactions.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

head(d1)
d1$agg <- paste(d1$V1, d1$V2, d1$V3, d1$V4, d1$V5, d1$V6, d1$V7, sep="_")

library(stringr)
str_extract_all(string = d1$agg, pattern = "substrate id=[^/]+")

s1 <- str_extract_all(string = d1$agg, pattern = "substrate id=[^/]+")
p1 <- str_extract_all(string = d1$agg, pattern = "product id=[^/]+")
r1 <- str_extract_all(string = d1$agg, pattern = "rn:R[0-9]+")
length(s1)
length(p1)
length(r1)
r1
rsp <- data.frame(RN = rep(NA, count), Sub = rep(NA, count), Prod = rep(NA, count))
count <- 1
for (i in 1:length(r1)){
  rtmp <- r1[[i]]
  stmp <- s1[[i]]
  ptmp <- p1[[i]]
  nr <- length(rtmp)
  ns <- length(stmp)
  np <- length(ptmp)
  if(length(stmp) > 0 & length(ptmp) > 0 & length(rtmp) > 0 ){
  for (j in 1:nr){
    for (k in 1:ns){
      for(l in 1:np){
        rsp[count, "RN"] <- rtmp[j]
        rsp[count, "Sub"] <- stmp[k]
        rsp[count, "Prod"] <- ptmp[l]
        count = count + 1
      }
    }
  }
    }
}
count

head(rsp)

rsp$SubClean <- str_extract(rsp$Sub, "cpd:C[0-9]+")
rsp$ProdClean <- str_extract(rsp$Prod, "cpd:C[0-9]+")
nrow(rsp)


# Now combine the data into a network
# reactions are connected when the product of one is the substrate for another
head(kom)
head(rsp)

# do an example
prod <- "cpd:C00158"
#when is this a substrate?
rsp[rsp$SubClean %in% prod, ]


# for each reaction what are its subsequent reactions? => create a adj list from that
arn <- data.frame(FirstRn = rep(NA, count), SecondRn = rep(NA, count))
count <- 1
for (i in 1:nrow(rsp)){
  rtmp <- rsp[i, "RN"]
  ptmp <- rsp[i, "ProdClean"]
  print( unique( rsp[rsp$SubClean %in% ptmp, ]$RN ) )
  adjrxns <- unique( rsp[rsp$SubClean %in% ptmp, ]$RN )
  if (length(adjrxns) > 0){
    for (j in 1:length(adjrxns)){
      arn[count, "FirstRn"] <- rtmp
      arn[count, "SecondRn"] <- adjrxns[j]
      count <- count + 1
    }
  }
}
count

head(arn)

# now expand this out acccording to the KO
head(kom)
arn2 <- merge(arn, kom, by.x = "FirstRn", by.y = "RN", all.x=FALSE, all.y=FALSE)
head(arn2)
colnames(arn2)[3] <- "FirstKO"

arn3 <- arn2[!duplicated(arn2), ]


arn4 <- merge(arn3, kom, by.x = "SecondRn", by.y = "RN", all.x=FALSE, all.y=FALSE)
arn4 <- arn4[!duplicated(arn4), ]

# still half a million!!!
head(arn4)
colnames(arn4)[4] <- "SecondKO"

arn5 <- arn4[ ,c("FirstKO", "SecondKO")]
arn5 <- arn5[!duplicated(arn5), ]
arn5 <- arn5[arn5$FirstKO != arn5$SecondKO, ]
head(arn5)

write.table(arn5, "Test_KEGG_Network.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

setwd("/Users/arr47/Documents/Dev/KEGG_Dev")
write.table(arn5, "KEGG_Network.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


setwd("/Users/arr47/Documents/Databases")
write.table(arn5, "KEGG_Network.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


# length(l1)
# length(l2)
# kom <- data.frame(KO = rep(NA, count), RN = rep(NA, count))
# count = 1
# for (i in 1:length(l1)){
#   nk <- length(l1[[i]])
#   nr <- length(l2[[i]])
#   for (j in 1:nk){
#     for (k in 1:nr){
#       print(c( l1[[i]][j], l2[[i]][k]))
#       if(length( l2[[i]][k] ) > 0){
#         kom[count, "KO"] <- l1[[i]][j]
#         kom[count, "RN"] <- l2[[i]][k]
#         count = count + 1
#       }
#     }
#   }
# }
# 
# 
# head(kom)
# kom$X1 <- NA
# kom$X2 <- NA
# kom$X3 <- NA
# kom$X4 <- NA
# kom$X5 <- NA
# 
# i <- 1
# for (i in 1:nrow(kom)){
#   rtmp <- kom[i, "RN"]
#   tmp2 <- d1[grep(pattern = rtmp, x=d1$V1), ]
#   if (nrow(tmp2) == 1){
#     if (is.na(tmp2$V6)){
#       kom[i, "X1"] <- tmp2$V2
#       kom[i, "X2"] <- tmp2$V3
#       kom[i, "X3"] <- tmp2$V4
#       kom[i, "X4"] <- tmp2$V5
#     }
#     else if (is.na(tmp2$V5)){
#       kom[i, "X1"] <- tmp2$V2
#       kom[i, "X2"] <- tmp2$V3
#       kom[i, "X3"] <- tmp2$V4
#     }
#     
#     else if (is.na(tmp2$V4)){
#       kom[i, "X1"] <- tmp2$V2
#       kom[i, "X2"] <- tmp2$V3
#     }
#   }
#   print(i)
# }
# head(kom)
# 
# 
# 
# 
# 
# 
# 
# 
close(conn)
