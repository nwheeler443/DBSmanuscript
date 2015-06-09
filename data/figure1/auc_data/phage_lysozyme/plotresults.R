library(ROCR)


dbs <- read.table("40.scan")
dbs <- dbs[order(dbs$V4),]
dbs$dbs <- dbs[dbs$V4=="wild_type",8] - dbs[,8]
dbs <- dbs[-(grep("wild_type", dbs$V4)),]
labels <- rep(0, nrow(dbs))
labels[grep("-", dbs$V4)] = 1
preddbs <- prediction(dbs$dbs, labels)
perfdbs <- performance(preddbs, "auc")

pfam <- read.table("scanproc2")
pfam <- pfam[order(pfam$V1),]
pfam <- pfam[-(grep("wild_type", pfam$V1)),]
labels <- rep(0, nrow(pfam))
labels[grep("-", pfam$V1)] = 1
predpfam <- prediction(pfam$V2, labels)
perfpfam <- performance(predpfam, "auc")

provean <- read.delim("provean.tsv", header=FALSE)
provean <- provean[order(provean$V1),]
provpred <- prediction(provean[,2]*(-1), labels)
provperf <- performance(provpred, "auc")

polyphen <- read.delim("polyphen.txt")
polyphen$name <- do.call(paste0, polyphen[c(3,2,4)])
polyphen <- polyphen[order(polyphen$name),]
polypred <- prediction(polyphen$pph2_prob, labels)
polyperf <- performance(polypred, "auc")

