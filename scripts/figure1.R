library("ROCR")

lysozyme <- read.delim("lysozymeAUC.txt", header=F)
laci <- read.delim("lacIAUC.txt", header=F)
hiv <- read.delim("HIVAUC.txt", header=F)
human <- read.delim("humanAUC.txt", header=T)

bestdist <- read.table("custom_lacI.scan", comment.char="#", header=F)
bestdist$dbs <- (bestdist[bestdist$V4=="wild_type",8] - bestdist[,8])
bestdistpositive <- bestdist[grep("-", bestdist$V4),24]
bestdistpositive <- cbind(bestdistpositive, rep(1, length(bestdistpositive)))
bestdistnegative <- bestdist[-(grep("-", bestdist$V4)),24]
bestdistnegative <- cbind(bestdistnegative, rep(0, length(bestdistnegative)))
bestdistROC <- rbind(bestdistpositive, bestdistnegative)
bestpred <- prediction(bestdistROC[,1], labels=bestdistROC[,2])
bestperf <- performance(bestpred, "tpr", "fpr")

lacpfam <- read.delim("Pfam_lacI.txt", header=F, stringsAsFactors=F)
lacpfam <- lacpfam[order(lacpfam[,1]),]
lacpfam <- lacpfam[match(bestdist[,4], lacpfam[,1]),]
lacpfamlabels <- rep(0, nrow(lacpfam))
lacpfamlabels[(grep("-", lacpfam[,1]))] = 1
lacpfampred <- prediction(lacpfam[,2], labels=lacpfamlabels)
lacpfamperf <- performance(lacpfampred, "tpr", "fpr")
lacpfam$name <- lapply(lacpfam[,1], function(x) strsplit(x, "\\W")[[1]][1])

lacprovean <- read.delim("provean_lacI.tsv", header=F, stringsAsFactors=F)
lacprovean <- lacprovean[match(lacpfam$name, lacprovean[,1]),]
lacproveanpred <- prediction((lacprovean[,2]*(-1)), labels=lacpfamlabels)
lacproveanperf <- performance(lacproveanpred, "tpr", "fpr")

lacpolyphen <- read.delim("polyphen_lacI.tsv", header=F, stringsAsFactors=F)
lacpolyphen <- lacpolyphen[match(lacpfam$name, lacpolyphen[,1]),]
lacpolyphenpred <- prediction(lacpolyphen[,2], labels=lacpfamlabels)
lacpolyphenperf <- performance(lacpolyphenpred, "tpr", "fpr")

humanneg <- read.table("treefam_polymorphism_gene_information.txt", header=T)
humanpos <- read.table("treefam_disease_gene_information.txt", header=T)
hneg <- hist(humanneg[,3], breaks=50000)
hpos <- hist(humanpos[,3], breaks=50000)

png("fig1.png", width=900, height=700)

par(fig=c(0,1,0.5,1), mar=c(7,5,3,12), cex.main=1.4, cex.axis=1.2)

plot(c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), lysozyme[1:11], lwd=2.5, type="l", main="AUC values for different methods across four benchmarking data sets", col="cadetblue", ylim=c(0.52,0.92), xlim=c(0, 170), xaxt="n", xlab="", ylab="AUC", cex.lab=1.2)
lines(c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), laci[1:11], col="coral4", lwd=2.5)
lines(c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), hiv[1:11], col="antiquewhite4", lwd=2.5)
points(c(115, 145, 155), lysozyme[12:14], col="cadetblue", pch=16, cex=1.3)
points(c(115, 145, 155), laci[12:14], col="coral4", pch=16, cex=1.3)
points(c(115, 145, 155), hiv[12:14], col="antiquewhite4", pch=16, cex=1.3)
points(c(115, 125, 145, 155, 165), human[c(12,13,15,16,17)], col="navy", pch=16, cex=1.3)
points(c(135,134,135,136), y=c(lysozyme[5], laci[5], hiv[5], human[14]), col=c("cadetblue", "coral4", "antiquewhite4", "navy"), pch=16, cex=1.3)
points(165, human[18], col="navy", pch=21, cex=1.3, lwd=3)

text(-24,0.98, labels="A", xpd=T, cex=2)

axis(1, at=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), labels=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
axis(1, at=c(115, 125, 135, 145, 155, 165), labels=c("Pfam DBS", "Treefam DBS", "Custom DBS", "PROVEAN", "PolyPhen-2", "FATHMM"), las=2)

legend(x=180, y=0.85, cex=1.3, xpd=T,fill=c("cadetblue", "coral4", "antiquewhite4", "navy"), legend=c("phage lysozyme", "LacI", "HIV protease", "Human variants"), title=expression(bold("Data set")))

text(50, 0.35, labels="Percent identity\n(Custom-built HMMs)", xpd=T, cex=1.3)

par(fig=c(0,0.5, 0,0.5), new=T, mar=c(5,5,5,5))

plot(bestperf, lwd=2.5, main="ROC curve for LacI\n40%identity model", cex.lab=1.4)
plot(lacpfamperf, add=T, lwd=2.5, lty=4)
plot(lacproveanperf, add=T, lwd=2.5, lty=2)
plot(lacpolyphenperf, add=T, lwd=2.5, lty=3)
legend(x=0.4, y=0.5, legend=c("Custom", "PROVEAN", "PolyPhen-2", "Pfam"), lty=c(1,2,3,4), xpd=T, lwd=2.5, cex=1.3)

text(-0.26,1.25, labels="B", xpd=T, cex=2)
text(1.25,1.25, labels="C", xpd=T, cex=2)

par(fig=c(0.5,1,0,0.5), new=T, mar=c(5,5,5,5))

plot(density(humanneg[abs(humanneg[,3])<10,3], bw=0.2), xlim=c(-10,10), col="navy", lwd=2.5, main="DBS for human polymorphisms\nTreefam models", cex.lab=1.4, xlab="Delta-bitscore", ylab="Density")
lines(density(humanpos[abs(humanpos[,3])<10,3], bw=0.2), col="coral4", lwd=2.5)

legend(x=1, y=0.38, legend=c("Polymorphic", "Disease"), fill=c("navy", "coral4"), xpd=T, cex=1.3)

dist <- rnorm(n=nrow(humanneg), m=mean(humanneg[abs(humanneg[,3])<5,3]), sd=sd(humanneg[abs(humanneg$DBS)<5,3]))

dev.off()
