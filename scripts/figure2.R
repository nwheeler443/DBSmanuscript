dbs <- read.delim("ent-gal.dbs", comment.char="#", header=F, stringsAsFactors=F)

dbsgenes <- unique(dbs[,1])
scores = data.frame()
for (i in 1:length(dbsgenes)) {
    scores[i,1] = sum(dbs[!is.na(match(dbs[,1], dbsgenes[i])),10])
}

dnds <- read.table("dndsvalues.tsv", header=TRUE)

dist <- rnorm(n=50000, mean=0, sd=2.13290575206885)
distdens <- density(dist, bw=0.3)

combined <- cbind(dnds, as.numeric(scores[match(dnds[,1], dbsgenes),1]))
combined <- combined[!is.na(combined[,6]),]
combined <- combined[combined[,6]>0 & combined$Ka>0 & combined$Ks>0,]

redcolour <- rgb(170,0,0,70, max=255)
bluecolour <- rgb(70, 120, 190,70, max=255)

pathdata <- read.table("pathloss.genes.unique.txt", header = T, sep = "\t")
pathdata <- pathdata[c(1,2,3,4,5,6,8,11,16,17,19,21),]
pathdata <- pathdata[nrow(pathdata):1,]

pathdata$totalgenes <- pathdata[,1]*2+(pathdata[,4]+pathdata[,5])*2
pathdata[,7:10] <- pathdata[,2:5]/pathdata$totalgenes						# losses and uniques as fractions of total
pathdata[,11:12] <- (pathdata[,1]-pathdata[,2:3]+pathdata[,4:5])/pathdata[,6]			# functional proportion of genome

pathdata <- pathdata[pathdata[,12]!="Inf",]
pathdataT<-matrix(c(pathdata[,10], pathdata[,7], pathdata[,11], pathdata[,12], pathdata[,8], pathdata[,9]),nrow=nrow(pathdata),ncol=6,dimnames=list(c(as.character(rownames(pathdata))),c("unique to pathogen", "loss in commensal", "commensal", "pathogen", "loss in pathogen", "unique to commensal")))

simpleCap <- function(x) {
	s <- strsplit(x, " ")[[1]]
	paste(toupper(substring(s, 1,1)), substring(s, 2),
	sep="", collapse=" ")
}

dbshist <- hist(log(combined[,6]+0.01), breaks=15)
dndshist <- hist(log10(combined[,5]), breaks=50)

################################################################################

plotfig <- function() {
	
	par(fig=c(0,0.78,0.7,1), mar=c(3.5,9,2,1), mgp=c(2,0.5,0))
    hist(scores[is.na(match(dbsgenes, exclude[,1])),1], xlab="Delta-bitscore", main="",xlim=c(-8,8), breaks=300, freq=FALSE, ylim=c(0,0.5))
	lines(distdens$x, distdens$y, col="coral4", lwd=2)
	text(-13.2, 0.52, labels="A", font=2, xpd=T, cex=1.4)
	
	###################
	
	op<-par(fig=c(0,1,0.4,0.7), las=2, mar=c(0,14,1,16.5), xpd=T, new=TRUE)			# mfrow for mulit-panelled plot
	bp <- barplot(t(pathdataT) , col=c(rgb(30, 30, 100, max=255), rgb(70, 95, 190, max=255), rgb(130, 160, 255, 200,max=255), rgb(205, 82, 82,200, max=255), rgb(145, 42, 42,240, max=255), rgb(100, 0, 0, max=255)),horiz=TRUE, xaxt="n")
	
	text(x=c(rep(.40, nrow(pathdataT)), rep(.60, nrow(pathdataT))), y=bp, labels=c(pathdata[,6]/2-pathdata[,5]-pathdata[,2],pathdata[,6]/2-pathdata[,4]-pathdata[,3]))
	
	legend(title="Classification of genes", x=1.1, y=nrow(pathdata), c("absent in Enteritidis", "loss of function in Enteritidis", "functional in Enteritidis", "functional in Gallinarum", "loss of function in Gallinarum", "absent in Gallinarum"), fill=c(rgb(30, 30, 100, max=255), rgb(70, 95, 190, max=255), rgb(130, 160, 255, max=255), rgb(205, 82, 82, max=255), rgb(145, 42, 42, max=255), rgb(100, 0, 0, max=255)))
	text(x=c(.25,.75), y=nrow(pathdataT)+3, labels=c("Enteritidis", "Gallinarum"), font=2)
	text(y=nrow(pathdataT)+3, x=-0.7, labels="B", font=2, cex=1.4)
	
	####################
	
	par(fig=c(0,0.93,0,0.36), mar=c(4,5,2,3), new=T, mgp=c(2.5,0.7,0))
	plot(combined[,5], abs(combined[,6])+0.01, xlab=expression(bold(paste("d"[N],"/d"[S]))), ylab=expression(bold(paste("|DBS|"))), pch=20, log="xy", xaxt="n", yaxt="n", cex.lab=1.2)
	axis(side=1, at=c(0.01,0.1,1,10,100), labels=c(0,0.1,1,10,100), cex.axis=1, las=0)
	axis(side=2, at=c(0.01,0.1,1,10,100), labels=c(0,0.1,1,10,100), cex.axis=1)
	
	par(xpd=F)
	
	polygon(x=c(0.0001,0.0001,0.5,0.5), y=c(6.6,0.0001,0.0001,6.6), col=rgb(70, 120, 190,40, max=255), lty=2)
	polygon(x=c(0.5,0.5,1.5,1.5), y=c(6.6,0.0001,0.0001,6.6), col=rgb(70, 120, 190,60, max=255), lty=2)
	polygon(x=c(1.5,1.5,1000,1000), y=c(6.6,0.0001,0.0001,6.6), col=rgb(70, 120, 190,80, max=255), lty=2)
	polygon(x=c(0.0001,0.0001,0.5,0.5), y=c(200,6.6,6.6,200), col=rgb(170,30,0,40, max=255), lty=2)
	polygon(x=c(0.5,0.5,1.5,1.5), y=c(6.6,200,200,6.6), col=rgb(170,30,0,60, max=255), lty=2)
	polygon(x=c(1.5,1.5,1000,1000), y=c(6.6,200,200,6.6), col=rgb(170,30,0,80, max=255), lty=2)
	
	text(200, 0.1, labels="functionally\ninsignificant", xpd=T, srt=270, font=2)
	text(200, 25, labels="functionally\nsignificant", xpd=T, srt=270,font=2)
	text(0.04, 185, labels="negative", xpd=T, font=2)
	text(0.85, 185, labels="neutral", xpd=T, font=2)
	text(20, 185, labels="positive", xpd=T, font=2)
	
	text(20, 0.04, labels=paste("n=",nrow(combined[combined[,5]>1.5&abs(combined[,6])<6.6,]), sep=""), font=2)
	text(0.9, 0.04, labels=paste("n=",nrow(combined[combined[,5]<1.5&combined[,5]>0.5&abs(combined[,6])<6.6,]), sep=""), font=2)
	text(0.04, 0.04, labels=paste("n=",nrow(combined[combined[,5]<0.5&abs(combined[,6])<6.6,]), sep=""), font=2)
	text(20, 40, labels=paste("n=",nrow(combined[combined[,5]>1.5&abs(combined[,6])>6.6,]), sep=""), font=2)
	text(0.9, 40, labels=paste("n=",nrow(combined[combined[,5]<1.5&combined[,5]>0.5&abs(combined[,6])>6.6,]), sep=""), font=2)
	text(0.04, 40, labels=paste("n=",nrow(combined[combined[,5]<0.5&abs(combined[,6])>6.6,]), sep=""), font=2)
	
	par(fig=c(0,0.88,0.35,0.39), new=TRUE, mar=c(0,4,0,0))
	barplot(dndshist$counts, axes=FALSE, space=0)
	text(y=30, x=-4.3, labels="C", font=2, xpd=T, cex=1.4)
	
	par(fig=c(0.9,0.97,0.016,0.35),new=TRUE, mar=c(2.8,1,1,0))
	barplot(dbshist$counts, horiz=TRUE, axes=FALSE, space=0)

}

png("figure2.png", width=700, height=800)
plotfig()
dev.off()

pdf("figure2.pdf", width=10, height=12)
plotfig()
dev.off()

cor(abs(combined[,6]), combined$KaKs, method="spearman")
cor(abs(combined[,6]), combined$Ka, method="spearman")

