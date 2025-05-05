

#### Counts and expression matrices

library(dplyr)
library(tibble)

setwd("Desktop/Neurogenomics/Projectes/CortexMalformation/DATA/BulkLines/")

samples <- sapply(strsplit(as.character(dir("counts")), split="\\."), "[", 1) %>% unique
  
countdata <- c()
for (s in samples){
  c <- read.table(paste0("counts/",s,".counts.txt"), header=TRUE, row.names=1)
  c1 <- c[ ,6:ncol(c)]
  countdata <- cbind(countdata, c1)
}  
colnames(countdata) <- sapply(strsplit(as.character(samples), split="_"), "[", 1)
rownames(countdata) <- rownames(c)

dup1Alpha <- which(colnames(countdata)=="1alpha")
dup6 <- which(colnames(countdata)=="6")

sum1Alpha <- rowSums(countdata[,dup1Alpha])
sum6 <- rowSums(countdata[,dup6])

countdata_filt <- countdata[,-c(dup1Alpha, dup6)]
countdata_filt <- as.data.frame(countdata_filt)
colnames(countdata_filt) <- paste0("S_",colnames(countdata_filt))
countdata_filt$S_6 <- sum6
countdata_filt$S_1alpha <- sum1Alpha

geneinfo <- read.table("geneInfo.tab", sep="\t", header=FALSE)

#### Load metadata

meta <- read.table("Metadata.txt", header=T, sep="\t")
meta$sample.ID <-  paste0("S_",meta$sample.ID)

countdata_filt <- countdata_filt[,match( meta$sample.ID, colnames(countdata_filt) )]
countdata_filt <- add_column(countdata_filt, GENE=rownames(countdata_filt), .before=1)

saveRDS(meta, file="Metadata.rds")
saveRDS(countdata_filt, file="BulkLines_RAWcounts.rds")
saveRDS(geneinfo, file = "GeneInfo.rds")

library(scater)

length <- c$Length
tpm <- calculateTPM(x = countdata_filt[,-1], lengths = length)
tpm <- add_column(as.data.frame(tpm), GENE=rownames(tpm), .before=1)

saveRDS(tpm, file="BulkLines_TPM.rds")

library(EDASeq)
plotRLE(as.matrix(log2(tpm[,-1])))

dv <- grep("DIV", meta$Passage2)


#### RUN LDA

require(DESeq)
counts <- countdata_filt[,-1]
counts <- counts[,dv]
lib.size <- estimateSizeFactorsForMatrix(counts)
ed <- t(t(counts)/lib.size)
means <- rowMeans(ed)
vars <- apply(ed,1,var)
cv2 <- vars/means^2
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9)
smoothScatter(log(means),log(cv2))
require(statmod)
minMeanForFit <- unname( quantile( means[ which( cv2 > .3 ) ], .95 ) )
useForFit <- means >= minMeanForFit # & spikeins
fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
fit$coefficients
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2));
xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
vfit <- a1/xg + a0
# add fit line
lines( log(xg), log(vfit), col="black", lwd=3 )
df <- ncol(ed) - 1
# add confidence interval
lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")
afit <- a1/means+a0
varFitRatio <- vars/(afit*means^2)
varorder <- order(varFitRatio,decreasing=T)
oed_dv <- ed[varorder,]
# save for the next exercise
save(oed_dv,file="oed_dv.RData")

# repeat previous plot
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2)); lines( log(xg), log(vfit), col="black", lwd=3 ); lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black"); lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black");
# add top 100 genes
points(log(means[varorder[1:100]]),log(cv2[varorder[1:100]]),col=2)


topvariable_p <- rownames(oed_p[1:100,])
topvariable_dv <- rownames(oed_dv[1:100,])

library(MASS)
library(reshape2)
data <- t(log2(tpm[,-1] + 1))


#### Using Passage Disease 

data_p  <- data[-dv,] 
data_p <- data_p[,colnames(data_p) %in% topvariable_p]
data_p <- add_column(as.data.frame(data_p), Variable=paste0(meta[-dv,"Passage2"], "_", meta[-dv,"Disease"]), .before=1)
model_p <- lda(Variable~., data = data_p)
plda <- predict(object = model_p,
                newdata = data_p)
dataset = cbind(meta[-dv,], plda$x)
p1 <- ggplot(dataset) + geom_point(aes(LD1, LD2, colour = Passage2, shape = Disease), size = 2.5) + 
  theme_linedraw()

p2 <- ggplot(dataset) + geom_point(aes(LD1, LD3, colour = Passage2, shape = Disease), size = 2.5) + 
  theme_linedraw()

p3 <- ggplot(dataset) + geom_point(aes(LD2, LD3, colour = Passage2, shape = Disease), size = 2.5) + 
  theme_linedraw()

multiplot(p1, p2, p3, cols=3)

dataset2 <- melt(dataset, id.vars = colnames(dataset)[1:5])
p2 <- ggplot(dataset2, aes(x=Disease, y=value, fill=Disease)) + geom_violin() + 
  facet_grid(variable~.) + scale_fill_manual(values = c("tomato", "cadetblue")) + theme_linedraw()

p3 <- ggplot(dataset2, aes(x=Passage2, y=value, fill=Passage2)) + geom_violin() + 
  facet_grid(variable~.) + scale_fill_manual(values = c("tomato", "cadetblue", "tan")) + theme_linedraw()

p4 <- ggplot(dataset2, aes(x=X20ng.ml.FGF2, y=value, fill=X20ng.ml.FGF2)) + geom_violin() + 
  facet_grid(variable~.) + scale_fill_manual(values = c("tomato", "cadetblue")) + theme_linedraw()

multiplot(p2, p3, p4, cols=3)

loadings <- as.data.frame(model_p$scaling)
loadings <- merge(loadings, geneinfo, by.x=0, by.y="V1")
loadings <- loadings[order(loadings$LD3,decreasing=T), ]

gene2plot <- cbind(meta[-dv,], data_p[,"ENSG00000143858"])
colnames(gene2plot)[6] <- "gene"
p1 <- ggplot(gene2plot, aes(x=X20ng.ml.FGF2, y=gene, fill=X20ng.ml.FGF2)) + geom_boxplot() + 
   scale_fill_manual(values = c("tomato", "cadetblue")) + theme_linedraw()
p2 <- ggplot(gene2plot, aes(x=Disease, y=gene, fill=Disease)) + geom_boxplot() + 
  scale_fill_manual(values = c("tomato", "cadetblue")) + theme_linedraw()
p3 <- ggplot(gene2plot, aes(x=Passage2, y=gene, fill=Passage2)) + geom_boxplot() + 
 scale_fill_manual(values = c("tomato", "cadetblue", "tan")) + theme_linedraw()
multiplot(p2, p3, p1, cols=3)
#### Using Passage Disease and FGF

data_p  <- data[-dv,] 
data_p <- data_p[,colnames(data_p) %in% topvariable_p]
data_p <- add_column(as.data.frame(data_p), Variable=paste0(meta[-dv,"Passage2"], "_", meta[-dv,"Disease"], "_", meta[-dv,"X20ng.ml.FGF2"]), .before=1)
model_p <- lda(Variable~., data = data_p)
plda <- predict(object = model_p,
                newdata = data_p)
dataset = cbind(meta[-dv,], plda$x)

p1 <- ggplot(dataset) + geom_point(aes(LD1, LD2, colour = Passage2, shape = Disease, alpha = as.factor(X20ng.ml.FGF2)), size = 2.5) + 
  theme_linedraw() + scale_alpha_manual(values = c(0.5, 1))

p2 <- ggplot(dataset) + geom_point(aes(LD1, LD3, colour = Passage2, shape = Disease, alpha = as.factor(X20ng.ml.FGF2)), size = 2.5) + 
  theme_linedraw() + scale_alpha_manual(values = c(0.5, 1))

p3 <- ggplot(dataset) + geom_point(aes(LD2, LD3, colour = Passage2, shape = Disease, alpha = as.factor(X20ng.ml.FGF2)), size = 2.5) + 
  theme_linedraw() + scale_alpha_manual(values = c(0.5, 1))

p4 <- ggplot(dataset) + geom_point(aes(LD1, LD4, colour = Passage2, shape = Disease, alpha = as.factor(X20ng.ml.FGF2)), size = 2.5) + 
  theme_linedraw() + scale_alpha_manual(values = c(0.5, 1))

p5 <- ggplot(dataset) + geom_point(aes(LD3, LD4, colour = Passage2, shape = Disease, alpha = as.factor(X20ng.ml.FGF2)), size = 2.5) + 
  theme_linedraw() + scale_alpha_manual(values = c(0.5, 1))

p6 <- ggplot(dataset) + geom_point(aes(LD7, LD4, colour = Passage2, shape = Disease, alpha = as.factor(X20ng.ml.FGF2)), size = 2.5) + 
  theme_linedraw() + scale_alpha_manual(values = c(0.5, 1))

multiplot(p1, p2, p3, p4, p5, p6, cols=3)

dataset2 <- melt(dataset, id.vars = colnames(dataset)[1:5])
p2 <- ggplot(dataset2, aes(x=Disease, y=value, fill=Disease)) + geom_violin() + 
  facet_grid(variable~.) + scale_fill_manual(values = c("tomato", "cadetblue")) + theme_linedraw()

p3 <- ggplot(dataset2, aes(x=Passage2, y=value, fill=Passage2)) + geom_violin() + 
  facet_grid(variable~.) + scale_fill_manual(values = c("tomato", "cadetblue", "tan")) + theme_linedraw()

p4 <- ggplot(dataset2, aes(x=X20ng.ml.FGF2, y=value, fill=X20ng.ml.FGF2)) + geom_violin() + 
  facet_grid(variable~.) + scale_fill_manual(values = c("tomato", "cadetblue")) + theme_linedraw()

multiplot(p2, p3, p4, cols=3)

loadings <- as.data.frame(model_p$scaling)
loadings <- merge(loadings, geneinfo, by.x=0, by.y="V1")
loadings <- loadings[order(loadings$LD4,decreasing=T), ]
write.table(loadings, file = "Loadings_LD_withFGF.txt", quote=F, row.names = FALSE, sep="\t")

gene2plot <- cbind(meta[-dv,], data_p[,"ENSG00000162692"])
gene2plot$Var <- paste0(gene2plot$Passage2, "_",gene2plot$X20ng.ml.FGF2, "_",gene2plot$Disease)
colnames(gene2plot)[6] <- "gene"
p1 <- ggplot(gene2plot, aes(x=X20ng.ml.FGF2, y=gene, fill=X20ng.ml.FGF2)) + geom_boxplot() + 
  scale_fill_manual(values = c("tomato", "cadetblue")) + theme_linedraw()
p2 <- ggplot(gene2plot, aes(x=Disease, y=gene, fill=Disease)) + geom_boxplot() + 
  scale_fill_manual(values = c("tomato", "cadetblue")) + theme_linedraw()
p3 <- ggplot(gene2plot, aes(x=Passage2, y=gene, fill=Passage2)) + geom_boxplot() + 
  scale_fill_manual(values = c("tomato", "cadetblue", "tan")) + theme_linedraw()
p4 <- ggplot(gene2plot, aes(x=Disease, y=gene, fill=Disease)) + geom_boxplot() + 
  scale_fill_manual(values = c("tomato", "cadetblue")) + theme_linedraw() + facet_grid(X20ng.ml.FGF2~Passage2)
multiplot(p2, p3, p1, cols=3)

####### Dorsoventral

data_dv <- data[dv,]
data_dv <- data_dv[,colnames(data_dv) %in% topvariable_dv]
data_dv <- add_column(as.data.frame(data_dv), Variable=paste0(meta[dv,"Passage2"], "_", meta[dv,"Disease"]), .before=1)

model_p <- lda(Variable~., data = data_dv)
plda <- predict(object = model_p,
                newdata = data_dv)
dataset = cbind(meta[dv,], plda$x)
dataset$Passage2 <- factor(dataset$Passage2, levels = paste0("DIV ", c(8,17,30,38)))
p1 <- ggplot(dataset) + geom_point(aes(LD1, LD2, colour = Passage2, shape = Disease), size = 2.5) + 
  theme_linedraw()

p2 <- ggplot(dataset) + geom_point(aes(LD1, LD3, colour = Passage2, shape = Disease), size = 2.5) + 
  theme_linedraw()

p3 <- ggplot(dataset) + geom_point(aes(LD2, LD3, colour = Passage2, shape = Disease), size = 2.5) + 
  theme_linedraw()

multiplot(p1, p2, p3, cols=3)

dataset2 <- melt(dataset, id.vars = colnames(dataset)[1:5])
p2 <- ggplot(dataset2, aes(x=Disease, y=value, fill=Disease)) + geom_violin() + 
  facet_grid(variable~.) + scale_fill_manual(values = c("tomato", "cadetblue")) + theme_linedraw()

p3 <- ggplot(dataset2, aes(x=Passage2, y=value, fill=Passage2)) + geom_violin() + 
  facet_grid(variable~.) + scale_fill_manual(values = c("chocolate4", "chocolate3", "chocolate2", "chocolate1")) + theme_linedraw()

multiplot(p2, p3,  cols=2)

loadings <- as.data.frame(model_p$scaling)
loadings <- merge(loadings, geneinfo, by.x=0, by.y="V1")
loadings <- loadings[order(loadings$LD3,decreasing=T), ]
write.table(loadings, file = "Loadings_LD_DV.txt", quote=F, row.names = FALSE, sep="\t")

gene2plot <- cbind(meta[dv,], data_dv[,"ENSG00000134438"])
colnames(gene2plot)[6] <- "gene"
gene2plot$Passage2 <- factor(gene2plot$Passage2, levels = paste0("DIV ", c(8,17,30,38)))
p2 <- ggplot(gene2plot, aes(x=Disease, y=gene, fill=Disease)) + geom_boxplot() + 
  scale_fill_manual(values = c("tomato", "cadetblue")) + theme_linedraw()
p3 <- ggplot(gene2plot, aes(x=Passage2, y=gene, fill=Passage2)) + geom_boxplot() + 
  scale_fill_manual(values = c("chocolate4", "chocolate3", "chocolate2", "chocolate1")) + theme_linedraw()
multiplot(p2, p3, cols=2)

ggplot(gene2plot, aes(x=Disease, y=gene, fill=Disease)) + geom_boxplot() + 
  scale_fill_manual(values = c("tomato", "cadetblue")) + theme_linedraw() + facet_grid(.~Passage2)


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
