setwd("C:/Users/Jack/Desktop/R_Lab/final_project")

ndat <- read.table("GDS5826.dat",header=T,row.names=1,sep="\t")
dat <- ndat[,2:13]
dat.order <-c("GSM1692587","GSM1692588","GSM1692589","GSM1692590","GSM1692591","GSM1692592","GSM1692593","GSM1692594","GSM1692595","GSM1692596","GSM1692597","GSM1692598")
print(dat.order)
dat <- dat[dat.order]

#Testing for outlier samples 
# Correlation plot
library(gplots)

dat.cor <- cor(dat)
reorder_dat <- dat[dat.order]

layout(matrix(c(1,1,1,1,1,1,1,1,2,2), 5, 2, byrow = TRUE))
par(oma=c(5,7,1,1))

cx <- rev(colorpanel(25,"blue","white","red"))
leg <- seq(min(dat.cor,na.rm=T),max(dat.cor,na.rm=T),length=10)

image(dat.cor,main="Correlation plot of carfilzomib-resistant MM cell lines",axes=F,col=cx)
axis(1,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis=0.9,las=2)
axis(2,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis=0.9,las=2)
#Remove whitespace abouve second plot
par(mar=c(5,4,2,2)+0.1)
#plot gradient
image(as.matrix(leg),col=cx,axes=FALSE)
tmp <- round(leg,2)
axis(1,at=seq(0,1,length=length(leg)),labels=tmp,cex.axis=1)


dev.copy(png,'correlation_plot.png')
dev.off()

# hierarchical clustering plot
# transpose dat
dat.t <- t(dat)
# get pairwise distances
dat.dist <- dist(dat.t)
# calculate and plot hierarchical tree
plot(hclust(dat.dist), labels=dimnames(dat)[[2]],main="Hierarchical Clustering Dendrogram of carfilzomib-resistant MM cell lines",cex.main=0.75)
dev.copy(png,'hierarchical_clustering_plot.png')
dev.off()


dat.mean <- apply(dat,2,mean)
# calculate samples standard deviations
dat.sd <- apply(dat,2,sd)

# calculate cv of samples
dat.cv <- dat.sd/dat.mean

# create CV plot
plot(dat.mean,dat.cv,main="carfilzomib-resistant MM cell lines CV vs. Mean",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n")
text(dat.mean,dat.cv,label=dimnames(dat)[[2]])
dev.copy(png,'CV_plot.png')
dev.off()

# average correlation plot
dat.avg <- apply(dat.cor,1,mean)
par(oma=c(6,0.1,0.1,0.1))
plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avg r",main="Avg correlation of carfilzomib-resistant MM cell lines",axes=F)
points(dat.avg,bg="red",col=1,pch=21,cex=1.5)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat)[[2]],las=2,cex=0.5)
axis(2)
abline(v=seq(0.5,58.5,1),col="grey")
dev.copy(png,'avg_CV_plot.png')
dev.off()

#Based on the suit of plots
#These are the outilers identified

#Cv/Mean plot
#GSM1692593


library(impute)
# remove outlier samples
length(dat)
daty <- dat[,-which(names(dat) %in% c("GSM1692593"))]

row_mean <- rowMeans(daty)
print(mean(row_mean))
#average row mean of all the columns is 0.006494236
#Filter out average gene expression < 50
datz <- daty[rowMeans(daty)>0.006494236,]
length(datz)

print(colnames(datz))

t.test.all.genes <- function(x,s1,s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}

t.test.run <- apply(datz,1,t.test.all.genes,s1=c(1:6),s2=c(7:11))
hist(t.test.run,main="Histogram of p-values using a Student's t-test\nfrom MM Study",xlab="p-values",col="lightblue")
dev.copy(png,'hist_plot.png')
dev.off()

sum(t.test.run<.001)
#61 genes with p<0.001
#filter out genes with p>0.001
gene_list <- c()
for (i in 1:length(t.test.run)){
  if (t.test.run[i]<.001){
    gene_list = append(gene_list,names(t.test.run[i]))
  }
}

dat.filter <- datz[gene_list,]
#Number of genese retained:61

t.test.run <- apply(dat.filter,1,t.test.all.genes,s1=c(1:6),s2=c(7:11))
hist(t.test.run,main="Histogram of p-values using a Student's t-test\nfrom MM Study filtered data",xlab="p-values",col="lightblue")
dev.copy(png,'hist_plot2.png')
dev.off()

dat.pca <- prcomp(t(dat.filter),cor=F)
print(dat.pca)
dat.loadings <- dat.pca$x[,1:3]
plot(range(dat.loadings[,1]),range(dat.loadings[,2]),xlab='p1',ylab='p2',main='PCA plot of MM Study')
points(dat.loadings[,1][1:6], dat.loadings[,2][1:6],col='red',pch=16,cex=1.5)
points(dat.loadings[,1][7:11], dat.loadings[,2][7:11],col='blue',pch=16,cex=1.5)
leg.names<-c('resistant to 12 nM carfilzomib','Normal')
leg.col = c("red","blue")
legend(-1,.6,leg.names,leg.col,pch=15,cex=.7,horiz=F)
dev.copy(png,'pca_plot.png')
dev.off()


dat.pca.var <- round(dat.pca$sdev^2 / sum(dat.pca$sdev^2)*100,2)
plot(c(1:length(dat.pca.var)),dat.pca.var,type='b',xlab='# components',ylab='% variance',main='Scree plot')
dev.copy(png,'scree_plot.png')
dev.off()

#about 90% variance

# The weighted graph Laplacian
k.speClust2 <- function (X, qnt=NULL) {
  dist2full <- function(dis) {
    	n <- attr(dis, "Size")
            full <- matrix(0, n, n)
            full[lower.tri(full)] <- dis
            full + t(full)
    }
  dat.dis <- dist(t(X),"euc")^2
  if(!is.null(qnt)) {eps <- as.numeric(quantile(dat.dis,qnt))}
  if(is.null(qnt)) {eps <- min(dat.dis[dat.dis!=0])}
  kernal <- exp(-1 * dat.dis/(eps))
  K1 <- dist2full(kernal)
  diag(K1) <- 0
  D = matrix(0,ncol=ncol(K1),nrow=ncol(K1))
  tmpe <- apply(K1,1,sum)
  tmpe[tmpe>0] <- 1/sqrt(tmpe[tmpe>0])
  tmpe[tmpe<0] <- 0
  diag(D) <- tmpe
  L <- D%*% K1 %*% D
  X <- svd(L)$u
  Y <- X / sqrt(apply(X^2,1,sum))
}

temp <- t(dat.filter)
temp <- scale(temp,center=T,scale=T) 
phi <- k.speClust2(t(temp),qnt=NULL)
plot(range(phi[,1]),range(phi[,2]),xlab="phi1",ylab="phi2",main="Weighted Graph Laplacian plot")
points(phi[,1][1:6],phi[,2][1:6],col="red",pch=16,cex=1.5)
points(phi[,1][7:11],phi[,2][7:11],col="blue",pch=16,cex=1.5)
legend(-.3,.2,c('resistant to 12 nM carfilzomib','Normal'),col=c("red", "blue"),pch=15,cex=.7,horiz=F)
dev.copy(png,'laplacian_plot.png')
dev.off()

dat2 <- dat.filter
colnames(dat2)[1:6] <- "C"
colnames(dat2)[7:11] <- "N"
dat2 <- t(dat2)
dat.hca <- hclust(dist(dat2,"man"),method="median")
plot(dat.hca,main='HCA of MM Study')
dev.copy(png,'hca_plot.png')
dev.off()

heatmap(as.matrix(t(dat2)),labCol=ann,main='2D HCA MM data; 61 genes',cex.main=0.7)
dev.copy(png,'heatmat_plot.png')
dev.off()


dat.pca <- prcomp(t(dat2))
dat.loads <- dat.pca$x[,1:2]

# k-means clustering of samples
dat.k <- kmeans(dat.loads,2,61)$cluster

print(dat.k)

# plot samples in 2-dimensions with labels, colored from kmeans cluster membership
plot(range(dat.loads[,1]),range(dat.loads[,2]),type='n',xlab='p1',ylab='p2',main='PCA with K-means classification')

k.class <- c(1:2)
k.col <- c("red","blue")
species <- c("b","g")
ann <- c("b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","b","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g","g")
for(i in 1:length(species)) {
  text(dat.loads[ann==species[i],1],dat.loads[ann==species[i],2],label=ann[ann==species[i]])
  for(j in 1:length(k.class)) {
    if(sum(ann==species[i] & dat.k==k.class[j])>0) {
      text(dat.loads[(ann==species[i] & dat.k==k.class[j]),1],dat.loads[(ann==species[i] & dat.k==k.class[j]),2],
           label=ann[(ann==species[i] & dat.k==k.class[j])],col=k.col[j])
    }
  }
}

print(dat.loads[,1])
#Top 5 discriminate genes
tail(sort(dat.loads[,1]),5)
head(sort(dat.loads[,1]),5)
print(rownames(dat.loads))



















