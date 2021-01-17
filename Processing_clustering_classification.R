# LOAD CYTOF DATA
load("PDX_characterization_raw.RData")

# NORMALIZATION
cosnorm<-function(x) {
  norm<-sqrt(sum(x^2))
  x.n<-x/norm
  return(x.n)
}
dataset.c<-dataset
exprs(dataset.c)<-t(apply(exprs(dataset),1,cosnorm))

# QC AND DOWNSAMPLING
th.events<-1600
events<-aggregate(fData(dataset.c)$SampleID,by=list(fData(dataset.c)$SampleID),FUN="length")
good.samples<-events$Group.1[events$x>th.events]
length(good.samples)
dataset<-dataset.c[fData(dataset.c)$SampleID%in%good.samples,]

samples<-as.vector(unique(fData(dataset)$SampleID))
IDs<-numeric()
set.seed(123)
for(i in 1:length(samples)) {
  x<-exprs(dataset)[fData(dataset)$SampleID==samples[i],]
  sample<-sample(nrow(x),th.events,replace=F)
  IDs<-c(IDs,rownames(x)[sample])
}
downset<-dataset[IDs,]

# CLUSTERING
m<-exprs(downset)
pheno<-Rphenograph(exprs(downset),k=250)
pheno.lab<-cbind(cell=rownames(m)[as.numeric(pheno$names)],pheno.clust=pheno$membership)
all.cells<-data.frame(cell=rownames(m))
pheno.merge<-merge(all.cells,pheno.lab,all.x=T,by.x="cell",by.y="cell",sort=F)
rownames(pheno.merge)<-pheno.merge$cell
pheno.merge<-pheno.merge[rownames(downset),]
fData(downset)<-cbind(fData(downset),pheno.clust=pheno.merge$pheno.clust)

# CLASSIFIER
x<-t(exprs(downset))
rownames(x)<-downset$protein
y<-fData(downset)$pheno.clust
mydata=list(x=x,y=y,genenames=rownames(x))
mytrain<-pamr.train(mydata)

