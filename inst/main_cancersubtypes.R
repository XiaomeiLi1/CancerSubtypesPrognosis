###load library
library(CancerSubtypesPrognosis)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggrepel)

### Predict breast cancer subtypes on five multi-omic data
dn = c("TCGA","MBR","UK","HEL","GSE19783")

methods = c("PAM50","IntClust","CC","CNMF","iCluster","SNF","SNF-CC","WSNF","CIMLR","PINS","NEMO")

for(omics in c("miRNA", "mRNA", "multiomics")) {
  CancerSubtypes(dn, omics, methods,
             fileFolder=omics, logFile = omics)
}

###1. plot p-values of cancer subtyping methods (Figure 1)
x1 = read.csv("miRNA/pvalueTable.csv",header = T,stringsAsFactors = F)
rownames(x1) = gsub("MBR", "METABRIC", x1$X) 
x1 = x1[,-1]
x1 = x1[c(5,4,3,1,2),3:11]
colnames(x1) = gsub("\\.", "-", colnames(x1))

x2 = read.csv("mRNA/pvalueTable.csv",header = T,stringsAsFactors = F)
rownames(x2) = gsub("MBR", "METABRIC", x2$X)
x2 = x2[,-1]
x2 = x2[c(5,4,3,1,2),3:11]
colnames(x2) = gsub("\\.", "-", colnames(x2))

x3 = read.csv("multiomics/pvalueTable.csv",header = T,stringsAsFactors = F)
rownames(x3) = gsub("MBR", "METABRIC", x3$X)
x3 = x3[,-1]
x3 = x3[c(5,4,3,1,2),3:11]
colnames(x3) = gsub("\\.", "-", colnames(x3))


library(ggplot2)
library(ggpubr)
library(ggrepel)

xdata = cbind(x1,x2,x3)

rn = c("GSE19783 (99)","HEL (115)","UK (207)","TCGA (753)","METABRIC (1283)")


create.performance.plot.split <- function(r,xdata){
  
  df = data.frame(row.names = seq(1:27))
  df$x = t(xdata[r,])
  df$y = colnames(xdata)
  df$y = gsub("\\.", "-", df$y)
  df$y = factor(df$y, level=rev(colnames(x1)))
  df$group = c(rep("miRNA",ncol(x1)),rep("mRNA",ncol(x2)),rep("miRNA-mRNA",ncol(x3)))
  #df$methods = colnames(xdata)
  #df$methods = gsub("\\.", "-", df$methods)
  # df$shape = c(rep(17,ncol(x1)),rep(18,ncol(x2)),rep(20,ncol(x3)))
  # df$color = c(rep("red",ncol(x1)),rep("blue",ncol(x2)),rep("black",ncol(x3)))
  
  p <- ggplot(df,aes(x=-log10(x),y=y,shape=group,color=group))+
    scale_shape_manual(values = c(17,18,20)) +
    scale_color_manual(values = c("red", "blue", "black")) +
    xlab(expression(bold('-log'[10]*'(p-value)'))) + 
    ylab("Method")+
    geom_point(size = 3) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          #panel.background = element_blank(),
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.title = element_text(face = "bold", size = 10),
          legend.position="bottom",
          legend.box = "horizontal",
          #legend.box.background = element_rect(colour = "black",linetype="solid"),
          legend.text = element_text(face = "bold", size = 10))+
    expand_limits(x = 0, y=0)+
    guides(color = guide_legend(title="Data type: ",nrow = 1),shape = guide_legend(title="Data type: ", nrow = 1)) +
    geom_vline(xintercept= -log10(0.05),color = "red",size=0.8)
  
  return(p)
}

summary.performance.plot <- function(r,xdata){
  df = data.frame(row.names = seq(1:27))
  df$x = t(xdata[r,])
  df$y = colnames(xdata)
  df$y = gsub("\\.", "-", df$y)
  df$y = factor(df$y, level=rev(colnames(x1)))
  df$group = c(rep("miRNA",ncol(x1)),rep("mRNA",ncol(x2)),rep("miRNA-mRNA",ncol(x3)))
  #df$methods = colnames(xdata)
  #df$methods = gsub("\\.", "-", df$methods)
  # df$shape = c(rep(17,ncol(x1)),rep(18,ncol(x2)),rep(20,ncol(x3)))
  # df$color = c(rep("red",ncol(x1)),rep("blue",ncol(x2)),rep("black",ncol(x3)))
  
  p <- ggplot(df,aes(x=-log10(x),y=y,shape=group,color=group))+
    scale_shape_manual(values = c(17,18,20)) +
    scale_color_manual(values = c("red", "blue", "black")) +
    xlab(expression(bold('Average of -log'[10]*'(p-value)'))) + 
    ylab("Method")+
    geom_point(size = 3) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          #panel.background = element_blank(),
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.title = element_text(face = "bold", size = 10),
          legend.position="bottom",
          legend.box = "horizontal",
          #legend.box.background = element_rect(colour = "black",linetype="solid"),
          legend.text = element_text(face = "bold", size = 10))+
    expand_limits(x = 0, y=0)+
    guides(color = guide_legend(title="Data type: ",nrow = 1),shape = guide_legend(title="Data type: ", nrow = 1)) +
    geom_vline(xintercept= -log10(0.05),color = "red",size=0.8)
  
  return(p)
}

temp = colMeans(xdata)
xdata = rbind(xdata,temp)
rownames(xdata)[6]= "OVERALL"
rn = c(rn, "OVERALL")

gs = list(length = 6)
for(r in 1:5) {
  gs[[r]] = create.performance.plot.split(r,xdata)
}
gs[[6]] = summary.performance.plot(6,xdata)
pdf("Figure1-6.pdf",width=6,height=15,onefile = FALSE)
ggarrange(gs[[5]], gs[[4]],gs[[3]],
          gs[[2]],gs[[1]],gs[[6]],
          ncol = 1, nrow = 6,
          #labels = c("A", "B","C"),
          align = "hv",
          common.legend = TRUE, legend = "bottom")
dev.off()

data = list()
for(i in 1:5){
  temp = matrix(data=0,nrow = 9, ncol = 3)
  rownames(temp) = colnames(x1)
  colnames(temp) = c("miRNA","miRNA-mRNA","mRNA")
  for(j in 1:9){
    a = c(x1[i,j],x3[i,j],x2[i,j])
    temp[j,which.min(a)]=1
  }
  data[[i]]=temp
}
overall = colMeans(xdata)
temp = matrix(data=0,nrow = 9, ncol = 3)
rownames(temp) = colnames(x1)
colnames(temp) = c("miRNA","miRNA-mRNA","mRNA")
for(j in 1:9){
  a = c(overall[j],overall[18+j],overall[9+j])
  temp[j,which.min(a)]=1
}
data[[6]]=temp

thebest = matrix(data = NA,nrow = 9,ncol = 6)
for(i in 1:6){
  thebest[,i] = colnames(data[[i]])[apply(data[[i]],1,function(x) which(x==1))]
}
rownames(thebest) = rownames(data[[1]])

create.performance.plot.thebest <- function(r,data){
  df = data.frame(row.names = rownames(data))
  df$y = factor(rownames(data), level=rev(colnames(x1)))
  df$x = rep("Data type for best performance",9)
  df$value = data[,r]
  p <- ggplot(df,aes(x=x,y=y,shape=value,color=value, label=value))+
    scale_shape_manual(values = c(17,18,20)) +
    scale_color_manual(values = c("red", "blue", "black")) +
    geom_point(size = 3) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          #panel.background = element_blank(),
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.title = element_text(face = "bold", size = 10),
          legend.position="bottom",
          legend.box = "horizontal",
          #legend.box.background = element_rect(colour = "black",linetype="solid"),
          legend.text = element_text(face = "bold", size = 10))+
    scale_x_discrete(expand = c(0, 0.05))+
    geom_text(nudge_x = 0.25, nudge_y = 0)
  return(p)
}

create.performance.plot.thebest2 <- function(r,data){
  df = data.frame(row.names = rownames(data))
  df$y = factor(rownames(data), level=rev(colnames(x1)))
  df$x = rep("Data type for best performance",9)
  df$value = data[,r]
  p <- ggplot(df,aes(x=x,y=y,shape=value,color=value, label=value))+
    scale_shape_manual(values = c(18,20)) +
    scale_color_manual(values = c("blue", "black")) +
    geom_point(size = 3) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          #panel.background = element_blank(),
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.title = element_text(face = "bold", size = 10),
          legend.position="bottom",
          legend.box = "horizontal",
          #legend.box.background = element_rect(colour = "black",linetype="solid"),
          legend.text = element_text(face = "bold", size = 10))+
    scale_x_discrete(expand = c(0, 0.05))+
    geom_text(nudge_x = 0.25, nudge_y = 0)
  return(p)
}

gsn = list(length = 6)
for(r in 1:3) {
  gsn[[r]] = create.performance.plot.thebest(r,thebest)
}
for(r in 4:5) {
  gsn[[r]] = create.performance.plot.thebest2(r,thebest)
}
pdf("Figure2-6.pdf",width=2.8,height=15,onefile = FALSE)
ggarrange(gsn[[5]], gsn[[4]],gsn[[3]],
          gsn[[2]],gsn[[1]],gsn[[6]],
          ncol = 1, nrow = 6,
          #labels = c("A", "B","C"),
          align = "hv",
          common.legend = TRUE, legend = "bottom")
dev.off()


###2. plot Silhouette scores of cancer subtyping methods (Figure 1)
x1 = read.csv("miRNA/silTable.csv",header = T,stringsAsFactors = F)
rownames(x1) = gsub("MBR", "METABRIC", x1$X) 
x1 = x1[,-1]
x1 = x1[c(5,4,3,1,2),3:11]
colnames(x1) = gsub("\\.", "-", colnames(x1))

x2 = read.csv("mRNA/silTable.csv",header = T,stringsAsFactors = F)
rownames(x2) = gsub("MBR", "METABRIC", x2$X)
x2 = x2[,-1]
x2 = x2[c(5,4,3,1,2),3:11]
colnames(x2) = gsub("\\.", "-", colnames(x2))

x3 = read.csv("multiomics/silTable.csv",header = T,stringsAsFactors = F)
rownames(x3) = gsub("MBR", "METABRIC", x3$X)
x3 = x3[,-1]
x3 = x3[c(5,4,3,1,2),3:11]
colnames(x3) = gsub("\\.", "-", colnames(x3))


library(ggplot2)
library(ggpubr)
library(ggrepel)

xdata = cbind(x1,x2,x3)
rn = c("GSE19783 (99)","HEL (115)","UK (207)","TCGA (753)","METABRIC (1283)")

create.Silhouette.plot.split <- function(r,xdata){
  
  df = data.frame(row.names = seq(1:27))
  df$x = t(xdata[r,])
  df$y = factor(colnames(xdata), level=rev(colnames(x1)))
  #df$y = gsub("\\.", "-", df$y)
  df$group = c(rep("miRNA",ncol(x1)),rep("mRNA",ncol(x2)),rep("miRNA-mRNA",ncol(x3)))
  #df$methods = colnames(xdata)
  #df$methods = gsub("\\.", "-", df$methods)
  # df$shape = c(rep(17,ncol(x1)),rep(18,ncol(x2)),rep(20,ncol(x3)))
  # df$color = c(rep("red",ncol(x1)),rep("blue",ncol(x2)),rep("black",ncol(x3)))
  
  p <- ggplot(df,aes(x=x,y=y,shape=group,color=group))+
    scale_shape_manual(values = c(17,18,20)) +
    scale_color_manual(values = c("red", "blue", "black")) +
    xlab("Silhouette score") + 
    ylab("Method")+
    geom_point(size = 3) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          #panel.background = element_blank(),
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.title = element_text(face = "bold", size = 10),
          legend.position="bottom",
          legend.box = "horizontal",
          #legend.box.background = element_rect(colour = "black",linetype="solid"),
          legend.text = element_text(face = "bold", size = 10))+
    expand_limits(x = 0, y=0)+
    guides(color = guide_legend(title="Data type: ",nrow = 1),shape = guide_legend(title="Data type: ", nrow = 1))
  
  return(p)
}

summary.Silhouette.plot <- function(r,xdata){
  
  df = data.frame(row.names = seq(1:27))
  df$x = t(xdata[r,])
  df$y = factor(colnames(xdata), level=rev(colnames(x1)))
  #df$y = gsub("\\.", "-", df$y)
  df$group = c(rep("miRNA",ncol(x1)),rep("mRNA",ncol(x2)),rep("miRNA-mRNA",ncol(x3)))
  #df$methods = colnames(xdata)
  #df$methods = gsub("\\.", "-", df$methods)
  # df$shape = c(rep(17,ncol(x1)),rep(18,ncol(x2)),rep(20,ncol(x3)))
  # df$color = c(rep("red",ncol(x1)),rep("blue",ncol(x2)),rep("black",ncol(x3)))
  
  p <- ggplot(df,aes(x=x,y=y,shape=group,color=group))+
    scale_shape_manual(values = c(17,18,20)) +
    scale_color_manual(values = c("red", "blue", "black")) +
    xlab("Silhouette score") + 
    ylab("Method")+
    geom_point(size = 3) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          #panel.background = element_blank(),
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.title = element_text(face = "bold", size = 10),
          legend.position="bottom",
          legend.box = "horizontal",
          #legend.box.background = element_rect(colour = "black",linetype="solid"),
          legend.text = element_text(face = "bold", size = 10))+
    expand_limits(x = 0, y=0)+
    guides(color = guide_legend(title="Data type: ",nrow = 1),shape = guide_legend(title="Data type: ", nrow = 1))
  
  return(p)
}

temp = colMeans(xdata)
xdata = rbind(xdata,temp)
rownames(xdata)[6]= "OVERALL"
rn = c(rn, "OVERALL")

gs = list(length = 6)
for(r in 1:5) {
  gs[[r]] = create.Silhouette.plot.split(r,xdata)
}
gs[[6]] = summary.Silhouette.plot(6,xdata)
pdf("FigureS1-1.pdf",width=6,height=15,onefile = FALSE)
ggarrange(gs[[5]], gs[[4]],gs[[3]],
          gs[[2]],gs[[1]],gs[[6]],
          ncol = 1, nrow = 6,
          #labels = c("A", "B","C"),
          align = "hv",
          common.legend = TRUE, legend = "bottom")
dev.off()

data = list()
for(i in 1:5){
  temp = matrix(data=0,nrow = 9, ncol = 3)
  rownames(temp) = colnames(x1)
  colnames(temp) = c("miRNA","miRNA-mRNA","mRNA")
  for(j in 1:9){
    a = c(x1[i,j],x3[i,j],x2[i,j])
    temp[j,which.max(a)]=1
  }
  data[[i]]=temp
}
overall = colMeans(xdata)
temp = matrix(data=0,nrow = 9, ncol = 3)
rownames(temp) = colnames(x1)
colnames(temp) = c("miRNA","miRNA-mRNA","mRNA")
for(j in 1:9){
  a = c(overall[j],overall[18+j],overall[9+j])
  temp[j,which.max(a)]=1
}
data[[6]]=temp


thebest = matrix(data = NA,nrow = 9,ncol = 6)
for(i in 1:6){
  thebest[,i] = colnames(data[[i]])[apply(data[[i]],1,function(x) which(x==1))]
}
rownames(thebest) = rownames(data[[1]])

create.performance.plot.thebest <- function(r,data){
  df = data.frame(row.names = rownames(data))
  df$y = factor(rownames(data), level=rev(colnames(x1)))
  df$x = rep("Data type for best performance",9)
  df$value = data[,r]
  p <- ggplot(df,aes(x=x,y=y,shape=value,color=value, label=value))+
    scale_shape_manual(values = c(17,18,20)) +
    scale_color_manual(values = c("red", "blue", "black")) +
    geom_point(size = 3) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          #panel.background = element_blank(),
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.title = element_text(face = "bold", size = 10),
          legend.position="bottom",
          legend.box = "horizontal",
          #legend.box.background = element_rect(colour = "black",linetype="solid"),
          legend.text = element_text(face = "bold", size = 10))+
    scale_x_discrete(expand = c(0, 0.05))+
    geom_text(nudge_x = 0.25, nudge_y = 0)
  return(p)
}

create.performance.plot.thebest2 <- function(r,data){
  df = data.frame(row.names = rownames(data))
  df$y = factor(rownames(data), level=rev(colnames(x1)))
  df$x = rep("Data type for best performance",9)
  df$value = data[,r]
  p <- ggplot(df,aes(x=x,y=y,shape=value,color=value, label=value))+
    scale_shape_manual(values = c(17,20)) +
    scale_color_manual(values = c("red", "black")) +
    geom_point(size = 3) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          #panel.background = element_blank(),
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.title = element_text(face = "bold", size = 10),
          legend.position="bottom",
          legend.box = "horizontal",
          #legend.box.background = element_rect(colour = "black",linetype="solid"),
          legend.text = element_text(face = "bold", size = 10))+
    scale_x_discrete(expand = c(0, 0.05))+
    geom_text(nudge_x = 0.25, nudge_y = 0)
  return(p)
}

gsn = list(length = 6)
for(r in 1:6) {
  gsn[[r]] = create.performance.plot.thebest(r,thebest)
}

gsn[[3]] = create.performance.plot.thebest2(3,thebest)

pdf("FigureS1-2.pdf",width=2.8,height=15,onefile = FALSE)
ggarrange(gsn[[5]], gsn[[4]],gsn[[3]],
          gsn[[2]],gsn[[1]],gsn[[6]],
          ncol = 1, nrow = 6,
          #labels = c("A", "B","C"),
          align = "hv",
          common.legend = TRUE, legend = "bottom")
dev.off()

###3. Compare multi-omics data with benchmark (Figure S2)
x1 = read.csv("miRNA/pvalueTable.csv",header = T,stringsAsFactors = F)
rownames(x1) = gsub("MBR", "METABRIC", x1$X) 
x1 = x1[,-1]
x1 = x1[c(5,4,3,1,2),]
colnames(x1) = gsub("\\.", "-", colnames(x1))

x2 = read.csv("mRNA/pvalueTable.csv",header = T,stringsAsFactors = F)
rownames(x2) = gsub("MBR", "METABRIC", x2$X)
x2 = x2[,-1]
x2 = x2[c(5,4,3,1,2),]
colnames(x2) = gsub("\\.", "-", colnames(x2))

x3 = read.csv("multiomics/pvalueTable.csv",header = T,stringsAsFactors = F)
rownames(x3) = gsub("MBR", "METABRIC", x3$X)
x3 = x3[,-1]
x3 = x3[c(5,4,3,1,2),]
colnames(x3) = gsub("\\.", "-", colnames(x3))

x1$PAM50 = x2$PAM50
x1$IntClust = x2$IntClust

x3$PAM50 = x2$PAM50
x3$IntClust = x2$IntClust

library(ggplot2)
library(ggpubr)
library(ggrepel)

xdata = cbind(x1,x2,x3)

temp = colMeans(x1)
x1 = rbind(x1,temp)
rownames(x1)[6]= "OVERALL"

temp = colMeans(x2)
x2 = rbind(x2,temp)
rownames(x2)[6]= "OVERALL"

temp = colMeans(x3)
x3 = rbind(x3,temp)
rownames(x3)[6]= "OVERALL"

rn = c("GSE19783 (99)","HEL (115)","UK (207)","TCGA (753)","METABRIC (1283)")
rn = c(rn, "OVERALL")

x1data = list()
for (r in 1:6) {
  data = matrix(data = 0, nrow = 2, ncol = 9)
  rownames(data) = colnames(x1)[1:2]
  colnames(data) = colnames(x1)[3:11]
  for (i in 1:2) {
    for (j in 1:9) {
      if (x1[r,2+j] <= x1[r,i]) data[i,j] = 1
    }
  }
  x1data[[r]] = data
}

x2data = list()
for (r in 1:6) {
  data = matrix(data = 0, nrow = 2, ncol = 9)
  rownames(data) = colnames(x2)[1:2]
  colnames(data) = colnames(x2)[3:11]
  for (i in 1:2) {
    for (j in 1:9) {
      if (x2[r,2+j] <= x2[r,i]) data[i,j] = 1
    }
  }
  x2data[[r]] = data
}

x3data = list()
for (r in 1:6) {
  data = matrix(data = 0, nrow = 2, ncol = 9)
  rownames(data) = colnames(x3)[1:2]
  colnames(data) = colnames(x3)[3:11]
  for (i in 1:2) {
    for (j in 1:9) {
      if (x3[r,2+j] <= x3[r,i]) data[i,j] = 1
    }
  }
  x3data[[r]] = data
}

xdata = list()
for (r in 1:6) {
  data = x1data[[r]]+x3data[[r]]
  xdata[[r]] = data
}

ydata = list()
for (r in 1:6) {
  data = x1data[[r]]+x2data[[r]]+x3data[[r]]
  ydata[[r]] = data
}

create.benchmark.plot.split <- function(r,xdata){
  data = xdata[[r]]
  refine = melt(as.matrix(data))
  refine$Var1 = as.factor(refine$Var1)
  refine$Var2 = as.factor(refine$Var2)
  
  
  p <- ggplot(refine, aes(Var2, Var1)) +
    geom_tile(aes(fill=value), color="white")+
    ylab(rn[r])+
    geom_text(aes(label = value), size=4) +
    theme_bw() +
    theme(legend.position = "none")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(color="Black", size=8, face="bold", hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(size = 10, angle=45, vjust=1, hjust=1))+
    scale_fill_gradient(low = "grey", high = "steelblue") +
    coord_equal()
  return (p)
}

gs = list(length = 6)
for(r in 1:6) {
  gs[[r]] = create.benchmark.plot.split(r,xdata)
}

pdf("FigureS2-2.pdf",width=5,height=8,onefile = FALSE)
ggarrange(gs[[5]], gs[[4]],gs[[3]],
          gs[[2]],gs[[1]],gs[[6]],
          ncol = 1, nrow = 6,
          #label.x = 0.5,
          #labels = c("(A)","(B)","(C)","(D)","(E)","(F)"),
          #align = "hv",
          legend = "none")
dev.off()

data = list()
for(i in 1:5){
  temp = matrix(data=0,nrow = 9, ncol = 3)
  rownames(temp) = colnames(x1)
  colnames(temp) = c("miRNA","miRNA-mRNA","mRNA")
  for(j in 1:9){
    a = c(x1[i,j],x3[i,j],x2[i,j])
    temp[j,which.min(a)]=1
  }
  data[[i]]=temp
}
overall = colMeans(xdata)
temp = matrix(data=0,nrow = 9, ncol = 3)
rownames(temp) = colnames(x1)
colnames(temp) = c("miRNA","miRNA-mRNA","mRNA")
for(j in 1:9){
  a = c(overall[j],overall[18+j],overall[9+j])
  temp[j,which.min(a)]=1
}
data[[6]]=temp

create.performance.plot.number <- function(r,data){
  df = data[[r]]
  library(reshape2)
  df = setNames(melt(df), c('rows', 'vars', 'values'))
  df$rows = as.character(df$rows)
  p <- ggplot(df,aes(x=vars,y=rows, shape = as.character(values)))+
    scale_shape_manual(values = c(32,15)) +
    geom_point(size = 3) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          #panel.background = element_blank(),
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.title = element_text(face = "bold", size = 10),
          legend.position="bottom",
          legend.box = "horizontal",
          #legend.box.background = element_rect(colour = "black",linetype="solid"),
          legend.text = element_text(face = "bold", size = 10))
  return (p)
}

gsn = list(length = 6)
for(r in 1:6) {
  gsn[[r]] = create.performance.plot.number(r,data)
}
pdf("Figure2-1.pdf",width=3,height=15,onefile = FALSE)
ggarrange(gsn[[5]], gsn[[4]],gsn[[3]],
          gsn[[2]],gsn[[1]],gsn[[6]],
          ncol = 1, nrow = 6,
          #labels = c("A", "B","C"),
          align = "hv",
          common.legend = TRUE, legend = "bottom")
dev.off()

###compare with PAM50
A = NULL
for(i in 1:(ncol(x2)-1))
{
  A = c(A, wilcox.test(x2[,i+1], x2[,1], alternative = "greater")$p.value)
}
