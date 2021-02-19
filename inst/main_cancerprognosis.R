library(CancerSubtypesPrognosis)
library(reshape2)
library(ggplot2)
library(irr)
library(gridExtra)
library(Unicode)

###calculate the risk scores
dn = c("transbig", "unt", "upp", "mainz", "nki","GSE6532", "GEO", "TCGA753", "TCGA500",
       "METABRIC", "UK", "HEL", "GSE19783")
mRNA_methods <- c("AURKA", "ESR1", "ERBB2", "GGI", "GENIUS", "Endopredict", "OncotypeDx",
                  "TAMR13", "PIK3CAGS", "GENE70", "rorS", "Ensemble", "RNAmodel")
miRNA_methods <- c("hsa-miR-21", "hsa-miR-155", "hsa-miR-210", "miRNA10")
lncRNA_methods <- c("HOTAIR", "MALAT1", "DSCAM-AS1", "lncRNA12","lncRNA6","lncRNA5")

resMatrix <- as.list(NULL)

for (i in 1:length(dn)){
  print(dn[i])
  res = CancerPrognosis_RNAData(data=get(dn[i]), platform="custom", method=mRNA_methods)
  res = cbind(res,CancerPrognosis_miRNAData(data=get(dn[i]),method=miRNA_methods))
  res = cbind(res,CancerPrognosis_LncRNAData(data=get(dn[i]),method=lncRNA_methods))
  resMatrix[[i]] = res
}
names(resMatrix) = dn
save(resMatrix, file = "resMatrix.rda")

###calculate C-index
riskPList = c("AURKA", "ESR1", "ERBB2", "GGI", "GENIUS", "Endopredict", "OncotypeDx", "TAMR13",
              "PIK3CAGS", "GENE70", "rorS", "Ensemble","RNAmodel","miR-21","miR-155","miR-210",
              "miRNA10","HOTAIR","MALAT1","DSCAM-AS1","lncRNA12","lncRNA6","lncRNA5")
ciMatrix <- as.list(NULL)

for(i in 1:length(dn)){
  ddata=get(dn[i])
  sampleInfo= pData(ddata)
  if (dn[i] %in% c("transbig", "unt", "mainz", "nki")) {
    survival = data.frame(time=sampleInfo$t.dmfs,event=sampleInfo$e.dmfs, row.names=rownames(sampleInfo))
  }
  else if(dn[i] == "GSE19783") {
    survival=data.frame(time=sampleInfo$`disease free survival time (months):ch1`,event=sampleInfo$`death status:ch1`, row.names=rownames(sampleInfo))
  }
  else if(dn[i] == "HEL") {
    survival = data.frame(time=sampleInfo$`bddm followup time (months):ch1`,event=pd$`bddm status:ch1`, row.names=rownames(sampleInfo))
  }
  else {
    survival=data.frame(time=as.numeric(sampleInfo$t.rfs),event=as.numeric(sampleInfo$e.rfs), row.names=rownames(sampleInfo))
  }
  ciMatrix[[i]] = Cindex(data=resMatrix[[i]], survival, outputFolder=paste("./Cindex_output",dn[i],sep = "/"))
  colnames(ciMatrix[[i]]) = riskPList
}

names(ciMatrix) = dn
#save(ciMatrix, file = "ciMatrix.rda")

###calculate Overall scores for each dataset using mRNA-based methods, miRNA-based methods, lncRNA-based methods
for(i in names(ciMatrix)){
  #Get a meta-estimate
  ceData <- combine.est(x=ciMatrix[[i]]["cindex",1:13], x.se=ciMatrix[[i]]["cindex.se",1:13], hetero=TRUE, na.rm = T)
  cLower <- ceData$estimate + qnorm(0.025, lower.tail=TRUE) * ceData$se
  cUpper <- ceData$estimate + qnorm(0.025, lower.tail=FALSE) * ceData$se

  cindexO <- cbind("cindex"=ceData$estimate, "cindex.se"=ceData$se, "lower"=cLower, "upper"=cUpper)
  ciMatrix[[i]] <- cbind(ciMatrix[[i]], t(cindexO))

  ceData <- combine.est(x=ciMatrix[[i]]["cindex",14:17], x.se=ciMatrix[[i]]["cindex.se",14:17], hetero=TRUE, na.rm = T)
  cLower <- ceData$estimate + qnorm(0.025, lower.tail=TRUE) * ceData$se
  cUpper <- ceData$estimate + qnorm(0.025, lower.tail=FALSE) * ceData$se

  cindexO <- cbind("cindex"=ceData$estimate, "cindex.se"=ceData$se, "lower"=cLower, "upper"=cUpper)
  ciMatrix[[i]] <- cbind(ciMatrix[[i]], t(cindexO))

  ceData <- combine.est(x=ciMatrix[[i]]["cindex",18:23], x.se=ciMatrix[[i]]["cindex.se",18:23], hetero=TRUE, na.rm = T)
  cLower <- ceData$estimate + qnorm(0.025, lower.tail=TRUE) * ceData$se
  cUpper <- ceData$estimate + qnorm(0.025, lower.tail=FALSE) * ceData$se

  cindexO <- cbind("cindex"=ceData$estimate, "cindex.se"=ceData$se, "lower"=cLower, "upper"=cUpper)
  ciMatrix[[i]] <- cbind(ciMatrix[[i]], t(cindexO))
  colnames(ciMatrix[[i]]) <- c(riskPList, "mRNA_Overall","miRNA_Overall","lncRNA_Overall")
}

#save(ciMatrix, file = "ciMatrix.rda")

###calculate overall scores for each methods
method_overall = NULL
for(i in colnames(ciMatrix[[1]])){
  mp <- sapply(ciMatrix, function(x) { return(x[c("cindex","cindex.se"),i]) })

  ceData <- combine.est(x=mp["cindex",], x.se=mp["cindex.se",], hetero=TRUE, na.rm = T)
  cLower <- ceData$estimate + qnorm(0.025, lower.tail=TRUE) * ceData$se
  cUpper <- ceData$estimate + qnorm(0.025, lower.tail=FALSE) * ceData$se

  cindexO <- rbind("cindex"=ceData$estimate, "cindex.se"=ceData$se, "lower"=cLower, "upper"=cUpper)
  method_overall <- cbind(method_overall, cindexO)
}
colnames(method_overall) <- colnames(ciMatrix[[1]])
ciMatrix$method_overall <-method_overall
save(ciMatrix, file = "ciMatrix.rda")

###plot C-index matrix
data = NULL
for(i in colnames(ciMatrix[[1]])){
  mp <- sapply(ciMatrix, function(x) { return(x["cindex",i]) })
  data <- rbind(data, mp)
}
data = t(data)
colnames(data) <- colnames(ciMatrix[[1]])
rownames(data) <- names(ciMatrix)
rownames(data)[14] <- "Overall"
rownames(data) <- toupper(rownames(data))

data = data[,1:23]
refine = melt(as.matrix(data))
refine$Var1 = as.factor(refine$Var1)
refine$Var2 = as.factor(refine$Var2)
refine = na.omit(refine)
refine$value = round(refine$value, digits = 2)

pdf(file = "cindex.pdf",width = 9)
ggplot(refine, aes(Var2, Var1)) +
  geom_tile(aes(fill=value), color="white")+
  ggtitle("Concordance index") +
  xlab("Methods") + ylab("Datasets") +
  labs(fill = "value") +
  geom_text(aes(label = value), size=3 ) +
  theme_bw() +
  theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(color="Black", size=14, face="bold", hjust = 0.5),
        axis.title.x = element_text(color="Black", size=12, face="plain"),
        axis.title.y = element_text(color="Black", size=12, face="plain"),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  scale_fill_gradient(low = "white", high = "steelblue") +
  coord_equal()

dev.off()

###calculate the number of performace level for each methods
nres = matrix(data = NA, nrow = 4, ncol = ncol(data))
tdata = data[1:13,]
rownames(nres) = c("Poor", "Intermediate", "Good", "Total")
colnames(nres) = colnames(tdata)
nres["Poor",] = apply(tdata, 2, function(x) {return(sum(x<=0.5,na.rm = T))})
nres["Intermediate",] = apply(tdata, 2, function(x) {return(sum(x>0.5 & x<=0.7,na.rm = T))})
nres["Good",] = apply(tdata, 2, function(x) {return(sum(x>0.7,na.rm = T))})
nres["Total",] = apply(tdata, 2, function(x) {return(sum(x>=0,na.rm = T))})
write.csv(t(nres), file = "countCI.csv")

###get the index of the maximum C-index for each dataset
midx = apply(tdata, 1, function(x) {which.max(x)})
midx = colnames(tdata)[midx]
table(midx)

###compute PValues H0: cindex > 0.5
pv <- apply(ciMatrix$method_overall, 2, function(x) { return(pnorm((x[1] - 0.5) / x[2], lower.tail=x[1] < 0.5)) })
printPV <- matrix(pv,ncol=26)
rownames(printPV) <- "P-value"
colnames(printPV) <- names(pv)
printPV<-t(printPV)
## ----printPvalue,results="asis"------------------------------------------
write.csv(printPV,file = "Pvalue.csv",row.names = TRUE)

# Calculate cindex computation for each predictor
## cindex and p-value computation per algorithm
ciObjects <- as.list(NULL)

for(i in 1:length(dn)){
  ddata=get(dn[i])
  sampleInfo= pData(ddata)
  if (dn[i] %in% c("transbig", "unt", "mainz", "nki")) {
    survival = data.frame(time=sampleInfo$t.dmfs,event=sampleInfo$e.dmfs, row.names=rownames(sampleInfo))
  }
  else if(dn[i] == "GSE19783") {
    survival=data.frame(time=sampleInfo$`disease free survival time (months):ch1`,event=sampleInfo$`death status:ch1`, row.names=rownames(sampleInfo))
  }
  else if(dn[i] == "HEL") {
    survival = data.frame(time=sampleInfo$`bddm followup time (months):ch1`,event=pd$`bddm status:ch1`, row.names=rownames(sampleInfo))
  }
  else {
    survival=data.frame(time=as.numeric(sampleInfo$t.rfs),event=as.numeric(sampleInfo$e.rfs), row.names=rownames(sampleInfo))
  }
  ciObjects[[i]] = Cindex(data=resMatrix[[i]], survival, PValue = TRUE, outputFolder=paste("./Cindex_output",dn[i],sep = "/"))
}

names(ciObjects) = dn

tmp_ciobj = as.list(NULL)
for(i in 1:length(riskPList)){
  tmp_ciobj[[i]] <- sapply(ciObjects, function(x) { return(x[riskPList[i]]) })
}
names(tmp_ciobj) = riskPList

###Student t test
ccmData <- tt <- rr <- NULL
for(i in 1:length(tmp_ciobj)){
  tt <- NULL
  for(j in 1:length(tmp_ciobj)){
    if(i != j) {
      index <- c()
      for(k in 1:length(tmp_ciobj[[i]]))
        if (tmp_ciobj[[i]][[k]]$n != tmp_ciobj[[j]][[k]]$n)
          index <- c(index, k)

        index = unique(index)
        if (length(index) > 0) {
          rr <- cindex.comp.meta(list.cindex1=tmp_ciobj[[i]][-index],
                                 list.cindex2=tmp_ciobj[[j]][-index], hetero=TRUE)$p.value
        } else {
          rr <- cindex.comp.meta(list.cindex1=tmp_ciobj[[i]],
                                 list.cindex2=tmp_ciobj[[j]], hetero=TRUE)$p.value
        }
    }
    else { rr <- 1 }
    tt <- cbind(tt, rr)
  }
  ccmData <- rbind(ccmData, tt)
}
ccmData <- as.data.frame(ccmData)
colnames(ccmData) <- riskPList
rownames(ccmData) <- riskPList

# sapply(tmp_ciobj[[1]], function(x) { return(x$n)})

## ----printCCM,results="asis"---------------------------------------------
options(xtable.comment = FALSE)

## ----computeCCMPval------------------------------------------------------
ccmDataPval <- matrix(p.adjust(data.matrix(ccmData), method="holm"),ncol=length(riskPList),
                      dimnames=list(rownames(ccmData),colnames(ccmData)))

## ----printCCMPval,results='asis'-----------------------------------------
write.csv(ccmDataPval,file = "ccmDataPval.csv",row.names = TRUE)

####plot CCMPval
data = ccmDataPval

refine = melt(as.matrix(data))
refine$Var1 = as.factor(refine$Var1)
refine$Var2 = as.factor(refine$Var2)
refine = na.omit(refine)

# refine = rbind(refine, refine1)

refine$Rate = unlist(lapply(refine$value, function(x) {
  if (x <= 0.05) return("1")
  else return("2")
}))

pdf(file = "ccmDataPval.pdf",width = 9)
ggplot(refine, aes(Var1, Var2)) +
  geom_tile(aes(fill=Rate), color="white")+
  #geom_text(aes(label = value), size=1.5 ) +
  scale_fill_discrete(labels=c("Significant", "Not significant")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y= element_blank(), 
        axis.text.y = element_text(vjust=1, size=7, hjust=1),
        axis.text.x = element_text(angle=90, vjust=1, size=7, hjust=1))+
  # text = element_text(size=1))+
  coord_equal()
dev.off()


###calculate P-value of Log-rank test
pvMatrix <- NULL

for(i in 1:length(dn)){
  ddata=get(dn[i])
  sampleInfo= pData(ddata)
  if (dn[i] %in% c("transbig", "unt", "mainz", "nki")) {
    survival = data.frame(time=sampleInfo$t.dmfs,event=sampleInfo$e.dmfs, row.names=rownames(sampleInfo))
  }
  else if(dn[i] == "GSE19783") {
    survival=data.frame(time=sampleInfo$`disease free survival time (months):ch1`,event=sampleInfo$`death status:ch1`, row.names=rownames(sampleInfo))
  }
  else if(dn[i] == "HEL") {
    survival = data.frame(time=sampleInfo$`bddm followup time (months):ch1`,event=pd$`bddm status:ch1`, row.names=rownames(sampleInfo))
  }
  else {
    survival=data.frame(time=as.numeric(sampleInfo$t.rfs),event=as.numeric(sampleInfo$e.rfs), row.names=rownames(sampleInfo))
  }
  pvMatrix = rbind(pvMatrix, LogRank(data=resMatrix[[i]], survival))
}

rownames(pvMatrix) = dn

write.csv(pvMatrix,file = "LogrankPV.csv",row.names = TRUE)

### binarized the risk scores
groupMatrix = NULL
for(i in 1:length(resMatrix)) {
  groupMatrix = rbind(groupMatrix,apply(resMatrix[[i]], 2, function(x) {return(binarize(x, na.rm = FALSE))}))
}

###calculate kappa value
Kacoeff = NULL
PV = NULL
for(i in 1:ncol(groupMatrix)){
  temp_ka = NULL
  temp_pv = NULL
  for (j in 1:ncol(groupMatrix)) {
    temp = kappa2(cbind(groupMatrix[,i],groupMatrix[,j]))
    temp_ka = c(temp_ka,temp$value)
    temp_pv = c(temp_pv,temp$p.value)
  }
  Kacoeff = rbind(Kacoeff,temp_ka)
  PV = rbind(PV,temp_pv)
}
rownames(Kacoeff) = riskPList
colnames(Kacoeff) = riskPList
rownames(PV) = riskPList
colnames(PV) = riskPList

write.csv(Kacoeff,file = "kappa.csv",row.names = TRUE)
write.csv(PV,file = "kappa_PV.csv",row.names = TRUE)

###plot P-value of Log-rank test
#data = read.csv("logrank_data.csv", row.names = 1, header = T)
data = t(pvMatrix)
colnames(data) = toupper(colnames(data))
rownames(data) = gsub("\\.", "-", rownames(data))

refine = melt(as.matrix(data))
refine$Var1 = as.factor(refine$Var1)
refine$Var2 = as.factor(refine$Var2)
refine = na.omit(refine)
#refine$Rate = as.factor(ifelse(refine$value<=0.05, "Sig", "NoS"))
refine$value = ifelse(refine$value > 0.05, 0.051, refine$value)
refine$value = ifelse(refine$value <= 0.05 & refine$value > 0.02, 0.05, refine$value)
refine$value = ifelse(refine$value <= 0.02 & refine$value > 0.01, 0.02, refine$value)
refine$value = ifelse(refine$value <= 0.01, 0.01, refine$value)
head(refine)

p1 <- ggplot(refine, aes(Var1, Var2)) +
  ggtitle("a") + xlab("Method") + ylab("Dataset")+
  geom_point(aes(size=factor(value), colour=factor(value)),shape=20)+
  scale_colour_manual(name   = "P-value",
                      values=c("steelblue", "steelblue","steelblue", "#8E8A8B"),
                      breaks = c(0.01, 0.02, 0.05, max(as.numeric(refine$value))),
                      labels = c("\U2264 0.01", "\U2264 0.02", "\U2264 0.05", "> 0.05")) +
  scale_size_manual(name   = "P-value", 
                    values=c(5, 3.5, 2, 1),
                      breaks = c(0.01, 0.02, 0.05, max(as.numeric(refine$value))),
                      labels = c("\U2264 0.01", "\U2264 0.02", "\U2264 0.05", "> 0.05")) +
  theme_bw() +
  theme(
        plot.title = element_text(hjust = -0.1),
        legend.position = c("bottom"), # position the legend in the upper left 
        legend.direction = "horizontal",
        legend.title=element_text(size=10),
        legend.text = element_text(size = 9),
        axis.title.x = element_text(size = 10, face = "bold"), 
        axis.title.y= element_text(size = 10, face = "bold"), 
        axis.text.y = element_text(vjust=1, size=7, hjust=1), 
        axis.text.x = element_text(angle=90, vjust=1, size=7, hjust=1))

#ggsave("Fig01_discrete.pdf", useDingbats=FALSE, limitsize = FALSE)

###plot kappa value
temp = Kacoeff[lower.tri(Kacoeff)]
nres = rep(0, 6)
names(nres) = c("Not effective", "Slight", "Fair", "Moderate","Substantial","Almost perfect")
nres["Not effective"] = sum(temp<0)
nres["Slight"] = sum(temp>0&temp<=0.2)
nres["Fair"] = sum(temp>0.2&temp<=0.4)
nres["Moderate"] = sum(temp>0.4&temp<=0.6)
nres["Substantial"] = sum(temp>0.6&temp<=0.8)
nres["Almost perfect"] = sum(temp>0.8)

df <- data.frame(Concordance=names(nres),
                 Percentage=nres/length(temp))
df$Concordance = as.character(df$Concordance)

p2 <- ggplot(df, aes(x=factor(Concordance,level = c("Not effective", "Slight", "Fair", "Moderate","Substantial","Almost perfect")), y=Percentage))+
  ggtitle("b") + xlab("Concordance") + ylab("Percentage")+
  geom_bar(stat="identity",fill="steelblue",width=0.5) + 
  scale_y_continuous(labels=scales::percent)+
  theme_bw() +
  theme(
    plot.title = element_text(hjust = -0.1),
    axis.text.x = element_text(angle=45, size=7, vjust=1, hjust=1),
    axis.title.x = element_text(size = 10, face = "bold"), 
    axis.title.y= element_text(size = 10, face = "bold"), 
    axis.text.y = element_text(vjust=1, size=7, hjust=1)) +
  geom_text(aes(label=nres), size = 2.5,
            position=position_dodge(width=0.9), vjust=-0.25)

cairo_pdf(file = "Fig4.pdf",width = 8, height = 4,family="DejaVu Sans")
# layout.matrix <- matrix(c(1, 3, 2, 4), nrow = 1, ncol = 2)
# 
# layout(mat = layout.matrix,
#        heights = c(3, 2), # Heights of the two rows
#        widths = c(2, 2)) # Widths of the two columns

grid.arrange(p1, p2, nrow = 1, ncol = 2, widths = 3:2)
dev.off()

#data = read.csv("Kappa.csv", row.names = 1, header = T)
# rownames(data) = gsub("-", ".", rownames(data))

data = Kacoeff
data[upper.tri(data, diag = FALSE)] <- NA

refine = melt(as.matrix(data))
refine$Var1 = as.factor(refine$Var1)
refine$Var2 = as.factor(refine$Var2)
refine$value = round(refine$value, digits = 3)
refine = na.omit(refine)


# refine = rbind(refine, refine1)

# refine$Rate = unlist(lapply(refine$value, function(x) {
#   if (x <= 0) return("1")
#   if (x <= 0.2) return("2")
#   if (x <= 0.4) return("3")
#   if (x <= 0.6) return("4")
#   if (x <= 0.8) return("5")
#   if (x <= 1) return("6")
# }))

pdf(file = "FigS3.pdf",width = 10, height = 10)
ggplot(refine, aes(Var2, Var1)) +
  geom_tile(aes(fill=value), color="white")+
  xlab("Method") + ylab("Method") +
  labs(fill = "Kappa\ncoefficient")+
  geom_text(aes(label = value), size=3 ) +
  theme_bw() +
  theme(
    legend.position=c(0.8, 0.3),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(color="Black", size=14, face="bold"), 
        axis.title.y = element_text(color="Black", size=14, face="bold"), 
        axis.text.x = element_text(angle=45, vjust=1,hjust=1))+
  # text = element_text(size=1))+
  scale_fill_gradient(low = "grey91", high = "steelblue") +
  coord_equal()
dev.off()

###Overall forest plot
pdf(file = "forestPlot.pdf",width = 7, height = 7)
par(mfrow=c(2,2),mai=c(0.5,0.2,0.5,0.2))
layout.matrix <- matrix(c(1, 3, 2, 4), nrow = 2, ncol = 2)

layout(mat = layout.matrix,
       heights = c(3, 2), # Heights of the two rows
       widths = c(2, 2)) # Widths of the two columns
#for(i in c("method_overall","mRNA_Overall","miRNA_Overall","lncRNA_Overall")){
for(i in c("mRNA_Overall","miRNA_Overall","lncRNA_Overall")){
  if(i == "method_overall"){
    tt <- ciMatrix[[14]][,1:23]
    tt = t(tt)
    tt <- as.data.frame(tt)
    labeltext <- colnames(ciMatrix[[14]])[1:23]
    ttemp = cbind(labeltext, labeltext)

    r.mean <- c(tt$cindex)
    r.lower <- c(tt$lower)
    r.upper <- c(tt$upper)

    metaplot.surv(mn=r.mean, lower=r.lower, upper=r.upper, labels=labeltext, xlim=c(0.3,0.8),
                  boxsize=0.5, zero=0.5, cex = 0.8,
                  col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"),
                  main=paste(tolower(i)))
  }
  else {
    tt <- sapply(ciMatrix, function(x) { return(x[,i]) })
    tt = t(tt)
    tt <- as.data.frame(tt[1:13,])
    tt <- tt[complete.cases(tt), ]
    labeltext <- toupper(rownames(tt))

    r.mean <- c(tt$cindex)
    r.lower <- c(tt$lower)
    r.upper <- c(tt$upper)

    metaplot.surv(mn=r.mean, lower=r.lower, upper=r.upper, labels=labeltext, xlim=c(0.3,0.8),
                  boxsize=0.5, zero=0.5, cex = 0.8,
                  col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"),
                  main=paste(i))
  }
}
dev.off()

###Overall forest plot
pdf(file = "Fig3a.pdf",width = 8, height = 4,onefile = FALSE)
par(mfrow=c(1,2),mai=c(0.5,0.2,0.2,0))
layout.matrix <- matrix(c(2, 1), nrow = 1, ncol = 2)

layout(mat = layout.matrix,
       widths = c(5, 4)) # Widths of the two columns
layout.show()

#stPV = formatC(printPV, format = "e", digits = 2)
metaplot.surv(mn=r.mean, lower=r.lower, upper=r.upper, labels=labeltext, pvalue = stPV[1:23], xlim=c(0.3,0.8),
                 boxsize=0.5, zero=0.5, cex = 0.8,
                 col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"))
dev.off()
#countCI = read.csv(file = "Fig3.csv", header = T, stringsAsFactors = F)
# df = data.frame(row.names = countCI$X)
# df$x = countCI$Good
# df$y = countCI$Win

pdf(file = "Fig3b.pdf",width = 3.9, height = 3.9,onefile = FALSE)
ggplot(df, aes(x = x, y = y)) +
  ggtitle("b")+
  xlab("Good performance") + 
  ylab("Top performance")+
  geom_point(shape = 18, size = 4, color = "tomato4")+
  geom_text_repel(aes(label = rownames(df)), size = 2.5)+
  theme(
    plot.title = element_text(hjust = -0.1),
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 10),
    # panel.background = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.position = "none")
dev.off()

