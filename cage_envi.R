############cage environmental data analysis#######
setwd("D:/")
#__________________________________________________________________
######## packages###############
#__________________________________________________________________
library("mgcv")
library("dplyr")
library("ggplot2")
library(GGally)
library(ggcorrplot)
library("factoextra")
library(FactoMineR)
library(cowplot)
library(dendextend)
library (vegan)
library (cluster)
library(NbClust)
library(circlize)

#__________________________________________________________________
######## data ###############
#__________________________________________________________________
orig.data<-read.csv("cage_envi_data.csv")
summary(orig.data)
View(orig.data)
boxplot(orig.data[,15])

#__________________________________________________________________
######## correlation among variables ###############
#__________________________________________________________________

str(orig.data)
corr <- round(cor(orig.data[,c(5:15)]),2)

#jpeg("corr_variables.jpg", res = 600,height = 6,width = 7,units = "in")
ggcorrplot(corr,tl.cex = 12,lab_size = 5, lab = F, digits = 2,legend.title = "Correlation",
           type="lower",outline.col = "white")
dev.off()

str(orig.data)

#jpeg("corr_all1.jpg", res = 600,height = 11,width = 17,units = "in")
ggpairs(orig.data,
        columns = c("Temperature","Salinity","pH","DO","GPP","NPP","TSS","Ammonia",
                    "Phosphate","Nitrate","Chlorophyll","SDI"), 
        title = " ", upper = list(continuous = wrap("cor",size = 3)),
        lower = list(continuous = wrap("smooth",
                                       alpha = 0.3,
                                       size = 0.1))
,
 mapping = aes(color = Time))
dev.off()

#__________________________________________________________________
######## GAM model fitting ###############
#__________________________________________________________________
set.seed(123)
gam.tem<-gam(SDI~s(Temperature),data =orig.data)
set.seed(123)
gam.sal<-gam(SDI~s(Salinity),data =orig.data)
set.seed(123)
gam.ph<-gam(SDI~s(pH),data =orig.data)
set.seed(123)
gam.do<-gam(SDI~s(DO),data =orig.data)
set.seed(123)
gam.gpp<-gam(SDI~s(GPP),data =orig.data)
set.seed(123)
gam.npp<-gam(SDI~s(NPP),data =orig.data)
set.seed(123)
gam.tss<-gam(SDI~s(TSS),data =orig.data)
set.seed(123)
gam.amm<-gam(SDI~s(Ammonia),data =orig.data)
set.seed(123)
gam.pho<-gam(SDI~s(Phosphate),data =orig.data)
set.seed(123)
gam.nit<-gam(SDI~s(Nitrate),data =orig.data)
set.seed(123)
gam.chl<-gam(SDI~s(Chlorophyll),data =orig.data)

summary(gam.chl)
AIC(gam.tem,gam.sal,gam.ph,gam.do,gam.gpp,gam.npp,gam.tss,gam.amm,gam.pho,gam.nit,gam.chl)
BIC(gam.tem,gam.sal,gam.ph,gam.do,gam.gpp,gam.npp,gam.tss,gam.amm,gam.pho,gam.nit,gam.chl)
plot(gam.ph)

jpeg("GAM plot_all.jpg", res = 600,height = 8,width =11,units = "in")
par(mfrow=c(3,4))
plot(gam.tem,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.sal,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.ph,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.do,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.gpp,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.npp,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.tss,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.amm,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.pho,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.nit,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.chl,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
dev.off()

#__________________________________________________________________
######## multi GAM model###############
#__________________________________________________________________
set.seed(123)
gam1<-gam(SDI~s(TSS)+s(Chlorophyll),data =orig.data)
set.seed(123)
gam2<-gam(SDI~s(TSS)+s(Chlorophyll)+s(Nitrate),data =orig.data)
set.seed(123)
gam3<-gam(SDI~s(TSS)+s(Chlorophyll)+s(Nitrate)+s(DO),data =orig.data)
set.seed(123)
gam4<-gam(SDI~s(TSS)+s(Chlorophyll)+s(Nitrate)+s(DO)+s(Ammonia),data =orig.data)
set.seed(123)
gam5<-gam(SDI~s(TSS)+s(Chlorophyll)+s(Nitrate)+s(DO)+s(Ammonia)+s(NPP),
          data =orig.data)
set.seed(123)
gam6<-gam(SDI~s(TSS)+s(Chlorophyll)+s(Nitrate)+s(DO)+s(Ammonia)+s(NPP)+
            s(pH),data =orig.data)
set.seed(123)
gam7<-gam(SDI~s(TSS)+s(Chlorophyll)+s(Nitrate)+s(DO)+s(Ammonia)+s(NPP)+
            s(pH)+s(Temperature),data =orig.data)
set.seed(123)
gam8<-gam(SDI~s(TSS)+s(Chlorophyll)+s(Nitrate)+s(DO)+s(Ammonia)+s(NPP)+
            s(pH)+s(Temperature)+s(Salinity),data =orig.data)
set.seed(123)
gam9<-gam(SDI~s(TSS)+s(Chlorophyll)+s(Nitrate)+s(DO)+s(Ammonia)+s(NPP)+
            s(pH)+s(Temperature)+s(Salinity)+s(GPP),data =orig.data)
set.seed(123)
gam10<-gam(SDI~s(TSS)+s(Chlorophyll)+s(Nitrate)+s(DO)+s(Ammonia)+s(NPP)+
            s(pH)+s(Temperature)+s(Salinity)+s(GPP)+s(Phosphate),data =orig.data)

AIC(gam1,gam2,gam3,gam4,gam5,gam6,gam7,gam8,gam9,gam10)

summary(gam8)
plot(gam9)

#__________________________________________________________________
######## PCA based GAM ###############
#__________________________________________________________________
set.seed(123)
gam.do1<-gam(Shannon.diversity~s(DO),data =orig.data)
set.seed(123)
gam.gpp1<-gam(Shannon.diversity~s(GPP),data =orig.data)
set.seed(123)
gam.npp1<-gam(Shannon.diversity~s(NPP),data =orig.data)
set.seed(123)
gam.tss1<-gam(Shannon.diversity~s(TSS),data =orig.data)
set.seed(123)
gam.amm1<-gam(Shannon.diversity~s(Ammonia),data =orig.data)
set.seed(123)
gam.pho1<-gam(Shannon.diversity~s(Phosphate),data =orig.data)
set.seed(123)
gam.nit1<-gam(Shannon.diversity~s(Nitrate),data =orig.data)
summary(gam.nit1)
AIC(gam.do1,gam.gpp1,gam.npp1,gam.tss1,gam.amm1,gam.pho1,gam.nit1)

#jpeg("GAM plot_all_pca.jpg", res = 600,height = 8,width =11,units = "in")
par(mfrow=c(2,4))
plot(gam.tss1,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.amm1,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.do1,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.pho1,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.npp1,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.gpp1,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
plot(gam.nit1,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue" ,scale=0, cex.lab=1.5, cex.axis=1.5,ylim = c(-0.2,0.3))
dev.off()

set.seed(123)
gam.101<-gam(Shannon.diversity~s(TSS)+s(Ammonia),data =orig.data)
set.seed(123)
gam.102<-gam(Shannon.diversity~s(TSS)+s(Ammonia)+s(DO), data =orig.data)
set.seed(123)
gam.103<-gam(Shannon.diversity~s(TSS)+s(Ammonia)+s(DO)+s(Phosphate), data =orig.data)
set.seed(123)
gam.104<-gam(Shannon.diversity~s(TSS)+s(Ammonia)+s(DO)+s(Phosphate)+s(NPP), data =orig.data)
set.seed(123)
gam.105<-gam(Shannon.diversity~s(TSS)+s(Ammonia)+s(DO)+s(Phosphate)+s(NPP)
             +s(GPP), data =orig.data)
set.seed(123)
gam.106<-gam(Shannon.diversity~s(TSS)+s(Ammonia)+s(DO)+s(Phosphate)+s(NPP)
             +s(GPP)+s(Nitrate), data =orig.data)
AIC(gam.101,gam.102,gam.103,gam.104,gam.105,gam.106)

summary(gam.106)

##GAM ploting for best model
#jpeg("GAM plot_gam106_.jpg", res = 600,height = 8,width =10,units = "in")
par(mfrow=c(2,4))
plot(gam.106,shade=TRUE,seWithMean=TRUE,shade.col = "skyblue", scale=0, cex.lab=1.5, cex.axis=1.5, ylim = c(-0.7,0.4))
dev.off()

#GAM diagnostic plots
jpeg("GAM diagnostic plot_gam15.jpg", res = 600,height = 8,width =10,units = "in")
par(mfrow=c(2,2))
gam.check(gam.106, cex.lab=1.5, cex.axis=1.5)
dev.off()

#GAM Partial affect plot
jpeg("GAM_partial.gam106_tss_do.jpg", res = 600,height = 8,width =8,units = "in")

vis.gam(gam.106, view=c("TSS","DO"), color = "topo", type = 'response',
        plot.type = 'persp',ticktype="detailed",
        xlab = "TSS", ylab = "DO", zlab = "SDI", zlim = c(2.5,4.0), cex.lab=1.15, cex.axis=1.15, 
        cex.main =2, main="", phi=30, theta=45, n.grid = 500, border=NA)
dev.off()


#__________________________________________________________________
######## PCA data arrangement ###############
#__________________________________________________________________

#for all data   
pca.data.all<-orig.data[,c(2:3,5:16)]#data filter
pca.data.all$site<-c(rep("Cage",120),rep("Reference",20)) # add site vari. to the data
str(pca.data.all)
data.all.sca<-scale(pca.data.all[,c(-1,-2,-15)], center = T) 

#########filter time data
####Before culture
str(orig.data)
pca.data.b<- orig.data[,c(-1,-4)] %>% filter(Time=="Before Culture") 
pca.data.b$site<-c(rep("Cage",48),rep("Reference",8)) # add site vari. to the data
#View(pca.data.b)
str(pca.data.b)

####During culture
str(orig.data)
pca.data.d<- orig.data[,c(2:3,5:16)] %>% filter(Time=='During Culture') 
pca.data.d$site<-c(rep("Cage",48),rep("Reference",8)) # add site vari. to the data
#####View(pca.data.d)
str(pca.data.d)

####After culture
str(orig.data)
pca.data.a<- orig.data[,c(2:3,5:16)] %>% filter(Time=='After Culture') 
pca.data.a$site<-c(rep("Cage",24),rep("Reference",4)) # add site vari. to the data
#View(pca.data.a)
str(pca.data.a)

########### All data PCA
#add the site to data.all.sca
data.all.sca1<-cbind(data.all.sca,pca.data.all["site"])
#write.csv(data.sca1,"data.sca1.csv")
str(data.all.sca1)

#__________________________________________________________________
######## PCA ###############
#__________________________________________________________________
pca.all<-PCA(data.all.sca1[,-13], ncp = 12)
pca.all
eigenvalues1 <- pca.all$eig
head(eigenvalues1[, 1:2])

fviz_screeplot(pca.all, ncp=10)

#jpeg("varianceplot_all.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_screeplot(pca.all, addlabels = TRUE, ylim = c(0, 30),
               main="",barfill = "#00AFBB",barcolor = "#00AFBB")
dev.off()

# Contributions of variables to PC1
vc1<-fviz_contrib(pca.all, choice = "var", axes = 1, top = 10, ylim = c(0, 25),
                  fill = "#00AFBB",color = "#00AFBB")+ 
  geom_hline(yintercept = 10, linetype = "dashed", color = "brown") +
  labs(title = "", x = "Variable",y = "Contribution")

# Contributions of variables to PC2
vc2<-fviz_contrib(pca.all, choice = "var", axes = 2, top = 10,ylim = c(0, 25),
                  fill = "#00AFBB",color = "#00AFBB")+ 
  geom_hline(yintercept = 10, linetype = "dashed", color = "brown") +
  labs(title = "", x = "Variable",y = "Contribution")

#jpeg("VCplot_all.jpg", res = 600,height = 6,width = 9,units = "in")
plot_grid(vc1,vc2,
          labels=c("(a)","(b)"),label_size=10, ncol=2)
dev.off()

#jpeg("biplot_all.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_pca_var(pca.all, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             title="")
dev.off()


fviz_pca_biplot(pca.all,repel = T ,habillage=as.factor(data.all.sca1$site))

#jpeg("biplot_all1.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_pca_biplot(pca.all, repel = T,habillage=as.factor(data.all.sca1$site),
                addEllipses=TRUE, ellipse.level=0.95)
dev.off()

jpeg("biplot_all2.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_pca_ind(pca.all,
             label = "none", # hide individual labels
             habillage = as.factor (data.all.sca1$site), # color by groups
             #palette = c("#00AFBB", "#E7B800", "#FC4E07","#0000FF","blue1","red1"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.level=0.95
)
dev.off()

########### Before culture PCA
#scale transformation for PCA
str(pca.data.b)
data.b.sca<-scale(pca.data.b[,c(-1:-2,-15)], center = T) 
#add the site to data.all.sca
data.b.sca1<-cbind(data.b.sca,pca.data.b["site"])
#write.csv(data.sca1,"data.sca1.csv")
str(data.b.sca1)
set.seed(123)
pca.b<-PCA(data.b.sca1[,-13], ncp = 12)

eigenvalues1 <- pca.b$eig
head(eigenvalues1[, 1:2])

fviz_screeplot(pca.b, ncp=10)

jpeg("varianceplot_before.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_screeplot(pca.b, addlabels = TRUE, ylim = c(0, 50),
               main="",barfill = "#00AFBB",barcolor = "#00AFBB")
dev.off()

# Contributions of variables to PC1
vc1<-fviz_contrib(pca.b, choice = "var", axes = 1, top = 10, ylim = c(0, 25),
                  fill = "#00AFBB",color = "#00AFBB")+ 
  geom_hline(yintercept = 10, linetype = "dashed", color = "brown") +
  labs(title = "", x = "Variable",y = "Contribution")

# Contributions of variables to PC2
vc2<-fviz_contrib(pca.b, choice = "var", axes = 2, top = 10,ylim = c(0, 25),
                  fill = "#00AFBB",color = "#00AFBB")+ 
  geom_hline(yintercept = 10, linetype = "dashed", color = "brown") +
  labs(title = "", x = "Variable",y = "Contribution")

#jpeg("VCplot_before.jpg", res = 600,height = 6,width = 9,units = "in")
plot_grid(vc1,vc2,
          labels=c("(a)","(b)"),label_size=10, ncol=2)
dev.off()

jpeg("biplot_before.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_pca_var(pca.b, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             title="")
dev.off()

fviz_pca_biplot(pca.b,repel = T ,habillage=as.factor(data.b.sca1$site))

#jpeg("biplot_b1.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_pca_biplot(pca.b, repel = T,habillage=as.factor(data.b.sca1$site),
                addEllipses=TRUE, ellipse.level=0.95)
dev.off()

#jpeg("biplot_b2.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_pca_ind(pca.b,
             label = "none", # hide individual labels
             habillage = as.factor (data.b.sca1$site), # color by groups
             #palette = c("#00AFBB", "#E7B800", "#FC4E07","#0000FF","blue1","red1"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.level=0.95
)
dev.off()


###########During culture PCA
#scale transformation for PCA
str(pca.data.d)
data.d.sca<-scale(pca.data.d[,c(-1:-2,-15)], center = T) 
#add the site to data.all.sca
data.d.sca1<-cbind(data.d.sca,pca.data.d["site"])
#write.csv(data.sca1,"data.sca1.csv")
data.d.sca1<-data.d.sca1[c(-1,-2,-9,-17,-25,-26,-33,-41),]
##extreme row values in the pca

str(data.d.sca1)
set.seed(123)
pca.d<-PCA(data.d.sca1[,-13], ncp = 12)

eigenvalues1 <- pca.d$eig
head(eigenvalues1[, 1:2])

fviz_screeplot(pca.d, ncp=10)

#jpeg("varianceplot_during.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_screeplot(pca.d, addlabels = TRUE, ylim = c(0, 30),
               main="",barfill = "#00AFBB",barcolor = "#00AFBB")
dev.off()

# Contributions of variables to PC1
vc1<-fviz_contrib(pca.d, choice = "var", axes = 1, top = 10, ylim = c(0, 25),
                  fill = "#00AFBB",color = "#00AFBB")+ 
  geom_hline(yintercept = 10, linetype = "dashed", color = "brown") +
  labs(title = "", x = "Variable",y = "Contribution")

# Contributions of variables to PC2
vc2<-fviz_contrib(pca.d, choice = "var", axes = 2, top = 10,ylim = c(0, 25),
                  fill = "#00AFBB",color = "#00AFBB")+ 
  geom_hline(yintercept = 10, linetype = "dashed", color = "brown") +
  labs(title = "", x = "Variable",y = "Contribution")

#jpeg("VCplot_during.jpg", res = 600,height = 6,width = 9,units = "in")
plot_grid(vc1,vc2,
          labels=c("(a)","(b)"),label_size=10, ncol=2)
dev.off()

#jpeg("biplot_during.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_pca_var(pca.d, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             title="")
dev.off()


fviz_pca_biplot(pca.d,repel = T ,habillage=as.factor(data.d.sca1$site))

#jpeg("biplot_d1.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_pca_biplot(pca.d, repel = T,habillage=as.factor(data.d.sca1$site),
                addEllipses=TRUE, ellipse.level=0.95)
dev.off()
data.d.sca1$site
#jpeg("biplot_d2.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_pca_ind(pca.d,
             label = "none", # hide individual labels
             habillage = as.factor (data.d.sca1$site), # color by groups
             #palette = c("#00AFBB", "#E7B800", "#FC4E07","#0000FF","blue1","red1"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.level=0.95
)
dev.off()

#######After culture PCA
#scale transformation for PCA
str(pca.data.a)
data.a.sca<-scale(pca.data.a[,c(-1:-2,-15)], center = T) 
#add the site to data.all.sca
data.a.sca1<-cbind(data.a.sca,pca.data.a["site"])
#write.csv(data.sca1,"data.sca1.csv")
data.a.sca1<- data.a.sca1[c(-49,-50),]

str(data.a.sca1)
set.seed(123)
pca.a<-PCA(data.a.sca1[,-13], ncp = 12)

#eigenvalues1 <- pca.b$eig
head(eigenvalues1[, 1:2])

fviz_screeplot(pca.a, ncp=10)

#jpeg("varianceplot_after.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_screeplot(pca.a, addlabels = TRUE, ylim = c(0, 40),
               main="",barfill = "#00AFBB",barcolor = "#00AFBB")
dev.off()

# Contributions of variables to PC1
vc1<-fviz_contrib(pca.a, choice = "var", axes = 1, top = 10, ylim = c(0, 40),
                  fill = "#00AFBB",color = "#00AFBB")+ 
  geom_hline(yintercept = 10, linetype = "dashed", color = "brown") +
  labs(title = "", x = "Variable",y = "Contribution")

# Contributions of variables to PC2
vc2<-fviz_contrib(pca.a, choice = "var", axes = 2, top = 10,ylim = c(0, 40),
                  fill = "#00AFBB",color = "#00AFBB")+ 
  geom_hline(yintercept = 10, linetype = "dashed", color = "brown") +
  labs(title = "", x = "Variable",y = "Contribution")

#jpeg("VCplot_after.jpg", res = 600,height = 6,width = 9,units = "in")
plot_grid(vc1,vc2,
          labels=c("(a)","(b)"),label_size=10, ncol=2)
dev.off()

#jpeg("biplot_after.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_pca_var(pca.a, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             title="")
dev.off()

fviz_pca_biplot(pca.a,repel = T ,habillage=as.factor(data.a.sca1$site))

#jpeg("biplot_a1.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_pca_biplot(pca.a, repel = T,habillage=as.factor(data.a.sca1$site),
                addEllipses=TRUE, ellipse.level=0.95)
dev.off()

#jpeg("biplot_a2.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_pca_ind(pca.a,
             label = "none", # hide individual labels
             habillage = as.factor (data.a.sca1$site), # color by groups
             #palette = c("#00AFBB", "#E7B800", "#FC4E07","#0000FF","blue1","red1"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.level=0.95
)
dev.off()

#__________________________________________________________________
######## Cluster analysis ###############
#__________________________________________________________________

#####data prepared from pca data
den.data.b<-pca.data.b
den.data.b$group<-c(rep("BCC",48),rep("BCR",8)) # add site vari. to the data
str(den.data.b)

den.data.d<-pca.data.d
den.data.d$group<-c(rep("DCC",48),rep("DCR",8)) # add site vari. to the data
str(den.data.d)
den.data.d<-den.data.d[c(-1,-2,-9,-17,-25,-33,-41,-49),]

den.data.a<-pca.data.a
den.data.a$group<-c(rep("ACC",48),rep("ACR",8)) # add site vari. to the data
str(den.data.a)
den.data.a<-den.data.a[c(-49,-50),]

######join all the group
dend.data<-rbind(den.data.b,den.data.d,den.data.a)
dend.data<-dend.data[,c(-1,-2,-15)]
rownames(dend.data)<-c(paste0("BCC", 1:48),paste0("BCR", 1:8),
                       paste0("DCC", 1:41),paste0("DCR", 1:7),
                       paste0("ACC", 1:48),paste0("ACR", 1:6))

######mean data for each group
meandend.data<-aggregate(dend.data[,1:12], by = list(group = dend.data$group), mean)
row.names(meandend.data)<-c("ACC","ACR","BCC","BCR","DCC","DCR")
str(meandend.data)

#remove the last column from data frame
###### Create a dendrogram and plot it
str(dend.data)
dend <- dend.data[-c(13,13)] %>%  scale %>% 
  dist %>% hclust %>% as.dendrogram
dend %>% plot

s<-dend.data[-c(13,13)] 

####Distance table created with bray methos for square root transformed data are
dis <- vegdist (sqrt (s), method = 'bray') # percentage cover data are transformed by square root

#### cluster using different method namely single, complete and average
cluster.single <- hclust (d = dis, method = 'single')
cluster.complete <- hclust (dis, 'complete')
cluster.average <- hclust (dis, 'average')

#####Plot the dendrogram for different method
par (mfrow = c (1,1)) # will draw all dendrogram into one figure

plot (cluster.single, main = 'Single linkage')
plot (cluster.complete, hang=-1,labels = dend.data$group, main = 'Complete linkage')
plot (cluster.average, hang=-1,labels = dend.data$group, main = 'average')
plot_horiz.dendrogram(cluster.average,  main = 'average')
# abline used to check no of cluster at perticular height
abline(h=0.025)

##save the dendrogram plot in to jpeg and label setting

#jpeg("winterdendro1.jpg", res = 600, height = 7,width =10,units = "in")
plot (cluster.average, 
      hang=-1, 
      labels = dend.data$group, 
      axes=FALSE,
      ylab = "Similarity",
      main = '', 
      sub = "",
      cex.lab= 1.25, 
      xlab = "",
      lwd=2, 
      cex=1.25)
axis(side=2, at=seq(0, max(cluster.average$height), length=6), labels=seq(100,0,-20),lwd=2) #replace the height lable with similarity
rect.hclust (cluster.average, h = 0.71, which = c(7,8), border = c('blue','red')) #h for hieght, which for tree no.
#rect.hclust (cluster.average, h = 0.815, which = c(7), border = 'green')
#rect.hclust (cluster.average, h = 0.87, which = c(2), border = 'purple')
dev.off()


plot (cluster.average, main = 'Average linkage')
rect.hclust (cluster.complete, h= 0.75, which = c(3,6,7)) 


clusters <- cutree (cluster.average, k = 7)
clusters
############cluster#####

cluster.flexible <- agnes (x = dis, method = 'average', par.method = 0.7)

cluster.flexible.hclust <- as.hclust (cluster.flexible)
plot (cluster.flexible.hclust, hang = -1,labels = dend.data$group, 
      main = 'Average linkage', xlab = 'cage',lwd=2, cex=1.25)
lines(x = c(0,0), y = c(0,100), type = "n") # force extension of y axis
axis(side = 2, at = seq(0,100,10), labels = seq(100,0,-10))

plot(silhouette(cutree(cluster.average,7),dis))


#varies technique for clustering
# Standardize the data
df <- scale(s)
head(df)
#The R code below determine the optimal number of clusters for k-means clustering:

# Elbow method
fviz_nbclust(df, kmeans, method = "wss") +
geom_vline(xintercept = 7, linetype = 2)+
labs(subtitle = "Elbow method")

set.seed(31)
# function to compute total within-cluster sum of squares
#jpeg("clust_EM.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_nbclust(df, kmeans, method = "wss", k.max = 15) + theme_minimal() + 
  ggtitle("the Elbow Method")+
  geom_vline(xintercept = 7,linetype = "dashed", color = "brown")
dev.off()

# Silhouette method
fviz_nbclust(df, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

fviz_nbclust(df, kmeans, method = "silhouette", k.max = 15) +
  theme_minimal() + ggtitle("The Silhouette Plot")

# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", nboot = 500)+
  labs(subtitle = "Gap statistic method")

gap_stat <- clusGap(df, FUN = kmeans, nstart = 30, K.max = 15, B = 500)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")

res.nbclust <- NbClust(df, distance = "euclidean",
                       min.nc = 2, max.nc = 9, 
                       method = "complete", index ="all")
factoextra::fviz_nbclust(res.nbclust) + theme_minimal() +
  ggtitle("NbClust's optimal number of clusters")

tmp <- NULL
for (k in 1:11){
  tmp[k] <- kmeans(df, k, nstart = 30)
}
df1 <- data.frame(tmp)
# add a prefix to the column names
colnames(df1) <- seq(1:11)
colnames(df1) <- paste0("k",colnames(df1))
# get individual PCA
df.pca <- prcomp(df1, center = TRUE, scale. = FALSE)
ind.coord <- df.pca$x
ind.coord <- ind.coord[,1:2]
df1 <- bind_cols(as.data.frame(df1), as.data.frame(ind.coord))
clustree(df1, prefix = "k")

#jpeg("clustree.jpg", res = 600,height = 11,width = 9,units = "in")
clustree(df1, prefix = "k")
dev.off()

#######group col and circle dendogram##########
d_s <- dist(s) # method="man" # is a bit better
hc_c <- hclust(d_s, method = "average")
str(d_s)
dend <- as.dendrogram(hc_c)

# Requires that the circlize package will be installed
# Color the branches based on the clusters:
dend<-dend %>% color_branches(k=7)%>%color_labels(k=7) #, groupLabels=iris_species)

par(mar = rep(0,4))
circlize_dendrogram(dend)

#jpeg("circle_dendogram.jpg", res = 600,height = 9,width = 9,units = "in")
circlize_dendrogram(dend)
dev.off()

plot (dend, 
      hang=-1)
#jpeg("vert_dendogram.jpg", res = 600,height = 17,width = 13,units = "in")
plot_horiz.dendrogram(dend)
dev.off()

d_s <- dist(s) # method="man" # is a bit better
hc_c <- hclust(d_s, method = "average")
library(dendextend)
dend <- as.dendrogram(hc_c)

some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))

#install.packages("gplots")
library(gplots)
gplots::heatmap.2(as.matrix(s), 
                  main = "",
                  srtCol = 20,
                  dendrogram = "row",
                  Rowv = dend,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(5,0.1),      
                  key.xlab = "Cm",
                  denscol = "grey",
                  density.info = "density",
                  col = some_col_func
)

########mean group dendogram#########
str(meandend.data)
meandend <- meandend.data[-c(1,1)] %>%  scale %>% 
  dist %>% hclust %>% as.dendrogram
meandend %>% plot

s_mean<-meandend.data[-c(1,1)] 

#Distance table created with bray methos for square root transformed data are
dis_mean <- vegdist (sqrt (s_mean), method = 'bray') # percentage cover data are transformed by square root

# cluster using different method namely single, complete and average
cluster.single_mean <- hclust (d = dis_mean, method = 'single')
cluster.complete_mean <- hclust (dis_mean, 'complete')
cluster.average_mean <- hclust (dis_mean, 'average')

#Plot the dendrogram for different method
par (mfrow = c (1,1)) # will draw all dendrogram into one figure

plot (cluster.single_mean, main = 'Single linkage')
plot (cluster.complete_mean, hang=-1,labels = meandend.data$group, main = 'Complete linkage')
plot (cluster.average_mean, hang=-1,labels = meandend.data$group, main = 'average')
plot_horiz.dendrogram(cluster.average_mean,  main = 'average')

#cluster analysis
# Standardize the data
df_mean <- scale(s_mean)
head(df_mean)

#jpeg("clust_EM_mean.jpg", res = 600,height = 7,width = 9,units = "in")
fviz_nbclust(df_mean, kmeans, method = "wss", k.max = 5) + theme_minimal() + 
  ggtitle("the Elbow Method")+
  geom_vline(xintercept = 3,linetype = "dashed", color = "brown")
dev.off()

#clustree
tmp1 <- NULL
for (k in 1:5){
  tmp1[k] <- kmeans(df_mean, k, nstart = 30)
}
df2 <- data.frame(tmp1)
# add a prefix to the column names
colnames(df2) <- seq(1:5)
colnames(df2) <- paste0("k",colnames(df2))
# get individual PCA
df.pca1 <- prcomp(df2, center = TRUE, scale. = FALSE)
ind.coord1 <- df.pca1$x
ind.coord1 <- ind.coord1[,1:2]
df2 <- bind_cols(as.data.frame(df2), as.data.frame(ind.coord1))
clustree(df2, prefix = "k")
#
#jpeg("clustree.mean.jpg", res = 600,height = 9,width = 7,units = "in")
clustree(df2, prefix = "k")
dev.off()

#dendogram
ds_mean <- dist(s_mean) # method="man" # is a bit better
hc_mean <- hclust(ds_mean, method = "average")

dend_mean <- as.dendrogram(hc_mean)

# Requires that the circlize package will be installed
# Color the branches based on the clusters:
dend_mean<-dend_mean %>% color_branches(k=3)%>%color_labels(k=3) #, groupLabels=iris_species)

par(mar = rep(0,4))
circlize_dendrogram(dend_mean)
#?circlize_dendrogram

#jpeg("circle_dendogram.mean.jpg", res = 600,height = 7,width = 7,units = "in")
circlize_dendrogram(dend_mean)
dev.off()

plot (dend_mean, 
      hang=-1,label=meandend.data$group)

#jpeg("vert_dendogram.mean.jpg", res = 600,height = 5,width = 7,units = "in")
plot_horiz.dendrogram(dend_mean)
dev.off()
#__________________________________________________________________
######## THE END ###############
#__________________________________________________________________
