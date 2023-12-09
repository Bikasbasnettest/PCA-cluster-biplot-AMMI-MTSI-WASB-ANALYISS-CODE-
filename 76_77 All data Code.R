setwd("E:/NARC1")
Bikas<-read.csv("E:/NARC1/R-76-77.csv", header = TRUE)
print(Bikas)
head(Bikas)
library(metan)
Elimin<-scale(Bikas[-1 : -3], center = TRUE)
Elimin
all<-corr_coef(Bikas)
all
plot(all)
Ex<-all$cor
Ex<-as.data.frame(Ex)
class(Ex)
library(writexl)
write_xlsx(Ex, "Correlation value of Maize 76-77.xlsx")
getwd()
P1<-network_plot(
  all,
  min_cor = NULL,
  show = c("signif", "all"),
  p_val = 0.05,
  legend = c("full", "range"),
  colours = c("red", "white", "darkgreen"),
  legend_width = 1,
  legend_height = 15,
  legend_position = c("right", "left", "top", "bottom"),
  curved = TRUE,
  angle = 90,
  curvature = 0.5,
  expand_x = 0.25,
  expand_y = 0.25
)
P1
ggsave(filename = "Network Correlation Plot.jpg", plot = P1,
       width = 35, height = 30, dpi = 1500, units = "cm")

library(readxl)
CC76_77<-read_excel("Correlation value of Maize 76-77.xlsx")
View(CC76_77)
library(sjPlot)
head(CC76_77)
tab_corr(CC76_77)
#to show the lower triangle only
tab_corr(CC76_77,
         triangle = "lower", 
         p.numeric = TRUE,
         file = "76-77 maize significant Correlation Neumeric.docs")
require(ggplot2)
library(reshape2)
Ex1<-all$cor
Ex1<-as.data.frame(Ex1)
class(Ex1)
corr<-data.frame(Ex1)
corr_df<-as.data.frame(corr)
corr_df$variable1<-rownames(corr_df)
melted_corr<-melt(corr_df, id.vars = "variable1", variable.name = "variable2", value.name = "value")
P2<-ggplot(melted_corr, aes(x=variable1, y=variable2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", mid = "white", high = "green", 
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  coord_polar() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 9.5),
        axis.text.y = element_text(size = 9))

P2
library(ggplot2)
ggsave(filename = "Circular Correlation Plot.jpg", plot = P2,
       width = 35, height = 30, dpi = 1500, units = "cm")

#Note for creation of the Dendro gram mean value of the each Genotype is Requires
Bikas1<-read.csv("E:/NARC1/Cluster gRpah and Blup model for 76-77/MG_76-77.csv")
attach(Bikas1)
library(NbClust)
Elimin1<-scale(Bikas1[, -1])
Elimin1
?dist
distance<-dist(Elimin1)
distance
#elbow method of calculation of the Optimum number of clusters
library(factoextra)
library(ggplot2)
fviz_nbclust(Elimin1, kmeans, method ="wss")
fviz_nbclust(Elimin1, kmeans, method ="silhouette")
#Best Method is NB cluste method 
require(NbClust)
?NbClust
NbClust(data = Elimin1, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 15, 
        
        method = "kmeans", index = "all", alphaBeale = 0.1)
#now
setk<-kmeans(Elimin1, 2)
setk$cluster
plot(Elimin1,col=setk$cluster)
#for hiearchial clustering
NbClust(data = Elimin1, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 15, 
        
        method = "ward.D", index = "all", alphaBeale = 0.1)
#this indicates the D value 
NbClust(data = Elimin1, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 15, 
        
        method = "ward.D2", index = "all", alphaBeale = 0.1)
wardse<-hclust(distance, "ward.D2")
plot(wardse)
plot(wardse, hang = -1)
#after clearing the all Plots
plot(wardse,  labels = Elimin1$GEN)
P1<-plot(wardse,hang = -1, labels = Bikas1$GEN)
print(distance,digits=2)
hc<-hclust(distance)
P1<-plot(hc,labels=Bikas1$GEN)
rect.hclust(hc, k=11, border = "purple") 
rect.hclust(hc, k=11, border = 1:11)

#FOR Creation of Circulr Dendrogram
library(circlize)
library(dendextend)
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = 11)  
dend <- set(dend, "branches_lwd", 2)
par(mar = c(1))  
circlize_dendrogram(dend, labels_track_height = 0.1) 
genotype_names <- as.character(Bikas1$GEN)  
labels(dend) <- genotype_names
par(mar = c(3, 3, 3, 3))  
P2<-circlize_dendrogram(dend, labels_track_height = 0.1)
P2


#For PCA analysis
library(factoextra)
library(FactoMineR)
str(Bikas)

#singular value decomposition
pca2<-prcomp(Elimin)
pca2
pca2$x
var=get_pca_var(pca2)
corr_plot(var$cos2)
summary(pca2)
biplot(pca2)
Mpca<-PCA(Elimin, ncp=24)
Mpca
Mpca$eig
table1<-Mpca$eig
class(table1)
table1<-as.data.frame(table1)
library(writexl)
write_xlsx(table1, "elgien value of the PCA 76-77.xlsx")
plot(table1$`cumulative percentage of variance`)

fviz_pca(pca2)
PCA<-fviz_pca_biplot(
  pca2,
  axes = c(1, 2),
  geom = c("point", "text"),
  geom.var = c("arrow", "text"),
  col.ind = "darkgreen",
  fill.ind = "red",
  col.var = "black",
  fill.var = "red",
  gradient.cols = NULL,
  label = "var",
  invisible = "none",
  repel = TRUE,
  habillage = "none",
  palette = NULL,
  addEllipses = TRUE,
  title = "PCA - Biplot",
)
PCA
ggsave(filename = "76-77 PCA plot.jpg", plot = PCA,
       width = 35, height = 30, dpi = 1500, units = "cm")

#ggbiplot(pca2, ellipse= TRUE, choices= c(1,2),labels= rownames(GGE), groups=GGE$ENV)

library(factoextra)
library(cluster)  
library(ggplot2)
library(dplyr)
pca_result <- prcomp(Elimin, scale. = TRUE)  
num_clusters <- 3  
kmeans_result <- kmeans(pca_result$x[, 1:2], centers = num_clusters)
cluster_assignments <- kmeans_result$cluster
cluster_colors <- c("red", "blue", "green") 
individual_colors <- cluster_colors[cluster_assignments]
biplot<-fviz_pca_biplot(
  pca_result,
  geom = c("point", "text"),
  geom.var = c("arrow", "text"),
  col.ind = individual_colors,  
  fill.ind = "white",
  col.var = "black",
  fill.var = "white",
  gradient.cols = NULL,
  label = "var",
  invisible = "none",
  repel = TRUE,
  habillage = "none",
  palette = NULL,
  addEllipses = TRUE,
  title = "PCA - Biplot"
)
biplot


fviz_pca_var(pca2,repel = TRUE)
Contri<-fviz_pca_var(pca2, col.var = "contrib",
                     repel = TRUE,
                     gradient.cols = c("black", "darkgreen", "red"),
                     ggtheme = theme_minimal())

Contri
SP1<-fviz_eig(pca2, addlabels = TRUE, barfill="darkgreen")
SP1
SP2<-fviz_screeplot(pca2, choice="eigenvalue", ncp=10, addlabels = TRUE, barfill="brown")
SP2
library(gridExtra)
grid.arrange(SP1,SP2, ncol=1, top='76-77 Plot of Eigen value and % of explain variance')
#too see the contribution of Particular contribution by variables
a<-fviz_contrib(pca2, choice = "var",fill = "lightgreen", axes = 1, top = 12)
a
a2<-fviz_contrib(pca2, choice = "var", axes = 2,fill = "#40E0D0",top = 9)
a2
a3<-fviz_contrib(pca2, choice = "var", axes = 3,fill = "#007BFF",
                 color = "green", top = 6)
a3
a4<-fviz_contrib(pca2, choice = "var", axes = 4,fill = "blue",
                 color = "pink", top = 4)
a4

grid.arrange(a,a2, ncol=2, top='Contribution of variables to First 2PCs')


#for analyis of the AMMI model, stablity and other aspects
library(metan)
options(max.print = 500)
require(readxl)
Bikas$ENV<-factor(Bikas$ENV, levels=unique(Bikas$ENV))
Bikas$GEN <-factor(Bikas$GEN, levels=unique(Bikas$GEN))
Bikas$REP<-factor(Bikas$REP, levels=unique(Bikas$REP))
Bikas$YEAR<- factor(Bikas$YEAR, levels=unique(Bikas$YEAR))
str(Bikas)
inspect(Bikas, plot=TRUE)
desc_stat(Bikas, stats="all")
ds<-desc_stat(Bikas, stats="all") 
ds
tab1<-(ds)
class(tab1)
tab1<-as.data.frame(tab1)
library(writexl)
write_xlsx(tab1, "Discriptive statstics of Genotypes Traits 76-77.xlsx")
print(Bikas)
find_outliers(Bikas, var=Y.ha, plots=TRUE)
find_outliers(Bikas, var=DTS, plots=TRUE)
remove_rows_na(Bikas)
replace_zero(Bikas)
find_text_in_num(Bikas$DTA)
find_text_in_num(Bikas$Y.ha)
#mg indicates the mean of the Genotypes

mg<-means_by(Bikas, GEN)
mg
me<-means_by(Bikas, ENV)
me
head(Bikas)
mge<-Bikas %>% 
  group_by(ENV, GEN) %>%
  desc_stat(RotL,ShtLod,DTA,DTS,TSI,PHT,EHT,
            E.P,P.Ha,EH.Ha,ProLi,HC,TUR,MAY,INSc,PA,
            EA,CoD,CoL,SP,GP.R,GR.C,TGW,Y.ha, stats="mean")
mge
View(mge)
tabl1<-(mg)
tabl2<-(mge)
class(tabl1)
class(tabl2)
tab1<-as.data.frame(tabl1)
tabl2<-as.data.frame(tabl2)
library(writexl)
write_xlsx(tabl1, "Mean of Genotypes 76-77.xlsx")
write_xlsx(tabl2, "Mean of Genotypes across environment 76-77.xlsx")
pyld<-ge_plot(Bikas, ENV, GEN, Y.ha)
pyld
?ge_plot
pyld2<-ge_plot(Bikas, ENV, GEN, Y.ha, type=2)
pyld2
DTA<-ge_plot(Bikas, ENV, GEN, DTA, type=2)
DTA
DTS<-ge_plot(Bikas, ENV, GEN, DTS, type=2)
DTS
#alternative code
ge_plot(
  Bikas,
  ENV,
  GEN,
  Y.ha,
  type = 2,
  values = TRUE,
  text_col_pos = c("top", "bottom"),
  text_row_pos = c("left", "right"),
  average = TRUE,
  order_g = NULL,
  order_e = NULL,
  xlab = NULL,
  ylab = NULL,
  width_bar = 5,
  heigth_bar = 15,
  plot_theme = theme_metan(),
  colour = TRUE
)
ggsave(filename = "Yield performance of Genotype across the Env.jpg", plot = pyld2,
       width = 35, height = 30, dpi = 1500, units = "cm")
ggsave(filename = "Days to anthesis of Genotype across the Env.jpg", plot = DTA,
       width = 35, height = 30, dpi = 1500, units = "cm")
ggsave(filename = "Days to silking of Genotype across the Env.jpg", plot = DTS,
       width = 35, height = 30, dpi = 1500, units = "cm")
#winner Genotype across the Environment
win<-ge_winners(Bikas, ENV, GEN, resp = everything())
win
POLAR<-ge_polar(Bikas, ENV, GEN, resp = everything(),base = 10, verbose = TRUE)
POLAR


?ge_polar
win1<-ge_winners(Bikas, ENV, GEN, resp = everything(), type = "ranks")
win1
View(win)
tbl1<-(win)
tabes1<-(win1)
class(tbl1)
class(tabes1)
tab1<-as.data.frame(tbl1)
tabes1<-as.data.frame(tabes1)
library(writexl)
write_xlsx(tbl1, "all aspects winner Genotypes all env 76-77.xlsx")
write_xlsx(tbl1, "all aspects winner Gen Rank all env 76-77.xlsx")
#for fixed effect model
options(max.print = 10000)
inanova<-anova_ind(Bikas, ENV, GEN,REP, resp = everything())
inanova
#anova for Yield
iava<-inanova$Y.ha$individual
iava
View(iava)
#for polled anova analyissis
Panova<-anova_joint(Bikas, ENV, GEN,REP, resp = everything())
Panova
#for individual variables of the Pooled analysis 
PAV1<-Panova$Y.ha$anova
PAV1
tab2<-(PAV1)
tab2<-as.data.frame(tab2)
class(tab2)
write_xlsx(tab2, "Pooled anova analyis of yield all env 76-77.xlsx")
#for the stablity analyiss
#for the calculation of the anicardo anova analysis
anna<-Annicchiarico(Bikas, ENV,GEN,REP,Y.ha)
anna

#for the stablity analyiss
?Annicchiarico
View(anna$Y.ha$environments)
#for the calculation of the ecovalence value of the data
Eco<-ecovalence(Bikas, ENV,GEN,REP,Y.ha)
Eco
tab3<-(Eco)
tab3<-as.data.frame(tab3)
class(tab3)
write_xlsx(tab3, "Ecovalence method(Yield) all env 76-77 a.xlsx")
#for the stablity analyiss of sukla anaysisis
shu1<-Shukla(Bikas, ENV,GEN,REP,Y.ha)
shu1
#Reg_based stablity  of the genotypes for the yiled assocated traits.
reg1<-ge_reg(Bikas, ENV,GEN,REP,Y.ha)
reg1
plot(reg1)
View(reg1)
#anova for reganalyiss of data
reganv1<-reg1$Y.ha$anova
reganv1
View(reganv1)
write_xlsx(reganv1, "regression anova for the Genotypes.xlsx")
plotA<-plot(reg1)
plotA
ggsave(filename = "ReEnv for yield.jpg", plot = plotA,
       width = 35, height = 30, dpi = 1500, units = "cm")
#supperior by Lin e Binns' superiority index method of the non parametric test
head(Bikas)
?superiority
library(metan)
?Fox
out<-Fox(Bikas, ENV, GEN, Y.ha)
print(out)
#for factor bAsed
fact1<-ge_factanal(Bikas, ENV, GEN,REP, Y.ha)
fact1
####3wrap the stablity parameters
?ge_stats
stab1<-ge_stats(Bikas,ENV,GEN,REP, Y.ha, verbose = TRUE, prob = 0.05)
stab1
#3#ammi model analysisis##3
amod1<-performs_ammi(Bikas, ENV, GEN, REP, Y.ha)
print(amod1)
amiplot<-plot(amod1)
amiplot
ggsave(filename = "AMMI plot of Yield.jpg", plot = amiplot,
       width = 35, height = 30, dpi = 1500, units = "cm")
View(amod1$Y.ha$ANOVA)
write_xlsx(amod1$Y.ha$ANOVA,"ammianova1.xlsx")
#significance value of ipca
?get_model_data
get_model_data(amod1, "ipca_pval")
#to get ammiplot
a1<-plot_scores(amod1)
a1
A1<-plot_scores(amod1, x.lab = "Yield/ha")
A1
b1<-plot_scores(amod1, type = 2)
b1
Bb1<-plot_scores(amod1, type = 2, polygon = TRUE)
Bb1
B1<-plot_scores(amod1, 
                type = 2,
                col.env = "blue",
                col.gen = transparent_color(),
                col.segm.env = "orange",
                highlight = c("MBS-1144", "SUPER6768"),
                col.highlight = "darkcyan",
                axis.expand = 1.5)
B1
c1<-plot_scores(amod1, type = 4)
c1
C1<-plot_scores(amod1, type=4, repulsion = 2)
C1
P1<-arrange_ggplot(a1, b1, tag_levels = "a", nrow = 1)
P1
A1+Bb1
ggsave(filename = "AMMI combine of Yield.jpg", plot = P1,
       width = 35, height = 30, dpi = 1500, units = "cm")
ggsave(filename = "AMMI combine of Yield vs PC.jpg", plot = A1+Bb1,
       width = 35, height = 30, dpi = 1500, units = "cm")

######################### ammi based stability statistics ##################
abs1<-ammi_indexes(amod1)
print(abs1)
ggemodel1<-gge(Bikas, ENV, GEN, Y.ha)
ggemodel1
GGE<-plot(ggemodel1)
GGE
options(max.print = 100000)
predict(ggemodel1)
bbp1<- plot(ggemodel1, col.gen = "brown")
bbp1
ggsave(filename = "GGE biplot of the 76-77.jpg", plot = bbp1,
       width = 35, height = 30, dpi = 1500, units = "cm")


#calculation of Discritivenss and Representiveness
dvr1<-plot(ggemodel1, type = 4, plot_theme = grey(level = 1))
dvr1
dvr1 + theme(plot.background = element_rect(fill = "green"))+
  theme(plot.title = element_text(size = rel(2)))+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
ggsave(filename = "Discritivenss and Representivenesse 76-77.jpg", plot = dvr1,
       width = 35, height = 30, dpi = 1500, units = "cm")


?theme
#ranking over environment of the Genotypes
rel<-plot(ggemodel1, type = 6)
RE<-plot(rel)
#relationship among environments
rae1<-plot(ggemodel1, type = 10)
RA<-plot(rae1)
ggsave(filename = "Relationship and Ranking the Environments 76-77.jpg", plot = RE+RA,
       width = 35, height = 30, dpi = 1500, units = "cm")


###mean performance vs Stablity
GGE1<-gge(
  Bikas,
  ENV,
  GEN,
  Y.ha,
  centering = "environment",
  scaling = "none",
  svp = "environment",
  by = NULL,
)

plot(GGE1)

#need table extraction from here 
gpg1<-gge(Bikas, ENV, GEN, Y.ha, svp = "genotypes")
pgpg<- predict(gpg1)
View(pgpg$Y.ha)
?theme
?gge
?repel
mvs1<-plot(gpg1, type = 2,repel = TRUE, size.text.gen = 2, repulsion = 1)
plot(mvs1)
#ranking of the Genotypes
rg1<-plot(gpg1, type = 8)
plot(rg1)
#single value partation for Symmetrical
gp1<-gge(Bikas, ENV, GEN, Y.ha, svp = "symmetrical")
pgps1<-predict(gp1)
View(pgps1$Y.ha)
#which won where
ww1<-plot(gp1, type = 3)
WWW<-plot(ww1)
ggsave(filename = "which won where.jpg", plot = WWW,
       width = 35, height = 30, dpi = 1500, units = "cm")

#examining the Particular environment(need figure extraction)
E1<-plot(gp1, type = 5,col.gen = "darkred",size.text.gen = 2, repel = TRUE, sel_env = "Rampur")
E1
View(Bikas)
E2<-plot(gp1, type = 5,repulsion = 1,size.text.gen = 2,size.line = 1,repel = TRUE,size.shape.win = 8, size.text.gen = 1.8, sel_env = "Parwanipur")
E2
E3<-plot(gp1, type = 5,size.text.gen = 2, repel = TRUE, sel_env = "Tarahara")
E3
E4<-plot(gp1, type = 5, repel = TRUE,size.text.gen = 2, sel_env = "Nepalgunj")
E4
E1+E2+E3+E4
#for the Blup modeling
# Example for removing missing values
# Check for missing values in each column
missing_values <- colSums(is.na(Bikas))

# Display columns with missing values

library(metan)
model<-gamem(Bikas, GEN, REP, resp = everything())
BLUP=gmd(model, "blupg")
BLUP
aku<-mgidi(model)
#for the extraction of the all data sets
table2<-aku$PCA
table3<-aku$FA
table4<-aku$MGIDI
table4<-aku$sel_dif
table5<-aku$contri_fac_rank
table6<-aku$contri_fac_rank_sel
class(table2)
class(table3)
class(table4)
class(table5)
class(table6)
table2<-as.data.frame(table2)
table3<-as.data.frame(table3)
table4<-as.data.frame(table4)
table5<-as.data.frame(table5)
table6<-as.data.frame(table6)
class(table2)
class(table3)
class(table4)
class(table5)
class(table6)
library(writexl)
write_xlsx(table2, "76-77 eigen value of .xlsx")
write_xlsx(table3, "76-77 FAmy mung bean.xlsx")
write_xlsx(table4, "76-77 differential of the Mung bean.xlsx")
write_xlsx(table5, "77-78 factor rank of the genotypes.xlsx")
write_xlsx(table6, "76-77contribution factor rank of the selected genotypes.xlsx")
getwd()
p1<-plot(aku, type = "contribution")
p1
p2<-plot(aku)
p2
p1+p2
library(ggplot2)
ggsave(filename = "76-77 Contribututon and selected strength and weakness.jpg", plot = p1+p2,
       width = 35, height = 30, dpi = 1500, units = "cm")



# Fitting the WAAS index
AMMI <- waasb(Bikas,
              env = ENV,
              gen = GEN,
              rep = REP,
              resp = c(everything()))

# Getting the weighted average of absolute scores
wasb<-gmd(AMMI, what = "WAASB")
wasb
was<-(wasb)
class(was)
was<-as.data.frame(was)
class(was)
library(writexl)
write_xlsx(was, "77-78 weighted averag value for WAASB index.xlsx")
#new method
GMD<-get_model_data(model, "lrt")
?get_model_data
GMD
plot(GMD)

BLP1<-plot_waasby(AMMI,
               var = 1,
               which = "gen",
               ncol = NULL,
               nrow = NULL,
               prob = 0.05,
               export = FALSE,
               file.type = "pdf",
               file.name = ,
               plot_theme = theme_metan(),
               width = 6,
               height = 6,
               err.bar = TRUE,
               size.err.bar = 0.5,
               size.shape = 3.5,
               size.tex.lab = 8,
               height.err.bar = 0.3,
               x.lim = NULL,
               x.breaks = waiver(),
               col.shape = c("darkgreen", "red"),
               y.lab = "Genotypes",
               x.lab = "Yield/ha",  # Specify your custom x-axis label here
               n.dodge = 1,
               check.overlap = FALSE,
               panel.spacing = 1
)
BLP1
ggsave(filename = "all genotype 76 and 77 wasb plot.jpg", plot = BLP1,
       width = 30, height = 25, dpi = 1500, units = "cm")

?waasb()
VCO<-get_model_data(model, "vcomp")
print(VCO, n=5000)


# Genetic parameters
genpar<-get_model_data(model, "genpar")
genpar
# random effects
ranef<-get_model_data(model, "ranef")
ranef
# Predicted values
predict(model)
?mgidi
P1<-plot(aku,
         type = "contribution",
         genotypes = "all",
         x.lab = "Treatments",
         width = 1,
         title = "The strengths and weaknesses view of treatments",
         rotate = TRUE)

P1
ggsave(filename = "all genotype 77 & 78 srength and weakness.jpg", plot = P1,
       width = 35, height = 30, dpi = 1500, units = "cm")

P1+p2
ggsave(filename = "Cracking Contribn and strengtha and weakness.jpg", plot = P1+p2,
       width = 35, height = 30, dpi = 1500, units = "cm")

#Genotype-environment analysis by mixed-effect models

GEAMDE<-gamem_met(Bikas,
                  env = ENV,
                  gen = GEN,
                  rep = REP,
                  resp = everything())
#f# Distribution of random effects (first variable)
plot(model, type = "re")
## Distribution of random effects overll aspects
plot(GEAMDE, type = "re")
#for the linear mixed model effectvanalyis 
library(lme4)
library(lmerTest)
library(agricolae)
library(car)
library(multcompView)
library(emmeans)
#for homogenety of the variance
str(Bikas)
Bikas$ENV=as.factor(Bikas$ENV)
Bikas$REP=as.factor(Bikas$REP)
Bikas$YEAR=as.factor(Bikas$YEAR)
leveneTest(Y.ha~ENV, data = Bikas)
leveneTest(Y.ha~YEAR, data = Bikas)
?leveneTest
#BLUP plot for Genotypes
library(ggplot2)
library(metan)
head(Bikas)
Mode<-gamem_met(Bikas, ENV, GEN,REP, Y.ha)
BLUP1<-plot(Mode)
BLUP1
ggsave(filename = "BLUP plot components 76-77.jpg", plot = BLUP1,
       width = 35, height = 30, dpi = 1500, units = "cm")
BLUP2<-plot(Mode, which = c(1,2,7), ncol = 1)
BLUP2
ggsave(filename = "BLUP plot components2 76-77.jpg", plot = BLUP2,
       width = 35, height = 30, dpi = 1500, units = "cm")
options(repr.plot.width = 13, repr.plot.height = 11)
D<- plot_blup(Mode)
D
arrange_ggplot(D, tag_levels = list(c("D")), ncol = 1) + coord_flip()
BLP<-plot_blup(Mode,
               var = 1,
               which = "gen",
               ncol = NULL,
               nrow = NULL,
               prob = 0.05,
               export = FALSE,
               file.type = "pdf",
               file.name = ,
               plot_theme = theme_metan(),
               width = 6,
               height = 6,
               err.bar = TRUE,
               size.err.bar = 0.5,
               size.shape = 3.5,
               size.tex.lab = 8,
               height.err.bar = 0.3,
               x.lim = NULL,
               x.breaks = waiver(),
               col.shape = c("darkgreen", "red"),
               y.lab = "Genotypes",
               x.lab = "Yield/ha",  # Specify your custom x-axis label here
               n.dodge = 1,
               check.overlap = FALSE,
               panel.spacing = 1
)
BLP
ggsave(filename = "BLUP plot for yield 76-77.jpg", plot = BLP,
       width = 35, height = 30, dpi = 1500, units = "cm")
