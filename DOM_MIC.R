library(ggplot2)
library(RColorBrewer)


############# Shannon index 

dom<-read.csv(file = "heatmap_dom _all1.csv",header = T)

dom$sample <- factor(dom$sample, levels=c("MH","DXAL","WY", "CBS", "DLS","BTM","JGS","WYS","DHS","JFL","XSBN"), ordered=TRUE)
cols <-colorRampPalette(c("#8B8989", "#FF8C00", "#008B45", "#5D478B", "#436EEE", "#CDAD00"))


p1<-ggplot(dom,mapping = aes(x =sample, y =Shannon_dom,fill=sample,color=sample))+
  geom_boxplot(linetype="solid",size=1, width=0.5,alpha=0.6,outlier.shape = NA)+
  scale_fill_manual(name=" ",values = cols(11))+ 
  scale_color_manual(name=" ",values = cols(11))+# Fill colors           
  ylim(7.2,8.2)+
  labs(y=paste("DOM shannon index", sep=""))+
  geom_jitter(position=position_jitter(0.1),shape=21,alpha=0.8,size=3)+
  theme(plot.title = element_text(family = "serif",size = 20, colour = "black"),
        axis.text = element_text(family = "serif",size = 20, colour = "black"),
        axis.text.x=element_text(angle=30, hjust=1, vjust=1,family = "serif",  size = 18,color="black"),
        axis.title = element_text(family = "serif",size = 20 ,colour = "black"),
        axis.title.x = element_blank(),
        legend.text =  element_blank(),legend.position = "none",
        legend.title = element_text(family = "serif", size = 20,color="black"))
p1

############# Observed richness 

p2<-ggplot(dom,mapping = aes(x =sample, y =Richness_dom,fill=sample,color=sample))+
  #stat_boxplot(geom="errorbar",size=0.8,width=0.3,linetype="solid")+
  geom_boxplot(linetype="solid",size=1, width=0.5,alpha=0.6,outlier.shape = NA)+
  scale_fill_manual(name=" ",values = cols(11))+ 
  scale_color_manual(name=" ",values = cols(11))+# Fill colors           
  ylim(5000,7400)+
  labs(y=paste("DOM observed richness", sep=""))+
  
  geom_jitter(position=position_jitter(0.1),shape=21,alpha=0.8,size=3)+
  #geom_label(aes(label = "Bacteria"),size=6,x=10.5,y=8.0,family = "serif",label.size=NA,alpha=0,colour = "black")+
  theme_bw() +
  theme(plot.title = element_text(family = "serif",size =20, colour = "black"),
        axis.text = element_text(family = "serif",size = 20, colour = "black"),
        axis.text.x=element_text(angle=30, hjust=1, vjust=1,family = "serif",  size = 18,color="black"),
        axis.title = element_text(family = "serif",size = 20, colour = "black"),
        axis.title.x = element_blank(),
        legend.text =  element_blank(),legend.position = "none",
        legend.title = element_text(family = "serif", size = 20,color="black"))
p2
#############  NMDS analysis of DOM
library(vegan)
DOM <- read.delim(file = "otu_like_fill_intensity_NMDS.txt", header = TRUE,row.names = 1)
DOM <- data.frame(t(DOM))
hel <- decostand(DOM, method = 'hellinger')
group <- read.delim('design.txt', stringsAsFactors = FALSE)

DOM_Sample <- adonis(DOM~sample, group, distance = 'bray', permutations = 999)
DOM_Sample


bray_dis <- vegdist(DOM, method = 'bray')
nmds_dis <- metaMDS(bray_dis, k = 2)
nmds_dis$stress

nmds_dis_site <- data.frame(nmds_dis$points)
nmds_dis_site$ID <- rownames(nmds_dis_site)

site<-merge(nmds_dis_site,group,by="ID",all.x=TRUE)

library(ggplot2)
library(ggrepel)
library(ggthemes)
col1 <- colorRampPalette(c( "#FFE4C4","#B0E2FF", "#63B8FF",  "#FFAEB9","#FF3030" ))
p3<-ggplot(data =site,aes(x= MDS1, y=  MDS2,fill=pH ))+geom_tile(aes(fill =pH))+
  geom_point( shape=21,size=3,stroke = 0.5,alpha=0.8)+
  scale_fill_gradientn(name="",colours=col1(21),
                       limits=c(3,7),
                       guide = guide_legend())+
  
  geom_hline(yintercept=0,linetype=3,size=1)+ 
  geom_vline(xintercept=0,linetype=3,size=1)+
  labs(x=paste("NMDS1", sep=""),y=paste("NMDS2", sep=""))+
  annotate('text', label = paste('Stress =', round(nmds_dis$stress, 4)), x = 0.22, y = -0.22,family = "serif", size = 3.5, colour = 'black')+ 	#标注应力函数值
  geom_encircle(aes(fill=sample), alpha = 0.1, show.legend = F)+
  theme_bw()+
  theme(axis.title = element_text(family = "serif", face = "bold", size = 12,colour = "black"))+
  theme(axis.text = element_text(family = "serif", face = "bold", size = 12,color="black"))+
  theme(legend.text = element_text(family = "serif", face = "bold", size = 12,color="black"))+
  theme(legend.title = element_text(family = "serif", face = "bold", size = 12,color="black"))+
  theme(panel.grid=element_blank())
p3



############ Relative abundancce of DOM
data<- read.csv(file = 'MF_percent.csv',header=TRUE)
colors<-c("darkolivegreen3","gold","dodgerblue4","darkseagreen",
          
          "chartreuse4","darkorange","burlywood2","brown3","#984EA3","cyan3","grey50")

data$sample <- factor(data$sample, levels=c("MH","DXAL","WY", "CBS", "DLS","BTM","JGS","WYS","DHS","JFL","XSBN"), ordered=TRUE)
data$MF <- factor(data$MF, levels=c("Lignin","Tannin","Aromatic", "Protein", "Lipid","Carbohydrate","Others","UHCs"), ordered=TRUE)



p4<-ggplot(data, aes( x = sample, y = percent, fill = MF))+
  
  geom_bar(position = "fill", stat = "identity")+
  
  theme_bw()+
  
  scale_fill_manual(values=colors)+ 
  
  scale_y_continuous(expand = c(0,0))+
  
  labs(x="",y="Relative Abundance (%)",fill="")+
  
  theme(text=element_text(family = "serif", size = 20),
        
        axis.text.y=element_text(family = "serif", size = 20,color = "black"),
        
        axis.text.x=element_text(family = "serif", size = 16,color = "black",angle = 30, hjust =1, vjust = 1),
        
        legend.title=element_text(family = "serif", size = 20), 
        
        legend.text=element_text(family = "serif", size = 15))+ 
  
  theme(panel.grid = element_blank(),
        
        panel.background = element_rect(color = 'black', fill = 'transparent')) +
  
  guides(fill=guide_legend(keywidth = 0.8, keyheight = 1.5))+
  theme(plot.margin = unit(c(0.9,0.8,0.8,0.8),"cm"))
p4

##############  the multiple regression models
data=read.csv("DOM.csv", sep=",", header = T)
names(data)
str(data)


data1 <- data[,-1]
data1_z <- as.data.frame(scale(data1))
str(data1_z)
library(MuMIn)
library(performance)
library(tidyverse)
rm(list = ls())
 
model_B <- lm(Shannon_F~MAT+pH+CN+NP+Shannon_aromatic+
                Shannon_carbohydrate+Shannon_lignin+Shannon_lipid+
                Shannon_protein+Shannon_tannin,
              data=data)
summary(model_B)


plot(check_collinearity(model_B))

dre_model_SOC<- dredge(model_B, 
                       fixed = c("pH", "MAT"),
                       options(na.action = "na.fail"),
                       extra = c("R^2", F = function(x)
                         summary(x)$fstatistic[[1]]))
sub_model_SOC <- subset(dre_model_SOC, delta < 2)


aveg_model_SOC<-model.avg(dre_model_SOC, subset = delta < 2)
summary(aveg_model_SOC)

a1=summary(aveg_model_SOC)
a2 <- as.data.frame(a1$coefmat.full)
a3 <- a2[-1,c(1,2)]
a3$variables <- row.names(a3)
names(a3)[1] <- "coefficients"
names(a3)[2] <- "sd"
a3$proporation <- abs(a3$coefficients)/sum(abs(a3$coefficients))*100
a3
a3$sd <- a3$sd*1.96
a3

write.csv(a3,file = "LM_shannonF.csv")

mean_R2=mean(sub_model_SOC$`R^2`)
mean_AICc=mean(sub_model_SOC$`AICc`)
mean_delta_AICc=mean(sub_model_SOC$`delta`)
mean_R2
mean_AICc

library(tidyverse)
subgroup <- read.csv("subgroup.csv", header = T, sep = ",")


comb_date <- a3 %>% left_join (subgroup,
                               by = c('variables' = 'variables'))
comb_date <- comb_date %>% mutate(Year = "2006")
str(comb_date)


comb_date$variables <- factor(comb_date$variables,
                              levels=c("PCoA1_lipid","PCoA1_aromatic","PCoA1_lignin",
                                       "PCoA1_tannin","PCoA1_carbohydrate",
                                       "CN","NP",  "pH","MAT" ))

comb_date$variables <- factor(comb_date$variables,
                              levels=c("Shannon_lipid","Shannon_aromatic","Shannon_tannin",
                                       "Shannon_protein","Shannon_lignin","Shannon_carbohydrate",
                                       "CN", "pH","MAT" ))

figure3a <- ggplot(comb_date,aes(coefficients,variables,fill=subgroup))+
  
  geom_errorbar(aes(xmin=coefficients-sd, xmax=coefficients+sd,colour =subgroup ), size= 0.8, width=.01)+
  
  geom_point(size=8,shape=22,color="black",stroke = 1)+
  scale_fill_manual(values=c( "#91bcdd", "#e4babB", "#e5c185"))+
  scale_color_manual(values=c( "#91bcdd", "#e4babB", "#e5c185"))+
  geom_vline(aes(xintercept=0),linetype=2.5,col="black")+
  scale_x_continuous(limits = c(-1.4,1.2))+
  scale_y_discrete(position = 'right')+
  labs(y=NULL)+
  labs(x=NULL)+
  theme_classic()+theme( axis.line.x = element_line(size = 1, color = "black"),  # x 轴线宽度设置为 2，颜色为蓝色
                         
                         axis.text.x = element_text(color='black',size=20),
                         legend.text =  element_text(size = 20,color = "black"),
                         axis.text.y = element_text(color='black',size=20),
                         axis.line.y =element_blank(),
                         axis.ticks.y =element_blank()) +
  theme(legend.position="none")
figure3a



comb_date$subgroup <- factor(comb_date$subgroup,levels=c("Climate",
                                                         "Soil",
                                                         "DOM"))
figure3b <- ggplot(data = comb_date,aes(Year, proporation, fill=subgroup))+
  geom_bar(stat="identity", size=5)+
  scale_fill_manual(values=c( "#91bcdd", "#e5c185", "#e4babB"))+
  ggtitle(expression("R "^2==0.75))+
  theme_classic()+
  scale_y_continuous(name = "Relative effect of estimates (%)")+
  annotate("text", x = 1, y = 78, label = paste("Climate"),  color = "black") +
  annotate("text", x = 1, y = 50, label = paste("Soil Properties"), color = "black") +
  annotate("text", x = 1, y = 21, label = paste("DOM Chemodiverisity"),color = "black") +
  
  theme(axis.ticks.x = element_blank(),
        axis.line.y = element_line(size = 1, color = "black") ,
        axis.title.y = element_text( size = 20,colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(size = 1, color = "black"),  # 调整刻度线粗细和颜色
        axis.text.y =element_text(color='black',size=20),
        
        axis.line.x =element_blank())+
  theme(legend.position="none")+
  labs(x=NULL) +
  theme(
    plot.title = element_text(
      size = 20,      
      hjust = 0.5,    
      color = "black"  
    )
  )
figure3b


###################  Structural equation modeling

data=read.csv("DOM.csv", sep=",", header = T)
names(data)
str(data1)
data1 <- data
data1$Shannon_aromatic <- log(data1$Shannon_aromatic)
data1$Shannon_carbohydrate <- log(data1$Shannon_carbohydrate)
data1$Shannon_lipid <- log(data1$Shannon_lipid)
data1$Shannon_tannin <- log(data1$Shannon_tannin)
data1$Shannon_lignin <- log(data1$Shannon_lignin)

## 数据Z转化，标准化的数据用来进行线性模型选择，随机森林可以用原始数据
data1 <- data[,-1]
data1_z <- as.data.frame(scale(data1))
str(data1_z)



DATA.NA<-data1_z
DATA.NA<-data1

colnames(DATA.NA)

###################  shannonF ~~ DOM diveristy
mod1 <- psem(
  lm(pH~ MAT , data=DATA.NA),
  lm(CN~ MAT +pH, data=DATA.NA),
  lm(NP~ MAT, data=DATA.NA),
  lm(Shannon_aromatic ~MAT+ pH+NP , data=DATA.NA),
  lm(Shannon_lipid~ MAT+pH, data=DATA.NA),
  lm(Shannon_F~  Shannon_carbohydrate+Shannon_lipid+Shannon_aromatic+CN+MAT
     , data=DATA.NA),
   data = DATA.NA
)

################### Shannon_B ~~ DOM diveristy
mod1 <- psem(
  lm(pH~ MAT , data=DATA.NA),
  lm(CN~ MAT +pH, data=DATA.NA),
  lm(NP~ MAT, data=DATA.NA),
  lm(Shannon_aromatic ~MAT+ pH+NP , data=DATA.NA),
  lm(Shannon_lipid~ MAT+pH, data=DATA.NA),
  lm(Shannon_B~  Shannon_lipid+CN+NP+pH
     , data=DATA.NA),
  data = DATA.NA
)


############ Phylogenetic tree

library(pacman)
library(ggsci)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(ggplot2)
library(treeio)


tree <- read.newick("otus.treefile")
tax<- read.delim("Fotu_annon.txt",header = T)
group <-  read.delim("Fotu_annon.txt",header = T,row.names=1)
tmap<-read.delim("FD_heatmap.txt",header = T)

groupInfo <- split(row.names(group), group$Class )

tree <- groupOTU(tree, groupInfo)


col1 <- colorRampPalette(c("#0000CD", "#FFFFFF", "#EE2C2C"))#定义颜色
col2<-colorRampPalette(c("#CD5B45", "#6495ED", "#698B69", "#EEAD0E", "#C7C7C7", "#CD6090",
                         "#00868B", "#CDBE70", "#551A8B", "#EE9A49", "#006400", "#EEC900",
                         "#00008B", "#8B0000", "#4A708B", "#CD5555", "#00008B", "#8B7D7B",
                         "#FFAEB9", "#AB82FF", "#FA8072", "#7B68EE"))
p<-ggtree(tree, layout="fan", open.angle=15,size=0.1, linetype=1)+
  geom_tippoint(shape=21,alpha=0.6, starstroke=0.1,size=2,aes(fill= group))+
  geom_hilight(node=683, fill="gray",alpha=0.5)+
  geom_hilight(node=1300, fill="gray",alpha=0.5)+
  geom_hilight(node=1314, fill="gray",alpha=0.5)+
  geom_hilight(node=1337, fill="gray",alpha=0.5)+
  geom_hilight(node=869, fill="gray",alpha=0.5)+
  geom_hilight(node=806, fill="gray",alpha=0.5)+
  geom_hilight(node=1049, fill="gray",alpha=0.5)+
  geom_hilight(node=1153, fill="gray",alpha=0.5)+
  geom_hilight(node=1262, fill="gray",alpha=0.5)+
  geom_hilight(node=1267, fill="gray",alpha=0.5)+
  geom_hilight(node=1276, fill="gray",alpha=0.5)+
  geom_hilight(node=1284, fill="gray",alpha=0.5)+
  geom_hilight(node=735, fill="gray",alpha=0.5)

p 
p1 <- p+ new_scale_fill()+
  geom_fruit(data=tmap,
             geom=geom_tile, 
             mapping=aes(y=label, fill=r_dom), 
             alpha=0.8, 
             width=0.3,
             offset=0.08
  )+ 
  geom_fruit(data=tmap,
             geom=geom_tile,
             mapping=aes(y=label, fill=r_ar),
             alpha=0.8, 
             width=0.3, 
             offset=0.05
  )+
  
  geom_fruit(data=tmap,
             geom=geom_tile, 
             mapping=aes(y=label, fill=r_ca), 
             alpha=0.8, 
             width=0.3, 
             offset=0.05 
  )+ 
  geom_fruit(data=tmap,
             geom=geom_tile, 
             mapping=aes(y=label, fill=r_lig),
             alpha=0.8, 
             width=0.3,
             offset=0.05
  )+
  geom_fruit(data=tmap,
             geom=geom_tile, 
             mapping=aes(y=label, fill=r_lip), 
             alpha=0.8, 
             width=0.3, 
             offset=0.05
  )+
  geom_fruit(data=tmap,
             geom=geom_tile,
             mapping=aes(y=label, fill=r_pro), 
             alpha=0.8, 
             width=0.3, 
             offset=0.05 
  )+ 
  geom_fruit(data=tmap,
             geom=geom_tile, 
             mapping=aes(y=label, fill=r_tan), 
             alpha=0.8, 
             width=0.3, 
             offset=0.05 
  )+scale_fill_gradientn(limit=c(-1,1),breaks=c(-1,-0.5,0,0.5,1),colours=col1(21),
                         name="",guide=guide_legend(keywidth = 0.5, 
                                                    keyheight = 0.5, order=2))+
  theme( 
    legend.background=element_rect(fill=NA),
    legend.title=element_text(size=6.5),
    legend.text=element_text(size=4.5),
    legend.spacing.y = unit(0.01, "cm")
  )



p1
p2<-open_tree(p1, 30) %>% rotate_tree(105)

p2


p3<-p2+geom_cladelabel(node=1337,label="Rozellomycotina_cls_Incertae_sedis", hjust =1,angle = 0,fontsize=2, offset=5.8,color="black",barsize = 0.8)+
  geom_cladelabel(node=1300,label="Agaricomycetes",hjust = 1,angle = 0, offset=4.6,fontsize=2,color="black",barsize = 0.8)+
  geom_cladelabel(node=1314,label="Agaricomycetes",hjust = 1,angle = 0, offset=5.3,fontsize=2,color="black",barsize =0.8)+
  geom_cladelabel(node=685,label="Mortierellomycetes",hjust = 1,angle = 0, offset=4.75,fontsize=2,color="black",barsize = 0.8)+
  geom_cladelabel(node=735,label="Agaricomycetes",hjust = 0.8,angle = 0, offset=3.8,fontsize=2,color="black",barsize = 0.8)+
  geom_cladelabel(node=1284,label="Agaricomycetes", hjust = -0.2,angle = 0,offset=4.1,fontsize=2,color="black",barsize = 0.8)+
  geom_cladelabel(node=798,label="Sordariomycetes",hjust = 1,angle = 0, offset=5.25,fontsize=2,color="black",barsize = 0.8)+
  geom_cladelabel(node=1276,label="Sordariomycetes",hjust = 1,angle = 0, offset=5.25,fontsize=2,color="black",barsize = 0.8)+
  geom_cladelabel(node=1267,label="Eurotiomycetes",hjust = 1,angle = 0, offset=4.65,fontsize=2,color="black",barsize = 0.8)+
  geom_cladelabel(node=1259,label="Pezizomycetes",hjust = 1,angle = 0, offset=5.1,fontsize=2,color="black",barsize = 0.8)+
  geom_cladelabel(node=1262,label="Pezizomycetes",hjust = 1,angle = 0, offset=5.1,fontsize=2,color="black",barsize = 0.8)+
  geom_cladelabel(node=1153,label="Leotiomycetes",hjust = 1,angle = 0, offset=4.65,fontsize=2,color="black",barsize = 0.8)+
  geom_cladelabel(node=1049,label="Eurotiomycetes",hjust = 1,angle = 0, offset=3.3,fontsize=2,color="black",barsize = 0.8)+
  geom_cladelabel(node=806,label="Dothideomycetes",hjust = 1,angle = 0, offset=3.8,fontsize=2,color="black",barsize = 0.8)+
  geom_cladelabel(node=869,label="Sordariomycetes",hjust = 1,angle = 0, offset=2.5,fontsize=2,color="black",barsize = 0.8)


p3
