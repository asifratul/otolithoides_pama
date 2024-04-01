########## Import packages and fonts ##############

library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(HH)
library(likert)
library(ggrepel)
library(extrafont)
library(ggpubr)

font_import()
loadfonts(device = "win")
windowsFonts("Serif" = windowsFont("Serif"))


##################### Fig3: RSCU Usage ######################

dataset2 <- read_csv("Dataset/RSCU_Pama.csv")
dataset <- dataset2[c(1:14)] %>% column_to_rownames("Codons") %>% as.matrix

annotation_row = data.frame(Codon = dataset1$Codons, AA = dataset1$Aminoacids...16)

rownames(annotation_row) <- annotation_row$Codon
annotation_row1 <- annotation_row %>% dplyr::select(c(0,2)) 

H1 = pheatmap(dataset, cluster_rows =FALSE, main = "RSCU in 13 PCGs of Pama Pama",border_color=FALSE,
         annotation_row = as.data.frame(annotation_row1),display_numbers=T, number_format = "%.1f",
         treeheight_row = 0, treeheight_col = 0,fontsize_number = 12, angle_col = "0")


ggsave("RSCU_fig.jpg",height = 8, width = 11, units = "in", dpi=300)


##################### Fig 4(a) Conserved/Variable ######################

df3 <- read_csv("Dataset/Pama_ConserveVariables.csv")

df3 <- df3[,c(1,5,6)]

library(ggpubr)

P1 <- plot.likert(Gene~., df3,
                  auto.key=list(space="right", columns=1,reverse=TRUE, padding.text=3),
                  main='', xlab="Percentages", values = "show", show.prc.sign = TRUE,
                  xlim=c(-70,110), scales=list(x=list(at=seq(-70,100,10), labels=c(seq(70,0,-10),seq(10,100,10)))),
                  col=c("#E94E1B", "#F7AA4E"))

P1

png("ConservedPercentages_17June23.png",height = 11, width = 8, units = "in", res = 360)

P1
dev.off()

######################## 4(b) K2P distance ######################


df2 <- read_csv("Dataset/K2P_distance.csv")


P2 = df2 %>%
  ggplot(aes(x = Gene, y = Distance, fill = Type))+
  geom_bar(stat="identity") + scale_fill_brewer(palette="Set2")+
  geom_errorbar(aes(ymin = Distance - S.E., ymax = Distance + S.E.), width=0.3, position = position_dodge(0.9))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Different genes\n") + ylab("K2P Distance\n")+
  theme(axis.title = element_text(size = 16), legend.title = element_blank(),legend.text = element_text(size = 16))+
  theme(text=element_text(size=16, family="Serif"))+ 
  theme(axis.text.x = element_text(face = "italic"))+
  theme(legend.position = c(0.8, 0.85))

P2

ggsave("K2P.jpg",height = 7, width = 10, units = "in", dpi=600)


######################## 4(c) Ka/Ks###########################


df3 <- read_csv("Dataset/Ka_Ks.csv")


P3 = df3 %>%
  ggplot(aes(x = Ka, y = Ks,color = PCG))+
  geom_label_repel(aes(label=PCG), fontface="bold.italic", size=5)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 1)) +
  theme_bw() + xlab("Ka rates\n") + ylab("Ks rates\n")+
  theme(text=element_text(size=16, family="Serif"))+
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.title = element_blank())+ theme(legend.position = "none")+
  theme(axis.text=element_text(size=16), axis.title = element_text(size = 16,family="Serif"))

P3
ggsave("Ka_Ks_17June23.jpg",height = 7, width = 10, units = "in", dpi=600)

ggarrange(P1,ggarrange(P2, P3, ncol = 1, labels = c("B", "C")),
          nrow = 1, labels = "A")


############ Fig. S3, AT-GC Skew##################

df4 <- read_csv("Dataset/AT_GCSKew.csv")


df4 %>%
  ggplot(aes(x = Skew, y = Total, color=Type))+
  geom_point(size=3)+
  facet_wrap(~Total_Type, ncol =2)+
  scale_x_continuous(limits = c(-0.7, 0.7),breaks = c(-0.60, -0.4, -0.2, 0, 0.2, 0.4, 0.60))+
  theme_bw() + xlab("\nSkewness") + ylab("Content(%)\n")+
  theme(text=element_text(size=16, family="Serif"))+
  theme_bw()+
  theme(axis.text=element_text(size=16), axis.title = element_text(size = 16)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 12)) +
  theme(strip.text = element_text(size = 16)) +
  theme(strip.background = element_blank(), strip.placement = "outside")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text=element_text(family="serif"))+
  theme(strip.background = element_blank(), strip.placement = "outside")+
  theme(axis.text=element_text(size=14), axis.title = element_text(size = 16,family="Serif"))

ggsave("AT_GC_SKEW.jpg",height = 6, width = 10, units = "in", dpi=600)

############## Figure S5: Hierchial agglomerative clustering#################

library(ape)
library(cluster)

d1 <- read_csv("RSCU_PCGGenomes.csv")
d2 <- d1 %>% remove_rownames %>% column_to_rownames(var="CODONS")
distance_mat <- dist(d2, method = 'euclidean')
distance_mat


## average###
Hierar_cl <- hclust(distance_mat, method = "average")
Hierar_cl

jpeg("PCG_genomes_Algorithm-average.jpeg", width = 14, height = 8, units = 'in', res = 600)
plot(as.phylo(Hierar_cl), label.offset = 0.05)
dev.off()