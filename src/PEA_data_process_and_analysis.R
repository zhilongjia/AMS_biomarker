# setwd("/Users/diudiu/Desktop/Figure_v4_20210827/01_PEA_DEPs/src")
library(dplyr)

# rawdata without CV<0.3, missing freqency < 0.25 and duplication.  If NPX < LOD and NPX> 0, replaced NPX with LOD/ sqrt；if NPX < LOD and NPX<0 ,replace NPX with LOD

pro_raw <- read.csv("../data/10_olink_rawdata.csv", sep = ",")
g <- read.csv("../data/10_olink_group.csv",sep = ",")
pro_raw2 <- pro_raw %>% select(-UniprotID) %>% 
  tibble::column_to_rownames(var ='symbol')
colnames(pro_raw2) <- g$sample
t_raw <- t(pro_raw2)
t_raw <- data.frame(t_raw, stringsAsFactors = F)

#repalce QC unqualified symbol with min value multiply by 0.1

raw <- data.frame(apply(t_raw,2,function(x){
  x[is.na(x)]= min(x, na.rm = T) *0.1
  return(x)
}),stringsAsFactors = F  ) 

temp=NULL
  for (j in 1:ncol(raw))
  {
   # j=1
    test_1=raw[1:10,j]
    test_2=raw[11:20,j]
      var=var.test(test_1,test_2)$p.value>0.05
      test_palue=t.test(test_1,test_2,paired = T,var.equal =var)
      test_FC= mean(test_2)-mean(test_1)
      log2FC = sprintf("%.4f", test_FC)
      temp=rbind(temp,c(paste("PEA_AMS4k_AMS1k"),colnames(raw)[j],test_palue$p.value,
                        paste("t_test; var equal=",var,sep=""), log2FC))
  }

temp2 <- data.frame(temp, stringsAsFactors = F)    
colnames(temp2) <- c("group",'symbol','p_value','method','log2FC')

q_value <- p.adjust(temp2$p_value, method = "BH")

PEA_test <- cbind(temp2, q_value)


dep_PEA <- PEA_test %>% filter(q_value < 0.05)

#replace min*0.1 rawdata merge

traw <- t(raw)
traw <- data.frame(traw, stringsAsFactors = F)
traw <- traw %>% tibble::rownames_to_column( var = 'symbol')
DEPs_raw <- dep_PEA %>% left_join(traw, by = 'symbol')
PEA_comb <- PEA_test %>% left_join(traw, by = 'symbol')

# rawdata contains NA 

save(PEA_comb, file = "../results/PEA_rawdata_result_combine.RData")

DEPs_raw2 <- dep_PEA %>% left_join(pro_raw, by = 'symbol')
colnames(DEPs_raw2)[8:27] <- g$sample
PEA_comb2 <- PEA_test %>% left_join(pro_raw, by = 'symbol')

write.csv(dep_PEA, file = "../results/PEA_DEPs_results_v2.csv")

save(DEPs_raw2, file = "../results/PEA_DEPs_rawdata_result_withNA.RData")
save(PEA_comb2, file = "../results/PEA_rawdata_result_combine_withNA.RData")


save.image("../results/PEA_test.RData")

########## Volcano plot of 887 proteins ###########

library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(magrittr)
library(dplyr)
# setwd("/Users/diudiu/Desktop/Figure_v4_20210827/04_olink_volc/src")

load("../data/PEA_rawdata_result_combine.RData")
data <- PEA_comb %>% select(-c(group, method)) %>% 
  tibble::column_to_rownames( var ='symbol')

data <- apply(data,2,as.numeric)
data <- data.frame(data, stringsAsFactors = F)
rownames(data) <- PEA_comb$symbol
save(data, file = "../results/PEA_comb_numeric.RData")

all <- data %>% tibble::rownames_to_column(var ='symbol') %>%
  subset(select = c("symbol","log2FC","q_value")) %>%
  filter(abs(log2FC) > 0.5 & q_value < 0.05)

save(all, file = "../results/olink_top15_DEPs.RData")

vol_raw <- data %>% select(-c(p_value)) 

EnhancedVolcano (vol_raw,
                 lab = rownames(vol_raw),
                 col = c("grey30","grey30","grey30","red2"),
                 x = 'log2FC',
                 y = 'q_value',
                 selectLab = all$symbol,
                 xlim = c(-1.5, 1.5),
                 ylim = c(0,3),
                 xlab = bquote(~Log~'FC'),
                 ylab = bquote(~-Log[10]~(q-value)),
                 subtitle = 'PEA-identified DEPs between AMS4k and AMS1k',
                 title = "",
                subtitleLabSize = 14,
                 titleLabSize = 0,
                 pCutoff =  0.05,
                 FCcutoff = 0,
                 cutoffLineType = 'longdash',
                 cutoffLineCol = 'black',
                 cutoffLineWidth = 0.2,
                 pointSize = 4,
                 labSize = 4,
                 gridlines.major = TRUE,
                 drawConnectors = TRUE,
                 widthConnectors = 0.5,
                 colConnectors = 'grey10',
                 lengthConnectors = unit(0.01, 'npc'),
                 colAlpha = 0.9,
                 labCol = 'black',
                 labFace = 'bold',
                 labhjust = 0.8,
                 labvjust = 0.5,
                 legendLabels = c('','','Nonsignificance',"DEPs"),
                 legendPosition = 'top',
                 legendLabSize = 12,
                 legendIconSize = 4,
                shadeBins = 1,
                caption = ""
)


ggsave( "../results/PEA_volcanov.tiff", height = 16, width = 15, units = "cm", dpi = 300)
save.image("../results/PEA_volcano.RData")

####     Heatmap of 47 PEA-identified DEPs ############

library("heatmaply")
library("RColorBrewer")

group <- read.csv("../data/10_olink_group.csv", sep=",")
# load("../data/PEA_DEPs_rawdata_result_withNA.RData")

heat_raw <- DEPs_raw2[,c("symbol",group$sample)]
heat_raw <- heat_raw %>% tibble::column_to_rownames(var = "symbol")

dir.create("../results/heatmap")
heatmaply(heat_raw,
          file ="../results/47DEPs_heatmap.html",
          color = c("#285BBC","lightyellow","#AF0000"),
          dendrogram = "both",
          hclust_row  = "ward.D2",
          hclust_col  = "ward.D2",
          dist_method = "canberra",
          na.value = "grey50",
          width = 12000,
          height = 12000,
          branches_lwd = 0.3,
          scale = "row",
          fontsize_row = 12,
          fontsize_col = 12,
          key.title = "log2FC",
          k_col=2, k_row =2,
          colorbar_xanchor = "top"  
          
)	

browseURL("../results/heatmap/heatmaply_plot.png")

########    PCA plot of 42 DEPs identified by PEA which with no NA value ############

library(dplyr)

# load("../data/PEA_DEPs_rawdata_result_withNA.RData")
group <- read.csv("../data/10_olink_group.csv",sep = ",",header=T, stringsAsFactors = F, check.names = F)

g <- group[,c("sample","group","LLS")]
pca_raw <- DEPs_raw2[,c("symbol", group$sample )]
pca_raw <- pca_raw %>%
  tibble::column_to_rownames(var = "symbol")

pca_traw <- t(pca_raw)  
pca_traw <- data.frame(pca_traw, stringsAsFactors = F)
raw <- pca_traw %>% mutate(sample = rownames(pca_traw)) %>%
  left_join(g, by = "sample", keep = F) 
sym_del <- c('CA1', 'FCN2', 'TCN2', 'ANG', 'MEGF9')
sym1 <- rownames(pca_raw) 
sym2 <- setdiff(sym1, sym_del)

#is.numeric(pca_raw$FGF23)

ord <- prcomp(raw[,sym2])
summary(ord)
dt <- ord$x
df <- data.frame(dt,raw$group,raw$LLS)
rownames(df) <- pca_raw$sample
head(df)
summ <- summary(ord)
xlab <- paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab <- paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

library(ggplot2)

p1 <- ggplot(df,aes(PC1,PC2,color=raw$group))+
  stat_ellipse(aes(fill=raw.group),type="norm",
               geom="polygon",alpha=0.2,color= NA)+
  guides(fill=F)+
  geom_point(size = 2.5)+
  labs(x=xlab,y=ylab,color="") +
  theme_bw()+theme(panel.grid = element_blank()) +
  theme(axis.text= element_text(size=6,color="black", angle = 0,
                                vjust=0.5, hjust=1)) +
  theme(axis.title= element_text(size=16,color="black", 
                                 vjust=0.5, hjust=0.5)) +
  theme(legend.text = element_text(size=14,color="black", 
                                   vjust=0.5, hjust=0.5))+
  theme(legend.title = element_text(size=16,color="black", 
                                    vjust=0.5, hjust=0.5)) +
  # theme(legend.justification=c(0.01,0.999),
  #       legend.position=c(0.01,0.999)) +
  theme(legend.background = element_rect())
p1

ggsave(p1, file = "../results/PEA_PCA2.tiff", height = 8, width = 11, units = "cm", dpi=300 )


########### KEGG and GO enrichment for two clusters verified in the heatmap ############

library(AnnotationHub)	
library(org.Hs.eg.db) 
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(pathview)

# setwd("/Users/diudiu/Desktop/Figure_v4_20210827/03_olink_enri/src")
# load("../data/PEA_DEPs_rawdata_result_withNA.RData")

symbol <-  DEPs_raw2[c("symbol")]

g2 <- DEPs_raw2 %>% filter(log2FC<0) %>% subset(select = symbol)
g1 <- DEPs_raw2 %>% filter(log2FC>0) %>% subset(select = symbol)


# GO & KEGG enrichment
slist <- list( c1 = g1,
               c2 = g2
)

kegglist <- list()
golist <- list()
tifflist <- list()
for (go_var in names(slist)){
  #  go_var <- names(slist)[[2]] 
  go_i <- slist[[go_var]]
  print (go_var)
  
  g <- apply(go_i,2,as.character) 
  
  id2 <- bitr(g,fromType = "SYMBOL", toType = "ENTREZID" ,OrgDb = "org.Hs.eg.db" )
  id <- subset(id2, select = ENTREZID)
  id <- apply(id,2,as.character)
  
  ego <- enrichGO(gene= id, OrgDb="org.Hs.eg.db",  keyType = 'ENTREZID', ont = 'BP', 
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", readable= TRUE) 
  golist[[go_var]] = paste( "../results/olink_enrichGOBP-2C", go_var,"result.csv",sep="")
  write.csv(ego,file = golist[[go_var]])
  #simplify remove redundant Go
  egosimp <- simplify(ego,cutoff=0.3,by="p.adjust", select_fun = min,measure="Wang")
  
  golist[[go_var]] = paste("../results/olink_enrichGOBP_2C_", go_var,"_simlify0.3_Wang_padjust.csv",sep="")
  write.csv(egosimp,golist[[go_var]])
  
  p1 <- dotplot(egosimp,showCategory=40,color = "qvalue", 
                font.size = 12) +
    theme(axis.text.y= element_text(size=12,color="black", angle = 0,
                                    vjust=1, hjust=1)) 
  p1
  
  tifflist[[go_var]] = paste("../results/PEA_enrichGOBP_", go_var,"_simlify0.3_Wang_padjust.pdf",sep="")
  ggsave(p1,file =tifflist[[go_var]],height = 10, width = 20, units = "cm", dpi = 300)
  
  ekk <- enrichKEGG(gene= id,organism  = 'hsa', 
                    pvalueCutoff = 0.05,pAdjustMethod = "none")	
  
  p2 <- dotplot(ekk,font.size=14, color = "pvalue")	+
    theme(axis.text.y= element_text(size=14,color="black", angle = 0,
                                    vjust=0.5, hjust=0.5)) 
  
  p2
  kegglist[[go_var]] = paste("../results/olink_heatmap_enrichKEGG_dotplot_", go_var,".csv",sep="")
  tifflist[[go_var]] = paste("../results/olink_heatmap_enrichKEGG_dotplot_", go_var,".pdf",sep="")
  write.csv(ekk, file=kegglist[[go_var]])
  ggsave(p2,file =tifflist[[go_var]], 
         height = 9, width = 20, units = "cm", dpi = 300)
  
  genelist <- DEPs_raw2 %>% select(c(symbol,log2FC)) %>%
    arrange(desc(log2FC) ) %>% rename( SYMBOL = 'symbol') %>%
    left_join(id2, by = "SYMBOL") 
  
  FCgenelist <- as.numeric(as.factor(genelist$log2FC)) #numeric vector
  names(FCgenelist) <- as.character(genelist$ENTREZID) #named vector

  print(go_var)
  kk <- ekk %>% filter(pvalue < 0.05)
  keglist <- list()
  for (i in 1:length(kk$ID)){
    ekk_id <- kk$ID[i]
    print(ekk_id)
    
    ekk_id <- pathview(gene.data  = FCgenelist,
                       kegg.dir = "../results/pathview/",
                       pathway.id =  ekk_id,
                       species    = "hsa",
                       limit      = list(gene=max(abs(FCgenelist)), cpd=1)) ## cpd, compound
    # Info: Downloading xml files for hsa04110, 1/1 pathways..
    # Info: Downloading png files for hsa04110, 1/1 pathways..
    # 'select()' returned 1:1 mapping between keys and columns
    # Info: Working in directory /YOUR PATH/Project/clusterProfiler
    # Info: Writing image file hsa04110.pathview.png
  }
}
