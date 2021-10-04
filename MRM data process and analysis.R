######   comparasion among all the four groups(AMS1k, AMS4k,nAMS1k, nAMS4k)  ##########

library(dplyr)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("paste", "base")
conflict_prefer("rename", "dplyr")

rawmetadata <- read.csv("../data/mrm_sample_list_v4.csv", sep =",")
rawdata <- readxl::read_xlsx("../data/mrm_df_v3.xlsx")

load("../data/PEA_comb_numeric.RData")

 PEA_DEP <- data %>% tibble::rownames_to_column(var = 'symbol') %>% 
   filter(q_value < 0.05) %>% 
 select( c(symbol,log2FC)) %>% rename(PEA_log2FC = log2FC)

library(matrixTests)


comp2groups <- list(AMS4k_AMS1k=c("AMS1k", "AMS4k"),
                    nAMS4k_nAMS1k=c("nAMS1k", "nAMS4k"),
                    AMS1k_nAMS1k=c("nAMS1k", "AMS1k"),
                    AMS4k_nAMS4k=c("nAMS4k", "AMS4k"),
                    AMS4k_nAMS1k=c("nAMS1k", "AMS1k"),
                    AMS1k_nAMS4k=c("nAMS1k", "AMS1k")
)
DEA_list <- list()
result_comb <- rownames(raw_sym)
for (comp_var in names(comp2groups) ) {
  #comp_var <- names(comp2groups)[[1]]
  comp_i <- comp2groups[[comp_var]]
  print (comp_i)
  ####### data prepared #######
  load(paste("../results/",comp_i[1],"_MRM_symbol_rawdata_replaceNA.RData", sep =""))
  test1 <- raw_sym
  load(paste("../results/",comp_i[2],"_MRM_symbol_rawdata_replaceNA.RData", sep =""))
  test2 <- raw_sym
  test_raw <- test1 %>% left_join(test2, by = 'symbol') %>% 
  
                        tibble::column_to_rownames(var = "symbol")
  #save(test_raw, file = paste("../results/",comp_var,"_MRM_rawdata_replaceNA.RData", sep =""))
  # log(rawdata) for t-test
  test_raw_log <- log2(test_raw) 
  #save(test_raw_log, file = paste("../results/",comp_var,"_MRM_rawdata_log_replaceNA.RData", sep =""))

  # number of symbol validated
  print(nrow(test_raw))

  # group info 
  subsample_pheno <- dplyr::filter(rawmetadata, group %in% comp_i )
  group <- cbind(subsample_pheno$group, subsample_pheno$sampleID)
  group <- as.data.frame(group)
  colnames(group) <- c("group", "sampleID")
  #######test_raw for mean and logFC ########
  cols <- split(group$sampleID, group$group)
  
  sel_mean<- data.frame(lapply(cols, function(i) rowMeans(test_raw[i])))
  sel_FC <- dplyr::mutate(sel_mean, FC= sel_mean[,comp_i[2]]/sel_mean[,comp_i[1]], log2FC=log2(FC)) %>%
     tibble::rownames_to_column(var = 'symbol')
 # View(sel_FC)
  
  x1 <- as.matrix(unlist(cols[1]))[,1]
  x2 <- as.matrix(unlist(cols[2]))[,1]

  ######### test_raw_log for t.test  ########
  temp=NULL
  for (j in 1:nrow(test_raw_log))
  {
    #j=1
    test_1=apply(test_raw_log[j,x1],2, as.numeric)
    test_2=apply(test_raw_log[j,x2],2, as.numeric)
   #t.test for log rawdata, check var.test and paired 
    pair_c=length(test_1)==length(test_2)
    var=var.test(as.matrix(test_raw_log[j,x1]),as.matrix(test_raw_log[j,x2]))$p.value>0.05
    test_palue=t.test(test_1,test_2,paired = pair_c,var.equal =var)
    temp=rbind(temp,c(paste(comp_var),rownames(test_raw_log)[j],test_palue$p.value,
                        paste("t.test pair=",pair_c,";var equal=",var,sep="")))
  }
  temp2 <- data.frame(temp, stringsAsFactors = F)    
  colnames(temp2) <- c("groups",'symbol','p_value','method')
  #merge result and adjust pvalue
  test_result <- temp2 %>% left_join(sel_FC, by = "symbol") %>% 
      mutate(q_value = p.adjust(p_value, method = 'BH'))
  write.csv(test_result, file = paste("../results/",comp_var,"_MRM_test_result.csv", sep =""))
 # save(test_result, file = paste("../results/",comp_var,"_MRM_test_result.RData", sep =""))
  # rawdata with result 
  MRM_raw_re <- test_raw %>% tibble::rownames_to_column(var ='symbol') %>%
    left_join(test_result, by ='symbol') %>% select(-c(groups,method,p_value,FC))
 # save(MRM_raw_re, file = paste("../results/",comp_var,"_MRM_raw_result.RData", sep =""))
  #log raw with result
  MRM_log_re <- test_raw_log %>% tibble::rownames_to_column(var ='symbol') %>%
    left_join(test_result, by ='symbol') %>% select(-c(groups,method,p_value,FC))
 # save(MRM_log_re, file = paste("../results/",comp_var,"_MRM_lograw_result.RData", sep =""))
  
  result_comb <- cbind(result_comb, test_result, stringsAsFactors =F)
  
  dep_MRM <- test_result %>% filter(q_value < 0.05) %>% select(c(symbol,log2FC))
  #save(dep_MRM, file = paste("../results/",comp_var,"_MRM_DEPs.RData", sep =""))
  
  #num of DEPs
  print(nrow(dep_MRM))
 
  #########################   match with PEA      ########################

   mrm_match_PEA <- test_result %>% left_join(PEA_DEP, by = "symbol") %>%
    mutate(sig = sign(log2FC*PEA_log2FC) ) 
  
  n <- mrm_match_PEA %>% group_by(sig) %>% summarise(sum(sig))
  print(n)

}
  write.csv(result_comb, file ="../results/_MRM_allgroup_test_result_combined.csv")
save.image("MRM_allgroup_test.RData")

######## FUnctional enrichment of DEPs involved in AMS and nAMS ##############

library(AnnotationHub)	
library(org.Hs.eg.db)  
library(clusterProfiler)
library(dplyr)
library(ggplot2)

load("../data/AMS4k_AMS1k_MRM_test_result.RData")
sym <- test_result %>% filter(q_value < 0.05) %>% select(symbol)
#sym <- subsample_result
g <- apply(sym,2,as.character) 

id <- bitr(g,fromType = "SYMBOL", toType = "ENTREZID" ,OrgDb = "org.Hs.eg.db" )
id <- subset(id, select = ENTREZID)
id <- apply(id,2,as.character)

#KEGG#
ekk <- enrichKEGG(gene= id,organism  = 'hsa',
                  pvalueCutoff = 0.05,pAdjustMethod = "none")	
ekk_symbol <- setReadable(ekk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")


p1 <- dotplot(ekk_symbol,font.size=10, color = "pvalue") +	
 theme(axis.text= element_text(size=14,color="black", angle = 0,
                                  vjust=0.5, hjust=0.5)) +
  theme(axis.title= element_text(size=14,color="black", 
                                 vjust=0.5, hjust=0.5)) +
  theme(legend.text = element_text(size=14,color="black", 
                                   vjust=0.5, hjust=0.5))+
  theme(legend.title = element_text(size=14,color="black", 
                                    vjust=0.5, hjust=0.5)) 
p1
ggsave(p1,file = "../results/MRM_AMS4k_AMS1k_DEPs_KEGG_dotplot.tiff",
       width = 16.5, height = 9.5, units = "cm", dpi =300)

## Enrichemnt of DEPs in nAMS

g <- c('ADAM15', 'CD38', 'CST6','KITLG', 'THBD')

id <- bitr(g,fromType = "SYMBOL", toType = "ENTREZID" ,OrgDb = "org.Hs.eg.db" )
id <- subset(id, select = ENTREZID)
id <- apply(id,2,as.character)
  
ego <- enrichGO(gene= id, OrgDb="org.Hs.eg.db",  keyType = 'ENTREZID', ont = "BP", 
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", readable= TRUE)

egosimp <- simplify(ego,cutoff=0.3,by="p.adjust", select_fun = min,measure="Wang")
p <- dotplot(egosimp,showCategory=40,color = "p.adjust",font.size = 8)+  
  theme(axis.text= element_text(size=20,color="black", angle = 0,
                                vjust=0.5, hjust=0.5)) +
  theme(axis.title= element_text(size=14,color="black", 
                                 vjust=0.5, hjust=0.5)) +
  theme(legend.text = element_text(size=12,color="black", 
                                   vjust=0.5, hjust=0.5))+
  theme(legend.title = element_text(size=12,color="black", 
                                    vjust=0.5, hjust=0.5)) 
print(p)
write.csv(egosimp, file = "../results/nAMS_5P_GenegoPathway.csv")
ggsave(p,file = "../results/nAMS_5P_GenegoPathway.tiff",
         width = 16, height = 13, units = "cm")
         
ekk <- enrichKEGG(gene= id,organism  = 'hsa', 
                  pvalueCutoff = 0.05,pAdjustMethod = "none")	 

p1 <- dotplot(ekk,font.size=10, color = "pvalue") +	
  theme(axis.text= element_text(size=20,color="black", angle = 0,
                                vjust=0.5, hjust=0.5)) +
  theme(axis.title= element_text(size=12,color="black",
                                 vjust=0.5, hjust=0.5)) +
  theme(legend.text = element_text(size=12,color="black", 
                                   vjust=0.5, hjust=0.5))+
  theme(legend.title = element_text(size=12,color="black", 
                                    vjust=0.5, hjust=0.5)) 
p1
write.csv(ekk, file="../results/mrm_nAMS_5symbol_KEGG_enri.csv")
ggsave(p1,file = "../results/nAMS_5P_KEGG.tiff",
       width = 18, height = 10, units = "cm")


##########    Violin plot of proteins involved in carbohydrate metabolism #########
getwd()
library(tidyverse)
library(rstatix)
library(ggpubr)
library(dplyr)
library(ggplot2)

load("../data/mrm_sample_list_v4.RData")

my_comprisons=list(c("AMS4k","AMS1k"),c("nAMS4k","nAMS1k"),c("AMS1k","nAMS1k"),
                   c("AMS4k","nAMS4k" ) )

sym <- c( "ALDOA","PCK1","RBKS",
          "PDHA1",'PFKM','PGM1',"G6PC","LDHA") 

######## ggplot data process 4group####
data_comb <- data.frame(d2=1,symbol=1,group=1)
for (i in 1:length(my_comprisons)) {
 #i=3

  filename <- paste("../data/",paste(my_comprisons[[i]],collapse ="_"),"_MRM_lograw_result.RData", sep="")
 load(filename)
car_raw <- MRM_log_re %>% select(-c(log2FC,q_value,my_comprisons[[i]][1],my_comprisons[[i]][2])) %>%
  filter(symbol %in% sym) %>% tibble::column_to_rownames(var= 'symbol')
td <- t(car_raw)
td <- data.frame(td, stringsAsFactors = F,check.names = F)
sam <- data.frame(rownames(td))
colnames(sam) <- "sampleID"

cat <- sam %>% mutate(group = sapply(strsplit(sampleID, "\\_"), `[`, 2 ))

convert <- function(td) {
  lapply(1:ncol(td), function(i) {
   # i=1
    d2 <- td[, c(i)]
    d3 <- data.frame(d2) %>% mutate(symbol = rownames(car_raw)[i])
    
    d4 <- cbind(d3,cat)
  }) %>% do.call('rbind', .)
}
dd = convert(td)
print(nrow(dd)/ncol(td))

dd2 <- subset(dd, select = -sampleID)
data_comb <- rbind(data_comb,dd2)
} 

data <- data_comb %>% filter(symbol %in% sym) %>%
          subset(select = d2) %>% mutate_if(is.character, as.numeric) %>%
  rename(log2Intensity = d2)
####check ??
data <- cbind(data, data_comb[-1,-1])
data <- data.frame(data, stringsAsFactors = F, check.names = F)
########## ggplot geom_violin #######

col <- c('#F3574A',"#01569B",'#FAB70F',"#5FC6DB")

p <-  ggplot(data,aes(group,log2Intensity, fill = group))+
    geom_violin(alpha=0.8,width=1,trim = FALSE) +
      geom_jitter(fill = "black",width =0.2,shape = 21,size=0.8) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 1.5, color = "black",fill='white') +
  facet_wrap(~symbol, scales = "free" , ncol =2) +
    theme_bw() +
    theme(strip.background = element_rect(
      color = "white", fill = "white"),
      panel.border = element_rect(
        color = "black")) +theme(panel.grid=element_blank())+
    theme(axis.text= element_text(size=16,color="black", angle = 0,face='bold',
                                  vjust=0.5, hjust=0.5)) +
    theme(axis.title= element_text(size=16,color="black", face='bold',
                                   vjust=0.5, hjust=0.5)) +
    theme(strip.text = element_text(size=16,color="black", face='bold',
                                    vjust=0.5, hjust=0.5))+
    theme(legend.text = element_text(size=16,color="black", 
                                     vjust=0.5, hjust=0.5))+
    theme(legend.title = element_text(size=16,color="black",
                                      vjust=0.5, hjust=0.5)) 
p
p1 <- change_palette(p, palette = col)
p1

ggsave(p1, file = "../results/MRM_4group_metabolsim_violin_v3.tiff",
       width = 28, height =40, units = "cm", dpi = 300)

save.image(file = "mrm_4group_violinplot.RData") 

########  Violin plot of 5 key proteins ##########

library(tidyverse)
library(rstatix)
library(ggpubr)
library(dplyr)

load("../data/mrm_sample_list_v4.RData")

my_comprisons=list(c("AMS4k","AMS1k"),c("nAMS4k","nAMS1k"),c("AMS1k","nAMS1k"),
                   c("AMS4k","nAMS4k" ) )

sym <- c( "RET",'ADAM15',"PHGDH","TRAF2","S100A12") 

######## ggplot data process 4group ####
data_comb <- data.frame(d2=1,symbol=1,group=1)
for (i in 1:length(my_comprisons)) {
  #i=3
  
  filename <- paste("../data/",paste(my_comprisons[[i]],collapse ="_"),"_MRM_lograw_result.RData", sep="")
  load(filename)
  car_raw <- MRM_log_re %>% select(-c(log2FC,q_value,my_comprisons[[i]][1],my_comprisons[[i]][2])) %>%
    filter(symbol %in% sym) %>% tibble::column_to_rownames(var= 'symbol')
  td <- t(car_raw)
  td <- data.frame(td, stringsAsFactors = F,check.names = F)
  sam <- data.frame(rownames(td))
  colnames(sam) <- "sampleID"
  
  cat <- sam %>% mutate(group = sapply(strsplit(sampleID, "\\_"), `[`, 2 ))
  
  convert <- function(td) {
    lapply(1:ncol(td), function(i) {
      # i=1
      d2 <- td[, c(i)]
      d3 <- data.frame(d2) %>% mutate(symbol = rownames(car_raw)[i])
      
      d4 <- cbind(d3,cat)
    }) %>% do.call('rbind', .)
  }
  dd = convert(td)
  print(nrow(dd)/ncol(td))
  
  dd2 <- subset(dd, select = -sampleID)
  data_comb <- rbind(data_comb,dd2)
} 

data <- data_comb %>% filter(symbol %in% sym) %>%
  subset(select = d2) %>% mutate_if(is.character, as.numeric) %>%
  rename(log2Intensity = d2)
####check ??
data <- cbind(data, data_comb[-1,-1])
data <- data.frame(data, stringsAsFactors = F, check.names = F)
########## ggplot geom_violin #######

col <- c('#F3574A',"#01569B",'#FAB70F',"#5FC6DB")

p <-  ggplot(data,aes(group,log2Intensity, fill = group))+
  geom_violin(alpha=0.8,width=1,trim = FALSE) +
  geom_jitter(fill = "black",width =0.2,shape = 21,size=0.8) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 1.5, color = "black",fill='white') +
  facet_wrap(~symbol, scales = "free" , ncol =1) +
  theme_bw() +
  theme(strip.background = element_rect(
    color = "white", fill = "white"),
    panel.border = element_rect(
      color = "black")) +theme(panel.grid=element_blank())+
  theme(axis.text= element_text(size=16,color="black", angle = 0,face='bold',
                                  vjust=0.5, hjust=0.5)) +
  theme(axis.title= element_text(size=16,color="black", face='bold',
                                 vjust=0.5, hjust=0.5)) +
  theme(strip.text = element_text(size=16,color="black", face='bold',
                                  vjust=0.5, hjust=0.5))+
  theme(legend.text = element_text(size=14,color="black", 
                                   vjust=0.5, hjust=0.5))+
  theme(legend.title = element_text(size=14,color="black",
                                    vjust=0.5, hjust=0.5)) 
p
p1 <- change_palette(p, palette = col)
p1

ggsave(p1, file = "../results/MRM_5_key_biomarker_violin.tiff",
       width = 18, height = 45, units = "cm", dpi = 300)

save.image(file = "mrm_4group_key_biomarker.RData") 

###################   Heatmap of DEPs identified by MRM & PEA #############

library("pheatmap")
library("RColorBrewer")
library(dplyr)


load("../data/mrm_sample_listv2.RData")
load('../data/PEA_rawdata_result_combine.RData')

PEA_DEP <- read.csv("../data/PEA_DEPs_results_t.csv", sep=",", stringsAsFactors = F)

adjplist <- list(AMS4k_AMS1k = "../data/AMS4k_AMS1k_MRM_test_result.RData",
                 nAMS4k_nAMS1k = "../data/nAMS4k_nAMS1k_MRM_test_result.RData",
                 AMS1k_nAMS1k = "../data/AMS1k_nAMS1k_MRM_test_result.RData",
                 AMS4k_nAMS4k = "../data/AMS4k_nAMS4k_MRM_test_result.RData"
)

########## 4 compgroup DEPs overlap ######
sym <- data.frame(var_a =character())
for (raw_var in names(adjplist)) {
  #raw_var <- names(adjplist)[[3]]
  raw_i <- adjplist[[raw_var]]
  load(raw_i)
  data_extract <- function(data) {
    raw <- data %>% filter(q_value < 0.05) %>%
      subset(select = c("symbol"))
  } 
  sym <- rbind(sym,data_extract(test_result))
} 

mrm_sym <- unique(sym)
write.csv(mrm_sym, file = "../results/MRM_DEPs_4group_overlap.csv")
#merge PEA DEPs
PEA_sym <- PEA_DEP %>% subset(select = c("symbol")) #, "log2FC", "adjust_p_value"
symbols <- unique(rbind(sym, PEA_sym))

PEA_raw <- PEA_comb %>% filter(symbol %in% symbols$symbol) %>%
  select( c(symbol,log2FC))
colnames(PEA_raw)[2] <- "PEA_AMS4k_AMS1k"
PEA_q <- PEA_comb %>% filter(symbol %in% symbols$symbol) %>%
  select( c(symbol,q_value)) 
colnames(PEA_q)[2] <- "PEA_AMS4k_AMS1k"

# MRM logFC qvalue 
raw_list <- list()
q_list <- list()

for (raw_var in names(adjplist)) {
  #raw_var <- names(adjplist)[[2]]
  raw_i <- adjplist[[raw_var]]
  load(raw_i)
  
  heat_raw <- test_result %>% right_join(symbols, by = 'symbol' ) %>%
    subset(select = c("symbol", "log2FC")) %>% arrange( symbol)  %>% 
    tibble::column_to_rownames(var = "symbol") 
  colnames(heat_raw) <- raw_var
  raw_list[[raw_var]] <- heat_raw #q_list
  q_sub <- test_result %>% right_join(symbols, by = 'symbol' ) %>%
    subset(select = c("symbol", "q_value")) %>% arrange( symbol)  %>% 
    tibble::column_to_rownames(var = "symbol") 
  colnames(q_sub) <- raw_var
  q_list[[raw_var]] <- q_sub
  
}

mrm_raw <- do.call(cbind,raw_list)
q_raw <- do.call(cbind,q_list)
# merge PEA *MRM log2FC qvalue 
raw1 <- mrm_raw %>% tibble::rownames_to_column(var = 'symbol') %>%
  left_join(PEA_raw, by = 'symbol', keep = F) %>% 
  tibble::column_to_rownames(var = "symbol") %>% select(5,1:4)
raw <- apply(raw1,2, as.numeric)
rownames(raw) <- rownames(raw1)
q1 <- q_raw %>% tibble::rownames_to_column(var = 'symbol') %>%
  left_join(PEA_q, by = 'symbol', keep = F) %>% 
  tibble::column_to_rownames(var = "symbol") %>% select(5,1:4)
q <- apply(q1,2, as.numeric)
rownames(q) <- rownames(q1)

# NA replace with 1 to represent nonsig
q[is.na(q)] <- 1

# sig replace with * 
if (!is.null(q)){
  ssmt <- q< 0.05
  q[ssmt] <-'*'
  q[!ssmt]<- ''
} else {
  q <- F
}

######## pheatmap plot #####
bk <- c(seq(-2,-0.1,by=0.1),seq(0,2,by=0.1))

# save image function
save_pheatmap_pdf <- function(x, filename, width=6, height=13) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
p <- pheatmap(raw,
              color = c(colorRampPalette(colors = c("#285BBC","lightyellow"))(length(bk)/2),
                        colorRampPalette(colors = c("lightyellow","#AF0000"))(length(bk)/2)),
              legend_breaks=seq(-1.5,1.5,1),
              breaks = bk,
              cluster_rows  = FALSE,
              cluster_cols  = FALSE,
              display_numbers = q,
              width = 100,
              height = 100,
              fontsize_row = 10,
              fontsize_col = 10,
              legend = TRUE,
              legend.position= "left",
              cellwidth = 22,
              cellheight = 10,
              fontsize_number = 12,
              # show_colnames = F,
              angle_col = 45,
              border=NA,
              na_col = "white"
              #  display_numbers = matrix(ifelse(heat_raw != 0, "*", ""))
)	
print(p)
save_pheatmap_pdf(p, filename = "../results/MRM_overlap_PEA_DEPs_heatmap.pdf")

save.image("PEA_mrm_overlap_proteins_heatmap.RData")

