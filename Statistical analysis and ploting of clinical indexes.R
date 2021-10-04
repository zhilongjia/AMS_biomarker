library(caret)
library(dplyr)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("paste", "base")
conflict_prefer("as.matrix", "base")
conflict_prefer("rename", "dplyr")

clinic_info=read.csv("../data/MRM_clinical_indexes_AL_YC_53_v2.csv",header = T,stringsAsFactors = F,sep=",",check.names = F)
clinic_group=read.csv("../data/MRM_sampleID_group.csv",header = T,sep=",", stringsAsFactors = F)
clinic_info=merge(clinic_group,clinic_info,by="sampleID")

# remove all undetected clinical indexes
clinic_nna <- clinic_info %>% tibble::column_to_rownames(var = 'sampleID') %>%
  select(-c(group)) %>% filter_all(any_vars(!is.na(.)))
del_sample <- setdiff(clinic_info$sampleID, rownames(clinic_nna))
print(paste("all NA sample:",del_sample,sep = ""))

############## data process for each group #####
group_list <- list(A4 = "AMS4k",
                   A1 ="AMS1k",
                   nA4 = 'nAMS4k',
                   nA1 = 'nAMS1k'
)
for (i in 1:length(group_list)) {
  #i=2
  #select rawdata for each group 
  sub_raw <- clinic_nna %>% tibble::rownames_to_column( var ='sampleID') %>% 
    left_join(clinic_group, by = 'sampleID') %>% filter(group %in% group_list[[i]]) %>%
    tibble::column_to_rownames(var='sampleID') %>% select(-c(group))
  
  # remove 75% undetected index
  sub_traw <- t(sub_raw)
  sub_traw <- data.frame(sub_traw, stringsAsFactors = F, check.names = F) 
  # either NA or 0 represent not detected
  sub_traw[is.na(sub_traw)] <- 0.00
  sub_traw2 <- sub_traw %>% mutate(n = rowSums((sub_traw) != 0),
                                   identified_ratio = n/ncol(sub_traw))%>%
    filter(identified_ratio > 0.75) %>% select(-c(n, identified_ratio))
  omi_sample <- setdiff(colnames(sub_traw),colnames(sub_traw2))
  print(paste(group_list[[i]],"_delete sample:",omi_sample, sep = ""))
  #replace 0 with NA f
  sub_traw2[sub_traw2 == 0] <- NA
  
  final_raw <- t(sub_traw2)
  final_raw <- data.frame(final_raw,stringsAsFactors = F,check.names = F)
  
  repl_raw <- data.frame(apply(final_raw,2,function(x){
    x[is.na(x)]= min(x, na.rm = T) *0.1
    return(x)
  }),stringsAsFactors = F ,check.names = F ) 
  #test 
  # mean(sub_traw1$SCGB3A1_LLLSSLGIPVNHLIEGSQK_y8)
  save(repl_raw, file =  paste("../results/",group_list[[i]],"_clinic_rawdata.RData", sep =""))
}  

######### clinic test ###########
comp2groups <- list(AMS4k_AMS1k=c("AMS1k", "AMS4k"),
                    nAMS4k_nAMS1k=c("nAMS1k", "nAMS4k"),
                    AMS1k_nAMS1k=c("nAMS1k", "AMS1k"),
                    AMS4k_nAMS4k=c("nAMS4k", "AMS4k")
)

for (comp_var in names(comp2groups) ) {
  #comp_var <- names(comp2groups)[[1]]
  comp_i <- comp2groups[[comp_var]]
  print (comp_i)
  ####### data prepared #######
  load(paste("../results/",comp_i[1],"_clinic_rawdata.RData", sep =""))
  test1 <- repl_raw %>% tibble::rownames_to_column(var = 'sampleID')
  load(paste("../results/",comp_i[2],"_clinic_rawdata.RData", sep =""))
  test2 <- repl_raw %>% tibble::rownames_to_column(var = 'sampleID')
  #############  remove index that only exists in one group #######
  n1 <- colnames(test1)
  n2 <- colnames(test2)
  if (length(n1) > length(n2)){
    omi <- setdiff(n1, n2)
    test1 <- test1 %>% select(-c(omi))
    n = comp_i[1]
  }else{
    omi <- setdiff(n2, n1)
    test2 <- test2 %>% select(-c(omi))
    n = comp_i[2]
  }
  # show which group has the omi index
  print(paste(comp_var,"_group_",n, "_delet_clinicl:",omi ,sep=""))
  # merge the 2 group 
  if(comp_var == 'AMS4k_AMS1k'|comp_var== 'nAMS4k_nAMS1k'){
  test_raw <- test1 %>% add_row(test2) %>% filter(sampleID != 'AL1904P73') %>% #filter 成对删除未检出样本及配对样本信息
    tibble::column_to_rownames(var = "sampleID")
  save(test_raw, file = paste("../results/",comp_var,"_clinic_rawdata_replaceNA_del2sample.RData", sep =""))
  
  # log(rawdata)
  test_raw_log <- log2(test_raw) 
  save(test_raw_log, file = paste("../results/",comp_var,"_clinic_rawdata_log_replaceNA_del2sample.RData", sep =""))
  }else{
    test_raw <- test1 %>% add_row(test2) %>% 
      tibble::column_to_rownames(var = "sampleID")
    save(test_raw, file = paste("../results/",comp_var,"_clinic_rawdata_replaceNA_del1sample.RData", sep =""))

    test_raw_log <- log2(test_raw) 
    save(test_raw_log, file = paste("../results/",comp_var,"_clinic_rawdata_log_replaceNA_del1sample.RData", sep =""))
  }
  # number of index
  print(ncol(test_raw))
  
  # group info 
  subsample_pheno <- dplyr::filter(clinic_group, group %in% comp_i )
  g <- cbind(subsample_pheno$group, subsample_pheno$sampleID)
  g <- as.data.frame(g)
  colnames(g) <- c("group", "sampleID")
  # find omi sampleID and remove it 
  group <- g %>% filter(sampleID %in% rownames(test_raw))
 # View(group)                       
  
  #######test_raw for mean and logFC ########
  cols <- split(group$sampleID, group$group)
  
  test_traw <- t(test_raw)
  test_traw <- data.frame(test_traw, stringsAsFactors = F, check.names = F)
 
  print(paste("num clin:",nrow(test_traw),"; num sample:", ncol(test_traw),sep = ""))
  
  sel_mean<- data.frame(lapply(cols, function(i) rowMeans(test_traw[i])))
  sel_FC <- dplyr::mutate(sel_mean, FC= sel_mean[,comp_i[2]]/sel_mean[,comp_i[1]], log2FC=log2(FC)) %>%
    tibble::rownames_to_column(var = 'clinic')
  
  ### raw for clinc test
  x1 <- as.matrix(unlist(cols[1]))[,1]
  x2 <- as.matrix(unlist(cols[2]))[,1]
  
  temp=NULL
  for (j in 1:nrow(test_traw)){
    # j=1
    test_1=apply(test_traw[j,x1],2, as.numeric)
    test_2=apply(test_traw[j,x2],2, as.numeric)
    if(length(unique(test_1))==1|length(unique(test_2))==1)
    {
      group_1=FALSE
      group_2=FALSE
    }else
    {
      group_1=shapiro.test(test_1)$p.value>0.05
      group_2=shapiro.test(test_2)$p.value>0.05
    }
    
    if(group_1&group_2)
    {
      pair_c=length(test_1)==length(test_2)
      var=var.test(test_1,test_2)$p.value>0.05
      test_palue=t.test(test_1,test_2,paired = pair_c,var.equal =var)
      test_mean= mean(test_1)/mean(test_2)
      log2FC = log2(test_mean)
      temp=rbind(temp,c(paste(comp_var),rownames(test_traw)[j],test_palue$p.value,
                        paste("t.test pair=",pair_c,";var equal=",var,sep=""), log2FC))
    }else
    {
      if(length(test_1)==length(test_2))
      {
        test_palue=wilcox.test(test_1,test_2,paired =  TRUE)
        test_mean= mean(test_1)/mean(test_2)
        log2FC = log2(test_mean)
        temp=rbind(temp,c(paste(comp_var),rownames(test_traw)[j],test_palue$p.value,"wilcox.test-paird",log2FC))
      }else
      {
        test_palue=wilcox.test(test_1,test_2)
        test_mean= mean(test_1)/mean(test_2)
        log2FC = log2(test_mean)
        temp=rbind(temp,c(paste(comp_var),rownames(test_traw)[j],test_palue$p.value,"wilcox.test-no-paird",log2FC))
      }
      
    }
  }
  
  colnames(temp)=c("groups","clinic","p_value","Method","log2FC")
  q_value <- p.adjust(temp[,'p_value'],"bonferroni")
  re = cbind(temp, q_value)
  re <- data.frame(re, stringsAsFactors = F,check.names = F)
  write.table(re,file = paste("../results/",comp_var,"_clinic_pq_FC.csv",sep = ""),sep=",",col.names = T,row.names = F,quote = F)
  
  clin_raw_re <- test_traw %>% tibble::rownames_to_column(var ='clinic') %>%
    left_join(re, by = 'clinic') %>% 
  select(c(clinic,colnames(test_traw), log2FC, p_value,q_value))
  
  save(clin_raw_re,file= paste("../results/",comp_var,"_clinic_raw_result.RData",sep = ""))
}


#########   violin plots of differential clinical indexes between each comparision #################
library(dplyr)
library(ggplot2)
library(ggpubr)

library(tidyr)

group=read.csv("../data/mrm_sample_list_v4.csv")
clinic_group <- group %>% select(c(sample_id, group))

my_comprisons=list(c("AMS4k","AMS1k"),c("nAMS4k","nAMS1k"),c("AMS1k","nAMS1k"),
                   c("AMS4k","nAMS4k" ) )

col_list <- list(AMS4k_vs_AMS1k = c('#F3574A',"#01569B"),
                 nAMS4k_vs_nAMS1k = c('#FAB70F',"#5FC6DB"),
                 AMS1k_vs_nAMS1k = c("#F3574A",'#FAB70F'),
                 AMS4k_vs_nAMS4k= c("#01569B","#5FC6DB")
                 
)


for (i in 1:length(my_comprisons) )
{
 i= 3
  print(paste(my_comprisons[[i]],collapse = "_vs_"))
  col <- col_list[[paste(my_comprisons[[i]],collapse = "_vs_")]]
  load(paste("../results/",paste(my_comprisons[[i]],collapse = "_"),"_clinic_raw_result.RData",sep = ""))
  clinc_meth <- read.csv(paste("../results/",paste(my_comprisons[[i]],collapse = "_"),"_clinic_pq_FC.csv",sep = ""),
                        stringsAsFactors = F,check.names = F)
   # clinical indexes in the pathogenesis and protection grups were selected with q-value less than 0.05, while the other two groups were selected with p-value less than 0.05
  if(i==1 |i==2){
    raw <- clin_raw_re %>% filter(q_value < 0.05) %>% select(-c(log2FC, p_value, q_value)) %>%
      tibble::column_to_rownames(var ='clinic')
    index <-  clinc_meth[which(clinc_meth[,6]<0.05),2]
  }else{
    raw <- clin_raw_re %>% filter(p_value < 0.05) %>% select(-c(log2FC, p_value, q_value)) %>%
      tibble::column_to_rownames(var ='clinic')
    index <-  clinc_meth[which(clinc_meth[,3]<0.05),2]
  }
  td <- t(raw)
  td <- data.frame(td, stringsAsFactors = F,check.names = F)
  
  sam <- data.frame(rownames(td))
  colnames(sam) <- "sample_id"
  
  cat <- sam %>% left_join(clinic_group, by ='sample_id')
  
  convert <- function(td) {
    lapply(1:ncol(td), function(i) {
      # i=1
      d2 <- td[, c(i)]
      d3 <- data.frame(d2) %>% mutate(symbol = rownames(raw)[i])
      
      d4 <- cbind(d3,cat)
    }) %>% do.call('rbind', .)
  }
  dd = convert(td)
  print(nrow(dd)/ncol(td))
  data <- subset(dd, select = -sample_id) %>% rename(value ='d2', feature ='symbol')
   
  
  if(i==1 |i==2){
    label <- clinc_meth[which(clinc_meth[,6]<0.05),4]
  }else{
    label <- clinc_meth[which(clinc_meth[,3]<0.05),4]
  }

  
  label[which(label == 'wilcox.test-paird'| label =='wilcox.test-no-paird')] <- 'Wilcox test'
  label[which(label == 't.test pair=TRUE;var equal=FALSE'| 
                label =='t.test pair=TRUE;var equal=TRUE'| 
                label =='t.test pair=FALSE;var equal=TRUE'| 
                label =='t.test pair=FALSE;var equal=FALSE')] <- 't test'
  
  labels <- data.frame(feature = index, labels = label)
  labels <- labels %>% mutate (q = 'q<0.05')
  lab <- tidyr::unite(labels, 'method', labels, q, sep = ' ')

  clinic_data <- data %>% left_join(lab, by = 'feature') 

  clinic_data2 <- tidyr::unite(clinic_data, 'clinic' , feature, method,sep = "\n")
  
  p <- ggplot(clinic_data2,aes(group,value, fill = group))+
    geom_violin(alpha=0.8,width=1,trim = FALSE) +
    geom_jitter(fill = "black",width =0.2,shape = 21,size=0.8) +
    # geom_label(data = labels, aes(label=label),
    #            x = -Inf, y = Inf, hjust=-0.8, vjust=1,
    #            inherit.aes = FALSE) +
    facet_wrap(~clinic, scales = "free" ,ncol=5) +
    theme_bw() +
    theme(strip.background = element_rect(
      color = "white", fill = "white"),
      panel.border = element_rect(
        color = "black")) +theme(panel.grid=element_blank())+
    theme(axis.text= element_text(size=12,color="black", angle = 0,
                                  vjust=0.5, hjust=0.5)) +
    theme(axis.title= element_text(size=12,color="black", 
                                   vjust=0.5, hjust=0.5)) +
    theme(strip.text = element_text(size=12,color="black", 
                                    vjust=0.5, hjust=0.5))+
    theme(legend.text = element_text(size=14,color="black", 
                                     vjust=0.5, hjust=0.5))+
    theme(legend.title = element_text(size=14,color="black", 
                                      vjust=0.5, hjust=0.5)) 
  
  p1 <- change_palette(p, palette = col)
  p1
  if (i ==1 |i==2 ){
    height <- 35
    width <- 35
  }else{
    height <- 7*ceiling(length(index)/5)
    width <- 35
  }
  ggsave(paste(paste(my_comprisons[[i]],collapse = "_vs_"),"_select_clinical_method_naomit.pdf",sep=""), p1,
         height = height, width =width, units = "cm", dpi = 300)

}
save.image("select_proteins_boxplot_method.RData")

############## Heatmap of differential clinical indexes in each comparasion ############
library("pheatmap")
library("RColorBrewer")
library(dplyr)
############ clinic overlap in 4 groups #######
adjplist <- list(AMS4k_AMS1k = "../data/AMS4k_AMS1k_clinic_pq_FC.csv",
                 nAMS4k_nAMS1k = "../data/nAMS4k_nAMS1k_clinic_pq_FC.csv",
                 AMS1k_nAMS1k = "../data/AMS1k_nAMS1k_clinic_pq_FC.csv",
                 AMS4k_nAMS4k = "../data/AMS4k_nAMS4k_clinic_pq_FC.csv"
)

sym <- data.frame(var_a =character())
for (raw_var in names(adjplist)) {
  #raw_var <- names(adjplist)[[3]]
  raw_i <- adjplist[[raw_var]]
  raw <- read.csv(raw_i, sep = ",",stringsAsFactors = F,check.names = F)
 
  data_extract <- function(data) {
    if(raw_var == "AMS4k_AMS1k" | raw_var == "nAMS4k_nAMS1k"){
      raw <- data %>% filter(q_value < 0.05) %>%
        subset(select = c("clinic"))
    }else{
      raw <- data %>% filter(p_value < 0.05) %>%
        subset(select = c("clinic"))
    }
  } 
  sym <- rbind(sym,data_extract(raw))
} 

clin_union <- unique(sym)

write.csv(clin_union , file = "../results/clinic_4groups_overlap.csv")

######## clinic log2FC qvalue extraction #####
raw_list <- list()
q_list <- list()

clin_union  <- read.csv("../results/clinic_4groups_overlap.csv", sep = ",", check.names = F, row.names = 1)

for (raw_var in names(adjplist)) {
  # raw_var <- names(adjplist)[[3]]
  raw_i <- adjplist[[raw_var]]
  raw <- read.csv(raw_i, sep = "," ,stringsAsFactors = F,check.names = F)
  
  heat_raw <- raw %>% right_join(clin_union, by = 'clinic' ) %>%
    subset(select = c("clinic", "log2FC")) %>% arrange( clinic)  %>% # adjusted_p_value
    tibble::column_to_rownames(var = "clinic") 
  colnames(heat_raw) <- raw_var
  raw_list[[raw_var]] <- heat_raw #q_list
  if(raw_var == "AMS4k_AMS1k" | raw_var == "nAMS4k_nAMS1k"){
    q_sub <- raw %>% right_join(clin_union, by = 'clinic' ) %>%
      subset(select = c("clinic", "q_value")) %>% arrange( clinic)  %>% # adjusted_p_value
      tibble::column_to_rownames(var = "clinic") 
  }else{
    q_sub <- raw %>% right_join(clin_union, by = 'clinic' ) %>%
      subset(select = c("clinic", "p_value")) %>% arrange( clinic)  %>% # adjusted_p_value
      tibble::column_to_rownames(var = "clinic") 
  }
  colnames(q_sub) <- raw_var
  q_list[[raw_var]] <- q_sub
  
}

raw <- do.call(cbind,raw_list)
q <- do.call(cbind,q_list)

q[is.na(q)] <- 1
q1 <- q %>% select(c(AMS4k_AMS1k,nAMS4k_nAMS1k))

if (!is.null(q1)){
  ssmt <- q1< 0.05
  q1[ssmt] <-'*'
  q1[!ssmt]<- ''
} else {
  q1 <- F
}
q2 <- q %>% select(c(AMS1k_nAMS1k,AMS4k_nAMS4k))
if (!is.null(q2)){
  ssmt <- q2< 0.05
  q2[ssmt] <-'+'
  q2[!ssmt]<- ''
} else {
  q2 <- F
}


q_final <- cbind(q1,q2)

save_pheatmap_pdf <- function(x, filename, width=4.5, height=8) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

bk <- c(seq(-1.2,-0.1,by=0.1),seq(0,1.2,by=0.1))
p <- pheatmap(raw,
              color = c(colorRampPalette(colors = c("#285BBC","lightyellow"))(length(bk)/2),
                        colorRampPalette(colors = c("lightyellow","#AF0000"))(length(bk)/2)),
              legend_breaks=seq(-1,1,0.5),
              breaks = bk,
              cluster_rows  = FALSE,
              cluster_cols  = FALSE,
              display_numbers = q_final,
              width = 100,
              height = 100,
              fontsize_row = 12,
              fontsize_col = 12,
              legend = TRUE,
              legend.position= "left",
              cellwidth = 22,
              cellheight = 12,
              fontsize_number = 14,
              # show_colnames = F,
              angle_col = 45,
              border=NA,
              na_col = "white"
              
)	
print(p)
save_pheatmap_pdf(p, file = "../results/MRM_4groups_clinical_indexes_heatmap_qp.pdf")

save.image("MRM_4groups_clinc_index_heatmap_qp.RData")
