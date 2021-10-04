library(dplyr)
library(psych)
library(ComplexHeatmap)

rawmetadata <- read.csv('../data/mrm_sample_list_v4.csv',sep = ',')
LLS=read.csv("../data/LLQ_score_1k_4k.csv",header = T,sep=",",check.names = F)
               
comp_i <- c("nAMS1k", "AMS1k")

####### data prepared #######
#MRM all symbol
load(paste("../data/",comp_var,"_MRM_rawdata_replaceNA.RData", sep =""))
load(paste("../data/",comp_var,"_MRM_raw_result.RData", sep =""))
symbol <- rownames(test_raw)
sub_symbol <- MRM_raw_re %>% filter(q_value <0.05) 
sub_symbol <- sub_symbol$symbol
load(paste("../data/",comp_var,"_clinic_raw_result.RData", sep =""))
clinic <- clin_raw_re$clinic
sub_clinic <- clin_raw_re %>% filter(p_value <0.05)
sub_clinic <- sub_clinic$clinic

raw <- read.csv(paste("../data/",comp_var,"_MRM_clinic_raw.csv", sep =""),check.names = F,stringsAsFactors = F)

colnames(raw)[1] <- 'sample_id'
raw[raw==0] <- 0.001
corr_raw <- raw %>% left_join(LLS, by = 'sample_id') 

#remove phen 75% unrecorded 
phen_data <- corr_raw %>%
  tibble::column_to_rownames(var='sample_id')  %>% select(-(colnames(raw)[-1])) 

t_phen <- t(phen_data)
t_phen <- data.frame(t_phen, stringsAsFactors = F,check.names = F)

phen_raw <- t_phen %>% mutate(n = rowSums(t_phen != 0),
                              identified_ratio = n/ncol(t_phen)) %>%
  filter(identified_ratio > 0.25) %>%
  select(- c("n", "identified_ratio")) 

t_phen_raw <- t( phen_raw)
t_phen_raw <- data.frame(t_phen_raw, stringsAsFactors = F,check.names = F)

#nAMS本身症状少，因此仅保留75% 记录有的
t_phen_raw2 <- t_phen_raw %>% tibble::rownames_to_column(var='sample_id') %>%
  left_join(corr_raw[,c('Vomiting','sample_id')], by = 'sample_id') %>%
  tibble::column_to_rownames(var='sample_id') %>%
  arrange(rownames(t_phen_raw)) 

all <- corr_raw  %>% tibble::column_to_rownames(var='sample_id') %>%
  subset(select=c(sub_symbol,sub_clinic)) %>% arrange(rownames(corr_raw)) %>%
  bind_cols(t_phen_raw2) 

#########correlation analysis #####

sp_corr <- corr.test(all, method  = "sp", adjust = "BH")

#test
#t<- cor.test(all$LLS, all$Thyrotropin, method = "sp")

r_mat <- sp_corr$r 
p_mat <- sp_corr$p

sym1 <- sub_symbol
sym2 <- sub_clinic
sym3 <- colnames(t_phen_raw2)

feature_list <- list(p = sym3,
                     d = sym1,
                     c = sym2
                     
)
# select q-value for ploting (BH)

comp_list_up <- list( pd = c('d','p'),
                       pc = c('c','p'),
                       dc = c('d','c')                      
)

tit_list <- list(pd = c("Proteins","AMS phenotypes"),
                 pc = c("Clinical indexs", "AMS phenotypes"),
                 dc = c("Clinical indexs",  "Proteins")               
)
r_result <- data.frame()
p_result <- data.frame()

comp_list <- comp_list_up
######### heatmap 
for(i in 1:length(comp_list) ){
  #i=2
  comp_var <- names(comp_list)[[i]]
  fea_i <- feature_list[[comp_list[[i]][1]]]
  fea_j <- feature_list[[comp_list[[i]][2]]]
  tit_i <- tit_list[[comp_var]]
  r_i <- r_mat[which(rownames(r_mat) %in% fea_i), which(colnames(r_mat) %in%  fea_j)]
  p_i <- p_mat[which(rownames(p_mat) %in% fea_i), which(colnames(p_mat) %in%  fea_j)]
  
  r_i <- t(r_i)
  p_i <- t(p_i)
  col_fun = circlize::colorRamp2(c(1, 0, -1), c("#AF0000", "lightyellow", "#285BBC"))

  print(i)
  p <- Heatmap(r_i, name = "Coefficient", col = col_fun, 
               rect_gp = gpar(type = "none"),
               row_title = tit_i[1], column_title = tit_i[2],
               row_title_gp = gpar(fontsize = 14, fontface = 2),
               column_title_gp = gpar(fontsize = 14, fontface = 2),
               row_names_gp = gpar(fontsize = 12, fontface = 2),
               column_names_gp = gpar(fontsize = 12, fontface = 2),
               cell_fun=function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "grey", fill = NA))
                 grid.circle(x = x, y = y, r = abs(r_i[i, j]) *1.5* min(unit.c(width, height)), #pc *1.5 
                             gp = gpar(fill = col_fun(r_i[i, j]), col = NA))
                 grid.text(sprintf(ifelse(p_i[i, j]<0.05, "*", "") ), x, y, 
                           gp = gpar(fontsize = 12,col = "black"))
               }  )
  print(p)
  r_frame <- data.frame(r_i, stringsAsFactors = F, check.names = F)
  p_frame <- data.frame(p_i, stringsAsFactors = F, check.names = F)
  
  convert <- function(d) {
    lapply(1:ncol(d), function(i) {
      d1 <- d %>% tibble::rownames_to_column(var = "feature1")
      d2 <- d[, c(i)]
      d3 <- data.frame(d1[c("feature1")],d2) %>% mutate(symbol = colnames(d)[i])
    }) %>% do.call('rbind', .)
  }
  
  r_raw = convert(r_frame) %>% rename(r = d2, feature2 = symbol)
  r_result <- rbind(r_result,r_raw)
  p_raw = convert(p_frame) %>% rename(p_value = d2, feature2 = symbol)
  p_result <- rbind(p_result,p_raw)
  print(i)
}

# output the correlation result
corr_results <- r_result %>% left_join(p_result, by = c("feature1","feature2"), keep = F)

write.csv(corr_results, file = "../results/prediction_3group_merge_q0.05_corr_BH_results.csv")

######## cytoscape inputdata 
data <- corr_results

data <- data %>% filter(p_value <= 0.05 )

r_label = c()
for(i in 1:nrow(data)){
  if(data[i, 2] < 0)
  {
    r_label = c(r_label, "neg")
  }
  else
  {
    r_label = c(r_label, "pos")
  }
}

data$r_label = r_label

write.table(data, file="../results/prediction_3feature_merge_inputada_BH_p0.05.txt", sep="\t", row.names=F, quote=F)

net_data <- data %>% subset(select = -c(p_value))
net_data[,'r'] <- abs(net_data[,'r'])

write.table(net_data, file="../results/prediction_3feature_merge_input_network_BH_p0.05.txt", sep="\t", row.names=F, quote=F)

save.image('prediction_BH_sig_corr.RData')
