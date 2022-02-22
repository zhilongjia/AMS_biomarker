library(dplyr)
library(pheatmap)

load("../data/PEA_comb_numeric.RData")

PEA_DEP <- data %>% tibble::rownames_to_column(var = 'symbol') %>% 
  filter(q_value < 0.05) %>% 
  select( c(symbol,log2FC)) %>% rename(PEA_log2FC = log2FC)
load(paste("../results/AMS4k_AMS1k_MRM_lograw_result.RData", sep =""))

heat_raw_MRM <- MRM_log_re %>% select(-c( "q_value","AMS1k","AMS4k")) %>%
  left_join(PEA_DEP, by = "symbol") %>% mutate(sig_PEA_A = sign(log2FC*PEA_log2FC)) %>% filter(sig_PEA_A == 1) %>% 
  select(-c("log2FC","PEA_log2FC","log2FC","sig_PEA_A")) %>%
  tibble::column_to_rownames(var = "symbol")

########### 23 proteins heatmap from PEA and MRM of 60 samples #######
sam_ID <- read.csv("../data/PEA_sampleID_list.csv")

heat_raw_PEA <- data %>% tibble::rownames_to_column(var = "symbol") %>%
  filter(symbol %in% rownames(heat_raw_MRM)) %>% select(-c( "p_value","q_value","log2FC")) %>%
  tibble::column_to_rownames(var = "symbol")
theat_raw_PEA <- t(heat_raw_PEA)
theat_raw_PEA <- data.frame(theat_raw_PEA, stringsAsFactors = F)

# the whole 60 sample
sample_mrm <- read.csv("../data/mrm_sample_list_v4.csv")
sample_60 <- sample_mrm %>% filter(group %in% c("AMS4k","AMS1k")) %>% select(c("sampleID"))

theat_raw_MRM <- t(heat_raw_MRM)
theat_raw_MRM <- data.frame(theat_raw_MRM, stringsAsFactors = F)
theat_raw_MRM_sample <- theat_raw_MRM %>% tibble::rownames_to_column(var ="sampleID") %>%
  select(c("sampleID"))
theat_raw_MRM_data <-theat_raw_MRM %>% tibble::rownames_to_column(var ="sampleID") 

t_raw_60_sample_PEA_MRM <- theat_raw_PEA %>% tibble::rownames_to_column(var = "sample") %>% 
  left_join(sam_ID, by = "sample", keep =F) %>% select(-c("sample"))%>% right_join(sample_60, by = "sampleID", keep =F) %>%
  right_join(theat_raw_MRM_sample, by ="sampleID", keep=F) %>%
  arrange(sampleID) %>% left_join(theat_raw_MRM_data, by ="sampleID", keep =F) %>%
  tibble::column_to_rownames(var = "sampleID") 

raw_60_PEA_MRM <- t(t_raw_60_sample_PEA_MRM)
raw_60_PEA_MRM <- data.frame(raw_60_PEA_MRM, stringsAsFactors = F)

arr_raw <- raw_60_PEA_MRM %>% tibble::rownames_to_column(var ="symbol") %>% 
  arrange(symbol) %>% tibble::column_to_rownames(var ="symbol")
  
save_pheatmap_pdf <- function(x, filename, width=18, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
p <-pheatmap(
  arr_raw,
  # legend_breaks=seq(-1.5,1.5,1),
  width = 250,
  height = 80,
  fontsize_row = 12,
  fontsize_col = 12,
  scale = "row",
  cluster_cols  = FALSE,
  cluster_rows = FALSE,
  gaps_col = c(30),
  border_color = "black"
)
print(p)

save_pheatmap_pdf(p, filename = paste("../results/MRM_same_trend_PEA_AMS4k_AMS1k_23proteins_heatmap.pdf"))

