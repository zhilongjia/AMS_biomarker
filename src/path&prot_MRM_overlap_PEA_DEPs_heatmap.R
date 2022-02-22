library(dplyr)
library(pheatmap)

load("../data/PEA_comb_numeric.RData")

PEA_DEP <- data %>% tibble::rownames_to_column(var = 'symbol') %>% 
  filter(q_value < 0.05) %>% 
  select( c(symbol,log2FC)) %>% rename(PEA_log2FC = log2FC)

###########prot_heatmap ####
load(paste("../results/nAMS4k_nAMS1k_MRM_lograw_result.RData", sep =""))
nAMS_log_raw <- MRM_log_re
load(paste("../results/AMS4k_AMS1k_MRM_lograw_result.RData", sep =""))
#筛选机制组与保护组趋势相反蛋白，不需要过滤显著性，只看趋势
AMS_logFC <- MRM_log_re %>% select(c("symbol","log2FC")) %>% rename(AMS_log2FC = log2FC)

heat_raw <- nAMS_log_raw %>% filter(q_value < 0.05) %>% select(-c( "q_value","nAMS1k","nAMS4k")) %>%
  left_join(AMS_logFC, by ="symbol") %>% mutate(sig_A_nA = sign(log2FC*AMS_log2FC)) %>% filter(sig_A_nA == -1) %>%
  left_join(PEA_DEP, by = "symbol") %>% mutate(sig_PEA_A = sign(AMS_log2FC*PEA_log2FC)) %>% filter(sig_PEA_A == 1) %>% 
  select(-c("log2FC","PEA_log2FC","AMS_log2FC","sig_A_nA","sig_PEA_A")) %>%
  tibble::column_to_rownames(var = "symbol")

print(nrow(heat_raw))


save_pheatmap_pdf <- function(x, filename, width=18, height=3.5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
p <-pheatmap(
  heat_raw,
  # legend_breaks=seq(-1.5,1.5,1),
  width = 200,
  height = 200,
  fontsize_row = 12,
  fontsize_col = 12,
  scale = "row",
  cluster_cols  = FALSE,
  gaps_col = c(23),
  border_color = "black"
)
print(p)

save_pheatmap_pdf(p, filename = paste("../results/MRM_nAMS4k_nAMS1k_5DEPs_heatmap.pdf"))

##########  path_heatmap  #####################
load(paste("../results/AMS4k_AMS1k_MRM_lograw_result.RData", sep =""))


heat_raw <- MRM_log_re %>% filter(q_value < 0.05) %>% select(-c( "q_value","AMS1k","AMS4k")) %>%
  left_join(PEA_DEP, by = "symbol") %>% mutate(sig_PEA_A = sign(log2FC*PEA_log2FC)) %>% filter(sig_PEA_A == 1) %>% 
  select(-c("log2FC","PEA_log2FC","log2FC","sig_PEA_A")) %>%
  tibble::column_to_rownames(var = "symbol")

print(nrow(heat_raw))


save_pheatmap_pdf <- function(x, filename, width=18, height=3.5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
p <-pheatmap(
  heat_raw,
  # legend_breaks=seq(-1.5,1.5,1),
  width = 200,
  height = 200,
  fontsize_row = 12,
  fontsize_col = 12,
  scale = "row",
  cluster_cols  = FALSE,
  gaps_col = c(23),
  border_color = "black"
)
print(p)

save_pheatmap_pdf(p, filename = paste("../results/MRM_AMS4k_AMS1k_4DEPs_heatmap.pdf"))


