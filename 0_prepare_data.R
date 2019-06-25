library(Hmisc)
library(tidyr)
library(dplyr)
library(qtl)
library(gplots)
library(tidyverse)
library(corrr)
library(igraph)
library(ggraph)


geno <- read.csv("raw01/B73xPT_genotype.csv") %>% arrange(id)

set1 <- read.csv("raw01/B73xPT_phenotype_set1.csv") %>% arrange(ROW)

set2 <- read.csv("raw01/B73xPT_phenotype_set2.csv") %>% arrange(ROW)


### set1:  flowering time and 1st set2s measurent ###########################

base_pheno <- geno %>% select(id) %>% filter(id !="") %>%
  left_join(set1, by = "id") %>%
  select(-ROW, -lipid_batch_number) %>% # remove row for ensuing summary calculations
  filter(grepl("LANML",id)) %>%
  filter(!is.na(FIELD)) %>%
  arrange(FIELD,id)

colnames(base_pheno)

group <- c("DG", "LPC", "PC", "SM", "TG")

for (g in group){
  base_pheno[paste0(g,"s")] <- base_pheno %>%
    ungroup() %>%
    select(starts_with(paste0(g,"_"))) %>% 
    rowSums(na.rm = TRUE)
}
 ncol( base_pheno)
 base_pheno$LPC_PC <- log2(base_pheno$LPCs/base_pheno$PCs)


mc <-cor(base_pheno[,7:73],use="pairwise.complete.obs")

pdf(file = "metabolite_modules.pdf")
gplots::heatmap.2(mc, density.info="none", trace="none")


flattenCorrMatrix <- function(cormat, pmat, nmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    n  = (nmat)[ut],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

#just the metabolites along all fields 
cor_m <- Hmisc::rcorr(as.matrix(base_pheno[,7:73])) 
flat_cor <- flattenCorrMatrix(cor_m$r, cor_m$P,cor_m$n)
n <- nrow(flat_cor)
p_thresh <-  0.05 / (n * (n - 1 ) / 2)

is_significant <- (flat_cor$p < p_thresh) & !is.na(flat_cor$cor) & flat_cor$n > 299
selected_pair <- flat_cor[is_significant & flat_cor$p ==0 & flat_cor$cor <0,]
write.csv(file = "quotients.csv", selected_pair)  

# adj <- flat_cor[,c("row","column","cor")]
# adj <- adj[is_significant & adj$cor > 0.5,]

# colnames(adj)[3] <- "weight"

#adj$weight <- 1 - adj$weight

# Convert correlations  to an undirected graph object

# g <- adj %>% 
#   graph_from_data_frame(directed = FALSE)
# 
# c <- components(g)
# modules <- list()
# modules$modLPC <- names(c$membership[c$membership ==1])
# modules$modPC1 <- names(c$membership[c$membership ==2])
# modules$modPC2 <- names(c$membership[c$membership ==3])
# modules$modTG  <- names(c$membership[c$membership ==4])
 
# Plot
# ggraph(g) +
#   geom_edge_link() +
#   geom_node_point() +
#   geom_node_text(aes(label = name), repel = TRUE) +
#   theme_graph()
# 
# dev.off()

# for (comp in names(modules)){
#   base_pheno[comp] <- base_pheno %>%
#     ungroup() %>%
#     select(one_of(modules[[comp]])) %>% 
#     rowSums(na.rm = TRUE)
# }
# qcomb <- combn(names(modules),2)
# 
# idx <- NULL
# for (idx in 1:ncol(qcomb)){
#   qname <- paste0(qcomb[1,idx],"_",qcomb[2,idx])
#   q <- base_pheno[[qcomb[1,idx]]] / base_pheno[[qcomb[2,idx]]] 
#   base_pheno[qname] <- log2(q)
# }
# 
# base_pheno$modLPC_modPC1



idx <- NULL
for (idx in 1:nrow(selected_pair)){
  num <- as.character(selected_pair$row[idx])
  dem <- as.character(selected_pair$column[idx])
  qname <- paste0(rownames(selected_pair)[idx], "_",
                  gsub("_.*","", num ),
                  gsub("_.*","", dem))
  q <- base_pheno[num]/base_pheno[dem]
  base_pheno[qname] <- log2(q)
}





combo_pheno <- base_pheno %>% 
  group_by(id) %>%
  select(-FIELD) %>%
  summarise_all(.funs = mean, na.rm = TRUE) %>%
  arrange(id)

write.csv(combo_pheno,"input/B73xPT_combined_phenotype_1.csv", row.names = FALSE)


highland_pheno <-  base_pheno %>%
  filter(FIELD == "high") %>%
  group_by(id) %>%
  select(-FIELD) %>%
  summarise_all(.funs = mean, na.rm = TRUE) %>%
  arrange(id)

write.csv(highland_pheno,"input/B73xPT_highland_phenotype_1.csv", row.names = FALSE)

lowland_pheno <-  base_pheno %>%
  group_by(id) %>%
  select(-FIELD) %>%
  summarise_all(.funs = mean, na.rm = TRUE) %>%
  arrange(id)

write.csv(lowland_pheno,"input/B73xPT_lowland_phenotype_1.csv", row.names = FALSE)


GxE_pheno <- base_pheno %>%
  group_by(id,FIELD) %>%
  # Means by id corresponding to multiple samples per field
  summarise_all(.funs = mean, na.rm = TRUE) %>%
  # Make a column per field (highland/lowland)
  # for Rqtl GxE interaction QTL analysis
  gather(metabolite, signal, -(id:FIELD)) %>%
  unite(met_field,metabolite,FIELD) %>%
  spread(met_field,signal) %>%
  arrange(id)

colnames(GxE_pheno)
write.csv(GxE_pheno,"input/B73xPT_GxE_phenotype_1.csv", row.names = FALSE)

