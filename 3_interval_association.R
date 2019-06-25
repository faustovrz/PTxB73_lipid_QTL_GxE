library(dplyr)
library(GenomicRanges)
library(ggbio)
library(ggpubr)

cob_wt <- read.csv("GWAS_GxE/Regions_Bare cob weightmeanTemp.csv")
field_wt <- read.csv("GWAS_GxE/Regions_Field weightmeanTemp.csv")
grain_wt <- read.csv("GWAS_GxE/Regions_Grain weight per hectaremeanTemp.csv")
DTA <- read.csv("GWAS_GxE/Regions_Days to anthesismeanTemp.csv")


cob_wt$trait <- "cob_wt"
field_wt$trait <- "field_wt"
grain_wt$trait <- "grain_wt"
DTA$trait <- "DTA"


cob_wt <- cob_wt %>% select(trait,chrom,start,end)
field_wt <- field_wt %>% select(trait,chrom,start,end)
grain_wt <- grain_wt %>% select(trait,chrom,start,end)
DTA <- DTA %>% select(trait,chrom,start,end)

GWAS_GxE <- rbind(cob_wt, field_wt, grain_wt,DTA) 
GWAS_GxE$chr <- as.integer(gsub("chr","",GWAS_GxE$chrom))

GWAS_GxE <- GWAS_GxE  %>% arrange(trait, chrom, start)

GWAS_GxE_gr <- GRanges(seqnames= GWAS_GxE$chr, 
                       IRanges(start=GWAS_GxE$start, end = GWAS_GxE$end),
                       chr = GWAS_GxE$chr,
                       trait = GWAS_GxE$trait,
                       lod = -3)
seqlevels(GWAS_GxE_gr) <- as.character(1:10)
GWAS_GxE_gr <- sort(GWAS_GxE_gr)
genome(GWAS_GxE_gr) <- "Zea_mays.AGPv3.22"

biparental_GxE<- read.csv("biparental_GxE/GxE_B73xPT_QTL_intervals_1.csv")
colnames(biparental_GxE)
biparental_GxE$chr

biparental_GxE_gr <- GRanges(seqnames= biparental_GxE$chr, 
                             IRanges(start=biparental_GxE$start, end = biparental_GxE$end),
                             chr = biparental_GxE$chr,
                             trait = biparental_GxE$trait,
                             lod = biparental_GxE$lod)
seqlevels(biparental_GxE_gr) <- as.character(1:10)
biparental_GxE_gr <- sort(biparental_GxE_gr)
genome(biparental_GxE_gr) <- "Zea_mays.AGPv3.22"

gr <- biparental_GxE_gr
df <- data.frame(seqnames=seqnames(gr),
                 starts=start(gr)-1,
                 ends=end(gr),
                 names=gr$trait)

write.table(df, file="interval_association/flower_lipid_QTL_GxE.bed", quote=F, sep="\t", row.names=F, col.names=F)


gr <- NULL

hits <- as.data.frame(GWAS_GxE_gr)[unique(olap$queryHits),]

GWAS_GxE_gr 


gr <- GWAS_GxE_gr
df <- data.frame(seqnames=seqnames(gr),
                 starts=start(gr)-1,
                 ends=end(gr),
                 names=gr$trait)

write.table(df, file="interval_association/GWAS_GxE.bed", quote=F, sep="\t", row.names=F, col.names=F)

gr <- NULL


hits_gr<- subsetByOverlaps(GWAS_GxE_gr, biparental_GxE_gr)
chr_GxE <- seqnames(hits_gr)
chr_GxE <- as.character(chr_GxE@values)


bgnd <-  ggbio::autoplot( GWAS_GxE_gr[seqnames(GWAS_GxE_gr) %in% chr_GxE],
                          aes(color = trait, fill = trait)) + xlim(0,350e6) +
  ggtitle("biparental_GxE Ross-Ibarra") +
  scale_x_continuous(labels=function(x)x/1000000) 

temp_plot <- ggbio::autoplot(hits_gr[seqnames(GWAS_GxE_gr) %in% chr_GxE],  
                             aes(color = trait, fill = trait)) +
  xlim(0,350e6) + ggtitle("Overlap") +
  scale_x_continuous(labels=function(x)x/1000000) 



gr <- sortSeqlevels(biparental_GxE_gr)
gr <- sort(gr)

components <- reduce(gr)

olap <- as.data.frame(findOverlaps(gr, components))

values(gr) <- cbind(values(gr), DataFrame(component = olap$subjectHits))

toplot <- as.data.frame(gr)
toplot$bar <- as.numeric(rownames(toplot))
toplot <- toplot %>% group_by(chr,component) %>% mutate(nbar = row_number())


alt_plot <- ggbio::autoplot(gr,  aes(color = lod, fill = lod)) + xlim(0,350e6) +
  ggtitle("biparental_GxE Rellan Alvarez") +
  scale_x_continuous(labels=function(x)x/1000000) +
  geom_text(aes(x = end + 5000000, y = toplot$nbar, 
                hjust =0, label = substr(trait, 1, 8)), size = 2)


out <- ggarrange(alt_plot@ggplot, temp_plot@ggplot, bgnd@ggplot,
                 heights = c(3.5, 3,3.5), widths = c(15, 15,15),
                 ncol = 1, nrow = 3, align = "v")

pdf(file="interval_association/interval_overlap.pdf")
print(out)
dev.off()

library(regioneR)

zm3 <- toGRanges("interval_association/workspace.bed")
seqlengths(zm3)
seqlengths(zm3) <- width(zm3)
seqlengths(zm3)
genome(zm3) <- "Zea_mays.AGPv3.22"


ref_trait <- levels(biparental_GxE_gr$trait)
query_trait <- levels(as.factor(GWAS_GxE_gr$trait))

pt <- list()

pdf(file = "interval_association/regioneR/overlap_n.pdf")
for (q_trait in query_trait) {
  A <- GWAS_GxE_gr[GWAS_GxE_gr$trait == q_trait]
  for (r_trait in ref_trait) {
    B <- biparental_GxE_gr[biparental_GxE_gr$trait == r_trait]
    
    pair <- paste(q_trait,r_trait, sep = "-")
    print(pair)
    
    permutation <- overlapPermTest(
      A = A,
      B = B,
      genome = zm3,
      ntimes = 10000,
      mc.cores = 8,
      alternative="greater",
      force.parallel=TRUE)
    
    pt[[pair]] <- permutation$numOverlaps
    
    plot(pt[[pair]], 
         plotType = "Tailed", 
         main = paste(q_trait,r_trait),
         xlab  = "overlap count")
  }
}
dev.off()

save(pt, file = "interval_association/regioneR/overlap_n.Rdata")

overlap.width <- function(A,B,...){sum(width(intersect(A,B)))}

no.overlap.randomize <- 
  function(A,...){
    randomizeRegions(A, allow.overlaps = FALSE, ...)
  }


overlap_width_PT <- function(A, B, per.chromosome = TRUE, alternative="auto", ...) {
  
  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(B)) stop("B is missing")
  alternative<-match.arg(alternative,c("less","greater", "auto"))
  
  
  B <- toGRanges(B)
  A <- toGRanges(A)
  
  return(permTest(A, B=B, 
                  randomize.function = no.overlap.randomize, 
                  evaluate.function = overlap.width,
                  evaluate.function.name = "overlap.width", 
                  alternative=alternative, ...))
  
}

pt <- list()

pdf(file = "interval_association/regioneR/overlap_width_GAT.pdf")
for (q_trait in query_trait) {
  A <- reduce(GWAS_GxE_gr[GWAS_GxE_gr$trait == q_trait])
  for (r_trait in ref_trait) {
    B <- biparental_GxE_gr[biparental_GxE_gr$trait == r_trait]
    
    pair <- paste(q_trait,r_trait, sep = "-")
    print(pair)
    
    permutation <- overlap_width_PT(
      A = A,
      B = B,
      genome = zm3,
      ntimes = 10000,
      mc.cores = 8,
      alternative="greater",
      force.parallel=TRUE)
    
    pt[[pair]] <- permutation$overlap.width
    
    plot(pt[[pair]], 
         plotType = "Tailed", 
         main = paste(q_trait,r_trait),
         xlab  = "overlap length")
  }
}

dev.off()

save(pt, file = "interval_association/regioneR/overlap_width_GAT.Rdata")

pt.summary <-function(pt) {
  m <- stringr::str_split_fixed(names(pt), "-", 2)
  
  pt_summary  <- data.frame(
    GxE_Ross_Ibarra = m[,1],
    GxE_Rellan_Alvarez = m[,2],
    observed = unlist(lapply(pt, function(x){x$observed})),
    permuted = unlist(lapply(pt, function(x){mean(x$permuted)})),
    pval = unlist(lapply(pt, function(x){x$pval}))
  )
  
  
  pt_summary$FDR <- p.adjust(pt_summary$pval, method = "fdr" )
  pt_summary$significant <- ifelse(pt_summary$FDR<0.05, "*", NA)
  pt_summary <- pt_summary[order(pt_summary$pval),]
  pt_summary
}

load(file = "interval_association/regioneR/overlap_n.Rdata")
write.csv(pt.summary(pt),file="interval_association/regioneR/overlap_n_summary.csv")

load(file = "interval_association/regioneR/overlap_width_GAT.Rdata")
write.csv(pt.summary(pt),file="interval_association/regioneR/overlap_width_GAT_summary.csv")

