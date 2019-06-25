library(qtl)
library(dplyr)
library(GenomicRanges)
library(ggbio) 

set.seed(1234567890) 

### First measurement  lipids from flower file ###########################

sets <- c(1)
radix <- "B73xPT_QTL"
for(set in sets)  {
cross <- read.cross("csvs",".","input/B73xPT_genotype.csv",
                    na.strings="NA",
                    phefile = paste0("input/B73xPT_GxE_phenotype_",set,".csv"),
                    genotypes=c("A","B"))

cross <- convert2riself(cross)
cross <- calc.genoprob(cross, step=1)

thresh <- c()
lod_intervals <- list()
significant <-list()

n_pheno <- length(colnames(cross$pheno))

for (i in seq(2, n_pheno , 2)){
  trait <- gsub("_high$","", colnames(cross$pheno)[i])
  trait_low <- cross$pheno[i+1]
  print(trait)
  
  possibleError <- tryCatch(
    out_hk  <- scanone(cross,, pheno.col = i,
                    addcov = trait_low , method="hk"),
    error=function(e) e
  )
  
  if(inherits(possibleError, "error")){
    print(paste(trait,possibleError))
  }
  
  if(!inherits(possibleError, "error")){
    
    #REAL WORK
    
    #Estimate LOD cut-off using permutation; add 5% cut off to plots
    #perm <- scanone(cross, pheno.col=i, addcovar = cross$pheno[2], n.perm=100, method="hk")
    perm <- scanone(cross,, pheno.col = i,
                    addcov = trait_low,
                    n.perm=1000, method="hk", n.cluster = 4)
    cut_off <- summary(perm, alpha=0.05)[1]
    
    if(any(out_hk$lod > cut_off)) {
      
      thresh[trait] <- cut_off
      significant[[trait]] <- out_hk
      
      peaks <- out_hk %>% 
        filter(lod > cut_off) %>% 
        group_by(chr) %>% 
        summarise(max(lod))
      for (c in peaks$chr) {
        intvl <- lodint(out_hk, chr = c, drop = 1.5)
        lod_intervals[[trait]][[c]] <- intvl
      }
    }
  }
}

save(significant, file = paste0("GxE/", radix,"_scan_", set, ".Rdata"))
save(lod_intervals, file = paste0("GxE/", radix,"_intervals_", set, ".Rdata"))
save(thresh, file = paste0("GxE/", radix,"_threshold_", set, ".Rdata"))

# Plot scans ##################################################################
pdf(file = paste0("GxE/", radix,"_scan_", set, ".pdf"))

for (trait in  names(significant)) {
  out_hk <- significant[[trait]]
  cut_off <- thresh[trait]
  plot(out_hk, col=c("black"), lwd = 0.7, main=trait)
  abline(h=cut_off,col="red",lty=2) #5% cut off for hk
}

dev.off()

##### Output QTL intervals ##########################################################


qtls <- data.frame()
gQTL <- GRanges()

for(trait in names(significant)){
  groups <- names(lod_intervals[[trait]])
  
  for (g in groups){
    r <- significant[[trait]]
    r <- r[r$chr==g,]
    i <- lod_intervals[[trait]][[g]]
    
    markers <- r[grepl("_",rownames(r)),]
    markers$id <- rownames(markers)
    markers$nt <- as.numeric(gsub("\\d+_","",markers$id))
    
    left  <- tail(markers[markers$pos <= min(i$pos),], n = 1)
    right <- head(markers[markers$pos >= max(i$pos),],  n = 1)
    
    left$trait <- trait
    right$trait <- trait
    
    pair <- rbind(left,right)
    qtls <- rbind(qtls, pair)
    
    gr <- GRanges(seqnames= g, 
                  IRanges(start=left$nt, end = right$nt),
                  chr = g,
                  trait = trait,
                  lod = max(i$lod))
    gQTL <- c(gQTL,gr)
  }
}

gr <- sortSeqlevels(gQTL)
gr <- sort(gr)

components <- reduce(gr)

olap <- as.data.frame(findOverlaps(gr, components))

values(gr) <- cbind(values(gr), DataFrame(component = olap$subjectHits))




toplot <- as.data.frame(gr)
toplot$bar <- as.numeric(rownames(toplot))
toplot <- toplot %>% group_by(chr,component) %>% mutate(nbar = row_number())

write.csv(toplot, row.names=FALSE,
          file = paste0("GxE/", radix,"_intervals_", set, ".csv"))

toplot <- read.csv(paste0("GxE/", radix,"_intervals_", set, ".csv"))

pdf(file = paste0("GxE/", radix,"_intervals_", set, ".pdf"),
    height = 3, width = 15)

p <- ggbio::autoplot(gr,  aes(color = lod, fill = lod)) + xlim(0,350e6) +
  ggtitle("PTxB73 significant QTLs (1.5 LOD interval, n = 87)",
          paste0("for flowering time and phospholipids: ", 
                 "GxE (elevation) interaction, set ", set)) +
  scale_x_continuous(labels=function(x)x/1000000) +
  geom_text(aes(x = end + 5000000, y = toplot$nbar, 
                hjust =0, label = substr(trait, 1, 5)), size = 2)
print(p)
dev.off()
}

