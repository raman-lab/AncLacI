setwd("//research.drive.wisc.edu/sraman4/General/Tony/ancestral/dna-binding/genomic_gfp/ngs/Zenodo/analysis")
s <- read.csv("sample_names", header=F, stringsAsFactors=F)

#LGF data
#----------------------------------------------------------------------------------------------------------------------
s.LGF <- subset(s, substr(s[,1],1,3) == "LGF")
lib.LGF <- read.csv("LGF_ref_lib.csv", header=F, stringsAsFactors=F)
colnames(lib.LGF) <- c("lib_id","aa_seq")

for(i in 1:length(lib.LGF[,1]))
{
  ss <- subset(lib.LGF, lib.LGF[,2]==lib.LGF[i,2])
  lib.LGF[i,"repeats"] <- length(ss[,1])
}

for(i in 1:length(s.LGF[,1]))
{
  temp <- read.csv(paste0(s.LGF[i,1],"/library_counts.csv"), header=F, stringsAsFactors=F)
  lib.LGF[,paste0(s.LGF[i,1])] <- round(temp[,2]/lib.LGF[,"repeats"],0)
}
write.csv(lib.LGF, file="LGF_raw_read_counts.csv", row.names=F, quote=F)

for(i in 4:length(lib.LGF))
{
  lib.LGF[,i] <- lib.LGF[,i]/sum(lib.LGF[,i])
}
write.csv(lib.LGF, file="LGF_relative_read_counts.csv", row.names=F, quote=F)

lib.LGF[,"enrichment_replicate1"] <- lib.LGF[,"LGF_repressed_sorted_replicate1"]/lib.LGF[,"LGF_presorted_replicate1"]
lib.LGF[,"enrichment_replicate2"] <- lib.LGF[,"LGF_repressed_sorted_replicate2"]/lib.LGF[,"LGF_presorted_replicate2"]
lib.LGF <- lib.LGF[,c("lib_id","aa_seq","enrichment_replicate1","enrichment_replicate2")]
write.csv(lib.LGF, file="LGF_fold_enrichment.csv", row.names=F, quote=F)

for(i in 1:length(lib.LGF[,1]))
{
  lib.LGF[i,"log2_enrichment_replicate1"] <- log((lib.LGF[i,"enrichment_replicate1"]+0.001),2)
  lib.LGF[i,"log2_enrichment_replicate2"] <- log((lib.LGF[i,"enrichment_replicate2"]+0.001),2)
  lib.LGF[i,"avg_log2_enrichment"] <- mean(lib.LGF[i,"log2_enrichment_replicate1"],lib.LGF[i,"log2_enrichment_replicate2"])
  lib.LGF[i,"sd_log2_enrichment"] <- sd(c(lib.LGF[i,"log2_enrichment_replicate1"],lib.LGF[i,"log2_enrichment_replicate2"]))
}
lib.LGF <- lib.LGF[,c("lib_id","aa_seq","log2_enrichment_replicate1","log2_enrichment_replicate2",
                      "avg_log2_enrichment","sd_log2_enrichment")]
write.csv(lib.LGF, file="LGF_log2_enrichment.csv", row.names=F, quote=F)
cor.test(x=lib.LGF[,"log2_enrichment_replicate1"], y=lib.LGF[,"log2_enrichment_replicate2"], method="spearman")
#----------------------------------------------------------------------------------------------------------------------


#DMS data
#----------------------------------------------------------------------------------------------------------------------
s.DMS <- subset(s, substr(s[,1],1,3) == "DMS")
lib.DMS <- read.csv("DMS_ref_lib.csv", header=F, stringsAsFactors=F)
colnames(lib.DMS) <- c("lib_id","aa_seq")

for(i in 1:length(lib.DMS[,1]))
{
  ss <- subset(lib.DMS, lib.DMS[,2]==lib.DMS[i,2])
  lib.DMS[i,"repeats"] <- length(ss[,1])
}

for(i in 1:length(s.DMS[,1]))
{
  temp <- read.csv(paste0(s.DMS[i,1],"/library_counts.csv"), header=F, stringsAsFactors=F)
  lib.DMS[,paste0(s.DMS[i,1])] <- round(temp[,2]/lib.DMS[,"repeats"],0)
}
write.csv(lib.DMS, file="DMS_raw_read_counts.csv", row.names=F, quote=F)

for(i in 4:length(lib.DMS))
{
  lib.DMS[,i] <- lib.DMS[,i]/sum(lib.DMS[,i])
}
write.csv(lib.DMS, file="DMS_relative_read_counts.csv", row.names=F, quote=F)

lib.DMS[,"enrichment_replicate1"] <- lib.DMS[,"DMS_repressed_sorted_replicate1"]/lib.DMS[,"DMS_presorted_replicate1"]
lib.DMS[,"enrichment_replicate2"] <- lib.DMS[,"DMS_repressed_sorted_replicate2"]/lib.DMS[,"DMS_presorted_replicate2"]
lib.DMS[,"normalized_enrichment_replicate1"] <- lib.DMS[,"enrichment_replicate1"]/lib.DMS[1,"enrichment_replicate1"]
lib.DMS[,"normalized_enrichment_replicate2"] <- lib.DMS[,"enrichment_replicate2"]/lib.DMS[1,"enrichment_replicate2"]
lib.DMS <- lib.DMS[,c("lib_id","aa_seq","enrichment_replicate1","enrichment_replicate2",
                      "normalized_enrichment_replicate1","normalized_enrichment_replicate2")]
write.csv(lib.DMS, file="DMS_fold_enrichment.csv", row.names=F, quote=F)

lib.DMS <- lib.DMS[,c("lib_id","aa_seq","normalized_enrichment_replicate1","normalized_enrichment_replicate2")]
for(i in 1:length(lib.DMS[,1]))
{
  lib.DMS[i,"log2_normalized_enrichment_replicate1"] <- log((lib.DMS[i,"normalized_enrichment_replicate1"]+0.001),2)
  lib.DMS[i,"log2_normalized_enrichment_replicate2"] <- log((lib.DMS[i,"normalized_enrichment_replicate2"]+0.001),2)
  lib.DMS[i,"avg_log2_normalized_enrichment"] <- mean(lib.DMS[i,"log2_normalized_enrichment_replicate1"],
                                                      lib.DMS[i,"log2_normalized_enrichment_replicate2"])
  lib.DMS[i,"sd_log2_normalized_enrichment"] <- sd(c(lib.DMS[i,"log2_normalized_enrichment_replicate1"],
                                                   lib.DMS[i,"log2_normalized_enrichment_replicate2"]))
}
lib.DMS <- lib.DMS[,c("lib_id","aa_seq","log2_normalized_enrichment_replicate1","log2_normalized_enrichment_replicate2",
                      "avg_log2_normalized_enrichment","sd_log2_normalized_enrichment")]
write.csv(lib.DMS, file="DMS_log2_normalized_enrichment.csv", row.names=F, quote=F)
cor.test(x=lib.DMS[,"log2_normalized_enrichment_replicate1"], y=lib.DMS[,"log2_normalized_enrichment_replicate2"], 
         method="spearman")
#----------------------------------------------------------------------------------------------------------------------

