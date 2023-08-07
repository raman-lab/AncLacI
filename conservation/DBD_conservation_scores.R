library(msa)
library(ggplot2)
d <- read.csv("csvs/enriched_LGF_DBDs.csv", header=F, stringsAsFactors=F)
align <- msa(d[,1], type="protein", method="ClustalOmega")
data(BLOSUM62)
con.sc <- as.data.frame(msaConservationScore(align, BLOSUM62)) 

top <- read.csv("csvs/laci_topology.csv", header=T, stringsAsFactors=F) 
top[,"Conservation_score"] <- con.sc[,1]
hel <- subset(top, substr(top[,"Topology"],1,1) == "H")

pdf(file="pdfs/conservation_scores.pdf", width=12, height=4, useDingbats=F)
ggplot(top, aes(x=Position, y=Conservation_score)) + 
  geom_point(color="black", size=1.5) +
  labs(y="Conservation score", x="Position") + 
  theme_classic() + theme(axis.ticks=element_blank(), text=element_text(size=16), axis.text.x=element_text(colour="black"), 
                          axis.text.y=element_text(colour="black"))
dev.off()

pdf(file="pdfs/boxplot_conservation_scores.pdf", width=5, height=4, useDingbats=F)
ggplot(hel, aes(x=factor(Topology, levels=c("H1","H2","H3","HH")), y=Conservation_score)) + 
  geom_boxplot() +
  labs(y="Conservation score", x="") + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6) + 
  theme_classic() + theme(axis.ticks=element_blank(), text=element_text(size=16), axis.text.x=element_text(colour="black"), 
                          axis.text.y=element_text(colour="black"))
dev.off()

