library(ggplot2)
dms <- read.csv("csvs/DMS_ngs_read_counts.csv", header=T, stringsAsFactors=F)
for(i in 1:length(dms[,1]))
{
  dms[i,"position"] <- as.numeric(substr(dms[i,1],2,nchar(dms[i,1])-1))
}
wt <- dms[1,]
dms <- dms[2:length(dms[,1]),]

dms.pe <- as.data.frame(matrix(0,60,3))
dms.pe[,1] <- seq(1,60,by=1)
colnames(dms.pe) <- c("position","avg_entropy","sd_entropy")

for(i in 2:60)
{
  ss <- subset(dms, dms[,"position"] == i)
  ss <- rbind(wt,ss)
  for(j in 3:6)
  {
    ss[,j] <- ss[,j] / sum(ss[,j])
  }
  ss[,"e1"] <- ss[,"R1_RepSorted2"] * log((ss[,"R1_RepSorted2"] / ss[,"reporter_strain_1"]), 2)
  ss[,"e2"] <- ss[,"R2_RepSorted2"] * log((ss[,"R2_RepSorted2"] / ss[,"reporter_strain_2"]), 2)
  dms.pe[i,2] <- mean(c(sum(ss[,"e1"], na.rm=T),sum(ss[,"e2"], na.rm=T)))
  dms.pe[i,3] <- sd(c(sum(ss[,"e1"], na.rm=T),sum(ss[,"e2"], na.rm=T)))
}

anc.ngs <- read.csv("csvs/LGF_ngs_read_counts.csv", header=T, stringsAsFactors=F)
aa  <- c("R","H","K","D","E","S","T","N","Q","C","G","P","A","V","I","L","M","F","Y","W","_")
anc <- as.data.frame(matrix(0,59*21,6))
colnames(anc) <- c("position","identity","reporter_strain_1","reporter_strain_2","R1_RepSorted2","R2_RepSorted2")
c <- 0
for(i in 2:60)
{
  for(j in aa)
  {
    c <- c + 1
    anc[c,1] <- i
    anc[c,2] <- j
    ss <- subset(anc.ngs, substr(anc.ngs[,2],i-1,i-1) == j)
    for(k in 3:6)
    {
      anc[c,k] <- sum(ss[,k])
    }
  }
}

anc.pe <- as.data.frame(matrix(0,60,3))
anc.pe[,1] <- seq(1,60,by=1)
colnames(anc.pe) <- c("position","avg_entropy","sd_entropy")

for(i in 2:60)
{
  ss <- subset(anc, anc[,"position"] == i)
  for(j in 3:6)
  {
    ss[,j] <- ss[,j] / sum(ss[,j])
  }
  ss[,"e1"] <- ss[,"R1_RepSorted2"] * log((ss[,"R1_RepSorted2"] / ss[,"reporter_strain_1"]), 2)
  ss[,"e2"] <- ss[,"R2_RepSorted2"] * log((ss[,"R2_RepSorted2"] / ss[,"reporter_strain_2"]), 2)
  anc.pe[i,2] <- mean(c(sum(ss[,"e1"], na.rm=T),sum(ss[,"e2"], na.rm=T)))
  anc.pe[i,3] <- sd(c(sum(ss[,"e1"], na.rm=T),sum(ss[,"e2"], na.rm=T)))
}

dms.pe[,"screen"] <- "dms"
anc.pe[,"screen"] <- "anc"

pdf(file="pdfs/dms_KLD_lineplot.pdf", width=12, height=1.25, useDingbats=F)
ggplot(pe, aes(x=position, y=avg_entropy)) + 
  geom_line(color="#D73027", size=0.7) +
  geom_point(color="#D73027", size=1.2) +
  geom_errorbar(aes(ymin=avg_entropy-sd_entropy, ymax=avg_entropy+sd_entropy), width=0.4) +
  geom_hline(yintercept=0.13, color="black", size=0.7, linetype="dashed") +
  labs(y="Kullback-Leibler divergence", x="Position") + 
  theme_classic() + theme(axis.ticks=element_blank(), text=element_text(size=10), axis.text.x=element_text(colour="black"), 
                          axis.text.y=element_text(colour="black"))
dev.off()

pdf(file="pdfs/anc_KLD_lineplot.pdf", width=12, height=1.25, useDingbats=F)
ggplot(anc.pe, aes(x=position, y=avg_entropy)) + 
  geom_line(color="#4575B4", size=0.7) +
  geom_point(color="#4575B4", size=1.2) +
  geom_errorbar(aes(ymin=avg_entropy-sd_entropy, ymax=avg_entropy+sd_entropy), width=0.4) +
  geom_hline(yintercept=0.75, color="black", size=0.7, linetype="dashed") +
  labs(y="Kullback-Leibler divergence", x="Position") + 
  theme_classic() + theme(axis.ticks=element_blank(), text=element_text(size=10), axis.text.x=element_text(colour="black"), 
                          axis.text.y=element_text(colour="black"))
dev.off()

#K-means clustering
#----------------------------------------------------------------------------------------------------
pe <- cbind(dms.pe[,1:3],anc.pe[,2:3])
colnames(pe) <- c("position","dms_avg_entropy","dms_sd_entropy","anc_avg_entropy","anc_sd_entropy")

pe.norm <- pe[,c(2,4)]
pe.norm[,1] <- scale(pe.norm[,1])
pe.norm[,2] <- scale(pe.norm[,2])

library(factoextra)
pdf(file="pdfs/kmeans_KLD_optimal_clusters.pdf", width=4, height=4, useDingbats=F)
fviz_nbclust(pe.norm, kmeans, method="wss") +
  geom_vline(xintercept=4, linetype=2)+
  labs(subtitle="Elbow method")
dev.off()

set.seed(123)
pe.norm.km <- kmeans(pe.norm, 4, nstart=100)
pe.norm.km$cluster

pdf(file="pdfs/kmeans_KLD.pdf", width=4, height=4, useDingbats=F)
fviz_cluster(pe.norm.km, data=pe.norm, palette=c("#F9F39A","#A6D499","#9EC1CF","#F26664"), geom="point", ellipse.type="convex", ggtheme=theme_bw(), 
             ylab="Phylogenetic screen standardized KLD", xlab="DMS screen standardized KLD")
dev.off()

quads <- as.data.frame(cbind(pe$position,pe.norm.km$cluster))
colnames(quads) <- c("position","cluster")
quads[,"group"] <- "group"

library(ggplot2)
pdf(file="pdfs/kmeans_KLD_quadrants_heatmap.pdf", width=12, height=1, useDingbats=F)
ggplot(quads, aes(as.factor(position), group)) +
  geom_tile(aes(fill=as.factor(cluster)), color="white") + 
  scale_fill_manual(values=c("#F9F39A","#A6D499","#9EC1CF","#F26664"), guide="legend", name="", na.value="white") +
  theme_gray(base_size=9) + 
  labs(y="", x="") + 
  theme_classic() + theme(axis.ticks=element_blank(), text=element_text(size=10), axis.text.x=element_text(colour="black"), 
                          axis.text.y=element_text(colour="black"))
dev.off()
#----------------------------------------------------------------------------------------------------


