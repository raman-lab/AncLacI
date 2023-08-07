library(ggplot2)
library(minpack.lm)
library(grid)
clonal <- read.csv("csvs/sequenced_clonal_data.csv", header=T, stringsAsFactors=F)

rep <- subset(clonal, clonal[,"avg_Log2_RepSorted2_enrich"] > 1)
ind <- subset(rep, rep[,"X1000"] > rep[,"X0"] * 1.1)
uni.ind <- as.data.frame(unique(ind[,"anc_id"]))
colnames(uni.ind)[1] <- "anc_id"

noind <- subset(rep, rep[,"X1000"] < rep[,"X0"] * 1.1)
uni.noind <- as.data.frame(unique(noind[,"anc_id"]))
colnames(uni.noind)[1] <- "anc_id"

norep <- subset(clonal, clonal[,"avg_Log2_RepSorted2_enrich"] < 1)
uni.norep <- as.data.frame(unique(norep[,"anc_id"]))
colnames(uni.norep)[1] <- "anc_id"

for(i in 1:length(uni.norep[,1]))
{
  ss <- subset(clonal, uni.norep[i,"anc_id"] == clonal[,"anc_id"])
  
  m <- as.data.frame(matrix(0,8,3))
  colnames(m) <- c("conc.IPTG","avg_Fl","sd_Fl")
  m[,1] <- c(0,0.1,0.5,1,5,25,100,1000)
  m[,1] <- m[,1] + 0.001
  
  for(j in 1:8)
  {
    m[j,2] <- mean(ss[,3+j]) 
    m[j,3] <- sd(ss[,3+j])
  }
  
  p <- ggplot(m, aes(x=conc.IPTG, y=avg_Fl)) +
    geom_point(fill="#D1D3D4", size=3, color="#D1D3D4", pch=21) +
    labs(x=expression(paste("[IPTG] (",mu,"M)")), y=expression(paste("RFU OD"^"-1"))) +
    geom_errorbar(aes(ymin=avg_Fl-sd_Fl, ymax=avg_Fl+sd_Fl), width=0.2, color="#D1D3D4") +
    ylim(5000,55000) +
    theme_classic() + theme(text=element_text(size=20, colour="black"), axis.text.x=element_text(colour="black"), 
                            axis.text.y=element_text(colour="black"), plot.title=element_text(size=20, colour="black", hjust=0.5)) +
    scale_x_log10() +
    ggtitle(paste0(uni.norep[i,1]))
  ggsave(paste0("pdfs/dose_response_curves/",ss[1,"anc_id_no_slash"],".eps"), plot=p, device=cairo_ps, width=12, height=8, units="cm")
}

for(i in 1:length(uni.noind[,1]))
{
  ss <- subset(clonal, uni.noind[i,"anc_id"] == clonal[,"anc_id"])
  
  m <- as.data.frame(matrix(0,8,3))
  colnames(m) <- c("conc.IPTG","avg_Fl","sd_Fl")
  m[,1] <- c(0,0.1,0.5,1,5,25,100,1000)
  m[,1] <- m[,1] + 0.001
  
  for(j in 1:8)
  {
    m[j,2] <- mean(ss[,3+j]) 
    m[j,3] <- sd(ss[,3+j])
  }
  
  p <- ggplot(m, aes(x=conc.IPTG, y=avg_Fl)) +
    geom_point(fill="#FDEA61", size=3, color="#FDEA61", pch=21) +
    labs(x=expression(paste("[IPTG] (",mu,"M)")), y=expression(paste("RFU OD"^"-1"))) +
    geom_errorbar(aes(ymin=avg_Fl-sd_Fl, ymax=avg_Fl+sd_Fl), width=0.2, color="#FDEA61") +
    ylim(5000,55000) +
    theme_classic() + theme(text=element_text(size=20, colour="black"), axis.text.x=element_text(colour="black"), 
                            axis.text.y=element_text(colour="black"), plot.title=element_text(size=20, colour="black", hjust=0.5)) +
    scale_x_log10() +
    ggtitle(paste0(uni.noind[i,1]))
  ggsave(paste0("pdfs/dose_response_curves/",ss[1,"anc_id_no_slash"],".eps"), plot=p, device=cairo_ps, width=12, height=8, units="cm")
}

for(i in 1:length(uni.ind[,1]))
{
  ss <- subset(clonal, uni.ind[i,"anc_id"] == clonal[,"anc_id"])
  
  m <- as.data.frame(matrix(0,8,3))
  colnames(m) <- c("conc.IPTG","avg_Fl","sd_Fl")
  m[,1] <- c(0,0.1,0.5,1,5,25,100,1000)
  m[,1] <- m[,1] + 0.001
  
  for(j in 1:8)
  {
    m[j,2] <- mean(ss[,3+j]) 
    m[j,3] <- sd(ss[,3+j])
  }
  
  n.init  <- 0.5
  k.init  <- 10
  f.init  <- m[8,2]
  b.init  <- m[1,2]
  x <- m[,1]
  y <- m[,2]
  
  y.nls  <- nlsLM(y ~ b+(f-b)*x^n/(k^n+x^n), start=list(n=n.init, k=f.init, f=f.init, b=b.init))
  b <- round(as.numeric(sprintf("%.3f",summary(y.nls)$coefficients[4])), digits=0)
  f <- round(as.numeric(sprintf("%.3f",summary(y.nls)$coefficients[3])), digits=0)
  k <- round(as.numeric(sprintf("%.3f",summary(y.nls)$coefficients[2])), digits=1)
  n <- round(as.numeric(sprintf("%.3f",summary(y.nls)$coefficients[1])), digits=2)
  
  x <- seq(0.001,1020.001,by=0.1)
  y.calc <- b+(f-b)*x^n/(k^n+x^n) 
  fit <- as.data.frame(cbind(x,y.calc))
  
  v.b <- grobTree(textGrob(paste0("Basal = ",b), x=0.1, y=0.95, hjust=0, gp=gpar(col="black", fontsize=14)))
  v.f <- grobTree(textGrob(paste0("Vmax = ",f), x=0.1, y=0.85, hjust=0, gp=gpar(col="black", fontsize=14)))
  v.k <- grobTree(textGrob(paste0("Kd = ",k), x=0.1, y=0.75, hjust=0, gp=gpar(col="black", fontsize=14)))
  v.n <- grobTree(textGrob(paste0("n = ",n), x=0.1, y=0.65, hjust=0, gp=gpar(col="black", fontsize=14)))
  
  p <- ggplot(m, aes(x=conc.IPTG, y=avg_Fl)) +
    geom_point(fill="#D73027", size=3, color="#D73027", pch=21) +
    geom_line(data= fit, aes(x=x, y=y.calc), color="#D73027") +
    labs(x=expression(paste("[IPTG] (",mu,"M)")), y=expression(paste("RFU OD"^"-1"))) +
    geom_errorbar(aes(ymin=avg_Fl-sd_Fl, ymax=avg_Fl+sd_Fl), width=0.2, color="#D73027") +
    ylim(5000,55000) +
    theme_classic() + theme(text=element_text(size=20, colour="black"), axis.text.x=element_text(colour="black"), 
                            axis.text.y=element_text(colour="black"), plot.title=element_text(size=20, colour="black", hjust=0.5)) +
    scale_x_log10() +
    ggtitle(paste0(uni.ind[i,1])) +
    annotation_custom(v.b) +
    annotation_custom(v.f) +
    annotation_custom(v.k) +
    annotation_custom(v.n) 
  ggsave(paste0("pdfs/dose_response_curves/",ss[1,"anc_id_no_slash"],".eps"), plot=p, device=cairo_ps, width=12, height=8, units="cm")
}

