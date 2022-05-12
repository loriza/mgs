
args <- commandArgs(trailingOnly = TRUE)

prof_infile <- args[1]
group_file <- args[2]
outpath <- args[3]
envfile <- args[4]
chosen_group <- as.numeric(args[5])
dir.create(outpath,showWarnings=F)

library(vegan)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(ggsci)

sampledata <- read.table(prof_infile, header = TRUE, row.names=1,sep='\t')
env <- read.table(envfile, header=TRUE, row.names=1,sep='\t')
group <-  read.table(group_file, header = TRUE, colClasses=c("character"),sep='\t')

env <- t(env)

print(dim(sampledata))
print(dim(group))
print(dim(env))

names(group)[chosen_group] <- 'chosen_group'
sampledata <- t(sampledata)
sampledata <- decostand(sampledata,method = "hellinger")

dca_calc <- function(sampledata,group,env){
  group <- as.list(group)
  dca <- decorana(veg = sampledata)
  dca1 <- max(dca$rproj[,1])
  dca2 <- max(dca$rproj[,2])
  dca3 <- max(dca$rproj[,3])
  dca4 <- max(dca$rproj[,4])
  GL <- data.frame(DCA1 = c(dca1), DCA2 = c(dca2), DCA3 = c(dca3), DCA4 = c(dca4))
  GL
  rownames(GL) <- c("Gradient length")
  write.table(GL, file = paste(outpath,"dca.tsv",sep = '/'),sep='\t',quote=F)
}
envfit_calc <- function(test_obj,env){
  envfit <- envfit(test_obj,env,permutations  = 999)  
  r <- as.matrix(envfit$vectors$r) 
  p <- as.matrix(envfit$vectors$pvals)
  env.p <- cbind(r,p)
  colnames(env.p) <- c("r2","p-value")
  K <- as.data.frame(env.p)
  obj <- deparse(substitute(test_obj))
  filename <- paste(obj, 'envfit.tsv',sep='.')
  write.table(as.data.frame(env.p), file = paste(outpath,filename,sep = '/'),sep = '\t',quote=F)
  return(envfit)
}
cca_calc <- function(sampledata,group,env){

  cca <- cca(sampledata, env, scale = TRUE)
  ccascore <- scores(cca)
  #envfit_calc(cca,env)
  write.table(ccascore$sites, file = paste(outpath, "cca.sample.tsv",sep = '/'), quote=F)
  write.table(cca$CCA$biplot, file = paste(outpath, "cca.env.tsv",sep = '/'),quote=F)
  write.table(ccascore$species, file = paste(outpath, "cca.species.tsv", sep = '/'),quote=F)
  print('hello')
  CCAE <- as.data.frame(cca$CCA$biplot[,1:2])
  
  CCAS1 <- ccascore$sites[,1]*0.3
  CCAS2 <- ccascore$sites[,2]*0.3
  
  plotdata <- data.frame(rownames(ccascore$sites), CCAS1, CCAS2, group$chosen_group)
  colnames(plotdata) <- c("sample","CCAS1","CCAS2","group")
  
  cca1 <- round(cca$CCA$eig[1]/sum(cca$CCA$eig)*100,2)
  cca2 <- round(cca$CCA$eig[2]/sum(cca$CCA$eig)*100,2)
 
  p <- ggplot(plotdata, aes(CCAS1, CCAS2)) +
    geom_point(aes(color=group),size = 4) + 
    scale_color_locuszoom()+
    #stat_chull(geom = "polygon", aes(group = group, color = group), alpha = 0.1) +
    xlab(paste("CCA1 ( ",cca1,"%"," )", sep = "")) + 
    ylab(paste("CCA2 ( ",cca2,"%"," )", sep = "")) +
    geom_segment(data = CCAE, 
                 aes(x = 0, y = 0, xend = CCAE[,1], yend = CCAE[,2]),
                 colour = "black", size = 0.8,
                 arrow = arrow(angle = 30, length = unit(0.4, "cm"))) +
    geom_text_repel(data = CCAE, segment.colour = "black",
                    aes(x = CCAE[,1], y = CCAE[,2], 
                        label = rownames(CCAE)),
                    size=3)+
    geom_vline(aes(xintercept = 0), linetype = "dotted") +
    geom_hline(aes(yintercept = 0), linetype = "dotted") +
    theme(panel.background = element_rect(fill = "white", colour = "black"), 
          panel.grid = element_blank(),
          axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(colour = "black", size = 10, face = "bold", ),
          axis.title.y = element_text(colour="black", size = 10, face = "bold", ),
          axis.text = element_text(colour = "black", size = 8),
          legend.title = element_blank(),
          legend.text = element_text(size = 6), 
          legend.key = element_blank(),
          plot.title = element_text(size = 4, colour = "black", 
                                    face = "bold", hjust = 0.5))
  
    
  p
  ggsave(p,filename = paste(outpath,'cca_plot.png',sep='/') ,device='png',width = 110,height = 110,units = 'mm',dpi=600)
  return(cca)
}


rda_calc <- function(sampledata,group,env){
  rda <- rda(sampledata, env, scale = TRUE)
  rdascore <- scores(rda)
  envf <- envfit_calc(rda,env)
  write.table(rdascore$sites,file=paste(outpath,"rda.sample.tsv",sep='/'),sep='\t',quote=F)
  write.table(rdascore$species,file=paste(outpath,"rda.species.tsv",sep='/'),quote=F)

  #环境矢量
  env_adj <- envf
  #env_adj$vectors$pvals <- p.adjust(env_adj$vectors$pvals, method='bonferroni')
  env_score <- as.data.frame(scores(envf, display = "vectors"))
  env_adj <- data.frame(cbind(env_score,env_adj$vectors$r,env_adj$vectors$pvals))
  names(env_adj) <- c('RDA1','RDA2','r2','p.adj')
  write.table(env_adj,file=paste(outpath,"rda.env.tsv",sep='/'),quote=F)
  RDAE <- env_adj[which(env_adj$p.adj < 0.05),1:2]


  plotdata <- data.frame(rownames(rdascore$sites), (rdascore$sites)*0.2, group$chosen_group)
  colnames(plotdata) <- c("sample","RDAS1","RDAS2","group")
  
  rda1 <- round(rda$CCA$eig[1]/sum(rda$CCA$eig)*100,2)
  rda2 <- round(rda$CCA$eig[2]/sum(rda$CCA$eig)*100,2)
  
  print(plotdata)
  #RDA plot        
  p <-ggplot(plotdata, aes(RDAS1, RDAS2)) +
    geom_point(aes(color = group),size = 5) + 
    scale_color_locuszoom()+
    #stat_chull(geom = "polygon", aes(group = group, color = group, fill = group), alpha = 0.1)+ 
    xlab(paste("RDA1 ( ",rda1,"%"," )", sep = "")) + 
    ylab(paste("RDA2 ( ",rda2,"%"," )", sep = "")) +
    geom_segment(data = RDAE, aes(x = 0, y = 0, xend = RDAE[,1], yend = RDAE[,2]),
                 colour = "black", size = 0.8,
                 arrow = arrow(angle = 30, length = unit(1, "cm"))) +
    geom_text_repel(data = RDAE, segment.colour = "black",
                    aes(x = RDAE[,1], y = RDAE[,2], label = rownames(RDAE)),
                    size=3) +
    geom_vline(aes(xintercept = 0), linetype = "dotted") +
    geom_hline(aes(yintercept = 0), linetype = "dotted") +
    theme(panel.background = element_rect(fill = "white", colour = "black"), 
          panel.grid = element_blank(),
          axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(colour = "black", size = 10, face = "bold", ),
          axis.title.y = element_text(colour="black", size = 10, face = "bold", ),
          axis.text = element_text(colour = "black", size = 8),
          legend.title = element_blank(),
          legend.text = element_text(size = 6), 
          legend.key = element_blank(),
          plot.title = element_text(size = 4, colour = "black", 
                                    face = "bold", hjust = 0.5))
    
  p
  ggsave(p,filename = paste(outpath,'rda_plot.png',sep='/') ,device = 'png',width=110,height = 110,units='mm',dpi=600)
  return(rda)
}	

cca_calc(sampledata,group,env)
rda_calc(sampledata,group,env)
