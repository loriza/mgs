library(UpSetR)        
library(VennDiagram) 



feature <- read.delim('prof/genUReads.tsv',row.names = 1,check.names = F) 
mapping <- read.delim('metadata/cat_var_clean2.tsv',row.names=1,encoding='UTF-8',check.names = F)



upset_list <- list()
for (ele in unique(mapping[[target]]){
  ele_samples <- rownames(mapping[which(mapping$target == ele),])
  ele_features <- feature[,which(colnames(feature) %in% ele_samples)]
  ele_features <- row.names(ele_features[rowSums(ele_features[])>0,])
  upset_list[[ele]] <- ele_features
}



png(filename = outfile,width = 480, height = 300, units = "mm",res=600)
upset(fromList(upset_list),  
      nsets = 100,     
      nintersects = 40,
      order.by = "freq",
      keep.order = F,
      mb.ratio = c(0.6,0.4),   
      text.scale = 2
)

dev.off()
inter <- get.venn.partitions(upset_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '|')
inter <- subset(inter, select = -..values.. )
inter <- subset(inter, select = -..set.. )
write.table(inter, outtab, row.names = FALSE, sep = ',', quote = FALSE)
