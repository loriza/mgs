library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(tidyverse)

dir1 = 'met_sp_ir_corr.txt'
dir2 = 'met_sp_ir_pval.txt'
c_list_dir = '../../../phenotype/ir_filter_list.txt'
anno_dir = 'lefse_sp.txt'

ca_comp_corr = readfile('caeca/ca_compound_corr.txt')

readfile <- function(path){
  df <- read_delim(path) %>% as.data.frame()
  row.names(df) = df[,1]
  df = df[,-1]
  return(df)
}

corr = read.table(dir1,row.names = 1,sep='\t',check.names = F,header = T) 
pval = read.table(dir2,row.names = 1,sep='\t',check.names = F,header = T)
collist = readLines(c_list_dir)
annot = read.delim(anno_dir,header = F,row.names = 1)


hm_corr_filter <- function(corr,pval,prefix,width,height,rcutoff=0,
                           anno=NULL,chosen_list=NULL,
                           collist=NULL,rowlist=NULL,
                           keyword=NULL,col_filt_na=TRUE,
                           row_filt_na = TRUE){
  
  
  
  if (!is.null(keyword)) { 
    s_corr <- corr[grep(keyword, names(corr), value=TRUE)] 
    s_pval <- pval[grep(keyword, names(pval), value=TRUE)] 
  }
  else{
    s_corr <- corr
    s_pval <- pval
  }
  
  if(!is.null(collist)) {
    subset_col <- function(list,df){
      filtered_df<- df[ ,which((names(df) %in% list)==TRUE)]
      return(filtered_df)
    }
    s_pval <- subset_col(collist,s_pval)
    s_corr <- subset_col(collist,s_corr)
  }
  if(!is.null(rowlist)) {
    
    subset_row <- function(list,df){
      filtered_df <- df[which((row.names(df) %in% list)==TRUE),]
      return(filtered_df)
    }
    s_pval <- subset_row(rowlist,s_pval)
    s_corr <- subset_row(rowlist,s_corr)
  }
  
  
  
  #s_corr = cont_clean(s_corr,0)
  #s_pval = cont_clean(s_pval,0)
  s_corr$id <- row.names(s_corr)
  s_corr <- reshape2::melt(s_corr)
  s_pval$id <- row.names(s_pval)
  s_pval <- reshape2::melt(s_pval)

  data <- cbind(s_corr,s_pval)
  data <- data[,-c(4,5)]
  names(data)[3] <- 'r'
  names(data)[4] <- 'p'
  
  
  
  datP <- data %>%  filter(abs(r) > rcutoff) %>%
    mutate(p = cut(x = p,  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                          labels = c("***", "**","*", ''))) %>% 
    dplyr::select(id, variable, p) %>% 
    spread(key = variable, value = p) %>% 
    column_to_rownames(var = "id") %>% 
    as.matrix() 

  datP <- apply(datP, 2, function(x) gsub("^$|^ $", NA, x))
  
  datCor <- data %>% select(id:r) %>% 
    spread(key = variable, value = r) %>% 
    column_to_rownames(var = "id") %>%
    as.matrix()
  
  
  if (isTRUE(col_filt_na)){
    datP <- datP %>% 
      as.matrix() %>% 
      .[,!apply(is.na(.), 2, all)]
    datCor <- datCor[,(colnames(datCor)%in% colnames(datP))]
  }
  
  
  if (isTRUE(row_filt_na)){
    datP <- datP %>% .[!apply(is.na(.), 1, all),]
    datCor <- datCor %>% .[(rownames(.) %in% rownames(datP)),]
  }
  
  
  datP[is.na(datP)] <- " "
  
  TextFunc <- function(dat, col = "black", fontsize = 15, numdat = TRUE,digit = 2){
    if(numdat == TRUE){
      function(j, i, x, y, width, height, fill){
        grid.text(round(dat, digit)[i, j], x, y, gp = gpar(fontsize = fontsize, col  = col))
      }
    }else{
      function(j, i, x, y, width, height, fill){
        grid.text(dat[i, j], x, y, gp = gpar(fontsize = fontsize, col  = col))
      }
    }}
  
  col_fun <- colorRamp2(breaks = c(-1, 0, 1), 
                        colors = c("#2b8cbe", "white", "#e41a1c"))
  
  datCor.sorted  <- datCor[order(rownames(datCor)),]
  if (!is.null(keyword)) { 
    pdf(paste0(prefix,"_",keyword,'.pdf'), width=width, height=height)
  }
  else{
    pdf(paste0(prefix,'.pdf'), width=width, height=height)
  }
  
  textLeg <- Legend(labels = c("P < 0.05", "P < 0.01", "P < 0.001"), title = "P.value", type = "points", 
                    pch = c("*", "**", "***"), legend_gp = gpar(col = "black"), background = "white")
  
  if (!is.null(anno)) { 
    anno <- data.frame(anno[chosen_list])
    names(anno)[1] <- 'group'
    anno$group <- as.factor(anno$group)
    
    anno <- anno[match(rownames(datCor.sorted),rownames(anno)),]
    
    row_anno <- HeatmapAnnotation(group=anno,which='row') 
    
    p1 <- Heatmap(datCor, name = "R2", col = col_fun,
                  rect_gp = gpar(col = "white", lwd = 1), 
                  right_annotation = row_anno,
                  cell_fun = TextFunc(datP, numdat = F))
    draw(p1, merge_legend = TRUE, heatmap_legend_side = "right", 
         annotation_legend_side = "right", annotation_legend_list = textLeg)
  }
  else{
    p1 <- Heatmap(datCor, name = "R2", col = col_fun,
                  rect_gp = gpar(col = "white", lwd = 1), 
                  cell_fun = TextFunc(datP, numdat = F)
    )
    draw(p1, merge_legend = TRUE, heatmap_legend_side = "right")
  }
  dev.off()

}

hm_corr_filter(corr = all_scfa_corr,pval = all_scfa_pval,
               prefix = 'all_scfa_0.6',
               chosen_list = 1,
               rcutoff = 0.6,
               #anno=annot,
              
               #collist = collist,
               
               width =5,height = 12)


