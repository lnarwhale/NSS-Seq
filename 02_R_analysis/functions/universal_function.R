import_strgtf_to_TPMmatrix <- function(gtf_dir){
  gtf_sample <- list.files(gtf_dir,full.names = T)
  gtf_sample <- gtf_sample[grep(pattern = "*_str_quant.gtf", gtf_sample)]
  gtf_sample <- sub("_str_quant.gtf","",basename(gtf_sample))
  gtf_lis <- list()
  for (i in 1:length(gtf_sample)) {
    gtf_lis[[i]] <- rtracklayer::import(paste0(gtf_dir,gtf_sample[i],"_str_quant.gtf"))
  }
  names(gtf_lis) <- gtf_sample
  gtf_dflis <- lapply(gtf_lis, function(x){as.data.frame(x)})
  gtf_dflis_tpm <- lapply(gtf_dflis, function(x){
    df <- dplyr::select(x,"ref_gene_name","TPM")
    df <- na.omit(df)
    df$TPM <- as.numeric(df$TPM)
    df <- sum_duplicates_dplyr(df) %>% as.data.frame()
  })
  for (i in 1:length(gtf_dflis_tpm)) {
    colnames(gtf_dflis_tpm[[i]])[1] <- "gene"
  }
  all_gene_uni_df <- 
    unique(c(
      unlist(lapply(gtf_dflis_tpm, function(df) df$gene))  # 
    ) %>% 
      as.data.frame() %>% 
      `colnames<-`("gene"))  # 
  gtf_dflis_tpm_lis <- lapply(gtf_dflis_tpm, function(x){
    df <- merge(all_gene_uni_df,x,by="gene",all.x=T)
    df <- df[order(df$gene),]
    df[is.na(df)] <- 0
    return(df)
  })
  for (i in 1:length(gtf_dflis_tpm_lis)) {
    gtf_dflis_tpm_lis[[i]]$sample <- names(gtf_dflis_tpm_lis)[i]
  }
  all_gtf_dflis_tpm_r <- do.call(rbind,gtf_dflis_tpm_lis)
  all_gtf_dflis_matrix <- all_gtf_dflis_tpm_r %>%
    pivot_wider(names_from = sample, values_from = TPM) %>% as.data.frame()
  return(all_gtf_dflis_matrix)
}



#df_change
dis_mer_df <- function(df,uni){
  n_col <- (ncol(df)-1)/2+1
  df1 <- df[,1:n_col]
  df2 <- df[,c(1,(n_col+1):ncol(df))]
  colnames(df2) <- colnames(df1)
  re_df <- rbind(df1,df2)
  if (uni) {
    re_df <- unique(re_df)  
  }
  return(re_df)
} #20230626
name_merge_df <- function(df,name){
  df_name <- as.data.frame(rownames(df))
  colnames(df_name) <- name
  all <- cbind(df_name,df)
  return(all)
}
sum_duplicates_dplyr <- function(df) {
  if (!is.data.frame(df)) {
    stop("输入必须是数据框")
  }
  # 确保第一列是基因名称
  if (ncol(df) < 2) {
    stop("数据框至少需要两列")
  }
  library(dplyr)
  result <- df %>%
    group_by(`[[`(., 1)) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = 'drop')
  return(result)
}

#peakanno
gtf_anno <- function(gtf_lis){
  tmp_gr_lis <- lapply(gtf_lis, function(tmp){
    filtered_tmp <- tmp[!(is.na(tmp$TPM)), ]
    filtered_tmp$TPM <- as.numeric(filtered_tmp$TPM)
    filtered_tmp <- filtered_tmp[filtered_tmp$TPM>0,]
    tmp_gr <- GRanges(seqnames = filtered_tmp$seqnames,
                      ranges = IRanges(start = filtered_tmp$start,
                                       end = filtered_tmp$end),
                      strand = filtered_tmp$strand)
    return(tmp_gr)
  })
  tmp_peakAnno_df <- lapply(tmp_gr_lis, function(tmp_gr){
    tmp_peakAnno <- annotatePeak(tmp_gr, TxDb = txdb,  tssRegion =c(-0,0),
                                 annoDb = org, sameStrand = TRUE, verbose = FALSE,
                                 genomicAnnotationPriority = c("Exon", "Intron","Intergenic"))
    tmp_peakAnno_df <- as.data.frame(tmp_peakAnno)
    tmp_peakAnno_df$type <- sub("\\s*\\(.*", "", tmp_peakAnno_df$annotation)
    return(tmp_peakAnno_df)
  })
  tmp_peakAnno_df_selec <- lapply(tmp_peakAnno_df, function(x){
    res <- dplyr::select(x,"seqnames","start","end","strand","type")
    return(res)
  })
  for (i in 1:length(tmp_peakAnno_df_selec)) {
    tmp_peakAnno_df_selec[[i]]$sample <- rep(names(tmp_peakAnno_df_selec)[i],nrow(tmp_peakAnno_df_selec[[i]]))
  }
  freq_df_lis <- list()
  for (i in 1:length(tmp_peakAnno_df_selec)) {
    t_df <- tmp_peakAnno_df_selec[[i]]
    type_counts <- table(t_df$type)
    total <- sum(type_counts)
    percentage <- round((type_counts / total) * 100, 2)
    freq_df_lis[[i]] <- data.frame(
      symbol = names(percentage),
      percentage = percentage
    )
    freq_df_lis[[i]] <- dplyr::select(freq_df_lis[[i]],"symbol","percentage.Freq")
    freq_df_lis[[i]]$group <- names(tmp_peakAnno_df_selec)[i]
  }
  names(freq_df_lis) <- names(tmp_peakAnno_df_selec)
  names <- names(freq_df_lis)
  freq_df_lis_rbind <- do.call("rbind", freq_df_lis)
  freq_df_lis_rbind$group <- factor(freq_df_lis_rbind$group,
                                    levels = names[length(names):1])
  cols <- c("#FFE35E", "#98694E", "#9BCCE3","#C85C19")
  feature <- c("Exon","Distal Intergenic","Intron","Downstream")
  p <- ggplot(freq_df_lis_rbind, aes(x = group, y = percentage.Freq, fill = symbol), alpha = 0.7) +
    geom_bar(stat = "identity", width = 0.75, colour = "black", lwd = 0.125) +
    coord_flip() +
    scale_fill_manual("Features", values = cols) +
    theme_bw() +
    theme(text = element_text(size = 14, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(y = "Percentage", x = NULL)
  return(p)
}
bed_anno <- function(bed_lis){
  tmp_gr_lis <- lapply(bed_lis, function(tmp){
    tmp_gr <- GRanges(seqnames = tmp$seqnames,
                      ranges = IRanges(start = tmp$start,
                                       end = tmp$end),
                      strand = tmp$strand)
    return(tmp_gr)
  })
  tmp_peakAnno_df <- lapply(tmp_gr_lis, function(tmp_gr){
    tmp_peakAnno <- annotatePeak(tmp_gr, TxDb = txdb,  tssRegion =c(-0,0),
                                 annoDb = org, sameStrand = TRUE, verbose = FALSE,
                                 genomicAnnotationPriority = c("Exon", "Intron","Intergenic"))
    tmp_peakAnno_df <- as.data.frame(tmp_peakAnno)
    tmp_peakAnno_df$type <- sub("\\s*\\(.*", "", tmp_peakAnno_df$annotation)
    return(tmp_peakAnno_df)
  })
  tmp_peakAnno_df_selec <- lapply(tmp_peakAnno_df, function(x){
    res <- dplyr::select(x,"seqnames","start","end","strand","type")
    return(res)
  })
  for (i in 1:length(tmp_peakAnno_df_selec)) {
    tmp_peakAnno_df_selec[[i]]$sample <- rep(names(tmp_peakAnno_df_selec)[i],nrow(tmp_peakAnno_df_selec[[i]]))
  }
  freq_df_lis <- list()
  for (i in 1:length(tmp_peakAnno_df_selec)) {
    t_df <- tmp_peakAnno_df_selec[[i]]
    type_counts <- table(t_df$type)
    total <- sum(type_counts)
    percentage <- round((type_counts / total) * 100, 2)
    freq_df_lis[[i]] <- data.frame(
      symbol = names(percentage),
      percentage = percentage
    )
    freq_df_lis[[i]] <- dplyr::select(freq_df_lis[[i]],"symbol","percentage.Freq")
    freq_df_lis[[i]]$group <- names(tmp_peakAnno_df_selec)[i]
  }
  names(freq_df_lis) <- names(tmp_peakAnno_df_selec)
  names <- names(freq_df_lis)
  freq_df_lis_rbind <- do.call("rbind", freq_df_lis)
  freq_df_lis_rbind$group <- factor(freq_df_lis_rbind$group,
                                    levels = names[length(names):1])
  cols <- c("#FFE35E", "#98694E", "#9BCCE3","#C85C19")
  feature <- c("Exon","Distal Intergenic","Intron","Downstream")
  p <- ggplot(freq_df_lis_rbind, aes(x = group, y = percentage.Freq, fill = symbol), alpha = 0.7) +
    geom_bar(stat = "identity", width = 0.75, colour = "black", lwd = 0.125) +
    coord_flip() +
    scale_fill_manual("Features", values = cols) +
    theme_bw() +
    theme(text = element_text(size = 14, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(y = "Percentage", x = NULL)
  return(p)
}

#matrix
fpkm2tpm <- function(data1){
  data1_name <- as.data.frame(data1[,1])
  colnames(data1_name) <- colnames(data1)[1]
  fpkm <- data1[,-1]
  fpkm <- exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  data1 <- cbind(data1_name,fpkm)
  return(data1)
}
de_point <- function(data1,col){
  library(stringi)
  for (i in 1:nrow(data1)) {
    data1[i,col] <- str_sub(data1[i,col],1,str_locate(data1[i,col],"[.]")[1]-1)
  }
  return(data1)
}
de_1 <- function(data1,col){
  library(stringi)
  library(stringr)
  for (i in 1:nrow(data1)) {
    if (str_detect(data1[i,col],"[|]")) {
      data1[i,col] <- str_sub(data1[i,col],1,str_locate(data1[i,col],"[|]")[1]-1)
    }else if (str_detect(data1[i,col],"[|]")==F) {
      data1[i,col] <- data1[i,col]
    }
  }
  return(data1)
}
de_1_af <- function(data1,col){
  library(stringi)
  library(stringr)
  for (i in 1:nrow(data1)) {
    if (str_detect(data1[i,col],"[|]")) {
      data1[i,col] <- str_sub(data1[i,col],str_locate(data1[i,col],"[|]")[1]+1,str_length(data1[i,col]))
    }else if (str_detect(data1[i,col],"[|]")==F) {
      data1[i,col] <- data1[i,col]
    }
  }
  return(data1)
}
de_rowall_zero <- function(data1){
  for (i in 2:ncol(data1)) {
    data1[,i] <- as.numeric(data1[,i])
  }
  data2 <- data1[which(rowSums(data1[2:ncol(data1)])>0),]
  return(data2)
}
tss2loc <- function(species,typ,data1){
  options(bedtools.path="/data1/shenluemou/biosoft/anaconda3/bin/")
  library(bedtoolsr)
  if (species=="mouse") {
    if (typ=="ENSEMBL") {
      ref <- read.csv("/data1/shenluemou/database/gene_local/R/mouse_grcm38_ens_geneloc.gtf",header = F)
    }else if (typ=="SYMBOL") {
      ref <- read.csv("/data1/shenluemou/database/gene_local/R/mouse_grcm38_geneloc.gtf",header = F)
    }
  }else if (species=="macacaf") {
    if (typ=="ENSEMBL") {
      ref <- read.csv("/data1/shenluemou/database/gene_local/R/macaca_fac_ens_geneloc.gtf",header = F)
    }else if (typ=="SYMBOL") {
      ref <- read.csv("/data1/shenluemou/database/gene_local/R/macaca_fac_sys_geneloc.gtf",header = F)
    }
  }else if (species=="human") {
    if (typ=="ENSEMBL") {
      ref <- read.csv("/data1/shenluemou/database/gene_local/R/human_hg38_ens_geneloc.gtf",header = F) 
    }else if (typ=="SYMBOL") {
      ref <- read.csv("/data1/shenluemou/database/gene_local/R/human_hg38_sys_geneloc.gtf",header = F)
    }
  }
  tmp1 <-select_(data1,"seqnames","start","end")
  tmp2 <- bt.intersect(a = tmp1,b = ref,wa = T,wb = T)
  tmp3 <- select_(tmp2,"V1","V2","V7")
  colnames(tmp3) <- c("seqnames","dominant_ctss","gene")
  return(tmp3)
}
genesymol_all <- function(type,data1){
  library(clusterProfiler)
  library(AnnotationDbi)
  library(dplyr)
  colnames(data1) <- "gene"
  if (type=="ENSEMBL") {
    tmp1 <- bitr(data1$gene,fromType = "ENSEMBL",
                 toType = c("ENTREZID","SYMBOL"),org,drop = TRUE)
  }else if (type=="ENTREZID") {
    tmp1 <- bitr(data1$gene,fromType = "ENTREZID",
                 toType = c("ENSEMBL","SYMBOL"),org,drop = TRUE)
  }else if (type=="SYMBOL") {
    tmp1 <- bitr(data1$gene,fromType = "SYMBOL",
                                 toType = c("ENSEMBL","ENTREZID"),org,drop = TRUE)
    # if (species=="macacaf") {
    #   tmp1 <- bitr(data1$gene,fromType = "SYMBOL",
    #                toType = c("ENTREZID"),org,drop = TRUE)
    # }else if (species!="macacaf") {
    #   tmp1 <- bitr(data1$gene,fromType = "SYMBOL",
    #                toType = c("ENSEMBL","ENTREZID"),org,drop = TRUE)     
    # }
  }
  data2 <- select_(tmp1,"ENSEMBL","SYMBOL","ENTREZID")
  # if (species=="macacaf") {
  #   data2 <- select_(tmp1,"SYMBOL","ENTREZID")
  # }else if (species!="macacaf") {
  #   data2 <- select_(tmp1,"ENSEMBL","SYMBOL","ENTREZID")
  # }
  return(data2)
}
enrich_go <- function(species,list1){
  library(AnnotationDbi)
  library(clusterProfiler)
  if (species=="mouse") {
    library(org.Mm.eg.db)
    org <- org.Mm.eg.db
  }else if (species=="human") {
    library(org.Hs.eg.db)
    org <- org.Hs.eg.db
  }else if (species=="macacaf") {
    org <- loadDb("/data1/luofc/hanx/database/Enrichment/GO.db/org.Macaca_fascicularis.eg.sqlite")
  }
  go_re <- enrichGO(gene=list1,
                    OrgDb = org,
                    keyType = "ENTREZID",
                    ont="ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = T)
  return(go_re)
}
enrich_kegg <- function(species,list1){
  library(clusterProfiler)
  if (species=="mouse") {
    org <- "mmu"
  }else if (species=="human") {
    org <- "hsa"
  }
  re <- enrichKEGG(list1, 
             organism=org, 
             pvalueCutoff=0.01, 
             pAdjustMethod="BH",
             keyType="kegg")
  return(re)
}
tolog10 <- function(data1){
  data1_name <- as.data.frame(data1[,1])
  colnames(data1_name) <- colnames(data1)[1]
  data1 <- data1[,-1]
  data1 <- log10(data1+1)
  data1 <- cbind(data1_name,data1)
  return(data1)
}
de_rowatleast <- function(data1,num){
  library(dplyr)
  data1_name <- as.data.frame(data1[,1])
  colnames(data1_name) <- colnames(data1)[1]
  data1 <- as.data.frame(data1[,-1])
  data1_name <- data1_name[apply(data1>num, 1, any),]
  data1 <- data1[apply(data1>num, 1, any),]
  data1 <- cbind(data1_name,data1)
  return(data1)
}
de_rowall <- function(data1,num){
  library(dplyr)
  data1_name <- as.data.frame(data1[,1])
  colnames(data1_name) <- colnames(data1)[1]
  data1 <- as.data.frame(data1[,-1])
  data1_name <- data1_name[apply(data1>num, 1, all),]
  data1 <- data1[apply(data1>num, 1, all),]
  data1 <- cbind(data1_name,data1)
  return(data1)
}
select_col <- function(data1,colnum){
  for (i in 1:length(colnum)) {
    tmp_num <- colnum[i]
    tmp <- as.data.frame(data1[,tmp_num])
    colnames(tmp) <- colnames(data1)[tmp_num]
    if (i==1) {
      data2 <- tmp
    }else if (i>1) {
      data2 <- cbind(data2,tmp)
    }
  }
  return(data2)
}
cluster_re <- function(data,num){
  library(stats)
  colnames(data)[1] <- "gene_id"
  rownames(data) <- make.unique(data$gene_id)
  data$gene_id <- NULL
  data <- data[apply(data>1,1,any),]
  data_cl <- data %>% t() %>% scale(center = FALSE,scale = TRUE) %>%
    t() %>% hcut(method="pearson")
  x <- c(1.7, 2.7, 3.7, -3.7, -5.7,3.7, 5.7, 4.8, 10.3, -12.9, 13.8, 12.3)
  y <- c(1.2, 2.3, 3.2, -3.5, -3.2, 2.1,4.7, .8, 1.2, 11.5, 1.3, 3.2)
  plot(x, y, cex = 1, pch = 3,xlab ="x", ylab ="y",col ="black")
  data_cluster <- rect.hclust(data_cl,num)
  return(data_cluster)
}
edger_single_compare <- function(data1,pv,fc,qua){
  library(edgeR)
  group <- 1:2
  data1[,2] <- as.double(data1[,2])
  data1[,3] <- as.double(data1[,3])
  y <- DGEList(counts = data1[,2:3],genes = data1[,1],group = group)
  keep <- rowSums(cpm(y)>1) >= 1
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y_bcv <- y
  bcv <- 0.1
  et <- exactTest(y_bcv, dispersion = bcv ^ 2)
  if (qua==TRUE) {
    gene1 <- decideTestsDGE(et, p.value = pv, lfc = fc)  
    colnames(gene1) <- "Signifi"
    results <- cbind(y$genes,y$counts,et$table,gene1)
  }else {
    results <- cbind(y$genes,y$counts,et$table)
    results$thres <- 0
    results$thres <- ifelse(results$PValue < pv & results$logFC > fc,"1",
                            ifelse(results$PValue < pv & results$logFC < -fc,"-1", "0"))
  }
  colnames(results)[4] <- "log2FoldChange"
  colnames(results)[6] <- "padj"
  colnames(results)[7] <- "threshold"
  results <- results[order(results$padj),]
  rownames(results) <- make.unique(results$genes)
  results <- results[,-1]
  results$threshold <- as.character(results$threshold)
  results$threshold <- ifelse(results$threshold=="1","up",ifelse(results$threshold=="-1",
                                                                 "down", "no_dif"))
  re_up <- results[which(results$threshold=="up"),]
  re_down <- results[which(results$threshold=="down"),]
  re_change <- results[which(results$threshold!="no_dif"),]
  all_res <- list(results,re_up,re_down,re_change)
  names(all_res) <- c("all","up","down","change")
  return(all_res)
}
deseq2_multi_compare <- function(df,group,pv,fc){
  library(DESeq2)
  rownames(df) <- make.unique(df[,1])
  df <- df[,-1]
  df <- df[rowMeans(df)>1,]
  colData <- data.frame(row.names = colnames(df),group)
  dds <- DESeqDataSetFromMatrix(countData = df, 
                                colData = colData, 
                                design = ~ group)
  dds1 <- DESeq(dds, fitType = 'mean', 
                minReplicatesForReplace = 7, 
                parallel = FALSE) 
  res <- results(dds1)
  res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  res1$threshold <- ifelse(res1$padj>pv,"no",ifelse(abs(res1$log2FoldChange)<fc,"no",
                                                          ifelse(res1$log2FoldChange>fc,"up","down")))
  res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
  res1_up<- res1[which(res1$log2FoldChange >= fc & res1$padj < pv),]
  res1_down<- res1[which(res1$log2FoldChange <= -fc & res1$padj < pv),]
  res1_change <- rbind(res1_up,res1_down)
  all_res <- list(res1,res1_up,res1_down,res1_change)
  names(all_res) <- c("all","up","down","change")
  return(all_res)
}
ens_2_sym <- function(df,col){
  df_ens <- as.data.frame(df$gene_id)
  colnames(df_ens) <- "ENSEMBL"
  df_allsym <- genesymol_all("ENSEMBL",df_ens)
  colnames(df)[col] <- "ENSEMBL"
  df_res <- merge(df_allsym,df,by="ENSEMBL")
  return(df_res)
}
sym_2_all <- function(df,col){
  colnames(df)[col] <- "gene_id"
  df_ens <- as.data.frame(df$gene_id)
  colnames(df_ens) <- "SYMBOL"
  df_allsym <- genesymol_all("SYMBOL",df_ens)
  colnames(df)[col] <- "SYMBOL"
  df_res <- merge(df_allsym,df,by="SYMBOL")
  return(df_res)
}
enz_2_all <- function(df,col){
  colnames(df)[col] <- "gene_id"
  df_ens <- as.data.frame(df$gene_id)
  colnames(df_ens) <- "ENTREZID"
  df_allsym <- genesymol_all("ENTREZID",df_ens)
  colnames(df)[col] <- "ENTREZID"
  df_res <- merge(df_allsym,df,by="ENTREZID")
  return(df_res)
}
col_replace <- function(df1,col,species,fromtype,totype){
  colnames(df1)[col] <- fromtype
  tmp <- as.data.frame(df1[,col])
  colnames(tmp) <- "gene"
  all_sys <- genesymol_all(species,fromtype,tmp)
  chose_sys <- select(all_sys,fromtype,totype)
  all <- merge(chose_sys,df1,by=fromtype)
  all <- all[,-1]
  return(all)
}

#tran 
featurec2tpmAfpkm <- function(mycounts){
  colnames(mycounts)[1] <- "Geneid"
  mycounts_gene <- select_(mycounts,"Geneid")
  mycounts_gene$Geneid <- make.unique(mycounts_gene$Geneid)
  mycounts <- mycounts[,6:ncol(mycounts)]
  mycounts <- cbind(mycounts_gene,mycounts)
  rownames(mycounts)<-mycounts[,1]
  mycounts<-mycounts[,-1]
  kb <- mycounts[,1] / 1000
  countdata <- mycounts[,2:ncol(mycounts)]
  rpk <- countdata / kb
  tpm <- t(t(rpk)/colSums(rpk) * 1000000)
  fpkm <- t(t(rpk)/colSums(countdata) * 10^6)
  all <- list(tpm,fpkm)
  names(all) <- c("TPM","FPKM")
  return(all)
}
lis2df <- function(total_map_txt){
  combined_df <- data.frame()
  for (i in 1:length(total_map_txt)) {
    df <- total_map_txt[[i]]
    df_name <- names(total_map_txt)[i]
    df_with_name <- cbind(df,Sample = rep(df_name,nrow(df)))
    combined_df <- bind_rows(combined_df,df_with_name)
  }
  return(combined_df)
}
split_col <- function(df,split_symbol,coln,split_res){
  split_columns <- strsplit(df[, coln],split = split_symbol, fixed = TRUE)
  split_data <- as.data.frame(do.call(rbind, split_columns))
  names(split_data) <- split_res
  res <- cbind(split_data, df[-coln])
  return(res)
}


#plot
plot_clu2heatmap <- function(list1,raw_matrix,rown){
  colnames(raw_matrix)[1] <- "gene"
  for (i in 1:length(list1)) {
    tmp <- as.data.frame(rownames(as.data.frame(list1[[1]])))
    colnames(tmp) <- "gene"
    tmp2 <- merge(tmp,raw_matrix,by="gene")
    if (i==1) {
      all <- tmp2
    }else if (i>1) {
      all <- rbind(all,tmp2)
    }
  }
  p <- heat_map(all,F,F,rown)
  return(list(all,p))
}
plot_venn <- function(list1){
  library(ggvenn)
  p <- ggvenn(list1, 
         show_elements = FALSE,
         label_sep = "\n", stroke_size = 1.5,set_name_size = 5,
         text_size = 3)
  return(p)
}
violin_plot <- function(data1,xname,yname){
  library(reshape2)
  colnames(data1)[1] <- "gene"
  data2 <- melt(data1)
  a <- ggplot(data2,aes(x=variable,y=value))+
    geom_violin(aes(color=variable),trim=T)+
    geom_boxplot(width=0.01)+
    labs(x=xname,y=yname)
  return(a)
}
tpm2pca_plot_group <- function(data1,group){
  library(factoextra)
  data1 <- de_rowall_zero(data1)
  rownames(data1) <- make.unique(data1[,1])
  data1 <- data1[,-1]
  data1 <- t(data1)
  data.pca <- prcomp(data1, scale = T)
  p <- fviz_pca_ind(data.pca, col.ind=group, 
                    fill = group,
                    pointshape = 21,
                    mean.point=F,  
                    label = "none", 
                    addEllipses = T, 
                    legend.title="Groups",
                    pointsize = 5,
                    ellipse.type="none", 
                    ellipse.level=0.9)
  p <- p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  return(p)
}
tpm2pca_plot_group_eli <- function(data1,group){
  library(factoextra)
  data1 <- de_rowall_zero(data1)
  rownames(data1) <- make.unique(data1[,1])
  data1 <- data1[,-1]
  data1 <- t(data1)
  data.pca <- prcomp(data1, scale = T)
  p <- fviz_pca_ind(data.pca, col.ind=group, 
                    fill = group,
                    pointshape = 21,
                    mean.point=F,  #
                    label = "none", # 
                    addEllipses = T, # 
                    legend.title="Groups",
                    #ellipse.type="none", #
                    pointsize = 5,
                    ellipse.level=0.8)
  p <- p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  return(p)
}
tpm2pca_plot_group_mark <- function(data1,group){
  library(factoextra)
  data1 <- de_rowall_zero(data1)
  rownames(data1) <- make.unique(data1[,1])
  data1 <- data1[,-1]
  data1 <- t(data1)
  data.pca <- prcomp(data1, scale = T)
  p <- fviz_pca_ind(data.pca, col.ind=group, 
                    fill = group,
                    pointshape = 21,
                    mean.point=F,  #
                    label = "none", # 
                    #addEllipses = T, 
                    legend.title="Groups",
                    pointsize = 5)
                    #ellipse.type="confidence", # 
                    #ellipse.level=0.9
  #p <- p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  return(p)
}
PM2plot_box <- function(df,xname,yname){
  library(ggplot2)
  p <- ggplot(df,aes(variable,value))+
    stat_boxplot(geom="errorbar",width=0.1,size=0.5,position=position_dodge(0.6),color="black")+
    geom_boxplot(position=position_dodge(0.6),
                 size=0.5,
                 width=0.3,
                 fill="light blue",
                 color="black",
                 outlier.color = "black",
                 outlier.fill = "black",
                 outlier.shape = 19,
                 outlier.size = 1.5,
                 outlier.stroke = 0.5,
                 outlier.alpha = 45,
                 notch = F,
                 notchwidth = 0.5)+
    theme(axis.title = element_text(size=18),
          axis.text = element_text(size=14)) +
    labs(x = xname, y = yname)
  return(p)
}
heat_map <- function(data1,clu_col,clu_row,rown){
  library(pheatmap)
  rownames(data1) <- make.unique(data1[,1])
  data1 <- data1[,-1]
  if (ncol(data1)<=2) {
    data1 <- t(data1) %>% scale(center = FALSE) %>% t()
  }else if (ncol(data1)>2) {
    data1 <- t(data1) %>% scale() %>% t()
  }
  p <- pheatmap::pheatmap(data1,cluster_rows = clu_row,
                          cluster_cols = clu_col,border=FALSE,
                          show_rownames = rown)
  return(p)
}
plot_vac <- function(Dat,pv,fc){
  library(ggplot2)
  library(ggrepel)
  Dat <- na.omit(Dat)
  Gene <- rownames(Dat[1:20,])
  ggplot(Dat,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
    geom_point()+
    scale_color_manual(values=c("#00008B","#808080","#DC143C"),
                       labels=c(paste0("Down:",nrow(Dat[Dat$threshold=="down",])),
                                paste0("No:",nrow(Dat[Dat$threshold=="no",])),
                                paste0("Up:",nrow(Dat[Dat$threshold=="up",]))))+#确定点的颜色
    geom_text_repel(
      data = Dat[1:20,],
      aes(label = Gene))+
    theme_bw()+
    theme(
      legend.title = element_blank(),#
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank()
    )+
    ylab('-log10 (p-adj)')+
    xlab('log2 (FoldChange)') + 
    geom_vline(xintercept=c(-fc,fc),lty=3,col="black",lwd=0.5) +#
    geom_hline(yintercept = -log10(pv),lty=3,col="black",lwd=0.5) 
}
plot_MA <- function(df,fdr,fc){
  library(ggpubr)
  p <- ggmaplot(df, fdr=fdr, fc = fc, size = 0.4,
           palette = c("#B31B21", "#1465AC", "darkgray"),
           legend = "top", top = 20,
           font.label = c("bold", 11),
           font.legend = "bold",
           font.main = "bold",
           ggtheme = ggplot2::theme_minimal())+
    theme_bw()+
    theme(
      legend.title = element_blank(),#
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank()
    )
  return(p)
}
plot_go <- function(S4,outdir,name){
  df <- as.data.frame(S4)
  n1 <- nrow(df[which(df[,1]=="BP"),])
  n2 <- nrow(df[which(df[,1]=="CC"),])
  n3 <- nrow(df[which(df[,1]=="MF"),])
  n <- max(n1,n2,n3)
  p_dot <- dotplot(S4,split="ONTOLOGY")+facet_grid(ONTOLOGY~.,scale = "free")
  p_bar <- barplot(S4,split="ONTOLOGY")+facet_grid(ONTOLOGY~.,scale = "free")
  if (n>10) {
    ggsave(paste0(outdir,name,"_dot.png"),p_dot,width = 10,height = 15)
    ggsave(paste0(outdir,name,"_bar.png"),p_bar,width = 10,height = 15)
  }else if (n<10) {
    ggsave(paste0(outdir,name,"_dot.png"),p_dot,width = 10,height = 1.5*n)
    ggsave(paste0(outdir,name,"_bar.png"),p_bar,width = 10,height = 1.5*n)
  }
}
plot_kegg <- function(S4,outdir,name){
  df <- as.data.frame(S4)
  n <- nrow(df)
  p_dot <- dotplot(S4)
  p_bar <- barplot(S4)
  if (n >= 10) {
    ggsave(paste0(outdir,name,"_dot.png"),p_dot,width = 10,height = 10)
    ggsave(paste0(outdir,name,"_bar.png"),p_bar,width = 10,height = 10)
  }else if (n < 10) {
    ggsave(paste0(outdir,name,"_dot.png"),p_dot,width = 10,height = 1*n)
    ggsave(paste0(outdir,name,"_bar.png"),p_bar,width = 10,height = 1*n)
  }
}
plot_tpm2cor <- function(df){
  library(ggcorrplot)
  rownames(df) <- df[,1]
  df <- df[,-1]
  corr<-cor(df)
  p <- ggcorrplot(corr,
             outline.color = "white",
             ggtheme = ggplot2::theme_bw,
             colors = c("#6D9EC1", "white", "#E46726"),lab = T,title = "Corralation"
  )
  return(p)
}
plot_clustertree <- function(df1){
  rownames(df1) <- make.unique(df1[,1])
  df1 <- df1[,-1]
  df2 <- t(df1)
  dists <- dist(df2,method = "euclidean")
  hc <- hclust(dists,method = "ave")
  dend1 <- as.dendrogram(hc)
  # png("test.png",width = 1000,height = 1300)
  p <- plot(dend1, type = "rectangle", 
       ylab="Height",
       main="Cluster Dendrogram")
  # dev.off()
  return(p)
}
plot_boxplot <- function(df,xname,yname){
  colnames(df) <- c("symbol","value")
  stats <- df %>% group_by(symbol) %>% summarize(
    median = median(value),
    q1=quantile(value,0.25),
    q3=quantile(value,0.75)
  )
  p <- ggplot(df,aes(x=symbol,y=value))+
    stat_boxplot(geom="errorbar",width=0.1,size=0.5,position=position_dodge(0.6),color="blue")+
    geom_boxplot(position=position_dodge(0.6),
                 size=0.5,
                 width=0.3,
                 fill="light blue",
                 color="blue",
                 outlier.color = "light blue",
                 outlier.fill = "red",
                 outlier.shape = 19,
                 outlier.size = 1.5,
                 outlier.stroke = 0.5,
                 outlier.alpha = 45,
                 notch = F,
                 notchwidth = 0.5)+
    theme(axis.title = element_text(size=18),
          axis.text = element_text(size=14)) +
    labs(x = xname, y = yname)+
    geom_text(data=stats,aes(x=symbol,y=median,label=paste("median：",round(median,2))),vjust=-1.5,
              size=2,fontface="bold",color="black")+
    geom_text(data=stats,aes(x=symbol,y=q3,label=paste("q3：",round(q3,2))),vjust=-1.5,
              size=2,fontface="bold",color="red")+
    geom_text(data=stats,aes(x=symbol,y=q1,label=paste("q1：",round(q1,2))),vjust=-1.5,
              size=2,fontface="bold",color="blue")
  return(p)
}
plot_kegg_df <- function(df,pv){
  library(dplyr)
  df <- df[which(df$pvalue<pv),]
  df2 <- dplyr::select(df,c("Description","Count","pvalue","qvalue"))
  colnames(df2) <- c("PATHWAY","Gene","Pvalue","Qvalue")
  p <- ggplot(df2,aes(Pvalue,PATHWAY))+
    geom_point(aes(size=Gene,color=-1*log10(Qvalue))) +
    scale_color_gradient(low="green",high="red") +
    labs(color=expression(-log[10](Qvalue)),size="Gene",
         x="Pvalue",y="Pathway name",title = "Pathway enrichment") 
  return(p)
}
plot_stacked_chart_per_nor <- function(combined_df,xname,yname){
  colnames(combined_df) <- c("treatment","value","sample")
  combined_df$proportion <- combined_df$value/tapply(combined_df$value,combined_df$sample,
                                                     sum)[combined_df$sample]
  
  combined_df$sample <- factor(combined_df$sample, 
                               levels = unique(rev(combined_df$sample)))
  p <- ggplot(combined_df,aes(x=sample,y=proportion,fill=treatment))+
    geom_bar(stat = "identity",position = "fill") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(y=xname,fill=yname)+
    coord_flip()
  return(p)
}

plot_boxplot_limit <- function(df,xname,yname){
  colnames(df) <- c("symbol","value")
  stats <- df %>% group_by(symbol) %>% summarize(
    median = median(value),
    q1=quantile(value,0.25),
    q3=quantile(value,0.75)
  )
  p <- ggplot(df,aes(x=symbol,y=value))+
    stat_boxplot(geom="errorbar",width=0.1,size=0.5,position=position_dodge(0.6),color="blue")+
    geom_boxplot(position=position_dodge(0.6),
                 size=0.5,
                 width=0.3,
                 fill="light blue",
                 color="blue",
                 outlier.color = "light blue",
                 outlier.fill = "red",
                 outlier.shape = 19,
                 outlier.size = 1.5,
                 outlier.stroke = 0.5,
                 outlier.alpha = 45,
                 notch = F,
                 notchwidth = 0.5)+
    theme(axis.title = element_text(size=18),
          axis.text = element_text(size=14)) +
    labs(x = xname, y = yname)+
    geom_text(data=stats,aes(x=symbol,y=q3,label=paste("q3：",round(q3,2))),vjust=-1.5,
              size=2,fontface="bold",color="red")+
    coord_flip()
  return(p)
}

plot_boxplot_noq3 <- function(df,xname,yname){
  colnames(df) <- c("symbol","value")
  p <- ggplot(df,aes(x=symbol,y=value))+
    stat_boxplot(geom="errorbar",width=0.1,size=0.5,position=position_dodge(0.6),color="blue")+
    geom_boxplot(position=position_dodge(0.6),
                 size=0.5,
                 width=0.3,
                 fill="light blue",
                 color="blue",
                 outlier.color = "light blue",
                 outlier.fill = "red",
                 outlier.shape = 19,
                 outlier.size = 1.5,
                 outlier.stroke = 0.5,
                 outlier.alpha = 45,
                 notch = F,
                 notchwidth = 0.5)+
    theme(axis.title = element_text(size=18),
          axis.text = element_text(size=14)) +
    labs(x = xname, y = yname)+
    coord_flip()
  return(p)
}

plot_smooth_lis <- function(df,xname,yname){
  colnames(df)[1] <- "symbol"
  p <- list()
  for (i in 1:nrow(df)) {
    tmp <- data.frame(gene=colnames(df[,-1]),
                      expression=as.numeric(df[,-1][i,]))
    tmp$order <- 1:nrow(tmp)
    p1 <- ggplot(tmp,aes(x=order,y=expression)) +
      geom_smooth(method = "loess",se=FALSE) +
      scale_x_continuous(breaks = tmp$order, labels = tmp$gene) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #plot.background = element_rect(fill = "white"),
        #panel.background = element_rect(fill = "white")
      )+
      labs(x=xname,y=yname)+
      ggtitle(paste0("( ",df$symbol[i]," )"))
    p[[i]] <- p1
  }
  names(p) <- df$symbol
  return(p)
}
plot_deseq2_multi_compare <- function(df,group,pv,fc){
  library(DESeq2)
  rownames(df) <- make.unique(df[,1])
  df <- df[,-1]
  df <- df[rowMeans(df)>1,]
  colData <- data.frame(row.names = colnames(df),
                        group = factor(group, levels = unique(group)))
  dds <- DESeqDataSetFromMatrix(countData = df, 
                                colData = colData, 
                                design = ~ group)
  dds1 <- DESeq(dds, fitType = 'mean', 
                minReplicatesForReplace = 7, 
                parallel = FALSE) 
  res <- results(dds1)
  res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  res1 <- na.omit(res1)
  res1$threshold <- ifelse(res1$padj>pv,"no",ifelse(abs(res1$log2FoldChange)<fc,"no",
                                                    ifelse(res1$log2FoldChange>fc,"up","down")))
  res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
  res1_up<- res1[which(res1$log2FoldChange >= fc & res1$padj < pv),]
  res1_down<- res1[which(res1$log2FoldChange <= -fc & res1$padj < pv),]
  res1_change <- rbind(res1_up,res1_down)
  the_res <- list(res1,res1_up,res1_down,res1_change)
  names(the_res) <- c("all","up","down","change")
  if (nrow(res1_up)>15) {
    Dat_up <- res1_up[1:15,]
  }else{
    Dat_up <- res1_up
  }
  if (nrow(res1_down)>15) {
    Dat_down <- res1_down[1:15,]
  }else{
    Dat_down <- res1_down
  }
  Dat <- rbind(Dat_up,Dat_down)
  Dat$symbol <- rownames(Dat)
  res1$symbol <- rownames(res1)
  Dat[1,]
  p <- ggplot(res1,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
    geom_point()+
    scale_color_manual(values=c("#00008B","#808080","#DC143C"),
                       labels=c(paste0("Down:",nrow(res1[res1$threshold=="down",])),
                                paste0("No:",nrow(res1[res1$threshold=="no",])),
                                paste0("Up:",nrow(res1[res1$threshold=="up",]))))+#确定点的颜色
    geom_text_repel(
      data = Dat,
      aes(label = symbol))+
    theme_bw()+
    theme(
      legend.title = element_blank(),#
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank()
    )+
    ylab('-log10 (p-adj)')+
    xlab('log2 (FoldChange)') + 
    geom_vline(xintercept=c(-fc,fc),lty=3,col="black",lwd=0.5) +#
    geom_hline(yintercept = -log10(pv),lty=3,col="black",lwd=0.5) 
  p1 <- ggplot(res1,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
    geom_point()+
    scale_color_manual(values=c("#00008B","#808080","#DC143C"),
                       labels=c(paste0("Down:",nrow(res1[res1$threshold=="down",])),
                                paste0("No:",nrow(res1[res1$threshold=="no",])),
                                paste0("Up:",nrow(res1[res1$threshold=="up",]))))+#
    theme_bw()+#修改图片背景
    theme(
      legend.title = element_blank(),#
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank()
    )+
    ylab('-log10 (p-adj)')+#
    xlab('log2 (FoldChange)') + #
    geom_vline(xintercept=c(-fc,fc),lty=3,col="black",lwd=0.5) +#
    geom_hline(yintercept = -log10(pv),lty=3,col="black",lwd=0.5) 
  rlog <- assay(rlogTransformation(dds1)) %>% as.data.frame()
  rlog$symbol <- rownames(rlog)
  select_gene <- Dat$symbol %>% as.data.frame()
  colnames(select_gene) <- "symbol"
  select_gene$order <- 1:nrow(select_gene)
  toheat <- merge(select_gene,rlog,by="symbol")
  toheat <- toheat[order(toheat$order),]
  toheat$order <- NULL
  p_heat <- heat_map(toheat,FALSE,FALSE,TRUE)
  all_res <- list(the_res,p,p1,p_heat)
  names(all_res) <- c("res_df","valco_plot_mark","valco_plot","heat_map")
  return(all_res)
}
plot_fa_2_logo <- function(infile){
  seqs <- readDNAStringSet(infile)
  freq_matrix <- consensusMatrix(seqs, as.prob = TRUE)[1:4, ]
  seq_length <- ncol(freq_matrix)
  p <- ggseqlogo(freq_matrix, method='prob', seq_type='dna') +
    theme_logo() +
    ggtitle("Sequence Conservation Logo") +
    xlab("Position") +
    ylab("Information Content (bits)") +
    scale_x_continuous(
      breaks = seq(5, seq_length, by = 5),  # 
      labels = seq(5, seq_length, by = 5)    # 
    ) +
    theme(
      axis.text.x = element_text(
        size = 10,                   # 
        angle = 45,                   #
        hjust = 1,                    # 
        vjust = 1                     # 
      ),
      plot.title = element_text(hjust=0.5, face="bold", size=14)
    )
  return(p)
}

