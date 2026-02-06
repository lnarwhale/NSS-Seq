#---color choose---
#80B1DA blue dark
#F69A8F red
#7FC9B0 green
#7CCCED blue light light
#F69247 organe
#42B0D1 blue light
#AB87C4 purple
#F09BB8 pink
#999798 gray
#-------plot-TSS-corr-------
tocorrPlot_nomark <- function(data, cor, xlab, ylab) {
  p <- ggplot(data, aes(data[,2], y = data[,1])) +
    geom_hex(bins = 100, show.legend = FALSE) +
    scale_fill_gradient(low = "gray24", high = "gray56") +
    theme(text = element_text(size = 20,
                              family = "Helvetica"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.subtitle = element_text(size = 24, hjust = -0.125, vjust = 1, face="bold", color="black"),
          axis.text = element_text(size = 18, family = "Helvetica", color = "black"),
          legend.title = element_text(size = 16, family = "Helvetica", colour = "black"),
          legend.text = element_text(size = 14, family = "Helvetica", colour = "black"),
          axis.title = element_blank()) +
    xlab(xlab) +
    ylab(ylab) +
    xlim(c(0, 4.5)) +
    ylim(c(0, 4.5)) +
    annotate("text", x = 0, y = 4.5, size = 6, hjust = 0, label = paste("R = ", cor))  +
    coord_fixed(ratio = 1)
  return(p)
}
plot_cager_cor_nomark <- function(LICAGE,control,treatment){
  tag <- CTSSnormalizedTpmDF(LICAGE) %>% as.data.frame()
  filtr <- tag[,c(control,treatment)]
  filtr <- filtr[(filtr[, 1] >= 1 | filtr[, 2] >= 1), ]
  cor.l <- round(cor(filtr), digits = 3)
  filtr.data.log.l <- log10(filtr + 1)
  names(filtr) <- treatment
  names(filtr.data.log.l) <- treatment
  labels_y <- treatment
  labels_x <- control
  p <- tocorrPlot_nomark(data = filtr.data.log.l, xlab = labels_x, ylab = labels_y, 
                  cor = cor.l[2])
  return(p)
}
plot_all_correlations <- function(LICAGE,order_sam) {
  sampleLabels <- order_sam
  sampleLabels_x <- rev(rev(sampleLabels)[-1])
  sampleLabels_y <- rev(sampleLabels[-1])
  plots <- list()
  for (i in 1:length(sampleLabels_x)) {
    for (j in 1:(length(sampleLabels_y) - i + 1)) {
      control <- sampleLabels_y[i]
      treatment <- sampleLabels_x[j]
      p <- plot_cager_cor_nomark(LICAGE, control, treatment)
      plots[[length(plots) + 1]] <- p
    }
  }
  n <- length(sampleLabels_x)
  pushViewport(viewport(layout = grid.layout(n + 1, n + 1, 
                                             widths = unit(c(0.1, rep(1, n)), "null"),
                                             heights = unit(c(0.1, rep(1, n)), "null"))))
  i <- 1
  for (row in 1:n) {
    for (col in 1:(n - row + 1)) {
      print(paste0("Placing plot at row ", row, ", col ", col))
      v <- viewport(layout.pos.row = row + 1, layout.pos.col = col + 1)
      print(plots[[i]], vp = v)
      i <- i + 1
    }
  }
  i <- 1
  for (row in 1:n) {
    for (col in 1:(n - row + 1)) {
      v <- viewport(layout.pos.row = row + 1, layout.pos.col = col + 1)
      print(plots[[i]], vp = v)
      if (col == 1) {
        grid.text(sampleLabels_y[row], x = 0.5, y = 0.5, 
                  just = c("left", "centre"), gp = gpar(cex = 2.0), rot = 90,
                  vp = viewport(layout.pos.row = row + 1, layout.pos.col = 1))
      }
      if (row == 1) {
        grid.text(sampleLabels_x[col], x = 0.5, y = 0.5, 
                  just = c("centre", "top"), gp = gpar(cex = 2.0),
                  vp = viewport(layout.pos.row = 1, layout.pos.col = col + 1))
      }
      i <- i + 1
    }
  }
  upViewport()
}
plot_all_correlations_dui <- function(LICAGE,order_sam) {
  sampleLabels <- order_sam
  sampleLabels_x <- order_sam
  sampleLabels_y <- rev(order_sam)
  plots <- list()
  for (i in 1:length(sampleLabels_x)) {
    for (j in 1:(length(sampleLabels_y) - i + 1)) {
      control <- sampleLabels_y[i]
      treatment <- sampleLabels_x[j]
      p <- plot_cager_cor_nomark(LICAGE, control, treatment)
      plots[[length(plots) + 1]] <- p
    }
  }
  n <- length(sampleLabels_x)
  pushViewport(viewport(layout = grid.layout(n + 1, n + 1, 
                                             widths = unit(c(0.1, rep(1, n)), "null"),
                                             heights = unit(c(0.1, rep(1, n)), "null"))))
  i <- 1
  for (row in 1:n) {
    for (col in 1:(n - row + 1)) {
      print(paste0("Placing plot at row ", row, ", col ", col))
      v <- viewport(layout.pos.row = row + 1, layout.pos.col = col + 1)
      print(plots[[i]], vp = v)
      i <- i + 1
    }
  }
  i <- 1
  for (row in 1:n) {
    for (col in 1:(n - row + 1)) {
      v <- viewport(layout.pos.row = row + 1, layout.pos.col = col + 1)
      print(plots[[i]], vp = v)
      if (col == 1) {
        grid.text(sampleLabels_y[row], x = 0.5, y = 0.5, 
                  just = c("left", "centre"), gp = gpar(cex = 2.0), rot = 90,
                  vp = viewport(layout.pos.row = row + 1, layout.pos.col = 1))
      }
      if (row == 1) {
        grid.text(sampleLabels_x[col], x = 0.5, y = 0.5, 
                  just = c("centre", "top"), gp = gpar(cex = 2.0),
                  vp = viewport(layout.pos.row = 1, layout.pos.col = col + 1))
      }
      i <- i + 1
    }
  }
  upViewport()
}

#------PCA------
plot_pca_fa_segment <- function(data1, group, order,ord) {
  library(factoextra)
  library(ggplot2)
  library(dplyr)
  data1 <- de_rowall_zero(data1)
  rownames(data1) <- make.unique(data1[, 1])
  data1 <- data1[,-1]
  data1 <- t(data1)
  data.pca <- prcomp(data1, scale = TRUE)
  pca_data <- as.data.frame(data.pca$x[, 1:2])
  combined_data <- cbind(pca_data, Group = group)
  group_means <- aggregate(. ~ Group, data = combined_data, FUN = mean)
  group_means_ordered <- group_means %>%
    filter(Group %in% order) %>%
    arrange(match(Group, order))
  arrows_data <- data.frame(
    Group = order[-length(order)],
    PC1 = group_means_ordered$PC1[-nrow(group_means_ordered)],
    PC2 = group_means_ordered$PC2[-nrow(group_means_ordered)],
    xend = group_means_ordered$PC1[-1],
    yend = group_means_ordered$PC2[-1]
  )
  p <- fviz_pca_ind(data.pca, col.ind = group, 
                    fill = group,
                    pointshape = 21,
                    mean.point = FALSE, 
                    label = "none", 
                    addEllipses = FALSE, 
                    legend.title = "Groups",
                    pointsize = 5,
                    geom = "none") +
    geom_point(data = group_means, 
               aes(x = PC1, y = PC2, color = Group), 
               size = 5, 
               shape = 16) +
    geom_segment(data = arrows_data, 
                 aes(x = PC1, y = PC2, xend = xend, yend = yend), 
                 #linetype = "dashed", 
                 color = "orange",     
                 arrow = arrow(length = unit(0.05, "inches")),  
                 size = 0.5) +
    geom_rect(aes(xmin = min(combined_data$PC1), xmax = max(combined_data$PC1), 
                  ymin = min(combined_data$PC2), ymax = max(combined_data$PC2)), 
              fill = NA, color = "black") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    scale_color_discrete(breaks = ord) +
    scale_fill_discrete(breaks = ord)
  return(p)
}
pap_pca_group <- function(data1,group,order = NULL ){
  library(factoextra)
  data1 <- de_rowall_zero(data1)
  rownames(data1) <- make.unique(data1[,1])
  data1 <- data1[,-1]
  data1 <- t(data1)
  data.pca <- prcomp(data1, scale = T)
  group <- factor(group, levels = order)
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
#------plot-GO-----
plot_GO_bar <- function(df_lis,order_lis,top_num){
  samples <- names(df_lis)
  for (i in 1:length(df_lis)) {
    df_lis[[i]]$log10pvalue <- -log10(df_lis[[i]]$p.adjust)
    df_lis[[i]]$sample <- samples[i]
  }
  df_lis <- lapply(df_lis, function(x){
    res <- x[order(x$log10pvalue,decreasing = T),]
    return(res)
  })
  df_chose_des_lis <- lapply(df_lis, function(x){
    gene <- x[1:top_num,]$Description %>% as.data.frame()
    colnames(gene) <- "Description"
    return(gene)
  })
  des_res <- data.frame()
  for (i in 1:length(order_lis)) {
    des_res <- rbind(des_res,df_chose_des_lis[[order_lis[i]]])
  }
  des_res$order <- 1:nrow(des_res)
  unique_rows <- !duplicated(des_res$Description)
  des_res_unique <- des_res[unique_rows, ]
  des_res_unique$order <- NULL
  expre_lis <- lapply(df_lis, function(lis){
    dplyr::select(lis,"Description","log10pvalue")
  })
  res_lis <- list()
  for (i in 1:length(expre_lis)) {
    res_lis[[i]] <- merge(des_res_unique,expre_lis[[i]],
                          by="Description",all.x=TRUE)
    res_lis[[i]]$log10pvalue[is.na(res_lis[[i]]$log10pvalue)] <- 0
    res_lis[[i]] <- res_lis[[i]][match(des_res_unique$Description,
                                       res_lis[[i]]$Description),]
  }
  names(res_lis) <- names(expre_lis)
  long_df <- data.frame()
  for (name in names(res_lis)) {
    current_df <- res_lis[[name]]
    current_df$Source <- name
    long_df <- bind_rows(long_df, current_df)
  }
  long_df <- na.omit(long_df)
  long_df$log10pvalue <- as.numeric(long_df$log10pvalue)
  wide_df <- long_df %>%
    pivot_wider(names_from = Source, values_from = log10pvalue) %>% as.data.frame()
  wide_df <- na.omit(wide_df)
  heat_map_tmp <- function(data1){
    library(pheatmap)
    data1[,1] <- as.character(data1[,1])
    rownames(data1) <- make.unique(data1[,1])
    data1 <- data1[,-1]
    data1_scaled <- as.data.frame(scale(data1, center = TRUE, scale = TRUE))
    p <- pheatmap::pheatmap(data1_scaled, 
                            cluster_rows = F, 
                            cluster_cols = F, 
                            border = TRUE,  
                            show_rownames = T,
                            cellwidth = 15,  
                            cellheight = 15,  
                            main = "GO term",
                            legend = T,
                            name = "Z-score of -log₁₀(padj)",
                            angle_col =45,
                            fontsize_col = 8,
                            color = colorRampPalette(c("#004A80","white","#C85C19"))(100))
    return(p)
  }
  p <- heat_map_tmp(wide_df)
  return(p)
}
plot_GO_bar_cluT <- function(df_lis,order_lis,top_num){
  samples <- names(df_lis)
  for (i in 1:length(df_lis)) {
    df_lis[[i]]$log10pvalue <- -log10(df_lis[[i]]$p.adjust)
    df_lis[[i]]$sample <- samples[i]
  }
  df_lis <- lapply(df_lis, function(x){
    res <- x[order(x$log10pvalue,decreasing = T),]
    return(res)
  })
  df_chose_des_lis <- lapply(df_lis, function(x){
    gene <- x[1:top_num,]$Description %>% as.data.frame()
    colnames(gene) <- "Description"
    return(gene)
  })
  des_res <- data.frame()
  for (i in 1:length(order_lis)) {
    des_res <- rbind(des_res,df_chose_des_lis[[order_lis[i]]])
  }
  des_res$order <- 1:nrow(des_res)
  unique_rows <- !duplicated(des_res$Description)
  des_res_unique <- des_res[unique_rows, ]
  des_res_unique$order <- NULL
  expre_lis <- lapply(df_lis, function(lis){
    dplyr::select(lis,"Description","log10pvalue")
  })
  res_lis <- list()
  for (i in 1:length(expre_lis)) {
    res_lis[[i]] <- merge(des_res_unique,expre_lis[[i]],
                          by="Description",all.x=TRUE)
    res_lis[[i]]$log10pvalue[is.na(res_lis[[i]]$log10pvalue)] <- 0
    res_lis[[i]] <- res_lis[[i]][match(des_res_unique$Description,
                                       res_lis[[i]]$Description),]
  }
  names(res_lis) <- names(expre_lis)
  long_df <- data.frame()
  for (name in names(res_lis)) {
    current_df <- res_lis[[name]]
    current_df$Source <- name
    long_df <- bind_rows(long_df, current_df)
  }
  long_df <- na.omit(long_df)
  long_df$log10pvalue <- as.numeric(long_df$log10pvalue)
  wide_df <- long_df %>%
    pivot_wider(names_from = Source, values_from = log10pvalue) %>% as.data.frame()
  wide_df <- na.omit(wide_df)
  heat_map_tmp <- function(data1){
    library(pheatmap)
    data1[,1] <- as.character(data1[,1])
    rownames(data1) <- make.unique(data1[,1])
    data1 <- data1[,-1]
    data1_scaled <- as.data.frame(scale(data1, center = TRUE, scale = TRUE))
    p <- pheatmap::pheatmap(data1_scaled, 
                            cluster_rows = T, 
                            cluster_cols = T, 
                            border = TRUE,  
                            show_rownames = T,
                            cellwidth = 15, 
                            cellheight = 15,  
                            main = "GO term",
                            legend = T,
                            name = "Z-score of -log₁₀(padj)",
                            angle_col =45,
                            fontsize_col = 8,
                            color = colorRampPalette(c("#004A80","white","#C85C19"))(100))
    return(p)
  }
  p <- heat_map_tmp(wide_df)
  return(p)
}
plot_GO_bar_orderx <- function(df_lis, order_lis, top_num){
  #X - order according to order_lis
  samples <- names(df_lis)
  for (i in 1:length(df_lis)) {
    df_lis[[i]]$log10pvalue <- -log10(df_lis[[i]]$p.adjust)
    df_lis[[i]]$sample <- samples[i]
  }
  df_lis <- lapply(df_lis, function(x){
    res <- x[order(x$log10pvalue, decreasing = T),]
    return(res)
  })
  df_chose_des_lis <- lapply(df_lis, function(x){
    gene <- x[1:top_num,]$Description %>% as.data.frame()
    colnames(gene) <- "Description"
    return(gene)
  })
  des_res <- data.frame()
  for (i in 1:length(order_lis)) {
    des_res <- rbind(des_res, df_chose_des_lis[[order_lis[i]]])
  }
  des_res$order <- 1:nrow(des_res)
  unique_rows <- !duplicated(des_res$Description)
  des_res_unique <- des_res[unique_rows, ]
  des_res_unique$order <- NULL
  expre_lis <- lapply(df_lis, function(lis){
    dplyr::select(lis, "Description", "log10pvalue")
  })
  res_lis <- list()
  for (i in 1:length(expre_lis)) {
    res_lis[[i]] <- merge(des_res_unique, expre_lis[[i]],
                          by = "Description", all.x = TRUE)
    res_lis[[i]]$log10pvalue[is.na(res_lis[[i]]$log10pvalue)] <- 0
    res_lis[[i]] <- res_lis[[i]][match(des_res_unique$Description,
                                       res_lis[[i]]$Description),]
  }
  names(res_lis) <- names(expre_lis)
  long_df <- data.frame()
  for (name in names(res_lis)) {
    current_df <- res_lis[[name]]
    current_df$Source <- name
    long_df <- bind_rows(long_df, current_df)
  }
  long_df <- na.omit(long_df)
  long_df$log10pvalue <- as.numeric(long_df$log10pvalue)
  wide_df <- long_df %>%
    pivot_wider(names_from = Source, values_from = log10pvalue) %>% as.data.frame()
  wide_df <- na.omit(wide_df)
  ordered_samples <- names(df_lis)[order_lis]
  ordered_samples <- ordered_samples[ordered_samples %in% colnames(wide_df)]
  wide_df <- wide_df[, c("Description", ordered_samples)]
  heat_map_tmp <- function(data1){
    library(pheatmap)
    data1[,1] <- as.character(data1[,1])
    rownames(data1) <- make.unique(data1[,1])
    data1 <- data1[,-1]
    data1_scaled <- as.data.frame(scale(data1, center = TRUE, scale = TRUE))
    p <- pheatmap::pheatmap(data1_scaled, 
                            cluster_rows = F, 
                            cluster_cols = F,  
                            border = TRUE,  
                            show_rownames = T,
                            cellwidth = 15,
                            cellheight = 15,
                            main = "GO term",
                            legend = T,
                            name = "Z-score of -log₁₀(padj)",
                            angle_col = 45,
                            fontsize_col = 8,
                            color = colorRampPalette(c("#004A80","white","#C85C19"))(100))
    return(p)
  }
  p <- heat_map_tmp(wide_df)
  return(p)
}


#------plot_TSS_enrichment------
pap_peakanno_lis <- function(domTSS){
  library(ggplot2)
  library(CAGEr)
  library(ChIPseeker)
  library(RColorBrewer)
  peakAnno_list <- lapply(domTSS, 
                          function(x) annotatePeak(x, TxDb = txdb,  tssRegion = c(-1000, 1000), annoDb = org, sameStrand = TRUE, verbose = FALSE))
  names <- names(domTSS)
  names(peakAnno_list) <- names
  feats.l <- lapply(peakAnno_list, function(x) x@annoStat)
  names(feats.l) <- names
  new_feats.l <- list()
  for (i in 1:length(names)) {
    new_feats.l[[i]] <- data.frame(Feature = c("Promoter","5' UTR","3' UTR","Exon","Intron","Distal Intergenic"),
                                   Frequency = c(feats.l[[i]]$Frequency[1],feats.l[[i]]$Frequency[2],feats.l[[i]]$Frequency[3],
                                                 feats.l[[i]]$Frequency[4]+feats.l[[i]]$Frequency[5],
                                                 feats.l[[i]]$Frequency[6]+feats.l[[i]]$Frequency[7],
                                                 feats.l[[i]]$Frequency[8]))
  }
  names(new_feats.l) <- names(feats.l)
  for (i in 1:length(names)) {
    new_feats.l[[i]]$sample <- rep(names[[i]], nrow(new_feats.l[[i]]))
  }
  feats.df <- do.call("rbind", new_feats.l)
  feats.df$sample <- factor(feats.df$sample, levels = names[length(names):1])
  col <- c("#C9E3EF","#4298C5", "#F58228", "#FFECA9", "#98694E", "#D2BEA2")
  feature_factors <- c("Promoter","5' UTR","3' UTR",
                       "Exon","Intron","Distal Intergenic")
  names(col) <- feature_factors
  col_sel <- col[names(col) %in% unique(feats.df$Feature)]
  p <- ggplot(feats.df, aes(x = sample, y = Frequency, fill = Feature), alpha = 0.7) +
    geom_bar(stat = "identity", width = 0.75, colour = "black", lwd = 0.125) +
    coord_flip() +
    scale_fill_manual("Features", values = col_sel) +
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
pap_actg_line_from_dom <- function(dom,up,down,name_pre,res_dir){
  seq_lis <- lapply(dom, function(x)
    cluter2seq(up,down,x))
  sample_name <- names(dom)
  names(seq_lis) <- sample_name
  base_colors <- c("A" = "#C85C19", "G" = "#FFECA9", "C" = "#004A80", "T" = "#C9E3EF")
  for (i in 1:length(seq_lis)) {
    png(paste0(res_dir,"/",name_pre,"_",sample_name[i],".png"))
    legend_elements <- gList(
      textGrob(sample_name[i], x = unit(0.53, "npc"), y = unit(0.9, "npc"), just = "centre"),
      do.call(gList, lapply(names(base_colors), function(label) {
        col <- base_colors[label]
        grobTree(
          segmentsGrob(x0 = unit(0.8, "npc"), x1 = unit(0.9, "npc"), 
                       y0 = unit(0.88 - 0.03* which(names(base_colors) == label), "npc"), 
                       y1 = unit(0.88 - 0.03* which(names(base_colors) == label), "npc"),
                       gp = gpar(col = col, lwd = 2)),
          textGrob(label, x = unit(0.91, "npc"), 
                   y = unit(0.88 - 0.03 * which(names(base_colors) == label), "npc"), just = "left")
        )
      }))
    )
    print(plotPatternOccurrenceAverage(regionsSeq = seq_lis[[i]],
                                       patterns = c("A","G","C","T"), flankUp = up, flankDown = down,
                                       smoothingWindow = 1, color = c("#C85C19", "#FFECA9","#004A80","#C9E3EF"),
                                       xLabel = "Distance to TSS (bp)",
                                       plotLegend = "F"))
    grid.draw(legend_elements)
    dev.off()
    pdf(paste0(res_dir,"/",name_pre,"_",sample_name[i],".pdf"))
    print(plotPatternOccurrenceAverage(regionsSeq = seq_lis[[i]],
                                       patterns = c("A","G","C","T"), flankUp = up, flankDown = down,
                                       smoothingWindow = 1, color = c("#C85C19", "#FFECA9","#004A80","#C9E3EF"),
                                       xLabel = "Distance to TSS (bp)",
                                       plotLegend = "F"))
    grid.draw(legend_elements)
    dev.off()
  }
}
pap_heatmap_ud_10_from_dom <- function(dom,res_dir,name_pre){
  plot_heatmap_ud_10 <- function(seq_df,name){
    tmp <- as.character(seq_df)
    tmp_data <- matrix(NA, nrow = length(tmp), ncol = 20)
    for (i in seq_along(tmp)) {
      tmp_data[i, ] <- unlist(strsplit(tmp[i], ""))
    }
    column_title <- c("-10", "", "", "", "", "", "", "", "","-1", 
                      "+1", "", "", "", "", "", "", "", "","+10")
    tmp2 <- as.data.frame(tmp_data)
    tmp2 <- tmp2[rowSums(tmp2 == "N") == 0, ]
    colnames(tmp2) <- column_title
    tmp2 <- as.matrix(tmp2)
    heat <- Heatmap(tmp2, 
                    name = "Base",
                    col = c("A" = "#C85C19", "G" = "#FFECA9", "C" = "#004A80", "T" = "#C9E3EF"),
                    column_split = rep(1:2,each=10),
                    column_title = list(bottom = paste0(name," Position relative to Tss (bp)")),
                    use_raster=FALSE,show_row_names = F)
    return(heat)
  }
  seq_1010 <- lapply(dom, function(x)
    cluter2seq(10,10,x))
  names(seq_1010) <- names(dom)
  sample_name <- names(dom)
  plot_heatmap_1010_lis <- list()
  for (i in 1:length(seq_1010)) {
    plot_heatmap_1010_lis[[i]] <- plot_heatmap_ud_10(seq_1010[[i]],
                                                          sample_name[i])
    png(paste0(res_dir,"/",name_pre,"_",sample_name[i],".png"))
    print(plot_heatmap_1010_lis[[i]])
    dev.off()
    pdf(paste0(res_dir,"/",name_pre,"_",sample_name[i],".pdf"))
    print(plot_heatmap_1010_lis[[i]])
    dev.off()
  }
}
pap_dflis_2_freq_din <- function(df_lis,sample_order){
  seq_lis <- list()
  for (i in 1:length(df_lis)) {
    seq_lis[[i]] <- cluter2seq(1,1,df_lis[[i]])
  }
  names(seq_lis) <- names(df_lis)
  seq.l <- lapply(seq_lis, table)
  seq_count.l <- lapply(seq.l, function(x) x/sum(x) * 100)
  seq_count_df.l <- lapply(seq_count.l, function(x) 
    as.data.frame(x))
  seq_count_df_tidy.l <- lapply(seq_count_df.l, function(x)
    x[order(x$Freq,decreasing = T),])
  seq_count_all.df <- lapply(seq_along(seq_count_df_tidy.l), function(i) {
    df <- seq_count_df_tidy.l[[i]]
    df$sample <- names(seq_count_df_tidy.l)[i]
    return(df)
  }) %>% 
    bind_rows()
  colnames(seq_count_all.df) <- c("dinucleotide","percentage","samples")
  max_freq <- max(seq_count_all.df$percentage)
  upper_limit <- ceiling(max_freq / 5) * 5
  tmp <- seq_count_all.df
  seq_count_all.df <- seq_count_all.df %>% 
    arrange(percentage)
  col <- c("#F8F2AB", "#B4D6A4", "#67BC9A", "#13B0A5", "#0071A7", "#3E606F", "#88F9D4", "#18C29C",
           "#0B877D", "#126872", "#031727","#CD5555")
  p <- ggplot(data = seq_count_all.df, aes(x = factor(dinucleotide, levels = unique(seq_count_all.df$dinucleotide)),
                                           y = percentage,
                                           fill = factor(samples, levels = rev(sample_order)))) + 
    geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25) + 
    coord_flip() +
    xlab("Initiator dinucleotide") + ylab("Percentage") + 
    scale_fill_manual(values = col) +
    scale_y_continuous(limits = c(0, upper_limit), breaks = seq(0, upper_limit, by = 5)) +
    ggtitle(NULL) + 
    theme_bw() +
    theme(text = element_text(size = 14, colour = "black"), 
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(fill = NULL)+
    guides(fill = guide_legend(reverse = TRUE))
  res <- list(p,seq_count_df_tidy.l)
  names(res) <- c("picture","frequency_df")
  return(res)
}
pap_dflis_2_freq_din_orderCT <- function(df_lis, sample_order, order_CT = NULL){
  seq_lis <- list()
  for (i in 1:length(df_lis)) {
    seq_lis[[i]] <- cluter2seq(1, 1, df_lis[[i]])
  }
  names(seq_lis) <- names(df_lis)
  
  seq.l <- lapply(seq_lis, table)
  seq_count.l <- lapply(seq.l, function(x) x/sum(x) * 100)
  seq_count_df.l <- lapply(seq_count.l, function(x) 
    as.data.frame(x))
  seq_count_df_tidy.l <- lapply(seq_count_df.l, function(x)
    x[order(x$Freq,decreasing = T),])
  
  seq_count_all.df <- lapply(seq_along(seq_count_df_tidy.l), function(i) {
    df <- seq_count_df_tidy.l[[i]]
    df$sample <- names(seq_count_df_tidy.l)[i]
    return(df)
  }) %>% 
    bind_rows()
  
  colnames(seq_count_all.df) <- c("dinucleotide","percentage","samples")
  max_freq <- max(seq_count_all.df$percentage)
  upper_limit <- ceiling(max_freq / 5) * 5
  
  tmp <- seq_count_all.df
  if (!is.null(order_CT)) {
    all_dinucleotides <- unique(seq_count_all.df$dinucleotide)
    missing_in_order <- setdiff(all_dinucleotides, order_CT)
    
    if (length(missing_in_order) > 0) {
      warning(paste("以下dinucleotide不在order_CT参数中，将被放在最后:", 
                    paste(missing_in_order, collapse = ", ")))
      order_CT <- c(order_CT, missing_in_order)
    }
    seq_count_all.df$dinucleotide <- factor(
      seq_count_all.df$dinucleotide,
      levels = rev(order_CT)
    )
    seq_count_all.df <- seq_count_all.df %>%
      arrange(dinucleotide) %>%
      mutate(dinucleotide = as.character(dinucleotide))
    seq_count_all.df$dinucleotide <- factor(
      seq_count_all.df$dinucleotide,
      levels = unique(seq_count_all.df$dinucleotide)
    )
    
  } else {
    seq_count_all.df <- seq_count_all.df %>% 
      arrange(percentage)
  }
  
  col <- c("#F8F2AB", "#B4D6A4", "#67BC9A", "#13B0A5", "#0071A7", "#3E606F", "#88F9D4", "#18C29C",
           "#0B877D", "#126872", "#031727","#CD5555")
  p <- ggplot(data = seq_count_all.df, aes(x = factor(dinucleotide, levels = unique(seq_count_all.df$dinucleotide)),
                                           y = percentage,
                                           fill = factor(samples, levels = rev(sample_order)))) + 
    geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25) + 
    coord_flip() +
    xlab("Initiator dinucleotide") + ylab("Percentage") + 
    scale_fill_manual(values = col) +
    scale_y_continuous(limits = c(0, upper_limit), breaks = seq(0, upper_limit, by = 5)) +
    ggtitle(NULL) + 
    theme_bw() +
    theme(text = element_text(size = 14, colour = "black"), 
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(fill = NULL)+
    guides(fill = guide_legend(reverse = TRUE))
  
  res <- list(p, seq_count_df_tidy.l)
  names(res) <- c("picture", "frequency_df")
  return(res)
}
pap_df_2_freq_din_single <- function(df, sample_name){
  seq <- cluter2seq(1, 1, df)
  seq_table <- table(seq)
  seq_count <- seq_table / sum(seq_table) * 100
  seq_count_df <- as.data.frame(seq_count)
  colnames(seq_count_df) <- c("dinucleotide", "percentage")
  seq_count_df_tidy <- seq_count_df[order(seq_count_df$percentage, decreasing = TRUE), ]
  seq_count_df_tidy.l <- list()
  seq_count_df_tidy.l[[sample_name]] <- seq_count_df_tidy
  seq_count_all.df <- seq_count_df_tidy
  seq_count_all.df$samples <- sample_name
  colnames(seq_count_all.df) <- c("dinucleotide", "percentage", "samples")
  max_freq <- max(seq_count_all.df$percentage)
  upper_limit <- ceiling(max_freq / 5) * 5
  tmp <- seq_count_all.df
  seq_count_all.df <- seq_count_all.df %>% 
    arrange(percentage)
  col <- c("#F8F2AB", "#B4D6A4", "#67BC9A", "#13B0A5", "#0071A7", "#3E606F", "#88F9D4", "#18C29C",
           "#0B877D", "#126872", "#031727","#CD5555")
  sample_order <- c(sample_name)
  p <- ggplot(data = seq_count_all.df, aes(x = factor(dinucleotide, levels = unique(seq_count_all.df$dinucleotide)),
                                           y = percentage,
                                           fill = factor(samples, levels = rev(sample_order)))) + 
    geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25) + 
    coord_flip() +
    xlab("Initiator dinucleotide") + ylab("Percentage") + 
    scale_fill_manual(values = col[1]) +  # 只使用第一个颜色
    scale_y_continuous(limits = c(0, upper_limit), breaks = seq(0, upper_limit, by = 5)) +
    ggtitle(NULL) + 
    theme_bw() +
    theme(text = element_text(size = 14, colour = "black"), 
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +  # 单样本不需要图例
    labs(fill = NULL)
  res <- list(p, seq_count_df_tidy.l)
  names(res) <- c("picture", "frequency_df")
  return(res)
}
pap_df_2_freq_din_single <- function(df, sample_name,order_ct){
  seq <- cluter2seq(1, 1, df)
  seq_table <- table(seq)
  seq_count <- seq_table / sum(seq_table) * 100
  seq_count_df <- as.data.frame(seq_count)
  colnames(seq_count_df) <- c("dinucleotide", "percentage")
  seq_count_df_tidy <- seq_count_df[order(seq_count_df$percentage, decreasing = TRUE), ]
  seq_count_df_tidy.l <- list()
  seq_count_df_tidy.l[[sample_name]] <- seq_count_df_tidy
  seq_count_all.df <- seq_count_df_tidy
  seq_count_all.df$samples <- sample_name
  colnames(seq_count_all.df) <- c("dinucleotide", "percentage", "samples")
  max_freq <- max(seq_count_all.df$percentage)
  upper_limit <- ceiling(max_freq / 5) * 5
  tmp <- seq_count_all.df
  seq_count_all.df <- seq_count_all.df %>% 
    arrange(percentage)
  col <- c("#F8F2AB", "#B4D6A4", "#67BC9A", "#13B0A5", "#0071A7", "#3E606F", "#88F9D4", "#18C29C",
           "#0B877D", "#126872", "#031727","#CD5555")
  sample_order <- c(sample_name)
  p <- ggplot(data = seq_count_all.df, aes(x = factor(dinucleotide, levels = unique(seq_count_all.df$dinucleotide)),
                                           y = percentage,
                                           fill = factor(samples, levels = rev(sample_order)))) + 
    geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25) + 
    coord_flip() +
    xlab("Initiator dinucleotide") + ylab("Percentage") + 
    scale_fill_manual(values = col[1]) +  # 只使用第一个颜色
    scale_y_continuous(limits = c(0, upper_limit), breaks = seq(0, upper_limit, by = 5)) +
    ggtitle(NULL) + 
    theme_bw() +
    theme(text = element_text(size = 14, colour = "black"), 
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +  # 单样本不需要图例
    labs(fill = NULL)
  res <- list(p, seq_count_df_tidy.l)
  names(res) <- c("picture", "frequency_df")
  return(res)
}
pap_df_2_freq_din_single_ordCT <- function(df, sample_name, order_ct = NULL){
  order_ct <- rev(order_ct)
  seq <- cluter2seq(1, 1, df)
  seq_table <- table(seq)
  seq_count <- seq_table / sum(seq_table) * 100
  seq_count_df <- as.data.frame(seq_count)
  colnames(seq_count_df) <- c("dinucleotide", "percentage")
  seq_count_df_tidy <- seq_count_df[order(seq_count_df$percentage, decreasing = TRUE), ]
  seq_count_df_tidy.l <- list()
  seq_count_df_tidy.l[[sample_name]] <- seq_count_df_tidy
  seq_count_all.df <- seq_count_df_tidy
  seq_count_all.df$samples <- sample_name
  colnames(seq_count_all.df) <- c("dinucleotide", "percentage", "samples")
  max_freq <- max(seq_count_all.df$percentage)
  upper_limit <- ceiling(max_freq / 5) * 5
  if (!is.null(order_ct)) {
    all_dinucleotides <- unique(seq_count_all.df$dinucleotide)
    missing_in_order <- setdiff(all_dinucleotides, order_ct)
    
    if (length(missing_in_order) > 0) {
      warning(paste("以下dinucleotide不在order_ct参数中，将被放在最后:", 
                    paste(missing_in_order, collapse = ", ")))
      order_ct <- c(order_ct, missing_in_order)
    }
    seq_count_all.df$dinucleotide <- factor(
      seq_count_all.df$dinucleotide, 
      levels = rev(order_ct) 
    )
    seq_count_all.df <- seq_count_all.df[order(seq_count_all.df$dinucleotide, decreasing = TRUE), ]
  } else {
    tmp <- seq_count_all.df
    seq_count_all.df <- seq_count_all.df %>% 
      arrange(percentage)
  }
  col <- c("#F8F2AB", "#B4D6A4", "#67BC9A", "#13B0A5", "#0071A7", "#3E606F", "#88F9D4", "#18C29C",
           "#0B877D", "#126872", "#031727","#CD5555")
  sample_order <- c(sample_name)
  p <- ggplot(data = seq_count_all.df, aes(x = factor(dinucleotide, levels = unique(seq_count_all.df$dinucleotide)),
                                           y = percentage,
                                           fill = factor(samples, levels = rev(sample_order)))) + 
    geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25) + 
    coord_flip() +
    xlab("Initiator dinucleotide") + ylab("Percentage") + 
    scale_fill_manual(values = col[1]) +  
    scale_y_continuous(limits = c(0, upper_limit), breaks = seq(0, upper_limit, by = 5)) +
    ggtitle(NULL) + 
    theme_bw() +
    theme(text = element_text(size = 14, colour = "black"), 
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +  
    labs(fill = NULL)
  res <- list(p, seq_count_df_tidy.l)
  names(res) <- c("picture", "frequency_df")
  return(res)
}
pap_logo_from_dom <- function(dom,up,down,name_pre,res_dir){
  seq2logo <- function(seq){
    library("ggseqlogo")
    all <- c()
    for (i in 1:length(seq)) {
      all <- c(all,as.character(seq[[i]]))
    }
    n <- length(seq[[1]])/2
    p <- ggseqlogo(all,method="bits") +
      scale_x_continuous(breaks = seq(1, length(seq[[1]])), 
                         labels = seq(-n, n)[seq(-n, n) != 0])
    return(p)
  }
  seq_lis <- lapply(dom, function(x){
    cluter2seq(up,down,x)
  })
  sample_name <- names(seq_lis)
  for (i in 1:length(seq_lis)) {
    logo_seq <- seq2logo(seq_lis[[i]])
    ggsave(paste0(res_dir,"/",name_pre,"_",sample_name[i],".png"),logo_seq)
    ggsave(paste0(res_dir,"/",name_pre,"_",sample_name[i],".pdf"),logo_seq)
  }
}

#-----heatmap------
heat_map_cut <- function(data1){
  library(pheatmap)
  library(circlize)
  rownames(data1) <- make.unique(data1[,1])
  data1 <- data1[,-1]
  if (ncol(data1)<=2) {
    data1 <- t(data1) %>% scale(center = FALSE) %>% t()
  }else if (ncol(data1)>2) {
    data1 <- t(data1) %>% scale() %>% t()
  }
  p <- Heatmap(data1,cluster_columns  = FALSE,
               cluster_rows = FALSE,
               show_row_names = T,
               row_names_gp = gpar(fontsize = 5),
               #row_split = split_factor,
               col = colorRamp2(c(-3, 0, 3), c("#004A80", "#FFFFFF", "#C85C19")),
               #row_title = split_group ,  
               row_title_rot = 0, 
               border =T,
               column_names_side = "top",
               column_names_rot = 45,
               use_raster=F,
               name = "Expression Level")
  return(p)
}
heat_map_cut_noname <- function(data1){
  library(pheatmap)
  library(circlize)
  rownames(data1) <- make.unique(data1[,1])
  data1 <- data1[,-1]
  if (ncol(data1)<=2) {
    data1 <- t(data1) %>% scale(center = FALSE) %>% t()
  }else if (ncol(data1)>2) {
    data1 <- t(data1) %>% scale() %>% t()
  }
  p <- Heatmap(data1,cluster_columns  = FALSE,
               cluster_rows = FALSE,
               show_row_names = F,
               row_names_gp = gpar(fontsize = 5),
               #row_split = split_factor,
               col = colorRamp2(c(-3, 0, 3), c("#004A80", "#FFFFFF", "#C85C19")),
               #row_title = split_group ,  
               row_title_rot = 0, 
               border =T,
               column_names_side = "top",
               column_names_rot = 45,
               use_raster=F,
               name = "Expression Level")
  return(p)
}

#-----boxplot----
plot_boxplot_enhanced <- function(df, xname, yname,sample_order = NULL) {
  colnames(df) <- c("symbol", "value")
  df$symbol <- factor(df$symbol, levels = rev(sample_order))
  n_symbols <- length(unique(df$symbol))
  color_palette <- colorRampPalette(brewer.pal(8, "Set2"))(n_symbols)
  p <- ggplot(df, aes(x = symbol, y = value, fill = symbol)) +
    # 箱线图主体
    geom_boxplot(
      width = 0.25, 
      alpha = 0.8,           
      color = "gray30",      # 
      # outlier.color = "firebrick3",  #
      # outlier.size = 2.5,            # 
      # outlier.alpha = 0.7,           # 
      # outlier.shape = 21,            # 
      outliers = F,
      notch = FALSE,
      lwd = 0.6                      # 
    ) +
    stat_boxplot(
      geom = "errorbar", 
      width = 0.15, 
      color = "gray30",
      size = 0.8
    ) +
    scale_fill_manual(values = color_palette) +
    labs(
      x = xname,
      y = yname,
      fill = "sample"  
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "grey92", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "gray50", size = 0.5),
      axis.text = element_text(color = "gray30", size = 12),
      axis.title = element_text(color = "gray20", size = 14, face = "bold"),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = "gray80"),
      legend.position = "bottom",
      plot.margin = margin(15, 15, 15, 15)
    ) +
    coord_flip() +
    geom_jitter(
      width = 0.1,
      alpha = 0.1,
      size = 0.1,
      color = "gray40"
    )
  if (n_symbols > 10) {
    p <- p + guides(fill = "none")
  }
  return(p)
}
plot_length_group_boxplot <- function(df, sample_order = NULL) {
  breaks <- c(0, 100, 300, 10000, Inf)
  labels <- c("<100", "100-300", "300-10000", ">10000")
  df <- df %>%
    mutate(length_group = cut(length, 
                              breaks = breaks, 
                              labels = labels,
                              include.lowest = TRUE))
  if (!is.null(sample_order)) {
    df$group <- factor(df$group, levels = sample_order)
  }
  p <- ggplot(df, aes(x = length_group, y = TPM, fill = group)) +
    geom_boxplot(
      width = 0.7,
      alpha = 0.8,
      outlier.shape = NA,  # 
      position = position_dodge(width = 0.85),
      coef = 1.5,  # 
      show.legend = TRUE
    ) +
    scale_y_log10() +
    scale_fill_brewer(palette = "Set2") +
    labs(
      x = "Transcript Length (nt)",
      y = "TPM (log10 scale)",
      fill = "Sample Group"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.line = element_line(color = "black"),  
      axis.ticks = element_line(color = "black"),  #
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "bottom"
    )
  
  return(p)
}


#-----bar----------------------------
plot_tss_nrctss_log2 <- function(data, range = 2000, line_width = 0.5, tit = NULL) {
  plot_data <- data %>%
    mutate(
      adjusted_distance = ifelse(strand == "+", distanceToTSS, -distanceToTSS),
      count_value = ifelse(strand == "+", nr_ctss, -nr_ctss),
      strand = ifelse(strand == "+", "+", "-")  # 简化为+/-符号
    ) %>%
    filter(abs(adjusted_distance) <= range) %>%
    group_by(adjusted_distance, strand) %>%
    summarise(total_count = sum(count_value, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      log2_total_count = sign(total_count) * log2(abs(total_count) + 1)
    )
  vibrant_colors <- c("+" = "#FF0000", "-" = "#0099FF") 
  p <- ggplot(plot_data, aes(x = adjusted_distance)) +
    geom_segment(
      aes(xend = adjusted_distance, 
          y = 0, 
          yend = log2_total_count,
          color = strand),
      linewidth = line_width,
      lineend = "round"
    ) +
    scale_x_continuous(
      breaks = seq(-range, range, by = range/2),
      labels = function(x) format(abs(x), big.mark = ",")
    ) +
    scale_y_continuous(
      labels = abs
    ) +
    scale_color_manual(values = vibrant_colors) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      x = "Distance to TSS (bp)",
      y = "log2(CTSS Count Sum + 1)",
      title = tit, 
      color = NULL  
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      text = element_text(size = 12, color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_blank(),  
      legend.text = element_text(face = "bold", size = 12)
    )
  return(p)
}
plot_tss_tpm_log2 <- function(data, range = 2000, line_width = 0.5, tit = NULL) {
  plot_data <- data %>%
    mutate(
      adjusted_distance = ifelse(strand == "+", distanceToTSS, -distanceToTSS),
      tpm_value = ifelse(strand == "+", tpm.dominant_ctss, -tpm.dominant_ctss),
      strand = ifelse(strand == "+", "+", "-") 
    ) %>%
    filter(abs(adjusted_distance) <= range) %>%
    group_by(adjusted_distance, strand) %>%
    summarise(
      total_tpm = sum(tpm_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      log2_total_tpm = sign(total_tpm) * log2(abs(total_tpm) + 1)
    )
  vibrant_colors <- c("+" = "#FF0000", "-" = "#0099FF")
  p <- ggplot(plot_data, aes(x = adjusted_distance)) +
    geom_segment(
      aes(xend = adjusted_distance, 
          y = 0, 
          yend = log2_total_tpm,
          color = strand),
      linewidth = line_width,
      lineend = "round"
    ) +
    scale_x_continuous(
      breaks = seq(-range, range, by = range/2),
      labels = function(x) format(abs(x), big.mark = ",")
    ) +
    scale_y_continuous(
      labels = abs
    ) +
    scale_color_manual(values = vibrant_colors) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      x = "Distance to TSS (bp)",
      y = "log2(Total CTSS TPM + 1)",
      title = tit,  
      color = NULL 
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      text = element_text(size = 12, color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_blank(),  
      legend.text = element_text(face = "bold", size = 12)
    )
  return(p)
}
plot_tss_nrctss <- function(data, range = 2000, line_width = 0.5, tit = NULL) {
  plot_data <- data %>%
    mutate(
      adjusted_distance = ifelse(strand == "+", distanceToTSS, -distanceToTSS),
      count_value = ifelse(strand == "+", nr_ctss, -nr_ctss),
      strand = ifelse(strand == "+", "+", "-")  # 简化为+/-符号
    ) %>%
    filter(abs(adjusted_distance) <= range) %>%
    group_by(adjusted_distance, strand) %>%
    summarise(total_count = sum(count_value, na.rm = TRUE), .groups = "drop") 
  vibrant_colors <- c("+" = "#FF0000", "-" = "#0099FF") 
  p <- ggplot(plot_data, aes(x = adjusted_distance)) +
    geom_segment(
      aes(xend = adjusted_distance, 
          y = 0, 
          yend = total_count,
          color = strand),
      linewidth = line_width,
      lineend = "round"
    ) +
    scale_x_continuous(
      breaks = seq(-range, range, by = range/2),
      labels = function(x) format(abs(x), big.mark = ",")
    ) +
    scale_y_continuous(
      labels = abs
    ) +
    scale_color_manual(values = vibrant_colors) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      x = "Distance to TSS (bp)",
      y = "CTSS Count Sum",
      title = tit, 
      color = NULL  
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      text = element_text(size = 12, color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_blank(),  
      legend.text = element_text(face = "bold", size = 12)
    )
  return(p)
}
plot_tss_tpm <- function(data, range = 2000, line_width = 0.5, tit = NULL) {
  plot_data <- data %>%
    mutate(
      adjusted_distance = ifelse(strand == "+", distanceToTSS, -distanceToTSS),
      tpm_value = ifelse(strand == "+", tpm.dominant_ctss, -tpm.dominant_ctss),
      strand = ifelse(strand == "+", "+", "-") 
    ) %>%
    filter(abs(adjusted_distance) <= range) %>%
    group_by(adjusted_distance, strand) %>%
    summarise(
      total_tpm = sum(tpm_value, na.rm = TRUE),
      .groups = "drop"
    ) 
  vibrant_colors <- c("+" = "#FF0000", "-" = "#0099FF")
  p <- ggplot(plot_data, aes(x = adjusted_distance)) +
    geom_segment(
      aes(xend = adjusted_distance, 
          y = 0, 
          yend = total_tpm,
          color = strand),
      linewidth = line_width,
      lineend = "round"
    ) +
    scale_x_continuous(
      breaks = seq(-range, range, by = range/2),
      labels = function(x) format(abs(x), big.mark = ",")
    ) +
    scale_y_continuous(
      labels = abs
    ) +
    scale_color_manual(values = vibrant_colors) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      x = "Distance to TSS (bp)",
      y = "Total CTSS TPM",
      title = tit,  
      color = NULL  
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      text = element_text(size = 12, color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_blank(),  # 移除图例标题
      legend.text = element_text(face = "bold", size = 12)
    )
  return(p)
}
plot_tss_loc <- function(data,range = 2000, line_width = 0.5, tit = NULL) {
  plot_data <- data %>%
    mutate(
      adjusted_distance = ifelse(strand == "+", distanceToTSS, -distanceToTSS),
      count_value = ifelse(strand == "+", 1, -1)
    ) %>%
    filter(abs(adjusted_distance) <= range) %>%
    group_by(adjusted_distance) %>%
    summarise(
      total_count = sum(count_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(strand = ifelse(total_count > 0, "+", "-"))
  vibrant_colors <- c("+" = "#FF0000", "-" = "#0099FF")
  
  p <- ggplot(plot_data, aes(x = adjusted_distance)) +
    geom_segment(
      aes(xend = adjusted_distance, 
          y = 0, 
          yend = total_count,
          color = strand), 
      linewidth = line_width,
      lineend = "round"
    ) +
    scale_x_continuous(
      breaks = seq(-range, range, by = range/2),
      labels = function(x) format(abs(x), big.mark = ",")
    ) +
    scale_y_continuous(labels = abs) +
    scale_color_manual(values = vibrant_colors) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      x = "Distance to TSS (bp)",
      y = "Number of Sites",
      title = tit,
      color = "Strand"
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      text = element_text(size = 12, color = "black"),
      axis.text = element_text(color = "black"),
      legend.title = element_blank(), 
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) 
  return(p)
}

plot_tss_loc_percent <- function(data, range = 2000, line_width = 0.5, tit = NULL) {
  plot_data <- data %>%
    mutate(
      adjusted_distance = ifelse(strand == "+", distanceToTSS, -distanceToTSS),
      count_value = ifelse(strand == "+", 1, -1)
    ) %>%
    filter(abs(adjusted_distance) <= range) %>%
    group_by(adjusted_distance) %>%
    summarise(
      total_count = sum(count_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(strand = ifelse(total_count > 0, "+", "-"))
  max_abs_pos <- plot_data %>% 
    filter(strand == "+") %>% 
    summarise(max_val = max(total_count, na.rm = TRUE)) %>% 
    pull(max_val) %>% 
    abs()
  
  max_abs_neg <- plot_data %>% 
    filter(strand == "-") %>% 
    summarise(max_val = min(total_count, na.rm = TRUE)) %>%  # 
    pull(max_val) %>% 
    abs()
  if (is.infinite(max_abs_pos)) max_abs_pos <- 1
  if (is.infinite(max_abs_neg)) max_abs_neg <- 1
  plot_data <- plot_data %>%
    mutate(
      percent_value = ifelse(
        strand == "+", 
        total_count / max_abs_pos * 100,
        total_count / max_abs_neg * 100
      )
    )
  
  vibrant_colors <- c("+" = "#FF0000", "-" = "#0099FF")
  
  p <- ggplot(plot_data, aes(x = adjusted_distance)) +
    geom_segment(
      aes(xend = adjusted_distance, 
          y = 0, 
          yend = percent_value,
          color = strand),
      linewidth = line_width,
      lineend = "round"
    ) +
    scale_x_continuous(
      breaks = seq(-range, range, by = range/2),
      labels = function(x) format(abs(x), big.mark = ",")
    ) +
    scale_y_continuous(
      labels = function(y) paste0(abs(y), "%")  # 
    ) +
    scale_color_manual(values = vibrant_colors) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      x = "Distance to TSS (bp)",
      y = "Percentage of number of Sites(%)",  # 修改y轴标签
      title = tit,
      color = "Strand"
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      text = element_text(size = 12, color = "black"),
      axis.text = element_text(color = "black"),
      legend.title = element_blank(), 
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) 
  return(p)
}
plot_tss_nrctss_percent <- function(data, range = 2000, line_width = 0.5, tit = NULL) {
  plot_data <- data %>%
    mutate(
      adjusted_distance = ifelse(strand == "+", distanceToTSS, -distanceToTSS),
      count_value = ifelse(strand == "+", nr_ctss, -nr_ctss),
      strand = ifelse(strand == "+", "+", "-")
    ) %>%
    filter(abs(adjusted_distance) <= range)
  position_counts <- plot_data %>%
    group_by(adjusted_distance, strand) %>%
    summarise(
      total_count = sum(count_value, na.rm = TRUE),
      .groups = "drop"
    )
  max_abs_pos <- position_counts %>%
    filter(strand == "+") %>%
    summarise(max_val = max(total_count, na.rm = TRUE)) %>%
    pull(max_val) %>%
    abs()
  
  max_abs_neg <- position_counts %>%
    filter(strand == "-") %>%
    summarise(max_val = min(total_count, na.rm = TRUE)) %>%  # 
    pull(max_val) %>%
    abs()
  if (is.infinite(max_abs_pos)) max_abs_pos <- 1
  if (is.infinite(max_abs_neg)) max_abs_neg <- 1
  plot_data <- position_counts %>%
    mutate(
      percent_value = ifelse(
        strand == "+",
        total_count / max_abs_pos * 100,
        total_count / max_abs_neg * 100
      )
    )
  
  vibrant_colors <- c("+" = "#FF0000", "-" = "#0099FF") 
  
  p <- ggplot(plot_data, aes(x = adjusted_distance)) +
    geom_segment(
      aes(xend = adjusted_distance, 
          y = 0, 
          yend = percent_value,
          color = strand),
      linewidth = line_width,
      lineend = "round"
    ) +
    scale_x_continuous(
      breaks = seq(-range, range, by = range/2),
      labels = function(x) format(abs(x), big.mark = ",")
    ) +
    scale_y_continuous(
      labels = function(y) paste0(abs(y), "%")
    ) +
    scale_color_manual(values = vibrant_colors) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      x = "Distance to TSS (bp)",
      y = "Percentage of CTSS count sum(%)",  
      title = tit, 
      color = NULL  
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      text = element_text(size = 12, color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_blank(),
      legend.text = element_text(face = "bold", size = 12)
    )
  
  return(p)
}
plot_tss_tpm_percent <- function(data, range = 2000, line_width = 0.5, tit = NULL) {
  plot_data <- data %>%
    mutate(
      adjusted_distance = ifelse(strand == "+", distanceToTSS, -distanceToTSS),
      tpm_value = ifelse(strand == "+", tpm.dominant_ctss, -tpm.dominant_ctss),
      strand = ifelse(strand == "+", "+", "-") 
    ) %>%
    filter(abs(adjusted_distance) <= range) %>%
    group_by(adjusted_distance, strand) %>%
    summarise(
      total_tpm = sum(tpm_value, na.rm = TRUE),
      .groups = "drop"
    )
  
  max_abs_pos <- plot_data %>%
    filter(strand == "+") %>%
    summarise(max_val = max(total_tpm, na.rm = TRUE)) %>%
    pull(max_val) %>%
    abs()
  
  max_abs_neg <- plot_data %>%
    filter(strand == "-") %>%
    summarise(max_val = min(total_tpm, na.rm = TRUE)) %>%  # 
    pull(max_val) %>%
    abs()
  
  if (is.infinite(max_abs_pos)) max_abs_pos <- 1
  if (is.infinite(max_abs_neg)) max_abs_neg <- 1
  
  plot_data <- plot_data %>%
    mutate(
      percent_value = ifelse(
        strand == "+",
        total_tpm / max_abs_pos * 100,
        total_tpm / max_abs_neg * 100
      )
    )
  
  vibrant_colors <- c("+" = "#FF0000", "-" = "#0099FF")
  
  p <- ggplot(plot_data, aes(x = adjusted_distance)) +
    geom_segment(
      aes(xend = adjusted_distance, 
          y = 0, 
          yend = percent_value,
          color = strand),
      linewidth = line_width,
      lineend = "round"
    ) +
    scale_x_continuous(
      breaks = seq(-range, range, by = range/2),
      labels = function(x) format(abs(x), big.mark = ",")
    ) +
    scale_y_continuous(
      labels = function(y) paste0(abs(y), "%")  # 
    ) + 
    scale_color_manual(values = vibrant_colors) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      x = "Distance to TSS (bp)",
      y = "Percentage of Peak TPM (%)",  #
      title = tit, 
      color = NULL  
    ) +
    theme_bw() +
    theme(
      legend.position = "top",  
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      text = element_text(size = 12, color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_blank(),
      legend.text = element_text(face = "bold", size = 12)
    )
  return(p)
}

plot_tss_loc_squ_percent <- function(data, range = 2000, line_width = 0.5, tit = NULL) {
  plot_data <- data %>%
    mutate(
      adjusted_distance = ifelse(strand == "+", distanceToTSS, -distanceToTSS),
      count_value = ifelse(strand == "+", 1, -1)
    ) %>%
    filter(abs(adjusted_distance) <= range) %>%
    group_by(adjusted_distance) %>%
    summarise(
      total_count = sum(count_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(strand = ifelse(total_count > 0, "+", "-"))
  
  max_abs_pos <- plot_data %>% 
    filter(strand == "+") %>% 
    summarise(max_val = max(total_count, na.rm = TRUE)) %>% 
    pull(max_val) %>% 
    abs()
  
  max_abs_neg <- plot_data %>% 
    filter(strand == "-") %>% 
    summarise(max_val = min(total_count, na.rm = TRUE)) %>%
    pull(max_val) %>% 
    abs()
  
  if (is.infinite(max_abs_pos)) max_abs_pos <- 1
  if (is.infinite(max_abs_neg)) max_abs_neg <- 1
  
  plot_data <- plot_data %>%
    mutate(
      sqrt_count = ifelse(
        strand == "+", 
        sqrt(pmax(0, total_count)),  #
        -sqrt(pmax(0, -total_count))  #
      ),
      percent_value = ifelse(
        strand == "+",
        sqrt_count / sqrt(max_abs_pos) * 100,
        sqrt_count / sqrt(max_abs_neg) * 100
      )
    )
  
  vibrant_colors <- c("+" = "#FF0000", "-" = "#0099FF")
  
  p <- ggplot(plot_data, aes(x = adjusted_distance)) +
    geom_segment(
      aes(xend = adjusted_distance, 
          y = 0, 
          yend = percent_value,
          color = strand),
      linewidth = line_width,
      lineend = "round"
    ) +
    scale_x_continuous(
      breaks = seq(-range, range, by = range/2),
      labels = function(x) format(abs(x), big.mark = ",")
    ) +
    scale_y_continuous(
      labels = function(y) paste0(abs(y), "%")  # 
    ) +
    scale_color_manual(values = vibrant_colors) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      x = "Distance to TSS (bp)",
      y = "Percentage of sqrt-scale number of Sites(%)",  #
      title = tit,
      color = "Strand"
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      text = element_text(size = 12, color = "black"),
      axis.text = element_text(color = "black"),
      legend.title = element_blank(), 
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) 
  return(p)
}
plot_tss_nrctss_squ_percent <- function(data, range = 2000, line_width = 0.5, tit = NULL) {
  plot_data <- data %>%
    mutate(
      adjusted_distance = ifelse(strand == "+", distanceToTSS, -distanceToTSS),
      count_value = ifelse(strand == "+", nr_ctss, -nr_ctss),
      strand = ifelse(strand == "+", "+", "-")
    ) %>%
    filter(abs(adjusted_distance) <= range)
  position_counts <- plot_data %>%
    group_by(adjusted_distance, strand) %>%
    summarise(
      total_count = sum(count_value, na.rm = TRUE),
      .groups = "drop"
    )
  max_abs_pos <- position_counts %>%
    filter(strand == "+") %>%
    summarise(max_val = max(total_count, na.rm = TRUE)) %>%
    pull(max_val) %>%
    abs()
  
  max_abs_neg <- position_counts %>%
    filter(strand == "-") %>%
    summarise(max_val = min(total_count, na.rm = TRUE)) %>%
    pull(max_val) %>%
    abs()
  
  if (is.infinite(max_abs_pos)) max_abs_pos <- 1
  if (is.infinite(max_abs_neg)) max_abs_neg <- 1
  
  plot_data <- position_counts %>%
    mutate(
      sqrt_count = ifelse(
        strand == "+", 
        sqrt(pmax(0, total_count)),  # 
        -sqrt(pmax(0, -total_count))  # 
      ),
      percent_value = ifelse(
        strand == "+",
        sqrt_count / sqrt(max_abs_pos) * 100,
        sqrt_count / sqrt(max_abs_neg) * 100
      )
    )
  
  vibrant_colors <- c("+" = "#FF0000", "-" = "#0099FF") 
  p <- ggplot(plot_data, aes(x = adjusted_distance)) +
    geom_segment(
      aes(xend = adjusted_distance, 
          y = 0, 
          yend = percent_value,
          color = strand),
      linewidth = line_width,
      lineend = "round"
    ) +
    scale_x_continuous(
      breaks = seq(-range, range, by = range/2),
      labels = function(x) format(abs(x), big.mark = ",")
    ) +
    scale_y_continuous(
      labels = function(y) paste0(abs(y), "%")  # 
    ) + 
    scale_color_manual(values = vibrant_colors) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      x = "Distance to TSS (bp)",
      y = "Percentage of sqrt-scale CTSS count sum(%)",  #
      title = tit, 
      color = NULL  
    ) +
    theme_bw() +
    theme(
      legend.position = "top",  
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      text = element_text(size = 12, color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_blank(),
      legend.text = element_text(face = "bold", size = 12)
    )
  return(p)
}
plot_tss_tpm_squ_percent <- function(data, range = 2000, line_width = 0.5, tit = NULL) {
  plot_data <- data %>%
    mutate(
      adjusted_distance = ifelse(strand == "+", distanceToTSS, -distanceToTSS),
      tpm_value = ifelse(strand == "+", tpm.dominant_ctss, -tpm.dominant_ctss),
      strand = ifelse(strand == "+", "+", "-")
    ) %>%
    filter(abs(adjusted_distance) <= range)
  
  position_tpm <- plot_data %>%
    group_by(adjusted_distance, strand) %>%
    summarise(
      total_tpm = sum(tpm_value, na.rm = TRUE),
      .groups = "drop"
    )
  
  max_abs_pos <- position_tpm %>%
    filter(strand == "+") %>%
    summarise(max_val = max(total_tpm, na.rm = TRUE)) %>%
    pull(max_val) %>%
    abs()
  
  max_abs_neg <- position_tpm %>%
    filter(strand == "-") %>%
    summarise(max_val = min(total_tpm, na.rm = TRUE)) %>%
    pull(max_val) %>%
    abs()
  
  if (is.infinite(max_abs_pos)) max_abs_pos <- 1
  if (is.infinite(max_abs_neg)) max_abs_neg <- 1
  
  plot_data <- position_tpm %>%
    mutate(
      sqrt_tpm = ifelse(
        strand == "+", 
        sqrt(pmax(0, total_tpm)),  # 
        -sqrt(pmax(0, -total_tpm))  # 
      ),
      percent_value = ifelse(
        strand == "+",
        sqrt_tpm / sqrt(max_abs_pos) * 100,
        sqrt_tpm / sqrt(max_abs_neg) * 100
      )
    )
  
  vibrant_colors <- c("+" = "#FF0000", "-" = "#0099FF")
  
  p <- ggplot(plot_data, aes(x = adjusted_distance)) +
    geom_segment(
      aes(xend = adjusted_distance, 
          y = 0, 
          yend = percent_value,
          color = strand),
      linewidth = line_width,
      lineend = "round"
    ) +
    scale_x_continuous(
      breaks = seq(-range, range, by = range/2),
      labels = function(x) format(abs(x), big.mark = ",")
    ) +
    scale_y_continuous(
      labels = function(y) paste0(abs(y), "%")  # 
    ) + 
    scale_color_manual(values = vibrant_colors) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      x = "Distance to TSS (bp)",
      y = "Percentage of sqrt-scale Peak TPM(%)",  # 
      title = tit, 
      color = NULL  
    ) +
    theme_bw() +
    theme(
      legend.position = "top",  
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      text = element_text(size = 12, color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_blank(),
      legend.text = element_text(face = "bold", size = 12)
    )
  return(p)
}



