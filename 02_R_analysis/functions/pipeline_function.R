#------pre end----
IFN_bulk_diff_pip <- function(count_dir,count_name,select_col,replace_col,
                          col_group,pv,fc,
                          res_dir,pre_name,top_num){
  #-----read in----
  #count_dir:count所在文件夹 ;count_name:count的名字 ;
  #select_col:选取count的列名 (control在前);replace_col:count的列名替代
  #----diff analysis----
  #col_group:差异分析的group; pv ; fc
  #res_dir:输出图片和文件的文件夹; pre_name:所有图片前置名称
  #----go analysis----
  #group_downup:down,up的bar的group; top_num:bar_go的前置go数目
  #----examples---
  # tmp <- IFN_bulk_diff_pip(count_dir = count_dir,count_name = "total_gene_count",
  #                          select_col = c("NSS_HCT_CA_0H_1_241115","NSS_HCT_CA_0H_2_241115","NSS_HCT_CA_2H_1_241115","NSS_HCT_CA_2H_2_241115"),
  #                          replace_col = c("nascent_CA_0h_1","nascent_CA_0h_2","nascent_CA_1h_1","nascent_CA_1h_2"),
  #                          col_group = c("nascent_CA_0h","nascent_CA_0h","nascent_CA_1h","nascent_CA_1h"),
  #                          pv=0.01,fc=1,res_dir=res_dir,pre_name  = "1_CA1h_CA0h",top_num=10)
  #-----read in-----
  raw_count <- read.csv(paste0(count_dir,"/",count_name,".txt"))
  raw_count_ens <- de_1(raw_count,1)
  raw_count_all <- ens_2_sym(raw_count_ens,1)
  raw_count_sym <- raw_count_all
  raw_count_sym$ENSEMBL <- NULL
  raw_count_sym$ENTREZID <- NULL
  raw_count_sym_se <- dplyr::select(raw_count_sym,"SYMBOL",select_col)
  colnames(raw_count_sym_se) <- c("SYMBOL",replace_col)
  #-----different analysis-------
  diff_raw_count_sym_se <- plot_deseq2_multi_compare(raw_count_sym_se,col_group,pv,fc)
  ggsave(paste0(res_dir,"/",pre_name,"_1_vac_mark.png"),diff_raw_count_sym_se$valco_plot_mark)
  ggsave(paste0(res_dir,"/",pre_name,"_1_vac_mark.pdf"),diff_raw_count_sym_se$valco_plot_mark)
  ggsave(paste0(res_dir,"/",pre_name,"_2_hmap.png"),diff_raw_count_sym_se$heat_map)
  ggsave(paste0(res_dir,"/",pre_name,"_2_hmap.pdf"),diff_raw_count_sym_se$heat_map)
  write.csv(diff_raw_count_sym_se$res_df$change,paste0(res_dir,"/",pre_name,"_25_gene.csv"))
  #-----GO analysis------
  geneup_diff_re <- diff_raw_count_sym_se$res_df$up
  genedown_diff_re <- diff_raw_count_sym_se$res_df$down
  GO_geneup_diff_re <- enrichGO(gene=rownames(geneup_diff_re),
                                           OrgDb = org,
                                           keyType = "SYMBOL",
                                           ont="ALL",
                                           pAdjustMethod = "BH",
                                           pvalueCutoff = 1,
                                           qvalueCutoff = 1,
                                           readable = T)
  GO_genedown_diff_re <- enrichGO(gene=rownames(genedown_diff_re),
                                             OrgDb = org,
                                             keyType = "SYMBOL",
                                             ont="ALL",
                                             pAdjustMethod = "BH",
                                             pvalueCutoff = 1,
                                             qvalueCutoff = 1,
                                             readable = T)
  GO_geneup_diff_re_df <- as.data.frame(GO_geneup_diff_re)
  GO_genedown_diff_re_df <- as.data.frame(GO_genedown_diff_re)
  #----up------
  GO_geneup_diff_re_df$log10pdj <- -log10(GO_geneup_diff_re_df$p.adjust)
  GO_geneup_diff_re_df_sorted <- GO_geneup_diff_re_df[order(GO_geneup_diff_re_df$p.adjust),]
  GO_geneup_diff_re_df_sorted$Description <- factor(GO_geneup_diff_re_df_sorted$Description, 
                                                    levels = GO_geneup_diff_re_df_sorted$Description)
  up_pattern <- c("type II interferon-mediated signaling pathway", "response to type II interferon",
                     "cellular response to type II interferon","regulation of response to type II interferon",
                     "receptor signaling pathway via JAK-STAT")
  p_GO_pattern_geneup_diff_re_df <- ggplot(data=GO_geneup_diff_re_df_sorted, aes(x=Description, y=log10pdj)) +
    geom_point(size = 0.01) + 
    geom_point(data = subset(GO_geneup_diff_re_df_sorted,Description %in% up_pattern),
               color = "red", size = 0.5) +
    geom_text(data = subset(GO_geneup_diff_re_df_sorted,Description %in% up_pattern),
              aes(label = Description),
              hjust = -0.02, vjust = -0.5, size = 2, color = "black")+
    theme_test() + 
      xlab("GO term") + 
      ylab("-log10 p.adjust") +
      labs(title = "Enriched GO Terms")+
      theme(
        axis.text.x = element_blank(),  
        axis.text.y = element_text(color = "black",size=5),
        plot.title = element_text(hjust = 0.5)
      ) +
      scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, by = 2),
                       expand = expansion(mult = c(0, 0.05))) +
      geom_hline(yintercept = 1.3, linetype = "dashed", color = "red", size = 0.2)
  ggsave(paste0(res_dir,"/",pre_name,"_3_p_GO_up_pattern.png"),
         p_GO_pattern_geneup_diff_re_df,width = 5,height=6)
  ggsave(paste0(res_dir,"/",pre_name,"_3_p_GO_up_pattern.pdf"),
         p_GO_pattern_geneup_diff_re_df,width = 5,height=6)
  #---down-----
  GO_genedown_diff_re_df$log10pdj <- -log10(GO_genedown_diff_re_df$p.adjust)
  GO_genedown_diff_re_df_sorted <- GO_genedown_diff_re_df[order(GO_genedown_diff_re_df$p.adjust),]
  GO_genedown_diff_re_df_sorted$Description <- factor(GO_genedown_diff_re_df_sorted$Description, 
                                                    levels = GO_genedown_diff_re_df_sorted$Description)
  down_pattern <- c("type II interferon-mediated signaling pathway", "response to type II interferon",
                  "cellular response to type II interferon","regulation of response to type II interferon",
                  "receptor signaling pathway via JAK-STAT")
  p_GO_pattern_genedown_diff_re_df <- ggplot(data=GO_genedown_diff_re_df_sorted, aes(x=Description, y=log10pdj)) +
    geom_point(size = 0.01) + 
    geom_point(data = subset(GO_genedown_diff_re_df_sorted,Description %in% down_pattern),
               color = "red", size = 0.5) +
    geom_text(data = subset(GO_genedown_diff_re_df_sorted,Description %in% down_pattern),
              aes(label = Description),
              hjust = -0.02, vjust = -0.5, size = 2, color = "black")+
    theme_test() + 
    xlab("GO term") + 
    ylab("-log10 p.adjust") +
    labs(title = "Enriched GO Terms")+
    theme(
      axis.text.x = element_blank(),  
      axis.text.y = element_text(color = "black",size=5),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, by = 2),
                       expand = expansion(mult = c(0, 0.05))) +
    geom_hline(yintercept = 1.3, linetype = "dashed", color = "red", size = 0.2)
  ggsave(paste0(res_dir,"/",pre_name,"_4_p_GO_down_pattern.png"),
         p_GO_pattern_genedown_diff_re_df,width = 5,height=6)
  ggsave(paste0(res_dir,"/",pre_name,"_4_p_GO_down_pattern.pdf"),
         p_GO_pattern_genedown_diff_re_df,width = 5,height=6)
  
  write.csv(GO_geneup_diff_re_df_sorted,paste0(res_dir,"/",pre_name,"_5_geneup_diff_res.csv"))
  write.csv(GO_genedown_diff_re_df_sorted,paste0(res_dir,"/",pre_name,"_5_genedown_diff_res.csv"))
  
  go_df <- list(GO_genedown_diff_re_df_sorted,GO_geneup_diff_re_df_sorted)
  names(go_df) <- unique(col_group)
  bar_go_df <- plot_GO_bar(go_df,c(1,2),top_num)
  ggsave(paste0(res_dir,"/",pre_name,"_5_bar_go_df.png"),bar_go_df)
  ggsave(paste0(res_dir,"/",pre_name,"_5_bar_go_df.pdf"),bar_go_df)

  bar_GO_geneup <- barplot(GO_geneup_diff_re)
  ggsave(paste0(res_dir,"/",pre_name,"_6_bar_GO_geneup.png"),bar_GO_geneup)
  ggsave(paste0(res_dir,"/",pre_name,"_6_bar_GO_geneup.pdf"),bar_GO_geneup)
  # dot_GO_geneup <- dotplot(GO_geneup_diff_re)
  # ggsave(paste0(res_dir,"/",pre_name,"_6_dot_GO_geneup.png"),dot_GO_geneup)
  # ggsave(paste0(res_dir,"/",pre_name,"_6_dot_GO_geneup.pdf"),dot_GO_geneup)

  bar_GO_genedown <- barplot(GO_genedown_diff_re)
  ggsave(paste0(res_dir,"/",pre_name,"_7_bar_GO_genedown.png"),bar_GO_genedown)
  ggsave(paste0(res_dir,"/",pre_name,"_7_bar_GO_genedown.pdf"),bar_GO_genedown)
  # dot_GO_genedown <- dotplot(GO_genedown_diff_re)
  # ggsave(paste0(res_dir,"/",pre_name,"_7_dot_GO_genedown.png"),dot_GO_genedown)
  # ggsave(paste0(res_dir,"/",pre_name,"_7_dot_GO_genedown.pdf"),dot_GO_genedown)
  the_res <- list(diff_raw_count_sym_se,GO_geneup_diff_re_df_sorted,GO_genedown_diff_re_df_sorted)
  names(the_res) <- c("diff_res","GO_up_df","GO_down_df")
  return(the_res)
}

IFN_gene_lis_to_GO <- function(gene_lis,res_dir,pre_name,top_num){
  #-----read in----
  #gene_lis <- list("a,b,c","d,e,f,g"...),names(gene_lis) <- name
  #-----GO analysis------
  group_name <- names(gene_lis)
  GO_re_df_lis <- list()
  for (i in 1:length(gene_lis)) {
    write.csv(gene_lis[[i]],paste0(res_dir,"/",pre_name,"_1_gene_",group_name[i],".csv"))
    GO_re <- enrichGO(gene=gene_lis[[i]],
                      OrgDb = org,
                      keyType = "SYMBOL",
                      ont="ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,
                      readable = T)
    GO_re_df <- as.data.frame(GO_re)
    GO_re_df$log10pdj <- -log10(GO_re_df$p.adjust)
    GO_re_df_sorted <- GO_re_df[order(GO_re_df$p.adjust),]
    GO_re_df_sorted$Description <- factor(GO_re_df_sorted$Description, 
                                          levels = GO_re_df_sorted$Description)
    write.csv(GO_re_df_sorted,paste0(res_dir,"/",pre_name,"_2_GO_",group_name[i],".csv"))
    up_pattern <- c("type II interferon-mediated signaling pathway", "response to type II interferon",
                    "cellular response to type II interferon","regulation of response to type II interferon",
                    "receptor signaling pathway via JAK-STAT")
    p_GO_pattern <- ggplot(data=GO_re_df_sorted, aes(x=Description, y=log10pdj)) +
      geom_point(size = 0.01) + 
      geom_point(data = subset(GO_re_df_sorted,Description %in% up_pattern),
                 color = "red", size = 0.5) +
      geom_text(data = subset(GO_re_df_sorted,Description %in% up_pattern),
                aes(label = Description),
                hjust = -0.02, vjust = -0.5, size = 2, color = "black")+
      theme_test() + 
      xlab("GO term") + 
      ylab("-log10 p.adjust") +
      labs(title = "Enriched GO Terms")+
      theme(
        axis.text.x = element_blank(),  
        axis.text.y = element_text(color = "black",size=5),
        plot.title = element_text(hjust = 0.5)
      ) +
      scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, by = 2),
                         expand = expansion(mult = c(0, 0.05))) +
      geom_hline(yintercept = 1.3, linetype = "dashed", color = "red", size = 0.2)
    ggsave(paste0(res_dir,"/",pre_name,"_3_p_GO_pattern_",group_name[i],".png"),
           p_GO_pattern,width = 5,height=6)
    ggsave(paste0(res_dir,"/",pre_name,"_3_p_GO_pattern_",group_name[i],".pdf"),
           p_GO_pattern,width = 5,height=6)
    
    bar_GO_res <- barplot(GO_re)
    ggsave(paste0(res_dir,"/",pre_name,"_4_bar_GO_",group_name[i],".png"),bar_GO_res)
    ggsave(paste0(res_dir,"/",pre_name,"_4_bar_GO_",group_name[i],".pdf"),bar_GO_res)
    dot_GO_res <- dotplot(GO_re)
    ggsave(paste0(res_dir,"/",pre_name,"_5_dot_GO_",group_name[i],".png"),dot_GO_res)
    ggsave(paste0(res_dir,"/",pre_name,"_5_dot_GO_",group_name[i],".pdf"),dot_GO_res)
    GO_re_df_lis[[i]] <- GO_re_df_sorted
  }
  names(GO_re_df_lis) <- group_name
  bar_go_df <- plot_GO_bar(GO_re_df_lis,c(1,2),top_num)
  #bar_go_df <- plot_GO_bar(GO_re_df_lis,c(1,2,3),top_num)
  ggsave(paste0(res_dir,"/",pre_name,"_6_bar_go_df.png"),bar_go_df)
  ggsave(paste0(res_dir,"/",pre_name,"_6_bar_go_df.pdf"),bar_go_df)
  the_res <- GO_re_df_lis
  return(the_res)
}
