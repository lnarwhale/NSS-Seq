#---setting----
work_dir <- "/project/code/02_r_analysis/analysis/"
bam_dir <- "/project/code/01_shell_pipeline/res/result/5_map/"
fuc_dir <- "/project/code/02_r_analysis/functions/"
res_dir <- "/project/code/02_r_analysis/analysis/res/"

#---library--------
library(AnnotationDbi)
library(bedtoolsr)
library(CAGEr)
library(ChIPseeker)
library(circlize)
library(clusterProfiler)
library(ComplexHeatmap)
library(DESeq2)
library(dplyr)
library(edgeR)
library(factoextra)
library(ggcorrplot)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggseqlogo)
library(ggvenn)
library(grid)
library(gridExtra)
library(hdf5r)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(seqPattern)
library(stats)
library(stringi)
library(stringr)
library(viridis)

#---pre-----
setwd(work_dir)

library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
seqnames(BSgenome.Hsapiens.UCSC.hg38) <- gsub("chr","",seqnames(BSgenome.Hsapiens.UCSC.hg38))
BSname <- "BSgenome.Hsapiens.UCSC.hg38"
BS <- BSgenome.Hsapiens.UCSC.hg38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
org <- "org.Hs.eg.db"
seqlevels(txdb) <- gsub("chr","",seqlevels(txdb))

source(paste0(fuc_dir,"/cager_function_set.R"))
source(paste0(fuc_dir,"/paper_plot.R"))
source(paste0(fuc_dir,"/pipeline_function.R"))
source(paste0(fuc_dir,"/universal_function.R"))

#---running----
cager_get_domTSS_order(bam_dir)
cager_pre <- cager_get_LICAGE_exp_pre(bam_dir,c(1,1,2,2),c("test_IFN_0h","test_IFN_1h"))
pdf(paste0(res_dir,"/1_rev_cagers_mer.pdf"))
plotReverseCumulatives(cager_pre,fitInRange = c(5,1000),onePlot = T)
dev.off()

cager_merge_af <- cager_get_LICAGE_exp_af(cager_pre,1.44,5*10^5,20)
cager_merge_agg <- cager_aggregate(cager_merge_af,100)
cager_merge_matrix <- cager_matrix(cager_merge_agg)
cager_merge_dom <- cager_get_domTSS(cager_merge_agg)
cager_merge_dom_df <- lapply(cager_merge_dom, function(x){
  as.data.frame(x)
})

cager_pre_nomer <- cager_get_LICAGE_exp_pre_notmer_tsscor(bam_dir)
pdf(paste0(res_dir,"/2_rev_cagers_nomer.pdf"))
plotReverseCumulatives(cager_pre_nomer,fitInRange = c(5,1000),onePlot = T)
dev.off()
cager_nomerge_af <- cager_get_LICAGE_exp_af(cager_pre_nomer,1.27,10^5,20)
cager_nomerge_agg <- cager_aggregate(cager_nomerge_af,100)
cager_nomerge_matrix <- cager_matrix(cager_nomerge_agg)
cager_nomerge_dom <- cager_get_domTSS(cager_nomerge_agg)
cager_nomerge_dom_df <- lapply(cager_nomerge_dom, function(x){
  as.data.frame(x)
})


pdf(paste0(res_dir,"/3_corr_sample.pdf"),width = 30,height = 30)
plot_all_correlations(cager_nomerge_agg,c("test_IFN_0h_1","test_IFN_0h_2","test_IFN_6h_1","test_IFN_6h_2"))
dev.off()

#------peak_bp------
plot_pap_din <- pap_dflis_2_freq_din_orderCT(cager_merge_dom_df,
                                             c("test_IFN_0h","test_IFN_1h"),
                                             c("CA","TA","TG","CG","CC","TC","GG","AG","CT","TT","GT","AA","AT","GA","GC","AC"))
ggsave(paste0(res_dir,'/4_pap_din_0h.pdf'),plot_pap_din$picture)


plot_pca <- pap_pca_group(cager_nomerge_matrix,group =c("test_IFN_0h_1","test_IFN_0h_2","test_IFN_6h_1","test_IFN_6h_2"),
                          order = c("test_IFN_0h_1","test_IFN_0h_2","test_IFN_6h_1","test_IFN_6h_2"))
ggsave(paste0(res_dir,"/5_pca_plot.pdf"),plot_pca,width = 10,height = 10)

plot_anno_peak <- pap_peakanno_lis(cager_merge_dom)
ggsave(paste0(res_dir,"/6_plot_anno_peak.pdf"),plot_anno_peak,width = 10,height = 5)

pap_actg_line_from_dom(cager_merge_dom,100,100,"7_atcg",res_dir)

pap_heatmap_ud_10_from_dom(cager_merge_dom,res_dir,"8_heatmap_ud_10")

pap_logo_from_dom(cager_merge_dom,5,5,"9_logo_ud10",res_dir)

#------clu-heatmap-----
cager_matrix_clu <- cluster_re(cager_merge_matrix,4)
cager_matrix_cctpm <- lapply(cager_matrix_clu, function(item){
  df <- data.frame(location = names(item))
  merge(df,cager_merge_matrix,by = "location")
})
for (i in 1:length(cager_matrix_cctpm)) {
  cager_matrix_cctpm[[i]] <- dplyr::select(cager_matrix_cctpm[[i]],"location",
                                           "test_IFN_0h","test_IFN_1h")
}
cager_matrix_cctpm_rbind <- rbind(cager_matrix_cctpm[[2]],cager_matrix_cctpm[[4]],
                                  cager_matrix_cctpm[[3]],cager_matrix_cctpm[[1]])
cager_matrix_cctpm_hmap <- heat_map_cut_noname(cager_matrix_cctpm_rbind)
pdf(paste0("res/10_heatmap_cut_",i,".pdf"))
cager_matrix_cctpm_hmap
dev.off()

cager_seq <- lapply(cager_merge_dom, function(x){
  cluter2seq(100,100,x)
})

plot_motifoccur_lis(cager_merge_dom,c("green3","red3"),100,20,70,3,"11_TATA_motif",TATA,"TATA",c(1,2),res_dir)




