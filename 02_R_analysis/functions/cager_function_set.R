meme_motif_in <- function(){
  #----human------
  TATA <- from_meme_get_ppm("/Share/user/shenlm/reference/motif_db/HOCOMOCOv11_full_HUMAN/TBP.meme")
  YY1 <- from_meme_get_ppm("/data1/shenluemou/reference/meme/YY1_hs.meme")
  ELK1 <- from_meme_get_ppm("/data1/shenluemou/reference/meme/ELK1_hs.meme")
  GABPA <- from_meme_get_ppm("/data1/shenluemou/reference/meme/GABPA_hs.meme")
  E2F1 <- from_meme_get_ppm("/data1/shenluemou/reference/meme/E2F_hs.meme")
  ARNT <- from_meme_get_ppm("/data1/shenluemou/reference/meme/ARNT_hs.meme")
  SP1 <- from_meme_get_ppm("/data1/shenluemou/reference/meme/SP1_hs.meme")
  NFYa <- from_meme_get_ppm("/data1/shenluemou/reference/meme/NFYA_hs.meme")
  NFYb <- from_meme_get_ppm("/data1/shenluemou/reference/meme/NFYB_hs.meme")
  NFYc <- from_meme_get_ppm("/data1/shenluemou/reference/meme/NFYC_hs.meme")
  NRF1 <- from_meme_get_ppm("/data1/shenluemou/reference/meme/NRF1_hs.meme")
  CREB <- from_meme_get_ppm("/data1/shenluemou/reference/meme/CREB1_hs.meme")
  ATF1 <- from_meme_get_ppm("/data1/shenluemou/reference/meme/ATF1_hs.meme")
}
# cager_ref_in <- function(){
#   #library
#   if (species == "human") {
#     library(BSgenome.Hsapiens.UCSC.hg38)
#     library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#     library(org.Hs.eg.db)
#     seqnames(BSgenome.Hsapiens.UCSC.hg38) <- gsub("chr","",seqnames(BSgenome.Hsapiens.UCSC.hg38))
#     BSname <- "BSgenome.Hsapiens.UCSC.hg38"
#     BS <- BSgenome.Hsapiens.UCSC.hg38
#     txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#     org <- "org.Hs.eg.db"
#     seqlevels(txdb) <- gsub("chr","",seqlevels(txdb))
#   }else if (species == "mouse") {
#     library(BSgenome.Mmusculus.UCSC.mm10)
#     library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#     library(org.Mm.eg.db)
#     seqnames(BSgenome.Mmusculus.UCSC.mm10) <- gsub("chr","",seqnames(BSgenome.Mmusculus.UCSC.mm10))
#     BSname <- "BSgenome.Mmusculus.UCSC.mm10"
#     BS <- BSgenome.Mmusculus.UCSC.mm10
#     txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#     org <- "org.Mm.eg.db"
#     seqlevels(txdb) <- gsub("chr","",seqlevels(txdb))
#   }else if (species== "macaca") {
#     library(BSgenome.Mfascicularis.NCBI.5.0)
#     seqnames(BSgenome.Mfascicularis.NCBI.5.0) <- gsub("MFA","",seqnames(BSgenome.Mfascicularis.NCBI.5.0))
#     BS <- BSgenome.Mfascicularis.NCBI.5.0
#     BSname <- "BSgenome.Mfascicularis.NCBI.5.0"
#     txdb <- makeTxDbFromGFF("/data1/shenluemou/reference/gtf/Macaca_fascicularis/MFA1912/macaca_MFA1912.gtf",format = "gtf")
#     org <- loadDb("/data1/luofc/hanx/database/Enrichment/GO.db/org.Macaca_fascicularis.eg.sqlite")
#     seqlevels(txdb) <- gsub("chr","",seqlevels(txdb))
#   }else if (species=="Dmelanogaster") {
#     library(BSgenome.Dmelanogaster.UCSC.dm6)
#     seqnames(BSgenome.Dmelanogaster.UCSC.dm6) <- gsub("chr","",seqnames(BSgenome.Dmelanogaster.UCSC.dm6))
#     BS <- BSgenome.Dmelanogaster.UCSC.dm6
#     BSname <- "BSgenome.Dmelanogaster.UCSC.dm6"  
#     library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
#     txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
#     library(org.Dm.eg.db)
#     org <- "org.Dm.eg.db"
#     seqlevels(txdb) <- gsub("chr","",seqlevels(txdb))
#   }else if (species=="zebrafish") {
#     library(BSgenome.Drerio.UCSC.danRer10)
#     seqnames(BSgenome.Drerio.UCSC.danRer10) <- gsub("chr","",seqnames(BSgenome.Drerio.UCSC.danRer10))
#     BS <- BSgenome.Drerio.UCSC.danRer10
#     BSname <- "BSgenome.Drerio.UCSC.danRer10"
#     library(TxDb.Drerio.UCSC.danRer10.refGene)
#     txdb <- TxDb.Drerio.UCSC.danRer10.refGene
#     library(org.Dr.eg.db)
#     org <- org.Dr.eg.db
#     seqlevels(txdb) <- gsub("chr","",seqlevels(txdb))
#   }
# }

#trans
write_pwm2meme <- function(pwm,name){
  #pwm的rowname为"A","C","G","T"
  pwm_to_meme <- function(pwm, name) {
    PPM <- 2^pwm*0.25
    ncol <- ncol(PPM)
    PPM <- t(PPM)
    lines <- c()
    for (i in 1:nrow(PPM)) {
      tmp <- paste0(paste(PPM[i,], collapse="  "))
      lines <- c(lines,tmp)
    }
    header <- paste0("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\nMOTIF ", name)
    footer <- "\n"
    header2 <- paste0("letter-probability matrix: alength= 4 w= ",ncol)
    paste(header, header2,paste(lines, collapse="\n"), footer, sep="\n")
  }
  tmp_meme <- pwm_to_meme(pwm,name)
  writeLines(tmp_meme,paste0(name,".meme"))
}
write_ppm2meme <- function(ppm,name){
  #pwm的rowname为"A","C","G","T"
  ppm_to_meme <- function(ppm, name) {
    ncol <- ncol(ppm)
    ppm <- t(ppm)
    lines <- c()
    for (i in 1:nrow(ppm)) {
      tmp <- paste0(paste(ppm[i,], collapse="  "))
      lines <- c(lines,tmp)
    }
    header <- paste0("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\nMOTIF ", name)
    footer <- "\n"
    header2 <- paste0("letter-probability matrix: alength= 4 w= ",ncol)
    paste(header, header2,paste(lines, collapse="\n"), footer, sep="\n")
  }
  tmp_meme <- ppm_to_meme(ppm,name)
  writeLines(tmp_meme,paste0(name,".meme"))
}
from_meme_get_pwm <- function(tmp_path){
  library(dplyr)
  tmp_line <- readLines(tmp_path)
  tmp_df <- as.data.frame(tmp_line)
  tmp_bac <- tmp_df[8,]
  background_frequencies <- strsplit(tmp_bac, " ")[[1]]
  df_background <- data.frame(matrix(0, ncol = 4, nrow = 1))
  colnames(df_background) <- c("A", "C", "G", "T")
  for (i in 1:4) {
    df_background[, i] <- background_frequencies[i*2]
  }
  tmp_start <- grep("etter-probability matrix", tmp_line,fixed = FALSE)
  tmp_end <- grep("http://jaspar.genereg.net/", tmp_line,fixed = FALSE)
  tmp_matrix <- tmp_df[(tmp_start+1):(tmp_end-1),] %>% na.omit() %>%as.data.frame()
  
  for (i in 1:nrow(tmp_matrix)) {
    tmp <- lapply(tmp_matrix, function(x) {
      strsplit(x, " ")[[i]]
    })
    tmp2 <- as.data.frame(tmp) %>% t() %>% as.data.frame()
    if (i==1) {
      tmp_matrix_df <- tmp2
    }else{
      tmp_matrix_df <- rbind(tmp_matrix_df,tmp2)
    }
  }
  tmp_matrix_df <- tmp_matrix_df %>%  dplyr::select(-V1, -V3, -V5, -V7)
  colnames(tmp_matrix_df) <- c("A","C","G","T")
  tmp_matrix_df <- sapply(tmp_matrix_df,as.numeric)
  df_background <- sapply(df_background, as.numeric)
  res_pwm <- tmp_matrix_df
  for (j in 1:ncol(tmp_matrix_df)) {
    for (k in 1:nrow(tmp_matrix_df)) {
      num <- tmp_matrix_df[k,j]
      num_trans <- log2(num/df_background[j])
      res_pwm[k,j] <- num_trans
    }
  } %>% as.data.frame()
  res_pwm <- t(res_pwm)
  return(res_pwm)
}
from_meme_get_ppm <- function(tmp_path){
  library(dplyr)
  tmp_line <- readLines(tmp_path)
  tmp_df <- as.data.frame(tmp_line)
  tmp_bac <- tmp_df[8,]
  background_frequencies <- strsplit(tmp_bac, " ")[[1]]
  df_background <- data.frame(matrix(0, ncol = 4, nrow = 1))
  colnames(df_background) <- c("A", "C", "G", "T")
  for (i in 1:4) {
    df_background[, i] <- background_frequencies[i*2]
  }
  tmp_start <- grep("etter-probability matrix", tmp_line,fixed = FALSE)
  tmp_end <- grep("URL http:", tmp_line,fixed = FALSE)
  tmp_matrix <- tmp_df[(tmp_start+1):(tmp_end-1),] %>% na.omit() %>%as.data.frame()
  
  for (i in 1:nrow(tmp_matrix)) {
    tmp <- strsplit(tmp_matrix[i,]," ") 
    tmp <- tmp[[1]]
    tmp <- tmp[tmp != ""]
    tmp2 <- as.data.frame(tmp) %>% t() %>% as.data.frame()
    if (i==1) {
      tmp_matrix_df <- tmp2
    }else{
      tmp_matrix_df <- rbind(tmp_matrix_df,tmp2)
    }
  }
  tmp_matrix_df <- tmp_matrix_df %>%  as.data.frame()
  colnames(tmp_matrix_df) <- c("A","C","G","T")
  tmp_matrix_df <- sapply(tmp_matrix_df,as.numeric)
  df_background <- sapply(df_background, as.numeric)
  res_ppm <- t(tmp_matrix_df)
  return(res_ppm)
}
cc_2_gr <- function(cc){
  res <- GRanges(seqnames = cc$seqnames,
                 ranges = IRanges(start = cc$start,
                                  end = cc$end),
                 strand = cc$strand,
                 tpm = cc$tpm,
                 interquantile_width = cc$interquantile_width)
  return(res)
}

#get sequence
location2seq <- function(df,up,down){
  #aim:根据df中location(chr:tss-end:strand)一列,获取up和down_bp的序列
  library(dplyr)
  library(stringr)
  df_chose <- dplyr::select(df,"location")
  df_re <- data.frame(str_split_fixed(df_chose$location,"[:-]",n=4))
  colnames(df_re) <- c("chr","strat","end","strand")
  df_re$strat <- as.numeric(df_re$strat)
  df_re$end <- as.numeric(df_re$end)
  df_re <- na.omit(df_re)
  gr <- GRanges(seqnames = df_re$chr,
                ranges = IRanges(start=df_re$strat,end=df_re$end),
                strand = df_re$strand,
                seqlengths = seqlengths(BS))
  range <- c(-up, down)
  win <- up + down
  gr_seq <- promoters(gr,up,down)
  gr_seq <- gr_seq[width(trim(gr_seq))==win]
  gr_seq <- getSeq(BS,gr_seq)
  return(gr_seq)
}


#get bam
cager_get_domTSS_order <- function(bam_dir){
  #这个function用于确认cager的bam文件读入顺序,以便确定名字,mergeindex
  paths1 <- list.files(bam_dir,full.names = T)
  pathsToInputFiles1 <- paths1[grep(pattern = "*_sorted.bam", paths1)]
  pathsToInputFiles <- c(pathsToInputFiles1)
  names <- sub("_sorted.bam","",basename(pathsToInputFiles))
  return(names)
}
cager_get_LICAGE_exp_pre <- function(bam_dir,mergeindex,mergesample){
    library(CAGEr)
    #function:from bam to get domTSS.grl
    #mergeIndex = c(1,2,3,4,4,4,5,5,5)
    #mergedSampleLabels = c("MII","onepre","oneaf","two","four"))
    paths1 <- list.files(bam_dir,full.names = T)
    pathsToInputFiles1 <- paths1[grep(pattern = "*_sorted.bam", paths1)]
    pathsToInputFiles <- c(pathsToInputFiles1)
    names <- sub("_sorted.bam","",basename(pathsToInputFiles))
    LICAGE <- CAGEexp(genomeName = BSname,
                      inputFiles = pathsToInputFiles,
                      inputFilesType = "bam",
                      sampleLabels = names)
    # LICAGE <- new("CAGEexp", genomeName = BSname, 
    #               inputFiles = pathsToInputFiles, inputFilesType = "bam", 
    #               sampleLabels = names)
    LICAGE <- getCTSS(LICAGE,removeFirstG = FALSE,correctSystematicG = FALSE,
                      sequencingQualityThreshold = 10,
                      mappingQualityThreshold = 20,
                      useMulticore = TRUE,
                      nrCores = 40)
    LICAGE <- mergeSamples(LICAGE,mergeIndex = mergeindex,mergedSampleLabels = mergesample)
    return(LICAGE)
}
cager_get_LICAGE_exp_pre_corr <- function(bam_dir,mergeindex,mergesample){
  library(CAGEr)
  #function:from bam to get domTSS.grl
  #mergeIndex = c(1,2,3,4,4,4,5,5,5)
  #mergedSampleLabels = c("MII","onepre","oneaf","two","four"))
  paths1 <- list.files(bam_dir,full.names = T)
  pathsToInputFiles1 <- paths1[grep(pattern = "*_sorted.bam", paths1)]
  pathsToInputFiles <- c(pathsToInputFiles1)
  names <- sub("_sorted.bam","",basename(pathsToInputFiles))
  LICAGE <- CAGEexp(genomeName = BSname,
                    inputFiles = pathsToInputFiles,
                    inputFilesType = "bam",
                    sampleLabels = names)
  # LICAGE <- new("CAGEexp", genomeName = BSname, 
  #               inputFiles = pathsToInputFiles, inputFilesType = "bam", 
  #               sampleLabels = names)
  LICAGE <- getCTSS(LICAGE,removeFirstG = TRUE,correctSystematicG = TRUE,
                    sequencingQualityThreshold = 10,
                    mappingQualityThreshold = 20,
                    useMulticore = TRUE,
                    nrCores = 40)
  LICAGE <- mergeSamples(LICAGE,mergeIndex = mergeindex,mergedSampleLabels = mergesample)
  return(LICAGE)
}
cager_get_LICAGE_exp_af <- function(licage,al,t,dist){
  licage <- normalizeTagCount(licage, method = "powerLaw", 
                                     fitInRange = c(5, 1000), 
                                     alpha = al, T = t)
  licage <- clusterCTSS(object = licage, threshold = 50, 
                               thresholdIsTpm = TRUE, nrPassThreshold = 1, 
                               method = "distclu", maxDist = dist, 
                               removeSingletons = TRUE, keepSingletonsAbove = 100,
                        useMulticore = TRUE,nrCores = 20)
  licage <- cumulativeCTSSdistribution(licage, clusters = "tagClusters",
                                       useMulticore = TRUE,nrCores = 20)
  licage <- quantilePositions(licage, clusters = "tagClusters", qLow = 0.1, qUp = 0.9,
                              useMulticore = TRUE,nrCores = 20)
  return(licage)
}
cager_get_LICAGE_exp_af_para <- function(licage,al,t,len){
  licage <- normalizeTagCount(licage, method = "powerLaw", 
                              fitInRange = c(5, 1000), 
                              alpha = al, T = t)
  licage <- clusterCTSS(object = licage, threshold = 50, 
                        thresholdIsTpm = TRUE, nrPassThreshold = 1, 
                        method = "paraclu",maxLength = len, minStability = 1, 
                        removeSingletons = TRUE, keepSingletonsAbove = 100,
                        reduceToNonoverlapping = TRUE)
  licage <- cumulativeCTSSdistribution(licage, clusters = "tagClusters")
  licage <- quantilePositions(licage, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
  return(licage)
}
cager_get_LICAGE_exp_pre_notmer_tsscor <- function(bam_dir){
  #st和mt一般取20
  library(CAGEr)
  #function:from bam to get domTSS.grl
  #mergeIndex = c(1,2,3,4,4,4,5,5,5)
  #mergedSampleLabels = c("MII","onepre","oneaf","two","four"))
  paths1 <- list.files(bam_dir,full.names = T)
  pathsToInputFiles1 <- paths1[grep(pattern = "*_sorted.bam", paths1)]
  pathsToInputFiles <- c(pathsToInputFiles1)
  names <- sub("_sorted.bam","",basename(pathsToInputFiles))
  # LICAGE <- CAGEexp(genomeName = BSname,
  #                   inputFiles = pathsToInputFiles,
  #                   inputFilesType = "bam",
  #                   sampleLabels = names)
  LICAGE <- new("CAGEexp", genomeName = BSname,
                inputFiles = pathsToInputFiles, inputFilesType = "bam",
                sampleLabels = names)
  LICAGE <- getCTSS(LICAGE,removeFirstG = FALSE,correctSystematicG = FALSE,
                    sequencingQualityThreshold = 10,#10
                    mappingQualityThreshold = 20,#20
                    useMulticore = TRUE,
                    nrCores = 1)
  return(LICAGE)
}
cager_get_LICAGE <- function(bam_dir,mergeindex,mergesample){
    library(CAGEr)
    #function:from bam to get domTSS.grl
    #mergeindex的格式为mergeIndex = c(1,2,3,4,4,4,5,5,5)
    #mergesample的格式为mergedSampleLabels = c("MII","onepre","oneaf","two","four"))
    paths1 <- list.files(bam_dir,full.names = T)
    pathsToInputFiles1 <- paths1[grep(pattern = "*_sorted.bam", paths1)]
    pathsToInputFiles <- c(pathsToInputFiles1)
    names <- sub("_sorted.bam","",basename(pathsToInputFiles))
    LICAGE <- new("CAGEset", genomeName = BSname, 
                  inputFiles = pathsToInputFiles, inputFilesType = "bam", 
                  sampleLabels = names)
    getCTSS(LICAGE)
    mergeSamples(LICAGE,mergeIndex = mergeindex,mergedSampleLabels = mergesample)
    normalizeTagCount(LICAGE, method = "powerLaw", fitInRange = c(100, 1000), alpha = 1.14, T = 1000000)
    #exportCTSStoBedGraph(LICAGE, values = "normalized", format = "bedGraph", oneFile = FALSE)
    clusterCTSS(object = LICAGE, threshold = 1, thresholdIsTpm = TRUE, nrPassThreshold = 1, 
                method = "distclu", maxDist = 5, removeSingletons = TRUE, keepSingletonsAbove = 5)
    cumulativeCTSSdistribution(LICAGE, clusters = "tagClusters")
    quantilePositions(LICAGE, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
    return(LICAGE)
}

#progress
cager_get_domTSS <- function(LICAGE){
  #sampleLabels <- unname(sampleLabels(LICAGE))
  sampleLabels <- sampleLabels(LICAGE)
  tc.l <- lapply(sampleLabels, 
                 function (x) tagClustersGR(LICAGE, sample = x, 
                                            returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9))
  names(tc.l) <- sampleLabels
  names <- sampleLabels
  tc.grl <- GRangesList()
  tc.grl <- lapply(tc.l, function(x) GRanges(seqnames = x@seqnames, 
                                             ranges = x@ranges,
                                             strand = x@strand,
                                             nr_ctss = x$nr_ctss,
                                             dominant_ctss = x$dominant_ctss,
                                             tpm = x$tpm.dominant_ctss,
                                             tpm.dominant_ctss = x$tpm.dominant_ctss,
                                             q_0.1 = x$q_0.1,
                                             q_0.9 = x$q_0.9,
                                             interquantile_width = x$interquantile_width,
                                             seqlengths = seqlengths(BS)))
  
  # add seqinfo information                
  for (i in 1:length(tc.grl)) {
    seqinfo(tc.grl[[i]]) <- seqinfo(BS) 
  }
  domTSS.grl <- list()
  for (i in 1:length(tc.grl)) {
    domTSS.grl[[i]] <- GRanges(seqnames = seqnames(tc.grl[[i]]),
                               ranges = IRanges(start = tc.grl[[i]]$dominant_ctss, 
                                                end = tc.grl[[i]]$dominant_ctss),
                               strand = strand(tc.grl[[i]]),
                               nr_ctss = tc.grl[[i]]$nr_ctss,
                               dominant_ctss = tc.grl[[i]]$dominant_ctss,
                               tpm = tc.grl[[i]]$tpm, 
                               tpm.dominant_ctss = tc.grl[[i]]$tpm.dominant_ctss, 
                               q_0.1 = tc.grl[[i]]$q_0.1,
                               q_0.9 = tc.grl[[i]]$q_0.9,
                               interquantile_width = tc.grl[[i]]$interquantile_width,
                               tc_start = tc.grl[[i]]$dominant_ctss,
                               tc_end = tc.grl[[i]]$dominant_ctss)
    
    seqlevels(domTSS.grl[[i]]) <- seqlevels(tc.grl[[i]])
    seqlengths(domTSS.grl[[i]]) <- seqlengths(tc.grl[[i]])
    genome(domTSS.grl[[i]]) <- genome(BS)
  }
  names(domTSS.grl) <- names
  return(domTSS.grl)
}
cager_aggregate <- function(licage,maxd){
  #function:begian from licage to aggregate,cumulate and quantile licage
  #from:licage
  #next progress:consensusClustersDESeq2 or get_CC_lis
  licage_aggqua <- aggregateTagClusters(licage,tpmThreshold = 5,
                                           qLow = 0.1,qUp = 0.9,maxDist = maxd)
  licage_aggqua <- cumulativeCTSSdistribution(licage_aggqua,clusters = "consensusClusters")
  licage_aggqua <- quantilePositions(licage_aggqua,clusters = "consensusClusters", 
                                        qLow = 0.1, qUp = 0.9)
  return(licage_aggqua)
}
cager_matrix <- function(cager_agg){
  cager_tagclu <- tagClustersGR(cager_agg,returnInterquantileWidth = TRUE,qLow=0.1,qUp=0.9) %>%
    as.data.frame()
  cager_tagclu$location <- sprintf("%s:%s-%s:%s",
                                   cager_tagclu$seqnames,
                                   cager_tagclu$start,
                                   cager_tagclu$end,
                                   cager_tagclu$strand)
  cager_mer_tag <- dplyr::select(cager_tagclu,"group_name","location",
                                 "tpm.dominant_ctss")
  cager_tag_matrix <- reshape2::dcast(cager_mer_tag,location ~ group_name,
                                          value.var = "tpm.dominant_ctss",fill = 0)
  return(cager_tag_matrix)
}
cager_cc_matrix <- function(cager_merge_mm_cc_fen) {
  # Convert list elements to data frames and add group_name column
  cager_merge_mm_cc_fen_df <- lapply(cager_merge_mm_cc_fen, as.data.frame)
  # Add group_name column to each dataframe
  for (i in seq_along(cager_merge_mm_cc_fen)) {
    cager_merge_mm_cc_fen_df[[i]]$group_name <- names(cager_merge_mm_cc_fen)[i]
  }
  # Bind rows of list elements into a single dataframe
  cager_merge_mm_cc_fm <- do.call(rbind, cager_merge_mm_cc_fen_df)
  # Select necessary columns
  cager_merge_mm_cc_fm_select <- dplyr::select(cager_merge_mm_cc_fm,
                                               "group_name", "consensus.cluster", "tpm")
  # Reshape data frame into matrix format
  cager_merge_mm_cc_matrix <- reshape2::dcast(cager_merge_mm_cc_fm_select,
                                              consensus.cluster ~ group_name,
                                              value.var = "tpm", fill = 0)
  new_name <- paste("x", cager_merge_mm_cc_matrix$consensus.cluster, sep="_")
  cager_merge_mm_cc_matrix$consensus.cluster <- new_name
  return(cager_merge_mm_cc_matrix)
}

#annotation
cc2anno <- function(cc){
  #function:输入cc,获取该cc的注释
  seqlevels(txdb) <- gsub("chr","",seqlevels(txdb))
  gene <- GenomicFeatures::genes(txdb)
  peakAnno.l <- annotatePeak(cc, TxDb = txdb,  annoDb = org,
                             sameStrand = TRUE, verbose = TRUE)
  peakan <- peakAnno.l@anno %>% as.data.frame()
  peakan_name <- as.data.frame(rownames(peakan))
  colnames(peakan_name) <- "location"
  all_anno <- cbind(peakan_name,peakan)
  return(all_anno)
}

#quality_control
dom_peakann_res <- function(domTSS){
  library(ggplot2)
  library(CAGEr)
  library(ChIPseeker)
  library(RColorBrewer)
  peakAnno_list <- lapply(domTSS, 
                          function(x) annotatePeak(x, TxDb = txdb,  tssRegion = c(-1000, 1000), annoDb = org, sameStrand = TRUE, verbose = FALSE))
  names <- names(domTSS)
  names(peakAnno_list) <- names
  return(peakAnno_list)
}
dom_peakann_plot <- function(domTSS){
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
  for (i in 1:length(names)) {
    feats.l[[i]]$sample <- rep(names[[i]], nrow(feats.l[[i]]))
  }
  feats.df <- do.call("rbind", feats.l)
  feats.df$sample <- factor(feats.df$sample, levels = names[length(names):1])
  col1 <- brewer.pal(5, "Greys")[3:5]
  col2 <- brewer.pal(8, "Blues")[8:6]
  col <- c(col2, col1)
  # new colour scheme
  col <- c("#F8F2AB", "#B4D6A4", "#67BC9A", "#13B0A5", "#0071A7", "#3E606F", "#88F9D4", "#18C29C",
           "#0B877D", "#126872", "#031727")
  # set feature factors
  feature_factors <- c("Promoter", "Promoter (<=1kb)", "Promoter (1-3kb)", "5' UTR", "1st Exon", "Other Exon", 
                       "1st Intron", "Other Intron", "3' UTR", "Downstream (<=3kb)", "Distal Intergenic")
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
dom_promoter_width <- function(domTSS,ncol){
  names <- names(domTSS)
  iq_all.l <- lapply(domTSS, function(x) data.frame(iq_width = x$interquantile_width))
  names(iq_all.l) <- names
  for (i in 1:length(names)) {
    iq_all.l[[i]]$sample <- rep(names[[i]], nrow(iq_all.l[[i]]))
  }
  iq_all.df <- do.call(rbind, iq_all.l)
  iq_all.df$sample <- factor(iq_all.df$sample, levels = names)
  # plotting
  p <- ggplot(iq_all.df, aes(iq_width)) +
    geom_histogram(aes(x = iq_width, y = ..ncount..* 100), binwidth = 1, 
                   fill = "gray60", col = "black", size = 0.1) +
    theme_bw() +
    theme(text = element_text(size = 16, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black", size = 14),
          axis.text.y = element_text(colour = "black", size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(colour = "black", fill = "gray87")) +
    labs(x = "Interquantile width", y = "percent") +
    xlim(0, 100)
  p1 <- p + facet_wrap(~ sample, ncol = ncol)
  return(p1)
}
plot_cager_TFctss <- function(licage,control){
  CTSSnorm <- CTSSnormalizedTpmDF(licage)
  CTSSnorm <- as.data.frame(CTSSnorm)
  samples <- setdiff(sampleLabels(licage),control)
  tpm_pairs.dfl <- lapply(samples, function(x) {
    df <- CTSSnorm[, c(x, control)]
    df <- df[rowSums(df) > 0, ]
    return(df) })
  names(tpm_pairs.dfl) <- samples
  TP_FP_l <- list()
  TP_FP.gg <- list()
  for(j in 1:length(tpm_pairs.dfl)){
    threshold <- seq(0, 10, 0.05) 
    for (i in 1:length(threshold)) {
      tpm_filt <- tpm_pairs.dfl[[j]][tpm_pairs.dfl[[j]][, 1] >= threshold[i], ]
      TP_FP_l[[i]] <- data.frame("TP" = tpm_filt[,1] > 0 & tpm_filt[,2] > 0,
                                 "FP" = tpm_filt[,1] > 0 & tpm_filt[,2] == 0)
      
    }
    TP_FP.gg[[j]] <- as.data.frame(do.call(rbind, lapply(TP_FP_l, colSums)))
    TP_FP.gg[[j]]$sample <- rep(names(tpm_pairs.dfl[j]), times = nrow(TP_FP.gg[[j]]))
  }
  # flatten for plotting
  TP_FP.gg <- as.data.frame(do.call(rbind, TP_FP.gg))
  TP_FP.gg$sample <- factor(TP_FP.gg$sample, levels = samples)
  
  library(viridis)
  col <- magma(10)
  
  p <- ggplot(TP_FP.gg, aes(x = FP, y = TP, fill = sample)) +
    geom_line(aes(col = sample), size = 1, alpha = 0.7) +
    geom_point(shape = 21, size = 1, alpha = 0.7) +
    scale_fill_manual(values = col) +
    scale_colour_manual(values = col) +
    theme(text = element_text(size = 12,
                              family = "Helvetica"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 12, family = "Helvetica", color = "black"),
          legend.title = element_text(size = 12, family = "Helvetica", colour = "black"),
          legend.text = element_text(size = 12, family = "Helvetica", colour = "black")) +
    xlab("False positive CTSSs") +
    ylab("True positive CTSSs") +
    xlim(c(0, 80000)) +
    ylim(c(0, 80000))
  return(p)
}
tocorrPlot <- function(data, cor, xlab, ylab) {
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
          legend.text = element_text(size = 14, family = "Helvetica", colour = "black")) +
    xlab(xlab) +
    ylab(ylab) +
    xlim(c(0, 4.5)) +
    ylim(c(0, 4.5)) +
    annotate("text", x = 0, y = 4.5, size = 6, hjust = 0, label = paste("R = ", cor))  +
    coord_fixed(ratio = 1)
  return(p)
}

plot_cager_cor <- function(LICAGE,control,treatment){
  tag <- CTSSnormalizedTpmDF(LICAGE) %>% as.data.frame()
  filtr <- tag[,c(control,treatment)]
  filtr <- filtr[(filtr[, 1] >= 1 | filtr[, 2] >= 1), ]
  cor.l <- round(cor(filtr), digits = 3)
  filtr.data.log.l <- log10(filtr + 1)
  names(filtr) <- treatment
  names(filtr.data.log.l) <- treatment
  labels_y <- treatment
  labels_x <- control
  p <- tocorrPlot(data = filtr.data.log.l, xlab = labels_x, ylab = labels_y, 
                  cor = cor.l[2])
  return(p)
}
#tag_cluster
plot_cager_cluster <- function(LICAGE,width,height){
  aggregateTagClusters(LICAGE, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)
  ce <- getExpressionProfiles(LICAGE, what = "consensusClusters", tpmThreshold = 10, nrPassThreshold = 1, method = "som", xDim = width, yDim = height)
  p <- plotExpressionProfiles(ce, what = "consensusClusters")
  return(p)
}
tag_cluster <- function(LICAGE,cluster,width,height){
  aggregateTagClusters(LICAGE, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)
  ce <- getExpressionProfiles(LICAGE, what = "consensusClusters", tpmThreshold = 10, nrPassThreshold = 1, method = "som", xDim = width, yDim = height)
  extract.cl <- consensusClustersGR(ce) |> subset(consensusClustersGR(ce)$exprClass ==  cluster)
  extract.grl <- makeGRangesFromDataFrame(extract.cl, keep.extra.columns = TRUE, 
                                          ignore.strand = FALSE,
                                          seqnames.field = "chr",
                                          start.field = "start",
                                          end.field = "end",
                                          strand.field = "strand")
  return(list(extract.cl,extract.grl))
}
tag_cluster_merge <- function(tag1,tag2){
  #该function用于将tag_cluter获取的两个S4 GRange结构进行合并
  tmp1 <- as.data.frame(tag1[[2]])
  tmp2 <- as.data.frame(tag2[[2]])
  tmp <- rbind(tmp1,tmp2)
  data2 <- GRanges(seqnames = tmp$seqnames,
                   ranges = IRanges(start = tmp$start, 
                                    end = tmp$end),
                   strand = tmp$strand,
                   consensus.cluster = tmp$consensus.cluster,
                   tpm = tmp$tpm,
                   G1 = tmp$G1,
                   G2 = tmp$G2,
                   S = tmp$S,
                   expression_class = tmp$expression_class)
  data3 <- rbind(tag1[[1]],tag2[[1]])
  return(list(data3,data2))
}
tag_cluster_width <- function(tag_cluster_licage,names){
  #该函数输入tag_cluster的LICAGE,获取该cluster号的promoter_width
  tag_cluster_licage$interquantile_width <- tag_cluster_licage$end+1 - tag_cluster_licage$start
  iq_all.l <- data.frame(iq_width = tag_cluster_licage$interquantile_width)
  iq_all.l <- list(iq_all.l)
  names(iq_all.l) <- names
  for (i in 1:length(names)) {
    iq_all.l[[i]]$sample <- rep(names[[i]], nrow(iq_all.l[[i]]))
  }
  iq_all.df <- do.call(rbind, iq_all.l)
  iq_all.df$sample <- factor(iq_all.df$sample, levels = names)
  # plotting
  p <- ggplot(iq_all.df, aes(iq_width)) +
    geom_histogram(aes(x = iq_width, y = ..ncount..* 100), binwidth = 1, 
                   fill = "gray60", col = "black", size = 0.1) +
    theme_bw() +
    theme(text = element_text(size = 16, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black", size = 14),
          axis.text.y = element_text(colour = "black", size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(colour = "black", fill = "gray87")) +
    labs(x = "Interquantile width", y = "percent") +
    xlim(0, 100)
  p1 <- p + facet_wrap(~ sample, ncol = 5)
  df <- data.frame(sample1=c(1,2),sample2=c(1,2))
  for (i in 1:length(names)) {
    tmp <- iq_all.l[[i]]
    sharp <- nrow(tmp[which(tmp$iq_width<10),])
    broad <- nrow(tmp[which(tmp$iq_width>=10),])
    df[1,i] <- sharp
    df[2,i] <- broad
  }
  colnames(df) <- names
  rownames(df) <- c("sharp","broad")
  all_re <- list(p1,df)
  return(all_re)
}
tag_from_cc <- function(tagcluster,ccluster){
  #根据consen_cluster,从对应的tagcluster中提取出对应tagcluster
  #例如从MII高表达的consen里面提取出对应的tagcluster
  overlaps <- findOverlaps(tagcluster, ccluster) 
  res <- tagcluster[queryHits(overlaps)]
  return(res)
}


#others
exp_to_deseq2 <- function(LICAGE,species,group,pv,fc){
  #该函数的前置函数为cager_get_LICAGE_exp,目的为获取差异表达TC
  if (species == "human") {
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  }else if (species == "mouse") {
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  }else if (species== "macaca") {
    txdb <- makeTxDbFromGFF("/data1/shenluemou/reference/gtf/Macaca_fascicularis/MFA1912/macaca_MFA1912.gtf",format = "gtf")
  }else if (species=="Dmelanogaster") {
    library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
    txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
  }
  seqlevels(txdb) <- gsub("chr","",seqlevels(txdb))
  aggregateTagClusters(LICAGE, 
                       tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)
  consensusClustersGR(LICAGE)
  gene <- genes(txdb)
  annotateConsensusClusters(LICAGE, gene)
  cumulativeCTSSdistribution(LICAGE, clusters = "consensusClusters", useMulticore = TRUE)
  quantilePositions(LICAGE, clusters = "consensusClusters", qLow = 0.1, qUp = 0.9, useMulticore = TRUE)
  LICAGE$group <- group#控制组,实验组,如c("g1","g1","g2","g2")
  dds <- consensusClustersDESeq2(LICAGE, ~group)
  dds1 <- DESeq(dds)
  res <- results(dds1)
  res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  res1 <- na.omit(res1)
  res1$threshold <- ifelse(res1$pvalue>pv,"no",ifelse(abs(res1$log2FoldChange)<fc,"no",
                                                      ifelse(res1$log2FoldChange>fc,"up","down")))
  res1$symbol <- rownames(res1)
  res2 <- res1 %>% 
      separate(symbol, sep = ":",into = c("char", "range","lian")) %>%
      separate(range,sep = "-",into = c("start","end"))
  res2$start <- as.numeric(res2$start)
  res2$end <- as.numeric(res2$end)
  res2$Interquantil_width <- res2$end - res2$start +1
  return(res2)
}
cager_diffres2dom <- function(df){
  #将cager_different_analysis_pad/cager_different_analysis的res结果中的某一个转换成为dom
  df$symbol <- rownames(df)
  df2 <- df %>% separate(symbol, sep = ":",into = c("char", "range","lian")) %>%
    separate(range,sep = "-",into = c("start","end"))
  #df2 <- na.omit(df2)
  df2$start <- as.numeric(df2$start)
  df2$end <- as.numeric(df2$end)
  df2$end[is.na(df2$end)] <- df2$start[is.na(df2$end)]
  dom1 <- GRanges(seqnames = df2$char,
                  ranges = IRanges(start = df2$start,end = df2$end),
                  strand = df2$lian)
  return(dom1)
}
cager_matrix2dom <- function(df){
  #将cager_matrix转换成为dom
  df$symbol <- df$location
  df2 <- df %>% separate(symbol, sep = ":",into = c("char", "range","lian")) %>%
    separate(range,sep = "-",into = c("start","end"))
  #df2 <- na.omit(df2)
  df2$start <- as.numeric(df2$start)
  df2$end <- as.numeric(df2$end)
  df2$end[is.na(df2$end)] <- df2$start[is.na(df2$end)]
  dom1 <- GRanges(seqnames = df2$char,
                  ranges = IRanges(start = df2$start,end = df2$end),
                  strand = df2$lian)
  return(dom1)
}
df2gr <- function(df_des){
  data2 <- GRanges(seqnames = df_des$seqnames,
                   ranges = IRanges(start = df_des$start,
                                    end = df_des$end),
                   strand = df_des$strand,
                   interquantile_width = df_des$Interquantil_width)
  return(data2)
}
dfdeseq2_to_width <- function(df_des,names){
  iq_all.l <- data.frame(iq_width = df_des$Interquantil_width)
  iq_all.l <- list(iq_all.l)
  names(iq_all.l) <- names
  for (i in 1:length(names)) {
    iq_all.l[[i]]$sample <- rep(names[[i]], nrow(iq_all.l[[i]]))
  }
  iq_all.df <- do.call(rbind, iq_all.l)
  iq_all.df$sample <- factor(iq_all.df$sample, levels = names)
  # plotting
  p <- ggplot(iq_all.df, aes(iq_width)) +
    geom_histogram(aes(x = iq_width, y = ..ncount..* 100), binwidth = 1, 
                   fill = "gray60", col = "black", size = 0.1) +
    theme_bw() +
    theme(text = element_text(size = 16, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black", size = 14),
          axis.text.y = element_text(colour = "black", size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(colour = "black", fill = "gray87")) +
    labs(x = "Interquantile width", y = "percent") +
    xlim(0, 100)
  p1 <- p + facet_wrap(~ sample, ncol = 5)
  df <- data.frame(sample1=c(1,2),sample2=c(1,2))
  for (i in 1:length(names)) {
    tmp <- iq_all.l[[i]]
    sharp <- nrow(tmp[which(tmp$iq_width<10),])
    broad <- nrow(tmp[which(tmp$iq_width>=10),])
    df[1,i] <- sharp
    df[2,i] <- broad
  }
  colnames(df) <- names
  rownames(df) <- c("sharp","broad")
  all_re <- list(p1,df)
  return(all_re)
}
upset_interse <- function(seq1,motif1,motif1_name,seq2,motif2,motif2_name){
  #----motif1_seq1------
  simple=1
  motif1_sy <- c()
  all_sy <- c()
  for (t in 1:length(seq1)) {
    tmp_seq <- DNAString(as.character(seq1[t]))
    tmp_count <- countPWM(motif1,tmp_seq,min.score = "80%")
    if (tmp_count > 0) {
      motif1_sy <- c(motif1_sy,simple)
    }
    simple=simple+1
    all_sy <- c(all_sy,simple)
  }
  #----motif1_seq2------
  motif2_sy <- c()
  simple=1
  for (t in 1:length(seq2)) {
    tmp_seq <- DNAString(as.character(seq2[t]))
    tmp_count <- countPWM(motif2,tmp_seq,min.score = "80%")
    if (tmp_count > 0) {
      motif2_sy <- c(motif2_sy,simple)
    }
    simple=simple+1
  }
  motif_all <- union(motif1_sy,motif2_sy)
  others_sy <- setdiff(all_sy,motif_all)
  all <- list(motif1 = motif1_sy,motif2=motif2_sy,others=others_sy)
  names(all) <- c(motif1_name,motif2_name,"others")
  p1 <- upset(fromList(all),text.scale = c(2,5,2,2,4,5))
  return(p1)
}
cager_diffclu_seq <- function(species,diff_clu,up,down){
  #该功能是根据clu的dominat_ctss的矩阵结果获取序列
  #1位:diff_clu；2位:domdf_list ; 3位:上游; 4位:下游
  if (species=="mouse") {
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
    seqnames(BSgenome.Mmusculus.UCSC.mm10) <- gsub("chr","",seqnames(BSgenome.Mmusculus.UCSC.mm10))
    BSname <- "BSgenome.Mmusculus.UCSC.mm10"
    BS <- BSgenome.Mmusculus.UCSC.mm10
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    org <- "org.Mm.eg.db"
  }else if (species== "macacaf") {
    library(BSgenome.Mfascicularis.NCBI.5.0)
    seqnames(BSgenome.Mfascicularis.NCBI.5.0) <- gsub("MFA","",seqnames(BSgenome.Mfascicularis.NCBI.5.0))
    BS <- BSgenome.Mfascicularis.NCBI.5.0
    BSname <- "BSgenome.Mfascicularis.NCBI.5.0"
  }else if (species=="human") {
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    seqnames(BSgenome.Hsapiens.UCSC.hg38) <- gsub("chr","",seqnames(BSgenome.Hsapiens.UCSC.hg38))
    BSname <- "BSgenome.Hsapiens.UCSC.hg38"
    BS <- BSgenome.Hsapiens.UCSC.hg38
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    org <- "org.Hs.eg.db"
  }
  seq_li <- list()
  for (i in 1:length(diff_clu)) {
    tmp_loc <- diff_clu[[i]]
    tmp <- GRanges(seqnames = tmp_loc$seqnames, 
                   ranges = IRanges(start = tmp_loc$start, end = tmp_loc$end),
                   strand = tmp_loc$strand,
                   nr_ctss = tmp_loc$nr_ctss,
                   dominant_ctss = tmp_loc$dominant_ctss,
                   tpm = tmp_loc$tpm,
                   tpm.dominant_ctss = tmp_loc$tpm.dominant_ctss,
                   q_0.1 = tmp_loc$q_0.1,
                   q_0.9 = tmp_loc$q_0.9,
                   interquantile_width = tmp_loc$interquantile_width,
                   tc_start = tmp_loc$tc_start,
                   tc_end = tmp_loc$tc_end,
                   seqlengths = seqlengths(BS))
    range <- c(-up, down)
    win <- up + down
    dom_tmp <- promoters(tmp,up=up,down=down)
    dom_tmp_win <- dom_tmp[width(trim(dom_tmp))==win]
    seq_li[[i]] <- getSeq(BS,dom_tmp_win)
  }
  return(seq_li)
}
cager_seq_scan_genekind_plot <- function(seq_li,motif,motif_name,pwm_score){
  motif_pic <- list()
  the_name <- names(seq_li)
  for (i in 1:length(seq_li)) {
    sum=0
    for (t in 1:length(seq_li[[i]])) {
      tmp_seq_file <- seq_li[[i]]
      tmp_seq <- DNAString(as.character(tmp_seq_file[t]))
      tmp_count <- countPWM(motif,tmp_seq,min.score = paste0(pwm_score,"%"))
      if (tmp_count > 0) {
        sum=sum+1
      }
    }
    motif_per <- sum/length(seq_li[[i]])
    motif_topie_name <- c(motif_name,paste0("not ",motif_name))
    motif_topie_percent <- c(motif_per,1-motif_per)
    motif_topie_percent <- round(motif_topie_percent*100,2)
    motif_topie_name <- paste(motif_topie_name," | ",motif_topie_percent,"%",sep = "")
    png(paste0(the_name[i],"_",motif_name,"_",pwm_score,"_seqkind.png"),height = 800,width = 800)
    motif_pic[[i]] <- pie(motif_topie_percent,labels = motif_topie_name,main=paste0(the_name[i],"_",motif_name),
                          radius = 0.8,cex=1.0,clockwise = TRUE)
    dev.off()
  }
}
cager_seq_scan_transnumber_plot <- function(seq_li,diff_clu,motif,motif_name){
  motif_pic <- list()
  the_name <- names(seq_li)
  for (i in 1:length(seq_li)) {
    sum=0
    for (t in 1:length(seq_li[[i]])) {
      tmp_seq_file <- seq_li[[i]]
      tmp_seq <- DNAString(as.character(tmp_seq_file[t]))
      tmp_count <- countPWM(motif,tmp_seq,min.score = "80%")
      if (tmp_count > 0) {
        sum=sum+diff_clu[[i]]$nr_ctss[t]
      }
    }
    all <- sum(diff_clu[[i]]$nr_ctss)
    motif_per <- sum/all
    motif_topie_name <- c(motif_name,paste0("not ",motif_name))
    motif_topie_percent <- c(motif_per,1-motif_per)
    motif_topie_percent <- round(motif_topie_percent*100,2)
    motif_topie_name <- paste(motif_topie_name," | ",motif_topie_percent,"%",sep = "")
    png(paste0(the_name[i],"_",motif_name,"_transnumber.png"),height = 800,width = 800)
    motif_pic[[i]] <- pie(motif_topie_percent,labels = motif_topie_name,main=paste0(the_name[i],"_",motif_name),
                          radius = 0.8,cex=1.0,clockwise = TRUE)
    dev.off()
  }
}
cager_seq_scanpattern_genekind_plot <- function(seq_li,motif,motif_name){
  motif_pic <- list()
  the_name <- names(seq_li)
  for (i in 1:length(seq_li)) {
    sum=0
    for (t in 1:length(seq_li[[i]])) {
      tmp_seq_file <- seq_li[[i]]
      tmp_seq <- DNAString(as.character(tmp_seq_file[t]))
      tmp_count <- countPattern(motif,tmp_seq)
      if (tmp_count > 0) {
        sum=sum+1
      }
    }
    motif_per <- sum/length(seq_li[[i]])
    motif_topie_name <- c(motif_name,paste0("not ",motif_name))
    motif_topie_percent <- c(motif_per,1-motif_per)
    motif_topie_percent <- round(motif_topie_percent*100,2)
    motif_topie_name <- paste(motif_topie_name," | ",motif_topie_percent,"%",sep = "")
    png(paste0(the_name[i],"_",motif_name,"_",pwm_score,"_seqkind.png"),height = 800,width = 800)
    motif_pic[[i]] <- pie(motif_topie_percent,labels = motif_topie_name,main=paste0(the_name[i],"_",motif_name),
                          radius = 0.8,cex=1.0,clockwise = TRUE)
    dev.off()
  }
}
df_to_S4 <- function(data1){
  data2 <- GRanges(seqnames = data1$seqnames,
                   ranges = IRanges(start = data1$dominant_ctss, 
                                    end = data1$dominant_ctss),
                   strand = data1$strand,
                   nr_ctss = data1$nr_ctss,
                   dominant_ctss = data1$dominant_ctss,
                   tpm = data1$tpm, 
                   tpm.dominant_ctss = data1$tpm.dominant_ctss, 
                   q_0.1 = data1$q_0.1,
                   q_0.9 = data1$q_0.9,
                   interquantile_width = data1$interquantile_width,
                   tc_start = data1$tc_start,
                   tc_end = data1$tc_end)
  return(data2)
}
clulis_to_gene <- function(cluster){
  datafram <- list()
  for (i in 1:length(cluster)) {
    datafram[[i]] <- as.data.frame(rownames(as.data.frame(cluster[[i]])))
    colnames(datafram[[i]]) <- "gene_id"
  }
  return(datafram)
}
cager_tss_dom <- function(tss,dom){
  dom <- merge(dom,tss,by="dominant_ctss")
  colnames(dom)[2] <- "seqnames"
  return(dom)
}
cager_te_count <- function(test_yuan,test_x2){
  #function:根据cager结果与te的结果统计不同te的count
  #input:1.cager表格 2.te结果
  test_x2 <- na.omit(test_x2)
  test_yuan <- unite(test_yuan,"sym",c("seqnames","start"),sep = "_",remove = FALSE)
  test_x2 <- unite(test_x2,"sym",c("V1","V2"),sep = "_",remove = FALSE)
  test <- merge(test_yuan,test_x2,by="sym")
  test_te_all <- select_(test,"seqnames","start","end","nr_ctss","V5")
  test_te_nr <- select_(test_te_all,"V5","nr_ctss")
  test_count <- aggregate(nr_ctss ~ V5, data = test_te_nr, sum)
  new_row <- nrow(test_count)+1
  gene_all <- colSums(test_yuan[8])
  gene_only <- gene_all-colSums(test_count[2])
  test_count[new_row,]<-c("no_TE",gene_only)
  return(test_count)
}
cluter2seq <- function(up,down,cluster_sp){
  library(seqPattern)
  cluster_sp_range <- as.data.frame(cluster_sp)
  tmp <- GRanges(seqnames = cluster_sp_range$seqnames, 
                 ranges = IRanges(start = cluster_sp_range$start, end = cluster_sp_range$end),
                 strand = cluster_sp_range$strand,
                 interquantile_width = cluster_sp_range$width,
                 seqlengths = seqlengths(BS))
  range <- c(-up, down)
  win <- up + down
  dom_tmp <- promoters(tmp,up=up,down=down)
  #dom_tmp_win <- dom_tmp[width(trim(dom_tmp))==win]
  seq_li <- getSeq(BS,dom_tmp)
  return(seq_li)
}
seq2logo <- function(seq){
  library("ggseqlogo")
  all <- c()
  for (i in 1:length(seq)) {
    all <- c(all,as.character(seq[[i]]))
  }
  p <- ggseqlogo(all)
  return(p)
}
select_gene_range <- function(GRange,genedf,species){
  library(ChIPseeker)
  #该函数是根据输入的genedf(SYMBOL),从指定的Granges获取这些SYMBOL的起始终止等信息
  #input:GRanges, genelist,species
  #output:df1(select_gene的信息,not_select_gene的信息)
  #add:txdb,org等至此开始,一开始就load,之后的函数不再进行load
  colnames(genedf) <- "SYMBOL"
  peakAnno.l <- annotatePeak(GRange, TxDb = txdb,  annoDb = NULL, sameStrand = TRUE, verbose = FALSE)
  clusters_anno.grl <- as.data.frame(peakAnno.l@anno)
  gene_enz <- select(clusters_anno.grl,"geneId")
  gene_all <- genesymol_all(species,"ENTREZID",gene_enz)
  colnames(gene_all)[3] <- "geneId"
  clusters_anno.grl_sym <- merge(clusters_anno.grl,gene_all,by="geneId")
  select_gene <- merge(clusters_anno.grl_sym,genedf,by="SYMBOL")
  all_gene_list <- clusters_anno.grl_sym$SYMBOL
  select_gene_list <- genedf$SYMBOL
  other_gene_list <- setdiff(all_gene_list,select_gene_list)
  other_gene_df <- as.data.frame(other_gene_list)
  colnames(other_gene_df) <- "SYMBOL"
  other_gene <- merge(clusters_anno.grl_sym,other_gene_df,by="SYMBOL")
  the_list <- list(select_gene,other_gene)
  names(the_list) <- c("select_gene","other_gene")
  return(the_list)
}

#different_analysis
cager_different_analysis <- function(nss_2_licage_fen,group1,pv,fc){
  #以licage_agg为起始,对licage进行差异分析,返回火山图和差异矩阵
  nss_2_licage_fen$group <- group1
  nss_2_th_fen_dds <- consensusClustersDESeq2(nss_2_licage_fen,~group)
  nss_2_th_fen_dds1 <- DESeq(nss_2_th_fen_dds, fitType = 'mean', 
                             minReplicatesForReplace = 7, 
                             parallel = FALSE) 
  ment <- c("group",group1[length(group1)],group1[1])
  nss_2_th_fen_res <- DESeq2::results(nss_2_th_fen_dds1,contrast = ment)
  nss_2_th_fen_res1 <- data.frame(nss_2_th_fen_res, stringsAsFactors = FALSE, check.names = FALSE)
  print(colnames(nss_2_th_fen_res1))
  nss_2_th_fen_res1$threshold <- ifelse(nss_2_th_fen_res1$pvalue>pv,"no",
                                        ifelse(abs(nss_2_th_fen_res1$log2FoldChange)<fc,"no",
                                               ifelse(nss_2_th_fen_res1$log2FoldChange>fc,"up","down")))
  
  nss_2_th_fen_res1 <- nss_2_th_fen_res1[order(nss_2_th_fen_res1$pvalue, 
                                               nss_2_th_fen_res1$log2FoldChange, 
                                               decreasing = c(FALSE, TRUE)), ]
  nss_2_th_fen_res1_up<- nss_2_th_fen_res1[which(nss_2_th_fen_res1$log2FoldChange >= fc 
                                                 & nss_2_th_fen_res1$pvalue < pv),]
  nss_2_th_fen_res1_down<- nss_2_th_fen_res1[which(nss_2_th_fen_res1$log2FoldChange <= -fc 
                                                   & nss_2_th_fen_res1$pvalue < pv),]
  nss_2_th_fen_res1_change <- rbind(nss_2_th_fen_res1_up,nss_2_th_fen_res1_down)
  rlog <- assay(rlogTransformation(nss_2_th_fen_dds)) %>% as.data.frame()
  rlog$symbol <- rownames(rlog) 
  select_all_sam_up <- nss_2_th_fen_res1_up[1:15,]
  select_all_sam_up$symbol <- rownames(select_all_sam_up)
  select_all_sam_up_sym <- as.data.frame(select_all_sam_up$symbol)
  colnames(select_all_sam_up_sym) <- "symbol"
  toheat_up <- merge(select_all_sam_up_sym,rlog,by="symbol")
  select_all_sam_down <- nss_2_th_fen_res1_down[1:15,]
  select_all_sam_down$symbol <- rownames(select_all_sam_down)
  select_all_sam_down_symbol <- as.data.frame(select_all_sam_down$symbol)
  colnames(select_all_sam_down_symbol) <- "symbol"
  toheat_down <- merge(select_all_sam_down_symbol,rlog,by="symbol")
  toheat <- rbind(toheat_up,toheat_down)
  p_heat <- heat_map(toheat,FALSE,FALSE,TRUE)
  all_nss_2_th_fen_res <- list(nss_2_th_fen_res1,
                               nss_2_th_fen_res1_up,
                               nss_2_th_fen_res1_down,nss_2_th_fen_res1_change)
  names(all_nss_2_th_fen_res) <- c("all","up","down","change")
  Dat <- all_nss_2_th_fen_res$all
  Dat <- na.omit(Dat)
  Dat$symbol <- rownames(Dat)
  Dat[1:2,]
  p <- ggplot(Dat,aes(x=log2FoldChange,y=-log10(pvalue),color=threshold))+
    geom_point()+
    scale_color_manual(values=c("#00008B","#808080","#DC143C"),
                       labels=c(paste0("Down:",nrow(Dat[Dat$threshold=="down",])),
                                paste0("No:",nrow(Dat[Dat$threshold=="no",])),
                                paste0("Up:",nrow(Dat[Dat$threshold=="up",]))))+#确定点的颜色
   geom_text_repel(
       data = Dat[1:20,],
       aes(label = symbol))+#添加关注的点的基因名
    theme_bw()+#修改图片背景
    theme(
      legend.title = element_blank(),#不显示图例标题
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank()
      #不显示图例标题
    )+
    ylab('-log10 ()')+#修改y轴名称
    xlab('log2 (FoldChange)') + #修改x轴名称
    geom_vline(xintercept=c(-fc,fc),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange
    geom_hline(yintercept = -log10(pv),lty=3,col="black",lwd=0.5) 
  p1 <- ggplot(Dat,aes(x=log2FoldChange,y=-log10(pvalue),color=threshold))+
    geom_point()+
    scale_color_manual(values=c("#00008B","#808080","#DC143C"),
                       labels=c(paste0("Down:",nrow(Dat[Dat$threshold=="down",])),
                                paste0("No:",nrow(Dat[Dat$threshold=="no",])),
                                paste0("Up:",nrow(Dat[Dat$threshold=="up",]))))+#确定点的颜色
    # geom_text_repel(
    #   data = Dat[1:20,],
    #   aes(label = symbol))+#添加关注的点的基因名
    theme_bw()+#修改图片背景
    theme(
      legend.title = element_blank(),#不显示图例标题
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank()
      #不显示图例标题
    )+
    ylab('-log10 ()')+#修改y轴名称
    xlab('log2 (FoldChange)') + #修改x轴名称
    geom_vline(xintercept=c(-fc,fc),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange
    geom_hline(yintercept = -log10(pv),lty=3,col="black",lwd=0.5) 
  tmp <- list(p,p1,all_nss_2_th_fen_res,p_heat)
  names(tmp) <- c("valco_plot_mark","valco_plot","res_df","heat_map")
  return(tmp)
}
cager_different_analysis_pad <- function(nss_2_licage_fen,group1,pv,fc){
  #以licage_agg为起始,对licage进行差异分析,返回火山图和差异矩阵
  nss_2_licage_fen$group <- factor(group1)
  nss_2_th_fen_dds <- consensusClustersDESeq2(nss_2_licage_fen,~group)
  nss_2_th_fen_dds1 <- DESeq(nss_2_th_fen_dds, fitType = 'mean', 
                             minReplicatesForReplace = 7, 
                             parallel = FALSE) 
  ment <- c("group",group1[length(group1)],group1[1])
  nss_2_th_fen_res <- DESeq2::results(nss_2_th_fen_dds1,contrast = ment)
  nss_2_th_fen_res1 <- data.frame(nss_2_th_fen_res, stringsAsFactors = FALSE, check.names = FALSE)
  nss_2_th_fen_res1 <- na.omit(nss_2_th_fen_res1)
  nss_2_th_fen_res1$threshold <- ifelse(nss_2_th_fen_res1$padj>pv,"no",
                                        ifelse(abs(nss_2_th_fen_res1$log2FoldChange)<fc,"no",
                                               ifelse(nss_2_th_fen_res1$log2FoldChange>fc,"up","down")))
  
  nss_2_th_fen_res1 <- nss_2_th_fen_res1[order(nss_2_th_fen_res1$padj, 
                                               nss_2_th_fen_res1$log2FoldChange, 
                                               decreasing = c(FALSE, TRUE)), ]
  nss_2_th_fen_res1_up<- nss_2_th_fen_res1[which(nss_2_th_fen_res1$log2FoldChange >= fc 
                                                 & nss_2_th_fen_res1$padj < pv),]
  nss_2_th_fen_res1_down<- nss_2_th_fen_res1[which(nss_2_th_fen_res1$log2FoldChange <= -fc 
                                                   & nss_2_th_fen_res1$padj < pv),]
  nss_2_th_fen_res1_change <- rbind(nss_2_th_fen_res1_up,nss_2_th_fen_res1_down)
  rlog <- assay(rlogTransformation(nss_2_th_fen_dds)) %>% as.data.frame()
  rlog$symbol <- rownames(rlog) 
  select_all_sam_up <- nss_2_th_fen_res1_up[1:15,]
  select_all_sam_up$symbol <- rownames(select_all_sam_up)
  select_all_sam_up_sym <- as.data.frame(select_all_sam_up$symbol)
  colnames(select_all_sam_up_sym) <- "symbol"
  toheat_up <- merge(select_all_sam_up_sym,rlog,by="symbol")
  select_all_sam_down <- nss_2_th_fen_res1_down[1:15,]
  select_all_sam_down$symbol <- rownames(select_all_sam_down)
  select_all_sam_down_symbol <- as.data.frame(select_all_sam_down$symbol)
  colnames(select_all_sam_down_symbol) <- "symbol"
  toheat_down <- merge(select_all_sam_down_symbol,rlog,by="symbol")
  toheat <- rbind(toheat_up,toheat_down)
  p_heat <- heat_map(toheat,FALSE,FALSE,TRUE)
  all_nss_2_th_fen_res <- list(nss_2_th_fen_res1,
                               nss_2_th_fen_res1_up,
                               nss_2_th_fen_res1_down,nss_2_th_fen_res1_change)
  names(all_nss_2_th_fen_res) <- c("all","up","down","change")
  Dat <- all_nss_2_th_fen_res$all
  Dat <- na.omit(Dat)
  Dat$symbol <- rownames(Dat)
  Dat[1,]
  p <- ggplot(Dat,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
    geom_point()+
    scale_color_manual(values=c("#00008B","#808080","#DC143C"),
                       labels=c(paste0("Down:",nrow(Dat[Dat$threshold=="down",])),
                                paste0("No:",nrow(Dat[Dat$threshold=="no",])),
                                paste0("Up:",nrow(Dat[Dat$threshold=="up",]))))+#确定点的颜色
    geom_text_repel(
      data = Dat[1:20,],
      aes(label = symbol))+#添加关注的点的基因名
    theme_bw()+#修改图片背景
    theme(
      legend.title = element_blank(),#不显示图例标题
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank()
      #不显示图例标题
    )+
    ylab('-log10 (p-adj)')+#修改y轴名称
    xlab('log2 (FoldChange)') + #修改x轴名称
    geom_vline(xintercept=c(-fc,fc),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange
    geom_hline(yintercept = -log10(pv),lty=3,col="black",lwd=0.5) 
  p1 <- ggplot(Dat,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
    geom_point()+
    scale_color_manual(values=c("#00008B","#808080","#DC143C"),
                       labels=c(paste0("Down:",nrow(Dat[Dat$threshold=="down",])),
                                paste0("No:",nrow(Dat[Dat$threshold=="no",])),
                                paste0("Up:",nrow(Dat[Dat$threshold=="up",]))))+#确定点的颜色
    # geom_text_repel(
    #   data = Dat[1:20,],
    #   aes(label = symbol))+#添加关注的点的基因名
    theme_bw()+#修改图片背景
    theme(
      legend.title = element_blank(),#不显示图例标题
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank()
      #不显示图例标题
    )+
    ylab('-log10 (p-adj)')+#修改y轴名称
    xlab('log2 (FoldChange)') + #修改x轴名称
    geom_vline(xintercept=c(-fc,fc),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange
    geom_hline(yintercept = -log10(pv),lty=3,col="black",lwd=0.5) 
  tmp <- list(p,p1,all_nss_2_th_fen_res,p_heat)
  names(tmp) <- c("valco_plot_mark","valco_plot","res_df","heat_map")
  return(tmp)
}


cager_diff2seq <- function(Dat,up,down){
  #从cager_different_analysis获取的差异矩阵list起始,获取上下调TC的seq
  Dat_ud <- list(Dat$up,Dat$down)
  names(Dat_ud) <- c("up","down")
  seq_li <- list()
  for (i in 1:length(Dat_ud)) {
    tmp <- Dat_ud[[i]]
    tmp$symbol <- rownames(tmp)
    tmp <- na.omit(tmp)
    split_df <- strsplit(as.character(tmp$symbol),"[:-]")
    split_df <- lapply(split_df, function(x){
      x[x==""] <- "-"
      x
    })
    split_df <- lapply(split_df, function(x){
      if (length(x)==3) {
        x[4] <- x[3]
        x[3] <- x[2]
      }
      return(x)
    })
    split_df <- lapply(split_df, function(x) 
      setNames(x,c("seqnames","start","end","strand")))
    df_split <- as.data.frame(do.call(rbind,split_df))
    tmp[,c("seqnames","start","end","strand")] <- df_split
    tmp$strand[tmp$strand==""] <- "-"
    tmp$start <- as.numeric(tmp$start)
    tmp$end <- as.numeric(tmp$end)
    gran <- GRanges(seqnames = tmp$seqnames,
                    ranges = IRanges(start = tmp$start, 
                                     end = tmp$end),
                    strand = tmp$strand)
    range <- c(-up,down)
    win <- up+down
    dom_tmp <- promoters(gran,up=up,down=down)
    dom_tmp_win <- dom_tmp[width(trim(dom_tmp))==win]
    seq_li[[i]] <- getSeq(BS,dom_tmp_win)
  }
  names(seq_li) <- c("up","down")
  return(seq_li)
}
cager_diffdf2seq <- function(Dat_ud,up,down){
  #从cager_different_analysis获取的差异矩阵res_dif集合,获取上下调TC的seq
  seq_li <- list()
  for (i in 1:length(Dat_ud)) {
    tmp <- Dat_ud[[i]]
    tmp$symbol <- rownames(tmp)
    tmp <- na.omit(tmp)
    split_df <- strsplit(as.character(tmp$symbol),"[:-]")
    split_df <- lapply(split_df, function(x){
      x[x==""] <- "-"
      x
    })
    split_df <- lapply(split_df, function(x){
      if (length(x)==3) {
        x[4] <- x[3]
        x[3] <- x[2]
      }
      return(x)
    })
    split_df <- lapply(split_df, function(x) 
      setNames(x,c("seqnames","start","end","strand")))
    df_split <- as.data.frame(do.call(rbind,split_df))
    tmp[,c("seqnames","start","end","strand")] <- df_split
    tmp$strand[tmp$strand==""] <- "-"
    tmp$start <- as.numeric(tmp$start)
    tmp$end <- as.numeric(tmp$end)
    gran <- GRanges(seqnames = tmp$seqnames,
                    ranges = IRanges(start = tmp$start, 
                                     end = tmp$end),
                    strand = tmp$strand)
    range <- c(-up,down)
    win <- up+down
    dom_tmp <- promoters(gran,up=up,down=down)
    dom_tmp_win <- dom_tmp[width(trim(dom_tmp))==win]
    seq_li[[i]] <- getSeq(BS,dom_tmp_win)
  }
  names(seq_li) <- names(Dat_ud)
  return(seq_li)
}
cager_diff2seq_na <- function(Dat,up,down){
  #从cager_different_analysis获取的差异矩阵list起始,获取上下调TC的seq
  Dat_na <- Dat$all[which(Dat$all$threshold=="no"),]
  Dat_na$symbol <- rownames(Dat_na)
  Dat_na <- na.omit(Dat_na)
  split_df <- strsplit(as.character(Dat_na$symbol),"[:-]")
  split_df <- lapply(split_df, function(x){
      x[x==""] <- "-"
      x
  })
  split_df <- lapply(split_df, function(x){
    if (length(x)==3) {
      x[4] <- x[3]
      x[3] <- x[2]
    }
    return(x)
  })
  split_df <- lapply(split_df, function(x) 
    setNames(x,c("seqnames","start","end","strand")))
  df_split <- as.data.frame(do.call(rbind,split_df))
  Dat_na[,c("seqnames","start","end","strand")] <- df_split
  Dat_na$start <- as.numeric(Dat_na$start)
  Dat_na$end <- as.numeric(Dat_na$end)
  gran <- GRanges(seqnames = Dat_na$seqnames,
                  ranges = IRanges(start = Dat_na$start, 
                                     end = Dat_na$end),
                  strand = Dat_na$strand)
  range <- c(-up,down)
  win <- up+down
  dom_tmp <- promoters(gran,up=up,down=down)
  dom_tmp_win <- dom_tmp[width(trim(dom_tmp))==win]
  seq_li <- getSeq(BS,dom_tmp_win)
  return(seq_li)
}
cager_diff2pca <- function(licage_agg,group1){
  #function:从aggre/cumu/quali后的licage开始,获取vlog和rlog之后的count进行pca
  licage_agg$group <- group1 
  licage_agg_dds <- consensusClustersDESeq2(licage_agg,~group)
  licage_agg_dds1 <- DESeq(licage_agg_dds, fitType = 'mean', 
                           minReplicatesForReplace = 7,parallel = FALSE) 
  licage_vst <- vst(licage_agg_dds1)
  #licage_rlog <- rlog(licage_agg_dds1)
  licage_raw <- SummarizedExperiment(counts(licage_agg_dds, normalized=FALSE),
                                     colData=colData(licage_agg_dds))
  p1 <- plotPCA(DESeqTransform(licage_raw), intgroup="sampleLabels")
  p2 <- plotPCA(DESeqTransform(licage_raw), intgroup="group")
  p3 <- plotPCA(licage_vst,intgroup="sampleLabels")
  p4 <- plotPCA(licage_vst,intgroup="group")
  #p5 <- plotPCA(licage_rlog,intgroup="sampleLabels")
  #p6 <- plotPCA(licage_rlog,intgroup="group")
  all_p <- list(p1,p2,p3,p4)
  names(all_p) <- c("raw_samlab","raw_group","vst_samlab","vst_group")
  return(all_p)
}

#CC
get_CC_lis <- function(licage_agg,mode){
  #mode is "feb" or "all"
  library(dplyr)
  # Licage_exp <- aggregateTagClusters(Licage_exp,tpmThreshold = 1, qLow = 0.1, qUp = 0.9, maxDist = maxd)
  all_name <- sampleLabels(licage_agg)
  # Licage_exp <- cumulativeCTSSdistribution(Licage_exp, clusters = "consensusClusters", useMulticore = TRUE)
  # Licage_exp <- quantilePositions(Licage_exp, clusters = "consensusClusters", qLow = 0.1, qUp = 0.9, useMulticore = TRUE)
  CC <- list()
  for (i in 1:length(all_name)) {
    CC[[i]] <- consensusClustersGR(licage_agg,sample=all_name[i],returnInterquantileWidth = T,
                                   qLow = 0.1, qUp = 0.9)
  }
  names(CC) <- all_name
  CC_all <- consensusClustersGR(licage_agg,returnInterquantileWidth = T,qLow = 0.1, qUp = 0.9)
  if (mode=="fen") {
    re <- CC
  }else if (mode=="all"){
    re <- CC_all
  }
  return(re)
}
get_tagcluster_lis <- function(Licage_exp){
  all_name <- sampleLabels(Licage_exp)
  tagclu <- list()
  for (i in 1:length(all_name)) {
    tagclu[[i]] <-  as.data.frame(Licage_exp@metadata$tagClusters[[i]])
  }
  names(tagclu) <- all_name
  return(tagclu)
}
get_CC_count <- function(CC_cager,tagcluster){
  library(bedtoolsr)
  options(bedtools.path = "/data1/shenluemou/biosoft/anaconda3/bin")
  #将consensusClustersGR获取的df和tagcluster进行合并,获取CC的count
  CC_cager_sym <- dplyr::select(CC_cager,c("seqnames","start","end","consensus.cluster"))
  tagcluster_sym <- dplyr::select(tagcluster,c("seqnames","start","end","nr_ctss"))
  CC_nr <- bedtoolsr::bt.intersect(CC_cager_sym,tagcluster_sym,wa = T,wb = T)
  CC_nr <- dplyr::select(CC_nr,"V1","V2","V3","V4","V8")
  CC_nr_mer <- aggregate(CC_nr[5], CC_nr[-5], FUN = function(x) sum(x))
  colnames(CC_nr_mer) <- c("seqnames","start","end","consensus.cluster","nr_ctss")
  CC_nr_mer <- dplyr::select(CC_nr_mer,"consensus.cluster","nr_ctss")
  res <- merge(CC_cager,CC_nr_mer,by="consensus.cluster",all.x=T)
  res[is.na(res)] <- 0
  return(res)
}
get_CC_count_lis <- function(CC_lis,tag_lis){
  CC_lis_df <- list()
  for (i in 1:length(CC_lis)) {
    CC_lis_df[[i]] <- as.data.frame(CC_lis[[i]])
  }
  names(CC_lis_df) <- names(CC_lis)
  tag_lis_df <- list()
  for (i in 1:length(tag_lis)) {
    tag_lis_df[[i]] <- as.data.frame(tag_lis[[i]])
  }
  names(tag_lis_df) <- names(tag_lis)
  licage_cc_count_lis <- list()
  for (i in 1:length(tag_lis)) {
    licage_cc_count_lis[[i]] <- get_CC_count(CC_lis_df[[i]],tag_lis_df[[i]])
  }
  names(licage_cc_count_lis) <- names(tag_lis)
  return(licage_cc_count_lis)
}
get_C_count_only <- function(CC_count_lis){
  all_name <- names(CC_count_lis)
  sym <- list()
  for (i in 1:length(all_name)) {
    sym[[i]] <- select(CC_count_lis[[i]],"consensus.cluster","nr_ctss")
  }
  for (i in 1:length(all_name)-1) {
    if (i==1) {
      C_count <- merge(sym[[1]],sym[[2]],by="consensus.cluster")
    }else if (i>1) {
      C_count <- merge(C_count,sym[[i+1]],by="consensus.cluster")
    }
  }
  colnames(C_count) <- c("consensus.cluster",all_name)
  return(C_count)
}
CC_to_tpmfc_choose <- function(consen_lis,minTPM,timeTPM){
  #该功能是根据TPM的阈值和fc获取各时期的特异CC
  #输入的consen_lis为consensusClustersGR获取的df所组成的list,需要names
  all_name <- names(consen_lis)
  con_tpm <- data.frame(consen_lis[[1]]$consensus.cluster)
  for (i in 1:length(consen_lis)) {
    con_tpm <- data.frame(con_tpm,consen_lis[[i]]$tpm)
  }
  colnames(con_tpm) <- c("consensus.cluster",all_name)
  con_yuan <- con_tpm
  con <- list()
  con_sym <- list()
  for (i in 2:ncol(con_tpm)) {
    tmp <- data.frame(con_tpm[,i])
    con_tpm[,i] <- con_tpm[,2]
    con_tpm[,2] <- tmp
    mark <- i-1
    con[[mark]] <- con_tpm
    con[[mark]] <- subset(con[[mark]],con[[mark]][,2]>minTPM)
    for (t in 3:ncol(con_tpm)) {
      con[[mark]] <- subset(con[[mark]],con[[mark]][,2]>timeTPM*con[[mark]][,t])
    }
    con_sym[[mark]] <- data.frame(con[[mark]]$consensus.cluster)
    colnames(con_sym[[mark]]) <- "consensus.cluster"
  }
  names(con_sym) <- all_name
  con_re <- list()
  for (i in 1:length(con_sym)) {
    con_re[[i]] <- merge(con_sym[[i]],consen_lis[[i]],by="consensus.cluster")
  }
  names(con_re) <- all_name
  return(con_re)
}
CC_to_tpmbase_choose <- function(consen_lis,maxTPM,minTPM){
  all_name <- names(consen_lis)
  con_tpm <- data.frame(consen_lis[[1]]$consensus.cluster)
  for (i in 1:length(consen_lis)) {
    con_tpm <- data.frame(con_tpm,consen_lis[[i]]$tpm)
  }
  colnames(con_tpm) <- c("consensus.cluster",all_name)
  con_yuan <- con_tpm
  con <- list()
  con_sym <- list()
  for (i in 2:ncol(con_tpm)) {
    tmp <- data.frame(con_tpm[,i])
    con_tpm[,i] <- con_tpm[,2]
    con_tpm[,2] <- tmp
    mark <- i-1
    con[[mark]] <- con_tpm
    con[[mark]] <- subset(con[[mark]],con[[mark]][,2]>maxTPM)
    for (t in 3:ncol(con_tpm)) {
      con[[mark]] <- subset(con[[mark]],con[[mark]][,t]<minTPM)
    }
    con_sym[[mark]] <- data.frame(con[[mark]]$consensus.cluster)
    colnames(con_sym[[mark]]) <- "consensus.cluster"
  }
  names(con_sym) <- all_name
  con_re <- list()
  for (i in 1:length(con_sym)) {
    con_re[[i]] <- merge(con_sym[[i]],consen_lis[[i]],by="consensus.cluster")
  }
  names(con_re) <- all_name
  return(con_re)
}
CC_to_tpm_fcbe_choose <- function(consen_lis,minTPM,timeTPM){
  all_name <- names(consen_lis)
  con_tpm <- data.frame(consen_lis[[1]]$consensus.cluster)
  for (i in 1:length(consen_lis)) {
    con_tpm <- data.frame(con_tpm,consen_lis[[i]]$tpm)
  }
  colnames(con_tpm) <- c("consensus.cluster",all_name)
  con <- list()
  con_sym <- list()
  for (i in 2:ncol(con_tpm)) {
    if (i==2) {
      con[[i-1]] <- subset(con_tpm,con_tpm[,2]>minTPM)
    }else if (i>2) {
      con[[i-1]] <- subset(con_tpm,con_tpm[,i]>minTPM)
      con[[i-1]] <- subset(con[[i-1]],con[[i-1]][,i]>timeTPM*con[[i-1]][,i-1])
    }
    con_sym[[i-1]] <- data.frame(con[[i-1]]$consensus.cluster)
    colnames(con_sym[[i-1]]) <- "consensus.cluster"
  }
  names(con_sym) <- all_name
  con_re <- list()
  for (i in 1:length(con_sym)) {
    con_re[[i]] <- merge(con_sym[[i]],consen_lis[[i]],by="consensus.cluster")
  }
  names(con_re) <- all_name
  return(con_re)
}
CC_to_sharp_broad <- function(CC){
  sharp <- CC[which(CC$interquantile_width<10),]
  broad <- CC[which(CC$interquantile_width>=10),]
  lis <- list(sharp,broad)
  names(lis) <- c("sharp","broad")
  return(lis)
}

#Grange_scan
Grange_2_scan <- function(Grange,up,down,motif,percent){
  range <- c(-up, down)
  win <- up + down
  dom_tmp <- promoters(Grange,up=up,down=down)
  dom_tmp_win <- dom_tmp[width(trim(dom_tmp))==win]
  seq_li <- getSeq(BS,dom_tmp_win)
  p1 <- plotMotifOccurrenceAverage(seq_li,motif,flankUp = up,flankDown = down,
                                   smoothingWindow = 3,color = "red",
                                   minScore = percent)
  return(p1)
}
cluter2plot_pattern <- function(up,down,cluster_sp,pattern){
  cluster_sp_range <- as.data.frame(cluster_sp@ranges)
  cluster_sp_range <-cluster_sp_range %>%
    separate(names,into = c("chr","ranges","strand"), sep = ":")
  tmp <- GRanges(seqnames = cluster_sp_range$chr, 
                 ranges = IRanges(start = cluster_sp_range$start, end = cluster_sp_range$end),
                 strand = cluster_sp_range$strand,
                 interquantile_width = cluster_sp_range$width,
                 seqlengths = seqlengths(BS))
  range <- c(-up, down)
  win <- up + down
  dom_tmp <- promoters(tmp,up=up,down=down)
  dom_tmp_win <- dom_tmp[width(trim(dom_tmp))==win]
  seq_li <- getSeq(BS,dom_tmp_win)
  p1 <- plotPatternOccurrenceAverage(seq_li,pattern,flankUp = up,flankDown = down,
                                     smoothingWindow = 3,color = "blue")
  return(list(p1))
}
cluter2plot <- function(up,down,cluster_sp,motif,percent){
  cluster_sp_range <- as.data.frame(cluster_sp)
  # cluster_sp_range <- as.data.frame(cluster_sp@ranges)
  # cluster_sp_range <-cluster_sp_range %>%
  #   separate(names,into = c("chr","ranges","strand"), sep = ":")
  tmp <- GRanges(seqnames = cluster_sp_range$seqnames, 
                 ranges = IRanges(start = cluster_sp_range$start, end = cluster_sp_range$end),
                 strand = cluster_sp_range$strand,
                 interquantile_width = cluster_sp_range$width,
                 seqlengths = seqlengths(BS))
  range <- c(-up, down)
  win <- up + down
  dom_tmp <- promoters(tmp,up=up,down=down)
  dom_tmp_win <- dom_tmp[width(trim(dom_tmp))==win]
  seq_li <- getSeq(BS,dom_tmp_win)
  p1 <- plotMotifOccurrenceAverage(seq_li,motif,flankUp = up,flankDown = down,
                                   smoothingWindow = 3,color = "red",
                                   minScore = percent)
  return(p1)
}
GRange2plot_dou <- function(GRange,up,down,motif1,motif1_name,percent1,color1,
                            motif2,motif2_name,percent2,color2){
  cluster_sp_range <- as.data.frame(GRange)
  tmp <- GRanges(seqnames = cluster_sp_range$seqnames, 
                 ranges = IRanges(start = cluster_sp_range$start, end = cluster_sp_range$end),
                 strand = cluster_sp_range$strand,
                 interquantile_width = cluster_sp_range$width,
                 seqlengths = seqlengths(BS))
  range <- c(-up, down)
  win <- up + down
  dom_tmp <- promoters(tmp,up=up,down=down)
  dom_tmp_win <- dom_tmp[width(trim(dom_tmp))==win]
  seq_li <- getSeq(BS,dom_tmp_win)
  plotMotifOccurrenceAverage(seq_li,motif1,flankUp = up,flankDown = down,
                             smoothingWindow = 3,color = color1,
                             minScore = percent1)
  plotMotifOccurrenceAverage(seq_li,motif2,flankUp = up,flankDown = down,
                             smoothingWindow = 3,color = color2,
                             minScore = percent2,add=TRUE)
  legend("topright",legend = c(paste0(motif1_name,"_",percent1),
                               paste0(motif2_name,"_",percent2)),
         col = c(color1,color2),bty = "n",lwd = 1)
}
GRange2width <- function(tag_cluster_licage,names){
  library(ggplot2)
  iq_all.l <- data.frame(iq_width = tag_cluster_licage$interquantile_width)
  iq_all.l <- list(iq_all.l)
  names(iq_all.l) <- names
  for (i in 1:length(names)) {
    iq_all.l[[i]]$sample <- rep(names[[i]], nrow(iq_all.l[[i]]))
  }
  iq_all.df <- do.call(rbind, iq_all.l)
  iq_all.df$sample <- factor(iq_all.df$sample, levels = names)
  # plotting
  p <- ggplot(iq_all.df, aes(iq_width)) +
    geom_histogram(aes(x = iq_width, y = ..ncount..* 100), binwidth = 1, 
                   fill = "gray60", col = "black", size = 0.1) +
    theme_bw() +
    theme(text = element_text(size = 16, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black", size = 14),
          axis.text.y = element_text(colour = "black", size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(colour = "black", fill = "gray87")) +
    labs(x = "Interquantile width", y = "percent") +
    xlim(0, 100)
  p1 <- p + facet_wrap(~ sample, ncol = 5)
  df <- data.frame(sample1=c(1,2),sample2=c(1,2))
  for (i in 1:length(names)) {
    tmp <- iq_all.l[[i]]
    sharp <- nrow(tmp[which(tmp$iq_width<10),])
    broad <- nrow(tmp[which(tmp$iq_width>=10),])
    df[1,i] <- sharp
    df[2,i] <- broad
  }
  colnames(df) <- names
  rownames(df) <- c("sharp","broad")
  all_re <- list(p1,df)
  return(all_re)
}
GRange_muti_motif_scan <- function(GRange_lis,up,down,motif,percent){
  all_color <- c("red","blue","black","green","purple","yellow","gray")
  all_name <- names(GRange_lis)
  all_legend <- c()
  all_col <- c()
  for (i in 1:length(GRange_lis)) {
    gr_df <- as.data.frame(GRange_lis[[i]])
    tmp <- GRanges(seqnames = gr_df$seqnames, 
                   ranges = IRanges(start = gr_df$start, end = gr_df$end),
                   strand = gr_df$strand,
                   interquantile_width = gr_df$width,
                   seqlengths = seqlengths(BS))
    range <- c(-up, down)
    win <- up + down
    dom_tmp <- promoters(tmp,up=up,down=down)
    dom_tmp_win <- dom_tmp[width(trim(dom_tmp))==win]
    seq_li <- getSeq(BS,dom_tmp_win)
    if (i==1) {
      plotMotifOccurrenceAverage(seq_li,motif,flankUp = up,flankDown = down,
                                   smoothingWindow = 3,minScore = percent,color = all_color[i],cex.axis = 0.9)
    }else if (i>1) {
      plotMotifOccurrenceAverage(seq_li,motif,flankUp = up,flankDown = down,
                                   smoothingWindow = 3,minScore = percent,color = all_color[i],cex.axis = 0.9,
                                   add=TRUE)
    }
    all_legend <- c(all_legend,paste0(all_name[i]))
    all_col <- c(all_col,all_color[i])
  }
  legend("topright",legend = all_legend,col = all_col,bty = "n",lwd = 1)
}
GRange_muti_pattern_scan <- function(GRange_lis,up,down,pattern){
  all_color <- c("red","blue","green","black","purple","yellow","gray")
  all_name <- names(GRange_lis)
  all_legend <- c()
  all_col <- c()
  for (i in 1:length(GRange_lis)) {
    gr_df <- as.data.frame(GRange_lis[[i]])
    tmp <- GRanges(seqnames = gr_df$seqnames, 
                   ranges = IRanges(start = gr_df$start, end = gr_df$end),
                   strand = gr_df$strand,
                   interquantile_width = gr_df$width,
                   seqlengths = seqlengths(BS))
    range <- c(-up, down)
    win <- up + down
    dom_tmp <- promoters(tmp,up=up,down=down)
    dom_tmp_win <- dom_tmp[width(trim(dom_tmp))==win]
    seq_li <- getSeq(BS,dom_tmp_win)
    if (i==1) {
      plotPatternOccurrenceAverage(seq_li,pattern,flankUp = up,flankDown = down,
                                   smoothingWindow = 10,color = all_color[i])
    }else if (i>1) {
      plotPatternOccurrenceAverage(seq_li,pattern,flankUp = up,flankDown = down,
                                   smoothingWindow = 10,color = all_color[i],add=TRUE)
    }
    all_legend <- c(all_legend,paste0(all_name[i]))
    all_col <- c(all_col,all_color[i])
  }
  legend("topright",legend = all_legend,col = all_col,bty = "n",lwd = 1)
}
centrimo_scan <- function(seq,name,meme_name){
  names(seq) <- paste0("seq",seq_along(seq))
  writeXStringSet(seq,paste0(name,".fa"),format = "fasta")
  system(paste0("/data1/shenluemou/biosoft/meme-5.5.3/src/centrimo --oc ",
                name,"_",meme_name," ",name,".fa ",meme_name,".meme"))
  system("mkdir html")
  system(paste0("mv ",name,"_",meme_name,"/centrimo.html ","html/",
                name,"_",meme_name,".html"))
}
centrimo_scan_motif_neg <- function(seq_neg,seq_pos,name_neg,name_pos,motif){
  seqlis <- list(seq_neg,seq_pos)
  names(seqlis) <- c(name_neg,name_pos)
  for (i in 1:length(seqlis)) {
    names(seqlis[[i]]) <- paste0("seq",seq_along(seqlis[[i]]))
    writeXStringSet(seqlis[[i]],paste0(names(seqlis)[i],".fa"),format = "fasta")
  }
  system(paste0("/data1/shenluemou/biosoft/meme-5.5.3/src/centrimo --oc ",
                name_pos,"_",name_neg,"_",motif," ",name_pos,".fa ",
                "--neg ",name_neg,".fa ",paste0(motif,".meme")))
  system("mkdir html")
  system(paste0("mv ",name_pos,"_",name_neg,"_",motif,"/centrimo.html ",
                "html/",name_pos,"_",name_neg,"_",motif,".html"))
}
centrimo_scan_db <- function(seq,name,db){
  names(seq) <- paste0("seq",seq_along(seq))
  writeXStringSet(seq,paste0(name,".fa"),format = "fasta")
  system(paste0("/Share/user/shenlm/soft/meme-5.5.3/src/centrimo --oc ",
                name,"_db"," ",name,".fa ",db))
  #system(paste0("mkdir html"))
  #system(paste0("mv ",name,"_db","/centrimo.html ","html/",name,"_db.html"))
}
centrimo_scan_db_neg <- function(seq_neg,seq_pos,name_neg,name_pos,db){
  seqlis <- list(seq_neg,seq_pos)
  names(seqlis) <- c(name_neg,name_pos)
  for (i in 1:length(seqlis)) {
    names(seqlis[[i]]) <- paste0("seq",seq_along(seqlis[[i]]))
    writeXStringSet(seqlis[[i]],paste0(names(seqlis)[i],".fa"),format = "fasta")
  }
  system(paste0("/data1/shenluemou/biosoft/meme-5.5.3/src/centrimo --oc ",
                name_pos,"_",name_neg,"_db"," ",name_pos,".fa ",
                "--neg ",name_neg,".fa ",db))
  system("mkdir html")
  system(paste0("mv ",name_pos,"_",name_neg,"_db/centrimo.html ",
                "html/",name_pos,"_",name_neg,"_db.html"))
}

#output
agg_write2be <- function(agg,pre_name){
  #input the agg and output the bw
  the_name <- sampleLabels(agg)
  for (i in 1:length(the_name)) {
    tmp_name <- the_name[i]
    tmp_trk <- exportToTrack(CTSSnormalizedTpmGR(agg,tmp_name))
    tmp_split <- split(tmp_trk,strand(tmp_trk),drop = TRUE)
    tmp_split_pos <- tmp_split$`+`
    tmp_split_neg <- tmp_split$`-`
    rtracklayer::export.bw(tmp_split_pos,paste0(pre_name,"_",tmp_name,"_pos.bw"))
    rtracklayer::export.bw(tmp_split_neg,paste0(pre_name,"_",tmp_name,"_neg.bw")) 
  }
}

#----Y(R)/Y(C)-----
from_tss_tpm_to_charkdf <- function(ctss_df,tpm_df,name){
  tmp_lis <- list(ctss_df,tpm_df)
  names(tmp_lis) <- c("ctss","tpm")
  tmp_tpm <- tapply(tmp_lis$tpm, substring(tmp_lis$ctss,1,1), sum)
  tmp_tpm_freq <- tmp_tpm / sum(tmp_tpm)
  tmp_tpm_freq <- data.frame(letter = names(tmp_tpm_freq),
                             tpm = tmp_tpm_freq)
  tss_freq <- as.data.frame(table(substring(tmp_lis$ctss,1,1))/sum(table(substring(tmp_lis$ctss,1,1))))
  names(tss_freq) <- c("letter","number")
  tmp_tpm_freq$name <- name
  tss_freq$name <- name
  tmp_lis <- list(tmp_tpm_freq,tss_freq)
  names(tmp_lis) <- c("tpm","frequency")
  return(tmp_lis)
}
from_tss_tpm_to_charkdf_pos1 <- function(ctss_df,tpm_df,name){
  tmp_lis <- list(ctss_df,tpm_df)
  names(tmp_lis) <- c("ctss","tpm")
  tmp_tpm <- tapply(tmp_lis$tpm, substring(tmp_lis$ctss,2,2), sum)
  tmp_tpm_freq <- tmp_tpm / sum(tmp_tpm)
  tmp_tpm_freq <- data.frame(letter = names(tmp_tpm_freq),
                             tpm = tmp_tpm_freq)
  tss_freq <- as.data.frame(table(substring(tmp_lis$ctss,2,2))/sum(table(substring(tmp_lis$ctss,2,2))))
  names(tss_freq) <- c("letter","number")
  tmp_tpm_freq$name <- name
  tss_freq$name <- name
  tmp_lis <- list(tmp_tpm_freq,tss_freq)
  names(tmp_lis) <- c("tpm","frequency")
  return(tmp_lis)
}
df2YR_YC <- function(df,letter_num,number,name){
  A=df[which(df[,letter_num]=="A"),][number]
  G=df[which(df[,letter_num]=="G"),][number]
  C=df[which(df[,letter_num]=="C"),][number]
  R=A+G
  per=R/C
  res <- as.data.frame(c("Y(R)","Y(C)","Y(R)/Y(C)"))
  res$value <- c(R,C,per)
  res$name <- name
  colnames(res)[1] <- "description"
  return(res)
}
cluter_Yr_yc <- function(df,num,name){
  df1 <- df
  colnames(df1)[1] <- "location"
  df1 <- separate(df1, location, 
                  into = c("seqnames", "start", "end", "strand"), 
                  sep = "[:-]", remove = FALSE, convert = TRUE)
  df1$strand[df1$strand==""] <- "-"
  df_seq <- cluter2seq(1,1,df1)
  df_tpm <- df[,num]
  df_neg <- from_tss_tpm_to_charkdf(df_seq,df_tpm,name)
  df_ctss_seq <- df_seq[grep('^[CT]',substring(df_seq,1,1))]
  df_ctss_tpm <- df_tpm[grep('^[CT]',substring(df_seq,1,1))]
  df_ctss_pos <- from_tss_tpm_to_charkdf_pos1(df_ctss_seq,df_ctss_tpm,name)
  YR_YC_freq <- df2YR_YC(df_ctss_pos$frequency,1,2,name)
  YR_YC_tpm <- df2YR_YC(df_ctss_pos$tpm,1,2,name)
  res <- list(df_neg,df_ctss_pos,YR_YC_freq,YR_YC_tpm)
  names(res) <- c("neg","pos","yr_yc_freq","yr_yc_tpm")
  return(res)
}

#------plot-------------------
plot_freq_dincle <- function(df){
  colnames(df) <- c("Var1","Freq")
  max_freq <- max(df$Freq)
  upper_limit <- ceiling(max_freq / 5) * 5
  p <- ggplot(data = df, aes(x = reorder(Var1, Freq),
                             y = Freq,
                             fill = Freq)) + 
    geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25) + 
    coord_flip() +
    scale_y_continuous(limits = c(0, upper_limit), breaks = seq(0, upper_limit, by = 5)) +
    xlab("Initiator dinucleotide") + ylab("Percentage") + 
    ggtitle(NULL) + 
    theme_bw() +
    scale_fill_gradient(low = "#0B877D", high = "#88F9D4") +
    theme(text = element_text(size = 14, colour = "black"), 
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(fill = FALSE)
  return(p)
}
plot_heatmap_up_down_10 <- function(seq_df,name){
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
                  col = c("A" = "#0B877D", "C" = "#0071A7", "G" = "orange", "T" = "red"),
                  column_split = rep(1:2,each=10),
                  show_row_names = F,
                  column_title = list(bottom = paste0(name," Position relative to Tss (bp)")),
                  use_raster=FALSE)
  return(heat)
}
plot_dflis_2_freq_din <- function(df_lis){
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
                                           fill = samples)) + 
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
    labs(fill = NULL)
  res <- list(p,seq_count_df_tidy.l)
  names(res) <- c("picture","frequency_df")
  return(res)
}

#----plot_pipline-----
plot_heatmap_up_down_10_lis_fr_dom <- function(dom,res_dir,name_pre){
  seq_1010 <- lapply(dom, function(x)
    cluter2seq(10,10,x))
  names(seq_1010) <- names(dom)
  sample_name <- names(dom)
  plot_heatmap_1010_lis <- list()
  for (i in 1:length(seq_1010)) {
    plot_heatmap_1010_lis[[i]] <- plot_heatmap_up_down_10(seq_1010[[i]],
                                                          sample_name[i])
    png(paste0(res_dir,"/",name_pre,"_",sample_name[i],".png"))
    print(plot_heatmap_1010_lis[[i]])
    dev.off()
    pdf(paste0(res_dir,"/",name_pre,"_",sample_name[i],".pdf"))
    print(plot_heatmap_1010_lis[[i]])
    dev.off()
  }
}
plot_atcg_lis_from_dom <- function(dom,up,down,name_pre,res_dir){
  seq_lis <- lapply(dom, function(x)
    cluter2seq(up,down,x))
  sample_name <- names(dom)
  names(seq_lis) <- sample_name
  for (i in 1:length(seq_lis)) {
    png(paste0(res_dir,"/",name_pre,"_",sample_name[i],".png"))
    print(plotPatternOccurrenceAverage(regionsSeq = seq_lis[[i]],
                                       patterns = c("A", "T","C","G"), flankUp = up, flankDown = down,
                                       smoothingWindow = 1, color = c("red3","green","yellow3","blue3"),plotLegend = F))
    dev.off()
    pdf(paste0(res_dir,"/",name_pre,"_",sample_name[i],".pdf"))
    print(plotPatternOccurrenceAverage(regionsSeq = seq_lis[[i]],
                                       patterns = c("A", "T","C","G"), flankUp = up, flankDown = down,
                                       smoothingWindow = 1, color = c("red3","green","yellow3","blue3"),plotLegend = F))
    dev.off()
  }
}
plot_logo_from_dom <- function(dom,up,down,name_pre,res_dir){
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
plot_motifoccur_lis <- function(dom,color_lis,up,down,minscore,smmothwin,name_pre,motif,motif_name,order_lis,res_dir){
  seq_lis <- lapply(dom, function(x)
    cluter2seq(up,down,x))
  sample_name <- names(dom)
  pdf(paste0(res_dir,"/",name_pre,"_",motif_name,"_",minscore,".pdf"))
  print(plotMotifOccurrenceAverage(seq_lis[[order_lis[1]]],motif,smoothingWindow = smmothwin,minScore = paste0(minscore,"%"),
                                   flankUp = up,flankDown = down,color = color_lis[order_lis[1]], cex.axis = 0.9))
  for (i in 2:length(order_lis)) {
    print(plotMotifOccurrenceAverage(seq_lis[[order_lis[i]]],motif,smoothingWindow = smmothwin,minScore = paste0(minscore,"%"),
                                     flankUp = up,flankDown = down,color = color_lis[order_lis[i]], cex.axis = 0.9,add = T))
  }
  print(legend("topright", 
               legend = sample_name, 
               col = color_lis,
               bty="n", lwd=1))
  dev.off()
  png(paste0(res_dir,"/",name_pre,"_",motif_name,"_",minscore,".png"))
  print(plotMotifOccurrenceAverage(seq_lis[[order_lis[1]]],motif,smoothingWindow = smmothwin,minScore = paste0(minscore,"%"),
                                   flankUp = up,flankDown = down,color = color_lis[order_lis[1]], cex.axis = 0.9))
  for (i in 2:length(order_lis)) {
    print(plotMotifOccurrenceAverage(seq_lis[[order_lis[i]]],motif,smoothingWindow = smmothwin,minScore = paste0(minscore,"%"),
                                     flankUp = up,flankDown = down,color = color_lis[order_lis[i]], cex.axis = 0.9,add = T))
  }
  print(legend("topright", 
               legend = sample_name, 
               col = color_lis,
               bty="n", lwd=1))
  dev.off()
}



