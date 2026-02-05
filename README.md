NSS-seq
=========
Nascent Strand-Specific RNA sequencing (NSS-seq) is a robust and streamlined method for genome-wide profiling of the capped 5′ ends of nascent RNAs.
This repository contains the computational pipeline for generating transcript-quantified GTF and BAM files from FASTQ data, alongside code for downstream analysis. The downstream component implements a subset of Transcription Start Site (TSS) analysis methods.

Abstract
---------
Transcriptional regulation is a highly dynamic process in which nascent RNAs provide the most immediate readout of transcriptional activity. Precise mapping of transcription start sites (TSSs) is therefore critical for understanding promoter architecture and gene regulation, yet remains technically challenging. Here, we introduce Nascent Strand-Specific RNA sequencing (NSS-seq), a robust and streamlined method for genome-wide profiling of the capped 5′ ends of nascent RNAs. By directly capturing transcription initiation events, NSS-seq overcomes the temporal delay inherent to conventional RNA-seq and enables time-resolved interrogation of transcriptional dynamics. Applied to interferon-gamma (IFN-γ)-stimulation, NSS-seq uncovered previously unrecognized IFN-γ-responsive genes and transient transcription factor activation patterns underlying interferon-mediated tumor-suppressive functions. Together, NSS-seq provides a cost-effective and technically accessible platform for dissecting promoter-level regulatory dynamics during cellular responses.

*01_shell_pipeline
----------
Shell scripts to process FASTQ files into alignment BAMs and quantitative GTF files.

*02_R_analysis
---------
R code that carries out the core TSS analysis workflow presented in our study, utilizing the CAGEr package for critical steps such as TSS clustering analysis.
