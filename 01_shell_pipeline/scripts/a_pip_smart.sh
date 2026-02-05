#!/bin/bash
#SBATCH -J test   
#SBATCH -o job-%j.log    
#SBATCH -e job-%j.err    
#SBATCH -N 2 -n 30                

#------setting------
work_dir="/project/code/01_shell_pipeline/scripts/" #scriptes_dir
all_dir="/project/code/01_shell_pipeline/res/" #Output folder
data_dir="/project/data/" 
species="human"
sample="test_IFN_0h_1 test_IFN_0h_2 test_IFN_6h_1 test_IFN_6h_2" #sample".fq"
reference_rna_dir="/hisat2/human_hg38_rrna" #rRNA-ref
reference_dir="/hisat2/human_hg38" #genome-ref
gtf_dir="/annotation/Homo_sapiens.GRCh38.113.gtf" #gtf
bed_dir="/annotation/grch38.bed" #bed

#------running----------
source $work_dir/bulk.sh
cd $work_dir
mkdir $all_dir/result
mkdir $all_dir/result/log
eval log=$all_dir/result/log
eval pair="single"

#-----smart-cut------
mkdir $all_dir/result/1_cut
for name in ${sample};do
	smartseq3_cut $data_dir $name $all_dir/result/1_cut
	echo "the smart-cut of ${name} finished"
done

#-------fastp---------------
mkdir $all_dir/result/2_fastp
for name in ${sample};do
	fastp_control ${all_dir}/result/1_cut ${name} ${all_dir}/result/2_fastp
	echo "the fastp of ${name} finished"
done

#---------trim----------------
mkdir $all_dir/result/3_trim
for name in ${sample};do
	trim_g ${all_dir}/result/2_fastp ${name} ${all_dir}/result/3_trim 50
	echo "the trim_galore of "${name}" finished"
done

#---------delete rrna-------
mkdir $all_dir/result/4_rrna
for name in ${sample};do
	de_rna ${all_dir}/result/3_trim ${name} ${all_dir}/result/4_rrna ${reference_rna_dir}
	echo "the de_rrna of "${name}" finished"
done

#-------bw--------------------
mkdir $all_dir/result/5_map
for name in ${sample};do
	fastqtobw ${all_dir}/result/4_rrna ${name} ${reference_dir} ${all_dir}/result/5_map
	echo "the mapping of "${name}" finished"
done

#--------------rnaqc-------------
mkdir $all_dir/result/6_qc
for name in ${sample};do
	rnaqc ${all_dir}/result/5_map ${name} ${all_dir}/result/6_qc ${gtf_dir} ${bed_dir}
	echo "the bam-qc of "${name}" finished"
done

#---------------abundance----------
mkdir $all_dir/result/7_abundance
for name in ${sample};do
	estimate_abundance ${all_dir}/result/5_map ${name} ${gtf_dir} ${all_dir}/result/7_abundance
	echo "the abundance of "${name}" finished"
done
mkdir ../result/8_matrix
total_abundance ${all_dir}/result/7_abundance ${all_dir}/result/8_matrix
read_matrix ${all_dir}/result/8_matrix ${all_dir}/result/7_abundance

