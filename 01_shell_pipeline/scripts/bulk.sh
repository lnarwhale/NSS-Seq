#!/bin/bash
function help_bulk(){
	echo "smartseq3_cut ｜
	echo "---------------------------"
	echo "fastp_control | 
	echo "---------------------------"
	echo "fastqtobw | 
	echo "---------------------------"
	echo "rnaqc | 
	echo "---------------------------"
	echo "estimate_abundance |
	echo "---------------------------"
	echo "total_abundance | 
	echo "---------------------------"
	echo "read_matrix | 
	echo "---------------------------"
	echo "trim_g | 
	echo "---------------------------"
	echo "de_rna | 
	echo "---------------------------"
	echo "saturation | 
	echo "---------------------------"
}
function smartseq3_cut(){
	#--function:cut the fastq of smart-seq3
	#--input: 1.fastq_dir; 2.name; 3.out_dir 
	#--output:fastq after cut
	#--need: precut.py; smartcut.py; (eval log)
	echo "the fastq is in "$1
	echo "the name is "$2
	echo "the output_dir is "$3
	fastq_dir=$1
	name=$2
	out_dir=$3
	mkdir $out_dir/tmp
	raw=$(cat $fastq_dir/$name".fq" | wc -l | awk '{print $1}')
	let raw=$raw/4
	python precut.py -i $fastq_dir -o $out_dir/tmp -n $name
	prec=$(cat $out_dir/tmp/umi.fq | wc -l | awk '{print $1}')
	let prec=$prec/4
	seqkit rmdup --by-seq --ignore-case $out_dir/tmp/umi.fq -D $out_dir/$name"_umi.txt"> $out_dir/tmp/deumi.fq
	all=$(wc -l $out_dir/tmp/deumi.fq | awk '{print $1}')
	let all=$all/4	
	cat $out_dir/$name"_umi.txt" | awk '{print $1}' > $out_dir/$name"_umi_1.txt"
	chongfu=$(cat $out_dir/$name"_umi_1.txt" | awk '{sum+=$1}END{print sum}')
	let dan=$all-$chongfu
	uniq -c $out_dir/$name"_umi_1.txt" > $out_dir/$name"_umi_2.txt"
	echo " "$dan" 1" >> $out_dir/$name"_umi_2.txt"
	sed -i 's/ \+/ /g' $out_dir/$name"_umi_2.txt"
	sed -i 's/ /,/g' $out_dir/$name"_umi_2.txt"
	sed -i  's/^.//g' $out_dir/$name"_umi_2.txt"
	mv $out_dir/$name"_umi_2.txt" $log/$name"_umi_sum.txt"
	nodu=$(cat $out_dir/tmp/deumi.fq | wc -l | awk '{print $1}')
	let nodu=$nodu/4
	mv  $out_dir/tmp/deumi.fq  $out_dir/tmp/$name.fq
	python smartcut.py -i $out_dir/tmp -o $out_dir -n $name
	rm -rf $out_dir/tmp
	touch $log/$name"_record.txt"
	echo $name"_raw "$raw >> $log/$name"_record.txt"
	echo $name"_prec "$prec >> $log/$name"_record.txt"
	echo $name"_nodu "$nodu >> $log/$name"_record.txt"
}
#----------------------------------------------------------------------
function fastp_control(){
	#--function:fastp quality control
	#--input:1.fastq_dir 2.name 3.out_dir
	#--need:fastp
	echo "the fastp of $1/$2 begin"
	fastp \
	-i $1/$2".fq" \
	-o $3/$2".fq" \
	--html $3/$2".html" \
	--trim_poly_g --poly_g_min_len 5 \
	--trim_poly_x --poly_x_min_len 5 \
	--cut_front --cut_tail --cut_window_size 4 \
	--qualified_quality_phred 15 \
	--low_complexity_filter \
        --complexity_threshold 30 \
        --length_required 30 \
        --thread 20
}
#-------------------------------------------------------------------
function fastqtobw(){
	#--function:to get bw and mapping_picture from fastq
	#--input:1.fastq_dir; 2.name; 3.reference_dir&name(index); 4.out_dir 
	#--output:sam bam bw mapping_picture
	#--need: (eval log); (eval pair); (eval work_dir(include plot_function.R));  map_stat.R
	echo "the fastq is in "$1
	echo "the name is "$2
	echo "the reference is "$3
	echo "the output_dir is "$4
	hisat2 -p 20 -x $3 -U $1/$2".fq" -S $4/$2".sam" > $4/$2".log" 2>&1
	map_g=$(samtools view -h $4/$2".sam" | grep -E 'NH:i' | awk '{print $1}' | sort | uniq | wc -l)
	echo $2"_map_genome "$map_g >> $log/$2"_record.txt" 
	samtools view -bS $4/$2".sam" > $4/$2".bam"
	samtools sort -@ 8 $4/$2".bam" -o $4/$2"_sorted.bam" 
	samtools index $4/$2"_sorted.bam"
	rm -rf $4/$2".bam"
	rm -rf $4/$2".sam"
	bamCoverage -b $4/$2"_sorted.bam" -o $4/$2".bw"
	mkdir $4/tmp
	if [ $pair = "double" ];then
		cat $4/$2".log" | grep -E ') aligned concordantly' > $4/tmp/tmp_mapstat.txt
	fi
	if [ $pair = "single" ];then
		cat $4/$2".log" | grep -E 'aligned' > $4/tmp/tmp_mapstat.txt
	fi
	sed -i 's/[a-zA-Z]//g' $4/tmp/tmp_mapstat.txt
	sed -i 's/[\(\)\%]//g' $4/tmp/tmp_mapstat.txt
	sed -i 's/[ ][ ]*/ /g' $4/tmp/tmp_mapstat.txt
	Rscript map_stat.R -i $4/tmp/tmp_mapstat.txt -w $4/tmp -f $work_dir/"plot_function.R" -s $2 > $4/tmp/$2"_forstat.log" 2>&1
	mv $4/tmp/tmp_map.txt $4/$2"_map.txt"
	mv $4/tmp/tmp_mappie.pdf $4/$2"_mappie.pdf"
	rm -rf $4/tmp
}
#-----------------------------------------------------------------
function rnaqc(){
	#--function:to do the rna_qc for bam
	#--input:1.bam_dir; 2.name; 3.out_dir; 4.reference_dir&name(gtf);	5.reference_dir&name(bed)
	#--output:picture of distribution and picture of deletation
	#--need: (eval log); (eval work_dir(include plot_function.R)); qualimap; quali.R; read_distribution.py; rseqc_stat.R
	mkdir $3/log
	mkdir $3/tmp
	mkdir $3/read_distribution/
	echo "the bam's name is "$2
	echo "the bam is in "$1
	echo "the output_dir is "$3
	echo "gtf is "$4
	echo "bed is "$5
	./qualimap/qualimap --java-mem-size=12G rnaseq -bam $1/$2"_sorted.bam" -gtf $4 -outdir $3/tmp/ -pe -outformat PDF:HTML > $3/$2".log" 2>&1
	cat $3/tmp/rnaseq_qc_results.txt | grep -E 'exonic|intronic|intergenic' > $3/tmp/tmp_readdis.txt
	sed -i 's/[\=\(\)\%\,]//g' $3/tmp/tmp_readdis.txt
	sed -i 's/[ ][ ]*/ /g' $3/tmp/tmp_readdis.txt
	sed -i 's/[\#]//g' $3/tmp/raw_data_qualimapReport/"coverage_profile_along_genes_(total).txt"
	sed -i 's/[\#]//g' $3/tmp/raw_data_qualimapReport/"coverage_profile_along_genes_(high).txt"
	sed -i 's/[\#]//g' $3/tmp/raw_data_qualimapReport/"coverage_profile_along_genes_(low).txt"	
	Rscript quali.R -d $3/tmp/tmp_readdis.txt -t $3/tmp/raw_data_qualimapReport/"coverage_profile_along_genes_(total).txt" -h $3/tmp/raw_data_qualimapReport/"coverage_profile_along_genes_(high).txt" -l $3/tmp/raw_data_qualimapReport/"coverage_profile_along_genes_(low).txt" -w $3/tmp -f ${work_dir}/"plot_function.R" -s $2 > $3/tmp/$2"_quali.log" 2>&1
	mv $3/tmp/tmp_quali.txt  $3/read_distribution/$2"_quali.txt"
	mv $3/tmp/tmp_readdis_pie.pdf $3/read_distribution/$2"_readdis_pie.pdf"
	mv $3/tmp/tmp_readdis.pdf $3/read_distribution/$2"_readdis_line.pdf"
	mkdir $3/read_distribution/tmp/
	read_distribution.py -i $1/$2"_sorted.bam" -r $5 > $3/read_distribution/tmp/tmp_tag.txt
	cp $3/read_distribution/tmp/tmp_tag.txt $3/read_distribution/$2"_taq.txt"
	cat $3/read_distribution/$2"_taq.txt" | grep -E 'Group|Exons|Introns|TSS|TES' > $3/read_distribution/tmp/tag.txt
	sed -i 's/[ ][ ]*/ /g' $3/read_distribution/tmp/tag.txt
	sed -i  "s/[\']//g" $3/read_distribution/tmp/tag.txt
	Rscript rseqc_stat.R -i $3/read_distribution/tmp/tag.txt -w $3/read_distribution/tmp/ -f ${work_dir}/"plot_function.R" -s $2 > $3/read_distribution/tmp/$2"reddis.log" 2>&1
	mv $3/read_distribution/tmp/tmp_rseqc.txt $3/read_distribution/$2"_rseqc.txt"
	mv $3/read_distribution/tmp/tmp_readdis.pdf $3/read_distribution/$2"_readdis_bar.pdf"
	rm -rf $3/read_distribution/tmp/
}
#-----------------------------------------------------------
function estimate_abundance(){
	#--function: using stringtie to calculate the rpkm or fpkm
	#--input:1.bam_dir; 2.name; 3.reference_dir&name(gtf) 4.output_dir
	#--output: fpkm/rpkm.txt
	#--need: (eval log); (eval work_dir(include plot_function.R))
	echo "the bam is"$1"/"$2"_sorted.bam"
	echo "the gtf is "$3
	echo "the output_dir is "$4
	mkdir $4/ballgown/
	mkdir $4/ballgown/$2
	stringtie -e -B -A $4/$2 -p 8 -G $3 -o $4/ballgown/$2/$2"_str_quant.gtf" $1/$2"_sorted.bam"
	awk '{print $1,$8}' $4/$2 > $4/$2"_fpkm_tmp.txt"
	sed -i '1d' $4/$2"_fpkm_tmp.txt"
	awk '$0=$0" '$2'"' $4/$2"_fpkm_tmp.txt" > $4/$2"_fpkm.txt"
	rm -rf $4/$2"_fpkm_tmp.txt"
}
#----------------------------------------------------------
function total_abundance(){
	#--function: merge fpkm and box-plot
	#--input: 1.fpkm_dir 2.working_dir
	#--output: fpkm_all.csv ; box-plot
	echo "-------------------------------------"
	echo "the fpkm is in "$1
	echo "the working_dir/output-dir is in "$2
	cat $1/*.txt > $2/total_expression.txt
	Rscript fpkm_dis.R -i $2/total_expression.txt -w $2 -f ${work_dir}/plot_function.R > $2/fpkm_dis.log 2>&1
	echo "the total_abundance finished"
	echo "-------------------------------------"
}
#----------------------------------------------------------
function read_matrix(){
	#--funtion:generate box_plot and read_matrix
	#--input: 1.working_dir 2.ballgown_dir 
	#--output: gene_count.csv ; fpkm.csv
	#--need : eval work_dir ; split_countmatrix.R
	echo "the working_dir is "$1
	echo "the ballgown_dir is "$2
	cd $1
	cp -r $2/ballgown $1
	/Share/user/shenlm/soft/prepDE.py -i ./ballgown > ./read_matrix.log 2>&1
	cd ${work_dir}
	Rscript split_countmatrix.R -i $1/gene_count_matrix.csv -w $1 >> ./read_matrix.log 2>&1
	echo "the ballgown has finished"
	echo "---------------------------------------------------------------------------------"
}
#----------------------------------------------------------
function trim_g(){
	#--function:trim_galore
	#--input: 1.fastq_dir; 2.name; 3.output_dir(must different from fastq_dir) 4.trim_length
	#--output:fastq after trim
	#--need:(eval pair)
	echo "the fastq is "$1"/"$2".fq"
	echo "the output is "$3
	echo "the pair is "${pair}
	if [ $pair = "single" ];then
		trim_galore -q 25 --phred33 --stringency 3 --length $4 -e 0.1 --no_report_file $1/$2".fq" -o $3/tmp/ > $3/$2".log" 2>&1
		mv $3/tmp/$2"_trimmed.fq" $3/$2".fq"
		rm -rf $3/tmp/
		tri=$(cat $3/$2".fq"|wc -l|awk '{print $1}')
		let tri=$tri/4
		echo $2"_trim "$tri >> $log/$2"_record.txt"
    fi
    if [ $pair = "double" ];then
        trim_galore -q 25 --phred33 --stringency 3 --length $4 -e 0.1 --no_report_file --paired $1/$2"_R1.fq" $1/$2"_R2.fq" -o $3/tmp/ > $3/$2".log" 2>&1
		mv $3/tmp/$2"_R1_trimmed.fq" $3/$2"_R1.fq"
		mv $3/tmp/$2"_R2_trimmed.fq" $3/$2"_R2.fq"
		rm -rf $3/tmp/
    fi
}
#---------------------------------------------------------
function de_rna(){
	#--function:delete the rna 
	#--input:1.fastq_dir; 2.name; 3.out_dir(different from input); 4.rna_index_dir&name(index)
	#--output:plot of rna%; fastq afer delete rna
	#--need:(eval work_dir);map_stat.R
	hisat2 -p 10 -x $4 -U $1/$2".fq" -S $3/$2"_rna.sam" > $3/$2".log" 2>&1
	cat $3/$2"_rna.sam" | grep -v "NH:i" > $3/$2"_norna.sam" 
	samtools fastq $3/$2"_norna.sam" > $3/$2".fq"
	de_rrna=$(wc -l $3/$2".fq"| awk '{print $1}')
	let de_rrna=$de_rrna/4
	echo  $2"_map_rrna "$de_rrna >> $log/$2"_record.txt"
	rm -rf $3/$2"_norna.sam" $3/$2"_rna.sam"
	mkdir $3/tmp
	cat $3/$2".log" | grep -E 'aligned' > $3/tmp/tmp_mapstat.txt
	sed -i 's/[a-zA-Z]//g' $3/tmp/tmp_mapstat.txt
	sed -i 's/[\(\)\%]//g' $3/tmp/tmp_mapstat.txt
	sed -i 's/[ ][ ]*/ /g' $3/tmp/tmp_mapstat.txt
	Rscript map_stat.R -i $3/tmp/tmp_mapstat.txt -w $3/tmp -f $work_dir/"plot_function.R" -s $2 > $3/tmp/$2"_forstat.log" 2>&1
	mv $3/tmp/tmp_map.txt $3/$2"_map.txt"
	mv $3/tmp/tmp_mappie.pdf $3/$2"_mappie.pdf"
}



