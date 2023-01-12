#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process PREPARE_SINGLE_LINKAGE{

	label "gcc"
	publishDir "${params.folder_template}/src", mode: 'copy'
	
	input:
	path cpp
	
	output:
	val "done"
	path "single_linkage"

	script:
	"""
	g++ ${cpp} \
		-o single_linkage
	"""	
}

process PREPARE_IN_SLC {

	label "bash"
	publishDir "${params.folder_template}/result/${specific_sample}", mode: 'copy'
	
	input:
	val specific_sample
	val species
	val prepare_singleLinkage

	output:
	path "in_slc"
	path "in_slc_PWM"
	path "in_slc_motif"
	path "in_slc_num_polyA"
	// pass as "val" and then when you pass full data as input use "path". This is to avoid scope issue. 
	// (you didnt receive full data as an input here but you are trying to use it as a file)
	tuple val("${params.folder_template}/data/${species}/${specific_sample}.bam"), val("${params.folder_template}/data/${species}/${specific_sample}.bam.bai"), val("${specific_sample}")

	script:
	"""
	mkdir in_slc
	mkdir in_slc_PWM
	mkdir in_slc_motif
	mkdir in_slc_num_polyA
	"""
}

process TRIM_NEG_CONTROL {

	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase1_output", mode: 'copy'

	input:
	tuple path(bam), path(bai), val(specific_sample)
	path python_script
	
	output:
	tuple path ("${specific_sample}_fixed.bam"), val("${specific_sample}")

	script:
	"""
	python3 ${python_script}\
		--bam_input ${bam}\
		--bam_out ${specific_sample}_fixed.bam
	"""	
}

process SPLIT_PHASE1 {
	
	echo true
	label "samtools"
	publishDir "${params.folder_template}/result/${specific_sample}/phase1_output", mode: 'copy'
	tag "$chrom"

	input:
	tuple path(bam), path(bai), val(specific_sample), val(chrom)

	output:
	tuple path("possorted_*.bam"), path("possorted_*.bam.bai"), val("${specific_sample}")
	tuple val("${specific_sample}"), path("bam_folder_${chrom}")

	script:
	"""
	# echo ${bam}
	# echo ${chrom}
	samtools view -bh ${bam} chr${chrom} -o intermediate_${chrom} 

	samtools sort intermediate_${chrom} -o possorted_${chrom}.bam
		
	samtools index possorted_${chrom}.bam

	mkdir bam_folder_${chrom}
	# copy possorted_number.bam, bai to bamfolder_number folder
	cp possorted_* bam_folder_${chrom}
	
	"""         
}

process DEDUP {

	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase1_output", mode: 'copy'

	input:
	tuple val(specific_sample), path(bam_folder)
	path python_script
	
	output:
	tuple path ("deduplicated_chr*"), val("${specific_sample}")
	tuple path ("splitted_chr*"), val("${specific_sample}")

	script:
	"""
	# copy everything from bam folders (1~19, X, Y) to the current working directory.
	cp ${bam_folder}/* .
	# * is for placeholder not "all"
	python3 ${python_script} \
		--bam_template possorted_*.bam \
		--out_template ${params.dedup_out_name_prefix}\
		--span_threshold ${params.span_threshold}\
		--split_read_template ${params.split_out_name_prefix}	
	"""
}

process SORT_PHASE1{
	
	echo true
	label "samtools"
	publishDir "${params.folder_template}/result/${specific_sample}/phase1_output/", mode: 'copy'
	
	input:
	tuple path(bams), val(specific_sample)

	output:
	tuple path("*_sorted.bam"), path("*_sorted.bam.bai"), val("${specific_sample}")

	script:
	"""	
	input=\$(basename ${bams})
	prefix=\$(echo \$input | cut -d '.' -f 1)

	samtools sort ${bams} -o \${prefix}_sorted.bam
	samtools index \${prefix}_sorted.bam
	"""
}

process CONVERT_TUPLE{
	label "bash"
	input:
	tuple path(bams), path(bais), val(specific_sample)
	
	output:
	tuple path("${bams}"), val("${specific_sample}")
	tuple path("${bais}"), val("${specific_sample}")

	script:
	"""
	"""
}

process SPLIT_TUPLE{
	label "bash"
	input:
	tuple path(bams), path(bais), val(specific_sample)
	
	output:
	path("${bams}")
	path("${bais}")
	val("${specific_sample}")

	script:
	"""
	"""
}

process MERGE{

	label "samtools"
	publishDir "${params.folder_template}/result/${specific_sample}/phase1_output", mode: 'copy'
	
	// input channel looks like: [([bam1, bam2...... bam6.....bamY], sample1), ([bam1, bam2...... bam6.....bamY], sample2)]
	// this is to merge once and only once rather than merging several times
	input:
	tuple path(bams), val(specific_sample) 
	tuple path(bais), val(specific_sample)

	output:
	tuple path ("dedup_full_sorted.bam"), path ("dedup_full_sorted.bam.bai"), val("${specific_sample}")
		
	script:
	"""
	samtools merge -f -o dedup_full.bam deduplicated_chr*_output_sorted.bam
	samtools sort dedup_full.bam -o dedup_full_sorted.bam
	samtools index dedup_full_sorted.bam
	"""
}

process GET_SPAN_RAW{

	label "custom_python"
	label "heavy_computation"
	publishDir "${params.folder_template}/result/${specific_sample}/phase1_output/", mode: 'copy'

	input:
	tuple path(raw_bams), path(raw_bais), val(specific_sample)
	val raw
	val sample_or_control
	path python_script

	output:
	tuple path("${params.raw_span}_*.csv"), val("${specific_sample}")
	
	script:
	"""
	if [[ ${sample_or_control} == "sample" ]]
	then
		raw_input=\$(basename ${raw_bams})
		raw_prefix=\$(echo \$raw_input | cut -d '_' -f 2)
		raw_number=\$(echo \$raw_prefix | cut -d '.' -f 1)

		python3 ${python_script}\
		--input ${raw_bams}\
		--output ${params.raw_span}_\${raw_number}.csv\
		--action_type ${raw}

	elif [[ ${sample_or_control} == "control" ]]
	then
		raw_input=\$(basename ${raw_bams})
		raw_prefix=\$(echo \$raw_input | cut -d '_' -f 8)
		raw_number=\$(echo \$raw_prefix | cut -d '.' -f 1)

		python3 ${python_script}\
		--input ${raw_bams}\
		--output ${params.raw_span}_\${raw_number}.csv\
		--action_type ${raw}
	fi	
	"""
}

process GET_SPAN_DEDUP{

	label "custom_python"
	label "heavy_computation"
	publishDir "${params.folder_template}/result/${specific_sample}/phase1_output/dedup_span", mode: 'copy'

	input:
	tuple path(dedup_bams), path(dedup_bais), val(specific_sample)
	val dedup
	val sample_or_control
	path python_script

	output:
	tuple path("${params.dedup_span}_*.csv"), val("${specific_sample}")
	
	script:
	"""
	if [[ ${sample_or_control} == "sample" ]]
	then
		dedup_input=\$(basename ${dedup_bams})
		dedup_prefix=\$(echo \$dedup_input | cut -d '_' -f 2)
		dedup_number=\$(echo \$dedup_prefix | cut -d 'r' -f 2)

		python3 ${python_script}\
		--input ${dedup_bams}\
		--output ${params.dedup_span}_\${dedup_number}.csv\
		--action_type ${dedup}

	elif [[ ${sample_or_control} == "control" ]]
	then
		dedup_input=\$(basename ${dedup_bams})
		dedup_prefix=\$(echo \$dedup_input | cut -d '_' -f 8)
		dedup_number=\$(echo \$dedup_prefix | cut -d '.' -f 1)

		python3 ${python_script}\
		--input ${dedup_bams}\
		--output ${params.dedup_span}_\${dedup_number}.csv\
		--action_type ${dedup}	

	fi		
	"""
}

process GET_SPAN_SPLIT{

	label "custom_python"
	label "heavy_computation"
	publishDir "${params.folder_template}/result/${specific_sample}/phase1_output/split_span", mode: 'copy'

	input:
	tuple path(split_bams), path(split_bais), val(specific_sample)
	val split
	path python_script

	output:
	tuple path("${params.split_span}_*.csv"), val("${specific_sample}")

	script:
	"""
	split_input=\$(basename ${split_bams})
	split_prefix=\$(echo \$split_input | cut -d '_' -f 2)
	split_number=\$(echo \$split_prefix | cut -d 'r' -f 2)

	python3 ${python_script}\
	--input splitted_chr\${split_number}_output_sorted.bam\
	--output ${params.split_span}_\${split_number}.csv\
	--action_type ${split}
	"""
}

process CONCAT_SPAN{

	label "bash"
	publishDir "${params.folder_template}/result/${specific_sample}/phase1_output", mode: 'copy'

	input:
	// channel looks like: [([csv1, csv2...... csv6.....csvY], sample1), ([csv1, csv2...... csv6.....csvY], sample2)]
	tuple path(csv_files), val(specific_sample)
	val(file_name)

	output:
	tuple path("${specific_sample}_combined_${file_name}.csv"), val("${specific_sample}")

	script:
	"""
	cat ${file_name}_*.csv > ${specific_sample}_combined_${file_name}.csv
	# cat ${params.dedup_span}_*.csv > combined_dedup_span.csv
	"""
}

process GET_SPAN_HIST{

	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase1_output", mode: 'copy'

	input:
	tuple path(combined_spans), val(specific_sample)
	path python_script
	val(file_name)

	output:
	path "${specific_sample}_${file_name}"

	script:
	"""
	python3 ${python_script}\
	--input_csv ${combined_spans}\
	--out_hist ${specific_sample}_${file_name}
	"""
}

process FIND_TERMINAL_EXONS{

	label "custom_python"
	publishDir "${params.folder_template}/result/common", mode: 'copy'

	input:
	path python_script

	output:
	path "${params.terminal_exons_out}"

	script:
	"""
	python3 ${python_script} --gtf_file ${params.gtf_fasta_location_template}/${params.sample_type}/${params.gtf} --bed_out ${params.terminal_exons_out}
	"""
}

process FIND_EXONS_GENES_BED{
	
	label "custom_python"
	publishDir "${params.folder_template}/result/common", mode: 'copy'

	input:
	path python_script

	output:
	path("${params.exons_bed_out}")
	path("${params.genes_bed_out}")

	script:
	"""
	python3 ${python_script}\
	--gtf_dir ${params.gtf_fasta_location_template}/${params.sample_type}/${params.gtf}\
	--genes_bed_out ${params.genes_bed_out}\
	--exons_bed_out ${params.exons_bed_out}
	"""
}

process GET_POLYA{

	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase2_output", mode: 'copy'
	
	input:
	// for each specific sample, run python script
	tuple path(dedup_full_bam), path(dedup_full_bai), val(specific_sample)
	val polyA_out
	val nonpolyA_out
	path python_script

	output:
	tuple path("${polyA_out}"), val("${specific_sample}")
	tuple path("${nonpolyA_out}"), val("${specific_sample}")
	
	script:
	"""
	python3 ${python_script} --bam_input ${dedup_full_bam}\
	--o_polyA ${polyA_out}\
	--o_nonpolyA ${nonpolyA_out}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--percentage_threshold ${params.polyA_percentage_threshold}\
	--length_threshold ${params.length_threshold}\
	--use_fc ${params.use_fc}
	"""
}

// sort on single full data for all samples 
process SORT_PHASE2{

	echo true
	label "samtools"
	publishDir "${params.folder_template}/result/${specific_sample}/phase2_output", mode: 'copy'
	
	input:
	tuple path(bam), val(specific_sample)

	output:
	tuple path("*_sorted.bam"), path("*_sorted.bam.bai"), val("${specific_sample}")

	script:
	"""	
	input=\$(basename ${bam})
	prefix=\$(echo \$input | cut -d '.' -f 1)

	samtools sort ${bam} -o \${prefix}_sorted.bam
	samtools index \${prefix}_sorted.bam
	"""
}

process SORT_PHASE2_ALTER{

	echo true
	label "samtools"
	publishDir "${params.folder_template}/result/${specific_sample}/phase2_output", mode: 'copy'
	
	input:
	tuple path(bam), val(specific_sample)

	output:
	tuple path("*_sorted.bam"), path("*_sorted.bam.bai"), val("${specific_sample}")

	script:
	"""	
	input=\$(basename ${bam})
	prefix=\$(echo \$input | cut -d '.' -f 1)

	samtools sort ${bam} -o ${specific_sample}_\${prefix}_sorted.bam
	samtools index ${specific_sample}_\${prefix}_sorted.bam
	"""
}

process INTERSECT{

	label "bedtools"
	publishDir "${params.folder_template}/result/${specific_sample}/phase2_output", mode: 'copy'

	input:
	// here it is ok to have bed file and bam, bai separate because the length are the same (i.e. number of full_bams == number of samples). 
	// but if two channels have different lengths, then you need to concatenate two channels into one channel
	// if channel of length = 1, it does not matter if two channels have different lengths or not (bed file is length of 1 so it is ok).
	tuple path(bam), path(bai), val(specific_sample)
	path(bed)
	val(out_name)

	output:
	tuple path("${out_name}.bam"), val("${specific_sample}")

	script:
	"""
	bedtools intersect -s -nonamecheck -abam ${bam} -b ${bed} > ${out_name}.bam
	"""
}

process NO_INTERSECT{
	
	label "bedtools"
	publishDir "${params.folder_template}/result/${specific_sample}/phase2_output", mode: 'copy'

	input:
	// here it is ok to have bed file and bams_bais separate because the length are the same (i.e. number of full_bams == number of samples). 
	// but if two channels have different lengths, then you need to concatenate two channels into one channel
	// if channel of length = 1, it does not matter if two channels have different lengths or not (bed file is length of 1 so it is ok).
	tuple path(bam), path(bai), val(specific_sample)
	path(bed)
	val(out_name)

	output:
	tuple path("${out_name}.bam"), val("${specific_sample}")

	script:
	"""
	bedtools intersect -s -v -nonamecheck -abam ${bam} -b ${bed} > ${out_name}.bam
	"""
}

// sort on single full data for all samples (after intersect)
process SORT_AFTER_INTERSECT{
	
	label "samtools"
	publishDir "${params.folder_template}/result/${specific_sample}/phase2_output", mode: 'copy'

	input:
	tuple path(bam), val(specific_sample)
	val(out_name)

	output:
	tuple path("${out_name}_sorted.bam"), path("${out_name}_sorted.bam.bai"), val("${specific_sample}")

	script:
	"""
	samtools sort ${bam} -o ${out_name}_sorted.bam
	samtools index ${out_name}_sorted.bam
	"""
}

process SPLIT_PHASE2 {

	echo true
	label "samtools"
	publishDir "${params.folder_template}/result/${specific_sample}/phase2_output", mode: 'copy'
	
	input:
	tuple path(bam), path(bai), val(specific_sample), val(chrom)
	val(out_name)

	output:
	tuple path("${specific_sample}_${out_name}_*.bam"), path ("${specific_sample}_${out_name}_*.bam.bai"), val("${specific_sample}")

	script:
	"""	
	samtools view -bh ${bam} chr${chrom} | samtools sort -o  ${specific_sample}_${out_name}_sorted_${chrom}.bam
	samtools index ${specific_sample}_${out_name}_sorted_${chrom}.bam
	"""         
}

process GET_PRO_INTRONIC_INTERGENIC{
	
	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase2_output", mode: 'copy'

	input:
	tuple path(polyA_bam), path(polyA_bai), val(specific_sample)
	path python_script

	output:
	tuple path("${params.out_pro_intronic_intergenic}"), val("${specific_sample}") 

	script:
	"""
	python3 ${python_script} --bam_input ${polyA_bam}\
	--bam_out ${params.out_pro_intronic_intergenic}\
	"""
}

process ISOLATE_POLYA_TERMINAL{

	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase2_output", mode: 'copy'

	input:
	// if channel of length = 1, it does not matter if two channels have different lengths or not (bed file is length of 1). so did not fuse them
	tuple path(bams), path(bais), val(specific_sample)
	path(terminal_exons)
	path python_script
	
	output:
	tuple path ("isolated_reads_above_chr*"), val("${specific_sample}")
	tuple path("isolated_reads_below_chr*"), val("${specific_sample}")

	script:
	"""
	python3 ${python_script} --bam_input ${bams}\
	--bed_dir ${terminal_exons}\
	--unannotated_template ${params.unannotated_in_template}\
	--annotated_template ${params.annotated_in_template}\
	--distance_threshold ${params.distance_threshold}\
	--use_fc ${params.use_fc}
	"""
}

process GET_COUNTS{

	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase3_output", mode: 'copy'

	input:
	tuple path(full_bams), path(full_bais), val(specific_sample)
	val(csv_name)
	val(bam_type)
	path python_script
	
	output:
	path("${specific_sample}_${csv_name}")
	
	script:
	"""
	python3 ${python_script} \
	--bam_input ${full_bams}\
	--csv_output ${specific_sample}_${csv_name}\
	--sample_name ${specific_sample}\
	--bam_type ${bam_type}
	"""
}

process CONCAT_COUNTS{
	label "bash"
	publishDir "${params.folder_template}/result/common", mode: 'copy'
	
	input:
	path (csvs)
	val(bam_type)

	output:
	path("merged_${bam_type}.csv")

	script:
	"""
	cat *.csv > merged_${bam_type}.csv
	"""	
}

process MERGE_ALL_COUNTS{
	label "bash"
	publishDir "${params.folder_template}/result/common", mode: 'copy'
	
	input:
	path(raw_concat_csv)
	path(dedup_concat_csv)
	path(nonpolyA_concat_csv)
	path(polyA_concat_csv)
	path(polyA_T_concat_csv)
	path(intronic_concat_csv)
	path(intergenic_concat_csv)
	path(exonic_concat_csv)

	output:
	path("total_merged.csv")

	script:
	"""
	cat *.csv > total_merged.csv
	"""
}

process SPLIT_PER_SAMPLE{

	label "custom_python"
	publishDir "${params.folder_template}/result/separate_merged_counts", mode: 'copy'

	input:
	path(totalCounts_csv)
	path(python_script)

	output:
	path("*X*.csv")

	script:
	"""
	python3 ${python_script}\
	--total_csv_input ${totalCounts_csv}
	"""
}

process ONE_BY_ONE{
	label "bash"
	
	input:
	each(path(split_count_csvs))

	output:
	path("${split_count_csvs}")
	
	script:
	"""
	"""	
}

process MAKE_COUNT_SAMPLE{
	label "bash"
	
	input:
	tuple path(count_csvs), val(specific_samples)

	output:
	tuple path("${count_csvs}"), val("${specific_samples}")
	
	script:
	"""
	"""	
}

process GET_BAR_CHART{

	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase3_output", mode: 'copy'
	
	input:
	tuple path(csvs), val(specific_sample)
	path python_script
	
	output:
	path "${specific_sample}_${params.barchart_output}"

	script:
	"""
	python3 ${python_script}\
	--input_dir ${csvs}\
	--out_file ${specific_sample}_${params.barchart_output}
	"""
}

process GET_DISTANCE_ALTER_SAMPLE{

	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase4_output", mode: 'copy'

	input:
	tuple path(polyA_t_full_bam), path(polyA_t_full_bai), val(specific_sample)
	path terminal_exons
	val sample_or_control
	path python_script

	output:
	tuple path("distances.csv"), val("${specific_sample}")
	path "outlier*.csv"
	
	script:
	"""
	python3 ${python_script}\
	--bam_input ${polyA_t_full_bam}\
	--bed_input ${terminal_exons}\
	--distance_out ${params.distance_output}\
	--use_fc ${params.use_fc}\
	--sample_or_control ${sample_or_control}
	"""
}

process GET_DISTANCE_ALTER_CONTROL{

	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase4_output", mode: 'copy'

	input:
	tuple path(polyA_t_full_bam), path(polyA_t_full_bai), val(specific_sample)
	path terminal_exons
	val sample_or_control
	path python_script

	output:
	tuple path("distances*.csv"), val("${specific_sample}")
	path "outlier*.csv"
	
	script:
	"""
	python3 ${python_script}\
	--bam_input ${polyA_t_full_bam}\
	--bed_input ${terminal_exons}\
	--distance_out ${params.distance_output}\
	--use_fc ${params.use_fc}\
	--sample_or_control ${sample_or_control}
	"""
}

process COMBINE_CONTROL_DISTANCES{

	label "bash"
	publishDir "${params.folder_template}/result/${specific_sample}/phase4_output", mode: 'copy'

	input:
	tuple path(csvs), val(specific_sample)
	val out_name

	output:
	path("${out_name}")

	script:
	"""
	cat *.csv > ${out_name}
	"""
}

process GET_D_HISTOGRAM_ALTER{

	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase4_output", mode: 'copy'

	input:
	tuple path(distance_csv), val(specific_sample), path(negative_control_distance_csv)
	path python_script

	output:
	path "${params.distance_histogram_file_name}"
	
	script:
	"""
	python3 ${python_script}\
	--csv_input ${distance_csv}\
	--cntrl_csv ${negative_control_distance_csv}\
	--file_name ${params.distance_histogram_file_name}\
	--use_fc ${params.use_fc}
	"""
}

process PLOT_ATGC{

	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase5_output", mode: 'copy'
	// pass all bam files at once (not one by one)
	// and then within python script it uses bam file one by one.
	input:
	tuple path(bams), path(bais), val(specific_sample)
	val bam_template
	val output_name
	val annotated
	path single_linkage
	val species
	path python_script

	output:
	path "${specific_sample}_${output_name}"

	script:
	"""
	python3 ${python_script}\
	--bam_template ${specific_sample}_${bam_template}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--slc_distance ${params.slc_distance_parameter}\
	--out_name ${specific_sample}_${output_name}\
	--annotated ${annotated}\
	--in_slc ${params.folder_template}/result/${specific_sample}/in_slc\
	--use_fc ${params.use_fc}\
	--species ${species}
	"""
}

process FIX_SOFTCLIPPED_REGION{
	
	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase2_output", mode: 'copy'

	input:
	tuple path(bam), path(bai), val(specific_sample)
	path python_script

	output:
	tuple path("${params.bam_after_fixation}"), val("${specific_sample}")
	
	script:
	"""
	python3 ${python_script}\
	--bam_file ${bam}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--bam_out ${params.bam_after_fixation}

	"""
}

process GET_MOTIF_FREQ_PLOT{

	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase6_output", mode: 'copy'

	input:
	tuple path(bams), path(bais), val(specific_sample)
	val bam_template
	val output_name
	val annotated
	path peaks_csv
	path single_linkage
	val species
	path python_script

	output:
	tuple path ("${specific_sample}_${output_name}*.png"), val("${specific_sample}")
	tuple path ("${specific_sample}_${output_name}_overlaid.png"), val("${specific_sample}")
	path ("${specific_sample}_${output_name}_ordered_motives.csv")
	path ("${specific_sample}_${output_name}_scores.csv")

	script:
	"""
	python3 ${python_script}\
	--bam_template ${specific_sample}_${bam_template}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--slc_distance ${params.slc_distance_parameter}\
	--out_name ${specific_sample}_${output_name}\
	--annotated ${annotated}\
	--motif_in_slc ${params.folder_template}/result/${specific_sample}/in_slc_motif\
	--window ${params.motif_search_window}\
	--motif_info_dir ${params.motif_info_dir}\
	--downstream ${params.down_limit}\
	--peaks ${peaks_csv}\
	--species ${species}
	"""
}

process GET_MOTIF_FREQ_GTF_PLOT{

	label "custom_python"
	publishDir "${params.folder_template}/result/common", mode: 'copy'

	input:
	path terminal_exons
	val output_name
	path single_linkage
	path python_script

	output:
	path "${output_name}_*.png"
	path "${output_name}_overlaid.png"
	path "${output_name}_ordered_*.csv"
	path "${output_name}_peaks.csv"

	script:
	"""
	python3 ${python_script}\
	--terminal_exons ${terminal_exons}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--slc_distance ${params.slc_distance_parameter}\
	--out_name ${output_name}\
	--motif_in_slc ${params.folder_template}/result/common\
	--window ${params.motif_search_window}\
	--motif_info_dir ${params.motif_info_dir}\
	--downstream ${params.down_limit}
	"""
}

process CONCAT_MOTIF_ORDERS{
	label "bash"
	publishDir "${params.folder_template}/result/common", mode: 'copy'
	
	input:
	path (motif_orders_csvs)
	path (motif_scores_csvs)
	val (bam_type)

	output:
	path("merged_${bam_type}_motif_order.csv")
	path("merged_${bam_type}_scores.csv")

	script:
	"""
	cat *_ordered_motives.csv > merged_${bam_type}_motif_order.csv
	cat *_scores.csv > merged_${bam_type}_scores.csv
	"""	
}

process MERGE_ALL_MOTIF_CSVS{
	label "bash"
	publishDir "${params.folder_template}/result/common", mode: 'copy'
	
	input:
	path(annotated_concat_csv)
	path(unannotated_concat_csv)
	path(intergenic_concat_csv)
	path(intronic_concat_csv)
	path(exonic_concat_csv)
	val(out_name)

	output:
	path("${out_name}")

	script:
	"""
	cat *.csv > ${out_name}
	"""
}

process SPLIT_MOTIF_ORDER_PER_SAMPLE{

	label "custom_python"
	publishDir "${params.folder_template}/result/separate_merged_motif_orders", mode: 'copy'

	input:
	path(total_motif_orders_csv)
	path(python_script)

	output:
	path("*10X*.csv")

	script:
	"""
	python3 ${python_script}\
	--total_csv_input ${total_motif_orders_csv}
	"""
}

process MOTIF_ORDERS_ONE_BY_ONE{
	label "bash"
	
	input:
	each(path(split_motif_orders_csvs))

	output:
	path("${split_motif_orders_csvs}")
	
	script:
	"""
	"""	
}

process COMBINE_MOTIF_ORDERS{

	label "bash"
	publishDir "${params.folder_template}/result/${specific_sample}/phase6_output", mode: 'copy'

	input:
	tuple path(csvs), val(specific_sample)
	path gtf_motif
	val out_name

	output:
	tuple path("${specific_sample}_${out_name}"), val("${specific_sample}")

	script:
	"""
	cat *.csv > ${specific_sample}_${out_name}
	"""
}

process COMPUTE_SCORES_FOR_MOTIF_ALTER{

	label "custom_python"
	label "heavy_computation"
	publishDir "${params.folder_template}/result/${specific_sample}/phase6_output", mode: 'copy'

	input:
	tuple path(bams), path(bais), val(specific_sample)
	val bam_template
	val output_name
	path peaks_csv
	path single_linkage
	val species
	path python_script

	output:
	tuple path ("${specific_sample}_${output_name}*.png"), val("${specific_sample}")
	tuple path ("${specific_sample}_${output_name}_overlaid.png"), val("${specific_sample}")
	path "${specific_sample}_${output_name}_*.csv"
	
	script:
	"""
	python3 ${python_script}\
	--bam_template ${specific_sample}_${bam_template}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--slc_distance ${params.slc_distance_parameter}\
	--out_name ${specific_sample}_${output_name}\
	--motif_in_slc ${params.folder_template}/result/${specific_sample}/in_slc_motif\
	--window ${params.motif_search_window}\
	--motif_info_dir ${params.motif_info_dir}\
	--downstream ${params.down_limit}\
	--peaks ${peaks_csv}\
	--species ${species}
	"""
}

process COMBINE_CSVS{

	label "bash"
	publishDir "${params.folder_template}/result/common", mode: 'copy'

	input:
	path csvs
	val out_name

	output:
	path "${out_name}"

	script:
	"""
	cat *.csv > ${out_name}
	"""
}

process GET_BAR_CHART_SCORE_MOTIF_POLYA{

	label "custom_python"
	publishDir "${params.folder_template}/result/common", mode: 'copy'
	
	input:
	path(csvs)
	path python_script
	
	output:
	path "${params.polyA_motif_score}"

	script:
	"""
	python3 ${python_script}\
	--input_dir ${csvs}\
	--out_file ${params.polyA_motif_score}
	"""
}

process GET_SOFTCLIPPED_DISTRIBUTION{

	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase6_output", mode: 'copy'

	input:
	tuple path(bams), path(bais), val(specific_sample)
	val output_name
	path python_script

	output:
	tuple path ("${specific_sample}_${output_name}"), val("${specific_sample}")
	
	script:
	"""
	python3 ${python_script}\
	--dedup_bam ${bams}\
	--file_name ${specific_sample}_${output_name}
	"""
}

process GET_NUM_POLYA_SITES{
	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase7_output", mode: 'copy'

	input:
	tuple path(bams), path(bais), val(specific_sample)
	val csv_output_name
	val barchart_output_name
	path single_linkage
	path python_script

	output:
	path "${specific_sample}_${csv_output_name}"
	path "${specific_sample}_${barchart_output_name}"

	script:
	"""
	python3 ${python_script}\
	--bam ${bams}\
	--slc_distance ${params.slc_distance_parameter}\
	--num_polyA_in_slc ${params.folder_template}/result/${specific_sample}/in_slc_num_polyA\
	--csv_out_name ${specific_sample}_${csv_output_name}\
	--barchart_out_name ${specific_sample}_${barchart_output_name}\
	--cluster_threshold ${params.cluster_threshold}
	"""	
}

process GET_BAR_CHART_NUM_POLYA{

	label "custom_python"
	publishDir "${params.folder_template}/result/common", mode: 'copy'
	
	input:
	path(csvs)
	val output_name
	path python_script
	
	output:
	path "${output_name}"

	script:
	"""
	python3 ${python_script}\
	--input_dir ${csvs}\
	--out_file ${output_name}
	"""
}

process GET_GENE_COVERAGE{
	label "custom_python"
	publishDir "${params.folder_template}/result/${specific_sample}/phase8_output", mode: 'copy'

	input:
	tuple path(bams), path(bais), val(specific_sample)
	val output_name
	path terminal_exons
	val type_read
	path python_script

	output:
	path "${specific_sample}_${output_name}"

	script:
	"""
	python3 ${python_script}\
	--bam ${bams}\
	--out_name ${specific_sample}_${output_name}\
	--read_count_threshold ${params.read_count_threshold}\
	--terminal_exons ${terminal_exons}\
	--type_read ${type_read}
	"""	
}

process GET_TOTAL_NUM_GENES{

	label "custom_python"
	publishDir "${params.folder_template}/result/common", mode: 'copy'
	
	input:
	path terminal_exons
	val out_name
	val species
	path python_script
	
	output:
	path "${out_name}"

	script:
	"""
	python3 ${python_script}\
	--terminal_exons ${terminal_exons}\
	--out_name ${out_name}\
	--species ${species}
	"""
}

process GET_BAR_CHART_GENE_COVERAGE{

	label "custom_python"
	publishDir "${params.folder_template}/result/common", mode: 'copy'
	
	input:
	path dedup_csv
	path polyA_csv
	path total_gene_csv
	val coverage_output_name
	path python_script
	
	output:
	path "${coverage_output_name}"
	
	script:
	"""
	python3 ${python_script}\
	--dedup_csv_dir ${dedup_csv}\
	--polyA_csv_dir ${polyA_csv}\
	--total_gene_dir ${total_gene_csv}\
	--coverage_out_file ${coverage_output_name}
	"""
}