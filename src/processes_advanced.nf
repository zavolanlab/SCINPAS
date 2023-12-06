#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process PREPARE_IN_SLC {

	label "bash"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}", mode: 'copy'
	
	input:
	val specific_sample
	val species

	output:
	// pass as "val" and then when you pass full data as input use "path". This is to avoid scope issue. 
	// (you didnt receive full data as an input here but you are trying to use it as a file)
	tuple val("${params.folder_template}/data/${species}/${params.result_folder}/${specific_sample}.bam"), val("${params.folder_template}/data/${species}/${params.result_folder}/${specific_sample}.bam.bai"), val("${specific_sample}")

	script:
	"""
	"""
}

process FIND_EXONS_GENES_BED{
	
	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/common", mode: 'copy'

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

process FIND_TERMINAL_EXONS{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/common", mode: 'copy'

	input:
	path python_script

	output:
	path "${params.terminal_exons_out}"

	script:
	"""
	python3 ${python_script} --gtf_file ${params.gtf_fasta_location_template}/${params.sample_type}/${params.gtf} --bed_out ${params.terminal_exons_out}
	"""
}

process FILTER_FURTHER_TERMINAL_EXONS{
	
	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/common", mode: 'copy'

	input:
	path terminal_exons
	path python_script

	output:
	path("${params.further_filtered_terminal_exons_out}")

	script:
	"""
	python3 ${python_script}\
	--terminal_exons_bed_input ${terminal_exons}\
	--bed_out ${params.further_filtered_terminal_exons_out}
	"""
}

process SPLIT_PHASE1 {
	
	echo true
	label "samtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/trimming_and_deduplication", mode: 'copy'

	input:
	tuple path(bam), path(bai), val(specific_sample), val(chrom)

	output:
	tuple path("${specific_sample}_possorted_*.bam"), path("${specific_sample}_possorted_*.bam.bai"), val("${specific_sample}")
	tuple val("${specific_sample}"), path("bam_folder_${chrom}")

	script:
	"""
	# echo ${bam}
	# echo ${chrom}
	samtools view -bh ${bam} chr${chrom} -o ${specific_sample}_intermediate_${chrom} 

	samtools sort ${specific_sample}_intermediate_${chrom} -o ${specific_sample}_possorted_${chrom}.bam
		
	samtools index ${specific_sample}_possorted_${chrom}.bam

	mkdir bam_folder_${chrom}
	# copy possorted_number.bam, bai to bamfolder_number folder
	cp ${specific_sample}_possorted_* bam_folder_${chrom}
	
	"""         
}

process DEDUP {

	label "custom_python"
	label "heavy_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/trimming_and_deduplication", mode: 'copy'

	input:
	tuple val(specific_sample), path(bam_folder)
	path python_script
	
	output:
	tuple path ("${specific_sample}_deduplicated_chr*"), val("${specific_sample}")
	tuple path ("${specific_sample}_splitted_chr*"), val("${specific_sample}")

	script:
	"""
	# copy everything from bam folders (1~19, X, Y) to the current working directory.
	cp ${bam_folder}/* .
	# * is for placeholder not "all"
	python3 ${python_script} \
		--bam_template ${specific_sample}_possorted_*.bam \
		--out_template ${specific_sample}_${params.dedup_out_name_prefix}\
		--span_threshold ${params.span_threshold}\
		--split_read_template ${specific_sample}_${params.split_out_name_prefix}
	"""
}

process SORT_PHASE1{
	
	echo true
	label "samtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/trimming_and_deduplication/", mode: 'copy'
	
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

process CHANGE_TUPLE_ORDER{
	label "bash"
	input:
	tuple path(bams), path(bais), val(specific_sample)
	
	output:
	tuple path("${bams}"), val("${specific_sample}"), path("${bais}")

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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/trimming_and_deduplication", mode: 'copy'
	
	// input channel looks like: [([bam1, bam2...... bam6.....bamY], sample1), ([bam1, bam2...... bam6.....bamY], sample2)]
	// this is to merge once and only once rather than merging several times
	input:
	tuple path(bams), path(bais), val(specific_sample)

	output:
	tuple path ("${specific_sample}_dedup_full_sorted.bam"), path ("${specific_sample}_dedup_full_sorted.bam.bai"), val("${specific_sample}")
		
	script:
	"""
	samtools merge -f -o ${specific_sample}_dedup_full.bam ${specific_sample}_deduplicated_chr*_output_sorted.bam
	samtools sort ${specific_sample}_dedup_full.bam -o ${specific_sample}_dedup_full_sorted.bam
	samtools index ${specific_sample}_dedup_full_sorted.bam
	"""
}

process MERGE_CONTROL{

	label "samtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/get_polyA", mode: 'copy'
	
	// input channel looks like: [([bam1, bam2...... bam6.....bamY], sample1), ([bam1, bam2...... bam6.....bamY], sample2)]
	// this is to merge once and only once rather than merging several times
	input:
	tuple path(bams), path(bais), val(specific_sample)

	output:
	tuple path ("${specific_sample}_dedup_fixed_full_sorted.bam"), path ("${specific_sample}_dedup_fixed_full_sorted.bam.bai"), val("${specific_sample}")
		
	script:
	"""
	samtools merge -f -o ${specific_sample}_dedup_fixed_full.bam ${specific_sample}_DedupNegative_ctrl_reads_sorted_*.bam
	samtools sort ${specific_sample}_dedup_fixed_full.bam -o ${specific_sample}_dedup_fixed_full_sorted.bam
	samtools index ${specific_sample}_dedup_fixed_full_sorted.bam
	"""
}

process TRIM_NEG_CONTROL {

	label "custom_python"
	label "super_heavy_computation"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/trimming_and_deduplication", mode: 'copy'

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

process FIX_SOFTCLIPPED_REGION{
	
	label "custom_python"
	label "super_heavy_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/trimming_and_deduplication", mode: 'copy'

	input:
	tuple path(bam), path(bai), val(specific_sample)
	path python_script

	output:
	tuple path("${params.bam_after_fixation}"), val("${specific_sample}")
	tuple path("num_fixed_unfixed.csv"), val("${specific_sample}")

	script:
	"""
	python3 ${python_script}\
	--bam_file ${bam}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--bam_out ${params.bam_after_fixation}

	"""
}

process GET_SPAN_RAW{

	label "custom_python"
	label "heavy_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/span_distribution", mode: 'copy'

	input:
	tuple path(raw_bams), path(raw_bais), val(specific_sample)
	val raw
	val sample_or_control
	path python_script

	output:
	tuple path("${specific_sample}_${params.raw_span}_*.csv"), val("${specific_sample}")
	
	script:
	"""
	if [[ ${sample_or_control} == "sample" ]]
	then
		raw_input=\$(basename ${raw_bams})
		raw_prefix=\$(echo \$raw_input | cut -d '_' -f 5)
		raw_number=\$(echo \$raw_prefix | cut -d '.' -f 1)

		python3 ${python_script}\
		--input ${raw_bams}\
		--output ${specific_sample}_${params.raw_span}_\${raw_number}.csv\
		--action_type ${raw}

	elif [[ ${sample_or_control} == "control" ]]
	then
		raw_input=\$(basename ${raw_bams})
		raw_prefix=\$(echo \$raw_input | cut -d '_' -f 8)
		raw_number=\$(echo \$raw_prefix | cut -d '.' -f 1)

		python3 ${python_script}\
		--input ${raw_bams}\
		--output ${specific_sample}_${params.raw_span}_\${raw_number}.csv\
		--action_type ${raw}
	fi	
	"""
}

process GET_SPAN_DEDUP{

	label "custom_python"
	label "middle_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/span_distribution", mode: 'copy'

	input:
	tuple path(dedup_bams), path(dedup_bais), val(specific_sample)
	val dedup
	val sample_or_control
	path python_script

	output:
	tuple path("${specific_sample}_${params.dedup_span}_*.csv"), val("${specific_sample}")
	
	script:
	"""
	if [[ ${sample_or_control} == "sample" ]]
	then
		dedup_input=\$(basename ${dedup_bams})
		dedup_prefix=\$(echo \$dedup_input | cut -d '_' -f 5)
		dedup_number=\$(echo \$dedup_prefix | cut -d 'r' -f 2)

		python3 ${python_script}\
		--input ${dedup_bams}\
		--output ${specific_sample}_${params.dedup_span}_\${dedup_number}.csv\
		--action_type ${dedup}

	elif [[ ${sample_or_control} == "control" ]]
	then
		dedup_input=\$(basename ${dedup_bams})
		dedup_prefix=\$(echo \$dedup_input | cut -d '_' -f 8)
		dedup_number=\$(echo \$dedup_prefix | cut -d '.' -f 1)

		python3 ${python_script}\
		--input ${dedup_bams}\
		--output ${specific_sample}_${params.dedup_span}_\${dedup_number}.csv\
		--action_type ${dedup}	

	fi		
	"""
}

process GET_SPAN_SPLIT{

	label "custom_python"
	label "heavy_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/span_distribution", mode: 'copy'

	input:
	tuple path(split_bams), path(split_bais), val(specific_sample)
	val split
	path python_script

	output:
	tuple path("${specific_sample}_${params.split_span}_*.csv"), val("${specific_sample}")

	script:
	"""
	split_input=\$(basename ${split_bams})
	split_prefix=\$(echo \$split_input | cut -d '_' -f 5)
	split_number=\$(echo \$split_prefix | cut -d 'r' -f 2)

	python3 ${python_script}\
	--input ${specific_sample}_splitted_chr\${split_number}_output_sorted.bam\
	--output ${specific_sample}_${params.split_span}_\${split_number}.csv\
	--action_type ${split}
	"""
}

process CONCAT_SPAN{

	label "bash"
	label "heavy_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/span_distribution", mode: 'copy'

	input:
	// channel looks like: [([csv1, csv2...... csv6.....csvY], sample1), ([csv1, csv2...... csv6.....csvY], sample2)]
	tuple path(csv_files), val(specific_sample)
	val(file_name)

	output:
	tuple path("${specific_sample}_combined_${file_name}.csv"), val("${specific_sample}")

	script:
	"""
	cat ${specific_sample}_${file_name}_*.csv > ${specific_sample}_combined_${file_name}.csv
	# cat ${params.dedup_span}_*.csv > combined_dedup_span.csv
	"""
}

process GET_SPAN_HIST{

	label "custom_python"
	label "heavy_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/span_distribution", mode: 'copy'

	input:
	tuple path(combined_spans), val(specific_sample)
	path python_script
	val(file_name)

	output:
	path "${specific_sample}_${file_name}*"

	script:
	"""
	python3 ${python_script}\
	--input_csv ${combined_spans}\
	--out_hist ${specific_sample}_${file_name}
	"""
}

process GET_POLYA{

	label "custom_python"
	label "super_heavy_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/get_polyA", mode: 'copy'
	
	input:
	// for each specific sample, run python script
	tuple path(dedup_full_bam), path(dedup_full_bai), val(specific_sample)
	val polyA_out
	val nonpolyA_out
	path python_script

	output:
	tuple path("${polyA_out}.bam"), val("${specific_sample}")
	tuple path("${nonpolyA_out}.bam"), val("${specific_sample}")
	
	script:
	"""
	python3 ${python_script} --bam_input ${dedup_full_bam}\
	--o_polyA ${polyA_out}.bam\
	--o_nonpolyA ${nonpolyA_out}.bam\
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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/get_polyA", mode: 'copy'
	
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

process GET_POLYA_UNIQUE_CLEAVAGE_SITES_BED{

	label "custom_python"
	label "middle_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/split_bed_clustering", mode: 'copy'
	
	input:
	tuple path(bams), path(bais), val(specific_sample)
	val(polyA_unique_cs_bed)
	path python_script

	output:
	tuple path("${specific_sample}_*.bed"), val("${specific_sample}")

	script:
	"""
	input=\$(basename ${bams})
	prefix=\$(echo \$input | cut -d '.' -f 1)
	number=\$(echo \$prefix | cut -d '_' -f 7)
	
	python3 ${python_script}\
	--bam ${bams}\
	--bed_out ${specific_sample}_${polyA_unique_cs_bed}\
	--use_fc ${params.use_fc}\
	--split \${number}
	"""
}

process MERGED_BED{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/merged_bed_clustering", mode: 'copy'
	
	input:
	tuple path(beds), val(specific_sample)
	val(out_name)
	path python_script

	output:
	tuple path("${specific_sample}_${out_name}"), val("${specific_sample}")

	script:
	"""
	python3 ${python_script}\
	in_bed ${beds}\
	--bed_template ${specific_sample}_${out_name}
	"""
}

process PERFORM_CLUSTERING_TUPLE{

	label "custom_python"
	label "heavy_computation"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/merged_bed_clustering", mode: 'copy'
	
	input:
	tuple path(cleavage_site_bed), val(specific_sample)
	val(clustered_bed)
	path python_script

	output:
	tuple path("${specific_sample}_${clustered_bed}"), val("${specific_sample}")

	script:
	"""
	python3 ${python_script}\
	in_bed ${cleavage_site_bed}\
	--out ${specific_sample}_${clustered_bed}\
	--du ${params.cluter_up}\
	--dd ${params.cluster_down}\
	--c ${params.num_cores_clustering}
	"""
}

process BED_INTERSECT_TUPLE{

	label "bedtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/bed_classification", mode: 'copy'

	input:
	// if channel of length = 1, it does not matter if two channels have different lengths or not (bed file is length of 1 so it is ok).
	tuple path(clustered_bed_file), val(specific_sample)
	// either terminal_exons, exons or genes.bed
	path(bed)
	val(intersect_out_name)
	
	output:
	tuple path("*_${intersect_out_name}.bed"), val("${specific_sample}")
	
	script:
	"""
	input=\$(basename ${clustered_bed_file})
	prefix=\$(echo \$input | cut -d '_' -f 1-3)
	# bedtools intersect -s -u -nonamecheck -a ${clustered_bed_file} -b ${bed} -wa -bed > \${prefix}_${intersect_out_name}.bed
	
	# window by default reports both A and B if "extended" entry of A overlap with entry of B.
	# -w 1 extend 1bp on both direction
	# -sm means report overlap on the same strand.
	# -u means report A only once if extended entry of A overlap with entry of B.
	bedtools window -w 1 -sm -u -a ${clustered_bed_file} -b ${bed} > \${prefix}_${intersect_out_name}.bed
	"""
}

process BED_NO_INTERSECT_TUPLE{

	label "bedtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/bed_classification", mode: 'copy'

	input:
	// if channel of length = 1, it does not matter if two channels have different lengths or not (bed file is length of 1 so it is ok).
	tuple path(clustered_bed_file), val(specific_sample)
	// either terminal_exons, exons or genes.bed
	path(bed)
	val(no_intersect_out_name)
	
	output:
	tuple path("*_${no_intersect_out_name}.bed"), val("${specific_sample}")

	script:
	"""
	input=\$(basename ${clustered_bed_file})
	prefix=\$(echo \$input | cut -d '_' -f 1-3)
	# bedtools intersect -s -v -nonamecheck -a ${clustered_bed_file} -b ${bed} -wa -bed > \${prefix}_${no_intersect_out_name}.bed
	bedtools window -w 1 -sm -v -a ${clustered_bed_file} -b ${bed} > \${prefix}_${no_intersect_out_name}.bed
	"""
}

process CLASSIFY_BAM{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/split_bam_classification", mode: 'copy'
	
	input:
	// all_files contain polyA.bam, polyA.bai, TE.bed, E.bed, Intronic.bed and Intergenic.bed in a mixed order. It does not matter. we can access them correctly by their names
	tuple val(specific_sample), path(bam), path(bai), path(input_bed)
	val(bam_out)
	val(remaining_bam)
	val(class_reads)
	path python_script

	output:
	tuple path("${specific_sample}_${bam_out}_*.bam"), val("${specific_sample}")
	tuple path("${specific_sample}_${remaining_bam}_*.bam"), val("${specific_sample}")

	script:
	"""
	python3 ${python_script}\
	--bam ${bam}\
	--input_bed ${input_bed}\
	--use_fc ${params.use_fc}\
	--out_bam ${specific_sample}_${bam_out}\
	--remaining_bam ${specific_sample}_${remaining_bam}\
	--class_reads ${class_reads}
	"""
}

process MERGE_CLASSIFIED{

	label "samtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/bam_classification", mode: 'copy'
	
	// input channel looks like: [[bam1, bam2...... bam6.....bamY], sample1, [bai1, bai2...... bai6.....baiY]]
	// this is to merge once and only once rather than merging several times
	input:
	tuple path(bams), val(specific_sample), path(bais)
	val(output_template)

	output:
	tuple path ("${specific_sample}_${output_template}_full_sorted.bam"), path ("${specific_sample}_${output_template}_full_sorted.bam.bai"), val("${specific_sample}")
		
	script:
	"""
	samtools merge -f -o ${specific_sample}_${output_template}_full.bam ${specific_sample}_${output_template}_*.bam
	samtools sort ${specific_sample}_${output_template}_full.bam -o ${specific_sample}_${output_template}_full_sorted.bam
	samtools index ${specific_sample}_${output_template}_full_sorted.bam
	"""
}

process ISOLATE_PA_T{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/bed_classification", mode: 'copy'
	
	input:
	tuple path(bed), val(specific_sample)
	path(terminal_exons)
	path python_script

	output:
	tuple path("${specific_sample}_${params.below_template}.bed"), val("${specific_sample}")
	tuple path("${specific_sample}_${params.above_template}.bed"), val("${specific_sample}")
	
	script:
	"""
	python3 ${python_script}\
	--bed_input ${bed}\
	--terminal_exons_gtf ${terminal_exons}\
	--annotated_out ${specific_sample}_${params.below_template}.bed\
	--unannotated_out ${specific_sample}_${params.above_template}.bed\
	--distance_threshold ${params.distance_threshold}\
	"""
}

process CLASSIFY_BELOW_ABOVE{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/split_bam_classification", mode: 'copy'
	
	input:
	tuple val(specific_sample), path(bam), path(bai), path(input_bed)
	val(bam_out)
	val(remaining_bam)
	val(class_reads)
	path python_script

	output:
	tuple path("${specific_sample}_${bam_out}_*.bam"), val("${specific_sample}")
	tuple path("${specific_sample}_${remaining_bam}_*.bam"), val("${specific_sample}")

	script:
	"""
	python3 ${python_script}\
	--bam ${bam}\
	--input_bed ${input_bed}\
	--use_fc ${params.use_fc}\
	--out_bam ${specific_sample}_${bam_out}\
	--remaining_bam ${specific_sample}_${remaining_bam}\
	--class_reads ${class_reads}
	"""
}

process SORT_CLASSIFIED{

	echo true
	label "samtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/split_bam_classification", mode: 'copy'
	
	input:
	tuple path(bam), val(specific_sample)

	output:
	tuple path("*_sorted.bam"), val("${specific_sample}"), path("*_sorted.bam.bai")

	script:
	"""	
	input=\$(basename ${bam})
	prefix=\$(echo \$input | cut -d '.' -f 1)

	samtools sort ${bam} -o \${prefix}_sorted.bam
	samtools index \${prefix}_sorted.bam
	"""
}

process SPLIT_PHASE2{

	echo true
	label "samtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/split_bam_classification", mode: 'copy'
	
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

process GET_NON_SOFTCLIP{

	label "custom_python"
	label "super_heavy_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/get_A_freq_in_nonsoftclipped_reads", mode: 'copy'
	
	input:
	// for each specific sample, run python script
	tuple path(dedup_full_bam), path(dedup_full_bai), val(specific_sample)
	path python_script

	output:
	tuple path("${specific_sample}_${params.no_softclipped_reads}"), val("${specific_sample}")

	
	script:
	"""
	python3 ${python_script}\
	--bam_input ${dedup_full_bam}\
	--o_file ${specific_sample}_${params.no_softclipped_reads}
	"""

}

process SORT_PHASE3{

	echo true
	label "samtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/get_A_freq_in_nonsoftclipped_reads", mode: 'copy'
	
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

process GET_A_FREQ_IN_NONSOFTCLIP{

	label "custom_python"
	label "middle_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/get_A_freq_in_nonsoftclipped_reads", mode: 'copy'
	
	input:
	// for each specific sample, run python script
	tuple path(bam), path(bai), val(specific_sample)
	path python_script

	output:
	path "${specific_sample}_${params.avg_a_freq_csv}"

	script:
	"""
	python3 ${python_script}\
	--bam_input ${bam}\
	--o_file ${specific_sample}_${params.avg_a_freq_csv}
	"""
}

process GET_COUNTS{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/num_reads_per_class", mode: 'copy'

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
	publishDir "${params.folder_template}/result/${params.sample_type}/common", mode: 'copy'
	
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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/separate_merged_counts", mode: 'copy'
	
	input:
	path(raw_concat_csv)
	path(dedup_concat_csv)
	path(dedup_softclipped_concat_csv)
	path(nonpolyA_concat_csv)
	path(polyA_concat_csv)
	path(polyA_T_concat_csv)
	path(annotated_concat_csv)
	path(unannotated_concat_csv)
	path(intronic_concat_csv)
	path(intergenic_concat_csv)
	path(exonic_concat_csv)
	val(out_name)
	
	output:
	path("${out_name}.csv")

	script:
	"""
	cat *.csv > ${out_name}.csv
	"""
}

process SPLIT_PER_SAMPLE{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/separate_merged_counts", mode: 'copy'

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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/num_reads_per_class", mode: 'copy'
	
	input:
	tuple path(csvs), val(specific_sample), path(control_csv)
	path python_script
	
	output:
	path "${specific_sample}_${params.barchart_output}.png"
	path "${specific_sample}_${params.barchart_output}.svg"

	script:
	"""
	python3 ${python_script}\
	--input_dir ${csvs}\
	--control_input_dir ${control_csv}\
	--out_file ${specific_sample}_${params.barchart_output}
	"""
}

process GET_DISTANCE_ALTER{

	label "custom_python"
	label "heavy_memory"
	label "long_time"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/distance", mode: 'copy'

	input:
	tuple path(polyA_t_full_bed), val(specific_sample)
	path terminal_exons
	path python_script

	output:
	tuple path("${specific_sample}_${params.distance_output}"), val("${specific_sample}")
	
	script:
	"""
	python3 ${python_script}\
	--bed_input ${polyA_t_full_bed}\
	--gtf_input ${terminal_exons}\
	--distance_out ${specific_sample}_${params.distance_output}
	"""
}

process UNCOUPLE_CSV_WITH_SAMPLE{

	label "bash"

	input:
	tuple path(csvs), val(specific_sample)

	output:
	path("${csvs}")

	script:
	"""
	"""
}

process GET_D_HISTOGRAM_ALTER{

	label "custom_python"
	label "heavy_memory"

	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/distance", mode: 'copy'

	input:
	tuple path(distance_csv), val(specific_sample), path(negative_control_distance_csv)
	path python_script

	output:
	path "${params.distance_histogram_file_name}.png"
	path "${params.distance_histogram_file_name}.svg"

	script:
	"""
	python3 ${python_script}\
	--csv_input ${distance_csv}\
	--cntrl_csv ${negative_control_distance_csv}\
	--file_name ${params.distance_histogram_file_name}
	"""
}

process PLOT_ATGC{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/ATGC_plot", mode: 'copy'
	// pass all bam files at once (not one by one)
	// and then within python script it uses bam file one by one.
	input:
	tuple path(beds), val(specific_sample)
	val output_name
	path python_script

	output:
	path "${specific_sample}_${output_name}*"

	script:
	"""
	python3 ${python_script}\
	--bed ${beds}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--out_name ${specific_sample}_${output_name}\
	--window_size ${params.window_size}
	"""
}

process GET_MOTIF_FREQ_PLOT{

	label "custom_python"
	label "heavy_computation"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/motif_freq_plot", mode: 'copy'

	input:
	tuple path(bed), val(specific_sample)
	val output_name
	val annotated
	path peaks_csv
	path python_script

	output:
	path ("${specific_sample}_${output_name}*.png")
	path ("${specific_sample}_${output_name}_overlaid.png")
	path ("${specific_sample}_${output_name}_overlaid.svg")
	path ("${specific_sample}_${output_name}_ordered_motives.csv")
	path ("${specific_sample}_${output_name}_scores.csv")

	script:
	"""
	python3 ${python_script}\
	--bed ${bed}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--out_name ${specific_sample}_${output_name}\
	--annotated ${annotated}\
	--window ${params.motif_search_window}\
	--motif_info_dir ${params.motif_info_dir}\
	--downstream ${params.down_limit}\
	--peaks ${peaks_csv}\
	--n_cores ${params.heavy_cores}
	"""
}

process GET_MOTIF_FREQ_GTF_PLOT{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/gtf_motif_freq_plots", mode: 'copy'

	input:
	path terminal_exons
	val output_name
	path python_script

	output:
	path "${output_name}_*.png"
	path "${output_name}_overlaid.png"
	path "${output_name}_overlaid.svg"
	path "${output_name}_ordered_*.csv"
	path "${output_name}_peaks.csv"

	script:
	"""
	python3 ${python_script}\
	--terminal_exons ${terminal_exons}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--out_name ${output_name}\
	--window ${params.motif_search_window}\
	--motif_info_dir ${params.motif_info_dir}\
	--downstream ${params.down_limit}
	"""
}

process CONCAT_MOTIF_ORDERS{
	label "bash"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/motif_orders", mode: 'copy'
	
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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/motif_orders", mode: 'copy'
	
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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/separate_merged_motif_orders", mode: 'copy'

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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/motif_freq_plot", mode: 'copy'

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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/motif_freq_plot", mode: 'copy'

	input:
	tuple path(beds), val(specific_sample)
	val output_name
	path peaks_csv
	path python_script

	output:
	tuple path ("${specific_sample}_${output_name}_overlaid.png"), val("${specific_sample}")
	tuple path ("${specific_sample}_${output_name}_overlaid.svg"), val("${specific_sample}")
	path "${specific_sample}_${output_name}_*.csv"
	
	script:
	"""
	python3 ${python_script}\
	--bed ${beds}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--out_name ${specific_sample}_${output_name}\
	--window ${params.motif_search_window}\
	--motif_info_dir ${params.motif_info_dir}\
	--downstream ${params.down_limit}\
	--peaks ${peaks_csv}\
	--n_cores ${params.heavy_cores}
	"""
}

process COMBINE_CSVS{

	label "bash"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/common", mode: 'copy'

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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/polyA_motif_score", mode: 'copy'
	
	input:
	path(csvs)
	path python_script
	
	output:
	path "${params.polyA_motif_score}.png"
	path "${params.polyA_motif_score}.svg"

	script:
	"""
	python3 ${python_script}\
	--input_dir ${csvs}\
	--out_file ${params.polyA_motif_score}
	"""
}

process GET_SOFTCLIPPED_DISTRIBUTION{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/softclipped_distribution", mode: 'copy'

	input:
	tuple path(bams), path(bais), val(specific_sample)
	val output_name
	path python_script

	output:
	tuple path ("${specific_sample}_${output_name}.png"), val("${specific_sample}")
	tuple path ("${specific_sample}_${output_name}.svg"), val("${specific_sample}")
	tuple path ("${specific_sample}_*.csv"), val("${specific_sample}")

	script:
	"""
	python3 ${python_script}\
	--dedup_bam ${bams}\
	--file_name ${specific_sample}_${output_name}
	"""
}

process GET_NUM_POLYA_SITES{
	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/num_polyA", mode: 'copy'

	input:
	tuple path(beds), val(specific_sample)
	val csv_output_name
	val annotated
	path python_script

	output:
	path "${specific_sample}_${csv_output_name}"
	
	script:
	"""
	python3 ${python_script}\
	--bed ${beds}\
	--csv_out_name ${specific_sample}_${csv_output_name}\
	--annotated ${annotated}\
	--cluster_threshold ${params.cluster_threshold}\
	--n_cores ${params.basic_cores}
	"""	
}

process MERGE_ALL_NUM_PAS_CSVS{
	label "bash"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/num_polyA", mode: 'copy'
	
	input:
	path(unannotated_concat_csv)
	path(annotated_concat_csv)
	path(intergenic_concat_csv)
	path(intronic_concat_csv)
	path(exonic_concat_csv)
	path(control_csv)
	val(out_name)

	output:
	path("${out_name}")

	script:
	"""
	cat *.csv > ${out_name}
	"""
}

process FUSE_MOTIF_SCORE_NUM_PAS{
	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/num_polyA", mode: 'copy'
	
	input:
	path(total_motif_scores)
	path(polyA_motif_scores)
	path(total_num_polyA_csv)
	val(out_name)
	path python_script

	output:
	path("${out_name}")

	script:
	"""
	python3 ${python_script}\
	--total_motif_scores_csv ${total_motif_scores}\
	--polyA_csv ${polyA_motif_scores}\
	--total_num_polyA_csv ${total_num_polyA_csv}\
	--csv_out_name ${out_name}	
	"""
}
process GET_BAR_CHART_NUM_POLYA{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/num_polyA", mode: 'copy'
	
	input:
	path(csvs)
	val output_name
	path python_script
	
	output:
	path "${output_name}.png"
	path "${output_name}.svg"

	script:
	"""
	python3 ${python_script}\
	--input_dir ${csvs}\
	--out_file ${output_name}
	"""
}

process GET_GENE_COVERAGE{
	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/gene_coverage", mode: 'copy'

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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/gene_coverage", mode: 'copy'
	
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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/gene_coverage", mode: 'copy'
	
	input:
	path dedup_csv
	path polyA_csv
	path total_gene_csv
	val coverage_output_name
	path python_script
	
	output:
	path "${coverage_output_name}.png"
	path "${coverage_output_name}.svg"

	script:
	"""
	python3 ${python_script}\
	--dedup_csv_dir ${dedup_csv}\
	--polyA_csv_dir ${polyA_csv}\
	--total_gene_dir ${total_gene_csv}\
	--coverage_out_file ${coverage_output_name}
	"""
}

process SOFTCLIP_PWM_LOGO{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/PWM_logo", mode: 'copy'
	
	input:
	tuple path(bams), path(bais), val(specific_sample)
	val length_softclip
	path python_script
	
	output:
	path "${specific_sample}_${params.softclip_pwm_logo}*"
	
	script:
	"""
	python3 ${python_script}\
	--bam_file ${bams}\
	--pwm_signature ${specific_sample}_${params.softclip_pwm_logo}\
	--length_softclipped ${length_softclip}\
	--use_FC ${params.use_fc}
	"""
}

process SPLIT_CSV_BY_CELL_TYPE{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/cell_type_specific", mode: 'copy'
	
	input:
	tuple val(specific_sample), path(metadata)
	path python_script

	output:
	tuple path("*_${params.type1}.txt"), path("*_${params.type2}.txt"), val("${specific_sample}")

	script:
	"""
	python3 ${python_script}\
		--cell_type_dir ${metadata}\
		--txt_output_template ${params.out_template}\
		--type1 ${params.type1}\
		--type2 ${params.type2}\
		--sample_name ${specific_sample}\
		--sample_type ${params.result_folder}
	"""
}

process SPLIT_SAMPLE_BY_TYPE{

	label "samtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${specific_sample}/cell_type_specific", mode: 'copy'
	
	input:
	tuple val(specific_sample), path(dedup_bam), path(dedup_bai), path(type1_txt), path(type2_txt)

	output:
	path("${specific_sample}_${params.type1}_sorted.bam")
	path("${specific_sample}_${params.type1}_sorted.bam.bai")
	path("${specific_sample}_${params.type2}_sorted.bam")
	path("${specific_sample}_${params.type2}_sorted.bam.bai")
	
	script:
	"""
	samtools view -D CB:${type1_txt} -bh ${dedup_bam} -o ${specific_sample}_${params.type1}.bam  
	samtools sort ${specific_sample}_${params.type1}.bam -o ${specific_sample}_${params.type1}_sorted.bam  
	samtools index ${specific_sample}_${params.type1}_sorted.bam

	samtools view -D CB:${type2_txt} -bh ${dedup_bam} -o ${specific_sample}_${params.type2}.bam  
	samtools sort ${specific_sample}_${params.type2}.bam -o ${specific_sample}_${params.type2}_sorted.bam  
	samtools index ${specific_sample}_${params.type2}_sorted.bam
	"""
}

process MERGE_BY_CELL_TYPE{

	label "samtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/merged_cell_type", mode: 'copy'
	
	input:
	path(type1_bam)
	path(type1_bai)
	path(type2_bam)
	path(type2_bai)

	output:
	tuple path("merged_${params.type1}_sorted.bam"), path("merged_${params.type1}_sorted.bam.bai") 
	tuple path("merged_${params.type2}_sorted.bam"), path("merged_${params.type2}_sorted.bam.bai") 

	script:
	"""
	samtools merge -f -o merged_${params.type1}.bam *_${params.type1}_sorted.bam
	samtools sort merged_${params.type1}.bam -o merged_${params.type1}_sorted.bam
	samtools index merged_${params.type1}_sorted.bam

	samtools merge -f -o merged_${params.type2}.bam *_${params.type2}_sorted.bam
	samtools sort merged_${params.type2}.bam -o merged_${params.type2}_sorted.bam
	samtools index merged_${params.type2}_sorted.bam
	"""
}

process GET_POLYA_CELL_TYPE{

	label "custom_python"
	label "super_heavy_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/merged_cell_type", mode: 'copy'
	
	input:
	tuple path(dedup_full_bam), path(dedup_full_bai)
	val polyA_out
	val nonpolyA_out
	path python_script

	output:
	path("${polyA_out}")
	path("${nonpolyA_out}")
	
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

process SORT_CELL_TYPE{

	echo true
	label "samtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/merged_cell_type", mode: 'copy'
	
	input:
	path(bam)

	output:
	tuple path("*_sorted.bam"), path("*_sorted.bam.bai")

	script:
	"""	
	input=\$(basename ${bam})
	prefix=\$(echo \$input | cut -d '.' -f 1)

	samtools sort ${bam} -o \${prefix}_sorted.bam
	samtools index \${prefix}_sorted.bam
	"""
}

process SPLIT_CELL_TYPE{

	echo true
	label "samtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/merged_cell_type", mode: 'copy'
	
	input:
	tuple path(bam), path(bai), val(chrom)
	val(cell_type)
	val(out_name)

	output:
	tuple path("${out_name}_*.bam"), path ("${out_name}_*.bam.bai")

	script:
	"""	
	samtools view -bh ${bam} chr${chrom} | samtools sort -o  ${out_name}_sorted_${chrom}_${cell_type}.bam
	samtools index ${out_name}_sorted_${chrom}_${cell_type}.bam
	"""         
}

process GET_CELL_TYPE_POLYA_UNIQUE_CS_BED{

	label "custom_python"
	label "middle_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/merged_cell_type_by_chrom", mode: 'copy'
	
	input:
	tuple path(bams), path(bais)
	val(polyA_unique_cs_bed)
	path python_script

	output:
	path("all_polyA_cs_*.bed")

	script:
	"""
	input=\$(basename ${bams})
	prefix=\$(echo \$input | cut -d '.' -f 1)
	number=\$(echo \$prefix | cut -d '_' -f 4)
		
	python3 ${python_script}\
	--bam ${bams}\
	--bed_out ${polyA_unique_cs_bed}\
	--use_fc ${params.use_fc}\
	--split \${number}
	"""
	
}

process MERGED_BED_CELL_TYPE{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/${params.result_folder}/merged_cell_type", mode: 'copy'
	
	input:
	path(beds)
	val(out_name)
	path python_script

	output:
	path("${out_name}")

	script:
	"""
	python3 ${python_script}\
	in_bed ${beds}\
	--bed_template ${out_name}
	"""
}

process PERFORM_CLUSTERING_CELL_TYPE{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/merged_cell_type", mode: 'copy'
	
	input:
	path(cleavage_site_bed)
	val(clustered_bed)
	path python_script

	output:
	path("${clustered_bed}")

	script:
	"""
	python3 ${python_script}\
	in_bed ${cleavage_site_bed}\
	--out ${clustered_bed}\
	--du ${params.cluter_up}\
	--dd ${params.cluster_down}\
	--c ${params.num_cores_clustering}
	"""
}

process BED_INTERSECT_CELL_TYPE{

	label "bedtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/merged_cell_type", mode: 'copy'

	input:
	path(clustered_bed_file)
	// either terminal_exons, exons or genes.bed
	path(bed)
	val(type)
	val(intersect_out_name)
	
	output:
	path("${type}_${intersect_out_name}.bed")
	
	script:
	"""
	# bedtools intersect -s -u -nonamecheck -a ${clustered_bed_file} -b ${bed} -wa -bed > ${type}_${intersect_out_name}.bed
	bedtools window -w 1 -sm -u -a ${clustered_bed_file} -b ${bed} > ${type}_${intersect_out_name}.bed
	"""
}

process BED_NO_INTERSECT_CELL_TYPE{

	label "bedtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/merged_cell_type", mode: 'copy'

	input:
	path(clustered_bed_file)
	// either terminal_exons, exons or genes.bed
	path(bed)
	val(type)
	val(no_intersect_out_name)
	
	output:
	path("${type}_${no_intersect_out_name}.bed")
	
	script:
	"""
	# bedtools intersect -s -v -nonamecheck -a ${clustered_bed_file} -b ${bed} -wa -bed > ${type}_${no_intersect_out_name}.bed
	bedtools window -w 1 -sm -v -a ${clustered_bed_file} -b ${bed} > ${type}_${no_intersect_out_name}.bed
	"""
}

process DISTANCE_SCATTER{

	label "custom_python"
	label "heavy_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/merged_cell_type", mode: 'copy'
	
	input:
	path(polyA_type1_bed)
	path(polyA_type2_bed)
	path terminal_exons
	path python_script

	output:
	path "${params.type1_type2_distance_csv}"

	script:
	"""
	python3 ${python_script}\
	--type1_polyA_bed ${polyA_type1_bed}\
	--type2_polyA_bed ${polyA_type2_bed}\
	--bed_input ${terminal_exons}\
	--distance_out ${params.type1_type2_distance_csv}
	"""
}

process T1_T2_SCATTER{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/merged_cell_type", mode: 'copy'
	
	input:
	path t1_t2_distance_csv
	path python_script

	output:
	path "${params.type1_type2_scatter}*"

	script:
	"""
	python3 ${python_script}\
	--distance_csv_dir ${t1_t2_distance_csv}\
	--output_name ${params.type1_type2_scatter}\
	--type1 ${params.type1}\
	--type2 ${params.type2}
	"""
}

process GET_CELL_TYPE_MOTIF_FREQ_PLOT{

	label "custom_python"
	label "heavy_computation"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/merged_cell_type", mode: 'copy'

	input:
	path (bed)
	val output_name
	val annotated
	path peaks_csv
	val type
	path python_script

	output:
	path ("${type}_${output_name}*.png")
	path ("${type}_${output_name}_overlaid.png")
	path ("${type}_${output_name}_overlaid.svg")
	path ("${type}_${output_name}_ordered_motives.csv")
	path ("${type}_${output_name}_scores.csv")

	script:
	"""
	python3 ${python_script}\
	--bed ${bed}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--out_name ${type}_${output_name}\
	--annotated ${annotated}\
	--window ${params.motif_search_window}\
	--motif_info_dir ${params.motif_info_dir}\
	--downstream ${params.down_limit}\
	--peaks ${peaks_csv}\
	--n_cores ${params.heavy_cores}
	"""
}

process GET_CELL_TYPE_NUM_POLYA_SITES{
	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/merged_cell_type", mode: 'copy'

	input:
	path(bed)
	val csv_output_name
	val annotated
	val type
	path python_script

	output:
	path "${type}_${csv_output_name}"
	
	script:
	"""
	python3 ${python_script}\
	--bed ${bed}\
	--csv_out_name ${type}_${csv_output_name}\
	--annotated ${annotated}\
	--cluster_threshold ${params.cluster_threshold}\
	--n_cores ${params.basic_cores}
	"""	
}

process UNCOUPLE_BAM_WITH_SAMPLE{

	label "bash"

	input:
	tuple path(bam), path(bai), val(specific_sample)

	output:
	path("${bam}")
	path("${bai}")

	script:
	"""
	"""
}

process MERGE_ALL_SAMPLES{

	label "samtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/overlap_with_catalog", mode: 'copy'
	
	input:
	path(bams)
	path(bais)

	output:
	tuple path ("all_samples_polyA_full_sorted.bam"), path ("all_samples_polyA_full_sorted.bam.bai")
		
	script:
	"""
	samtools merge -f -o all_samples_polyA_full.bam 10X*sorted.bam
	samtools sort all_samples_polyA_full.bam -o all_samples_polyA_full_sorted.bam
	samtools index all_samples_polyA_full_sorted.bam
	"""
}

process SPLIT_ALL_POLYA_ALL_SAMPLES{

	echo true
	label "samtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/overlap_with_catalog", mode: 'copy'
	
	input:
	tuple path(bam), path(bai), val(chrom)
	val(out_name)

	output:
	tuple path("${out_name}_*.bam"), path ("${out_name}_*.bam.bai")

	script:
	"""	
	samtools view -bh ${bam} chr${chrom} | samtools sort -o  ${out_name}_sorted_${chrom}.bam
	samtools index ${out_name}_sorted_${chrom}.bam
	"""         
}

process GET_POLYA_UNIQUE_CLEAVAGE_SITES_ALL_SAMPLES{

	label "custom_python"
	label "middle_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/overlap_with_catalog", mode: 'copy'
	
	input:
	tuple path(bam), path(bai)
	val(polyA_unique_cs_bed)
	path python_script

	output:
	path("*.bed")

	script:
	"""
	input=\$(basename ${bam})
	prefix=\$(echo \$input | cut -d '.' -f 1)
	number=\$(echo \$prefix | cut -d '_' -f 4)

	python3 ${python_script}\
	--bam ${bam}\
	--bed_out ${polyA_unique_cs_bed}\
	--use_fc ${params.use_fc}\
	--split \${number}
	"""
}

process MERGED_BED_ALL_SAMPLES{

	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/overlap_with_catalog", mode: 'copy'
	
	input:
	path(beds)
	val(out_name)
	path python_script

	output:
	path("${out_name}")

	script:
	"""
	python3 ${python_script}\
	in_bed ${beds}\
	--bed_template ${out_name}
	"""
}

process PERFORM_CLUSTERING_ALL_SAMPLES{

	label "custom_python"
	label "heavy_computation"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/overlap_with_catalog", mode: 'copy'
	
	input:
	path(cleavage_site_bed)
	val(clustered_bed)
	path python_script

	output:
	path("all_${clustered_bed}")

	script:
	"""
	python3 ${python_script}\
	in_bed ${cleavage_site_bed}\
	--out all_${clustered_bed}\
	--du ${params.cluter_up}\
	--dd ${params.cluster_down}\
	--c ${params.num_cores_clustering}
	"""
}

process GET_OVERLAP{
	label "custom_python"
	label "heavy_computation"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/overlap_with_catalog", mode: 'copy'

	input:
	path(pas_bed_dir)
	path python_script

	output:
	path "${params.result_folder}_${params.overlap_catalog_figure}.png"
	path "${params.result_folder}_${params.overlap_catalog_figure}.svg"

	script:
	"""
	python3 ${python_script}\
	--pas_bed_dir ${pas_bed_dir}\
	--catalog_dir ${params.gtf_fasta_location_template}/${params.sample_type}/${params.catalog}\
	--overlap_barchart ${params.result_folder}_${params.overlap_catalog_figure}\
	--species ${params.sample_type}\
	--n_cores ${params.num_cores_overlap}\
	--threshold ${params.cluster_threshold}
	"""	
}

process GET_TOP_OVERLAP{
	label "custom_python"
	label "heavy_computation"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.result_folder}/overlap_with_catalog", mode: 'copy'

	input:
	path(pas_bed_dir)
	path python_script

	output:
	path "${params.result_folder}_${params.top_overlap_catalog_figure}.png"
	path "${params.result_folder}_${params.top_overlap_catalog_figure}.svg"
	
	script:
	"""
	python3 ${python_script}\
	--pas_bed_dir ${pas_bed_dir}\
	--catalog_dir ${params.gtf_fasta_location_template}/${params.sample_type}/${params.catalog}\
	--overlap_barchart ${params.result_folder}_${params.top_overlap_catalog_figure}\
	--species ${params.sample_type}\
	--n_cores ${params.num_cores_overlap}\
	--threshold ${params.cluster_threshold}
	"""	
}

// Extended version (for catalog)
process PREPARE_IN_SLC_CATALOG{

	label "bash"

	input:
	tuple val(scinpas), val(organ)
	val species

	output:
	// pass as "val" and then when you pass full data as input use "path". This is to avoid scope issue. 
	// (you didnt receive full data as an input here but you are trying to use it as a file)
	tuple val("${params.folder_template}/data/${species}/${organ}/${scinpas}.bam"), val("${params.folder_template}/data/${species}/${organ}/${scinpas}.bam.bai"), val("${scinpas}"), val("${organ}")

	script:
	"""
	"""
}

process FILTER_MULTIMAPPING_CATALOG{
	
	echo true
	label "samtools"
	label "short_time"
	memory {10.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
	
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas}-${organs}/multimap_filtered", mode: 'copy'

	input:
	tuple path(bam), path(bai), val(scinpas), val(organs)

	output:
	tuple path("*_unique_mapped_sorted.bam"), path("*_unique_mapped_sorted.bam.bai"), val("${scinpas}"), val("${organs}")
	path ("${scinpas}-${organs}_readInfo.csv")

	script:
	"""	
	# samtools view ${bam} -bq 255 -o ${scinpas}-${organs}_mapqfiltered.bam
	# samtools sort ${scinpas}-${organs}_mapqfiltered.bam -o ${scinpas}-${organs}_mapqfiltered_sorted.bam
	# samtools index ${scinpas}-${organs}_mapqfiltered_sorted.bam	

	samtools view ${bam} -F 0x904 -bq 255 | samtools sort - -o ${scinpas}-${organs}_unique_mapped_sorted.bam; samtools index ${scinpas}-${organs}_unique_mapped_sorted.bam

	# -F 0X900 means exclude secondary and supplementary alignments = total reads
	total=\$(samtools view ${bam} -F 0X900 -c)
	
	# -F 0X900 means exclude secondary, supplementary alignments and unmapped reads = mapped reads
	mapped=\$(samtools view ${bam} -F 0X904 -c)
	
	# uniquely mapped reads
	unique=\$(samtools view ${scinpas}-${organs}_unique_mapped_sorted.bam -c)
	
	# unmapped reads
	unmapped=\$(samtools view -f 4 ${bam} -c)

	csvfile="${scinpas}-${organs}_readInfo.csv"
	echo "\${total}, \${mapped}, \${unmapped}, \${unique}" >> \${csvfile}

	"""         
}

process CONCAT_CSVS_CATALOG{
	label "bash"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas}-${organs}/multimap_filtered", mode: 'copy'
	
	input:
	path(csvs)
	val(out_name)
	
	output:
	path("${out_name}.csv")

	script:
	"""
	cat *.csv > ${out_name}.csv
	"""
}

process SPLIT_PHASE1_CATALOG{
	
	echo true
	label "samtools"
	label "short_time"
	memory {10.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
	
	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas}-${organs}/multimap_filtered", mode: 'copy'

	input:
	tuple path(bam), path(bai), val(scinpas), val(organs), val(chrom)

	output:
	tuple path("*_possorted_*.bam"), path("*_possorted_*.bam.bai"), val("${scinpas}"), val("${organs}")

	script:
	"""
	# echo ${bam}
	# echo ${chrom}
	
	samtools view -bh ${bam} chr${chrom} -o ${scinpas}-${organs}_intermediate_${chrom}.bam

	samtools sort ${scinpas}-${organs}_intermediate_${chrom}.bam -o ${scinpas}-${organs}_possorted_${chrom}.bam
		
	samtools index ${scinpas}-${organs}_possorted_${chrom}.bam	
	"""         
}

// the memory and execution time limits are defined dynamically. 
// The first time the process is executed the task.attempt is set to 1, thus it will request a two GB of memory and one hour of maximum execution time.
// If the task execution stops reporting an exit status in the range between 137 and 140, the task is re-submitted (otherwise terminates immediately). 
// This time the value of task.attempt is 2, thus increasing the amount of the memory to four GB and the time to 2 hours, and so on.
// The directive maxRetries set the maximum number of time the same task can be re-executed.

process DEDUP_CATALOG{

	label "custom_python"
	label "short_time"
	memory {250.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5

	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas}-${organs}/trimming_and_deduplication", mode: 'copy'

	input:
	tuple path(bam), path(bai), val(scinpas), val(organs)
	path python_script
	
	output:
	tuple path ("*-${organs}_deduplicated_chr*"), val("${scinpas}"), val("${organs}")

	script:
	"""

	# * is for placeholder not "all"
	python3 ${python_script} \
		--bam_template ${scinpas}-${organs}_possorted_*.bam \
		--out_template ${scinpas}-${organs}_${params.dedup_out_name_prefix}\
		--span_threshold ${params.span_threshold}\
		--split_read_template ${scinpas}-${organs}_${params.split_out_name_prefix}
	"""
}

process SORT_PHASE1_CATALOG{
	
	echo true
	label "samtools"
	label "short_time"
	memory {40.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5

	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas}-${organs}/trimming_and_deduplication", mode: 'copy'
	
	input:
	tuple path(bams), val(scinpas), val(organs)

	output:
	tuple path("*_sorted.bam"), path("*_sorted.bam.bai"), val("${scinpas}-${organs}")

	script:
	"""	
	input=\$(basename ${bams})
	prefix=\$(echo \$input | cut -d '.' -f 1)
	
	samtools sort ${bams} -o \${prefix}_sorted.bam
	samtools index \${prefix}_sorted.bam
	"""
}

process MERGE_DEDUP_CATALOG{

	label "samtools"
	label "middle_time"
	cpus = params.heavy_cores
	memory {40.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5

	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/trimming_and_deduplication", mode: 'copy'
	
	// input channel looks like: [([bam1, bam2...... bam6.....bamY], sample1), ([bam1, bam2...... bam6.....bamY], sample2)]
	// this is to merge once and only once rather than merging several times
	input:
	tuple path(bams), path(bais), val(scinpas_organs)
	
	output:
	tuple path ("${scinpas_organs}_deduplicated_full_sorted.bam"), path ("${scinpas_organs}_deduplicated_full_sorted.bam.bai"), val("${scinpas_organs}")
		
	script:
	"""
	samtools merge -f -o ${scinpas_organs}_deduplicated_full.bam ${scinpas_organs}_deduplicated_chr*_sorted.bam
	samtools sort ${scinpas_organs}_deduplicated_full.bam -o ${scinpas_organs}_deduplicated_full_sorted.bam
	samtools index ${scinpas_organs}_deduplicated_full_sorted.bam
	"""
}

process FASTQC_CATALOG{
	
	echo true
	label "fastqc"
	label "short_time"
	memory {10.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
	
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/fastQC", mode: 'copy'

	input:
	tuple path(bam), path(bai), val(scinpas_organs)

	output:
	path("*_fastqc.zip")

	script:
	"""	
	fastqc ${bam} -o ./
	"""         
}

process FASTQC_SWARM_PLOT_CATALOG{
	
	echo true
	label "custom_python"
	label "short_time"
	memory {10.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
	
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/fastQC_swarm", mode: 'copy'

	input:
	path(fastqc_zips)

	output:
	path("${params.sequence_quality_csv}")
	path("${params.swarm_plot_out}")
	
	script:
	"""	
	python3 ${python_script}\
	in_bed ${fastqc_zips}\
	--swarm_out ${params.swarm_plot_out}\
	--csv_out ${params.sequence_quality_csv}
	"""         
}

process SPLIT_FILTERED_DEDUP_CATALOG{
	
	echo true
	label "samtools"
	label "short_time"
	memory {10.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
	
	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/trimming_and_deduplication", mode: 'copy'

	input:
	tuple path(bam), path(bai), val(scinpas_organs), val(chrom)

	output:
	tuple path("${scinpas_organs}_selected_deduplicated_chr*_sorted.bam"), path("${scinpas_organs}_selected_deduplicated_chr*_sorted.bam.bai"), val("${scinpas_organs}")

	script:
	"""
	samtools view -bh ${bam} chr${chrom} -o ${scinpas_organs}_selected_deduplicated_chr${chrom}.bam

	samtools sort ${scinpas_organs}_selected_deduplicated_chr${chrom}.bam -o ${scinpas_organs}_selected_deduplicated_chr${chrom}_sorted.bam
		
	samtools index ${scinpas_organs}_selected_deduplicated_chr${chrom}_sorted.bam
	"""         
}

process FIX_SOFTCLIPPED_REGION_CATALOG{
	
	label "custom_python"
	label "short_time"
	memory {60.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 6
	
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/trimming_and_deduplication", mode: 'copy'

	input:
	tuple path(bam), path(bai), val(scinpas_organs)
	path python_script

	output:
	tuple path("${scinpas_organs}_${params.bam_after_fixation_alter}_*"), val("${scinpas_organs}")

	script:
	"""
	python3 ${python_script}\
	--bam_file ${bam}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--bam_out ${scinpas_organs}_${params.bam_after_fixation_alter}
	"""
}

// sort on single full data for all samples 
process SORT_PHASE2_CATALOG{

	echo true
	label "samtools"
	label "short_time"
	memory {40.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
		
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/get_polyA", mode: 'copy'
	
	input:
	tuple path(bam), val(scinpas_organs)

	output:
	tuple path("*_sorted.bam"), path("*_sorted.bam.bai"), val("${scinpas_organs}")

	script:
	"""	
	input=\$(basename ${bam})
	prefix=\$(echo \$input | cut -d '.' -f 1)

	samtools sort ${bam} -o \${prefix}_sorted.bam
	samtools index \${prefix}_sorted.bam
	"""
}

process GET_POLYA_CATALOG{

	label "custom_python"
	label "short_time"
	memory {50.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5

	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/get_polyA", mode: 'copy'
	
	input:
	// for each specific sample, run python script
	tuple path(dedup_full_bam), path(dedup_full_bai), val(scinpas_organs)
	val polyA_out
	val nonpolyA_out
	path python_script

	output:
	tuple path("${scinpas_organs}_${polyA_out}_*.bam"), val("${scinpas_organs}")
	tuple path("${scinpas_organs}_${nonpolyA_out}_*.bam"), val("${scinpas_organs}")
	
	script:
	"""
	python3 ${python_script} --bam_input ${dedup_full_bam}\
	--o_polyA ${scinpas_organs}_${polyA_out}.bam\
	--o_nonpolyA ${scinpas_organs}_${nonpolyA_out}.bam\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--percentage_threshold ${params.polyA_percentage_threshold}\
	--length_threshold ${params.length_threshold}\
	--use_fc ${params.use_fc}
	"""
}

process MERGE_POLYA_CATALOG{

	label "samtools"
	label "middle_time"
	cpus = params.heavy_cores
	memory {40.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5

	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/get_polyA", mode: 'copy'
	
	// input channel looks like: [([bam1, bam2...... bam6.....bamY], sample1), ([bam1, bam2...... bam6.....bamY], sample2)]
	// this is to merge once and only once rather than merging several times
	input:
	tuple path(bams), path(bais), val(scinpas_organs)
	val polyA_out
	
	output:
	tuple path ("${scinpas_organs}_${polyA_out}_full_sorted.bam"), path ("${scinpas_organs}_${polyA_out}_full_sorted.bam.bai"), val("${scinpas_organs}")
		
	script:
	"""
	samtools merge -f -o ${scinpas_organs}_${polyA_out}_full.bam ${scinpas_organs}_${polyA_out}_chr*_sorted.bam
	samtools sort ${scinpas_organs}_${polyA_out}_full.bam -o ${scinpas_organs}_${polyA_out}_full_sorted.bam
	samtools index ${scinpas_organs}_${polyA_out}_full_sorted.bam
	"""
}

process GET_COUNTS_CATALOG{

	label "custom_python"
	label "short_time"

	input:
	tuple path(bam), path(bai), val(scinpas_organs)
	val(csv_name)
	val(bam_type)
	path python_script
	
	output:
	tuple path("${bam}"), path("${bai}"), path("${scinpas_organs}_${csv_name}"), val("${scinpas_organs}")
	
	script:
	"""
	python3 ${python_script} \
	--bam_input ${bam}\
	--csv_output ${scinpas_organs}_${csv_name}\
	--sample_name ${scinpas_organs}\
	--bam_type ${bam_type}
	"""
}

process SPLIT_PHASE2_CATALOG{

	echo true
	label "samtools"
	label "short_time"
	memory {10.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5

	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/split_bam_classification", mode: 'copy'
	
	input:
	tuple path(bam), path(bai), path(polyA_count_csv), val(scinpas_organs), val(chrom)
	val(out_name)

	output:
	tuple path("${scinpas_organs}_${out_name}_*.bam"), path("${scinpas_organs}_${out_name}_*.bam.bai"), path("${polyA_count_csv}"), val("${scinpas_organs}"), val("${chrom}")

	script:
	"""	
	samtools view -bh ${bam} chr${chrom} | samtools sort -o  ${scinpas_organs}_${out_name}_sorted_${chrom}.bam
	samtools index ${scinpas_organs}_${out_name}_sorted_${chrom}.bam
	"""         
}

process GET_POLYA_UNIQUE_CLEAVAGE_SITES_CATALOG{

	label "custom_python"
	label "short_time"
	memory {40.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
	
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/split_bed_clustering", mode: 'copy'
	
	input:
	tuple path(bams), path(bais), path(polyA_count_csv), val(scinpas_organs), val(chrom)
	val(polyA_unique_cs_bed)
	val(multiple_samples)
	path python_script

	output:
	// tuple path("${scinpas_organs}_*.bed"), val("${scinpas_organs}")
	tuple path("${scinpas_organs}_*.bed"), val("${scinpas_organs}"), val("${chrom}")
	
	script:
	"""
	input=\$(basename ${bams})
	prefix=\$(echo \$input | cut -d '.' -f 1)
	prefix2=\$(echo \$prefix | cut -d '_' -f 4)
	if [[ \$prefix2 == "immune" ]]
	then	
		number=\$(echo \$prefix | cut -d '_' -f 8)
	else
		number=\$(echo \$prefix | cut -d '_' -f 7)
	fi

	# multiple samples: whether use RPM and the number of experiments that support a particular cleavage site
	# for catalog purpose, this value should be 1. otherwise set to 0.
	python3 ${python_script}\
	--bam ${bams}\
	--bed_out ${scinpas_organs}_${polyA_unique_cs_bed}\
	--use_fc ${params.use_fc}\
	--split \${number}\
	--multiple_samples ${multiple_samples}\
	--count_dir ${polyA_count_csv}
	"""
}

process SPLIT_BY_DIRECTION{

	label "custom_python"
	label "short_time"
	memory {40.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
	
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/all_organs_split_bed_clustering", mode: 'copy'
	
	input:
	tuple path(beds), val(scinpas_organs), val(chromosome), val(direction)
	path python_script

	output:
	// tuple path("${scinpas_organs}_*.bed"), val("${scinpas_organs}")
	tuple path("*ForClustering_${chromosome}_${direction}.bed"), val("${chromosome}"), val("${direction}")
	tuple path("*ForOrganGrouping_${chromosome}_${direction}.bed"), val("${chromosome}"), val("${direction}"), val("${scinpas_organs}")
	
	script:
	"""
	python3 ${python_script}\
	--bed_in ${beds}\
	--direction ${direction}\
	--chromo ${chromosome}
	"""
}

process GROUPBY_BED_CATALOG{

	label "custom_python"
	label "very_long_time"
	cpus = 5
	memory {150.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5	
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/all_organs_merged_bed_clustering", mode: 'copy'
	
	input:
	tuple path(beds), val(chromosome), val(direction)
	val(out_name)
	val(num_samples)
	path python_script

	output:
	tuple path("${out_name}_${chromosome}_${direction}.bed"), val("${chromosome}"), val("${direction}")

	script:
	"""
	# Pass all bed files of a given chromosome and direction at once.
	python3 ${python_script}\
	in_bed ${beds}\
	--bed_out ${out_name}\
	--chrom ${chromosome}\
	--n_cores ${params.basic_cores}\
	--num_samples ${num_samples}\
	--direction ${direction}
	"""
}

process PERFORM_CLUSTERING_CATALOG{

	label "custom_python"
	label "very_long_time"
	cpus = 2
	memory {50.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
	
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/all_organs_merged_bed_clustering", mode: 'copy'
	
	input:
	tuple path(cleavage_site_bed), val(chromosome), val(direction)
	val(clustered_bed)
	val(modified_cs_out)
	path python_script

	output:
	path("${clustered_bed}*.bed")
	path("${modified_cs_out}*.bed")

	script:
	"""
	python3 ${python_script}\
	--in_bed ${cleavage_site_bed}\
	--out ${clustered_bed}\
	--du ${params.cluter_up}\
	--dd ${params.cluster_down}\
	--c ${params.basic_cores}\
	--num ${chromosome}\
	--strand ${direction}\
	--cs_out ${modified_cs_out}
	"""
}

process LEFT_JOIN_CATALOG{
	label "custom_python"
	label "middle_time"
	cpus = params.basic_cores
	memory {50.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
	
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/per_organ/${organs}/PAS_bed_split", mode: 'copy'
	
	input:
	tuple path(unique_cs_beds), val(chromosome), val(direction), val(organs)
	path (modified_unique_cs_all_samples)
	path (python_script)

	output:
	tuple path("${organs}_${params.organ_specific_cs_out}_${chromosome}_${direction}.bed"), val("${organs}")
	tuple path("${organs}_${params.organ_specific_cs_out}_filtered_${chromosome}_${direction}.bed"), val("${organs}")
	tuple path("${organs}_${params.organ_specific_pas_out}_${chromosome}_${direction}.bed"), val("${organs}")

	script:
	"""
	python3 ${python_script}\
	in_bed ${unique_cs_beds}\
	--modified_unique_cs ${params.modified_unique_cs}_${chromosome}_${direction}.bed\
	--out_cs ${organs}_${params.organ_specific_cs_out}\
	--out_pas ${organs}_${params.organ_specific_pas_out}\
	--chrom ${chromosome}\
	--direction ${direction}\
	--organ ${organs}\
	--n_cores ${params.basic_cores}
	"""
}

process MERGED_BED_CATALOG{

	label "custom_python"
	label "short_time"
	memory {40.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5

	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/per_organ/${scinpas_organs}/PAS_bed", mode: 'copy'
	
	input:
	tuple path(beds), val(scinpas_organs)
	val(out_name)
	path python_script

	output:
	tuple path("${scinpas_organs}_${out_name}.bed"), val("${scinpas_organs}")

	script:
	"""
	python3 ${python_script}\
	in_bed ${beds}\
	--bed_template ${scinpas_organs}_${out_name}
	"""
}

process BED_INTERSECT_CATALOG{

	label "bedtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/per_organ/${organs}/PAS_bed", mode: 'copy'

	input:
	// if channel of length = 1, it does not matter if two channels have different lengths or not (bed file is length of 1 so it is ok).
	tuple path(clustered_bed_file), val(organs)
	// either terminal_exons, exons or genes.bed
	path(bed)
	val(intersect_out_name)
	
	output:
	tuple path("${organs}_${intersect_out_name}.bed"), val("${organs}")
	
	script:
	"""
	# window by default reports both A and B if "extended" entry of A overlap with entry of B.
	# -w 1 extend 1bp on both direction
	# -sm means report overlap on the same strand.
	# -u means report A only once if extended entry of A overlap with entry of B.
	bedtools window -w 1 -sm -u -a ${clustered_bed_file} -b ${bed} > ${organs}_${intersect_out_name}.bed
	"""
}

process BED_NO_INTERSECT_CATALOG{

	label "bedtools"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/per_organ/${organs}/PAS_bed", mode: 'copy'

	input:
	// if channel of length = 1, it does not matter if two channels have different lengths or not (bed file is length of 1 so it is ok).
	tuple path(clustered_bed_file), val(organs)
	// either terminal_exons, exons or genes.bed
	path(bed)
	val(no_intersect_out_name)
	
	output:
	tuple path("${organs}_${no_intersect_out_name}.bed"), val("${organs}")

	script:
	"""
	bedtools window -w 1 -sm -v -a ${clustered_bed_file} -b ${bed} > ${organs}_${no_intersect_out_name}.bed
	"""
}

process ISOLATE_PA_T_CATALOG{

	label "custom_python"
	label "middle_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/per_organ/${organs}/PAS_bed", mode: 'copy'
	
	input:
	tuple path(te_bed), val(organs)
	path(terminal_exons)
	path python_script

	output:
	tuple path("${organs}_${params.below_template}.bed"), val("${organs}")
	tuple path("${organs}_${params.above_template}.bed"), val("${organs}")
	
	script:
	"""
	python3 ${python_script}\
	--bed_input ${te_bed}\
	--terminal_exons_gtf ${terminal_exons}\
	--annotated_out ${organs}_${params.below_template}.bed\
	--unannotated_out ${organs}_${params.above_template}.bed\
	--distance_threshold ${params.distance_threshold}
	"""
}

process PLOT_ATGC_CATALOG{

	label "custom_python"
	label "short_time"
	memory {20.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 10

	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/per_organ/${organs}/ATGC_plot", mode: 'copy'

	input:
	tuple path(beds), val(organs)
	val output_name
	path python_script

	output:
	tuple path("${organs}_${output_name}*"), val("${organs}")

	script:
	"""
	python3 ${python_script}\
	--bed ${beds}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--out_name ${organs}_${output_name}\
	--window_size ${params.window_size}
	"""
}

process GET_MOTIF_FREQ_PLOT_CATALOG{

	label "custom_python"
	label "middle_time"
	cpus = params.heavy_cores
	memory {40.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
	
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/per_organ/${organs}/motif_freq_plot", mode: 'copy'

	input:
	tuple path(bed), val(organs)
	val output_name
	val annotated
	path peaks_csv
	path python_script

	output:
	tuple path ("${organs}_${output_name}*.png"), val("${organs}")
	tuple path ("${organs}_${output_name}_overlaid.png"), val("${organs}")
	tuple path ("${organs}_${output_name}_overlaid.svg"), val("${organs}")
	tuple path ("${organs}_${output_name}_ordered_motives.csv"), val("${organs}")
	tuple path ("${organs}_${output_name}_scores.csv"), val("${organs}")

	script:
	"""
	python3 ${python_script}\
	--bed ${bed}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--out_name ${organs}_${output_name}\
	--annotated ${annotated}\
	--window ${params.motif_search_window}\
	--motif_info_dir ${params.motif_info_dir}\
	--downstream ${params.down_limit}\
	--peaks ${peaks_csv}\
	--n_cores ${params.heavy_cores}
	"""
}

process COMPUTE_SCORES_FOR_MOTIF_ALTER_ALL_SAMPLES{

	label "custom_python"
	label "long_time"
	cpus = params.super_heavy_cores
	memory {50.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5

	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/per_organ/${organs}/motif_score", mode: 'copy'

	input:
	tuple path(bed), val(organs)
	val output_name
	path peaks_csv
	path python_script

	output:
	tuple path ("${organs}_${output_name}_overlaid.png"), val("${organs}")
	tuple path ("${organs}_${output_name}_overlaid.svg"), val("${organs}")
	tuple path ("${organs}_${output_name}_*.csv"), val("${organs}")
	
	script:
	"""
	python3 ${python_script}\
	--bed ${bed}\
	--fasta ${params.gtf_fasta_location_template}/${params.sample_type}/${params.genome_fasta}\
	--out_name ${organs}_${output_name}\
	--window ${params.motif_search_window}\
	--motif_info_dir ${params.motif_info_dir}\
	--downstream ${params.down_limit}\
	--peaks ${peaks_csv}\
	--n_cores ${params.super_heavy_cores}
	"""
}

process GET_NUM_POLYA_SITES_CATALOG{
	label "custom_python"
	label "middle_time"
	cpus = params.heavy_cores
	memory {50.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
	
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/per_organ/${organs}/num_polyA", mode: 'copy'

	input:
	tuple path(bed), val(organs)
	val csv_output_name
	val annotated
	path python_script

	output:
	tuple path ("${organs}_${csv_output_name}"), val("${organs}")
	
	script:
	"""
	python3 ${python_script}\
	--bed ${bed}\
	--csv_out_name ${organs}_${csv_output_name}\
	--annotated ${annotated}\
	--cluster_threshold ${params.cluster_threshold}\
	--n_cores ${params.heavy_cores}
	"""	
}
