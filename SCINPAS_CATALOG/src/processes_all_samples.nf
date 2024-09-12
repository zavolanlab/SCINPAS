#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FIND_EXONS_GENES_BED{
	
	label "custom_python"
	label "short_time"
	label "middle_memory"
	cpus = 1
	publishDir "${params.folder_template}/result/${params.sample_type}/common", mode: 'copy'

	input:
	path (input_gtf)
	path python_script

	output:
	path("${params.exons_bed_out}")
	path("${params.genes_bed_out}")

	script:
	"""
	python3 ${python_script}\
	--gtf_dir ${input_gtf}\
	--genes_bed_out ${params.genes_bed_out}\
	--exons_bed_out ${params.exons_bed_out}\
	--species ${params.sample_type}
	"""
}

process FIND_TERMINAL_EXONS{

	label "custom_python"
	label "short_time"
	label "middle_memory"
	cpus = params.heavy_cores
	publishDir "${params.folder_template}/result/${params.sample_type}/common", mode: 'copy'

	input:
	path (input_gtf)
	path python_script

	output:
	path "${params.terminal_exons_out}"

	script:
	"""
	python3 ${python_script} --gtf_file ${input_gtf} --bed_out ${params.terminal_exons_out} --n ${params.heavy_cores} --species ${params.sample_type}
	"""
}

// Extended version (for catalog)
process PREPARE_IN_SLC_CATALOG{

	label "bash"

	input:
	tuple val(directory), val(scinpas), val(organ)
	val species

	output:
	// pass as "val" and then when you pass full data as input use "path". This is to avoid scope issue. 
	// (you didnt receive full data as an input here but you are trying to use it as a file)
	tuple val("${directory}/${scinpas}.bam"), val("${directory}/${scinpas}.bam.bai"), val("${scinpas}"), val("${organ}")

	script:
	"""
	"""
}

process FILTER_MULTIMAPPING_CATALOG{
	
	echo true
	label "samtools"
	label "long_time"
	memory {50.GB * task.attempt}
	cpus = params.middle_cores
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
	
	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas}-${organs}/multimap_filtered", mode: 'copy'

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

	# commands are executed consecutively. First, samtools view, then samtools sort and samtools index
	samtools view -@ ${params.middle_cores} ${bam} -F 0x904 -bq 255 | samtools sort -@ ${params.middle_cores} - -o ${scinpas}-${organs}_unique_mapped_sorted.bam; samtools index -@ ${params.middle_cores} ${scinpas}-${organs}_unique_mapped_sorted.bam

	# total, mapped, unique, unmapped are executed sequentially but within total, within mapped etc, they run in parallel.
	# -F 0X900 means exclude secondary and supplementary alignments = total reads
	total=\$(samtools view -@ ${params.middle_cores} ${bam} -F 0X900 -c)
	
	# -F 0X900 means exclude secondary, supplementary alignments and unmapped reads = mapped reads
	mapped=\$(samtools view -@ ${params.middle_cores} ${bam} -F 0X904 -c)
	
	# uniquely mapped reads
	unique=\$(samtools view -@ ${params.middle_cores} ${scinpas}-${organs}_unique_mapped_sorted.bam -c)
	
	# unmapped reads
	unmapped=\$(samtools view -@ ${params.middle_cores} -f 4 ${bam} -c)

	csvfile="${scinpas}-${organs}_readInfo.csv"
	echo "\${total}, \${mapped}, \${unmapped}, \${unique}" >> \${csvfile}

	"""         
}

process CONCAT_CSVS_CATALOG{
	label "custom_python"
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/mapped_unmapped", mode: 'copy'
	
	input:
	path(csvs)
	val(out_name)
	path(python_script)

	output:
	path("${out_name}.csv")
	path("${out_name}_good.csv")
	path("${out_name}_bad.csv")
	path("${params.first_filtered_samples_out}")

	script:
	"""
	python3 ${python_script}\
	csv_inputs ${csvs}\
	--out_csv ${out_name}\
	--original_input ${params.samples_list}\
	--modified_input_out ${params.first_filtered_samples_out}
	"""
}

process MODIFY_TUPLE{

	label "bash"

	input:
	tuple val(scinpas), val(organ), path(bam), path(bai), val(percentage), val(directory)

	output:
	tuple path("${bam}"), path("${bai}"), val("${scinpas}"), val("${organ}")

	script:
	"""
	"""
}

process SPLIT_PHASE1_CATALOG{
	
	echo true
	label "samtools"
	label "short_time"
	label "middle_memory"
	memory {40.GB * task.attempt}
	cpus = params.middle_cores
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
	
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas}-${organs}/multimap_filtered", mode: 'copy'

	input:
	tuple path(bam), path(bai), val(scinpas), val(organs), val(chrom)
	
	output:
	tuple path("*_possorted_*.bam"), path("*_possorted_*.bam.bai"), val("${scinpas}"), val("${organs}")

	script:
	"""
	# samtools view -bh ${bam} chr${chrom} -o ${scinpas}-${organs}_intermediate_${chrom}.bam
	# samtools sort ${scinpas}-${organs}_intermediate_${chrom}.bam -o ${scinpas}-${organs}_possorted_${chrom}.bam
	# samtools index ${scinpas}-${organs}_possorted_${chrom}.bam

	# Determine chromosome format based on species
	if [[ "${params.sample_type}" == "worm" ]]
	then
		CHROM_FORMAT=${chrom}
	else
		CHROM_FORMAT="chr${chrom}"
	fi

	# echo "Using CHROM_FORMAT: \${CHROM_FORMAT}"

	samtools view -@ ${params.middle_cores} -bh ${bam} \${CHROM_FORMAT} | samtools sort -@ ${params.middle_cores} - -o ${scinpas}-${organs}_possorted_${chrom}.bam
	samtools index -@ ${params.middle_cores} ${scinpas}-${organs}_possorted_${chrom}.bam
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

	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas}-${organs}/trimming_and_deduplication", mode: 'copy'

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
	cpus = params.middle_cores
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5

	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas}-${organs}/trimming_and_deduplication", mode: 'copy'
	
	input:
	tuple path(bams), val(scinpas), val(organs)

	output:
	tuple path("*_sorted.bam"), path("*_sorted.bam.bai"), val("${scinpas}-${organs}")

	script:
	"""	
	input=\$(basename ${bams})
	prefix=\$(echo \$input | cut -d '.' -f 1)
	
	samtools sort -@ ${params.middle_cores} ${bams} -o \${prefix}_sorted.bam
	samtools index -@ ${params.middle_cores} \${prefix}_sorted.bam
	"""
}

process MERGE_DEDUP_CATALOG{

	label "samtools"
	label "middle_time"
	cpus = params.heavy_cores
	memory {40.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5

	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/trimming_and_deduplication", mode: 'copy'
	
	// input channel looks like: [([bam1, bam2...... bam6.....bamY], sample1), ([bam1, bam2...... bam6.....bamY], sample2)]
	// this is to merge once and only once rather than merging several times
	input:
	tuple path(bams), path(bais), val(scinpas_organs)
	
	output:
	tuple path ("${scinpas_organs}_deduplicated_full_sorted.bam"), path ("${scinpas_organs}_deduplicated_full_sorted.bam.bai"), val("${scinpas_organs}")
		
	script:
	"""
	# samtools merge -f -o ${scinpas_organs}_deduplicated_full.bam ${scinpas_organs}_deduplicated_chr*_sorted.bam
	# samtools sort ${scinpas_organs}_deduplicated_full.bam -o ${scinpas_organs}_deduplicated_full_sorted.bam
	# samtools index ${scinpas_organs}_deduplicated_full_sorted.bam
	# rm ${scinpas_organs}_deduplicated_full.bam

	samtools merge -@ ${params.heavy_cores} -f - ${scinpas_organs}_deduplicated_chr*_sorted.bam | samtools sort -@ ${params.heavy_cores} - -o ${scinpas_organs}_deduplicated_full_sorted.bam
	samtools index -@ ${params.heavy_cores} ${scinpas_organs}_deduplicated_full_sorted.bam
	"""
}

process FASTQC_CATALOG{
	
	echo true
	label "fastqc"
	label "short_time"
	memory {10.GB * task.attempt}
	cpus = params.heavy_cores
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
	
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/fastQC", mode: 'copy'

	input:
	tuple path(bam), path(bai), val(scinpas_organs)

	output:
	path("*_fastqc.zip")

	script:
	"""	
	fastqc -t ${params.heavy_cores} ${bam} -o ./
	"""         
}

process FASTQC_SWARM_PLOT_CATALOG{
	
	echo true
	label "custom_python"
	label "short_time"
	memory {10.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5
	
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/fastQC_swarm", mode: 'copy'

	input:
	path(fastqc_zips)
	path(first_filtered)
	path(python_script)

	output:
	path("${params.sequence_quality_csv}_full.csv")
	path("${params.swarm_plot_out}")
	path("${params.sequence_quality_csv}_filtered.csv")
	path("${params.second_filtered_samples_out}")

	script:
	"""	
	python3 ${python_script}\
	in_fastqc ${fastqc_zips}\
	--swarm_out ${params.swarm_plot_out}\
	--csv_out ${params.sequence_quality_csv}\
	--threshold ${params.seq_qual_threshold}\
	--first_filtered_dir ${first_filtered}\
	--second_filtered_out ${params.second_filtered_samples_out}\
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
	path(genome_fasta)
	path python_script

	output:
	tuple path("${scinpas_organs}_${params.bam_after_fixation_alter}_*"), val("${scinpas_organs}")

	script:
	"""
	python3 ${python_script}\
	--bam_file ${bam}\
	--fasta ${genome_fasta}\
	--bam_out ${scinpas_organs}_${params.bam_after_fixation_alter}
	"""
}

// sort on single full data for all samples 
process SORT_PHASE2_CATALOG{

	echo true
	label "samtools"
	label "short_time"
	memory {40.GB * task.attempt}
	cpus = params.middle_cores
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

	samtools sort -@ ${params.middle_cores} ${bam} -o \${prefix}_sorted.bam
	samtools index -@ ${params.middle_cores} \${prefix}_sorted.bam
	"""
}

process GET_POLYA_CATALOG{

	label "custom_python"
	label "short_time"
	memory {50.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 5

	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/get_polyA", mode: 'copy'
	
	input:
	// for each specific sample, run python script
	tuple path(dedup_full_bam), path(dedup_full_bai), val(scinpas_organs)
	val polyA_out
	val nonpolyA_out
	path (genome_fasta)
	path python_script

	output:
	tuple path("${scinpas_organs}_${polyA_out}_chr*.bam"), val("${scinpas_organs}")
	
	script:
	"""
	python3 ${python_script} --bam_input ${dedup_full_bam}\
	--o_polyA ${scinpas_organs}_${polyA_out}.bam\
	--o_nonpolyA ${scinpas_organs}_${nonpolyA_out}.bam\
	--fasta ${genome_fasta}\
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

	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/get_polyA", mode: 'copy'
	
	// input channel looks like: [([bam1, bam2...... bam6.....bamY], sample1), ([bam1, bam2...... bam6.....bamY], sample2)]
	// this is to merge once and only once rather than merging several times
	input:
	tuple path(bams), path(bais), val(scinpas_organs)
	val polyA_out
	
	output:
	tuple path ("${scinpas_organs}_${polyA_out}_full_sorted.bam"), path ("${scinpas_organs}_${polyA_out}_full_sorted.bam.bai"), val("${scinpas_organs}")
		
	script:
	"""
	# samtools merge -f -o ${scinpas_organs}_${polyA_out}_full.bam ${scinpas_organs}_${polyA_out}_chr*_sorted.bam
	# samtools sort ${scinpas_organs}_${polyA_out}_full.bam -o ${scinpas_organs}_${polyA_out}_full_sorted.bam
	# samtools index ${scinpas_organs}_${polyA_out}_full_sorted.bam
	# rm ${scinpas_organs}_${polyA_out}_full.bam

	samtools merge -@ ${params.heavy_cores} -f - ${scinpas_organs}_${polyA_out}_chr*_sorted.bam | samtools sort -@ ${params.heavy_cores} - -o ${scinpas_organs}_${polyA_out}_full_sorted.bam
	samtools index -@ ${params.heavy_cores} ${scinpas_organs}_${polyA_out}_full_sorted.bam
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
	python3 ${python_script}\
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
	label "middle_memory"
	memory {40.GB * task.attempt}
	cpus = params.middle_cores
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
	# Determine chromosome format based on species
	if [[ "${params.sample_type}" == "worm" ]]
	then
		CHROM_FORMAT=${chrom}
	else
		CHROM_FORMAT="chr${chrom}"
	fi

	samtools view -@ ${params.middle_cores} -bh ${bam} \${CHROM_FORMAT} | samtools sort -@ ${params.middle_cores} - -o  ${scinpas_organs}_${out_name}_sorted_${chrom}.bam
	samtools index -@ ${params.middle_cores} ${scinpas_organs}_${out_name}_sorted_${chrom}.bam
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
	
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/further_split_bed_clustering", mode: 'copy'
	
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
	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/all_organs_merged_bed_clustering", mode: 'copy'
	
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
	tuple val(chromosome), val(direction), path("${clustered_bed}_${chromosome}_${direction}.bed")
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

process GET_INTRONIC_BED{

	label "custom_python"
	label "short_time"
	label "super_heavy_memory"
	cpus = params.heavier_cores
	memory {250.GB * task.attempt}
	errorStrategy {task.exitStatus in 137..140 ? 'retry' : 'terminate'}
	maxRetries 6

	publishDir "${params.folder_template}/result/${params.sample_type}/common", mode: 'copy'

	input:
	path (exons)
	val (out_name)
	path (python_script)

	output:
	path "${out_name}"

	script:
	"""
	python3 ${python_script}\
	--exon_bed ${exons}\
	--out_name ${out_name}\
	--n ${params.heavier_cores}
	"""
}

process CHANGE_BED{

	label "custom_python"
	cpus = 1
	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/modified_PAS", mode: 'copy'

	input:
	tuple val(seqid), val(strand), path(clustered_bed_file)
	path python_script

	output:
	tuple val(seqid), val(strand), path("modified*.bed")

	script:
	"""
	python3 ${python_script} --bed_in ${clustered_bed_file}
	"""
}

process BED_INTERSECT_CATALOG{

	label "bedtools"
	label "middle_memory"
	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/classification", mode: 'copy'

	input:
	// if channel of length = 1, it does not matter if two channels have different lengths or not (bed file is length of 1 so it is ok).
	tuple val(seqid), val(strand), path(clustered_bed_file)
	// either terminal_exons, exons or genes.bed
	path(bed)
	val(intersect_out_name)
	
	output:
	tuple val(seqid), val(strand), path("${intersect_out_name}_${seqid}_${strand}.bed")
	
	script:
	"""
	# window by default reports both A and B if "extended" entry of A overlap with entry of B.
	# -w 1 extend 1bp on both direction
	# -sm means report overlap on the same strand.
	# -u means report A only once if extended entry of A overlap with entry of B.
	bedtools window -w 1 -sm -u -a ${clustered_bed_file} -b ${bed} > ${intersect_out_name}_${seqid}_${strand}.bed
	"""
}

process BED_NO_INTERSECT_CATALOG{

	label "bedtools"
	label "middle_memory"
	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/classification", mode: 'copy'

	input:
	// if channel of length = 1, it does not matter if two channels have different lengths or not (bed file is length of 1 so it is ok).
	tuple val(seqid), val(strand), path(clustered_bed_file)
	// either terminal_exons, exons or genes.bed
	path(bed)
	val(no_intersect_out_name)
	
	output:
	tuple val(seqid), val(strand), path("${no_intersect_out_name}_${seqid}_${strand}.bed")

	script:
	"""
	bedtools window -w 1 -sm -v -a ${clustered_bed_file} -b ${bed} > ${no_intersect_out_name}_${seqid}_${strand}.bed
	"""
}

process BED_INTERSECT_CATALOG_ANTISENSE{

	label "bedtools"
	label "middle_memory"
	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/antisense_classification", mode: 'copy'

	input:
	// if channel of length = 1, it does not matter if two channels have different lengths or not (bed file is length of 1 so it is ok).
	tuple val(seqid), val(strand), path(clustered_bed_file)
	// either terminal_exons, exons or genes.bed
	path(bed)
	val(intersect_out_name)
	
	output:
	tuple val(seqid), val(strand), path("${intersect_out_name}_${seqid}_${strand}.bed")
	
	script:
	"""
	# window by default reports both A and B if "extended" entry of A overlap with entry of B.
	# -w 1 extend 1bp on both direction
	# -sm means report overlap on the same strand.
	# -u means report A only once if extended entry of A overlap with entry of B.
	bedtools window -w 1 -u -a ${clustered_bed_file} -b ${bed} > ${intersect_out_name}_${seqid}_${strand}.bed
	"""
}

process BED_NO_INTERSECT_CATALOG_ANTISENSE{

	label "bedtools"
	label "middle_memory"
	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/antisense_classification", mode: 'copy'

	input:
	// if channel of length = 1, it does not matter if two channels have different lengths or not (bed file is length of 1 so it is ok).
	tuple val(seqid), val(strand), path(clustered_bed_file)
	// either terminal_exons, exons or genes.bed
	path(bed)
	val(no_intersect_out_name)
	
	output:
	tuple val(seqid), val(strand), path("${no_intersect_out_name}_${seqid}_${strand}.bed")

	script:
	"""
	bedtools window -w 1 -v -a ${clustered_bed_file} -b ${bed} > ${no_intersect_out_name}_${seqid}_${strand}.bed
	"""
}

process ADD_CLASS_COLUMN{

	label "custom_python"
	cpus = 1
	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/modified_PAS", mode: 'copy'

	input:
	tuple val(seqid), val(strand), path(beds)
	path python_script

	output:
	path("${params.add_class_col_out}_${seqid}_${strand}.bed")

	script:
	"""
	python3 ${python_script} in_bed ${beds} --bed_out ${params.add_class_col_out}_${seqid}_${strand}.bed
	"""
}

process MERGE_ADD_COL_BED{
	label "bash"
	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/modified_PAS", mode: 'copy'
	
	input:
	path(beds)
	val(out_name)

	output:
	path("${out_name}")

	script:
	"""
	cat *.bed > ${out_name}
	"""
}

process GENE_ID{
	label "custom_python"
	label "middle_memory"
	label "short_time"
	cpus = params.basic_cores
	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/modified_PAS", mode: 'copy'
	
	input:
	path(pas_bed)
	path(genes_bed)
	val(out_template)
	path (python_script)

	output:
	path("${out_template}_${params.version}_intersect_out.bed")

	script:
	"""
	python3 ${python_script}\
	--pas ${pas_bed}\
	--genes ${genes_bed}\
	--out ${out_template}_${params.version}\
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
	tuple path("${organs}_${params.organ_specific_cs_out}_filtered_${chromosome}_${direction}.bed"), val("${organs}")

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

process GET_ORGAN_SCORE{
	label "custom_python"
	label "middle_memory"
	label "short_time"
	cpus = params.heavy_cores
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/organ_score", mode: 'copy'
	
	input:
	tuple path(organ_beds), val(organs)
	val(out_name)
	path (python_script)

	output:
	path("${organs}_${out_name}")

	script:
	"""
	python3 ${python_script}\
	in_bed ${organ_beds}\
	--out_name ${organs}_${out_name}\
	--organ ${organs}\
	--n ${params.heavy_cores}\
	"""
}

process MERGE_ORGAN_SCORE_TO_PAS{
	label "custom_python"
	label "middle_memory"
	label "short_time"
	cpus = params.basic_cores
	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/modified_PAS", mode: 'copy'
	
	input:
	path(organ_scores)
	path(pas)
	val(out_template)
	path (python_script)

	output:
	path("${out_template}_${params.version}_${params.merge_pas_organ_score}")

	script:
	"""
	python3 ${python_script}\
	in_bed ${organ_scores}\
	--pas ${pas}\
	--out_name ${out_template}_${params.version}_${params.merge_pas_organ_score}\
	"""
}

process RCS_MOTIF_CHECK{

	label "custom_python"
	label "middle_memory"
	label "short_time"
	cpus = 1
	// publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/modified_PAS", mode: 'copy'

	input:
	path(pas_bed)
	val(out_template)
	val(out_template2)
	path(genome_fasta)
	path python_script

	output:
	path ("${params.sample_type}_${params.version}_${out_template}_${out_template2}")

	script:
	"""
	python3 ${python_script}\
	--rcs_dir ${pas_bed}\
	--motif_dir ${params.motif_dir}\
	--fasta_dir ${genome_fasta}\
	--rcs_out ${params.sample_type}_${params.version}_${out_template}_${out_template2}
	"""
}

process CONVERT_GZIP_UNIQUE_CS{

	label "bash"

	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/${scinpas_organs}/split_bed_clustering", mode: 'copy'

	input:
	tuple path(unique_cs_beds), val(scinpas_organs), val(chrom)

	output:
	path "*.gz"

	script:
	"""
	gzip -c ${unique_cs_beds} > ${unique_cs_beds}.gz
	"""
}

process CONVERT_GZIP_ALL_SAMPLES_UNIQUE_CS{

	label "bash"

	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/all_organs_merged_bed_clustering", mode: 'copy'

	input:
	tuple path(all_unique_cs_beds), val(chromosome), val(direction)

	output:
	path "*.gz"

	script:
	"""
	gzip -c ${all_unique_cs_beds} > ${all_unique_cs_beds}.gz
	"""
}

process CONVERT_GZIP_CATALOG{

	label "bash"

	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/catalog_out", mode: 'copy'

	input:
	path (pas_bed)

	output:
	path "*.gz"

	script:
	"""
	gzip -c ${pas_bed} > ${pas_bed}.gz
	"""
}