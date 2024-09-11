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

process GET_INTRONIC_BED{

	label "custom_python"
	label "short_time"
	label "middle_memory"
	cpus = params.heavier_cores
	memory {60.GB * task.attempt}
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
	--n ${params.mega_heavy_cores}
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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/classification", mode: 'copy'

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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/classification", mode: 'copy'

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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/antisense_classification", mode: 'copy'

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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/antisense_classification", mode: 'copy'

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
	publishDir "${params.folder_template}/result/${params.sample_type}/${params.version}/modified_PAS", mode: 'copy'

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
	publishDir "${params.folder_template}/result/${params.sample_type}/modified_PAS", mode: 'copy'
	
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
	publishDir "${params.folder_template}/result/${params.sample_type}/GENE_ID", mode: 'copy'
	
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
	
	publishDir "${params.folder_template}/result/${params.sample_type}/per_organ/${organs}/PAS_bed_split", mode: 'copy'
	
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
	publishDir "${params.folder_template}/result/${params.sample_type}/organ_score", mode: 'copy'
	
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
	publishDir "${params.folder_template}/result/${params.sample_type}/modified_PAS", mode: 'copy'
	
	input:
	path(organ_scores)
	path(pas)
	val(out_template)
	path (python_script)

	output:
	path("${out_template}_${params.merge_pas_organ_score}")

	script:
	"""
	python3 ${python_script}\
	in_bed ${organ_scores}\
	--pas ${pas}\
	--out_name ${out_template}_${params.merge_pas_organ_score}\
	"""
}

process RCS_MOTIF_CHECK{

	label "custom_python"
	label "heavy_memory"
	label "short_time"
	cpus = 1
	publishDir "${params.folder_template}/result/${params.sample_type}/modified_PAS", mode: 'copy'

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

process FURHTER_FILTER_GENES{
	
	label "custom_python"
	label "short_time"
	label "middle_memory"
	cpus = params.basic_cores
	publishDir "${params.folder_template}/result/${params.sample_type}/common", mode: 'copy'

	input:
	path(genes)
	val(out)
	path(python_script)

	output:
	path("${out}")

	script:
	"""
	python3 ${python_script}\
	--genes ${genes}\
	--out ${out}\
	"""
}

process NON_OVERLAPPING{
	
	label "custom_python"
	label "short_time"
	label "middle_memory"
	cpus = params.basic_cores
	disk = '200 GB'
	publishDir "${params.folder_template}/result/${params.sample_type}/common", mode: 'copy'

	input:
	tuple path(genes), val(chromosome), val(direction)
	val(out)
	path(python_script)

	output:
	path("${out}_${chromosome}_${direction}.bed")

	script:
	"""
	python3 ${python_script}\
	--genes_dir ${genes}\
	--out ${out}\
	--chrom ${chromosome}\
	--direction ${direction}\
	--species ${params.sample_type}
	"""
}

process MERGE_ALL_DF{
	
	label "custom_python"
	label "short_time"
	label "middle_memory"
	cpus = params.basic_cores
	publishDir "${params.folder_template}/result/${params.sample_type}/non_overlap", mode: 'copy'

	input:
	path(beds)
	val(out1)
	val(out2)
	path(python_script)

	output:
	path("${out1}_${out2}.bed")

	script:
	"""
	python3 ${python_script}\
	in_bed ${beds}\
	--out ${out1}_${out2}.bed
	"""
}

process CHANGE_BED2{

	label "custom_python"
	label "middle_memory"
	cpus = 1
	publishDir "${params.folder_template}/result/${params.sample_type}/modified_PAS", mode: 'copy'

	input:
	path(clustered_bed_file)
	path (python_script)

	output:
	path("modified_v2_full_*.bed")

	script:
	"""
	python3 ${python_script} --bed_in ${clustered_bed_file}
	"""
}

process BED_INTERSECT_R{

	label "bedtools"
	label "middle_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/non_overlap", mode: 'copy'

	input:
	path(clustered_bed_file)
	// either terminal_exons, exons or genes.bed
	path(bed)
	val(intersect_out_name)
	
	output:
	path("${intersect_out_name}.bed")
	
	script:
	"""
	bedtools window -w 1 -sm -u -a ${clustered_bed_file} -b ${bed} > ${intersect_out_name}.bed
	"""
}

process BED_NO_INTERSECT_R{

	label "bedtools"
	label "middle_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/non_overlap", mode: 'copy'

	input:
	path(clustered_bed_file)
	// either terminal_exons, exons or genes.bed
	path(bed)
	val(no_intersect_out_name)
	
	output:
	path("${no_intersect_out_name}.bed")

	script:
	"""
	bedtools window -w 1 -sm -v -a ${clustered_bed_file} -b ${bed} > ${no_intersect_out_name}.bed
	"""
}

process BED_SUBTRACT{

	label "bedtools"
	label "middle_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/non_overlap", mode: 'copy'

	input:
	path(genes_bed)
	path(nonoverlapping_bed)
	val(out_name)
	
	output:
	path("${out_name}")
	
	script:
	"""
	bedtools subtract -a ${genes_bed} -b ${nonoverlapping_bed} -s > ${out_name}
	"""
}

process ASSIGN_OVERLAP_COLUMN{
	
	label "custom_python"
	label "short_time"
	label "middle_memory"
	cpus = params.basic_cores
	publishDir "${params.folder_template}/result/${params.sample_type}/non_overlap", mode: 'copy'

	input:
	path(pas_in_nonoverlapping)
	path(pas_in_overlapping)
	path(intergenic)
	val(out)
	path(python_script)

	output:
	path("${out}")

	script:
	"""
	python3 ${python_script}\
	--non_overlapping ${pas_in_nonoverlapping}\
	--overlapping ${pas_in_overlapping}\
	--intergenic ${intergenic}\
	--out ${out}
	"""
}

process INTERGENIC_GENE_ID{
	
	label "custom_python"
	label "long_time"
	label "middle_memory"
	cpus = params.super_heavy_cores
	publishDir "${params.folder_template}/result/${params.sample_type}/reassign", mode: 'copy'

	input:
	path(pas)
	path(original)
	path(te)
	val(out)
	path(python_script)

	output:
	path("${out}")

	script:
	"""
	python3 ${python_script}\
	--pas ${pas}\
	--original ${original}\
	--te ${te}\
	--out ${out}\
	--n ${params.super_heavy_cores}
	"""
}

process REASSIGN_GENEID{
	
	label "custom_python"
	label "short_time"
	label "heavy_memory"
	cpus = params.super_heavy_cores
	publishDir "${params.folder_template}/result/${params.sample_type}/reassign", mode: 'copy'

	input:
	tuple path(pas), val(chromosome), val(direction)
	path(original)
	path(genes)
	val(out)
	val(method)
	path(te)
	path(python_script)

	output:
	tuple path("${method}_${out}_${chromosome}_${direction}.bed"), val("${chromosome}"), val("${direction}")

	script:
	"""
	python3 ${python_script}\
	--pas ${pas}\
	--original ${original}\
	--genes ${genes}\
	--out ${out}\
	--chrom ${chromosome}\
	--strand ${direction}\
	--species ${params.sample_type}\
	--method ${method}\
	--te ${te}\
	"""
}

process REASSIGN_CLASSID{
	
	label "custom_python"
	label "short_time"
	label "middle_memory"
	cpus = params.heavy_cores
	publishDir "${params.folder_template}/result/${params.sample_type}/reassign", mode: 'copy'

	input:
	tuple path(pas), val(chromosome), val(direction)
	path(exons)
	path(te)
	path(introns)
	val(out)
	path(python_script)

	output:
	path("${out}*")

	script:
	"""
	python3 ${python_script}\
	--pas ${pas}\
	--exons ${exons}\
	--te ${te}\
	--introns ${introns}\
	--out ${out}\
	--chrom ${chromosome}\
	--strand ${direction}\
	--species ${params.sample_type}
	"""
}

process MERGE_AGAIN{
	
	label "custom_python"
	label "short_time"
	label "middle_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/reassign_merged", mode: 'copy'

	input:
	path(genic_pases)
	path(intergenic_pas)
	val(out)
	path(python_script)

	output:
	path("${out}")

	script:
	"""
	python3 ${python_script}\
	in_bed ${genic_pases}\
	--intergenic ${intergenic_pas}\
	--out ${out}
	"""
}

process COMPUTE_USAGE{
	
	label "custom_python"
	label "short_time"
	label "heavy_memory"
	cpus = params.heavy_cores
	publishDir "${params.folder_template}/result/${params.sample_type}/filtered_pas", mode: 'copy'

	input:
	path(pas)
	val(out)
	path(python_script)

	output:
	path("*${out}.bed")

	script:
	"""
	python3 ${python_script}\
	--pas ${pas}\
	--species ${params.sample_type}\
	--out ${out}\
	"""
}

process FILTER_PAS{
	
	label "custom_python"
	label "long_time"
	label "heavy_memory"
	cpus = params.heavy_cores
	publishDir "${params.folder_template}/result/${params.sample_type}/filtered_pas", mode: 'copy'

	input:
	path(pas)
	val(out)
	val(thres)
	path(python_script)

	output:
	path("*${out}.bed")
	path("*_filtering_summary_results.bed")

	script:
	"""
	python3 ${python_script}\
	in_bed ${pas}\
	--species ${params.sample_type}\
	--out ${out}\
	--thres ${thres}\
	"""
}

process MAKE_FULL_TABLE{
	
	label "custom_python"
	label "short_time"
	label "middle_memory"
	publishDir "${params.folder_template}/result/${params.sample_type}/filtered_pas/", mode: 'copy'

	input:
	path(subsets)
	path(original_pas)
	val(out)
	val(which)
	path(python_script)

	output:
	path("${params.sample_type}_${which}_${out}")

	script:
	"""
	python3 ${python_script}\
	in_bed ${subsets}\
	--original_pas ${original_pas}\
	--out ${params.sample_type}_${which}_${out}\
	"""
}