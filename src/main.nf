#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {PREPARE_SINGLE_LINKAGE} from './processes_advanced'
include {FIND_TERMINAL_EXONS} from './processes_advanced'
include {FIND_EXONS_GENES_BED} from './processes_advanced'
include {polyA} from './polyA_workflow_advanced'

workflow{
	(prepare_singleL, sl_script) = PREPARE_SINGLE_LINKAGE(params.cpp)

	// To get different inputs, chromosomes and negative controls
	// according to the species.
	if(params.sample_type == "mouse"){
		inputs = Channel
			.fromPath(params.mouse_input_dir)
			.map { file -> file.baseName }

		chrs = Channel.fromList(params.mouse_chromosomes)

		negative_control_dedup = Channel
					.fromPath(params.mouse_negative_control_dir)
					.map { file -> file.baseName }

		negative_control_raw = Channel
					.fromPath(params.mouse_negative_control_raw)
					.map { file -> file.baseName }
		
		// 1 ~ 19 + X + Y
		num_chromosomes = 21		
	}

	else if(params.sample_type == "human"){
		inputs = Channel
			.fromPath(params.human_input_dir)
			.map { file -> file.baseName }	

		chrs = Channel.fromList(params.human_chromosomes)	

		negative_control_dedup = Channel
				.fromPath(params.human_negative_control_dir)
				.map { file -> file.baseName }

		negative_control_raw = Channel
				.fromPath(params.human_negative_control_raw)
				.map { file -> file.baseName }
		
		// 1 ~ 22 + X + Y
		num_chromosomes = 24
	}
	
	// Get different terminal exons.bed, exons.bed and genes.bed according to the sample type. (This is done automatically behind the scene. no need to do if phrase)
	// Get terminal exons.bed
	terminal_exons = FIND_TERMINAL_EXONS(params.find_terminal_exons_script)
	// Get exons.bed and genes.bed
	(exons, genes) = FIND_EXONS_GENES_BED(params.find_exons_genes_script)

	polyA(inputs, prepare_singleL, sl_script, chrs, terminal_exons, exons, genes, negative_control_dedup, negative_control_raw, num_chromosomes)
}

