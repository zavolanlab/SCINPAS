#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {FIND_TERMINAL_EXONS} from './processes_all_samples'
include {FIND_EXONS_GENES_BED} from './processes_all_samples'
include {polyA_all_samples} from './polyA_workflow_advanced_all_samples'

workflow{
	
	if(params.check == "yes"){
		inputs = Channel
			.fromPath(params.samples_list)
			.splitCsv(header: true)
			.map {row -> tuple(row.dir, row.sample, row.organ)}
	}

	else if(params.check == "no"){
		inputs = Channel
			.fromPath(params.filtered_samples_list)
			.splitCsv(header: true)
			.map {row -> tuple(row.dir, row.sample, row.organ)}
	}

	our_terminal_exons = FIND_TERMINAL_EXONS(params.extended_annotation, params.find_terminal_exons_script)
	(our_exons, our_genes) = FIND_EXONS_GENES_BED(params.extended_annotation, params.find_exons_genes_script)		
	
	// To get different chromosomes according to the species.
	if(params.sample_type == "mouse"){

		chrs = Channel.fromList(params.mouse_chromosomes)
		// 1 ~ 19 + X + Y
		num_chromosomes = 21		
	}

	else if(params.sample_type == "human"){

		chrs = Channel.fromList(params.human_chromosomes)	
		// 1 ~ 22 + X + Y
		num_chromosomes = 24
	}

	else if(params.sample_type == "worm"){

		chrs = Channel.fromList(params.worm_chromosomes)	
		// "I" ~ "V" + "X"
		num_chromosomes = 6
	}

	polyA_all_samples(inputs, chrs, our_terminal_exons, our_exons, our_genes, num_chromosomes)	

}

