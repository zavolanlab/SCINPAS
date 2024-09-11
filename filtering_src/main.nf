#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {FIND_TERMINAL_EXONS; FIND_EXONS_GENES_BED; GET_INTRONIC_BED; CHANGE_BED; ADD_CLASS_COLUMN; GENE_ID;\
LEFT_JOIN_CATALOG; GET_ORGAN_SCORE; MERGE_ORGAN_SCORE_TO_PAS; RCS_MOTIF_CHECK; FURHTER_FILTER_GENES; NON_OVERLAPPING; CHANGE_BED2;\
MERGE_ALL_DF; BED_SUBTRACT; ASSIGN_OVERLAP_COLUMN; INTERGENIC_GENE_ID; REASSIGN_GENEID; REASSIGN_CLASSID; MERGE_AGAIN; COMPUTE_USAGE; FILTER_PAS; MAKE_FULL_TABLE;} from './process'

include {BED_INTERSECT_CATALOG as GENES_INTERSECT_ALL_SAMPLES} from './process'
include {BED_NO_INTERSECT_CATALOG as GENES_NO_INTERSECT_ALL_SAMPLES} from './process'

include {BED_INTERSECT_CATALOG as INTRONS_INTERSECT_ALL_SAMPLES} from './process'
include {BED_NO_INTERSECT_CATALOG as INTRONS_NO_INTERSECT_ALL_SAMPLES} from './process'

include {BED_INTERSECT_CATALOG as TE_INTERSECT_ALL_SAMPLES} from './process'
include {BED_NO_INTERSECT_CATALOG as TE_NO_INTERSECT_ALL_SAMPLES} from './process'

include {BED_INTERSECT_CATALOG_ANTISENSE as EXTRA_STEP_ONE_INTERSECT_ALL_SAMPLES} from './process'
include {BED_NO_INTERSECT_CATALOG_ANTISENSE as EXTRA_STEP_ONE_NO_INTERSECT_ALL_SAMPLES} from './process'

include {BED_INTERSECT_CATALOG_ANTISENSE as EXTRA_STEP_TWO_INTERSECT_ALL_SAMPLES} from './process'
include {BED_NO_INTERSECT_CATALOG_ANTISENSE as EXTRA_STEP_TWO_NO_INTERSECT_ALL_SAMPLES} from './process'

include {BED_INTERSECT_CATALOG_ANTISENSE as EXTRA_STEP_THREE_INTERSECT_ALL_SAMPLES} from './process'
include {BED_NO_INTERSECT_CATALOG_ANTISENSE as EXTRA_STEP_THREE_NO_INTERSECT_ALL_SAMPLES} from './process'

include {MERGE_ADD_COL_BED as MERGE_ADD_COL_BEDS} from './process'

include {BED_INTERSECT_R as BED_INTERSECT_NONOVERLAPPING_R} from './process'
include {BED_NO_INTERSECT_R as BED_NO_INTERSECT_NONOVERLAPPING_R} from './process'

include {BED_INTERSECT_R as BED_INTERSECT_OVERLAPPING_R} from './process'
include {BED_NO_INTERSECT_R as BED_NO_INTERSECT_OVERLAPPING_R} from './process'


workflow{

	our_terminal_exons = FIND_TERMINAL_EXONS(params.extended_annotation, params.find_terminal_exons_script)
	(our_exons, our_genes) = FIND_EXONS_GENES_BED(params.extended_annotation, params.find_exons_genes_script)	
	our_introns = GET_INTRONIC_BED(our_exons, params.introns_out, params.get_intronic_script)

	if(params.sample_type == "human"){

		chromosomes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y"]
	}

	else if(params.sample_type == "mouse"){

		chromosomes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,"X","Y"]
	}

	else if(params.sample_type == "worm"){

		chromosomes = ["I", "II", "III", "IV", "V", "X"]
	}

	inputs = Channel
		.fromPath(params.catalog)
		.splitCsv(header: true)
		.map {row -> tuple(row.chrom, row.direction, row.path)}

	modified_pas = CHANGE_BED(inputs, params.change_to_bed_script)

	// classification of intergenic PAS
	polyA_intergenic_clustered_beds_all_samples = GENES_NO_INTERSECT_ALL_SAMPLES(modified_pas, our_genes, params.intergenic_all_samples)
	all_genic_clustered_beds_all_samples = GENES_INTERSECT_ALL_SAMPLES(modified_pas, our_genes, params.all_genic_all_samples)
	
	// prioritize intronic pas
	all_exonic_clustered_beds_all_samples = INTRONS_NO_INTERSECT_ALL_SAMPLES(all_genic_clustered_beds_all_samples, our_introns, params.all_exonic_all_samples)
	polyA_intronic_clustered_beds_all_samples = INTRONS_INTERSECT_ALL_SAMPLES(all_genic_clustered_beds_all_samples, our_introns, params.intronic_all_samples)
	
	// classification of TE or exonic PAS
	polyA_exonic_clustered_beds_all_samples = TE_NO_INTERSECT_ALL_SAMPLES(all_exonic_clustered_beds_all_samples, our_terminal_exons, params.exonic_all_samples)
	polyA_TE_clustered_beds_all_samples = TE_INTERSECT_ALL_SAMPLES(all_exonic_clustered_beds_all_samples, our_terminal_exons, params.te_all_samples)

	/////////////////////////////////////////////////////////////////////////
	///////////////////// Further Classification of PAS /////////////////////
	polyA_true_intergenic_beds_all_samples = EXTRA_STEP_ONE_NO_INTERSECT_ALL_SAMPLES(polyA_intergenic_clustered_beds_all_samples, our_genes, params.true_intergenic_all_samples)
	all_antisense_genic_beds_all_samples = EXTRA_STEP_ONE_INTERSECT_ALL_SAMPLES(polyA_intergenic_clustered_beds_all_samples, our_genes, params.all_antisense_genic_all_samples)	

	polyA_antisense_intronic_beds_all_samples = EXTRA_STEP_TWO_NO_INTERSECT_ALL_SAMPLES(all_antisense_genic_beds_all_samples, our_exons, params.antisense_intronic_all_samples)
	all_antisense_exonic_beds_all_samples = EXTRA_STEP_TWO_INTERSECT_ALL_SAMPLES(all_antisense_genic_beds_all_samples, our_exons, params.all_antisense_exonic_all_samples)	

	polyA_antisense_exonic_beds_all_samples = EXTRA_STEP_THREE_NO_INTERSECT_ALL_SAMPLES(all_antisense_exonic_beds_all_samples, our_terminal_exons, params.antisense_exonic_all_samples)
	polyA_antisense_te_beds_all_samples = EXTRA_STEP_THREE_INTERSECT_ALL_SAMPLES(all_antisense_exonic_beds_all_samples, our_terminal_exons, params.antisense_te_all_samples)	

	all_classes_beds = modified_pas.mix(polyA_intronic_clustered_beds_all_samples, polyA_exonic_clustered_beds_all_samples, polyA_TE_clustered_beds_all_samples, polyA_true_intergenic_beds_all_samples,\
	polyA_antisense_intronic_beds_all_samples, polyA_antisense_exonic_beds_all_samples, polyA_antisense_te_beds_all_samples)

	additional_col_all_pas_no_tuple = ADD_CLASS_COLUMN(all_classes_beds.groupTuple(by: [0,1], sort: true), params.add_class_column_script)
	additional_col_all_pas_merged = MERGE_ADD_COL_BEDS(additional_col_all_pas_no_tuple.collect(), params.merge_add_col_bed)

	pas_w_geneid = GENE_ID(additional_col_all_pas_merged, our_genes, params.gene_id_pas_out, params.gene_id_alter_nextflow_script)

	polyA_unique_cs_beds_chr_dir_organs_improved = Channel
		.fromPath(params.file_list)
		.splitCsv(header: true)
		.map {row -> tuple(row.filename, row.chromosome, row.direction, row.organ)}

	all_polyA_modified_unique_cs_beds = Channel
		.fromPath(params.modified_cs_list)
		.splitCsv(header: true)
		.map {row -> row.sample}

	all_polyA_modified_unique_cs_beds.view()

	organ_specific_cs_filtered = LEFT_JOIN_CATALOG(polyA_unique_cs_beds_chr_dir_organs_improved.groupTuple(by:[1,2,3]), all_polyA_modified_unique_cs_beds.collect(), params.pas_split_by_organ_script)
	organ_specific_scores = GET_ORGAN_SCORE(organ_specific_cs_filtered.groupTuple(by:1), params.organ_score, params.organ_score_script)

	pas_w_gene_and_organ_score = MERGE_ORGAN_SCORE_TO_PAS(organ_specific_scores.collect(), pas_w_geneid, params.gene_id_pas_out, params.merge_pas_organ_score_script)
	pas_w_gene_w_organ_score_w_motif = RCS_MOTIF_CHECK(pas_w_gene_and_organ_score, params.gene_id_pas_out, params.rcs_motif_check_out, params.genome_fasta, params.rcs_motif_check_script)
	
	chrs = Channel.fromList(chromosomes)
	direction = Channel
		.from(['+', '-'])

	further_filtered_genes = FURHTER_FILTER_GENES(our_genes, params.further_filtered_genes_out, params.furhter_filter_genes_script)
	non_overlapping_regions = NON_OVERLAPPING(further_filtered_genes.combine(chrs).combine(direction), params.non_overlapping_genes, params.non_overlap_script)

	full_changed_bed = CHANGE_BED2(pas_w_gene_w_organ_score_w_motif, params.change_to_bed2_script)
	merged_non_overlapping_regions = MERGE_ALL_DF(non_overlapping_regions.collect(), params.merge_all_df_out, params.non_overlapping_genes, params.merge_all_df_script)
	merged_overlapping_regions = BED_SUBTRACT(further_filtered_genes, merged_non_overlapping_regions, params.overlapping_genes)
	
	pas_in_nonoverlapping_regions = BED_INTERSECT_NONOVERLAPPING_R(full_changed_bed, merged_non_overlapping_regions, params.overlapping_with_nonoverlapping_regions_out)
	pas_nonoverlapping_with_nonoverlapping_regions = BED_NO_INTERSECT_NONOVERLAPPING_R(full_changed_bed, merged_non_overlapping_regions, params.non_overlapping_with_nonoverlapping_regions_out)

	pas_in_overlapping_regions = BED_INTERSECT_OVERLAPPING_R(pas_nonoverlapping_with_nonoverlapping_regions, merged_overlapping_regions, params.overlapping_with_overlapping_regions_out)
	pas_in_intergenic_regions = BED_NO_INTERSECT_OVERLAPPING_R(pas_nonoverlapping_with_nonoverlapping_regions, merged_overlapping_regions, params.non_overapping_with_overlapping_regions_out)

	pas_assigned_overlap = ASSIGN_OVERLAP_COLUMN(pas_in_nonoverlapping_regions, pas_in_overlapping_regions, pas_in_intergenic_regions, params.assign_overlap_out, params.assign_overlap_script)

	intergenic_pas = INTERGENIC_GENE_ID(pas_assigned_overlap, pas_w_gene_w_organ_score_w_motif, our_terminal_exons, params.get_intergenic_gene_out, params.get_intergenic_gene_script)

	pas_reassigned_genes = REASSIGN_GENEID(pas_assigned_overlap.combine(chrs).combine(direction), pas_w_gene_w_organ_score_w_motif, further_filtered_genes, params.reassign_gene_out, "closest", our_terminal_exons, params.reassign_gene_script)
	pas_reassigned_genes_class = REASSIGN_CLASSID(pas_reassigned_genes, our_exons, our_terminal_exons, our_introns, params.reassign_class_out, params.reassign_class_script)

	// total_pas_reassigned = MERGE_AGAIN(pas_reassigned_genes_class.collect(), params.intergenic_pas, params.merge_again_out, params.merge_again_script)
	total_pas_reassigned = MERGE_AGAIN(pas_reassigned_genes_class.collect(), intergenic_pas, params.merge_again_out, params.merge_again_script)
	total_pas_w_usage = COMPUTE_USAGE(total_pas_reassigned, params.compute_usage_out, params.compute_usage_script)
	thresholds = Channel.fromList([10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 65, 70, 75, 80, 85, 90, 95])

	(filtered_pas, filtered_summary) = FILTER_PAS(total_pas_w_usage.collect(), params.filter_pas_out, thresholds, params.filter_pas_script)
	// (filtered_pas, filtered_summary) = FILTER_PAS(total_pas_reassigned, params.filter_pas_out, thresholds, params.filter_pas_script)
	

	full_filtered_table = MAKE_FULL_TABLE(filtered_pas.collect(), total_pas_reassigned, params.make_full_out, "closest", params.make_full_table_script)
}

