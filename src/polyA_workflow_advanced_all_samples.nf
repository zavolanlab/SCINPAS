#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//import and make alias. need to make alias because process can be invoked once and only once
include {PREPARE_IN_SLC_CATALOG; FILTER_MULTIMAPPING_CATALOG; CONCAT_CSVS_CATALOG; SPLIT_PHASE1_CATALOG; DEDUP_CATALOG; FASTQC_CATALOG; FASTQC_SWARM_PLOT_CATALOG;\
SPLIT_FILTERED_DEDUP_CATALOG; FIX_SOFTCLIPPED_REGION_CATALOG; MERGE_DEDUP_CATALOG; MERGE_POLYA_CATALOG; GET_COUNTS_CATALOG; GET_POLYA_UNIQUE_CLEAVAGE_SITES_CATALOG;\
MERGED_BED_CATALOG; SPLIT_BY_DIRECTION; GROUPBY_BED_CATALOG; PERFORM_CLUSTERING_CATALOG; LEFT_JOIN_CATALOG; FILTER_FURTHER_TERMINAL_EXONS; GET_MOTIF_FREQ_GTF_PLOT;\
COMPUTE_SCORES_FOR_MOTIF_ALTER_ALL_SAMPLES;} from './processes_advanced'

include {SORT_PHASE1_CATALOG as SORT_DEDUP} from './processes_advanced'
include {SORT_PHASE2_CATALOG as SORT_FIXED_DEDUP} from './processes_advanced'

include {GET_POLYA_CATALOG as GET_POLYA_SAMPLE_ONLY} from './processes_advanced'

include {SORT_PHASE2_CATALOG as SORT_POLYA_SAMPLE_ONLY} from './processes_advanced'
include {SPLIT_PHASE2_CATALOG as SPLIT_ALL_POLYA} from './processes_advanced'

// sample (classify reads into different classes)
include {BED_INTERSECT_CATALOG as GENES_INTERSECT_ALL_SAMPLES} from './processes_advanced'
include {BED_NO_INTERSECT_CATALOG as GENES_NO_INTERSECT_ALL_SAMPLES} from './processes_advanced'

include {BED_INTERSECT_CATALOG as EXONS_INTERSECT_ALL_SAMPLES} from './processes_advanced'
include {BED_NO_INTERSECT_CATALOG as EXONS_NO_INTERSECT_ALL_SAMPLES} from './processes_advanced'

include {BED_INTERSECT_CATALOG as TE_INTERSECT_ALL_SAMPLES} from './processes_advanced'
include {BED_NO_INTERSECT_CATALOG as TE_NO_INTERSECT_ALL_SAMPLES} from './processes_advanced'

include {ISOLATE_PA_T_CATALOG as ISOLATE_PA_T_ALL_SAMPLES} from './processes_advanced'

// frequency plot
include {PLOT_ATGC_CATALOG as PLOT_ATGC_ALL_POLYA} from './processes_advanced'
include {PLOT_ATGC_CATALOG as PLOT_ATGC_ABOVE} from './processes_advanced'
include {PLOT_ATGC_CATALOG as PLOT_ATGC_BELOW} from './processes_advanced'
include {PLOT_ATGC_CATALOG as PLOT_ATGC_INTRONIC} from './processes_advanced'
include {PLOT_ATGC_CATALOG as PLOT_ATGC_INTERGENIC} from './processes_advanced'
include {PLOT_ATGC_CATALOG as PLOT_ATGC_EXONIC} from './processes_advanced'

// motif frequency plot
include {GET_MOTIF_FREQ_PLOT_CATALOG as GET_MOTIF_FREQ_PLOT_ABOVE} from './processes_advanced'
include {GET_MOTIF_FREQ_PLOT_CATALOG as GET_MOTIF_FREQ_PLOT_BELOW} from './processes_advanced'
include {GET_MOTIF_FREQ_PLOT_CATALOG as GET_MOTIF_FREQ_PLOT_INTRONIC} from './processes_advanced'
include {GET_MOTIF_FREQ_PLOT_CATALOG as GET_MOTIF_FREQ_PLOT_INTERGENIC} from './processes_advanced'
include {GET_MOTIF_FREQ_PLOT_CATALOG as GET_MOTIF_FREQ_PLOT_EXONIC} from './processes_advanced'

// number of PAS
include {GET_NUM_POLYA_SITES_CATALOG as GET_NUM_POLYA_SITES_CATALOG_ALL_POLYA} from './processes_advanced'
include {GET_NUM_POLYA_SITES_CATALOG as GET_NUM_POLYA_SITES_CATALOG_ABOVE} from './processes_advanced'
include {GET_NUM_POLYA_SITES_CATALOG as GET_NUM_POLYA_SITES_CATALOG_BELOW} from './processes_advanced'
include {GET_NUM_POLYA_SITES_CATALOG as GET_NUM_POLYA_SITES_CATALOG_INTRONIC} from './processes_advanced'
include {GET_NUM_POLYA_SITES_CATALOG as GET_NUM_POLYA_SITES_CATALOG_INTERGENIC} from './processes_advanced'
include {GET_NUM_POLYA_SITES_CATALOG as GET_NUM_POLYA_SITES_CATALOG_EXONIC} from './processes_advanced'

workflow polyA_all_samples{

	take: 
	specific_dataFolders
	chromosomes
	terminal_exons_bed
	exons_bed
	genes_bed
	num_chrom

	main:

	if(params.sample_type == "mouse"){
		organs = Channel
			.fromPath(params.mouse_organs)
			.splitCsv(header: true)
			.map {row -> tuple(row.sample, row.organ)}
			// .view()		
	}

	else if(params.sample_type == "human"){
		organs = Channel
			.fromPath(params.human_organs)
			.splitCsv(header: true)
			.map {row -> tuple(row.sample, row.organ)}
			// .view()	
	}	

	//////////////////////////// trimming_and_deduplication ////////////////////////////

	full_data_tuple = PREPARE_IN_SLC_CATALOG(organs, params.sample_type)
	
	(mapq_filtered_full_data_tuple, readInfo_csvs) = FILTER_MULTIMAPPING_CATALOG(full_data_tuple)
	readInfo_merged_csv = CONCAT_CSVS_CATALOG(readInfo_csvs.collect, params.mapped_unmapped)
	// combine: cartesian product channel between full_data_tuple and chromosomes
	// possorted_bams_bais = SPLIT_PHASE1_CATALOG(mapq_filtered_full_data_tuple.combine(chromosomes))
	
	// // deduplication of "samples"
	// dedup_bams = DEDUP_CATALOG(possorted_bams_bais, params.dedup_script)

	// sorted_dedup_bams_bais = SORT_DEDUP(dedup_bams)
	
	// dedup_bams_bais_full = MERGE_DEDUP_CATALOG(sorted_dedup_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom))
	
	// fastqc_out_zips = FASTQC_CATALOG(dedup_bams_bais_full)
	// (seq_quality_csv, swarm_plot) = FASTQC_SWARM_PLOT_CATALOG(fastqc_out_zips.collect())

	// // seq_quality_csv is already a channel.
	// filtered_samples = seq_quality_csv
	// 	.splitCsv(header: true)
	// 	.map{row -> tuple(row.organ, row.sequencing_quality, row.scinpas_organ)}
	// 	.filter{row.sequencing_quality >= params.seq_qual_threshold}
	// 	.view()
	// n_samples = fastqc_out_zips.count()
	// n_samples.view()

	// // fix softclipped region (Alignment fixation) in the dedup file.
	// fixed_dedup_bams = FIX_SOFTCLIPPED_REGION_CATALOG(sorted_dedup_bams_bais, params.fix_softclipped_alter_script)

	// // sort fixed deduplicated bam in sample
	// fixed_dedup_sorted_bams_bais = SORT_FIXED_DEDUP(fixed_dedup_bams)
	
	// //////////////////////////// get_polyA ////////////////////////////
	// // sort_phase1, sort_phase2 can both do sorting and indexing (for only 1 file) and (for all files) because of tuple
	// // if for all files: (sorted_bam1, sorted_bam1, specific_sample6), (sorted_bam2, sorted_bam2, specific_sample6), (sorted_bam5, sorted_bam5, specific_sample6), (sorted_bam6, sorted_bam6, specific_sample6)
	// // if for only 1 file: (sorted_full_bam, sorted_full_bai, specific_sample5), (sorted_full_bam', sorted_full_bai', specific_sample6)

	// // polyA and non polyA reads within certain chromosome (for all samples)
	// (polyA_sample_only, non_polyA_sample_only) = GET_POLYA_SAMPLE_ONLY(fixed_dedup_sorted_bams_bais, params.out_polyA, params.out_non_polyA, params.getPolyA_advanced_alter_script)

	// // sort polyA reads within certain chromosome (for all samples)
	// polyA_sorted_bams_bais_samples = SORT_POLYA_SAMPLE_ONLY(polyA_sample_only)

	// // groupTuple allows to group file paths by samples.
	// // by: 2 means group by 3rd element of the tuple (i.e. samples)
	// // channel looks like: [[bam1, bam2...... bam6.....bamY], [bai1, bai2...... bai6.....baiY], sample1]
	// // dont need to use .collect()
	// // sorted_dedup_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom).view()
	// polyA_sorted_full_bams_bais_sample_only = MERGE_POLYA_CATALOG(polyA_sorted_bams_bais_samples.groupTuple(by: 2, sort: true, size: num_chrom), params.out_polyA)
	
	// polyA_sorted_full_bams_bais_counts_samples = GET_COUNTS_CATALOG(polyA_sorted_full_bams_bais_sample_only, params.counts_polyA, params.polyA_type, params.get_count_script)

	// polyA_sorted_bams_bais_counts_samples = SPLIT_ALL_POLYA(polyA_sorted_full_bams_bais_counts_samples.combine(chromosomes), params.out_polyA_partial)
	// // polyA_sorted_bams_bais_counts_samples.view()

	// polyA_unique_cs_beds_organs_chrom = GET_POLYA_UNIQUE_CLEAVAGE_SITES_CATALOG(polyA_sorted_bams_bais_counts_samples, params.polyA_unique_cs_bed_sample, 1, params.get_polyA_unique_cs_script)

	// // polyA_unique_cs_beds_sample.groupTuple(by:1, sort: true).view()

	// direction = Channel
	// 	.from(['+', '-'])
	
	// // polyA_unique_cs_beds_sample.combine(direction).view()

	// // combine bed, chromosome tuple with direction and then split bed file by direction.
	// (polyA_unique_cs_beds_chrom_direction, polyA_unique_cs_beds_chr_dir_organs) = SPLIT_BY_DIRECTION(polyA_unique_cs_beds_organs_chrom.combine(direction), params.split_by_direction_script)
	// // all samples bed files are grouped by its chromosome and direction
	// // sort removed in groupTuple for performance. 
	// // if you dont specify the size, it will wait until all inputs are recieved.
	// all_samples_polyA_unique_cs_beds_chrom_strand = GROUPBY_BED_CATALOG(polyA_unique_cs_beds_chrom_direction.groupTuple(by:[1,2]), params.grouped_by_out, n_samples, params.grouped_by_script)	

	// (all_polyA_clustered_beds, all_polyA_modified_unique_cs_beds) = PERFORM_CLUSTERING_CATALOG(all_samples_polyA_unique_cs_beds_chrom_strand, params.all_samples_cluster_out, params.modified_unique_cs, params.cluster_pas_alter_script)

	// // Grouping with same chromosome, direction and organ
	// // to generate organ specific PAS and organ specific modified_unique_cs_beds for a particular chromosome and direction.

	// // get() gets each element of the tuple/list. 0-based index
	// // tokenize('-') means split string by '-' so that you can get organ out of scinpas-organ by again using get(1)
	// // This will enable group bed files by organ.
	// polyA_unique_cs_beds_chr_dir_organs_improved = polyA_unique_cs_beds_chr_dir_organs.map{it->[it.get(0), it.get(1), it.get(2), it.get(3).tokenize('-').get(1)]}
	// // polyA_unique_cs_beds_chr_dir_organs_improved.groupTuple(by:[1,2,3]).view()

	// (organ_specific_cs, organ_specific_cs_filtered, all_PAS_split_by_organs) = LEFT_JOIN_CATALOG(polyA_unique_cs_beds_chr_dir_organs_improved.groupTuple(by:[1,2,3]), all_polyA_modified_unique_cs_beds.collect(), params.pas_split_by_organ_script)

	// // merge PAS bed files of a particular organ for all chromosomes and direction
	// all_PAS_split_by_organs_merged = MERGED_BED_CATALOG(all_PAS_split_by_organs.groupTuple(by:1), params.all_samples_cluster_out, params.merge_bed_script)
	
	// // classification of bed (sample)
	// polyA_intergenic_clustered_beds_all_samples = GENES_NO_INTERSECT_ALL_SAMPLES(all_PAS_split_by_organs_merged, genes_bed, params.intergenic_all_samples)
	// all_genic_clustered_beds_all_samples = GENES_INTERSECT_ALL_SAMPLES(all_PAS_split_by_organs_merged, genes_bed, params.all_genic_all_samples)

	// polyA_intronic_clustered_beds_all_samples = EXONS_NO_INTERSECT_ALL_SAMPLES(all_genic_clustered_beds_all_samples, exons_bed, params.intronic_all_samples)
	// all_exonic_clustered_beds_all_samples = EXONS_INTERSECT_ALL_SAMPLES(all_genic_clustered_beds_all_samples, exons_bed, params.all_exonic_all_samples)

	// polyA_exonic_clustered_beds_all_samples = TE_NO_INTERSECT_ALL_SAMPLES(all_exonic_clustered_beds_all_samples, terminal_exons_bed, params.exonic_all_samples)
	// polyA_TE_clustered_beds_all_samples = TE_INTERSECT_ALL_SAMPLES(all_exonic_clustered_beds_all_samples, terminal_exons_bed, params.te_all_samples)

	// (polyA_TE_below_beds_all_samples, polyA_TE_above_beds_all_samples) = ISOLATE_PA_T_ALL_SAMPLES(polyA_TE_clustered_beds_all_samples, terminal_exons_bed, params.isolate_pA_T_bed_script)	

	// //////////////////////////// analysis ////////////////////////////
	// if (params.analysis == "yes"){			
	// 	//////////////////////////// ATGC_plot ////////////////////////////

	// 	// plot frequency plot in each class
	// 	all_polyA_atgc_plts = PLOT_ATGC_ALL_POLYA(all_PAS_split_by_organs_merged, params.all_polyA_ATGC_out, params.plot_frequencies_script)

	// 	above_atgc_plts = PLOT_ATGC_ABOVE(polyA_TE_above_beds_all_samples, params.unannotated_ATGC_out, params.plot_frequencies_script)

	// 	below_atgc_plts = PLOT_ATGC_BELOW(polyA_TE_below_beds_all_samples, params.annotated_ATGC_out, params.plot_frequencies_script)

	// 	intergenic_atgc_plts = PLOT_ATGC_INTERGENIC(polyA_intergenic_clustered_beds_all_samples, params.intergenic_ATGC_out, params.plot_frequencies_script)

	// 	intronic_atgc_plts = PLOT_ATGC_INTRONIC(polyA_intronic_clustered_beds_all_samples, params.intronic_ATGC_out, params.plot_frequencies_script)

	// 	exonic_atgc_plt = PLOT_ATGC_EXONIC(polyA_exonic_clustered_beds_all_samples, params.polyA_exonic_ATGC_out, params.plot_frequencies_script)
	// 	//////////////////////////// motif_freq_plot ////////////////////////////
	// 	// draw 12 motives frequency plots based on gtf file
	// 	further_filtered_terminal_exons = FILTER_FURTHER_TERMINAL_EXONS(terminal_exons_bed, params.filter_terminal_exons_script)
	// 	(gtf_motif_graphs, gtf_motif_graphs_overlaid_png, gtf_motif_graphs_overlaid_svg, gtf_motif_order, gtf_peaks) = GET_MOTIF_FREQ_GTF_PLOT(further_filtered_terminal_exons, params.motif_search_out_gtf, params.gtf_search_motif_script)
		
	// 	// draw 12 motives frequency plots from our data in each class in each sample
	// 	(annotated_motif_graphs, annotated_motif_graphs_overlaid_png, annotated_motif_graphs_overlaid_svg, annotated_motif_orders, annotated_motif_scores) = GET_MOTIF_FREQ_PLOT_BELOW(polyA_TE_below_beds_all_samples,\
	// 	params.motif_search_out_template_annotated, params.below_annotated, gtf_peaks, params.search_motif_script)

	// 	(unannotated_motif_graphs, unannotated_motif_graphs_overlaid_png, unannotated_motif_graphs_overlaid_svg, unannotated_motif_orders, unannotated_motif_scores) = GET_MOTIF_FREQ_PLOT_ABOVE(polyA_TE_above_beds_all_samples,\
	// 	params.motif_search_out_template_unannotated, params.above_annotated, gtf_peaks, params.search_motif_script)

	// 	(intergenic_motif_graphs, intergenic_motif_graphs_overlaid_png, intergenic_motif_graphs_overlaid_svg, intergenic_motif_orders, intergenic_motif_scores) = GET_MOTIF_FREQ_PLOT_INTERGENIC(polyA_intergenic_clustered_beds_all_samples,\
	// 	params.motif_search_out_template_intergenic, params.intergenic_annotated, gtf_peaks, params.search_motif_script)	

	// 	(intronic_motif_graphs, intronic_motif_graphs_overlaid_png, intronic_motif_graphs_overlaid_svg, intronic_motif_orders, intronic_motif_scores) = GET_MOTIF_FREQ_PLOT_INTRONIC(polyA_intronic_clustered_beds_all_samples,\
	// 	params.motif_search_out_template_intronic, params.intronic_annotated, gtf_peaks, params.search_motif_script)	

	// 	(exonic_motif_graphs, exonic_motif_graphs_overlaid_png, exonic_motif_graphs_overlaid_svg, exonic_motif_orders, exonic_motif_scores) = GET_MOTIF_FREQ_PLOT_EXONIC(polyA_exonic_clustered_beds_all_samples,\
	// 	params.motif_search_out_template_exonic, params.polyA_exonic_annotated, gtf_peaks, params.search_motif_script)	
		
	// 	(polyA_motif_graphs_overlaid_png, polyA_motif_graphs_overlaid_svg, polyA_scores_csv) = COMPUTE_SCORES_FOR_MOTIF_ALTER_ALL_SAMPLES(all_PAS_split_by_organs_merged,\
	// 	params.motif_search_out_polyA, gtf_peaks, params.compute_motif_score_advanced_script)

	// 	// number of PAS in each class
	// 	num_polyA_csvs_all_pA = GET_NUM_POLYA_SITES_CATALOG_ALL_POLYA(all_PAS_split_by_organs_merged,\
	// 	params.allPAS_num_polyA_csv, params.all_polyA_annotated, params.get_numA_polyA_script)

	// 	num_polyA_csvs_above = GET_NUM_POLYA_SITES_CATALOG_ABOVE(polyA_TE_above_beds_all_samples,\
	// 	params.above_num_polyA_csv, params.above_annotated, params.get_numA_polyA_script)

	// 	num_polyA_csvs_below = GET_NUM_POLYA_SITES_CATALOG_BELOW(polyA_TE_below_beds_all_samples,\
	// 	params.below_num_polyA_csv, params.below_annotated, params.get_numA_polyA_script)

	// 	num_polyA_csvs_intronic = GET_NUM_POLYA_SITES_CATALOG_INTRONIC(polyA_intronic_clustered_beds_all_samples,\
	// 	params.intronic_num_polyA_csv, params.intronic_annotated, params.get_numA_polyA_script)

	// 	num_polyA_csvs_intergenic = GET_NUM_POLYA_SITES_CATALOG_INTERGENIC(polyA_intergenic_clustered_beds_all_samples,\
	// 	params.intergenic_num_polyA_csv, params.intergenic_annotated, params.get_numA_polyA_script)
	// }	
}