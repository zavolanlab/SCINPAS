#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// import and make alias. need to make alias because process can be invoked once and only once
include {TRIM_NEG_CONTROL; SPLIT_PHASE1; DEDUP; MERGE; UNCOUPLE_CSV_WITH_SAMPLE; FILTER_FURTHER_TERMINAL_EXONS;\
SPLIT_PER_SAMPLE; ONE_BY_ONE; GET_NON_SOFTCLIP; GET_A_FREQ_IN_NONSOFTCLIP; MAKE_COUNT_SAMPLE; GET_BAR_CHART;\
GET_SPAN_SPLIT; GET_SPAN_HIST; GET_D_HISTOGRAM_ALTER; FIX_SOFTCLIPPED_REGION; MERGE_CONTROL;\
COMBINE_MOTIF_ORDERS; GET_MOTIF_FREQ_GTF_PLOT; SPLIT_MOTIF_ORDER_PER_SAMPLE; MOTIF_ORDERS_ONE_BY_ONE; COMPUTE_SCORES_FOR_MOTIF_ALTER;\
GET_BAR_CHART_SCORE_MOTIF_POLYA; GET_SOFTCLIPPED_DISTRIBUTION; MERGE_ALL_NUM_PAS_CSVS; FUSE_MOTIF_SCORE_NUM_PAS; GET_BAR_CHART_NUM_POLYA; GET_TOTAL_NUM_GENES;\
GET_BAR_CHART_GENE_COVERAGE; SOFTCLIP_PWM_LOGO; SPLIT_CSV_BY_CELL_TYPE; SPLIT_SAMPLE_BY_TYPE; MERGE_BY_CELL_TYPE; DISTANCE_SCATTER; T1_T2_SCATTER;\
UNCOUPLE_BAM_WITH_SAMPLE; MERGE_ALL_SAMPLES; GET_POLYA_UNIQUE_CLEAVAGE_SITES_ALL_SAMPLES;\
PERFORM_CLUSTERING_ALL_SAMPLES; GET_OVERLAP; GET_TOP_OVERLAP;} from './processes_advanced'

include {PREPARE_IN_SLC as PREPARE_IN_SLC_DATA} from './processes_advanced'
include {PREPARE_IN_SLC as PREPARE_IN_SLC_CONTROL_DEDUP} from './processes_advanced'
include {PREPARE_IN_SLC as PREPARE_IN_SLC_CONTROL_RAW} from './processes_advanced'

include {GET_SPAN_RAW as GET_SPAN_RAW_SAMPLE} from './processes_advanced'
include {GET_SPAN_RAW as GET_SPAN_RAW_CONTROL} from './processes_advanced'

include {GET_SPAN_DEDUP as GET_SPAN_DEDUP_SAMPLE} from './processes_advanced'
include {GET_SPAN_DEDUP as GET_SPAN_DEDUP_CONTROL} from './processes_advanced'

include {GET_POLYA as GET_POLYA_SAMPLE_ONLY} from './processes_advanced'
include {GET_POLYA as GET_POLYA_NEG_CONTROL} from './processes_advanced'
include {GET_POLYA_CELL_TYPE as GET_POLYA_TYPE1_SPECIFIC} from './processes_advanced'
include {GET_POLYA_CELL_TYPE as GET_POLYA_TYPE2_SPECIFIC} from './processes_advanced'

include {CONCAT_SPAN as CONCAT_SPAN_RAW} from './processes_advanced'
include {CONCAT_SPAN as CONCAT_SPAN_DEDUP} from './processes_advanced'
include {CONCAT_SPAN as CONCAT_SPAN_SPLIT} from './processes_advanced'
include {GET_SPAN_HIST as GET_SPAN_HIST_RAW} from './processes_advanced'
include {GET_SPAN_HIST as GET_SPAN_HIST_DEDUP} from './processes_advanced'
include {GET_SPAN_HIST as GET_SPAN_HIST_SPLIT} from './processes_advanced'

include {SORT_PHASE1 as SORT_DEDUP} from './processes_advanced'
include {SORT_PHASE1 as SORT_SPLITTED} from './processes_advanced'
include {SORT_PHASE2 as SORT_FIXED_DEDUP} from './processes_advanced'

include {SORT_PHASE2 as SORT_POLYA_NEG_CONTROL} from './processes_advanced'
include {SORT_PHASE2 as SORT_POLYA_SAMPLE_ONLY} from './processes_advanced'
include {SORT_PHASE2 as SORT_TRIM_NEG_CONTROL} from './processes_advanced'

include {SORT_PHASE2 as SORT_NONPOLYA} from './processes_advanced'
include {SORT_PHASE2 as SORT_NONPOLYA_CONTROL} from './processes_advanced'

// polyA clustering (sample) 
include {GET_POLYA_UNIQUE_CLEAVAGE_SITES_BED as GET_POLYA_UNIQUE_CLEAVAGE_SITES_BED_SAMPLE} from './processes_advanced'
include {PERFORM_CLUSTERING_TUPLE as PERFORM_CLUSTERING_SAMPLE_TUPLE} from './processes_advanced'

// non-polyA clustering (sample)
include {GET_POLYA_UNIQUE_CLEAVAGE_SITES_BED as GET_NON_POLYA_UNIQUE_CLEAVAGE_SITES_BED_SAMPLE} from './processes_advanced'
include {PERFORM_CLUSTERING_TUPLE as PERFORM_CLUSTERING_SAMPLE_TUPLE_NONPOLYA} from './processes_advanced'

// polyA clustering (control)
include {GET_POLYA_UNIQUE_CLEAVAGE_SITES_BED as GET_POLYA_UNIQUE_CLEAVAGE_SITES_BED_CONTROL} from './processes_advanced'
include {PERFORM_CLUSTERING_TUPLE as PERFORM_CLUSTERING_CONTROL_TUPLE} from './processes_advanced'

// dedup clustering (control)
include {GET_POLYA_UNIQUE_CLEAVAGE_SITES_BED as GET_DEDUP_UNIQUE_CLEAVAGE_SITES_BED_CONTROL} from './processes_advanced'
include {PERFORM_CLUSTERING_TUPLE as PERFORM_CLUSTERING_DEDUP_CONTROL_TUPLE} from './processes_advanced'

// non-polyA clustering (control)
include {GET_POLYA_UNIQUE_CLEAVAGE_SITES_BED as GET_NON_POLYA_UNIQUE_CLEAVAGE_SITES_BED_CONTROL} from './processes_advanced'
include {PERFORM_CLUSTERING_TUPLE as PERFORM_CLUSTERING_CONTROL_TUPLE_NONPOLYA} from './processes_advanced'

// sample (classify reads into different classes)
include {BED_INTERSECT_TUPLE as GENES_INTERSECT_TUPLE} from './processes_advanced'
include {BED_NO_INTERSECT_TUPLE as GENES_NO_INTERSECT_TUPLE} from './processes_advanced'
include {BED_INTERSECT_TUPLE as EXONS_INTERSECT_TUPLE} from './processes_advanced'
include {BED_NO_INTERSECT_TUPLE as EXONS_NO_INTERSECT_TUPLE} from './processes_advanced'
include {BED_INTERSECT_TUPLE as TE_INTERSECT_TUPLE} from './processes_advanced'
include {BED_NO_INTERSECT_TUPLE as TE_NO_INTERSECT_TUPLE} from './processes_advanced'

include {BED_INTERSECT_TUPLE as TE_FURTHER_INTERSECT_TUPLE} from './processes_advanced'

include {SPLIT_PHASE2 as SPLIT_ALL_POLYA} from './processes_advanced'
include {CHANGE_TUPLE_ORDER as CHANGE_TUPLE_ORDER_POLYA_SAMPLE} from './processes_advanced'
include {CLASSIFY_BAM as CLASSIFY_TE_SAMPLE} from './processes_advanced'
include {CLASSIFY_BAM as CLASSIFY_E_SAMPLE} from './processes_advanced'
include {CLASSIFY_BAM as CLASSIFY_INTRONIC_SAMPLE} from './processes_advanced'
include {CLASSIFY_BAM as CLASSIFY_INTERGENIC_SAMPLE} from './processes_advanced'
include {CLASSIFY_BELOW_ABOVE as CLASSIFY_BELOW_SAMPLE} from './processes_advanced'
include {CLASSIFY_BELOW_ABOVE as CLASSIFY_ABOVE_SAMPLE} from './processes_advanced'

include {SORT_CLASSIFIED as SORT_CLASSIFIED_TE} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_CLASSIFIED_NON_TE} from './processes_advanced'

include {SORT_CLASSIFIED as SORT_CLASSIFIED_EXONIC} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_CLASSIFIED_NON_TE_NON_E} from './processes_advanced'

include {SORT_CLASSIFIED as SORT_CLASSIFIED_INTRONIC} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_CLASSIFIED_NON_TE_NON_E_NON_I} from './processes_advanced'

include {SORT_CLASSIFIED as SORT_CLASSIFIED_INTERGENIC} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_CLASSIFIED_NON_TE_NON_E_NON_I_NON_IG} from './processes_advanced'

include {SORT_CLASSIFIED as SORT_CLASSIFIED_ANNOTATED} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_CLASSIFIED_NON_ANNOTATED} from './processes_advanced'

include {SORT_CLASSIFIED as SORT_CLASSIFIED_UNANNOTATED} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_CLASSIFIED_NON_ANNOTATED_UNANNOTATED} from './processes_advanced'

include {MERGE_CLASSIFIED as MERGE_CLASSIFIED_TE} from './processes_advanced'
include {MERGE_CLASSIFIED as MERGE_CLASSIFIED_E} from './processes_advanced'
include {MERGE_CLASSIFIED as MERGE_CLASSIFIED_INTRONIC} from './processes_advanced'
include {MERGE_CLASSIFIED as MERGE_CLASSIFIED_INTERGENIC} from './processes_advanced'
include {MERGE_CLASSIFIED as MERGE_CLASSIFIED_ANNOTATED} from './processes_advanced'
include {MERGE_CLASSIFIED as MERGE_CLASSIFIED_UNANNOTATED} from './processes_advanced'

include {CONVERT_TUPLE as CONVERT_TUPLE_POLYAT_SAMPLE} from './processes_advanced'
include {ISOLATE_PA_T as ISOLATE_PA_T_SAMPLE} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_ABOVE} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_BELOW} from './processes_advanced'

// control (classify reads into different classes)
include {BED_INTERSECT_TUPLE as GENES_INTERSECT_CONTROL_TUPLE} from './processes_advanced'
include {BED_NO_INTERSECT_TUPLE as GENES_NO_INTERSECT_CONTROL_TUPLE} from './processes_advanced'
include {BED_INTERSECT_TUPLE as EXONS_INTERSECT_CONTROL_TUPLE} from './processes_advanced'
include {BED_NO_INTERSECT_TUPLE as EXONS_NO_INTERSECT_CONTROL_TUPLE} from './processes_advanced'
include {BED_INTERSECT_TUPLE as TE_INTERSECT_CONTROL_TUPLE} from './processes_advanced'
include {BED_NO_INTERSECT_TUPLE as TE_NO_INTERSECT_CONTROL_TUPLE} from './processes_advanced'

include {SPLIT_PHASE2 as SPLIT_ALL_POLYA_CONTROL} from './processes_advanced'
include {CHANGE_TUPLE_ORDER as CHANGE_TUPLE_ORDER_POLYA_CONTROL} from './processes_advanced'
include {CLASSIFY_BAM as CLASSIFY_TE_CONTROL} from './processes_advanced'
include {CLASSIFY_BAM as CLASSIFY_E_CONTROL} from './processes_advanced'
include {CLASSIFY_BAM as CLASSIFY_INTRONIC_CONTROL} from './processes_advanced'
include {CLASSIFY_BAM as CLASSIFY_INTERGENIC_CONTROL} from './processes_advanced'
include {CLASSIFY_BELOW_ABOVE as CLASSIFY_BELOW_CONTROL} from './processes_advanced'
include {CLASSIFY_BELOW_ABOVE as CLASSIFY_ABOVE_CONTROL} from './processes_advanced'

include {SORT_CLASSIFIED as SORT_CLASSIFIED_TE_CONTROL} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_CLASSIFIED_NON_TE_CONTROL} from './processes_advanced'

include {SORT_CLASSIFIED as SORT_CLASSIFIED_EXONIC_CONTROL} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_CLASSIFIED_NON_TE_NON_E_CONTROL} from './processes_advanced'

include {SORT_CLASSIFIED as SORT_CLASSIFIED_INTRONIC_CONTROL} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_CLASSIFIED_NON_TE_NON_E_NON_I_CONTROL} from './processes_advanced'

include {SORT_CLASSIFIED as SORT_CLASSIFIED_INTERGENIC_CONTROL} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_CLASSIFIED_NON_TE_NON_E_NON_I_NON_IG_CONTROL} from './processes_advanced'

include {SORT_CLASSIFIED as SORT_CLASSIFIED_ANNOTATED_CONTROL} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_CLASSIFIED_NON_ANNOTATED_CONTROL} from './processes_advanced'

include {SORT_CLASSIFIED as SORT_CLASSIFIED_UNANNOTATED_CONTROL} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_CLASSIFIED_NON_ANNOTATED_UNANNOTATED_CONTROL} from './processes_advanced'

include {MERGE_CLASSIFIED as MERGE_CLASSIFIED_TE_CONTROL} from './processes_advanced'
include {MERGE_CLASSIFIED as MERGE_CLASSIFIED_E_CONTROL} from './processes_advanced'
include {MERGE_CLASSIFIED as MERGE_CLASSIFIED_INTRONIC_CONTROL} from './processes_advanced'
include {MERGE_CLASSIFIED as MERGE_CLASSIFIED_INTERGENIC_CONTROL} from './processes_advanced'
include {MERGE_CLASSIFIED as MERGE_CLASSIFIED_ANNOTATED_CONTROL} from './processes_advanced'
include {MERGE_CLASSIFIED as MERGE_CLASSIFIED_UNANNOTATED_CONTROL} from './processes_advanced'

include {CONVERT_TUPLE as CONVERT_TUPLE_POLYAT_CONTROL} from './processes_advanced'
include {ISOLATE_PA_T as ISOLATE_PA_T_CONTROL} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_ABOVE_CONTROL} from './processes_advanced'
include {SORT_CLASSIFIED as SORT_BELOW_CONTROL} from './processes_advanced'

include {SORT_PHASE3 as SORT_NO_SOFTCLIP} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_RAW_CONTROL} from './processes_advanced'
include {GET_COUNTS as GET_COUNTS_DEDUP_CONTROL} from './processes_advanced'
include {GET_COUNTS as GET_COUNTS_SOFTCLIPPED_CONTROL} from './processes_advanced'
include {GET_COUNTS as GET_COUNTS_NONPOLYA_CONTROL} from './processes_advanced'
include {GET_COUNTS as GET_COUNTS_POLYA_CONTROL} from './processes_advanced'
include {GET_COUNTS as GET_COUNTS_POLYA_T_CONTROL} from './processes_advanced'
include {GET_COUNTS as GET_COUNTS_ANNOTATED_CONTROL} from './processes_advanced'
include {GET_COUNTS as GET_COUNTS_UNANNOTATED_CONTROL} from './processes_advanced'
include {GET_COUNTS as GET_COUNTS_INTRONIC_CONTROL} from './processes_advanced'
include {GET_COUNTS as GET_COUNTS_INTERGENIC_CONTROL} from './processes_advanced'
include {GET_COUNTS as GET_COUNTS_EXONIC_CONTROL} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_RAW} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_RAW} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_DEDUP} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_DEDUP} from './processes_advanced'

include {GET_COUNTS as GET_SOFTCLIPPED_COUNTS_DEDUP} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_SOFTCLIPPED_COUNTS_DEDUP} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_NONPOLYA} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_NONPOLYA} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_INTERNAL_P} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_INTERNAL_P} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_POLYA} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_POLYA} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_POLYA_TERMINAL} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_POLYA_TERMINAL} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_ANNOTATED} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_ANNOTATED} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_UNANNOTATED} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_UNANNOTATED} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_INTRONIC} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_INTRONIC} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_INTERGENIC} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_INTERGENIC} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_EXONIC} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_EXONIC} from './processes_advanced'

include {MERGE_ALL_COUNTS as MERGE_ALL_COUNTS_SAMPLE} from './processes_advanced'
include {MERGE_ALL_COUNTS as MERGE_ALL_COUNTS_CONTROL} from './processes_advanced'

include {GET_DISTANCE_ALTER as GET_DISTANCE_ALTER_SAMPLE} from './processes_advanced'
include {GET_DISTANCE_ALTER as GET_DISTANCE_ALTER_CONTROL} from './processes_advanced'

include {SPLIT_PHASE2 as SPLIT_NEG_CONTROL_RAW} from './processes_advanced'
include {SPLIT_PHASE2 as SPLIT_NEG_CONTROL_DEDUP} from './processes_advanced'

include {PLOT_ATGC as PLOT_ATGC_ALL_POLYA} from './processes_advanced'
include {PLOT_ATGC as PLOT_ATGC_ABOVE} from './processes_advanced'
include {PLOT_ATGC as PLOT_ATGC_BELOW} from './processes_advanced'
include {PLOT_ATGC as PLOT_ATGC_NONPOLYA} from './processes_advanced'
include {PLOT_ATGC as PLOT_ATGC_INTRONIC} from './processes_advanced'
include {PLOT_ATGC as PLOT_ATGC_INTERGENIC} from './processes_advanced'
include {PLOT_ATGC as PLOT_ATGC_EXONIC} from './processes_advanced'
include {PLOT_ATGC as PLOT_ATGC_CONTROL} from './processes_advanced'
include {PLOT_ATGC as PLOT_ATGC_NOSOFTCLIP} from './processes_advanced'

include {GET_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_ABOVE} from './processes_advanced'
include {GET_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_BELOW} from './processes_advanced'
include {GET_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_INTRONIC} from './processes_advanced'
include {GET_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_INTERGENIC} from './processes_advanced'
include {GET_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_EXONIC} from './processes_advanced'

include {CONCAT_MOTIF_ORDERS as CONCAT_MOTIF_ORDERS_ANNOTATED} from './processes_advanced'
include {CONCAT_MOTIF_ORDERS as CONCAT_MOTIF_ORDERS_UNANNOTATED} from './processes_advanced'
include {CONCAT_MOTIF_ORDERS as CONCAT_MOTIF_ORDERS_INTERGENIC} from './processes_advanced'
include {CONCAT_MOTIF_ORDERS as CONCAT_MOTIF_ORDERS_INTRONIC} from './processes_advanced'
include {CONCAT_MOTIF_ORDERS as CONCAT_MOTIF_ORDERS_EXONIC} from './processes_advanced'

include {MERGE_ALL_MOTIF_CSVS as MERGE_ALL_MOTIF_ORDERS} from './processes_advanced'
include {MERGE_ALL_MOTIF_CSVS as MERGE_ALL_MOTIF_SCORES} from './processes_advanced'

include {GET_NUM_POLYA_SITES  as GET_NUM_POLYA_SITES_ALLPOLYA} from './processes_advanced'
include {GET_NUM_POLYA_SITES as GET_NUM_POLYA_SITES_ABOVE} from './processes_advanced'
include {GET_NUM_POLYA_SITES as GET_NUM_POLYA_SITES_BELOW} from './processes_advanced'
include {GET_NUM_POLYA_SITES as GET_NUM_POLYA_SITES_INTRONIC} from './processes_advanced'
include {GET_NUM_POLYA_SITES as GET_NUM_POLYA_SITES_INTERGENIC} from './processes_advanced'
include {GET_NUM_POLYA_SITES as GET_NUM_POLYA_SITES_EXONIC} from './processes_advanced'
include {GET_NUM_POLYA_SITES as GET_NUM_POLYA_SITES_CONTROL} from './processes_advanced'

include {GET_GENE_COVERAGE as GET_GENE_COVERAGE_DEDUP} from './processes_advanced'
include {GET_GENE_COVERAGE as GET_GENE_COVERAGE_POLYA} from './processes_advanced'

include {COMBINE_CSVS as COMBINE_SCORES_FOR_MOTIF} from './processes_advanced'

include {COMBINE_CSVS as COMBINE_NUM_POLYA_ALL_PAS} from './processes_advanced'
include {COMBINE_CSVS as COMBINE_NUM_POLYA_ABOVE} from './processes_advanced'
include {COMBINE_CSVS as COMBINE_NUM_POLYA_BELOW} from './processes_advanced'
include {COMBINE_CSVS as COMBINE_NUM_POLYA_INTRONIC} from './processes_advanced'
include {COMBINE_CSVS as COMBINE_NUM_POLYA_INTERGENIC} from './processes_advanced'
include {COMBINE_CSVS as COMBINE_NUM_POLYA_EXONIC} from './processes_advanced'
include {COMBINE_CSVS as COMBINE_NUM_POLYA_CONTROL} from './processes_advanced'

include {COMBINE_CSVS as COMBINE_GENE_COVERAGE_DEDUP} from './processes_advanced'
include {COMBINE_CSVS as COMBINE_GENE_COVERAGE_POLYA} from './processes_advanced'

include {SOFTCLIP_PWM_LOGO as SOFTCLIP_PWM_LOGO_LENGTH_FIVE} from './processes_advanced'
include {SOFTCLIP_PWM_LOGO as SOFTCLIP_PWM_LOGO_LENGTH_SIX} from './processes_advanced'
include {SOFTCLIP_PWM_LOGO as SOFTCLIP_PWM_LOGO_LENGTH_SEVEN} from './processes_advanced'
include {SOFTCLIP_PWM_LOGO as SOFTCLIP_PWM_LOGO_LENGTH_EIGHT} from './processes_advanced'
include {SOFTCLIP_PWM_LOGO as SOFTCLIP_PWM_LOGO_LENGTH_NINE} from './processes_advanced'
include {SOFTCLIP_PWM_LOGO as SOFTCLIP_PWM_LOGO_LENGTH_TEN} from './processes_advanced'

include {SORT_CELL_TYPE as SORT_POLYA_TYPE1_SPECIFIC} from './processes_advanced'
include {SORT_CELL_TYPE as SORT_POLYA_TYPE2_SPECIFIC} from './processes_advanced'

// cell type1 specific
include {GET_CELL_TYPE_POLYA_UNIQUE_CS_BED as POLYA_TYPE1_BAM_TO_BED} from './processes_advanced'
include {PERFORM_CLUSTERING_CELL_TYPE as PERFORM_CLUSTERING_CELL_TYPE1} from './processes_advanced'

include {BED_NO_INTERSECT_CELL_TYPE as GENE_NO_INTERSECT_TYPE1} from './processes_advanced'
include {BED_INTERSECT_CELL_TYPE as GENE_INTERSECT_TYPE1} from './processes_advanced'
include {BED_NO_INTERSECT_CELL_TYPE as EXON_NO_INTERSECT_TYPE1} from './processes_advanced'
include {BED_INTERSECT_CELL_TYPE as EXON_INTERSECT_TYPE1} from './processes_advanced'
include {BED_NO_INTERSECT_CELL_TYPE as TE_NO_INTERSECT_TYPE1} from './processes_advanced'
include {BED_INTERSECT_CELL_TYPE as TE_INTERSECT_TYPE1} from './processes_advanced'

include {BED_INTERSECT_CELL_TYPE as TE_FURTHER_INTERSECT_TYPE1} from './processes_advanced'

include {GET_CELL_TYPE_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_INTERGENIC_TYPE1} from './processes_advanced'
include {GET_CELL_TYPE_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_INTRONIC_TYPE1} from './processes_advanced'
include {GET_CELL_TYPE_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_EXONIC_TYPE1} from './processes_advanced'
include {GET_CELL_TYPE_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_TE_TYPE1} from './processes_advanced'

include {GET_CELL_TYPE_NUM_POLYA_SITES as GET_TYPE1_NUM_POLYA_SITES_ALLPAS} from './processes_advanced'
include {GET_CELL_TYPE_NUM_POLYA_SITES as GET_TYPE1_NUM_POLYA_SITES_INTRONIC} from './processes_advanced'
include {GET_CELL_TYPE_NUM_POLYA_SITES as GET_TYPE1_NUM_POLYA_SITES_EXONIC} from './processes_advanced'
include {GET_CELL_TYPE_NUM_POLYA_SITES as GET_TYPE1_NUM_POLYA_SITES_TE} from './processes_advanced'

// cell type2 specific
include {GET_CELL_TYPE_POLYA_UNIQUE_CS_BED as POLYA_TYPE2_BAM_TO_BED} from './processes_advanced'
include {PERFORM_CLUSTERING_CELL_TYPE as PERFORM_CLUSTERING_CELL_TYPE2} from './processes_advanced'

include {BED_NO_INTERSECT_CELL_TYPE as GENE_NO_INTERSECT_TYPE2} from './processes_advanced'
include {BED_INTERSECT_CELL_TYPE as GENE_INTERSECT_TYPE2} from './processes_advanced'
include {BED_NO_INTERSECT_CELL_TYPE as EXON_NO_INTERSECT_TYPE2} from './processes_advanced'
include {BED_INTERSECT_CELL_TYPE as EXON_INTERSECT_TYPE2} from './processes_advanced'
include {BED_NO_INTERSECT_CELL_TYPE as TE_NO_INTERSECT_TYPE2} from './processes_advanced'
include {BED_INTERSECT_CELL_TYPE as TE_INTERSECT_TYPE2} from './processes_advanced'

include {BED_INTERSECT_CELL_TYPE as TE_FURTHER_INTERSECT_TYPE2} from './processes_advanced'

include {GET_CELL_TYPE_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_INTERGENIC_TYPE2} from './processes_advanced'
include {GET_CELL_TYPE_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_INTRONIC_TYPE2} from './processes_advanced'
include {GET_CELL_TYPE_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_EXONIC_TYPE2} from './processes_advanced'
include {GET_CELL_TYPE_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_TE_TYPE2} from './processes_advanced'

include {GET_CELL_TYPE_NUM_POLYA_SITES as GET_TYPE2_NUM_POLYA_SITES_ALLPAS} from './processes_advanced'
include {GET_CELL_TYPE_NUM_POLYA_SITES as GET_TYPE2_NUM_POLYA_SITES_INTRONIC} from './processes_advanced'
include {GET_CELL_TYPE_NUM_POLYA_SITES as GET_TYPE2_NUM_POLYA_SITES_EXONIC} from './processes_advanced'
include {GET_CELL_TYPE_NUM_POLYA_SITES as GET_TYPE2_NUM_POLYA_SITES_TE} from './processes_advanced'

workflow polyA{

	take: 
	specific_dataFolders
	chromosomes
	terminal_exons_bed
	exons_bed
	genes_bed
	negative_cntrl_dedup
	negative_cntrl_raw
	num_chrom
	
	main:

	//chromosomes.view()
	//////////////////////////// trimming_and_deduplication ////////////////////////////

	// prepare data
	full_data_tuple = PREPARE_IN_SLC_DATA(specific_dataFolders, params.sample_type)

	// prepare folders for negative control
	negative_control_data_raw_tuple = PREPARE_IN_SLC_CONTROL_RAW(negative_cntrl_raw, params.sample_type)
	negative_control_data_dedup_tuple = PREPARE_IN_SLC_CONTROL_DEDUP(negative_cntrl_dedup, params.sample_type)

	// combine: cartesian product channel between full_data_tuple and chromosomes
	(possorted_bams_bais, bam_folders) = SPLIT_PHASE1(full_data_tuple.combine(chromosomes))
	
	// deduplication of "samples"
	(dedup_bams, splitted_bams) = DEDUP(bam_folders, params.dedup_script)

	sorted_dedup_bams_bais = SORT_DEDUP(dedup_bams)
		
	splitted_sorted_bams_bais = SORT_SPLITTED(splitted_bams)

	// groupTuple allows to group file paths by samples.
	// by: 2 means group by 3rd element of the tuple (i.e. samples)
	// channel looks like: [[bam1, bam2...... bam6.....bamY], [bai1, bai2...... bai6.....baiY], sample1]
	// dont need to use .collect()
	//sorted_dedup_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom).view()
	dedup_full_sorted_tuple = MERGE(sorted_dedup_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom))
	
	// fix softclipped region (Alignment fixation) in the dedup file.
	(fixed_dedup_full_bams, num_fixed_unfixed) = FIX_SOFTCLIPPED_REGION(dedup_full_sorted_tuple, params.fix_softclipped_script)
	// // add 'FC' tag to negative control.
	fixed_negative_control_dedup_bam = TRIM_NEG_CONTROL(negative_control_data_dedup_tuple, params.trim_neg_control_script)

	// sort fixed deduplicated bam in sample
	fixed_dedup_sorted_full_bams_bais = SORT_FIXED_DEDUP(fixed_dedup_full_bams)
	
	// sort fixed deduplicated bam in control
	fixed_negative_control_dedup_bam_bai = SORT_TRIM_NEG_CONTROL(fixed_negative_control_dedup_bam)

	// raw negative control reads per chr.
	raw_negative_control_sorted_bams_bais = \
	SPLIT_NEG_CONTROL_RAW(negative_control_data_raw_tuple.combine(chromosomes), params.raw_negative_ctrl_template)

	// UMI-deduplicated negative control reads per chr.
	dedup_negative_control_sorted_bams_bais = \
	SPLIT_NEG_CONTROL_DEDUP(fixed_negative_control_dedup_bam_bai.combine(chromosomes), params.dedup_negative_ctrl_template)

	// merge UMI-deduplicated negative control by all chromosomes. (needed because we discard chromosomes that are not normal chromosomes. e.g. mouse: 1~19, X, Y)
	fixed_negative_control_dedup_full_bam_bai = MERGE_CONTROL(dedup_negative_control_sorted_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom))

	//////////////////////////// get_polyA ////////////////////////////

	// sort_phase1, sort_phase2 can both do sorting and indexing (for only 1 file) and (for all files) because of tuple
	// if for all files: (sorted_bam1, sorted_bam1, specific_sample6), (sorted_bam2, sorted_bam2, specific_sample6), (sorted_bam5, sorted_bam5, specific_sample6), (sorted_bam6, sorted_bam6, specific_sample6)
	// if for only 1 file: (sorted_full_bam, sorted_full_bai, specific_sample5), (sorted_full_bam', sorted_full_bai', specific_sample6)

	// get all full polyA reads and full non polyA reads without a negative control.
	// This is for the motif scoring. You want to compare score between negative control (UMI-tool deduplicated) and those sample data which underwent the pipeline.
	(polyA_sample_only, non_polyA_sample_only) = GET_POLYA_SAMPLE_ONLY(fixed_dedup_sorted_full_bams_bais, params.out_polyA, params.out_non_polyA, params.getPolyA_advanced_script)
	// get all full polyA reads and full non polyA reads of negative control (UMI-tools deduplicated).
	(polyA_control, non_polyA_control) = GET_POLYA_NEG_CONTROL(fixed_negative_control_dedup_full_bam_bai, params.out_polyA_control, params.out_non_polyA_control, params.getPolyA_advanced_script)
	
	// sort full polyA reads for all samples
	polyA_sorted_full_bams_bais_sample_only = SORT_POLYA_SAMPLE_ONLY(polyA_sample_only)
	// sort nonpolyA reads for all samples
	nonpolyA_sorted_full_bams_bais = SORT_NONPOLYA(non_polyA_sample_only)	

	// sort full polyA reads for negative control
	polyA_sorted_full_bams_bais_control = SORT_POLYA_NEG_CONTROL(polyA_control)	
	// sort nonpolyA reads for negative control
	nonpolyA_sorted_full_bams_bais_control = SORT_NONPOLYA_CONTROL(non_polyA_control)

	// polyA_clustering (sample)
	polyA_unique_cs_beds_sample = GET_POLYA_UNIQUE_CLEAVAGE_SITES_BED_SAMPLE(polyA_sorted_full_bams_bais_sample_only, params.polyA_unique_cs_bed_sample, params.get_polyA_unique_cs_script)	
	polyA_clustered_beds_sample_tuple = PERFORM_CLUSTERING_SAMPLE_TUPLE(polyA_unique_cs_beds_sample, params.sample_cluster_out, params.cluster_pas_script)
	
	// non_polyA clustering (sample)
	non_polyA_unique_cs_beds_sample = GET_NON_POLYA_UNIQUE_CLEAVAGE_SITES_BED_SAMPLE(nonpolyA_sorted_full_bams_bais, params.non_polyA_unique_cs_bed_sample, params.get_polyA_unique_cs_script)
	non_polyA_clustered_beds_sample_tuple = PERFORM_CLUSTERING_SAMPLE_TUPLE_NONPOLYA(non_polyA_unique_cs_beds_sample, params.sample_nonpolyA_cluster_out, params.cluster_pas_script)

	// polyA_clustering (control)
	polyA_unique_cs_beds_control = GET_POLYA_UNIQUE_CLEAVAGE_SITES_BED_CONTROL(polyA_sorted_full_bams_bais_control, params.polyA_unique_cs_bed_control, params.get_polyA_unique_cs_script)
	polyA_clustered_bed_control_tuple = PERFORM_CLUSTERING_CONTROL_TUPLE(polyA_unique_cs_beds_control, params.control_cluster_out, params.cluster_pas_script)
	
	// dedup clustering (control)
	dedup_unique_cs_beds_control = GET_DEDUP_UNIQUE_CLEAVAGE_SITES_BED_CONTROL(fixed_negative_control_dedup_full_bam_bai, params.dedup_unique_cs_bed_control, params.get_polyA_unique_cs_script)
	dedup_clustered_bed_control_tuple = PERFORM_CLUSTERING_DEDUP_CONTROL_TUPLE(dedup_unique_cs_beds_control, params.control_dedup_cluster_out, params.cluster_pas_script)

	// non_polyA clustering (control)
	non_polyA_unique_cs_beds_control = GET_NON_POLYA_UNIQUE_CLEAVAGE_SITES_BED_CONTROL(nonpolyA_sorted_full_bams_bais_control, params.non_polyA_unique_cs_bed_control, params.get_polyA_unique_cs_script)
	non_polyA_clustered_beds_control_tuple = PERFORM_CLUSTERING_CONTROL_TUPLE_NONPOLYA(non_polyA_unique_cs_beds_control, params.control_nonpolyA_cluster_out, params.cluster_pas_script)
	
	// classification of bed (sample)
	polyA_intergenic_clustered_beds_tuple = GENES_NO_INTERSECT_TUPLE(polyA_clustered_beds_sample_tuple, genes_bed, params.intergenic_template)
	all_genic_clustered_beds_tuple = GENES_INTERSECT_TUPLE(polyA_clustered_beds_sample_tuple, genes_bed, params.all_genic_template)

	polyA_intronic_clustered_beds_tuple = EXONS_NO_INTERSECT_TUPLE(all_genic_clustered_beds_tuple, exons_bed, params.intronic_template)
	all_exonic_clustered_beds_tuple = EXONS_INTERSECT_TUPLE(all_genic_clustered_beds_tuple, exons_bed, params.all_exonic_template)

	polyA_exonic_clustered_beds_tuple = TE_NO_INTERSECT_TUPLE(all_exonic_clustered_beds_tuple, terminal_exons_bed, params.exonic_template)
	polyA_TE_clustered_beds_tuple = TE_INTERSECT_TUPLE(all_exonic_clustered_beds_tuple, terminal_exons_bed, params.te_template)
		
	// classification of bam (sample)	
	polyA_sorted_bams_bais_sample = SPLIT_ALL_POLYA(polyA_sorted_full_bams_bais_sample_only.combine(chromosomes), params.out_polyA_partial)
	polyA_sorted_bams_sample_bais = CHANGE_TUPLE_ORDER_POLYA_SAMPLE(polyA_sorted_bams_bais_sample)

	(te_bams, non_te_bams) = CLASSIFY_TE_SAMPLE(polyA_sorted_bams_sample_bais.combine(polyA_TE_clustered_beds_tuple, by:1), params.te_template, params.non_te_template, "TE", params.classify_script)
	te_sorted_bams_sample_bais = SORT_CLASSIFIED_TE(te_bams)
	non_te_sorted_bams_sample_bais = SORT_CLASSIFIED_NON_TE(non_te_bams)

	(e_bams, non_te_e_bams) = CLASSIFY_E_SAMPLE(non_te_sorted_bams_sample_bais.combine(polyA_exonic_clustered_beds_tuple, by:1), params.exonic_template, params.non_te_e_template, "exonic", params.classify_script)
	e_sorted_bams_sample_bais = SORT_CLASSIFIED_EXONIC(e_bams)
	non_te_e_sorted_bams_sample_bais = SORT_CLASSIFIED_NON_TE_NON_E(non_te_e_bams)

	(intronic_bams, non_te_e_i_bams) = CLASSIFY_INTRONIC_SAMPLE(non_te_e_sorted_bams_sample_bais.combine(polyA_intronic_clustered_beds_tuple, by:1), params.intronic_template, params.non_te_e_i_template, "intronic", params.classify_script)
	intronic_sorted_bams_sample_bais = SORT_CLASSIFIED_INTRONIC(intronic_bams)
	non_te_e_i_sorted_bams_sample_bais = SORT_CLASSIFIED_NON_TE_NON_E_NON_I(non_te_e_i_bams) 

	(intergenic_bams, non_te_e_i_intergenic_bams) = CLASSIFY_INTERGENIC_SAMPLE(non_te_e_i_sorted_bams_sample_bais.combine(polyA_intergenic_clustered_beds_tuple, by:1), params.intergenic_template, params.non_te_e_i_intergenic_template, "intergenic", params.classify_script)
	intergenic_sorted_bams_sample_bais = SORT_CLASSIFIED_INTERGENIC(intergenic_bams)
	non_te_e_i_intergenic_bams_sample_bais = SORT_CLASSIFIED_NON_TE_NON_E_NON_I_NON_IG(non_te_e_i_intergenic_bams)

	polyA_terminal_sorted_full_bams_bais = MERGE_CLASSIFIED_TE(te_sorted_bams_sample_bais.groupTuple(by:1, sort: true, size: num_chrom), params.te_template)
	polyA_exonic_sorted_full_bams_bais = MERGE_CLASSIFIED_E(e_sorted_bams_sample_bais.groupTuple(by:1, sort: true, size: num_chrom), params.exonic_template)
	intronic_sorted_full_bams_bais = MERGE_CLASSIFIED_INTRONIC(intronic_sorted_bams_sample_bais.groupTuple(by:1, sort: true, size: num_chrom), params.intronic_template)
	intergenic_sorted_full_bams_bais = MERGE_CLASSIFIED_INTERGENIC(intergenic_sorted_bams_sample_bais.groupTuple(by:1, sort: true, size: num_chrom), params.intergenic_template)
	
	(polyA_TE_below_beds_tuple, polyA_TE_above_beds_tuple) = ISOLATE_PA_T_SAMPLE(polyA_TE_clustered_beds_tuple, terminal_exons_bed, params.isolate_pA_T_bed_script)
	
	(annotated_bams, non_annotated_bams) = CLASSIFY_BELOW_SAMPLE(te_sorted_bams_sample_bais.combine(polyA_TE_below_beds_tuple, by:1), params.below_template, params.non_annotated_template, "annotated", params.classify_below_above_script)
	annotated_sorted_bams_sample_bais = SORT_CLASSIFIED_ANNOTATED(annotated_bams)
	non_annotated_sorted_bams_sample_bais = SORT_CLASSIFIED_NON_ANNOTATED(non_annotated_bams)

	(unannotated_bams, non_annotated_unannotated_bams) = CLASSIFY_ABOVE_SAMPLE(non_annotated_sorted_bams_sample_bais.combine(polyA_TE_above_beds_tuple, by:1), params.above_template, params.non_annotated_unannotated_template, "unannotated", params.classify_below_above_script)
	unannotated_sorted_bams_sample_bais = SORT_CLASSIFIED_UNANNOTATED(unannotated_bams)
	non_annotated_unannotated_sorted_bams_sample_bais = SORT_CLASSIFIED_NON_ANNOTATED_UNANNOTATED(non_annotated_unannotated_bams)

	isolated_below_full_bams_bais = MERGE_CLASSIFIED_ANNOTATED(annotated_sorted_bams_sample_bais.groupTuple(by:1, sort: true, size: num_chrom), params.below_template)
	isolated_above_full_bams_bais = MERGE_CLASSIFIED_UNANNOTATED(unannotated_sorted_bams_sample_bais.groupTuple(by:1, sort: true, size: num_chrom), params.above_template)

	// classification of bed (control)
	polyA_intergenic_clustered_beds_control_tuple = GENES_NO_INTERSECT_CONTROL_TUPLE(polyA_clustered_bed_control_tuple, genes_bed, params.intergenic_control_template)
	all_genic_clustered_beds_control_tuple = GENES_INTERSECT_CONTROL_TUPLE(polyA_clustered_bed_control_tuple, genes_bed, params.all_genic_control_template)

	polyA_intronic_clustered_beds_control_tuple = EXONS_NO_INTERSECT_CONTROL_TUPLE(all_genic_clustered_beds_control_tuple, exons_bed, params.intronic_control_template)
	all_exonic_clustered_beds_control_tuple = EXONS_INTERSECT_CONTROL_TUPLE(all_genic_clustered_beds_control_tuple, exons_bed, params.all_exonic_control_template)

	polyA_exonic_clustered_beds_control_tuple = TE_NO_INTERSECT_CONTROL_TUPLE(all_exonic_clustered_beds_control_tuple, terminal_exons_bed, params.exonic_control_template)
	polyA_TE_clustered_beds_control_tuple = TE_INTERSECT_CONTROL_TUPLE(all_exonic_clustered_beds_control_tuple, terminal_exons_bed, params.te_control_template)

	// classification of bam (control)
	polyA_sorted_bams_bais_control = SPLIT_ALL_POLYA_CONTROL(polyA_sorted_full_bams_bais_control.combine(chromosomes), params.out_polyA_partial_control)
	polyA_sorted_bams_sample_bais_control = CHANGE_TUPLE_ORDER_POLYA_CONTROL(polyA_sorted_bams_bais_control)

	(te_bams_control, non_te_bams_control) = CLASSIFY_TE_CONTROL(polyA_sorted_bams_sample_bais_control.combine(polyA_TE_clustered_beds_control_tuple, by:1), params.te_control_template, params.non_te_control_template, "TE", params.classify_script)
	te_sorted_bams_sample_bais_control = SORT_CLASSIFIED_TE_CONTROL(te_bams_control)
	non_te_sorted_bams_sample_bais_control = SORT_CLASSIFIED_NON_TE_CONTROL(non_te_bams_control)

	(e_bams_control, non_te_e_bams_control) = CLASSIFY_E_CONTROL(non_te_sorted_bams_sample_bais_control.combine(polyA_exonic_clustered_beds_control_tuple, by:1), params.exonic_control_template, params.non_te_e_control_template, "exonic", params.classify_script)
	e_sorted_bams_sample_bais_control = SORT_CLASSIFIED_EXONIC_CONTROL(e_bams_control)
	non_te_e_sorted_bams_sample_bais_control = SORT_CLASSIFIED_NON_TE_NON_E_CONTROL(non_te_e_bams_control)

	(intronic_bams_control, non_te_e_i_bams_control) = CLASSIFY_INTRONIC_CONTROL(non_te_e_sorted_bams_sample_bais_control.combine(polyA_intronic_clustered_beds_control_tuple, by:1), params.intronic_control_template, params.non_te_e_i_control_template, "intronic", params.classify_script)
	intronic_sorted_bams_sample_bais_control = SORT_CLASSIFIED_INTRONIC_CONTROL(intronic_bams_control)
	non_te_e_i_sorted_bams_sample_bais_control = SORT_CLASSIFIED_NON_TE_NON_E_NON_I_CONTROL(non_te_e_i_bams_control) 

	(intergenic_bams_control, non_te_e_i_intergenic_bams_control) = CLASSIFY_INTERGENIC_CONTROL(non_te_e_i_sorted_bams_sample_bais_control.combine(polyA_intergenic_clustered_beds_control_tuple, by:1), params.intergenic_control_template, params.non_te_e_i_intergenic_control_template, "intergenic", params.classify_script)
	intergenic_sorted_bams_sample_bais_control = SORT_CLASSIFIED_INTERGENIC_CONTROL(intergenic_bams_control)
	non_te_e_i_intergenic_bams_sample_bais_control = SORT_CLASSIFIED_NON_TE_NON_E_NON_I_NON_IG_CONTROL(non_te_e_i_intergenic_bams_control)

	polyA_terminal_sorted_full_control_bams_bais = MERGE_CLASSIFIED_TE_CONTROL(te_sorted_bams_sample_bais_control.groupTuple(by:1, sort: true, size: num_chrom), params.te_control_template)
	polyA_exonic_sorted_full_bams_bais_control = MERGE_CLASSIFIED_E_CONTROL(e_sorted_bams_sample_bais_control.groupTuple(by:1, sort: true, size: num_chrom), params.exonic_control_template)
	intronic_sorted_full_bams_bais_control = MERGE_CLASSIFIED_INTRONIC_CONTROL(intronic_sorted_bams_sample_bais_control.groupTuple(by:1, sort: true, size: num_chrom), params.intronic_control_template)
	intergenic_sorted_full_bams_bais_control = MERGE_CLASSIFIED_INTERGENIC_CONTROL(intergenic_sorted_bams_sample_bais_control.groupTuple(by:1, sort: true, size: num_chrom), params.intergenic_control_template)

	(polyA_TE_below_beds_control_tuple, polyA_TE_above_beds_control_tuple) = ISOLATE_PA_T_CONTROL(polyA_TE_clustered_beds_control_tuple, terminal_exons_bed, params.isolate_pA_T_bed_script)
	
	(annotated_bams_control, non_annotated_bams_control) = CLASSIFY_BELOW_CONTROL(te_sorted_bams_sample_bais_control.combine(polyA_TE_below_beds_control_tuple, by:1), params.below_template, params.non_annotated_template, "annotated", params.classify_below_above_script)
	annotated_sorted_bams_sample_bais_control = SORT_CLASSIFIED_ANNOTATED_CONTROL(annotated_bams_control)
	non_annotated_sorted_bams_sample_bais_control = SORT_CLASSIFIED_NON_ANNOTATED_CONTROL(non_annotated_bams_control)

	(unannotated_bams_control, non_annotated_unannotated_bams_control) = CLASSIFY_ABOVE_CONTROL(non_annotated_sorted_bams_sample_bais_control.combine(polyA_TE_above_beds_control_tuple, by:1), params.above_template, params.non_annotated_unannotated_template, "unannotated", params.classify_below_above_script)
	unannotated_sorted_bams_sample_bais_control = SORT_CLASSIFIED_UNANNOTATED_CONTROL(unannotated_bams_control)
	non_annotated_unannotated_sorted_bams_sample_bais_control = SORT_CLASSIFIED_NON_ANNOTATED_UNANNOTATED_CONTROL(non_annotated_unannotated_bams_control)

	isolated_below_full_bams_bais_control = MERGE_CLASSIFIED_ANNOTATED_CONTROL(annotated_sorted_bams_sample_bais_control.groupTuple(by:1, sort: true, size: num_chrom), params.below_template)
	isolated_above_full_bams_bais_control = MERGE_CLASSIFIED_UNANNOTATED_CONTROL(unannotated_sorted_bams_sample_bais_control.groupTuple(by:1, sort: true, size: num_chrom), params.above_template)
	
	//////////////////////////// analysis ////////////////////////////

	if (params.analysis == "yes"){
		//////////////////////////// draw span distribution to compare raw vs dedup vs split ////////////////////////////
	
		// raw span distribution of samples and control.
		raw_spans_samples = GET_SPAN_RAW_SAMPLE(possorted_bams_bais, "raw", "sample", params.span_script)
		raw_spans_control = GET_SPAN_RAW_CONTROL(raw_negative_control_sorted_bams_bais, "raw", "control", params.span_script)
		raw_spans = raw_spans_samples.mix(raw_spans_control)

		// span distribution after deduplication of samples and control.
		dedup_spans_samples = GET_SPAN_DEDUP_SAMPLE(sorted_dedup_bams_bais, "dedup", "sample", params.span_script)
		dedup_spans_control = GET_SPAN_DEDUP_CONTROL(dedup_negative_control_sorted_bams_bais, "dedup", "control", params.span_script)
		dedup_spans = dedup_spans_samples.mix(dedup_spans_control)
		
		// span after "split" step in samples.
		// This is only available for samples (not control).
		split_spans = GET_SPAN_SPLIT(splitted_sorted_bams_bais, "split", params.span_script)

		raw_span_merged = CONCAT_SPAN_RAW(raw_spans.groupTuple(by: 1, sort: true, size: num_chrom), params.raw_span)
		dedup_span_merged = CONCAT_SPAN_DEDUP(dedup_spans.groupTuple(by: 1, sort: true, size: num_chrom), params.dedup_span)
		split_span_merged = CONCAT_SPAN_SPLIT(split_spans.groupTuple(by: 1, sort: true, size: num_chrom), params.split_span)

		raw_span_hist = GET_SPAN_HIST_RAW(raw_span_merged, params.span_histogram_script, params.raw_span_hist)
		dedup_span_hist = GET_SPAN_HIST_DEDUP(dedup_span_merged, params.span_histogram_script, params.dedup_span_hist)
		split_span_hist = GET_SPAN_HIST_SPLIT(split_span_merged, params.span_histogram_script, params.split_span_hist)
		
		further_filtered_terminal_exons = FILTER_FURTHER_TERMINAL_EXONS(terminal_exons_bed, params.filter_terminal_exons_script)

		//////////////////////////// get_A_freq_in_nonsoftclipped_reads ////////////////////////////
		// get reads that do not have soft-clipped regions.
		no_softclip_full_bams = GET_NON_SOFTCLIP(fixed_dedup_sorted_full_bams_bais, params.get_no_softclipped_script)
		no_softclip_sorted_full_bams_bais = SORT_NO_SOFTCLIP(no_softclip_full_bams)

		// get frequency of "A" nucleotide appearing by chance (per each sample). it should be <= 0.3
		// if this "A" frequency <= 0.3, softclipped length >= 5 is useable as polyA reads. (probability of happening by chance is reasonably low)
		avg_A_freq_csvs = GET_A_FREQ_IN_NONSOFTCLIP(no_softclip_sorted_full_bams_bais, params.get_avg_a_freq_script)

		////////////////////////////  num_reads_per_class ////////////////////////////

		// get counts csvs per "raw" class and per sample name.
		raw_counts_csvs = GET_COUNTS_RAW(full_data_tuple, params.counts_raw, params.raw_type, params.get_count_script)
		raw_counts_csvs_control = GET_COUNTS_RAW_CONTROL(negative_control_data_raw_tuple, params.counts_raw, params.raw_type, params.get_count_script)
		// merges counts csvs per "raw" class across all samples.
		merged_raw_counts = CONCAT_COUNTS_RAW(raw_counts_csvs.collect(), params.raw_type)
		
		// get counts csvs per "dedup" class and per sample name.
		dedup_counts_csvs = GET_COUNTS_DEDUP(fixed_dedup_sorted_full_bams_bais, params.counts_dedup, params.dedup_type, params.get_count_script)
		dedup_counts_csvs_control = GET_COUNTS_DEDUP_CONTROL(fixed_negative_control_dedup_full_bam_bai, params.counts_dedup, params.dedup_type, params.get_count_script)
		// merges counts csvs per "dedup" across all samples.
		merged_dedup_counts = CONCAT_COUNTS_DEDUP(dedup_counts_csvs.collect(), params.dedup_type)

		// get counts csvs of # of softclipped reads per "dedup" class and per sample name.
		dedup_softclipped_counts_csvs = GET_SOFTCLIPPED_COUNTS_DEDUP(fixed_dedup_sorted_full_bams_bais,\
		params.counts_dedup_softclipped, params.dedup_softclipped_type, params.get_softclip_count_script)

		dedup_softclipped_counts_csvs_control = GET_COUNTS_SOFTCLIPPED_CONTROL(fixed_negative_control_dedup_full_bam_bai,\
		params.counts_dedup_softclipped, params.dedup_softclipped_type, params.get_softclip_count_script)
		// merges counts csvs of # of softclipped reads per "dedup" across all samples.
		merged_softclipped_dedup_counts = CONCAT_SOFTCLIPPED_COUNTS_DEDUP(dedup_softclipped_counts_csvs.collect(), params.dedup_softclipped_type)

		// get counts csvs per non-polyA class and per sample name.
		nonpolyA_counts_csvs = GET_COUNTS_NONPOLYA(nonpolyA_sorted_full_bams_bais, params.counts_nonpolyA, params.nonpolyA_type, params.get_count_script)
		nonpolyA_counts_csvs_control = GET_COUNTS_NONPOLYA_CONTROL(nonpolyA_sorted_full_bams_bais_control, params.counts_nonpolyA, params.nonpolyA_type, params.get_count_script)
		// merges counts csvs per non-polyA across all samples.
		merged_nonpolyA_counts = CONCAT_COUNTS_NONPOLYA(nonpolyA_counts_csvs.collect(), params.nonpolyA_type)

		// get counts csvs per polyA class and per sample name.
		polyA_counts_csvs = GET_COUNTS_POLYA(polyA_sorted_full_bams_bais_sample_only, params.counts_polyA, params.polyA_type, params.get_count_script)
		polyA_counts_csvs_control = GET_COUNTS_POLYA_CONTROL(polyA_sorted_full_bams_bais_control, params.counts_polyA, params.polyA_type, params.get_count_script)
		// merges counts csvs per polyA class across all samples.
		merged_polyA_counts = CONCAT_COUNTS_POLYA(polyA_counts_csvs.collect(), params.polyA_type)
		
		// get counts csvs per polyA + terminal class and per sample name.
		polyAT_counts_csvs = GET_COUNTS_POLYA_TERMINAL(polyA_terminal_sorted_full_bams_bais, params.counts_polyAT, params.polyAT_type, params.get_count_script)
		polyAT_counts_csvs_control = GET_COUNTS_POLYA_T_CONTROL(polyA_terminal_sorted_full_control_bams_bais, params.counts_polyAT, params.polyAT_type, params.get_count_script)
		// merges counts csvs per polyA + terminal class across all samples.
		merged_polyAT_counts = CONCAT_COUNTS_POLYA_TERMINAL(polyAT_counts_csvs.collect(), params.polyAT_type)
		
		// get counts csvs per annotated class and per sample name.
		annotated_counts_csvs = GET_COUNTS_ANNOTATED(isolated_below_full_bams_bais, params.counts_annotated, params.annotated_type, params.get_count_script)
		annotated_counts_csvs_control = GET_COUNTS_ANNOTATED_CONTROL(isolated_below_full_bams_bais_control, params.counts_annotated, params.annotated_type, params.get_count_script)
		// merges counts csvs per annotated class across all samples.
		merged_annotated_counts = CONCAT_COUNTS_ANNOTATED(annotated_counts_csvs.collect(), params.annotated_type)

		// get counts csvs per unannotated class and per sample name.
		unannotated_counts_csvs = GET_COUNTS_UNANNOTATED(isolated_above_full_bams_bais, params.counts_unannotated, params.unannotated_type, params.get_count_script)
		unannotated_counts_csvs_control = GET_COUNTS_UNANNOTATED_CONTROL(isolated_above_full_bams_bais_control, params.counts_unannotated, params.unannotated_type, params.get_count_script)
		// merges counts csvs per unannotated class across all samples.
		merged_unannotated_counts = CONCAT_COUNTS_UNANNOTATED(unannotated_counts_csvs.collect(), params.unannotated_type)

		// get counts csvs per intronic class and per sample name.
		intronic_counts_csvs = GET_COUNTS_INTRONIC(intronic_sorted_full_bams_bais, params.counts_intronic, params.intronic_type, params.get_count_script)
		intronic_counts_csvs_control = GET_COUNTS_INTRONIC_CONTROL(intronic_sorted_full_bams_bais_control, params.counts_intronic, params.intronic_type, params.get_count_script)
		// merges counts csvs per intronic class across all samples.
		merged_intronic_counts = CONCAT_COUNTS_INTRONIC(intronic_counts_csvs.collect(), params.intronic_type)
		
		// get counts csvs per intergenic class and per sample name.
		intergenic_counts_csvs = GET_COUNTS_INTERGENIC(intergenic_sorted_full_bams_bais, params.counts_intergenic, params.intergenic_type, params.get_count_script)
		intergenic_counts_csvs_control = GET_COUNTS_INTERGENIC_CONTROL(intergenic_sorted_full_bams_bais_control, params.counts_intergenic, params.intergenic_type, params.get_count_script)
		// merges counts csvs per intergenic class across all samples.
		merged_intergenic_counts = CONCAT_COUNTS_INTERGENIC(intergenic_counts_csvs.collect(), params.intergenic_type)
		
		// get counts csvs per exonic class and per sample name.
		exonic_counts_csvs = GET_COUNTS_EXONIC(polyA_exonic_sorted_full_bams_bais, params.counts_exonic, params.exonic_type, params.get_count_script)
		exonic_counts_csvs_control = GET_COUNTS_EXONIC_CONTROL(polyA_exonic_sorted_full_bams_bais_control, params.counts_exonic, params.exonic_type, params.get_count_script)
		// merges counts csvs per exonic class across all samples.
		merged_exonic_counts = CONCAT_COUNTS_EXONIC(exonic_counts_csvs.collect(), params.exonic_type)

		// merges counts csvs across all classes and across all samples
		// total_counts_csv = MERGE_ALL_COUNTS_SAMPLE(merged_raw_counts, merged_dedup_counts, merged_softclipped_dedup_counts, merged_nonpolyA_counts,\
		// merged_polyA_counts, merged_polyAT_counts, merged_intronic_counts, merged_intergenic_counts, merged_exonic_counts, "total_merged_sample")

		total_counts_csv = MERGE_ALL_COUNTS_SAMPLE(merged_raw_counts, merged_dedup_counts, merged_softclipped_dedup_counts, merged_nonpolyA_counts,\
		merged_polyA_counts, merged_polyAT_counts, merged_annotated_counts, merged_unannotated_counts, merged_intronic_counts, merged_intergenic_counts, merged_exonic_counts, "total_merged_sample")

		total_counts_csv_control = MERGE_ALL_COUNTS_CONTROL(raw_counts_csvs_control, dedup_counts_csvs_control, dedup_softclipped_counts_csvs_control, nonpolyA_counts_csvs_control, polyA_counts_csvs_control,\
		polyAT_counts_csvs_control, annotated_counts_csvs_control, unannotated_counts_csvs_control, intronic_counts_csvs_control, intergenic_counts_csvs_control, exonic_counts_csvs_control, "total_merged_control")
		
		// split counts csv by sample
		split_by_samples_counts_csvs_collected = SPLIT_PER_SAMPLE(total_counts_csv, params.get_count_split_script)

		split_by_samples_counts_csvs = ONE_BY_ONE(split_by_samples_counts_csvs_collected)
		
		// make tuple of file path and sample name
		counts_and_samples = split_by_samples_counts_csvs.map { file -> [file, file.baseName] }

		// make counts bar chart
		// use combine so that you use negative control as many times as the number of samples.
		(count_barcharts_png, count_barcharts_svg) = GET_BAR_CHART(counts_and_samples.combine(total_counts_csv_control), params.get_barchart_script)

		//////////////////////////// distance ////////////////////////////

		// for distance distribution and avg distance (type1 vs type2)
		polyA_TE_further_filtered_clustered_beds_tuple = TE_FURTHER_INTERSECT_TUPLE(polyA_TE_clustered_beds_tuple, further_filtered_terminal_exons, params.filtered_te_template)
		
		// get distance distribution between end of terminal exon and representative cleavage site of a pA cluster in the sample.
		distance_outputs_alter_tuple = GET_DISTANCE_ALTER_SAMPLE(polyA_TE_further_filtered_clustered_beds_tuple,\
		further_filtered_terminal_exons, params.get_distance_final_script)

		// get distance distribution between end of terminal exon and representative cleavage site of a pA cluster in the control.
		distance_outputs_alter_control_tuple = GET_DISTANCE_ALTER_CONTROL(dedup_clustered_bed_control_tuple,\
		further_filtered_terminal_exons, params.get_distance_final_script)
		
		distance_outputs_alter_control = UNCOUPLE_CSV_WITH_SAMPLE(distance_outputs_alter_control_tuple)
		
		//distance histogram (before and after fixing cleavage site)
		(distance_histograms_alter_png, distance_histograms_alter_svg) = GET_D_HISTOGRAM_ALTER(distance_outputs_alter_tuple.combine(distance_outputs_alter_control), params.get_histogram_final_script)

		//////////////////////////// ATGC_plot ////////////////////////////

		// plot frequency plot in each class
		all_polyA_atgc_plts = PLOT_ATGC_ALL_POLYA(polyA_clustered_beds_sample_tuple, params.all_polyA_ATGC_out, params.plot_frequencies_script)

		above_atgc_plts = PLOT_ATGC_ABOVE(polyA_TE_above_beds_tuple, params.unannotated_ATGC_out, params.plot_frequencies_script)

		below_atgc_plts = PLOT_ATGC_BELOW(polyA_TE_below_beds_tuple, params.annotated_ATGC_out, params.plot_frequencies_script)

		nonpolyA_atgc_plts = PLOT_ATGC_NONPOLYA(non_polyA_clustered_beds_sample_tuple, params.non_polyA_ATGC_out, params.plot_frequencies_script)

		intergenic_atgc_plts = PLOT_ATGC_INTERGENIC(polyA_intergenic_clustered_beds_tuple, params.intergenic_ATGC_out, params.plot_frequencies_script)

		intronic_atgc_plts = PLOT_ATGC_INTRONIC(polyA_intronic_clustered_beds_tuple, params.intronic_ATGC_out, params.plot_frequencies_script)

		exonic_atgc_plt = PLOT_ATGC_EXONIC(polyA_exonic_clustered_beds_tuple, params.polyA_exonic_ATGC_out, params.plot_frequencies_script)

		control_atgc_plt = PLOT_ATGC_CONTROL(dedup_clustered_bed_control_tuple, params.control_ATGC_out, params.plot_frequencies_script)

		//////////////////////////// motif_freq_plot ////////////////////////////

		// draw 12 motives frequency plots based on gtf file
		(gtf_motif_graphs, gtf_motif_graphs_overlaid_png, gtf_motif_graphs_overlaid_svg, gtf_motif_order, gtf_peaks) = GET_MOTIF_FREQ_GTF_PLOT(further_filtered_terminal_exons, params.motif_search_out_gtf, params.gtf_search_motif_script)
		
		// draw 12 motives frequency plots from our data in each class in each sample
		(annotated_motif_graphs, annotated_motif_graphs_overlaid_png, annotated_motif_graphs_overlaid_svg, annotated_motif_orders, annotated_motif_scores) = GET_MOTIF_FREQ_PLOT_BELOW(polyA_TE_below_beds_tuple,\
		params.motif_search_out_template_annotated, params.below_annotated, gtf_peaks, params.search_motif_script)

		(unannotated_motif_graphs, unannotated_motif_graphs_overlaid_png, unannotated_motif_graphs_overlaid_svg, unannotated_motif_orders, unannotated_motif_scores) = GET_MOTIF_FREQ_PLOT_ABOVE(polyA_TE_above_beds_tuple,\
		params.motif_search_out_template_unannotated, params.above_annotated, gtf_peaks, params.search_motif_script)

		(intergenic_motif_graphs, intergenic_motif_graphs_overlaid_png, intergenic_motif_graphs_overlaid_svg, intergenic_motif_orders, intergenic_motif_scores) = GET_MOTIF_FREQ_PLOT_INTERGENIC(polyA_intergenic_clustered_beds_tuple,\
		params.motif_search_out_template_intergenic, params.intergenic_annotated, gtf_peaks, params.search_motif_script)	

		(intronic_motif_graphs, intronic_motif_graphs_overlaid_png, intronic_motif_graphs_overlaid_svg, intronic_motif_orders, intronic_motif_scores) = GET_MOTIF_FREQ_PLOT_INTRONIC(polyA_intronic_clustered_beds_tuple,\
		params.motif_search_out_template_intronic, params.intronic_annotated, gtf_peaks, params.search_motif_script)	

		(exonic_motif_graphs, exonic_motif_graphs_overlaid_png, exonic_motif_graphs_overlaid_svg, exonic_motif_orders, exonic_motif_scores) = GET_MOTIF_FREQ_PLOT_EXONIC(polyA_exonic_clustered_beds_tuple,\
		params.motif_search_out_template_exonic, params.polyA_exonic_annotated, gtf_peaks, params.search_motif_script)

		// save motif scores/orders (min:0, max:12) in each class and in each sample
		(merged_annotated_motif_orders, merged_annotated_motif_scores) = CONCAT_MOTIF_ORDERS_ANNOTATED(annotated_motif_orders.collect(), annotated_motif_scores.collect(), params.motif_search_out_template_annotated)

		(merged_unannotated_motif_orders, merged_unannotated_motif_scores)  = CONCAT_MOTIF_ORDERS_UNANNOTATED(unannotated_motif_orders.collect(), unannotated_motif_scores.collect(), params.motif_search_out_template_unannotated)
		
		(merged_intergenic_motif_orders, merged_intergenic_motif_scores) = CONCAT_MOTIF_ORDERS_INTERGENIC(intergenic_motif_orders.collect(), intergenic_motif_scores.collect(), params.motif_search_out_template_intergenic)
		
		(merged_intronic_motif_orders, merged_intronic_motif_scores) = CONCAT_MOTIF_ORDERS_INTRONIC(intronic_motif_orders.collect(), intronic_motif_scores.collect(), params.motif_search_out_template_intronic)
		
		(merged_exonic_motif_orders, merged_exonic_motif_scores) = CONCAT_MOTIF_ORDERS_EXONIC(exonic_motif_orders.collect(), exonic_motif_scores.collect(), params.motif_search_out_template_exonic)

		// merge motif scores/orders across all samples
		total_motif_orders_csv = MERGE_ALL_MOTIF_ORDERS(merged_annotated_motif_orders, merged_unannotated_motif_orders, merged_intergenic_motif_orders, merged_intronic_motif_orders, merged_exonic_motif_orders, params.total_merged_motif_orders)
		total_motif_scores_csv = MERGE_ALL_MOTIF_SCORES(merged_annotated_motif_scores, merged_unannotated_motif_scores, merged_intergenic_motif_scores, merged_intronic_motif_scores, merged_exonic_motif_scores, params.total_merged_motif_scores)

		// split motif orders by each sample. Each sample has motif orders for annotated, unannotated, intronic, intergenic, exonic.
		split_by_samples_motif_orders_csv_collected = SPLIT_MOTIF_ORDER_PER_SAMPLE(total_motif_orders_csv, params.split_motif_order_script)
		split_by_samples_motif_orders_csv = MOTIF_ORDERS_ONE_BY_ONE(split_by_samples_motif_orders_csv_collected)
		
		// make tuple of file path and sample name
		motif_orders_and_samples = split_by_samples_motif_orders_csv.map { file -> [file, file.baseName] }

		// combine motif orders of 1 sample + motif order of gtf. Each sample has motif orders for annotated, unannotated, intronic, intergenic, exonic and gtf.
		combined_motif_order = COMBINE_MOTIF_ORDERS(motif_orders_and_samples, gtf_motif_order, params.motif_order_outname)

		// compute motif scores (min:0, max:12) at the pulled polyA reads level per each sample.
		(polyA_motif_graphs, polyA_motif_graphs_overlaid_png, polyA_motif_graphs_overlaid_svg, polyA_scores_csv) = COMPUTE_SCORES_FOR_MOTIF_ALTER(polyA_clustered_beds_sample_tuple.mix(non_polyA_clustered_beds_control_tuple),\
		params.motif_search_out_polyA, gtf_peaks, params.compute_motif_score_advanced_script)

		// combine motif scores (min:0, max:12) at the pulled polyA reads level across all samples.
		combined_motif_scores_polyA_level = COMBINE_SCORES_FOR_MOTIF(polyA_scores_csv.collect(), params.combined_motif_scores)
		(polyA_motif_scores_png, polyA_motif_scores_svg) = GET_BAR_CHART_SCORE_MOTIF_POLYA(combined_motif_scores_polyA_level, params.polyA_motif_score_script)

		// plot histogram of length of softclipped region (before and after fixation).
		(size_softclipped_dist_png, size_softclipped_dist_svg, percentage_csv) = GET_SOFTCLIPPED_DISTRIBUTION(fixed_dedup_sorted_full_bams_bais, params.sizeSoftclipped_out, params.get_softclipped_script)

		//////////////////////////// num_polyA (number of polyA sites identified) ////////////////////////////
		
		// get number of polyA sites in each sample (using all polyA reads). + get the number of polyA sites/gene distribution in each sample.
		num_polyA_csvs = GET_NUM_POLYA_SITES_ALLPOLYA(polyA_clustered_beds_sample_tuple,\
		params.allPAS_num_polyA_csv, params.all_polyA_annotated, params.get_numA_polyA_script)

		// get number of polyA sites in each sample and in each class. + get the number of polyA sites/gene distribution in each sample. (for isolated above/below reads that have GX tag rest are dummy figures) 
		above_num_polyA_csvs = GET_NUM_POLYA_SITES_ABOVE(polyA_TE_above_beds_tuple,\
		params.above_num_polyA_csv, params.above_annotated, params.get_numA_polyA_script)

		above_combined_numPolyA_csv = COMBINE_NUM_POLYA_ABOVE(above_num_polyA_csvs.collect(), params.combined_above_num_polyA_csv)

		below_num_polyA_csvs = GET_NUM_POLYA_SITES_BELOW(polyA_TE_below_beds_tuple,\
		params.below_num_polyA_csv, params.below_annotated, params.get_numA_polyA_script)

		below_combined_numPolyA_csv = COMBINE_NUM_POLYA_BELOW(below_num_polyA_csvs.collect(), params.combined_below_num_polyA_csv)

		intronic_num_polyA_csvs = GET_NUM_POLYA_SITES_INTRONIC(polyA_intronic_clustered_beds_tuple,\
		params.intronic_num_polyA_csv, params.intronic_annotated, params.get_numA_polyA_script)

		intronic_combined_numPolyA_csv = COMBINE_NUM_POLYA_INTRONIC(intronic_num_polyA_csvs.collect(), params.combined_intronic_num_polyA_csv)

		intergenic_num_polyA_csvs = GET_NUM_POLYA_SITES_INTERGENIC(polyA_intergenic_clustered_beds_tuple,\
		params.intergenic_num_polyA_csv, params.intergenic_annotated, params.get_numA_polyA_script)

		intergenic_combined_numPolyA_csv = COMBINE_NUM_POLYA_INTERGENIC(intergenic_num_polyA_csvs.collect(), params.combined_intergenic_num_polyA_csv)

		exonic_num_polyA_csvs = GET_NUM_POLYA_SITES_EXONIC(polyA_exonic_clustered_beds_tuple,\
		params.exonic_num_polyA_csv, params.polyA_exonic_annotated, params.get_numA_polyA_script)

		exonic_combined_numPolyA_csv = COMBINE_NUM_POLYA_EXONIC(exonic_num_polyA_csvs.collect(), params.combined_exonic_num_polyA_csv)

		control_num_polyA_csvs = GET_NUM_POLYA_SITES_CONTROL(polyA_clustered_bed_control_tuple,\
		params.control_num_polyA_csv, params.control_annotated, params.get_numA_polyA_script)

		control_combined_numPolyA_csv = COMBINE_NUM_POLYA_CONTROL(control_num_polyA_csvs.collect(), params.combined_control_num_polyA_csv)

		// combine all csvs containing number of pA sites (at all polyA reads) from each sample into a single csv file.
		combined_numPolyA_csv = COMBINE_NUM_POLYA_ALL_PAS(num_polyA_csvs.mix(control_combined_numPolyA_csv).collect(), params.combined_all_PAS_num_polyA_csv)
		// plot a barchart of number of polyA sites for all samples (using all polyA reads).
		(num_polyA_barchart_png, num_polyA_barchart_svg) = GET_BAR_CHART_NUM_POLYA(combined_numPolyA_csv, params.num_polyA_barchart, params.get_barchart_numpolyA_script)
		
		// merge all number of polyA sites from different sample and different class into a single csv file.
		all_combined_numPolyA_csv = MERGE_ALL_NUM_PAS_CSVS(above_combined_numPolyA_csv, below_combined_numPolyA_csv,\
		intronic_combined_numPolyA_csv, intergenic_combined_numPolyA_csv, exonic_combined_numPolyA_csv, control_combined_numPolyA_csv, params.all_combined_num_polyA_csv)
		
		// combine total motif score csv file, total number of polyA sites csv file and total number of genes csv file together and make a barchart.
		motif_score_and_numPAS_csv = FUSE_MOTIF_SCORE_NUM_PAS(total_motif_scores_csv, combined_motif_scores_polyA_level, all_combined_numPolyA_csv, params.numPAS_motifScore_csv, params.fuse_numPAS_motifscore_script)
		
		//////////////////////////// gene_coverage (number of genes covered by each sample) ////////////////////////////

		// compute number of genes expressed by each sample and save it in each csv file. (includes negative control)
		dedup_gene_coverage_csvs = GET_GENE_COVERAGE_DEDUP(fixed_dedup_sorted_full_bams_bais.mix(fixed_negative_control_dedup_full_bam_bai), params.dedup_gene_coverage_csv, terminal_exons_bed, "dedup", params.gene_coverage_script)
		// compute number of genes covered by each sample and save it in each csv file. (includes negative control)
		polyA_gene_coverage_csvs = GET_GENE_COVERAGE_POLYA(polyA_terminal_sorted_full_bams_bais.mix(polyA_terminal_sorted_full_control_bams_bais), params.polyA_gene_coverage_csv, terminal_exons_bed, "polyA", params.gene_coverage_script)
		
		// combine all csv files from each sample into a single csv file.
		dedup_combined_gene_coverage_csv = COMBINE_GENE_COVERAGE_DEDUP(dedup_gene_coverage_csvs.collect(), params.combined_dedup_gene_coverage_csv)
		polyA_combined_gene_coverage_csv = COMBINE_GENE_COVERAGE_POLYA(polyA_gene_coverage_csvs.collect(), params.combined_polyA_gene_coverage_csv)
		
		// total number of genes in a given species 
		total_num_genes_csv = GET_TOTAL_NUM_GENES(terminal_exons_bed, params.total_num_genes_csv, params.sample_type, params.total_num_genes_script)
		
		// draw bar chart of gene coverage for all samples.
		(coverage_barchart_png, coverage_barchart_svg) = GET_BAR_CHART_GENE_COVERAGE(dedup_combined_gene_coverage_csv, polyA_combined_gene_coverage_csv,\
		total_num_genes_csv, params.geneCoverage_barchart, params.get_barchart_gene_coverage_script)

		//////////////////////////// PWM_logo ////////////////////////////
		// draw a logo of PWM of softclipped reads of length = 5, 6, 7, 8, 9, 10 separately.
		// This is to see the frequency of "A" appearing in each position of softclipped region.
		length_five_PWM_logo = SOFTCLIP_PWM_LOGO_LENGTH_FIVE(fixed_dedup_sorted_full_bams_bais, 5, params.softclip_pwm_logo_script)
		length_six_PWM_logo = SOFTCLIP_PWM_LOGO_LENGTH_SIX(fixed_dedup_sorted_full_bams_bais, 6, params.softclip_pwm_logo_script)
		length_seven_PWM_logo = SOFTCLIP_PWM_LOGO_LENGTH_SEVEN(fixed_dedup_sorted_full_bams_bais, 7, params.softclip_pwm_logo_script)
		length_eight_PWM_logo = SOFTCLIP_PWM_LOGO_LENGTH_EIGHT(fixed_dedup_sorted_full_bams_bais, 8, params.softclip_pwm_logo_script)
		length_nine_PWM_logo = SOFTCLIP_PWM_LOGO_LENGTH_NINE(fixed_dedup_sorted_full_bams_bais, 9, params.softclip_pwm_logo_script)
		length_ten_PWM_logo = SOFTCLIP_PWM_LOGO_LENGTH_TEN(fixed_dedup_sorted_full_bams_bais, 10, params.softclip_pwm_logo_script)

		// merged_cell_type
		if(params.overlap == "yes"){
			(polyA_bams, polyA_bais) = UNCOUPLE_BAM_WITH_SAMPLE(polyA_sorted_full_bams_bais_sample_only)
			(all_samples_polyA_bam, all_samples_polyA_bai) = MERGE_ALL_SAMPLES(polyA_bams.collect(), polyA_bais.collect())

			all_samples_polyA_bed = GET_POLYA_UNIQUE_CLEAVAGE_SITES_ALL_SAMPLES(all_samples_polyA_bam, all_samples_polyA_bai, params.polyA_unique_cs_bed_sample, params.get_polyA_unique_cs_script)

			all_samples_polyA_clustered_bed = PERFORM_CLUSTERING_ALL_SAMPLES(all_samples_polyA_bed, params.sample_cluster_out, params.cluster_pas_script)

			// Get overlap between SCINPAS PAS with catalog PAS
			(all_samples_catalog_overlap_barchart_png, all_samples_catalog_overlap_barchart_svg) = GET_OVERLAP(all_samples_polyA_clustered_bed, params.get_overlap_script)
			(top_score_PAS_catalog_overlap_barchart_png, top_score_PAS_catalog_overlap_barchart_svg) = GET_TOP_OVERLAP(all_samples_polyA_clustered_bed, params.get_top_overlap_script)
		}		
		
		// cell_type_specific
		if(params.cell_type_analysis == "yes"){
			// prepare cell type annotation file (metadta)
			if(params.sample_type == "mouse"){
				cell_type_metadata = Channel
					.fromPath(params.mouse_cell_metadata)

			}

			else if(params.sample_type == "human"){
				cell_type_metadata = Channel
					.fromPath(params.human_cell_metadata)

			}

			// make a txt file from metadata containing all CBs of a specific sample and a cell type (type1 and type2)
			cell_type_txts = SPLIT_CSV_BY_CELL_TYPE(specific_dataFolders.combine(cell_type_metadata), params.split_cell_type_script)	

			// split fixed deduplicated bam file by type1 and type 2 cell type for each sample
			(type1_bams, type1_bais, type2_bams, type2_bais) = SPLIT_SAMPLE_BY_TYPE(fixed_dedup_sorted_full_bams_bais.join(cell_type_txts, by: 2).view())

			// merge type1 bam files into a single bam file and type 2 bam files into a single bam file.
			(merged_type1_bam_bai, merged_type2_bam_bai) = MERGE_BY_CELL_TYPE(type1_bams.collect(), type1_bais.collect(), type2_bams.collect(), type2_bais.collect())

			(type1_polyA_full, type1_non_polyA_full) = GET_POLYA_TYPE1_SPECIFIC(merged_type1_bam_bai, params.out_polyA_type1, params.out_non_polyA_type1, params.getPolyA_advanced_script)
			(type2_polyA_full, type2_non_polyA_full) = GET_POLYA_TYPE2_SPECIFIC(merged_type2_bam_bai, params.out_polyA_type2, params.out_non_polyA_type2, params.getPolyA_advanced_script)

			type1_polyA_full_bams_bais = SORT_POLYA_TYPE1_SPECIFIC(type1_polyA_full)
			type2_polyA_full_bams_bais = SORT_POLYA_TYPE2_SPECIFIC(type2_polyA_full)

			type1_polyA_full_beds = POLYA_TYPE1_BAM_TO_BED(type1_polyA_full_bams_bais, params.polyA_unique_cs_bed_type1, params.get_polyA_unique_cs_script)
			type2_polyA_full_beds = POLYA_TYPE2_BAM_TO_BED(type2_polyA_full_bams_bais, params.polyA_unique_cs_bed_type2, params.get_polyA_unique_cs_script)
			
			type1_polyA_clustered_beds = PERFORM_CLUSTERING_CELL_TYPE1(type1_polyA_full_beds, params.polyA_cluster_out_type1, params.cluster_pas_script)
			type2_polyA_clustered_beds = PERFORM_CLUSTERING_CELL_TYPE2(type2_polyA_full_beds, params.polyA_cluster_out_type2, params.cluster_pas_script)

			// type1 classification (bed)
			type1_polyA_intergenic_clustered_beds = GENE_NO_INTERSECT_TYPE1(type1_polyA_clustered_beds, genes_bed, params.type1, params.intergenic_template)
			type1_all_genic_clustered_beds = GENE_INTERSECT_TYPE1(type1_polyA_clustered_beds, genes_bed, params.type1, params.all_genic_template)
			
			type1_polyA_intronic_clustered_beds = EXON_NO_INTERSECT_TYPE1(type1_all_genic_clustered_beds, exons_bed, params.type1, params.intronic_template)
			type1_all_exonic_clustered_beds = EXON_INTERSECT_TYPE1(type1_all_genic_clustered_beds, exons_bed, params.type1, params.all_exonic_template)
			
			type1_polyA_exonic_clustered_beds = TE_NO_INTERSECT_TYPE1(type1_all_exonic_clustered_beds, terminal_exons_bed, params.type1, params.exonic_template)
			type1_polyA_TE_clustered_beds = TE_INTERSECT_TYPE1(type1_all_exonic_clustered_beds, terminal_exons_bed, params.type1, params.te_template)

			type1_polyA_filtered_TE_clustered_beds = TE_FURTHER_INTERSECT_TYPE1(type1_polyA_TE_clustered_beds, further_filtered_terminal_exons, params.type1, params.filtered_te_template)

			// type2 classification (bed)
			type2_polyA_intergenic_clustered_beds = GENE_NO_INTERSECT_TYPE2(type2_polyA_clustered_beds, genes_bed, params.type2, params.intergenic_template)
			type2_all_genic_clustered_beds = GENE_INTERSECT_TYPE2(type2_polyA_clustered_beds, genes_bed, params.type2, params.all_genic_template)
			
			type2_polyA_intronic_clustered_beds = EXON_NO_INTERSECT_TYPE2(type2_all_genic_clustered_beds, exons_bed, params.type2, params.intronic_template)
			type2_all_exonic_clustered_beds = EXON_INTERSECT_TYPE2(type2_all_genic_clustered_beds, exons_bed, params.type2, params.all_exonic_template)
			
			type2_polyA_exonic_clustered_beds = TE_NO_INTERSECT_TYPE2(type2_all_exonic_clustered_beds, terminal_exons_bed, params.type2, params.exonic_template)
			type2_polyA_TE_clustered_beds = TE_INTERSECT_TYPE2(type2_all_exonic_clustered_beds, terminal_exons_bed, params.type2, params.te_template)

			type2_polyA_filtered_TE_clustered_beds = TE_FURTHER_INTERSECT_TYPE2(type2_polyA_TE_clustered_beds, further_filtered_terminal_exons, params.type2, params.filtered_te_template)
			
			// type1 type 2 distance
			type1_type2_d_csv = DISTANCE_SCATTER(type1_polyA_filtered_TE_clustered_beds, type2_polyA_filtered_TE_clustered_beds, further_filtered_terminal_exons, params.get_type1_type2_d_script) 
			t1_t2_scatter_plots = T1_T2_SCATTER(type1_type2_d_csv, params.get_type1_type2_dScatter_script)

			// type 1 motif plot
			(intergenic_motif_graphs_type1, intergenic_m_overlaid_type1_png, intergenic_m_overlaid_type1_svg, intergenic_motif_orders_type1, intergenic_motif_scores_type1) = GET_MOTIF_FREQ_PLOT_INTERGENIC_TYPE1(type1_polyA_intergenic_clustered_beds,\
			params.motif_search_out_template_intergenic, params.intergenic_annotated, gtf_peaks, params.type1, params.search_motif_script)	

			(intronic_motif_graphs_type1, intronic_m_overlaid_type1_png, intronic_m_overlaid_type1_svg, intronic_motif_orders_type1, intronic_motif_scores_type1) = GET_MOTIF_FREQ_PLOT_INTRONIC_TYPE1(type1_polyA_intronic_clustered_beds,\
			params.motif_search_out_template_intronic, params.intronic_annotated, gtf_peaks, params.type1, params.search_motif_script)	

			(exonic_motif_graphs_type1, exonic_m_overlaid_type1_png, exonic_m_overlaid_type1_svg, exonic_motif_orders_type1, exonic_motif_scores_type1) = GET_MOTIF_FREQ_PLOT_EXONIC_TYPE1(type1_polyA_exonic_clustered_beds,\
			params.motif_search_out_template_exonic, params.polyA_exonic_annotated, gtf_peaks, params.type1, params.search_motif_script)

			(te_motif_graphs_type1, te_m_overlaid_type1_png, te_m_overlaid_type1_svg, te_motif_orders_type1, te_motif_scores_type1) = GET_MOTIF_FREQ_PLOT_TE_TYPE1(type1_polyA_TE_clustered_beds,\
			params.motif_search_out_template_te, params.te_annotated, gtf_peaks, params.type1, params.search_motif_script)
			
			// type 1 # of all PAS
			type1_all_pas_num_polyA_csv = GET_TYPE1_NUM_POLYA_SITES_ALLPAS(type1_polyA_clustered_beds,\
			params.allPAS_num_polyA_csv, params.all_polyA_annotated, params.type1, params.get_numA_polyA_script)

			// type 1 # of intronic PAS
			type1_intronic_num_polyA_csv = GET_TYPE1_NUM_POLYA_SITES_INTRONIC(type1_polyA_intronic_clustered_beds,\
			params.intronic_num_polyA_csv, params.intronic_annotated, params.type1, params.get_numA_polyA_script)

			// type 1 # of exonic PAS
			type1_exonic_num_polyA_csv = GET_TYPE1_NUM_POLYA_SITES_EXONIC(type1_polyA_exonic_clustered_beds,\
			params.exonic_num_polyA_csv, params.polyA_exonic_annotated, params.type1, params.get_numA_polyA_script)

			// type 1 # of TE PAS
			type1_TE_num_polyA_csv = GET_TYPE1_NUM_POLYA_SITES_TE(type1_polyA_TE_clustered_beds,\
			params.te_num_polyA_csv, params.te_annotated, params.type1, params.get_numA_polyA_script)

			// type 2 motif plot
			(intergenic_motif_graphs_type2, intergenic_m_overlaid_type2_png, intergenic_m_overlaid_type2_svg, intergenic_motif_orders_type2, intergenic_motif_scores_type2) = GET_MOTIF_FREQ_PLOT_INTERGENIC_TYPE2(type2_polyA_intergenic_clustered_beds,\
			params.motif_search_out_template_intergenic, params.intergenic_annotated, gtf_peaks, params.type2, params.search_motif_script)	

			(intronic_motif_graphs_type2, intronic_m_overlaid_type2_png, intronic_m_overlaid_type2_svg, intronic_motif_orders_type2, intronic_motif_scores_type2) = GET_MOTIF_FREQ_PLOT_INTRONIC_TYPE2(type2_polyA_intronic_clustered_beds,\
			params.motif_search_out_template_intronic, params.intronic_annotated, gtf_peaks, params.type2, params.search_motif_script)	

			(exonic_motif_graphs_type2, exonic_m_overlaid_type2_png, exonic_m_overlaid_type2_svg, exonic_motif_orders_type2, exonic_motif_scores_type2) = GET_MOTIF_FREQ_PLOT_EXONIC_TYPE2(type2_polyA_exonic_clustered_beds,\
			params.motif_search_out_template_exonic, params.polyA_exonic_annotated, gtf_peaks, params.type2, params.search_motif_script)

			(te_motif_graphs_type2, te_m_overlaid_type2_png, te_m_overlaid_type2_svg, te_motif_orders_type2, te_motif_scores_type2) = GET_MOTIF_FREQ_PLOT_TE_TYPE2(type2_polyA_TE_clustered_beds,\
			params.motif_search_out_template_te, params.te_annotated, gtf_peaks, params.type2, params.search_motif_script)
			
			// type 2 # of all PAS
			type2_all_pas_num_polyA_csv = GET_TYPE2_NUM_POLYA_SITES_ALLPAS(type2_polyA_clustered_beds,\
			params.allPAS_num_polyA_csv, params.all_polyA_annotated, params.type2, params.get_numA_polyA_script)

			// type 2 # of intronic PAS
			type2_intronic_num_polyA_csv = GET_TYPE2_NUM_POLYA_SITES_INTRONIC(type2_polyA_intronic_clustered_beds,\
			params.intronic_num_polyA_csv, params.intronic_annotated, params.type2, params.get_numA_polyA_script)

			// type 2 # of exonic PAS
			type2_exonic_num_polyA_csv = GET_TYPE2_NUM_POLYA_SITES_EXONIC(type2_polyA_exonic_clustered_beds,\
			params.exonic_num_polyA_csv, params.polyA_exonic_annotated, params.type2, params.get_numA_polyA_script)

			// type 2 # of TE PAS
			type2_TE_num_polyA_csv = GET_TYPE2_NUM_POLYA_SITES_TE(type2_polyA_TE_clustered_beds,\
			params.te_num_polyA_csv, params.te_annotated, params.type2, params.get_numA_polyA_script)	
		}

	}
	
}