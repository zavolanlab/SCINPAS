#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// import and make alias. need to make alias because process can be invoked once and only once
include {PREPARE_SINGLE_LINKAGE; TRIM_NEG_CONTROL; SPLIT_PHASE1; DEDUP; MERGE;\
GET_PRO_INTRONIC_INTERGENIC; ISOLATE_POLYA_TERMINAL; MERGE_ALL_COUNTS; SPLIT_PER_SAMPLE; ONE_BY_ONE; MAKE_COUNT_SAMPLE; GET_BAR_CHART;\
GET_SPAN_SPLIT; GET_SPAN_HIST; GET_DISTANCE_ALTER_SAMPLE; GET_DISTANCE_ALTER_CONTROL; COMBINE_CONTROL_DISTANCES; GET_D_HISTOGRAM_ALTER; FIX_SOFTCLIPPED_REGION;\
COMBINE_MOTIF_ORDERS; GET_MOTIF_FREQ_GTF_PLOT; SPLIT_MOTIF_ORDER_PER_SAMPLE; MOTIF_ORDERS_ONE_BY_ONE; COMPUTE_SCORES_FOR_MOTIF_ALTER;\
GET_BAR_CHART_SCORE_MOTIF_POLYA; GET_SOFTCLIPPED_DISTRIBUTION; GET_NUM_POLYA_SITES; GET_BAR_CHART_NUM_POLYA; GET_TOTAL_NUM_GENES; GET_BAR_CHART_GENE_COVERAGE} from './processes_advanced'

include {PREPARE_IN_SLC as PREPARE_IN_SLC_DATA} from './processes_advanced'
include {PREPARE_IN_SLC as PREPARE_IN_SLC_CONTROL_DEDUP} from './processes_advanced'
include {PREPARE_IN_SLC as PREPARE_IN_SLC_CONTROL_RAW} from './processes_advanced'

include {GET_SPAN_RAW as GET_SPAN_RAW_SAMPLE} from './processes_advanced'
include {GET_SPAN_RAW as GET_SPAN_RAW_CONTROL} from './processes_advanced'

include {GET_SPAN_DEDUP as GET_SPAN_DEDUP_SAMPLE} from './processes_advanced'
include {GET_SPAN_DEDUP as GET_SPAN_DEDUP_CONTROL} from './processes_advanced'

include {GET_POLYA as GET_POLYA_SAMPLE_ONLY} from './processes_advanced'
include {GET_POLYA as GET_POLYA_NEG_CONTROL} from './processes_advanced'

include {CONCAT_SPAN as CONCAT_SPAN_RAW} from './processes_advanced'
include {CONCAT_SPAN as CONCAT_SPAN_DEDUP} from './processes_advanced'
include {CONCAT_SPAN as CONCAT_SPAN_SPLIT} from './processes_advanced'
include {GET_SPAN_HIST as GET_SPAN_HIST_RAW} from './processes_advanced'
include {GET_SPAN_HIST as GET_SPAN_HIST_DEDUP} from './processes_advanced'
include {GET_SPAN_HIST as GET_SPAN_HIST_SPLIT} from './processes_advanced'
include {CONVERT_TUPLE as CONVERT_TUPLE_DEDUP} from './processes_advanced'
include {CONVERT_TUPLE as CONVERT_TUPLE_SPLIT} from './processes_advanced'

include {SORT_PHASE1 as SORT_DEDUP} from './processes_advanced'
include {SORT_PHASE1 as SORT_SPLITTED} from './processes_advanced'
include {SORT_PHASE2 as SORT_FIXED_DEDUP} from './processes_advanced'
include {SORT_PHASE2 as SORT_POLYA_NEG_CONTROL} from './processes_advanced'
include {SORT_PHASE2 as SORT_POLYA_SAMPLE_ONLY} from './processes_advanced'
include {SORT_PHASE2 as SORT_NOGX} from './processes_advanced'
include {SORT_PHASE2 as SORT_NONPOLYA} from './processes_advanced'
include {SORT_PHASE2 as SORT_INTERNAL_P} from './processes_advanced'
include {SORT_PHASE2_ALTER as SORT_ABOVE} from './processes_advanced'
include {SORT_PHASE2_ALTER as SORT_BELOW} from './processes_advanced'
include {SORT_PHASE2_ALTER as SORT_TRIM_NEG_CONTROL} from './processes_advanced'

// polyA -> polyA_terminal
include {INTERSECT as INTERSECT_POLYA} from './processes_advanced'
include {SORT_AFTER_INTERSECT as SORT_AFTER_INTERSECT_POLYA} from './processes_advanced'
// polyA -> polyA_terminal (samples only)
include {INTERSECT as INTERSECT_POLYA_SAMPLE_ONLY} from './processes_advanced'
include {SORT_AFTER_INTERSECT as SORT_AFTER_INTERSECT_POLYA_SAMPLE_ONLY} from './processes_advanced'

// noGX -> noGX_noExon
include {NO_INTERSECT as NO_INTERSECT_NOGX} from './processes_advanced'
include {SORT_AFTER_INTERSECT as SORT_AFTER_NO_INTERSECT_NOGX} from './processes_advanced'
// noGX_noExon -> intronic
include {INTERSECT as INTERSECT_NOGX_NOEXON} from './processes_advanced'
include {SORT_AFTER_INTERSECT as SORT_AFTER_INTERSECT_NOGX_NOEXON} from './processes_advanced'
// noGX_noExon -> intergenic
include {NO_INTERSECT as NO_INTERSECT_NOGX_NOEXON} from './processes_advanced'
include {SORT_AFTER_INTERSECT as SORT_AFTER_NO_INTERSECT_NOGX_NOEXON} from './processes_advanced'
// polyA -> polyA_nonTerminal
include {NO_INTERSECT as NO_INTERSECT_POLYA} from './processes_advanced'
include {SORT_AFTER_INTERSECT as SORT_AFTER_NO_INTERSECT_POLYA} from './processes_advanced'
// polyA_nonTerminal -> polyA_exonic
include {INTERSECT as INTERSECT_POLYA_NONTERMINAL} from './processes_advanced'
include {SORT_AFTER_INTERSECT as SORT_AFTER_INTERSECT_POLYA_NONTERMINAL} from './processes_advanced'


include {GET_COUNTS as GET_COUNTS_RAW} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_RAW} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_DEDUP} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_DEDUP} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_NONPOLYA} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_NONPOLYA} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_INTERNAL_P} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_INTERNAL_P} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_POLYA} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_POLYA} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_POLYA_TERMINAL} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_POLYA_TERMINAL} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_INTRONIC} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_INTRONIC} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_INTERGENIC} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_INTERGENIC} from './processes_advanced'

include {GET_COUNTS as GET_COUNTS_EXONIC} from './processes_advanced'
include {CONCAT_COUNTS as CONCAT_COUNTS_EXONIC} from './processes_advanced'

include {SPLIT_PHASE2 as SPLIT_NEG_CONTROL_RAW} from './processes_advanced'
include {SPLIT_PHASE2 as SPLIT_NEG_CONTROL_DEDUP} from './processes_advanced'
include {SPLIT_PHASE2 as SPLIT_POLYA_TERMINAL} from './processes_advanced'
include {SPLIT_PHASE2 as SPLIT_INTRONIC} from './processes_advanced'
include {SPLIT_PHASE2 as SPLIT_INTERGENIC} from './processes_advanced'
include {SPLIT_PHASE2 as SPLIT_NONPOLYA} from './processes_advanced'
include {SPLIT_PHASE2 as SPLIT_INTERNAL_P} from './processes_advanced'
include {SPLIT_PHASE2 as SPLIT_POLYA_EXONIC} from './processes_advanced'
include {SPLIT_PHASE2 as SPLIT_POLYA_WITH_NEGCONTROL} from './processes_advanced'
include {PLOT_ATGC as PLOT_ATGC_ABOVE} from './processes_advanced'
include {PLOT_ATGC as PLOT_ATGC_BELOW} from './processes_advanced'
include {PLOT_ATGC as PLOT_ATGC_NONPOLYA} from './processes_advanced'
include {PLOT_ATGC as PLOT_ATGC_INTRONIC} from './processes_advanced'
include {PLOT_ATGC as PLOT_ATGC_INTERGENIC} from './processes_advanced'
include {PLOT_ATGC as PLOT_ATGC_INTERNAL_P} from './processes_advanced'
include {PLOT_ATGC as PLOT_ATGC_EXONIC} from './processes_advanced'

include {GET_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_ABOVE} from './processes_advanced'
include {GET_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_BELOW} from './processes_advanced'
include {GET_MOTIF_FREQ_PLOT as GET_MOTIF_FREQ_PLOT_NONPOLYA} from './processes_advanced'
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

include {GET_GENE_COVERAGE as GET_GENE_COVERAGE_DEDUP} from './processes_advanced'
include {GET_GENE_COVERAGE as GET_GENE_COVERAGE_POLYA} from './processes_advanced'

include {COMBINE_CSVS as COMBINE_SCORES_FOR_MOTIF} from './processes_advanced'
include {COMBINE_CSVS as COMBINE_NUM_POLYA} from './processes_advanced'
include {COMBINE_CSVS as COMBINE_GENE_COVERAGE_DEDUP} from './processes_advanced'
include {COMBINE_CSVS as COMBINE_GENE_COVERAGE_POLYA} from './processes_advanced'

workflow polyA{

	take: 
	specific_dataFolders
	prepare_singleLinkage
	single_linkage
	chromosomes
	terminal_exons_bed
	exons_bed
	genes_bed
	negative_cntrl_dedup
	negative_cntrl_raw
	num_chrom
	
	main:

	chromosomes.view()
	// phase 1

	// prepare folders for sample data
	(in_slc_Folder, in_slc_PWM_Folder, in_slc_motif_Folder, in_slc_num_polyA_Folder, full_data_tuple) = PREPARE_IN_SLC_DATA(specific_dataFolders, params.sample_type, prepare_singleLinkage)

	// prepare folders for negative control
	(neg_in_slc_raw, neg_in_slc_PWM_raw, neg_in_slc_motif_raw, neg_in_slc_num_polyA_dedup, negative_control_data_raw_tuple) = PREPARE_IN_SLC_CONTROL_RAW(negative_cntrl_raw, params.sample_type, prepare_singleLinkage)
	(neg_in_slc_dedup, neg_in_slc_PWM_dedup, neg_in_slc_motif_dedup, neg_in_slc_num_polyA_dedup, negative_control_data_dedup_tuple) = PREPARE_IN_SLC_CONTROL_DEDUP(negative_cntrl_dedup, params.sample_type, prepare_singleLinkage)
	
	// combine: cartesian product channel between full_data_tuple and chromosomes
	(possorted_bams_bais, bam_folders) = SPLIT_PHASE1(full_data_tuple.combine(chromosomes))
	
	// deduplication
	(dedup_bams, splitted_bams) = DEDUP(bam_folders, params.dedup_script)

	sorted_dedup_bams_bais = SORT_DEDUP(dedup_bams)
	
	(sorted_dedup_bams, sorted_dedup_bais) = CONVERT_TUPLE_DEDUP(sorted_dedup_bams_bais)
	
	splitted_sorted_bams_bais = SORT_SPLITTED(splitted_bams)
	(splitted_sorted_bams, splitted_sorted_bais) = CONVERT_TUPLE_SPLIT(splitted_sorted_bams_bais)
	// groupTuple allows to group file paths by samples.
	// by: 1 means group by 2nd element of the tuple (i.e. samples)
	// channel looks like: [([bam1, bam2...... bam6.....bamY], sample1), ([bam1, bam2...... bam6.....bamY], sample2)]
	// dont need to use .collect()
	// sorted_dedup_bams.groupTuple(by: 1, sort: true, size: 21).view()
	
	dedup_full_sorted_tuple = MERGE(sorted_dedup_bams.groupTuple(by: 1, sort: true, size: num_chrom), sorted_dedup_bais.groupTuple(by: 1, sort: true, size: num_chrom))

	// add 'FC' tag to negative control.
	fixed_negative_control_dedup_bam = TRIM_NEG_CONTROL(negative_control_data_dedup_tuple, params.trim_neg_control_script)
	fixed_negative_control_dedup_bam_bai = SORT_TRIM_NEG_CONTROL(fixed_negative_control_dedup_bam)

	// raw negative control reads per chr
	raw_negative_control_sorted_bams_bais = \
	SPLIT_NEG_CONTROL_RAW(negative_control_data_raw_tuple.combine(chromosomes), params.raw_negative_ctrl_template)

	// raw negative control reads per chr
	dedup_negative_control_sorted_bams_bais = \
	SPLIT_NEG_CONTROL_DEDUP(fixed_negative_control_dedup_bam_bai.combine(chromosomes), params.dedup_negative_ctrl_template)

	// draw span distribution to compare raw vs dedup vs split
	raw_spans_samples = GET_SPAN_RAW_SAMPLE(possorted_bams_bais, "raw", "sample", params.span_script)
	raw_spans_control = GET_SPAN_RAW_CONTROL(raw_negative_control_sorted_bams_bais, "raw", "control", params.span_script)
	raw_spans = raw_spans_samples.mix(raw_spans_control)

	dedup_spans_samples = GET_SPAN_DEDUP_SAMPLE(sorted_dedup_bams_bais, "dedup", "sample", params.span_script)
	dedup_spans_control = GET_SPAN_DEDUP_CONTROL(dedup_negative_control_sorted_bams_bais, "dedup", "control", params.span_script)
	dedup_spans = dedup_spans_samples.mix(dedup_spans_control)

	split_spans = GET_SPAN_SPLIT(splitted_sorted_bams_bais, "split", params.span_script)

	raw_span_merged = CONCAT_SPAN_RAW(raw_spans.groupTuple(by: 1, sort: true, size: num_chrom), params.raw_span)
	dedup_span_merged = CONCAT_SPAN_DEDUP(dedup_spans.groupTuple(by: 1, sort: true, size: num_chrom), params.dedup_span)
	split_span_merged = CONCAT_SPAN_SPLIT(split_spans.groupTuple(by: 1, sort: true, size: num_chrom), params.split_span)

	raw_span_hist = GET_SPAN_HIST_RAW(raw_span_merged, params.span_histogram_script, params.raw_span_hist)
	dedup_span_hist = GET_SPAN_HIST_DEDUP(dedup_span_merged, params.span_histogram_script, params.dedup_span_hist)
	split_span_hist = GET_SPAN_HIST_SPLIT(split_span_merged, params.span_histogram_script, params.split_span_hist)
	
	// phase 2

	//fix softclipped region (Alignment fixation) in the dedup file
	fixed_dedup_full_bams = FIX_SOFTCLIPPED_REGION(dedup_full_sorted_tuple, params.fix_softclipped_script)
	fixed_dedup_sorted_full_bams_bais = SORT_FIXED_DEDUP(fixed_dedup_full_bams)
	
	// get all full polyA reads and full non polyA reads without negative control.
	// This is for the motif scoring. You want to compare score between negative control (UMI-tool deduplicated) and those sample data which underwent the pipeline.
	(polyA_sample_only, non_polyA_sample_only) = GET_POLYA_SAMPLE_ONLY(fixed_dedup_sorted_full_bams_bais, params.out_polyA, params.out_non_polyA, params.getPolyA_advanced_script)
	
	// get all full polyA reads and full non polyA reads of negative control.
	(polyA_control, non_polyA_control) = GET_POLYA_NEG_CONTROL(fixed_negative_control_dedup_bam_bai, params.out_polyA_control, params.out_non_polyA_control, params.getPolyA_advanced_script)

	// sort_phase1, sort_phase2 can both do sorting and indexing (for only 1 file) and (for all files) because of tuple
	// if for all files: (sorted_bam1, sorted_bam1, specific_sample6), (sorted_bam2, sorted_bam2, specific_sample6), (sorted_bam5, sorted_bam5, specific_sample6), (sorted_bam6, sorted_bam6, specific_sample6)
	// if for only 1 file: (sorted_full_bam, sorted_full_bai, specific_sample5), (sorted_full_bam', sorted_full_bai', specific_sample6)
	// sort full polyA reads for all samples
	polyA_sorted_full_bams_bais_sample_only = SORT_POLYA_SAMPLE_ONLY(polyA_sample_only)
	polyA_sorted_full_bams_bais_control = SORT_POLYA_NEG_CONTROL(polyA_control)

	// get polyA reads mapping to terminal exon (sample only)
	// This is to compute distance between read end and end of terminal exon.
	// For this purpose, samples should be polyA reads mapping to terminal exons.
	// For this purpose, negative control should not be filtered for polyA.
	polyA_terminal_full_sample_only_bams = \
	INTERSECT_POLYA_SAMPLE_ONLY(polyA_sorted_full_bams_bais_sample_only, terminal_exons_bed, "polyA_reads_terminal_fullSampleOnly")

	polyA_terminal_sorted_full_sample_only_bams_bais = \
	SORT_AFTER_INTERSECT_POLYA_SAMPLE_ONLY(polyA_terminal_full_sample_only_bams, "polyA_reads_terminal_fullSampleOnly")

	// get polyA reads mapping to terminal exon (sample + negative control)
	polyA_terminal_full_bams = \
	INTERSECT_POLYA(polyA_sorted_full_bams_bais_sample_only.mix(polyA_sorted_full_bams_bais_control), terminal_exons_bed, "polyA_reads_terminal_full")

	polyA_terminal_sorted_full_bams_bais = \
	SORT_AFTER_INTERSECT_POLYA(polyA_terminal_full_bams, "polyA_reads_terminal_full")

	// split polyA reads mapping to terminal exon by chr, sorting and indexing for all chromosomes. repeat this for all samples
	polyA_terminal_sorted_bams_bais = \
	SPLIT_POLYA_TERMINAL(polyA_terminal_sorted_full_bams_bais.combine(chromosomes), "polyA_reads_terminal")

	// split polyA + terminal reads (chr1~Y) by distance threshold.
	(isolated_above_bams, isolated_below_bams)= \
	ISOLATE_POLYA_TERMINAL(polyA_terminal_sorted_bams_bais, terminal_exons_bed, params.isolate_polyATerminal_script)
	
	// get full noGX polyA reads.
	noGX_polyA_full_bams = \
	GET_PRO_INTRONIC_INTERGENIC(polyA_sorted_full_bams_bais_sample_only.mix(polyA_sorted_full_bams_bais_control), params.get_pro_intronic_intergenic_script)
	
	// sort full noGX_polyA reads
	noGX_polyA_sorted_full_bams_bais = SORT_NOGX(noGX_polyA_full_bams)

	// noGX_reads that do not overlap with exons.bed will be either intronic or intergenic
	noGX_noExon_full_bams = \
	NO_INTERSECT_NOGX(noGX_polyA_sorted_full_bams_bais, exons_bed, "noGX_noExon_polyA")	

	noGX_noExon_sorted_full_bams_bais = \
	SORT_AFTER_NO_INTERSECT_NOGX(noGX_noExon_full_bams, "noGX_noExon_polyA")

	// if noGX_noExon reads overlap with genes.bed, then they will be intronic reads
	intronic_full_bams = \
	INTERSECT_NOGX_NOEXON(noGX_noExon_sorted_full_bams_bais, genes_bed, "intronic_reads_full")

	intronic_sorted_full_bams_bais = \
	SORT_AFTER_INTERSECT_NOGX_NOEXON(intronic_full_bams, "intronic_reads_full")

	// if noGX_noExon reads do not overlap with genes.bed, then they will be intergenic reads
	intergenic_full_bams = \
	NO_INTERSECT_NOGX_NOEXON(noGX_noExon_sorted_full_bams_bais, genes_bed, "intergenic_reads_full")

	intergenic_sorted_full_bams_bais = \
	SORT_AFTER_NO_INTERSECT_NOGX_NOEXON(intergenic_full_bams, "intergenic_reads_full")
	
	// sort nonpolyA reads
	nonpolyA_sorted_full_bams_bais = SORT_NONPOLYA(non_polyA_sample_only.mix(non_polyA_control))

	// polyA reads that do not overlap with terminal exons but overlap with exons will be exonic reads
	polyA_nonTerminal_full_bams = \
	NO_INTERSECT_POLYA(polyA_sorted_full_bams_bais_sample_only.mix(polyA_sorted_full_bams_bais_control), terminal_exons_bed, "polyA_non_terminal_full")	

	polyA_nonTerminal_sorted_full_bams_bais = \
	SORT_AFTER_NO_INTERSECT_POLYA(polyA_nonTerminal_full_bams, "polyA_non_terminal_full")	

	polyA_exonic_full_bams = \
	INTERSECT_POLYA_NONTERMINAL(polyA_nonTerminal_sorted_full_bams_bais, exons_bed, "polyA_exonic_full")

	polyA_exonic_sorted_full_bams_bais = \
	SORT_AFTER_INTERSECT_POLYA_NONTERMINAL(polyA_exonic_full_bams, "polyA_exonic_full")

	// sort isolated reads (above, below) for all chromosomes (1~19, X, Y)
	isolated_above_sorted_bams_bais = SORT_ABOVE(isolated_above_bams)
	isolated_below_sorted_bams_bais = SORT_BELOW(isolated_below_bams)
	
	// intronic reads per chr
	intronic_sorted_bams_bais = \
	SPLIT_INTRONIC(intronic_sorted_full_bams_bais.combine(chromosomes), params.intronic_template)

	// intergenic reads per chr
	intergenic_sorted_bams_bais = \
	SPLIT_INTERGENIC(intergenic_sorted_full_bams_bais.combine(chromosomes), params.intergenic_template)

	// nonpolyA reads per chr
	nonpolyA_sorted_bams_bais = \
	SPLIT_NONPOLYA(nonpolyA_sorted_full_bams_bais.combine(chromosomes), params.non_polyA_template)
	
	// exonic polyA reads per chr
	polyA_exonic_sorted_bams_bais = \
	SPLIT_POLYA_EXONIC(polyA_exonic_sorted_full_bams_bais.combine(chromosomes), params.polyA_exonic_template)

	// phase 3

	// get counts csvs per "raw" class and per sample name (including negative control).
	raw_counts_csvs = GET_COUNTS_RAW(full_data_tuple.mix(negative_control_data_raw_tuple), params.counts_raw, params.raw_type, params.get_count_script)
	
	// merges counts csvs per "raw" class across all samples (including negative control).
	merged_raw_counts = CONCAT_COUNTS_RAW(raw_counts_csvs.collect(), params.raw_type)
	
	// get counts csvs per "dedup" class and per sample name (including negative control).
	dedup_counts_csvs = GET_COUNTS_DEDUP(fixed_dedup_sorted_full_bams_bais.mix(fixed_negative_control_dedup_bam_bai), params.counts_dedup, params.dedup_type, params.get_count_script)
	// merges counts csvs per "dedup" across all samples (including negative control).
	merged_dedup_counts = CONCAT_COUNTS_DEDUP(dedup_counts_csvs.collect(), params.dedup_type)

	// get counts csvs per non-polyA class and per sample name (including negative control).
	nonpolyA_counts_csvs = GET_COUNTS_NONPOLYA(nonpolyA_sorted_full_bams_bais, params.counts_nonpolyA, params.nonpolyA_type, params.get_count_script)
	// merges counts csvs per non-polyA across all samples (including negative control).
	merged_nonpolyA_counts = CONCAT_COUNTS_NONPOLYA(nonpolyA_counts_csvs.collect(), params.nonpolyA_type)

	// get counts csvs per polyA class and per sample name (including negative control).
	polyA_counts_csvs = GET_COUNTS_POLYA(polyA_sorted_full_bams_bais_sample_only.mix(polyA_sorted_full_bams_bais_control), params.counts_polyA, params.polyA_type, params.get_count_script)
	// merges counts csvs per polyA class across all samples (including negative control).
	merged_polyA_counts = CONCAT_COUNTS_POLYA(polyA_counts_csvs.collect(), params.polyA_type)
	
	// get counts csvs per polyA + terminal class and per sample name (including negative control).
	polyAT_counts_csvs = GET_COUNTS_POLYA_TERMINAL(polyA_terminal_sorted_full_bams_bais, params.counts_polyAT, params.polyAT_type, params.get_count_script)
	// merges counts csvs per polyA + terminal class across all samples (including negative control).
	merged_polyAT_counts = CONCAT_COUNTS_POLYA_TERMINAL(polyAT_counts_csvs.collect(), params.polyAT_type)
	
	// get counts csvs per intronic class and per sample name (including negative control).
	intronic_counts_csvs = GET_COUNTS_INTRONIC(intronic_sorted_full_bams_bais, params.counts_intronic, params.intronic_type, params.get_count_script)
	// merges counts csvs per intronic class across all samples (including negative control).
	merged_intronic_counts = CONCAT_COUNTS_INTRONIC(intronic_counts_csvs.collect(), params.intronic_type)
	
	// get counts csvs per intergenic class and per sample name (including negative control).
	intergenic_counts_csvs = GET_COUNTS_INTERGENIC(intergenic_sorted_full_bams_bais, params.counts_intergenic, params.intergenic_type, params.get_count_script)
	// merges counts csvs per intergenic class across all samples (including negative control).
	merged_intergenic_counts = CONCAT_COUNTS_INTERGENIC(intergenic_counts_csvs.collect(), params.intergenic_type)
	
	// get counts csvs per exonic class and per sample name (including negative control).
	exonic_counts_csvs = GET_COUNTS_EXONIC(polyA_exonic_sorted_full_bams_bais, params.counts_exonic, params.exonic_type, params.get_count_script)
	// merges counts csvs per exonic class across all samples (including negative control).
	merged_exonic_counts = CONCAT_COUNTS_EXONIC(exonic_counts_csvs.collect(), params.exonic_type)

	// merges counts csvs across all classes and across all samples
	total_counts_csv = MERGE_ALL_COUNTS(merged_raw_counts, merged_dedup_counts, merged_nonpolyA_counts,\
	merged_polyA_counts, merged_polyAT_counts, merged_intronic_counts, merged_intergenic_counts, merged_exonic_counts)

	// split counts csv by sample
	split_by_samples_counts_csvs_collected = SPLIT_PER_SAMPLE(total_counts_csv, params.get_count_split_script)

	split_by_samples_counts_csvs = ONE_BY_ONE(split_by_samples_counts_csvs_collected)
	
	// make tuple of file path and sample name
	counts_and_samples = split_by_samples_counts_csvs.map { file -> [file, file.baseName] }

	// make counts bar chart
	barcharts = GET_BAR_CHART(counts_and_samples, params.get_barchart_script)
	
	// phase 4

	// get distance distribution between end of terminal exon and end of a read.
	(distance_outputs_alter_tuple, outliers) = GET_DISTANCE_ALTER_SAMPLE(polyA_terminal_sorted_full_sample_only_bams_bais,\
	terminal_exons_bed, "sample", params.get_distance_final_script)

	// get distance distribution between end of terminal exon and end of a control. (per chromosome basis)
	(distance_outputs_alter_control, outliers_control) = GET_DISTANCE_ALTER_CONTROL(dedup_negative_control_sorted_bams_bais,\
	terminal_exons_bed, "control", params.get_distance_final_script)

	// combine negative control distance chr1 ~ 19 + X + Y.
	combined_negative_control_output = COMBINE_CONTROL_DISTANCES(distance_outputs_alter_control.groupTuple(by: 1, sort: true, size: num_chrom), params.combined_negative_control_distance_output)
	
	//distance histogram (before and after fixing cleavage site)
	distance_histograms_alter = GET_D_HISTOGRAM_ALTER(distance_outputs_alter_tuple.combine(combined_negative_control_output), params.get_histogram_final_script)

	// phase 5

	// plot frequency plot in each class
	// sorted_bams_bais are tuple (bams, bais, samples). you need to group the "keys" by samples. hence it should be by: 2 (0-based indexing)
	// data structure: [[bam1~bamY], [bai1~baiY], sample_name]

	above_atgc_plts = PLOT_ATGC_ABOVE(isolated_above_sorted_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom), params.unannotated_in_template,\
	params.unannotated_ATGC_out, params.above_annotated, single_linkage, params.sample_type, params.plot_frequencies_script)
	// isolated_above_sorted_bams_bais.groupTuple(by: 2, sort: true, size: 21).view()

	below_atgc_plts = PLOT_ATGC_BELOW(isolated_below_sorted_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom), params.annotated_in_template,\
	params.annotated_ATGC_out, params.below_annotated, single_linkage, params.sample_type, params.plot_frequencies_script)

	nonpolyA_atgc_plts = PLOT_ATGC_NONPOLYA(nonpolyA_sorted_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom), params.non_polyA_template,\
	params.non_polyA_ATGC_out, params.non_polyA_annotated, single_linkage, params.sample_type, params.plot_frequencies_script)

	intergenic_atgc_plts = PLOT_ATGC_INTERGENIC(intergenic_sorted_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom), params.intergenic_template,\
	params.intergenic_ATGC_out, params.intergenic_annotated, single_linkage, params.sample_type, params.plot_frequencies_script)

	intronic_atgc_plts = PLOT_ATGC_INTRONIC(intronic_sorted_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom), params.intronic_template,\
	params.intronic_ATGC_out, params.intronic_annotated, single_linkage, params.sample_type, params.plot_frequencies_script)

	exonic_atgc_plt = PLOT_ATGC_EXONIC(polyA_exonic_sorted_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom), params.polyA_exonic_template,\
	params.polyA_exonic_ATGC_out, params.polyA_exonic_annotated, single_linkage, params.sample_type, params.plot_frequencies_script)

	// phase 6 (motif frequency plots)

	// draw 18 motives frequency plots based on gtf file
	(gtf_motif_graphs, gtf_motif_graphs_overlaid, gtf_motif_order, gtf_peaks) = GET_MOTIF_FREQ_GTF_PLOT(terminal_exons_bed, params.motif_search_out_gtf, single_linkage, params.gtf_search_motif_script)
	
	// draw 18 motives frequency plots from our data in each class in each sample

	(annotated_motif_graphs, annotated_motif_graphs_overlaid, annotated_motif_orders, annotated_motif_scores) = GET_MOTIF_FREQ_PLOT_BELOW(isolated_below_sorted_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom), params.annotated_in_template,\
	params.motif_search_out_template_annotated, params.below_annotated, gtf_peaks, single_linkage, params.sample_type, params.search_motif_script)

	(unannotated_motif_graphs, unannotated_motif_graphs_overlaid, unannotated_motif_orders, unannotated_motif_scores) = GET_MOTIF_FREQ_PLOT_ABOVE(isolated_above_sorted_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom), params.unannotated_in_template,\
	params.motif_search_out_template_unannotated, params.above_annotated, gtf_peaks, single_linkage, params.sample_type, params.search_motif_script)
		
	(intergenic_motif_graphs, intergenic_motif_graphs_overlaid, intergenic_motif_orders, intergenic_motif_scores) = GET_MOTIF_FREQ_PLOT_INTERGENIC(intergenic_sorted_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom), params.intergenic_template,\
	params.motif_search_out_template_intergenic, params.intergenic_annotated, gtf_peaks, single_linkage, params.sample_type, params.search_motif_script)	

	(intronic_motif_graphs, intronic_motif_graphs_overlaid, intronic_motif_orders, intronic_motif_scores) = GET_MOTIF_FREQ_PLOT_INTRONIC(intronic_sorted_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom), params.intronic_template,\
	params.motif_search_out_template_intronic, params.intronic_annotated, gtf_peaks, single_linkage, params.sample_type, params.search_motif_script)	

	(exonic_motif_graphs, exonic_motif_graphs_overlaid, exonic_motif_orders, exonic_motif_scores) = GET_MOTIF_FREQ_PLOT_EXONIC(polyA_exonic_sorted_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom), params.polyA_exonic_template,\
	params.motif_search_out_template_exonic, params.polyA_exonic_annotated, gtf_peaks, single_linkage, params.sample_type, params.search_motif_script)

	// save 18 motif scores in each class and in each sample
	// save 18 motif orders in each class and in each sample
	(merged_annotated_motif_orders, merged_annotated_motif_scores) = CONCAT_MOTIF_ORDERS_ANNOTATED(annotated_motif_orders.collect(), annotated_motif_scores.collect(), params.motif_search_out_template_annotated)

	(merged_unannotated_motif_orders, merged_unannotated_motif_scores)  = CONCAT_MOTIF_ORDERS_UNANNOTATED(unannotated_motif_orders.collect(), unannotated_motif_scores.collect(), params.motif_search_out_template_unannotated)
	
	(merged_intergenic_motif_orders, merged_intergenic_motif_scores) = CONCAT_MOTIF_ORDERS_INTERGENIC(intergenic_motif_orders.collect(), intergenic_motif_scores.collect(), params.motif_search_out_template_intergenic)
	
	(merged_intronic_motif_orders, merged_intronic_motif_scores) = CONCAT_MOTIF_ORDERS_INTRONIC(intronic_motif_orders.collect(), intronic_motif_scores.collect(), params.motif_search_out_template_intronic)
	
	(merged_exonic_motif_orders, merged_exonic_motif_scores) = CONCAT_MOTIF_ORDERS_EXONIC(exonic_motif_orders.collect(), exonic_motif_scores.collect(), params.motif_search_out_template_exonic)

	// merge motif orders across all samples
	total_motif_orders_csv = MERGE_ALL_MOTIF_ORDERS(merged_annotated_motif_orders, merged_unannotated_motif_orders, merged_intergenic_motif_orders, merged_intronic_motif_orders, merged_exonic_motif_orders, params.total_merged_motif_orders)
	total_motif_scores_csv = MERGE_ALL_MOTIF_SCORES(merged_annotated_motif_scores, merged_unannotated_motif_scores, merged_intergenic_motif_scores, merged_intronic_motif_scores, merged_exonic_motif_scores, params.total_merged_motif_scores)

	// split motif orders by each sample. Each sample has motif orders for annotated, unannotated, intronic, intergenic, exonic.
	split_by_samples_motif_orders_csv_collected = SPLIT_MOTIF_ORDER_PER_SAMPLE(total_motif_orders_csv, params.split_motif_order_script)
	split_by_samples_motif_orders_csv = MOTIF_ORDERS_ONE_BY_ONE(split_by_samples_motif_orders_csv_collected)
	
	// make tuple of file path and sample name
	motif_orders_and_samples = split_by_samples_motif_orders_csv.map { file -> [file, file.baseName] }

	// combine motif orders of 1 sample + motif order of gtf. Each sample has motif orders for annotated, unannotated, intronic, intergenic, exonic and gtf.
	combined_motif_order = COMBINE_MOTIF_ORDERS(motif_orders_and_samples, gtf_motif_order, params.motif_order_outname)
	
	// split polyA by chromosome (negative control also included). This is for efficient computing in process "COMPUTE_SCORES_FOR_MOTIF_ALTER"
	
	// The mix operator combines the items emitted by two (or more) channels into a single channel.
	// negative control is not filtered for polyA and directly fed into the rest of the pipeline.
	// name of negative control is e.g. A10X_P7_14UmiDedup_polyA_sorted_6.bam but they are not polyA reads.
	// This is because we want to compare out pipeline to really a negative control.
	polyA_with_negative_control_sorted_bams_bais = SPLIT_POLYA_WITH_NEGCONTROL\
	(polyA_sorted_full_bams_bais_sample_only.mix(fixed_negative_control_dedup_bam_bai).combine(chromosomes), params.polyA_with_control_template)
	
	// compute motif scores for 18 motives at the pulled polyA reads level.

	// polyA_with_negative_control_sorted_bams_bais.view()
	// groupTuple allows to group file paths by samples.
	// by: 2 means group by 3rd element of the tuple (i.e. samples) (0-based index)
	// channel looks like: [([bam1, bam2...... bam6.....bamY], sample1), ([bam1, bam2...... bam6.....bamY], sample2)]
	// polyA_with_negative_control_sorted_bams_bais.groupTuple(by: 2, sort: true, size: 21).view()
	(polyA_motif_graphs, polyA_motif_graphs_overlaid, polyA_scores_csv) = COMPUTE_SCORES_FOR_MOTIF_ALTER(polyA_with_negative_control_sorted_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom),\
	params.polyA_with_control_template, params.motif_search_out_polyA, gtf_peaks, single_linkage, params.sample_type, params.compute_motif_score_advanced_script)

	// combine 18 motif scores for each sample across all samples
	combined_motif_scores = COMBINE_SCORES_FOR_MOTIF(polyA_scores_csv.collect(), params.combined_motif_scores)
	polyA_motif_scores = GET_BAR_CHART_SCORE_MOTIF_POLYA(combined_motif_scores, params.polyA_motif_score_script)

	// plot histogram of length of softclipped region. (before and after fixation)	
	size_softclipped_dist = GET_SOFTCLIPPED_DISTRIBUTION(fixed_dedup_sorted_full_bams_bais, params.sizeSoftclipped_out, params.get_softclipped_script)

	// phase 7 (number of polyA sites identified)

	// get number of polyA sites in each sample.
	(num_polyA_csvs, num_polyA_distributions) = GET_NUM_POLYA_SITES(polyA_sorted_full_bams_bais_sample_only.mix(polyA_sorted_full_bams_bais_control), params.num_polyA_csv, params.num_polyA_barchart, single_linkage, params.get_numA_polyA_script)
	
	// combine all csv files from each sample into a single csv file.
	combined_numPolyA_csv = COMBINE_NUM_POLYA(num_polyA_csvs.collect(), params.combined_num_polyA_csv)
	
	// plot a barchart of number of polyA sites for all samples.
	num_polyA_barchart = GET_BAR_CHART_NUM_POLYA(combined_numPolyA_csv, params.num_polyA_barchart, params.get_barchart_numpolyA_script)

	// phase 8 (number of genes covered by each sample)

	// compute number of genes expressed by each sample and save it in each csv file. (includes negative control)
	dedup_gene_coverage_csvs = GET_GENE_COVERAGE_DEDUP(fixed_dedup_sorted_full_bams_bais.mix(fixed_negative_control_dedup_bam_bai), params.dedup_gene_coverage_csv, terminal_exons_bed, "dedup", params.gene_coverage_script)
	// compute number of genes covered by each sample and save it in each csv file. (includes negative control)
	polyA_gene_coverage_csvs = GET_GENE_COVERAGE_POLYA(polyA_terminal_sorted_full_bams_bais, params.polyA_gene_coverage_csv, terminal_exons_bed, "polyA", params.gene_coverage_script)
	
	// combine all csv files from each sample into a single csv file.
	dedup_combined_gene_coverage_csv = COMBINE_GENE_COVERAGE_DEDUP(dedup_gene_coverage_csvs.collect(), params.combined_dedup_gene_coverage_csv)
	polyA_combined_gene_coverage_csv = COMBINE_GENE_COVERAGE_POLYA(polyA_gene_coverage_csvs.collect(), params.combined_polyA_gene_coverage_csv)

	total_num_genes_csv = GET_TOTAL_NUM_GENES(terminal_exons_bed, params.total_num_genes_csv, params.sample_type, params.total_num_genes_script)
	
	// draw bar chart of gene coverage for all samples.
	coverage_barchart = GET_BAR_CHART_GENE_COVERAGE(dedup_combined_gene_coverage_csv, polyA_combined_gene_coverage_csv,\
	total_num_genes_csv, params.geneCoverage_barchart, params.get_barchart_gene_coverage_script)
}
