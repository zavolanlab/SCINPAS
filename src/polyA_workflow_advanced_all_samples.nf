#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//import and make alias. need to make alias because process can be invoked once and only once
include {PREPARE_IN_SLC_CATALOG; FILTER_MULTIMAPPING_CATALOG; CONCAT_CSVS_CATALOG; MODIFY_TUPLE;\
SPLIT_PHASE1_CATALOG; DEDUP_CATALOG; MERGE_DEDUP_CATALOG; FASTQC_CATALOG; FASTQC_SWARM_PLOT_CATALOG; FIX_SOFTCLIPPED_REGION_CATALOG;\
MERGE_POLYA_CATALOG; GET_COUNTS_CATALOG; GET_POLYA_UNIQUE_CLEAVAGE_SITES_CATALOG; SPLIT_BY_DIRECTION; GROUPBY_BED_CATALOG; PERFORM_CLUSTERING_CATALOG; GET_INTRONIC_BED;\
CHANGE_BED; ADD_CLASS_COLUMN; GENE_ID; LEFT_JOIN_CATALOG; GET_ORGAN_SCORE; MERGE_ORGAN_SCORE_TO_PAS; RCS_MOTIF_CHECK;
CONVERT_GZIP_UNIQUE_CS; CONVERT_GZIP_ALL_SAMPLES_UNIQUE_CS} from './processes_all_samples'

include {SORT_PHASE1_CATALOG as SORT_DEDUP} from './processes_all_samples'
include {SORT_PHASE2_CATALOG as SORT_FIXED_DEDUP} from './processes_all_samples'

include {GET_POLYA_CATALOG as GET_POLYA_SAMPLE_ONLY} from './processes_all_samples'

include {SORT_PHASE2_CATALOG as SORT_POLYA_SAMPLE_ONLY} from './processes_all_samples'
include {SPLIT_PHASE2_CATALOG as SPLIT_ALL_POLYA} from './processes_all_samples'

include {BED_INTERSECT_CATALOG as GENES_INTERSECT_ALL_SAMPLES} from './processes_all_samples'
include {BED_NO_INTERSECT_CATALOG as GENES_NO_INTERSECT_ALL_SAMPLES} from './processes_all_samples'

include {BED_INTERSECT_CATALOG as INTRONS_INTERSECT_ALL_SAMPLES} from './processes_all_samples'
include {BED_NO_INTERSECT_CATALOG as INTRONS_NO_INTERSECT_ALL_SAMPLES} from './processes_all_samples'

include {BED_INTERSECT_CATALOG as TE_INTERSECT_ALL_SAMPLES} from './processes_all_samples'
include {BED_NO_INTERSECT_CATALOG as TE_NO_INTERSECT_ALL_SAMPLES} from './processes_all_samples'

include {BED_INTERSECT_CATALOG_ANTISENSE as EXTRA_STEP_ONE_INTERSECT_ALL_SAMPLES} from './processes_all_samples'
include {BED_NO_INTERSECT_CATALOG_ANTISENSE as EXTRA_STEP_ONE_NO_INTERSECT_ALL_SAMPLES} from './processes_all_samples'

include {BED_INTERSECT_CATALOG_ANTISENSE as EXTRA_STEP_TWO_INTERSECT_ALL_SAMPLES} from './processes_all_samples'
include {BED_NO_INTERSECT_CATALOG_ANTISENSE as EXTRA_STEP_TWO_NO_INTERSECT_ALL_SAMPLES} from './processes_all_samples'

include {BED_INTERSECT_CATALOG_ANTISENSE as EXTRA_STEP_THREE_INTERSECT_ALL_SAMPLES} from './processes_all_samples'
include {BED_NO_INTERSECT_CATALOG_ANTISENSE as EXTRA_STEP_THREE_NO_INTERSECT_ALL_SAMPLES} from './processes_all_samples'

include {MERGE_ADD_COL_BED as MERGE_ADD_COL_BEDS} from './processes_all_samples'

include {CONVERT_GZIP_CATALOG as CONVERT_GZIP_CATALOG_GENEID} from './processes_all_samples'
include {CONVERT_GZIP_CATALOG as CONVERT_GZIP_CATALOG_GENEID_ORGAN_MOTIF} from './processes_all_samples'

workflow polyA_all_samples{

	take: 
	organs
	chromosomes
	our_terminal_exons
	our_exons
	our_genes
	num_chrom

	main:

	//////////////////////////// trimming_and_deduplication ////////////////////////////

	full_data_tuple = PREPARE_IN_SLC_CATALOG(organs, params.sample_type)
	// tuple = (bam, bai, sample, organ)
	(mapq_filtered_full_data_tuple, readInfo_csvs) = FILTER_MULTIMAPPING_CATALOG(full_data_tuple)
	// 1st filtering: filtering of samples manually by % uniquely mapped reads
	if(params.check == "yes"){
		(readInfo_merged_total, readInfo_merged_good, readInfo_merged_bad, first_filtered_input) = CONCAT_CSVS_CATALOG(readInfo_csvs.collect(), params.mapped_unmapped, params.merge_mapped_unmapped_script)

		first_filtered_samples = first_filtered_input
			.splitCsv(header: true)
			.map{row -> tuple(row.percentage, row.dir, row.sample, row.organ)}
		
		// resulting tuple = (sample, organ, bam, bai, percentage, dir)
		// join allows to only keep samples that match (and hence filtering) 
		first_filtered_full_data = mapq_filtered_full_data_tuple.join(first_filtered_samples, by: [2, 3])
		// first_filtered_full_data.view()
		
		n_first_filtered_samples = first_filtered_full_data.count()
		n_first_filtered_samples.view()

		// resulting tuple = (bam, bai, sample, organ)
		first_filtered_modified = MODIFY_TUPLE(first_filtered_full_data)

		// combine: cartesian product channel between full_data_tuple and chromosomes
		// first_filtered_modified.combine(chromosomes) is a channel of tuple:
		// ((bam, bai, sample, organ) + chromosome (cartesian joint)
		// resulting tuple = (bam, bai, sample, organ, chromosome)
		possorted_bams_bais = SPLIT_PHASE1_CATALOG(first_filtered_modified.combine(chromosomes))

		// deduplication of 1st filtered "samples"
		dedup_bams = DEDUP_CATALOG(possorted_bams_bais, params.dedup_script)
		sorted_dedup_bams_bais = SORT_DEDUP(dedup_bams)
		dedup_bams_bais_full = MERGE_DEDUP_CATALOG(sorted_dedup_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom))

		fastqc_out_zips = FASTQC_CATALOG(dedup_bams_bais_full)
		// 2nd filtering: filter samples by sequencing quality (fastqc)
		// you filtered good samples in seq_quality_filtered_csv within python script
		// the threshold in here should be checked per species before running the whole workflow. it is not universal.
		(seq_quality_full_csv, swarm_plot, seq_quality_filtered_csv, final_input_samples) = FASTQC_SWARM_PLOT_CATALOG(fastqc_out_zips.collect(), first_filtered_input, params.swarm_script)
	}

	else if(params.check == "no"){
		// mapq_filtered_full_data_tuple is already filtered at this point
		n_second_filtered_samples = mapq_filtered_full_data_tuple.count()
		n_second_filtered_samples.view()
		possorted_bams_bais = SPLIT_PHASE1_CATALOG(mapq_filtered_full_data_tuple.combine(chromosomes))

		// deduplication of 1st and 2nd filtered "samples"
		dedup_bams = DEDUP_CATALOG(possorted_bams_bais, params.dedup_script)
		sorted_dedup_bams_bais = SORT_DEDUP(dedup_bams)

		// fix softclipped region (Alignment fixation) in the dedup file.
		fixed_dedup_bams = FIX_SOFTCLIPPED_REGION_CATALOG(sorted_dedup_bams_bais, params.genome_fasta, params.fix_softclipped_alter_script)

		// sort fixed deduplicated bam in sample
		fixed_dedup_sorted_bams_bais = SORT_FIXED_DEDUP(fixed_dedup_bams)

		//////////////////////////// get_polyA ////////////////////////////

		// polyA and non polyA reads within certain chromosome (for all samples)
		polyA_sample_only = GET_POLYA_SAMPLE_ONLY(fixed_dedup_sorted_bams_bais, params.out_polyA, params.out_non_polyA, params.genome_fasta, params.getPolyA_catalog_script)		

		// sort polyA reads within certain chromosome (for all samples)
		polyA_sorted_bams_bais_samples = SORT_POLYA_SAMPLE_ONLY(polyA_sample_only)

		// groupTuple allows to group file paths by samples.
		// by: 2 means group by 3rd element of the tuple (i.e. samples)
		// channel looks like: [[bam1, bam2...... bam6.....bamY], [bai1, bai2...... bai6.....baiY], sample1]
		// dont need to use .collect()
		// sorted_dedup_bams_bais.groupTuple(by: 2, sort: true, size: num_chrom).view()
		polyA_sorted_full_bams_bais_sample_only = MERGE_POLYA_CATALOG(polyA_sorted_bams_bais_samples.groupTuple(by: 2, sort: true, size: num_chrom), params.out_polyA)

		polyA_sorted_full_bams_bais_counts_samples = GET_COUNTS_CATALOG(polyA_sorted_full_bams_bais_sample_only, params.counts_polyA, params.polyA_type, params.get_count_script)

		polyA_sorted_bams_bais_counts_samples = SPLIT_ALL_POLYA(polyA_sorted_full_bams_bais_counts_samples.combine(chromosomes), params.out_polyA_partial)

		//////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////// get_unique_cleavage_site ////////////////////////////////
		polyA_unique_cs_beds_organs_chrom = GET_POLYA_UNIQUE_CLEAVAGE_SITES_CATALOG(polyA_sorted_bams_bais_counts_samples, params.polyA_unique_cs_bed_sample, 1, params.get_polyA_unique_cs_script)

		direction = Channel
			.from(['+', '-'])

		// combine bed, chromosome tuple with direction and then split bed file by direction.
		(polyA_unique_cs_beds_chrom_direction, polyA_unique_cs_beds_chr_dir_organs) = SPLIT_BY_DIRECTION(polyA_unique_cs_beds_organs_chrom.combine(direction), params.split_by_direction_script)

		// all samples bed files are grouped by its chromosome and direction
		// sort removed in groupTuple for performance. 
		// if you dont specify the size, it will wait until all inputs are recieved.
		all_samples_polyA_unique_cs_beds_chrom_strand = GROUPBY_BED_CATALOG(polyA_unique_cs_beds_chrom_direction.groupTuple(by:[1,2]), params.grouped_by_out, n_second_filtered_samples, params.grouped_by_script)	

		/////////////////////////////////////////////////////////////////////////
		//////////////////////////// Clustering ////////////////////////////////

		(all_polyA_clustered_beds, all_polyA_modified_unique_cs_beds) = PERFORM_CLUSTERING_CATALOG(all_samples_polyA_unique_cs_beds_chrom_strand, params.all_samples_cluster_out, params.modified_unique_cs, params.cluster_pas_alter_script)

		/////////////////////////////////////////////////////////////////////////
		//////////////////////////// Classification of PAS /////////////////////
		our_introns = GET_INTRONIC_BED(our_exons, params.introns_out, params.get_intronic_script)
		
		// Only keep the columns relevant for bedtools intersect
		modified_pas = CHANGE_BED(all_polyA_clustered_beds, params.change_to_bed_script)

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
		
		/////////////////////////////////////////////////////////////////////////
		//////////////////////// Assign gene name to PAS ///////////////////////

		pas_w_geneid = GENE_ID(additional_col_all_pas_merged, our_genes, params.gene_id_pas_out, params.gene_id_alter_nextflow_script)

		/////////////////////////////////////////////////////////////////////////
		//////////////////////////// Group PAS by organs ///////////////////////

		// Grouping with same chromosome, direction and organ
		// to generate organ specific PAS and organ specific modified_unique_cs_beds for a particular chromosome and direction.

		// get() gets each element of the tuple/list. 0-based index
		// tokenize('-') means split string by '-' so that you can get organ out of scinpas-organ by again using get(1)
		// This will enable group bed files by organ.
		// tuple is now (bed, chr, dir, organ). After this you group by chr, dir and organ (collect all beds for same chr, same dir, same organ)
		polyA_unique_cs_beds_chr_dir_organs_improved = polyA_unique_cs_beds_chr_dir_organs.map{it->[it.get(0), it.get(1), it.get(2), it.get(3).tokenize('-').get(1)]}

		// left join with original sample unique cs (which has sample and organ info) and modified unique cs (which has cluster info)
		// you consider all_polyA_modified_unique_cs_beds by one pair of chr and direction at a time but send all them once for simplicity
		organ_specific_cs_filtered = LEFT_JOIN_CATALOG(polyA_unique_cs_beds_chr_dir_organs_improved.groupTuple(by:[1,2,3]), all_polyA_modified_unique_cs_beds.collect(), params.pas_split_by_organ_script)

		/////////////////////////////////////////////////////////////////////////
		/////////////////////////// Organ specific score ///////////////////////

		organ_specific_scores = GET_ORGAN_SCORE(organ_specific_cs_filtered.groupTuple(by:1), params.organ_score, params.organ_score_script)
		pas_w_gene_and_organ_score = MERGE_ORGAN_SCORE_TO_PAS(organ_specific_scores.collect(), pas_w_geneid, params.gene_id_pas_out, params.merge_pas_organ_score_script)
		
		/////////////////////////////////////////////////////////////////
		/////////////////////////// Motif Search ///////////////////////
		
		pas_w_gene_w_organ_score_w_motif = RCS_MOTIF_CHECK(pas_w_gene_and_organ_score, params.gene_id_pas_out, params.rcs_motif_check_out, params.genome_fasta, params.rcs_motif_check_script)

		// worm does not have organs

		///////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////// Convert necessary output to bed.gz ///////////////////////	

		gzip_polyA_unique_cs_beds_organs_chrom = CONVERT_GZIP_UNIQUE_CS(polyA_unique_cs_beds_organs_chrom)
		
		gzip_all_samples_polyA_unique_cs_beds_chrom_strand = CONVERT_GZIP_ALL_SAMPLES_UNIQUE_CS(all_samples_polyA_unique_cs_beds_chrom_strand)
		
		gzip_pas_w_geneid = CONVERT_GZIP_CATALOG_GENEID(pas_w_geneid)

		gzip_pas_w_gene_w_organ_score_w_motif = CONVERT_GZIP_CATALOG_GENEID_ORGAN_MOTIF(pas_w_gene_w_organ_score_w_motif)
		
	}


}