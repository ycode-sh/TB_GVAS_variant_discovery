
include {snpEff_annotation; custom_annotations; bedtools_clin_assay; bedtools_epi_assay; phylo_tree; drug_res_profiling; snp_typing} from "./bin/module_script.nf"
// Choices include: bio_assays, clin_assays, epi_assays
// Inputs: a vcf file per sample

if(params.ref_seq_path == "h37rv.fasta"){
    params.ref_seq_path_int = params.ref_seq_path_d
} else {
    params.ref_seq_path_int = params.ref_seq_path
}


params.vcf_path = params.vcf_file_path + "/sample_*.vcf"
params.dr_catalogue_files = params.dr_catalogue_file_path + "/sorted*tsv"

workflow variant_transformation {
    take:
        vcf_file_wf
        snpEff_annot_script_wf 
        snpeff_jar_path
        snpeff_config_path
        custom_annotations_script_wf
        gene_annotation_file_wf
        gene_annotation_header_wf
        lineage_annotation_file_wf
        lineage_annotation_header_wf
        amr_region_file_wf
        amr_region_header_wf
        variable_region_annotation_file_wf
        variable_region_annotation_header_wf
        bedtools_script_wf
        lineage_tsv_wf
        assay_type_wf
        generate_consensus_script_wf
        fastafile
        dr_profiling_scipt_wf
        snp_typing_script_wf
        dr_extr_scipt_2_wf
        dr_extr_scipt_3_wf
        dr_catalogue_file_path
        prepped_lineage_file
        variant_caller
        dp_cov
        dr_R_formatting_script_wf
        snp_typing_R_formating_script_wf
  
    main:
        snpEff_annotation(snpEff_annot_script_wf, vcf_file_wf, snpeff_jar_path, snpeff_config_path)
        if (params.assay_type == "bio_assays"){
            custom_annotations(custom_annotations_script_wf, gene_annotation_file_wf, gene_annotation_header_wf, lineage_annotation_file_wf, lineage_annotation_header_wf, amr_region_file_wf, amr_region_header_wf, variable_region_annotation_file_wf, variable_region_annotation_header_wf, snpEff_annotation.out.snpeff_annotated_vcf)
        } else if (params.assay_type == "clin_assays"){
            bedtools_clin_assay(bedtools_script_wf, snpEff_annotation.out.snpeff_annotated_vcf, amr_region_file_wf, lineage_tsv_wf, assay_type_wf)
            collected_preped_vcfs = bedtools_clin_assay.out | collect
            dr_res_catalogue  = dr_catalogue_file_path | collect
            drug_res_profiling(dr_profiling_scipt_wf, collected_preped_vcfs, dr_res_catalogue, variant_caller, dp_cov, dr_extr_scipt_2_wf, dr_extr_scipt_3_wf, dr_R_formatting_script_wf)
            snp_typing(snp_typing_script_wf, collected_preped_vcfs, prepped_lineage_file, variant_caller, dp_cov, dr_extr_scipt_2_wf, dr_extr_scipt_3_wf, snp_typing_R_formating_script_wf)
        } else if (params.assay_type == "epi_assays"){
            bedtools_epi_assay(bedtools_script_wf, vcf_file_wf, variable_region_annotation_file_wf, assay_type_wf)
            phylo_tree(generate_consensus_script_wf, fastafile, bedtools_epi_assay.out)
        } else if (params.assay_type == "all_assays"){
            custom_annotations(custom_annotations_script_wf, gene_annotation_file_wf, gene_annotation_header_wf, lineage_annotation_file_wf, lineage_annotation_header_wf, amr_region_file_wf, amr_region_header_wf, variable_region_annotation_file_wf, variable_region_annotation_header_wf, snpEff_annotation.out.snpeff_annotated_vcf)
            bedtools_clin_assay(bedtools_script_wf, snpEff_annotation.out.snpeff_annotated_vcf, amr_region_file_wf, lineage_tsv_wf, assay_type_wf)
            collected_preped_vcfs = bedtools_clin_assay.out | flatten | collect
            dr_res_catalogue  = dr_catalogue_file_path | collect
            drug_res_profiling(dr_profiling_scipt_wf, collected_preped_vcfs, dr_res_catalogue, variant_caller, dp_cov, dr_extr_scipt_2_wf, dr_extr_scipt_3_wf, dr_R_formatting_script_wf)
            snp_typing(snp_typing_script_wf, collected_preped_vcfs, prepped_lineage_file, variant_caller, dp_cov, dr_extr_scipt_2_wf, dr_extr_scipt_3_wf, snp_typing_R_formating_script_wf)
            bedtools_epi_assay(bedtools_script_wf, vcf_file_wf, variable_region_annotation_file_wf, assay_type_wf)
            phylo_tree(generate_consensus_script_wf, fastafile, bedtools_epi_assay.out)
        }

    emit:
        bedtools_clin_assay.out
}


workflow {
    variant_transformation(Channel.fromPath(params.vcf_path), Channel.value(params.snpEff_annot_script), Channel.value(params.snpeff_jar_path), Channel.value(params.snpeff_config_path),
    Channel.value(params.custom_annotations_script), Channel.value(params.gene_annotation_file), Channel.value(params.gene_annotation_header), Channel.value(params.lineage_annotation_file),
    Channel.value(params.lineage_annotation_header), Channel.value(params.amr_region_file), Channel.value(params.amr_region_header), Channel.value(params.variable_region_annotation_file),
    Channel.value(params.variable_region_annotation_header), Channel.value(params.bedtools_script), Channel.value(params.lineage_tsv), Channel.value(params.assay_type), 
    Channel.value(params.generate_consensus_script), Channel.value(params.ref_seq_path_int), Channel.value(params.dr_profiling_scipt), Channel.value(params.snp_typing_python_script), Channel.value(params.dr_extr_scipt_2), 
    Channel.value(params.dr_extr_scipt_3), Channel.fromPath(params.dr_catalogue_files),
    Channel.value(params.prepped_lineage_file), Channel.value(params.variant_caller), Channel.value(params.reads_parsing_depth), Channel.value(params.dr_R_formatting_script), 
    Channel.value(params.snp_typing_R_formating_script))
}