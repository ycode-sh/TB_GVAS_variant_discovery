
if(params.in_data_type == "pe_illumina_reads"){
    params.pe_reads = params.seq_path + "sample_*_{1,2}.fastq"
    params.pe_reads_channel = Channel.fromFilePairs(params.pe_reads)
    params.se_reads_channel = Channel.value("no_se_reads")
    
} else if(params.in_data_type == "se_illumina_reads" | params.in_data_type == "minion_ont_reads"){
    params.se_reads = params.seq_path + "sample_*.fastq"
    params.pe_reads_channel = Channel.value("no_pe_reads")
    params.se_reads_channel = Channel.fromPath(params.se_reads)
    
}


if(params.seq_path == "test_data_path"){
     exit 1, "No fastq files specified -- aborting"
}
if(params.adp_path == "nextera_pe_path"){
    params.adp_path_int = params.adp_path_d
}else {
    params.adp_path_int = params.adp_path
}
if(params.ref_seq_path == "h37rv.fasta"){
    params.ref_seq_path_int = params.ref_seq_path_d
} else {
    params.ref_seq_path_int = params.ref_seq_path
}
if(params.out_dir == "current_working_directory"){
    params.out_dir_int = params.out_dir_d
} else {
    params.out_dir_int = params.out_dir
}

include {trim_fastq; emit_sam; coordsort_sam; bcf_call; gatk_call; gatk_call_gvcf; prepare_manifest; adjudicate_var; minos_joint_genotyping; gatk_joint_genotyping; isec} from "./bin/module_script.nf"

workflow call_genotype_var {
    take:
        trim_sample_script_wf
        emit_sam_script_wf
        coordsort_sam_script_wf
        bcf_call_script_wf
        gatk_call_script_wf
        var_adj_script_wf
        gatk_joint_geno_script_wf
        minos_jg_script_wf
        regenotype_nf_script_wf
        isec_script_wf
        se_reads_wf
        pe_reads_wf
        fastafile
        adapter_seq
        input_read_type
        minos_config_wf
        minos_profile_wf
        min_BQ_wf
        min_MQ_wf
        ploidy_wf
        tandem_qual_wf
        indel_bias_wf
        gvcf_choice
        min_frs_wf
        min_gcp_wf

    main:
        trim_fastq(trim_sample_script_wf, pe_reads_wf, adapter_seq, input_read_type, se_reads_wf)
        emit_sam(emit_sam_script_wf, trim_fastq.out, fastafile, input_read_type)
        // Sort sam and convert to bam 
        coordsort_sam(coordsort_sam_script_wf, emit_sam.out)
        flattened_bam = coordsort_sam.out | flatten | collect
        // Call variants
        bcf_call(bcf_call_script_wf, emit_sam.out, fastafile, min_BQ_wf, min_MQ_wf, ploidy_wf, tandem_qual_wf, indel_bias_wf)
        flattened_bt_vcfs = bcf_call.out | flatten | collect
        if (params.in_data_type != "minion_ont_reads"){
            gatk_call(gatk_call_script_wf, coordsort_sam.out, fastafile, min_BQ_wf, ploidy_wf)
            // Adjudicate and Joint_genotype illumina pe or se calls
            if (params.sample_count == "gvcf"){
                adjudicate_var(bcf_call.out, gatk_call.out, trim_fastq.out, var_adj_script_wf, input_read_type, fastafile, min_frs_wf, min_gcp_wf)
                gatk_call_gvcf(gatk_call_script_wf, coordsort_sam.out, fastafile, min_BQ_wf, ploidy_wf, gvcf_choice)
                flattened_gvcfs = gatk_call_gvcf.out | flatten | collect 
                prepare_manifest(coordsort_sam.out, bcf_call.out)
                minos_joint_genotyping(minos_jg_script_wf, minos_config_wf, minos_profile_wf, regenotype_nf_script_wf, fastafile, prepare_manifest.out, flattened_bam, flattened_bt_vcfs)
                gatk_joint_genotyping(gatk_joint_geno_script_wf, flattened_gvcfs, fastafile)
                isec(isec_script_wf, minos_joint_genotyping.out, gatk_joint_genotyping.out)
                    
            } else if (params.sample_count == "no_gvcf"){
                adjudicate_var(bcf_call.out, gatk_call.out, trim_fastq.out, var_adj_script_wf, input_read_type, fastafile, min_frs_wf, min_gcp_wf)
            }
        
        } else if(params.in_data_type == "minion_ont_reads"){
            // Joint_genotype minioin calls
            if (params.sample_count == "gvcf"){
                prepare_manifest(coordsort_sam.out, bcf_call.out)
                minos_joint_genotyping(minos_jg_script_wf, minos_config_wf, minos_profile_wf, regenotype_nf_script_wf, fastafile, prepare_manifest.out, flattened_bam, flattened_bt_vcfs)
            }
        }
       
}


workflow {
    call_genotype_var(Channel.value(params.trim_sample_script), Channel.value(params.emit_sam_script), Channel.value(params.coordsort_sam_script), Channel.value(params.bcf_call_script), 
    Channel.value(params.gatk_call_script), Channel.value(params.var_adj_script), Channel.value(params.gatk_joint_geno_script), Channel.value(params.minos_jg_script), Channel.value(params.regenotype_nf_script), 
    Channel.value(params.isec_script_wf), params.se_reads_channel, params.pe_reads_channel, Channel.value(params.ref_seq_path_int), Channel.value(params.adp_path_int),
    Channel.value(params.in_data_type), Channel.value(params.minos_config), Channel.value(params.sample_size_profile), Channel.value(params.min_BQ), Channel.value(params.min_MQ),
    Channel.value(params.ploidy), Channel.value(params.tandem_qual) , Channel.value(params.indel_bias), Channel.value(params.sample_count), Channel.value(params.frs), Channel.value(params.gcp))    
}
