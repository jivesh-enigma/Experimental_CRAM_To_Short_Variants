version 1.0
import "https://raw.githubusercontent.com/broadinstitute/warp/develop/tasks/broad/GermlineVariantDiscovery.wdl" as Calling
import "https://raw.githubusercontent.com/broadinstitute/warp/develop/tasks/broad/Utilities.wdl" as Utils
import "https://raw.githubusercontent.com/broadinstitute/warp/WholeGenomeGermlineSingleSample_v3.1.21/pipelines/broad/dna_seq/germline/variant_calling/VariantCalling.wdl" as DragenCaller
import "https://raw.githubusercontent.com/jivesh-enigma/Experimental_CRAM_To_Short_Variants/main/MitochondriaPipeline.wdl" as MitochondrialPipeline

workflow Short_Variant_Pipeline {
    input {
        File ref_fasta
        File DoC_interval_list
        File ref_fasta_index
        File ref_fasta_dict
        File inputCram
        File inputCramIndex
        File master_gene_list
        File GOF_gene_list
        String samplename
        String clinical_bucket_path
        String TestType
        File geneList
        Int runtime_disk
        File WGS_interval_list
        File WES_interval_list
        Int diskSpace
        Int resource_log_interval
        Int runtime_cpus
        String runtime_docker
        Int runtime_preemptible
        Int memoryGb
        Int minBaseQuality
        Int minMappingQuality
        Int preemptible
        Int GATK_diskGb
        Int GATK_memoryGb
        Int RAM
        Int HDD
        Int boot_disk_gb
        File chain
        File Clinvar_genes
        Int cpu_cores
        File GCD_genes
        File OMIM_genes
        Int output_disk_gb
        Int ram_gb
        File Rscript_file
        Int bootDiskSizeGb_VEP
        Int cpu_VEP
        Int diskGb_VEP
        Int fork
        Int memoryGb_VEP
        Int nearestGeneDistance

        # # Set DRAGEN-related arguments according to the "functional equivalence" mode
        # Int dragen_scatter_count = 10
        # Boolean run_dragen_mode_variant_calling_ = true
        # Boolean use_spanning_event_genotyping_ = false
        # Boolean unmap_contaminant_reads_ = false
        # Boolean perform_bqsr_ = false
        # Boolean use_bwa_mem_ = false
        # Boolean use_gatk3_haplotype_caller_ = false
        # Boolean use_dragen_hard_filtering_ = true
        # File dbSNP_vcf
        # File dbSNP_vcf_index
        # File dragen_bam
        # File dragen_bai
    }

    call depthOfCov {
        input: 
            refFasta=ref_fasta,
            intervalList=DoC_interval_list,
            refFastaIndex=ref_fasta_index,
            refFastaDict=ref_fasta_dict,
            inputCram=inputCram,
            inputCramIndex=inputCramIndex,
            sampleName=samplename,
            geneList=geneList,
            memoryGb=memoryGb,
            minBaseQuality=minBaseQuality,
            minMappingQuality=minMappingQuality,
            preemptible=preemptible
            
    }

    # call depthOfCov as depthOfCov_chrM {
    #     input:
    #         refFasta=ref_fasta,
    #         intervalList="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/twist_chrM_panel/Twist_MitoPanel_chrM_all_hg38_target.bed",
    #         refFastaIndex=ref_fasta_index,
    #         refFastaDict=ref_fasta_dict,
    #         inputCram=inputCram,
    #         inputCramIndex=inputCramIndex,
    #         sampleName=samplename,
    #         geneList=geneList,
    #         memoryGb=memoryGb,
    #         minBaseQuality=minBaseQuality,
    #         minMappingQuality=minMappingQuality,
    #         preemptible=preemptible
    # }
    
    call ChooseBed {
        input:
            TestType=TestType, 
            WGS_interval_list=WGS_interval_list,
            WES_interval_list=WES_interval_list

    }

    # Break the calling interval_list into sub-intervals
    # Perform variant calling on the sub-intervals, and then gather the results
    call Utils.ScatterIntervalList as ScatterIntervalList {
        input:
            interval_list = ChooseBed.model_interval_list,
            scatter_count = 10,
            break_bands_at_multiples_of = 100000
    }

    # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
    # If we take the number we are scattering by and reduce by 20 we will have enough disk space
    # to account for the fact that the data is quite uneven across the shards.
    Int potential_hc_divisor = ScatterIntervalList.interval_count - 20
    Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

    # Call variants in parallel over WGS/WES calling intervals
    scatter (scattered_interval_list in ScatterIntervalList.out) {

        call interval_list_to_bed {
            input:
                interval_list=scattered_interval_list
        }

        call deep_variant {
            input:
                sample=samplename,
                capture_bed=interval_list_to_bed.bed,
                Cram=inputCram,
                crai=inputCramIndex,
                reference_fasta=ref_fasta,
                reference_fasta_fai=ref_fasta_index,
                model_type=TestType,
                resource_log_interval=resource_log_interval,
                runtime_cpus=runtime_cpus,
                runtime_docker=runtime_docker,
                runtime_preemptible=runtime_preemptible,
                dv_scatter = hc_divisor
        }

        call bgzip {
            input:
                sample=samplename,
                uncompressed_vcf=deep_variant.vcf,
                runtime_disk=runtime_disk
        }

        File vcfs_to_merge = bgzip.filtered_vcf
        File vcf_indices_to_merge = bgzip.filtered_vcf_index
        File gvcfs_to_merge = deep_variant.gvcf
        File gvcf_indices_to_merge = deep_variant.gvcf_index
        
    }

    # Combine by-interval (g)VCFs into a single sample (g)VCF file
    
    call Calling.MergeVCFs as DV_MergeVCFs {
        input:
            input_vcfs = vcfs_to_merge,
            input_vcfs_indexes = vcf_indices_to_merge,
            output_vcf_name = samplename + ".vcf.gz",
            preemptible_tries = runtime_preemptible
    }

    call Calling.MergeVCFs as DV_MergeGVCFs {
        input:
            input_vcfs_indexes = gvcf_indices_to_merge,
            input_vcfs = gvcfs_to_merge,
            output_vcf_name = samplename + ".gvcf.gz",
            preemptible_tries = runtime_preemptible
    }

    # call DragenCaller.VariantCalling as DragenVCF {
    #     input:
    #         run_dragen_mode_variant_calling = run_dragen_mode_variant_calling_,
    #         use_spanning_event_genotyping = use_spanning_event_genotyping_,
    #         calling_interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list",
    #         evaluation_interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_evaluation_regions.hg38.interval_list",
    #         haplotype_scatter_count = dragen_scatter_count,
    #         break_bands_at_multiples_of = 100000,
    #         # contamination = UnmappedBamToAlignedBam.contamination,
    #         input_bam = dragen_bam,
    #         input_bam_index = dragen_bai,
    #         ref_fasta = ref_fasta,
    #         ref_fasta_index = ref_fasta_index,
    #         ref_dict = ref_fasta_dict,
    #         ref_str = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.str",
    #         dbsnp_vcf = dbSNP_vcf,
    #         dbsnp_vcf_index = dbSNP_vcf_index,
    #         base_file_name = samplename + ".DRAGEN",
    #         final_vcf_base_name = samplename + ".DRAGEN",
    #         agg_preemptible_tries = 5,
    #         use_gatk3_haplotype_caller = use_gatk3_haplotype_caller_,
    #         use_dragen_hard_filtering = use_dragen_hard_filtering_,
    #         make_gvcf = false
    # }

    # call bgzip as Dragen_bgzip {
    #     input:
    #         sample=samplename + ".DRAGEN",
    #         uncompressed_vcf=DragenVCF.output_vcf,
    #         runtime_disk=runtime_disk
    # }

    # call MitochondrialPipeline.MitochondriaPipeline as chrM_calling {
    #     input:
    #         wgs_aligned_input_bam_or_cram = inputCram,
    #         wgs_aligned_input_bam_or_cram_index = inputCramIndex,

    #         ref_fasta = ref_fasta,
    #         ref_fasta_index = ref_fasta_index,
    #         ref_dict = ref_fasta_dict,

    #         mt_dict = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.dict",
    #         mt_fasta = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta",
    #         mt_fasta_index = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.fai",
    #         mt_amb = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.amb",
    #         mt_ann = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.ann",
    #         mt_bwt = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.bwt",
    #         mt_pac = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.pac",
    #         mt_sa = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.sa",
    #         blacklisted_sites = "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed",
    #         blacklisted_sites_index = "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed.idx",

    #         mt_shifted_dict = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict",
    #         mt_shifted_fasta = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta",
    #         mt_shifted_fasta_index = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai",
    #         mt_shifted_amb = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb",
    #         mt_shifted_ann = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann",
    #         mt_shifted_bwt = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt",
    #         mt_shifted_pac = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac",
    #         mt_shifted_sa = "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa",

    #         shift_back_chain = "gs://gcp-public-data--broad-references/hg38/v0/chrM/ShiftBack.chain",

    #         control_region_shifted_reference_interval_list = "gs://gcp-public-data--broad-references/hg38/v0/chrM/control_region_shifted.chrM.interval_list",
    #         non_control_region_interval_list = "gs://gcp-public-data--broad-references/hg38/v0/chrM/non_control_region.chrM.interval_list"
    # }

    # call bgzip as chrM_bgzip {
    #     input:
    #         sample=samplename + ".chrM",
    #         uncompressed_vcf=chrM_calling.split_vcf,
    #         runtime_disk=runtime_disk
    # }

    # call addVAF {
    #     input:
    #         sample=samplename,
    #         dragenVCF = Dragen_bgzip.filtered_vcf,
    #         chrmVCF=chrM_bgzip.filtered_vcf,
    #         runtime_disk = runtime_disk
    # }

    # call gatkVCF {
    #     input:
    #         sample_id=samplename,
    #         dragenVCF=addVAF.dragen_vcf,
    #         chrmVCF=addVAF.chrM_vcf,
    #         runtime_disk = runtime_disk
    # }

    # call mergeVCF {
    #     input:
    #         deepvcf = DV_MergeVCFs.output_vcf,
    #         gatkvcf = gatkVCF.gatk_vcf,
    #         gatkvcf_index = gatkVCF.gatk_vcf_index,
    #         sample_id = samplename,
    #         runtime_disk = runtime_disk
    # }

    call variantcount_vcf {
        input:
            vcf = DV_MergeVCFs.output_vcf,
            sampleId=samplename,
            HDD=HDD,
            RAM=RAM

    }   

    
    call UnZip { 
        input:
            vcfFileGz = DV_MergeVCFs.output_vcf,
            sampleId = samplename,
            RAM=RAM,
            HDD=HDD
    }
    
    call VTRecal {
        input:

            vcfFile=UnZip.vcfFile,
            refFasta=ref_fasta,
            refFastaDict=ref_fasta_dict,
            refFastaIdx=ref_fasta_index,
            sampleId=samplename, 
            RAM=RAM,
            HDD=HDD
    }
    
    call ScatterVCF {
        input:
            samplename = samplename,
            vcf = VTRecal.normalizedVCF,
            vcf_index = VTRecal.normalizedVCF_index,
            memoryGb = memoryGb,
            preemptible = preemptible
    }

    scatter (chrvcf in ScatterVCF.scatteredvcfs) {

        call GATKVariantsToTable {
            input:
                normalizedvcfFileGz=chrvcf,
                refFasta=ref_fasta,
                refFastaFai=ref_fasta_index,
                refFastaDict=ref_fasta_dict,
                samplesetId=basename(chrvcf, ".vcf.gz"),
                GATK_diskGb=GATK_diskGb,
                GATK_memoryGb=GATK_memoryGb
           
        }

        call vep_task {
            input:
        
                refFasta=ref_fasta,
                refFastaFai=ref_fasta_index,
                refFastaDict=ref_fasta_dict,
                samplesetId=basename(chrvcf, ".vcf.gz"),
                normalizedvcfFileGz=chrvcf,
                bootDiskSizeGb_VEP=bootDiskSizeGb_VEP,
                cpu_VEP=cpu_VEP,
                diskGb_VEP=diskGb_VEP,
                fork=fork,
                memoryGb_VEP=memoryGb_VEP,
                nearestGeneDistance=nearestGeneDistance
        }

        call combineOutputFiles {
            input:
                samplesetId=basename(chrvcf, ".vcf.gz"),
                vepOutputFile=vep_task.VEP_Output,
                gatkOutputFile=GATKVariantsToTable.GATK_output,
                diskSpace=diskSpace
            
        }

        File vep_genotypes_to_merge = combineOutputFiles.vepannotated_vcf
    }

    call Merge_VEP_Genotypes {
        input:
            input_veps = vep_genotypes_to_merge,
            samplename = samplename,
            memoryGb = memoryGb_VEP,
            preemptible = preemptible
    }

    call VariantFilter {
        input: 
            input_vcf=Merge_VEP_Genotypes.vepannotated_vcf_merged,
            master_gene_list=master_gene_list,
            GOF_gene_list=GOF_gene_list,
            samplename=samplename,
            boot_disk_gb=boot_disk_gb,
            chain=chain,
            Clinvar_genes=Clinvar_genes,
            cpu_cores=cpu_cores,
            GCD_genes=GCD_genes,
            OMIM_genes=OMIM_genes,
            output_disk_gb=output_disk_gb,
            ram_gb=ram_gb,
            Rscript_file=Rscript_file
    }

    call IGV_Snapshots {
        input: 
            inputCram=inputCram,
            inputCramIndex=inputCramIndex,
            path_var_HQ_IGV_bed=VariantFilter.path_var_HQ_IGV_bed,
            path_var_HQ_non_clinical_IGV_bed=VariantFilter.path_var_HQ_non_clinical_IGV_bed,
            path_var_LQ_IGV_bed=VariantFilter.path_var_LQ_IGV_bed,
            sampleID=samplename,
            memoryGb=memoryGb
    }

    call data_transfer_clinical{
        input:
            clinical_bucket_path=clinical_bucket_path,
            DOC_sampleSummary=depthOfCov.sampleSummary,
            DOC_sampleCumulativeCoverageProportions=depthOfCov.sampleCumulativeCoverageProportions,
            DOC_sampleGeneSummary=depthOfCov.sampleGeneSummary,
            path_var_HQ=path_var_HQ,
            path_var_HQ_non_clinical=path_var_HQ_non_clinical,
            path_var_LQ=path_var_LQ,
            total_variant_count=total_variant_count,
            variants_PGx=variants_PGx,
            variants_HQ_IGV_snapshots=variants_HQ_IGV_snapshots,
            variants_HQ_non_clinical_IGV_snapshots=variants_HQ_non_clinical_IGV_snapshots,
            variants_LQ_IGV_snapshots=variants_LQ_IGV_snapshots,
            sampleID=samplename
    }

  # Outputs that will be retained when execution is complete  
  output {
        #DOC:
        File sampleGeneSummary = depthOfCov.sampleGeneSummary
        File sampleSummary = depthOfCov.sampleSummary
        Float sampleMeanCoverage = depthOfCov.sampleMeanCoverage
        Float sample10XCoverage = depthOfCov.sample10XCoverage
        File sampleStatistics = depthOfCov.sampleStatistics
        File sampleIntervalSummary = depthOfCov.sampleIntervalSummary
        File sampleIntervalStatistics = depthOfCov.sampleIntervalStatistics
        File sampleCumulativeCoverageProportions = depthOfCov.sampleCumulativeCoverageProportions
        File sampleCumulativeCoverageCounts = depthOfCov.sampleCumulativeCoverageCounts
        #DOC chrM:
        # File chrM_sampleGeneSummary = depthOfCov_chrM.sampleGeneSummary
        # File chrM_sampleSummary = depthOfCov_chrM.sampleSummary
        # Float chrM_sampleMeanCoverage = depthOfCov_chrM.sampleMeanCoverage
        # Float chrM_sample10XCoverage = depthOfCov_chrM.sample10XCoverage
        # File chrM_sampleStatistics = depthOfCov_chrM.sampleStatistics
        # File chrM_sampleIntervalSummary = depthOfCov_chrM.sampleIntervalSummary
        # File chrM_sampleIntervalStatistics = depthOfCov_chrM.sampleIntervalStatistics
        # File chrM_sampleCumulativeCoverageProportions = depthOfCov_chrM.sampleCumulativeCoverageProportions
        # File chrM_sampleCumulativeCoverageCounts = depthOfCov_chrM.sampleCumulativeCoverageCounts
        #DV:
        File DV_gvcf = DV_MergeGVCFs.output_vcf
        File DV_gvcf_index = DV_MergeGVCFs.output_vcf_index
        # File DV_resource_log = deep_variant.resource_log
        # File DV_stats_report = deep_variant.stats_report
        File DV_filtered_vcf = DV_MergeVCFs.output_vcf
        File DV_filtered_vcf_index = DV_MergeVCFs.output_vcf_index
        #DRAGEN:
        # File dragen_vcf_summary_metrics = DragenVCF.vcf_summary_metrics
        # File dragen_vcf_detail_metrics = DragenVCF.vcf_detail_metrics
        # File dragen_output_vcf = DragenVCF.output_vcf
        # File dragen_output_vcf_index = DragenVCF.output_vcf_index
        #Mitochondrial:
        # File chrM_subset_bam = chrM_calling.subset_bam
        # File chrM_subset_bai = chrM_calling.subset_bai
        # File chrM_mt_aligned_bam = chrM_calling.mt_aligned_bam
        # File chrM_mt_aligned_bai = chrM_calling.mt_aligned_bai
        # File chrM_out_vcf = chrM_calling.out_vcf
        # File chrM_out_vcf_index = chrM_calling.out_vcf_index
        # File chrM_split_vcf = chrM_calling.split_vcf
        # File chrM_split_vcf_index = chrM_calling.split_vcf_index
        # File chrM_input_vcf_for_haplochecker = chrM_calling.input_vcf_for_haplochecker
        # File chrM_duplicate_metrics = chrM_calling.duplicate_metrics
        # File chrM_coverage_metrics = chrM_calling.coverage_metrics
        # File chrM_theoretical_sensitivity_metrics = chrM_calling.theoretical_sensitivity_metrics
        # File chrM_contamination_metrics = chrM_calling.contamination_metrics
        # File chrM_base_level_coverage_metrics = chrM_calling.base_level_coverage_metrics
        # Int chrM_mean_coverage = chrM_calling.mean_coverage
        # Float chrM_median_coverage = chrM_calling.median_coverage
        # String chrM_major_haplogroup = chrM_calling.major_haplogroup
        # Float chrM_contamination = chrM_calling.contamination
        # # GATK VCF:
        # File gatk_vcf = gatkVCF.gatk_vcf
        # File gatk_vcf_index = gatkVCF.gatk_vcf_index
        # # Merged VCF:
        # File merged_vcf = mergeVCF.merged_vcf
        # File merged_vcf_index = mergeVCF.merged_vcf_index
        # Variant Count
        String variantcount = variantcount_vcf.variantcount
        # VTRecal Normalized VCF
        File normalizedVCF = VTRecal.normalizedVCF
        File normalizedVCF_index = VTRecal.normalizedVCF_index
        # VEP:
        File vepannotated_vcf= Merge_VEP_Genotypes.vepannotated_vcf_merged
        #variant_filtering:
        File path_var_HQ = VariantFilter.path_var_HQ
        File path_var_HQ_non_clinical = VariantFilter.path_var_HQ_non_clinical        
        File path_var_LQ = VariantFilter.path_var_LQ
        File total_variant_count = VariantFilter.total_variant_count
        File variants_PGx = VariantFilter.variants_PGx
        #IGV:
        File variants_HQ_IGV_snapshots = IGV_Snapshots.variants_HQ_IGV_snapshots
        File variants_HQ_non_clinical_IGV_snapshots = IGV_Snapshots.variants_HQ_non_clinical_IGV_snapshots
        File variants_LQ_IGV_snapshots = IGV_Snapshots.variants_LQ_IGV_snapshots
  }
    
}





task depthOfCov {
    input {
    File inputCram
    File inputCramIndex
    Int minBaseQuality
    Int minMappingQuality
    String sampleName
    File intervalList
    File geneList
    File refFasta
    File refFastaDict
    File refFastaIndex
    Int memoryGb
    Int preemptible
    #parameters
    Int addtional_disk_space_gb = 10
    Int disk_space_gb = ceil(size(inputCram, "GB")  * 2 ) + addtional_disk_space_gb
    }
    command <<<
        ln -s ~{inputCram} Cramfile.cram
        ln -s ~{inputCramIndex} Cramfile.crai

        gatk --java-options "-Xmx16G" DepthOfCoverage -R ~{refFasta} -O ~{sampleName} --omit-depth-output-at-each-base true -pt sample -gene-list ~{geneList} -I Cramfile.cram -L ~{intervalList} --min-base-quality ~{minBaseQuality} --summary-coverage-threshold 10
        
        sed -n "2p" < "~{sampleName}.sample_summary" | cut -d',' -f3 > sample_mean_coverage.txt
        sed -n "2p" < "~{sampleName}.sample_summary" | cut -d',' -f7 > sample_10X_coverage.txt

        mv "~{sampleName}.sample_gene_summary" "~{sampleName}.sample_gene_summary.csv"
        mv "~{sampleName}.sample_summary" "~{sampleName}.sample_summary.csv"
        mv "~{sampleName}.sample_interval_summary" "~{sampleName}.sample_interval_summary.csv"
        mv "~{sampleName}.sample_statistics" "~{sampleName}.sample_statistics.csv"
        mv "~{sampleName}.sample_cumulative_coverage_proportions" "~{sampleName}.sample_cumulative_coverage_proportions.csv"
        mv "~{sampleName}.sample_interval_statistics" "~{sampleName}.sample_interval_statistics.csv"


    >>>

    output {
        File sampleGeneSummary = "~{sampleName}.sample_gene_summary.csv"
        File sampleSummary = "~{sampleName}.sample_summary.csv"
        Float sampleMeanCoverage = read_float("sample_mean_coverage.txt")
        Float sample10XCoverage = read_float("sample_10X_coverage.txt")
        File sampleStatistics = "~{sampleName}.sample_statistics.csv"
        File sampleIntervalSummary = "~{sampleName}.sample_interval_summary.csv"
        File sampleIntervalStatistics = "~{sampleName}.sample_interval_statistics.csv"
        File sampleCumulativeCoverageProportions = "~{sampleName}.sample_cumulative_coverage_proportions.csv"
        File sampleCumulativeCoverageCounts = "~{sampleName}.sample_cumulative_coverage_counts"
    }

    runtime {
        docker: "broadinstitute/gatk:latest"
        memory: "~{memoryGb} GB"
        cpu: "1"
        disks: "local-disk " + disk_space_gb + " HDD"
        preemptible: "~{preemptible}"
    }
}


task ChooseBed {
    input {
  File WES_interval_list
  File WGS_interval_list
  String TestType
}


   command <<<

# Check if the input variable is "WES"
if [ ~{TestType} = "WES" ]; then
    bed=~{WES_interval_list}
else
    bed=~{WGS_interval_list}
fi

    mv "$bed" "model_interval_list.bed"



    >>>


  output {
    File model_interval_list = "model_interval_list.bed"
  }
  
    runtime {
        memory: '2 GB'
        disks: 'local-disk 1 HDD'
        preemptible: 5
        docker: 'quay.io/biocontainers/bedtools:2.28.0--hdf88d34_0'
    }
}


task interval_list_to_bed {
    input {
    File interval_list
    String bed_path = sub(basename(interval_list), 'interval_list', 'bed')
}
    command <<<
    set -xeuo pipefail

    # interval lists have headers that need to be removed and are 1-indexed
    # see also https://www.biostars.org/p/84686/
    grep -v '^@' ~{interval_list} \
    | awk -v OFS='\t' '{print $1, $2 - 1, $3}' \
    | sort -k1,1 -k2,2n -k3,3n \
    | bedtools merge \
    > ~{bed_path}
    >>>

    output {
        File bed = '~{bed_path}'
    }

    runtime {
        memory: '1 GB'
        disks: 'local-disk 1 HDD'
        preemptible: 5
        docker: 'quay.io/biocontainers/bedtools:2.28.0--hdf88d34_0'
    }
}




task deep_variant {
    input {
        String sample
        File Cram
        File crai
        String model_type 
        File reference_fasta
        File reference_fasta_fai
        File capture_bed

        Int runtime_cpus
        String runtime_docker
        Int runtime_memory = ceil(1.1 * runtime_cpus)
        Int runtime_preemptible
        Int resource_log_interval
        Int dv_scatter
    }

    Float ref_size = ceil(size(reference_fasta, "GiB") + size(reference_fasta_fai, "GiB"))
    Int disk_space_gb = ceil(((size(Cram, "GiB") + 30) / dv_scatter) + ref_size) + 20

    command <<<
        # log resource usage for debugging purposes
        function runtimeInfo() {
            echo [$(date)]
            echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
            echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
            echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
        }
        while true;
            do runtimeInfo >> resource.log;
            sleep ~{resource_log_interval};
        done &
        lscpu

        set -xeuo pipefail

        # make symbolic links to ensure Cram and index are in expected structure even after localization
        ln -s ~{crai} reads.crai
        ln -s ~{Cram} reads.cram

        # make symbolic links to ensure reference and index are in expected structure even after localization
        ln -s ~{reference_fasta} reference.fa
        ln -s ~{reference_fasta_fai} reference.fa.fai

        mkdir deepvariant_tmp

        /opt/deepvariant/bin/run_deepvariant \
            --model_type=~{model_type} \
            --ref=reference.fa \
            --reads=reads.cram \
            --regions=~{capture_bed} \
            --intermediate_results_dir=deepvariant_tmp \
            --output_vcf=~{sample}.vcf \
            --output_gvcf=~{sample}.g.vcf.gz \
            --num_shards=~{runtime_cpus}
    >>>

    output {
        File vcf = '~{sample}.vcf'
        File gvcf = '~{sample}.g.vcf.gz'
        File gvcf_index = '~{sample}.g.vcf.gz.tbi'
        File resource_log = 'resource.log'
        File stats_report = '~{sample}.visual_report.html'
    }
    
    runtime {
        memory: '~{runtime_memory} GB'
        cpu: '~{runtime_cpus}'
        disks: "local-disk " + disk_space_gb + " HDD"
        preemptible: '~{runtime_preemptible}'
        docker: '~{runtime_docker}'
    }
}

task bgzip {
    input {
    String sample
    File uncompressed_vcf

    Int runtime_disk
}
    command <<<
    set -xeuo pipefail

    bcftools view -Oz -o ~{sample}.vcf.gz ~{uncompressed_vcf}
    bcftools index --tbi ~{sample}.vcf.gz

    # create version of VCF with only PASSing variants
    bcftools view -Oz -o ~{sample}.filtered_callset.vcf.gz -f PASS ~{uncompressed_vcf}
    bcftools index --tbi ~{sample}.filtered_callset.vcf.gz
    >>>

    output {
        File vcf = '~{sample}.vcf.gz'
        File vcf_index = '~{sample}.vcf.gz.tbi'
        File filtered_vcf = '~{sample}.filtered_callset.vcf.gz'
        File filtered_vcf_index = '~{sample}.filtered_callset.vcf.gz.tbi'
    }

    runtime {
        memory: '1 GB'
        disks: 'local-disk ~{runtime_disk} HDD'
        preemptible: 5
        docker: 'quay.io/biocontainers/bcftools:1.9--ha228f0b_3'
    }
}


task variantcount_vcf {
    input {
    File vcf
    String sampleId
    Int RAM
    Int HDD
}
    command <<<
    
        bcftools query -f 'pos=%POS\n' ~{vcf} -o temp.txt
        cat temp.txt | wc -l > ~{sampleId}.variantcount.txt
        
   >>>

    output {
        String variantcount = read_string("~{sampleId}.variantcount.txt")
        File variantcount_file= '~{sampleId}.variantcount.txt'
    }
    
    runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
        memory: RAM + " GB"
        disks: "local-disk " + HDD + " HDD"
        preemptible: 5
    }
}


task UnZip {
    input {
    File vcfFileGz
    String sampleId
    Int RAM
    Int HDD
}
    command <<<
    # Decompress bgzipped merged VCF file
    echo "bgzip -d -c ~{vcfFileGz} > vcfFile.vcf"
    bgzip -d -c ~{vcfFileGz} > ~{sampleId}.vcf

   >>>

    output {
        File vcfFile="~{sampleId}.vcf"
    }
    runtime {
        docker: "vanallenlab/vt:3.13.2018"
        memory: RAM + "GB"
        disks: "local-disk" + " " + HDD + " " + "HDD"
        preemptible: 5
    }
}


task VTRecal {
    input {
    File vcfFile 
    File refFasta
    File refFastaIdx
    File refFastaDict
    String sampleId
    Int RAM
    Int HDD
}
    command <<<
    # VCF:
        echo "########### decompose VCF"
        /software/vt/./vt decompose -s \
        -o ~{sampleId}.vt1.vcf \
        ~{vcfFile}

        echo "########### normalize VCF using ch38 genome build"
        /software/vt/./vt normalize \
        -r ~{refFasta} \
        -o ~{sampleId}.vt2.vcf \
        ~{sampleId}.vt1.vcf
        
        echo "########### normalizing the spanning alleles (*):"
        sed 's/*/-/g' ~{sampleId}.vt2.vcf > ~{sampleId}.vt2_normalized_spanning_alleles.vcf
        bgzip ~{sampleId}.vt2_normalized_spanning_alleles.vcf
        
        echo "########### creating an index for vcf.gz:"
        tabix -p vcf ~{sampleId}.vt2_normalized_spanning_alleles.vcf.gz 


   >>>

    output {
        File normalizedVCF="~{sampleId}.vt2_normalized_spanning_alleles.vcf.gz"
        File normalizedVCF_index="~{sampleId}.vt2_normalized_spanning_alleles.vcf.gz.tbi"

    }

    runtime {
        docker: "vanallenlab/vt:3.13.2018"
        memory: RAM + " GB"
        disks: "local-disk " + HDD + " HDD"
        preemptible: 5
    }
}


task ScatterVCF {
    input {
        String samplename
        File vcf
        File vcf_index
        Int memoryGb
        Int preemptible
    }

    Int scatterdisk = ceil(size(vcf,"GiB") * 4)

    command <<<
        mkdir chrvcf

        for i in {1..22} X Y
        do
            bcftools view ~{vcf} --regions chr${i} -o chrvcf/~{samplename}.chr${i}.vcf.gz -Oz
        done

        ls chrvcf

    >>>

    output {
        Array[File] scatteredvcfs = glob("chrvcf/*.vcf.gz")
    }

    runtime {
        docker: 'quay.io/biocontainers/bcftools:1.9--ha228f0b_3'
        memory: "~{memoryGb} GB"
        cpu: "1"
        disks: "local-disk " + scatterdisk + " HDD"
        preemptible: "~{preemptible}"
    }
}

task vep_task {
    input {
    File normalizedvcfFileGz
    String samplesetId
    # Customizations
    Int nearestGeneDistance
    # Optimizations
    Int fork

    File refFasta
    File refFastaDict
    File refFastaFai
    Int bootDiskSizeGb_VEP
    Int cpu_VEP
    Int diskGb_VEP
    Int memoryGb_VEP
    
    # Cache files
    File speciesCacheTarGzFile="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/VEP_files/homo_sapiens_merged_vep_109_GRCh38.tar.gz"
    String speciesCacheLabel="homo_sapiens_merged"
    String speciesCacheParameter="--merged"
    
    #dbnSFP
    File dbNSFPData="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/dbNSFPv4.2a/dbNSFPv4.2a_modified_header.gz"
    File dbNSFPDataTbi="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/dbNSFPv4.2a/dbNSFPv4.2a_modified_header.gz.tbi"
    File dbNSFPPlugin="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/dbNSFPv4.2a/dbNSFP.pm"
    File dbNSFPReadme="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/dbNSFPv4.2a/dbNSFP4.1a.readme.txt"
    
    #   MCAP
    #File mcap="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch37/annotation_files/VEP_files/mcap_v1.4_modified.vcf.gz"
    #File mcap_index="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch37/annotation_files/VEP_files/mcap_v1.4_modified.vcf.gz.tbi"

    #   S-CAP
    File scap="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/VEP_files/scap_COMBINED_v1.0_modified.vcf.gz"
    File scap_index="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/VEP_files/scap_COMBINED_v1.0_modified.vcf.gz.tbi"
    
    #gnomAD genome
    #File gnomAD_genome="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch37/annotation_files/VEP_files/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz"
    #File gnomAD_genome_index="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch37/annotation_files/VEP_files/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz.tbi"
    
   
    
     #   CLINVAR
    File clinvar="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/VEP_files/clinvar_20240805.vcf.gz"
    File clinvar_index="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/VEP_files/clinvar_20240805.vcf.gz.tbi"
}
    
    
    command <<<
        # Prepare directories for data
        echo "mkdir /opt/vep/.vep/~{speciesCacheLabel}/"
        mkdir /opt
        mkdir /opt/vep
        mkdir /opt/vep/.vep
        mkdir /opt/vep/.vep/~{speciesCacheLabel}/
        mkdir /opt/vep/.vep/Plugins
        mkdir /opt/vep/.vep/Plugins/data
        
        # Put Plugins in correct folder
        echo "symbolic links..."
       
      
        # dbNSFP
        ln -s ~{dbNSFPPlugin} /opt/vep/.vep/Plugins/dbNSFP.pm
        ln -s ~{dbNSFPData} /opt/vep/.vep/Plugins/data/dbNSFPv4.2a_modified_header.gz
        ln -s ~{dbNSFPDataTbi} /opt/vep/.vep/Plugins/data/dbNSFPv4.2a_modified_header.gz.tbi
        ln -s ~{dbNSFPReadme} /opt/vep/.vep/Plugins/data/dbNSFP4.1a.readme.txt
      
        
        # Uncompress the species cache to the .vep directory
        echo "# Uncompress the species cache to the .vep directory"
        echo "tar xzf ~{speciesCacheTarGzFile} -C ~"
        tar xzf ~{speciesCacheTarGzFile} -C .
        
        echo "ls -lh"
        ls -lh
        
        echo "mv ~{speciesCacheLabel}/* /opt/vep/.vep/~{speciesCacheLabel}/"
        mv ~{speciesCacheLabel}/* /opt/vep/.vep/~{speciesCacheLabel}
        echo "ls -lh /opt/vep/.vep/~{speciesCacheLabel}/109_GRCh38/*"
        ls -lh /opt/vep/.vep/~{speciesCacheLabel}/109_GRCh38/*
    
        # log progress to make sure that the VEP output is being generated
        set -xeuo pipefail
        function runtimeInfo() {            
            echo [$(date)]
            echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
            echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
            echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
        }
        while true;
            do runtimeInfo;
            sleep 30;
        done &
   
        
        
        # Run VEP
        echo "running VEP..."
        vep -v -i ~{normalizedvcfFileGz} -o ~{samplesetId}.vep.txt \
        --tab \
        --offline --cache ~{speciesCacheParameter} --dir /opt/vep/.vep --fasta ~{refFasta} \
        --force_overwrite --stats_text --symbol --everything \
        --regulatory --distance ~{nearestGeneDistance}  \
        --total_length --numbers --domains --pick --variant_class --hgvs --hgvsg --ccds  --fork ~{fork} \
        --custom ~{scap},scap_v1.0,vcf,exact,0,Allele_region,rawscore,sensscore,rawscore_dom,sensscore_dom,rawscore_rec,senscore_rec \
        --plugin dbNSFP,/opt/vep/.vep/Plugins/data/dbNSFPv4.2a_modified_header.gz,hg19_chr,hg19_pos,FATHMM_score,FATHMM_pred,PROVEAN_score,MetaSVM_score,MetaLR_score,MetaLR_pred,MetaRNN_score,MetaRNN_pred,M-CAP_score,M-CAP_pred,REVEL_score,MutPred_score,MVP_score,Aloft_pred,LINSIGHT,CADD_raw,GenoCanyon_score,integrated_fitCons_score,Interpro_domain,gnomAD_genomes_MID_AC,gnomAD_genomes_MID_AN,gnomAD_genomes_MID_AF,gnomAD_genomes_MID_nhomalt \
        --custom ~{clinvar},ClinVar_updated_2024Aug,vcf,exact,0,ID,ALLELEID,CLNDN,CLNDISDB,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNVI,DBVARID 
    
        
        
        echo "ls -lh"
        ls -lh
        
        echo "ls -lh opt/vep/.vep/*"
        ls -lh /opt/vep/.vep/*
        
        echo "Number of VEP variants (grep -v # ~{samplesetId}.vep.txt | wc -l):"
        grep -v "#" ~{samplesetId}.vep.txt | wc -l
        
        # gzip ~{samplesetId}.vep.vcf
    >>>
    
    runtime {
        docker: "ensemblorg/ensembl-vep:release_109.2"    
        bootDiskSizeGb : "~{bootDiskSizeGb_VEP}"
        preemptible    : 5
        cpu            : "~{cpu_VEP}"
        disks          : "local-disk ~{diskGb_VEP} SSD"
        memory         : "~{memoryGb_VEP} GB"
    }

    output {        
        File VEP_Output="~{samplesetId}.vep.txt"
        File VEP_Summary="~{samplesetId}.vep.txt_summary.txt"
    }
}

task GATKVariantsToTable {
    input {   
        File normalizedvcfFileGz
        File refFasta
        File refFastaDict
        File refFastaFai
        String samplesetId
        Int GATK_diskGb
        Int GATK_memoryGb
    }   
        command <<<      
        echo "ls -lh"
        ls -lh
        ls ~{normalizedvcfFileGz}
        
        mv ~{normalizedvcfFileGz} vcfFile.vcf.gz
        
        echo "ls -lh"
        ls -lh
        
        echo "bgzip decompressing vt recal VCF file"
        bgzip --decompress vcfFile.vcf.gz
        
        echo "ls -lh"
        ls -lh 
        
        echo "ls -lh vcfFile.vcf"
        ls -lh vcfFile.vcf
        
        echo "########### Using GATK to extract variants into a table format (GRCh38)"
        java -jar /usr/GenomeAnalysisTK.jar -R ~{refFasta} -T VariantsToTable \
        -V vcfFile.vcf -o ~{samplesetId}.vt2_GATK_annotations.txt \
        -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -GF AD -GF DP  -GF GQ  -GF VAF -GF PL -GF GT --allowMissingData --showFiltered 

        echo '########### Done extracting relevant fields using GATK'

        # count the number of GATK variants:
        echo "########### number of GATK variants: "
        cat ~{samplesetId}.vt2_GATK_annotations.txt | wc -l
        
        # gzip ~{samplesetId}.vt2_GATK_annotations.txt
        
        echo "ls -lh"
        ls -lh
        >>>
        
        output {
            File GATK_output="~{samplesetId}.vt2_GATK_annotations.txt"
        }
        
        runtime {
            docker: "vanallenlab/gatk3.7_with_htslib1.9:1.0"
            memory: "~{GATK_memoryGb} GB"
            cpu: "1"
            disks: "local-disk ~{GATK_diskGb} HDD"
            preemptible: 5
        }
}

task combineOutputFiles {
    input {
        File vepOutputFile
        File gatkOutputFile
        String samplesetId
        Int diskSpace
    }

    command <<<
      cp ~{vepOutputFile} vepOutputFile.txt
      cp ~{gatkOutputFile} gatkOutputFile.txt
    
      # remove the '#' from the first line before parsing:
      echo "########### removing the # sign from the first line of the VEP output file"
      sed  "s/#Uploaded_variation/Uploaded_variation/g" vepOutputFile.txt > ~{samplesetId}_vt2_VEP_temp2.txt

      # remove the excess header part:
      grep "#" ~{samplesetId}_vt2_VEP_temp2.txt > ~{samplesetId}_VEP_annotation_list.txt
      grep -v "##" ~{samplesetId}_vt2_VEP_temp2.txt > ~{samplesetId}_vt2_VEP.txt

      # count the number of VEP variants:
      echo "########### number of VEP variants"
      cat ~{samplesetId}_vt2_VEP.txt | wc -l
      
      echo "########### Combining VEP  output files"
      paste ~{samplesetId}_vt2_VEP.txt gatkOutputFile.txt > ~{samplesetId}_vt2_VEP_Genotypes.txt
    >>>
    
    output {
        File vepannotated_vcf="~{samplesetId}_vt2_VEP_Genotypes.txt"
        File annotationsList="~{samplesetId}_VEP_annotation_list.txt"
    }
    
    runtime {
        docker: "ensemblorg/ensembl-vep:release_109.2"    
        preemptible    : 5
        cpu            : "1"
        disks          : "local-disk ~{diskSpace} HDD"
        memory         : "10 GB"
    }
}

task Merge_VEP_Genotypes {
    input {
        Array[File] input_veps
        String samplename
        Int memoryGb
        Int preemptible
    }

    Int disk_size = ceil(size(input_veps,"GiB") * 4) + 10

    command <<<
        mkdir vep_files
        # counter=1
        
        for file in $(cat ~{write_lines(input_veps)}); do
            mv ${file} vep_files/
        done
        
        ls vep_files

        python3 <<CODE
        import pandas as pd
        import os
        import re

        # Set the directory containing the files
        directory = 'vep_files'

        # Define a function to extract the chromosome number from the filename
        def extract_chr_number(filename):
            # Match the pattern for chromosomes (numeric, X, Y, M)
            match = re.search(r'\.chr(\d+|X|Y|M)_vt2_VEP_Genotypes\.txt', filename)
            if match:
                chr_str = match.group(1)
                if chr_str == 'X':
                    return 23
                elif chr_str == 'Y':
                    return 24
                elif chr_str == 'M':
                    return 25
                else:
                    return int(chr_str)
            # Return infinity for unexpected formats to sort them at the end
            return float('inf')

        # List all files in the directory
        files = [f for f in os.listdir(directory) if f.endswith('_vt2_VEP_Genotypes.txt')]

        # Sort the files by chromosome number
        files.sort(key=extract_chr_number)

        # Initialize an empty list to hold the dataframes
        dataframes = []

        # Loop over the sorted files and read them into pandas dataframes
        for file in files:
            file_path = os.path.join(directory, file)
            df = pd.read_csv(file_path, sep='\t', low_memory=False)
            dataframes.append(df)

        # Concatenate all dataframes
        combined_df = pd.concat(dataframes, ignore_index=True)

        # Save the combined dataframe to a new file
        combined_df.to_csv('~{samplename}_vt2_VEP_Genotypes.txt', sep='\t', index=False)
        CODE

        # Compressing the output
        gzip ~{samplename}_vt2_VEP_Genotypes.txt
    >>>

    output {
        File vepannotated_vcf_merged = '~{samplename}_vt2_VEP_Genotypes.txt.gz'
    }

    runtime {
        docker: 'jiveshenigma/terra-pgs-env:v3'
        memory: "~{memoryGb} GB"
        cpu: "2"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: "~{preemptible}"
    }
}


task VariantFilter {
    input {
    File input_vcf
    File master_gene_list
    File GOF_gene_list
    File OMIM_genes
    File Clinvar_genes
    File GCD_genes
    File chain
    File Rscript_file
    Int cpu_cores
    Int ram_gb
    Int boot_disk_gb
    Int output_disk_gb
    String samplename
}
   command <<<
       Rscript ~{Rscript_file} ~{master_gene_list} ~{GOF_gene_list} ~{OMIM_genes} ~{Clinvar_genes} ~{GCD_genes} ~{input_vcf} ~{samplename} ~{chain}  
   >>>
   
   runtime { 
     docker : "lbwang/rocker-genome"
     bootDiskSizeGb: "~{boot_disk_gb}"
     preemtible : 0
     disks: "local-disk ~{output_disk_gb} HDD"
     cpu: "~{cpu_cores}"
     memory: "~{ram_gb}GB"
    }   
    output {
        File path_var_HQ="~{samplename}.pathogenic_variants_all_high_quality.csv"
        File path_var_HQ_non_clinical="~{samplename}.pathogenic_variants_all_high_quality_non_clinical.csv"        
        File path_var_LQ="~{samplename}.pathogenic_variants_all_low_quality.csv"
        File path_var_HQ_IGV_bed="~{samplename}_pathogenic_variants_all_high_quality_IGV.bed"
        File path_var_HQ_non_clinical_IGV_bed="~{samplename}_pathogenic_variants_all_high_quality_non_clinical_IGV.bed"
        File path_var_LQ_IGV_bed="~{samplename}_pathogenic_variants_all_low_quality_IGV.bed"
        File all_variants_cancer_genes="~{samplename}.all_variants_cancer_genes.csv.gz"
        File all_variants_carrier_genes="~{samplename}.all_variants_carrier_genes.csv.gz"
        File total_variant_count="~{samplename}.total_variant_count.csv"
        File variants_PGx="${samplename}.variants_PGx.csv"
       }  
}



task IGV_Snapshots {
    input {
    File inputCram
    File inputCramIndex
    File path_var_HQ_IGV_bed
    File path_var_HQ_non_clinical_IGV_bed
    File path_var_LQ_IGV_bed
    Int memoryGb
    Int addtional_disk_space_gb = 10
    Int disk_space_gb = ceil(size(inputCram, "GB")  * 2 ) + addtional_disk_space_gb
    Int bootDiskSizeGb = ceil(size(inputCram, "GB")  * 2 ) + addtional_disk_space_gb
    String sampleID
    String genome = "hg38"
}
    command <<<
    ls

## HQ:
    echo "moving input files to correct directory..."
    
    mv ~{path_var_HQ_IGV_bed} /IGV-snapshot-automator/regions.bed

    cd /IGV-snapshot-automator
    
    mv ~{inputCram} /cromwell_root/~{sampleID}.cram
    mv ~{inputCramIndex} /cromwell_root/~{sampleID}.cram.crai

    /bin/sh -c "if [ ! -d 'IGV_Snapshots' ]; then mkdir IGV_Snapshots; fi"
    
    echo "creating batch file..."
    
    python3 make_IGV_snapshots.py -g ~{genome} /cromwell_root/~{sampleID}.cram
    
    xvfb-run --auto-servernum --server-num=1 igv.sh -g ~{genome} -b IGV_Snapshots/IGV_snapshots.bat

    echo "---- finished running igv snapshot automator ----"
    echo "compressing IGV_Snapshots directory for output..."
    
    tar -zcf IGV_Snapshots.tar.gz IGV_Snapshots
    
    echo "renaming & moving IGV_Snapshots directory to cromwell_root..."
    
    mv IGV_Snapshots.tar.gz ~{sampleID}.IGV_Snapshots_HQ.tar.gz
    mv ~{sampleID}.IGV_Snapshots_HQ.tar.gz /cromwell_root
    
    echo "---- completed running IGV_Snapshots ----"




    ## HQ_non_clinical:
    echo "moving input files to correct directory..."
    
    mv ~{path_var_HQ_non_clinical_IGV_bed} /IGV-snapshot-automator/regions.bed

    cd /IGV-snapshot-automator
    
    mv ~{inputCram} /cromwell_root/~{sampleID}.cram
    mv ~{inputCramIndex} /cromwell_root/~{sampleID}.cram.crai

    /bin/sh -c "if [ ! -d 'IGV_Snapshots' ]; then mkdir IGV_Snapshots; fi"
    
    echo "creating batch file..."
    
    python3 make_IGV_snapshots.py -g ~{genome} /cromwell_root/~{sampleID}.cram
    
    xvfb-run --auto-servernum --server-num=1 igv.sh -g ~{genome} -b IGV_Snapshots/IGV_snapshots.bat

    echo "---- finished running igv snapshot automator ----"
    echo "compressing IGV_Snapshots directory for output..."
    
    tar -zcf IGV_Snapshots.tar.gz IGV_Snapshots
    
    echo "renaming & moving IGV_Snapshots directory to cromwell_root..."
    
    mv IGV_Snapshots.tar.gz ~{sampleID}.IGV_Snapshots_HQ_non_clinical.tar.gz
    mv ~{sampleID}.IGV_Snapshots_HQ_non_clinical.tar.gz /cromwell_root
    
    echo "---- completed running IGV_Snapshots ----"



## LQ:
    echo "moving input files to correct directory..."
    
    mv ~{path_var_LQ_IGV_bed} /IGV-snapshot-automator/regions.bed

    cd /IGV-snapshot-automator
    
    mv ~{inputCram} /cromwell_root/~{sampleID}.cram
    mv ~{inputCramIndex} /cromwell_root/~{sampleID}.cram.crai

    /bin/sh -c "if [ ! -d 'IGV_Snapshots' ]; then mkdir IGV_Snapshots; fi"
    
    echo "creating batch file..."
    
    python3 make_IGV_snapshots.py -g ~{genome} /cromwell_root/~{sampleID}.cram
    
    xvfb-run --auto-servernum --server-num=1 igv.sh -g ~{genome} -b IGV_Snapshots/IGV_snapshots.bat

    echo "---- finished running igv snapshot automator ----"
    echo "compressing IGV_Snapshots directory for output..."
    
    tar -zcf IGV_Snapshots.tar.gz IGV_Snapshots
    
    echo "renaming & moving IGV_Snapshots directory to cromwell_root..."
    
    mv IGV_Snapshots.tar.gz ~{sampleID}.IGV_Snapshots_LQ.tar.gz
    mv ~{sampleID}.IGV_Snapshots_LQ.tar.gz /cromwell_root
    
    echo "---- completed running IGV_Snapshots ----"


    >>>

    output {
        File variants_HQ_IGV_snapshots = "~{sampleID}.IGV_Snapshots_HQ.tar.gz"
        File variants_HQ_non_clinical_IGV_snapshots = "~{sampleID}.IGV_Snapshots_HQ_non_clinical.tar.gz"
        File variants_LQ_IGV_snapshots = "~{sampleID}.IGV_Snapshots_LQ.tar.gz"
    }

    runtime {
        docker: "tylerchinskydfci/igv_snapshot:0.1"
        memory: "~{memoryGb} GB"
        cpu: "1"
        disks: "local-disk " + disk_space_gb + " HDD"
    }
}



task data_transfer_clinical {
    input {
        File DOC_sampleSummary
        File DOC_sampleCumulativeCoverageProportions
        File DOC_sampleGeneSummary
        File path_var_HQ
        File path_var_HQ_non_clinical
        File path_var_LQ
        File total_variant_count
        File variants_PGx
        File variants_HQ_IGV_snapshots
        File variants_HQ_non_clinical_IGV_snapshots
        File variants_LQ_IGV_snapshots
        String clinical_bucket_path
        String sampleID
    }

    command <<<
    
        gsutil -m cp ~{DOC_sampleSummary} ~{clinical_bucket_path}/~{sampleID}.sample_summary.csv
        gsutil -m cp ~{DOC_sampleCumulativeCoverageProportions} ~{clinical_bucket_path}/~{sampleID}.sample_cumulative_coverage_proportions.csv
        gsutil -m cp ~{DOC_sampleGeneSummary} ~{clinical_bucket_path}/~{sampleID}.sample_gene_summary.csv
        gsutil -m cp ~{path_var_HQ} ~{clinical_bucket_path}/~{sampleID}.path_var_HQ.csv
        gsutil -m cp ~{path_var_HQ_non_clinical} ~{clinical_bucket_path}/~{sampleID}.path_var_HQ_non_clinical.csv
        gsutil -m cp ~{path_var_LQ} ~{clinical_bucket_path}/~{sampleID}.path_var_LQ.csv
        gsutil -m cp ~{total_variant_count} ~{clinical_bucket_path}/~{sampleID}.total_variant_count.csv
        gsutil -m cp ~{variants_PGx} ~{clinical_bucket_path}/~{sampleID}.variants_PGx.csv
        gsutil -m cp ~{variants_HQ_IGV_snapshots} ~{clinical_bucket_path}/~{sampleID}.variants_HQ_IGV_snapshots.tar.gz
        gsutil -m cp ~{variants_HQ_non_clinical_IGV_snapshots} ~{clinical_bucket_path}/~{sampleID}.variants_HQ_non_clinical_IGV_snapshots.tar.gz
        gsutil -m cp ~{variants_LQ_IGV_snapshots} ~{clinical_bucket_path}/~{sampleID}.variants_LQ_IGV_snapshots.tar.gz

   >>>

    output {

    }
    runtime {
        docker: "google/cloud-sdk"
        memory: "1GB"
        disks: 'local-disk 1 HDD' 
        preemptible: 5
    }
}

task addVAF {
    input {
        String sample
        File dragenVCF
        File chrmVCF
        Int runtime_disk
    }

    command <<<
        # Add VAF to the dragen and mitochondrial vcfs
        bcftools +fill-tags -Oz -o ~{sample}.DRAGEN.vcf.gz ~{dragenVCF} -- -t VAF
        
        bcftools +fill-tags -Oz -o ~{sample}.chrM.vcf.gz ~{chrmVCF} -- -t VAF
    >>>

    output {
        File dragen_vcf = "~{sample}.DRAGEN.vcf.gz"
        File chrM_vcf = "~{sample}.chrM.vcf.gz"
    }

    runtime {
        memory: '1 GB'
        disks: 'local-disk ~{runtime_disk} HDD'
        preemptible: 5
        docker: 'jiveshenigma/htslib-samtools-bcftools:v1'
    }
}

task gatkVCF {
    input {
        String sample_id
        File dragenVCF
        File chrmVCF
        Int runtime_disk
    }

    command <<<
        # Concatenate the vcfs
        bcftools concat ~{dragenVCF} ~{chrmVCF} -Oz -o ~{sample_id}.concat.vcf.gz

        # Create a renaming file to change sample name
        echo "~{sample_id} ~{sample_id}_GATK" > gatk_renaming.txt

        # Change sample name in the concatenated vcf
        bcftools reheader -s gatk_renaming.txt -o ~{sample_id}.gatk.vcf.gz ~{sample_id}.concat.vcf.gz

        # Create index
        bcftools index -t ~{sample_id}.gatk.vcf.gz

    >>>

    output {
        File gatk_vcf = "~{sample_id}.gatk.vcf.gz"
        File gatk_vcf_index = "~{sample_id}.gatk.vcf.gz.tbi"
    }

    runtime {
        memory: '1 GB'
        disks: 'local-disk ~{runtime_disk} HDD'
        preemptible: 5
        docker: 'jiveshenigma/htslib-samtools-bcftools:v1'
    }
}

task mergeVCF {
    input {
        File deepvcf
        File gatkvcf
        File gatkvcf_index
        String sample_id
        Int runtime_disk
    }

    command <<<
        # Create a renaming file to change sample name in DV vcf
        echo "~{sample_id} ~{sample_id}_DV" > dv_renaming.txt

        # Change sample name in the DV vcf
        bcftools reheader -s dv_renaming.txt -o ~{sample_id}.dv.vcf.gz ~{deepvcf}

        # Creating index for new dv vcf
        bcftools index -t ~{sample_id}.dv.vcf.gz

        # Merging the DV vcf and GATK VCF
        bcftools merge ~{sample_id}.dv.vcf.gz ~{gatkvcf} -Oz -o ~{sample_id}.short_variants_merged.vcf.gz

        # index for merged vcf
        bcftools index -t ~{sample_id}.short_variants_merged.vcf.gz

    >>>

    output {
        File merged_vcf = "~{sample_id}.short_variants_merged.vcf.gz"
        File merged_vcf_index = "~{sample_id}.short_variants_merged.vcf.gz.tbi"
    }

    runtime {
        memory: '1 GB'
        disks: 'local-disk ~{runtime_disk} HDD'
        preemptible: 5
        docker: 'jiveshenigma/htslib-samtools-bcftools:v1'
    }
    
}