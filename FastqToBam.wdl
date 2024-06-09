version 1.0

workflow FastqToBam {
  input {
    String sample_name
    File fastq_1
    File fastq_2
    Boolean clean_fastqs = false
    Int n_fastq_shards = 1
    #Boolean cram_out = true

    #String reference_name
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    File? reference_alt
    File reference_sa 
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_amb
    # File dbSNP_vcf
    # File dbSNP_vcf_index
    # Array[File] known_indels_sites_VCF
    # Array[File] known_indels_sites_index

    # File geneList

    # String permanent_bucket_path

    String fastq_pair_docker = "vanallenlab/fastq-pair:latest"
    String fastqsplitter_docker = "vanallenlab/fastqsplitter:latest"
    String gotc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    # String gatk_docker = "broadinstitute/gatk:latest"
    # String gatk_path = "/gatk/gatk"
    # String python_docker = "python:2.7"

    String fastq_suffix = ".fq.gz"
    Int clean_fastq_reads_per_hash_cell = 8
    Int clean_fastq_disk_scaling_factor = 6

    Float clean_fastq_mem_gb = 64
    # Int? sort_and_fix_tags_disk_size
    Float? bwa_mem_gb
    Int? bwa_num_cpu
    Int? bwa_disk_size
    # Int? gather_bam_files_disk_size
  }

  # Clean up FASTQs with fastq-pair
  if ( clean_fastqs ) {
    call CleanFastqs {
      input:
        fastq_1 = fastq_1,
        fastq_2 = fastq_2,
        fastq_suffix = fastq_suffix,
        reads_per_hash_cell = clean_fastq_reads_per_hash_cell,
        disk_scaling_factor = clean_fastq_disk_scaling_factor,
        mem_gb = clean_fastq_mem_gb,
        docker = fastq_pair_docker
    }
  }
  File fq1 = select_first([CleanFastqs.fq_paired_1, fastq_1])
  File fq2 = select_first([CleanFastqs.fq_paired_2, fastq_2])
  Array[File] fq_singles = select_first([CleanFastqs.fq_singles, []])

  # Split paired fastqs to parallelize alignment
  if ( n_fastq_shards > 1 ) {
    call SplitFastq as SplitFq1 {
      input:
        fastq = fq1,
        n_shards = n_fastq_shards,
        docker = fastqsplitter_docker
    }
    call SplitFastq as SplitFq2 {
      input:
        fastq = fq2,
        n_shards = n_fastq_shards,
        docker = fastqsplitter_docker
    }
  }
  Array[File] sharded_fq1 = select_first([SplitFq1.fq_shards, [fq1]])
  Array[File] sharded_fq2 = select_first([SplitFq2.fq_shards, [fq2]])
  Array[Pair[File, File]] sharded_fq_pairs = zip(sharded_fq1, sharded_fq2)

  # Align paired fastqs with BWA
  scatter ( fq_pair in sharded_fq_pairs ) {
    call Bwa as AlignPairs {
      input:
        input_fastqs = [fq_pair.left, fq_pair.right],
        output_basename = sample_name + "." + basename(fq_pair.left, ".fq.gz"),
        ref_fasta = reference_fasta,
        ref_fasta_index = reference_fasta_index,
        ref_dict = reference_dict,
        ref_alt = reference_alt,
        ref_amb = reference_amb,
        ref_ann = reference_ann,
        ref_bwt = reference_bwt,
        ref_pac = reference_pac,
        ref_sa = reference_sa,
        mem_gb = bwa_mem_gb,
        num_cpu = bwa_num_cpu,
        disk_size = bwa_disk_size,
        docker = gotc_docker
    }
  }

  scatter ( fastq in fq_singles ) {
    call Bwa as AlignSingles {
      input:
        input_fastqs = [fastq],
        output_basename = sample_name,
        ref_fasta = reference_fasta,
        ref_fasta_index = reference_fasta_index,
        ref_dict = reference_dict,
        ref_alt = reference_alt,
        ref_amb = reference_amb,
        ref_ann = reference_ann,
        ref_bwt = reference_bwt,
        ref_pac = reference_pac,
        ref_sa = reference_sa,
        mem_gb = bwa_mem_gb,
        num_cpu = bwa_num_cpu,
        disk_size = bwa_disk_size,
        docker = gotc_docker
    }
  }

  output {
    # Merge BAMs, mark duplicates
    Array[File] bams_for_gatk = flatten([AlignPairs.output_bam, AlignSingles.output_bam])
  }
}

# Clean up FASTQs with fastq-pair
task CleanFastqs {
  input {
    File fastq_1
    File fastq_2
    String docker
    String fastq_suffix = ".fq.gz"
    Int reads_per_hash_cell = 3
    Int disk_scaling_factor = 10
    Float mem_gb = 31.5
  }

  String f1_basename = basename(fastq_1, fastq_suffix)
  String f2_basename = basename(fastq_2, fastq_suffix)

  Int disk_gb_base = ceil(3 * size([fastq_1, fastq_2], "GB"))
  Int disk_gb = (disk_scaling_factor * disk_gb_base) + 10

  command <<<
    set -eu -o pipefail

    # Check whether fastqs are gzipped and uncompress if necessary
    if [ $( file ~{fastq_1} | fgrep gzip | wc -l ) -gt 0 ]; then
      zcat ~{fastq_1} > ~{f1_basename}.fastq
    else
      mv ~{fastq_1} > ~{f1_basename}.fastq
    fi
    if [ $( file ~{fastq_2} | fgrep gzip | wc -l ) -gt 0 ]; then
      zcat ~{fastq_1} > ~{f2_basename}.fastq
    else
      mv ~{fastq_1} > ~{f2_basename}.fastq
    fi

    # Get number of reads
    n_reads=$( cat ~{f1_basename}.fastq | wc -l | awk '{ print $1 / 4 }' | cut -f1 -d\. )
    echo -e "\nCounted $n_reads total reads\n"

    # Set hash size to n_reads / reads_per_hash_cell
    hash_size=$( echo $n_reads | awk -v denom=~{reads_per_hash_cell} '{ printf "%.0f\n", $1 / denom }' | cut -f1 -d\. )
    echo -e "\nSetting hash size to $hash_size\n"

    fastq_pair -p -t $hash_size \
      ~{f1_basename}.fastq \
      ~{f2_basename}.fastq
    for out in paired single; do
      gzip -f ~{f1_basename}.fastq.$out.fq
      gzip -f ~{f2_basename}.fastq.$out.fq
    done
  >>>

  output {
    File fq_paired_1 = "~{f1_basename}.fastq.paired.fq.gz"
    File fq_paired_2 = "~{f2_basename}.fastq.paired.fq.gz"
    Array[File] fq_singles = ["~{f1_basename}.fastq.single.fq.gz",
                              "~{f2_basename}.fastq.single.fq.gz"]
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: 4
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 5
  }
}

# Evenly divide paired fastqs to parallelize alignment
task SplitFastq {
  input {
    File fastq
    Int n_shards
    String docker
  }

  Int disk_gb = ceil(size(fastq, "GB") * 3) + 10
  String out_prefix = basename(fastq, ".fq.gz")

  command <<<
    set -eu -o pipefail

    # Build fastqsplitter command
    cmd="fastqsplitter -i ~{fastq}"
    for i in $( seq 1 ~{n_shards} ); do
      cmd="$cmd -o ~{out_prefix}.$i.fq.gz"
    done

    # Split fastqs
    echo -e "Now splitting fastqs with the following command:\n$cmd\n"
    eval $cmd
  >>>

  output {
    Array[File] fq_shards = glob("~{out_prefix}.*.fq.gz")
  }

  runtime {
    docker: docker
    memory: "3.5GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 5
  }
}

# Align single or paired fastqs with BWA MEM
task Bwa {
  input {
    Array[File] input_fastqs
    String output_basename
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa

    Float mem_gb = 32
    String num_cpu = 8

    Int preemptible_tries = 3
    Int? disk_size

    String docker
  }

  Int default_disk_size = ceil(3 * size(input_fastqs, "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Align fastq(s) with bwa
    /usr/gitc/bwa mem -K 100000000 -v 3 -t ~{num_cpu} -Y \
      ~{ref_fasta} \
      ~{sep=" " input_fastqs} \
    | samtools view -b --output-fmt-option level=2 - \
    > ~{output_basename}.bam
  >>>
  
  output {
    File output_bam = "~{output_basename}.bam"
  }

  runtime {
    preemptible: preemptible_tries
    docker: docker
    memory: "~{mem_gb} GiB"
    cpu: num_cpu
    disks: "local-disk " + select_first([disk_size, default_disk_size]) + " HDD"
  }
}