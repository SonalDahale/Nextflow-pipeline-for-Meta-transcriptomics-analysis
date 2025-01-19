#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "/*/SRR16378842_*{1,2}.fastq"
params.outdir="/*/results/"
params.threads = 2
params.slidingwindow = "SLIDINGWINDOW:4:15"
params.avgqual = "AVGQUAL:30"
params.genome_dir = "/*/ecoli_ref/ecoli_index_dir/"
params.container = "/*/singularity_trinity/trinityrnaseq.v2.15.2.simg"
params.busco_container = "/*/busco_v5.8.2_cv1.sif"
params.busco_db = "/*/bacteria_odb12/"
params.karken2_container = "/*/kraken2_latest.sif"
params.karken2_db = "/*/minikraken2_v1_8GB/"

//Define channels for raw reads from RNA seq
reads_ch = Channel
    .fromFilePairs(params.reads, checkIfExists:true)

// Define channels for genome index
genome_dir_ch = Channel
    .fromPath(params.genome_dir, checkIfExists:true)

// Process trimmomatic

process trimmomatic {
  
  publishDir "${params.outdir}/trimmed-reads-${sample}", mode: 'copy'
  
  // Same input as fastqc on raw reads, comes from the same channel.
  
  input:
  tuple val(sample), path(reads)
  
  output:
  tuple val("${sample}"), path("${sample}*_P.fq"), emit: paired_fq
  tuple val("${sample}"), path("${sample}*_U.fq"), emit: unpaired_fq

  script:
  """
  trimmomatic PE -threads ${params.threads} ${reads[0]} ${reads[1]} ${sample}1_P.fq ${sample}1_U.fq ${sample}2_P.fq ${sample}2_U.fq ${params.slidingwindow} ${params.avgqual}
  """
} 
//Here we are just confirming our bacteiral strain by running whole data against Karken2 db. We will not be submitting this data to STAR.
process kraken2_contamination_check {

  publishDir "${params.outdir}/kraken2_output", mode: 'copy'

  input:
  tuple val(sample), path(paired_fq)

  output:
  path("${sample}_contam_report.txt"), emit: kraken_report

  script:
  """
  # Run Kraken 2 for contamination detection
  singularity exec -e ${params.karken2_container} kraken2 --db ${params.karken2_db} --paired ${paired_fq[0]} ${paired_fq[1]} \
    --output ${sample}_contam_report.txt \
    --report ${sample}_contam_report.txt
  """
}   
//We are running STAR mapping using E. coli as reference based on our wet lab results. We cross checked this by running Karken2.
process star_genome_mapping {
  
  publishDir "${params.outdir}/umapped_&_mapped_read", mode: 'copy'
  
  input:
  path(genome_dir)
  tuple val(sample), path(paired_fq)

  output:
  path("${sample}_Aligned.sortedByCoord.out.bam"), emit: bacterial_bam
  path("${sample}_Unmapped.out.mate1"), emit: unmapped_left
  path("${sample}_Unmapped.out.mate2"), emit: unmapped_right

  script:
  """
  STAR --runThreadN 16 \
  --genomeDir ${genome_dir} \
  --readFilesIn ${paired_fq[0]} ${paired_fq[1]} \
  --outFileNamePrefix ${sample}_ \
  --outSAMtype BAM SortedByCoordinate \
  --outReadsUnmapped Fastx 
  """
}

process bam_to_fastq {
    
    publishDir "${params.outdir}/bam_to_fastq", mode: 'copy'

    input:
    path(bacterial_bam)

    output:
    path("reads_1.fastq"), emit: bact_left_read
    path("reads_2.fastq"), emit: bact_right_read

    script:
    """
    bedtools bamtofastq -i ${bacterial_bam} -fq reads_1.fastq -fq2 reads_2.fastq
    """
}

process trinity_assembly {
    
    publishDir "${params.outdir}/trinity_output", mode: 'copy'
    
    // Use Singularity container
 
    input:
    path(bact_left_read)
    path(bact_right_read)

    output:
    // The directory for Trinity output
    path ("trinity_out_dir"), emit: trinity_output

    script:
    """
    # Run Trinity assembly using Singularity container
    singularity exec -e ${params.container} Trinity \
    --seqType fq \
    --left ${bact_left_read} \
    --right ${bact_right_read} \
    --max_memory 2G --CPU 8 \
    --no_normalize_reads\
    --output trinity_out_dir
    
    mv *Trinity.fasta* trinity_out_dir/
    """
}

process trinity_bowtie2_indexing {

    publishDir "${params.outdir}/trinity_indexing", mode: 'copy'

    // Use Singularity container

    input:
    path(trinity_output)

    output:
    // The directory for Trinity output
    path ("Trinity.*.bt2"), emit: trinity_indexing

    script:
    """
    # Run Trinity assembly using Singularity container
    singularity exec -e ${params.container} bowtie2-build ${trinity_output}/trinity_out_dir.Trinity.fasta Trinity
   
    """
}
process trinity_bowtie2_assesment {

    publishDir "${params.outdir}/trinity_assesment", mode: 'copy'

    // Use Singularity container

    input:
    path(trinity_indexing)
    path(bact_left_read)
    path(bact_right_read)

    output:
    // The directory for Trinity output
    path ("bowtie2.bam"), emit: trinity_bam
    path("align_stats.txt"), emit: trinity_align_stats

    script:
    """
    index_base=\$(dirname ${trinity_indexing[0]})/Trinity
    # Run Trinity assembly using Singularity container
    singularity exec -e ${params.container} bowtie2 -p 10 -q --no-unal -k 20 -x \${index_base}\
      -1 ${bact_left_read} -2 ${bact_right_read} 2> align_stats.txt | \
      samtools view -@10 -Sb -o bowtie2.bam
    """
}
process busco_analysis {

    publishDir "${params.outdir}/busco_output", mode: 'copy'

    input:
    path(trinity_output)

    output:
    path("busco_output_dir"), emit: busco_results

    script:
    """
    singularity exec -e ${params.busco_container} busco -i ${trinity_output}/trinity_out_dir.Trinity.fasta \
          -o busco_output_dir \
          -m transcriptome \
          -l ${params.busco_db}
    """
}
process trinity_cele_assembly {

    publishDir "${params.outdir}/trinity_output_cele", mode: 'copy'

    // Use Singularity container

    input:
    path(unmapped_left)
    path(unmapped_right)

    output:
    // The directory for Trinity output
    path ("trinity_out_dir_cele"), emit: trinity_output

    script:
    """
    # Run Trinity assembly using Singularity container
    singularity exec -e ${params.container} Trinity \
    --seqType fq \
    --left ${unmapped_left} \
    --right ${unmapped_left} \
    --max_memory 2G --CPU 8 \
    --no_normalize_reads\
    --output trinity_out_dir_cele

    mv *Trinity.fasta* trinity_out_dir_cele/
    """
}
// Running a workflow with the defined processes here.  
workflow {
  trimmomatic(reads_ch)
  kraken2_contamination_check(trimmomatic.out.paired_fq)
  star_genome_mapping(genome_dir_ch, trimmomatic.out.paired_fq)
  bam_to_fastq(star_genome_mapping.out.bacterial_bam)
  trinity_assembly(bam_to_fastq.out.bact_left_read, bam_to_fastq.out.bact_right_read)
  trinity_bowtie2_indexing(trinity_assembly.out.trinity_output)
  trinity_bowtie2_assesment(trinity_bowtie2_indexing.out.trinity_indexing, bam_to_fastq.out.bact_left_read, bam_to_fastq.out.bact_right_read)
  busco_analysis(trinity_assembly.out.trinity_output)
  trinity_cele_assembly(star_genome_mapping.out.unmapped_left, star_genome_mapping.out.unmapped_right)
}

workflow.onComplete {
    println "Pipeline completed at: ${workflow.complete}"
    println "Time to complete workflow execution: ${workflow.duration}"
    println "Execution status: ${workflow.success ? 'Succesful' : 'Failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
