![Pipeline_workflow](https://github.com/SonalDahale/Nextflow-pipeline-for-Meta-transcriptomics-analysis/blob/main/Metatrasncritpomic_nextflow.png)

############# Gut Microbial Metatranscriptomics Pipeline #############

This metatranscriptomics pipeline processes RNA-seq data, specifically bacterial contaminants in organisms like nematodes. The pipeline includes quality control, contamination checks, transcriptome assembly, and quality assessment, and is tailored to identify and filter out known contaminants.

Key Features:

- Quality Control: Uses Trimmomatic to trim and filter raw RNA-seq reads.\
- Contamination Check: Runs Kraken2 to detect microbial contaminant abundance.\
- Genome Alignment: Maps reads to a reference genome using STAR.I am mapping reads to E. coli genome as we used its pure culture for infecting nematode.\
- Transcriptome Assembly: Performs de novo transcript assembly with Trinity.\
- Quality Assessment: Assesses transcriptome assembly quality using BUSCO.\
- Flexibility: Supports containerized tools via Singularity for reproducibility.

Installation:

Install Nextflow: Follow the Nextflow Installation Guide.

The pipeline uses the following tools. Ensure they are installed and accessible in your environment:

Trimmomatic  
STAR

Singularity Containers:

The following tools are used via Singularity containers to ensure reproducibility, portability, simplified dependencies, and efficient use of storage:

Trinity  
Kraken2  
BUSCO  

Usage:

Modify the nextflow.config file to specify your input data paths and parameters.

Run the pipeline with the following command:

nextflow run main.nf

or on your cluster using:

qsub Nextflow_job_submission.sh

Output files will be saved in the directory specified by params.outdir.

Parameters:

Update the following parameters in nextflow.config as needed:

- params.reads: Path to raw RNA-seq reads (supports paired-end reads).\
- params.outdir: Directory for pipeline outputs.\
- params.genome_dir: Path to the reference genome index directory.\
- params.container: Path to the Singularity container for Trinity.\
- params.busco_container: Path to the Singularity container for BUSCO.\
- params.busco_db: Path to the BUSCO database for quality assessment.\
- params.kraken2_container: Path to the Singularity container for Kraken2.\
- params.kraken2_db: Path to the Kraken2 database.

Outputs:

- Trimmed Reads: Quality-controlled FASTQ files.\
- Contamination Report: Kraken2 contamination summary highlighting lab-origin contaminants.\
- Mapped and Unmapped Reads: BAM files and unmapped FASTQ reads*.\
- Transcriptome Assembly: Assembled transcripts in FASTA format.\
- Quality Metrics: BUSCO results and Bowtie2 alignment statistics.\

Continuous Updates:

This pipeline is actively under development, with new features and analyses being continuously added. Stay tuned for updates, including functional annotation, taxonomic profiling, and downstream analysis.

Note on Unmapped FASTQ Reads:

Unmapped FASTQ reads are mRNA reads from the host (nematode in this case). Depending on the availability of a reference genome, these reads will be used for either reference-based or de novo assembly. These reads will also be employed for functional annotation, GSEA, and downstream analysis
