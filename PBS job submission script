#!/usr/bin/env bash

#PBS -lselect=1:ncpus=32:mem=256gb
#PBS -l walltime=48:00:00
#PBS -o nf_output4.txt
#PBS -e nf_error4.txt

module load anaconda3/personal
source activate nextflow_env

export JAVA_HOME=~/anaconda3/envs/nextflow_env
export PATH=$JAVA_HOME/bin:$PATH

~/anaconda3/envs/nextflow_env/bin/nextflow run /*/nextflow/bulk/main.nf -resume -with-dag main.png
