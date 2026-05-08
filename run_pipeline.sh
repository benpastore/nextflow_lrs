#!/bin/bash

input_file=$PWD/family_14.txt
output_directory=$PWD/results_family_14

#######################DONT MODIFY##############################################
script=/fs/ess/PAS0631/00_vicki/pipelines/nextflow_lrs/workflows/nanopore/run.sh

sh $script --design $input_file --results $output_directory --medaka