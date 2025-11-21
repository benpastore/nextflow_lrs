#!/usr/bin/env python3

import pandas as pd
import argparse
import os

def process_design(file) :

    df = pd.read_csv(file, sep = "\t")

    # make csv for fastq channel
    df['simple_name'] = df.apply(lambda row: os.path.basename(row['Reads']).split(".")[0], axis = 1)
    df['condition'] = df.apply(lambda row: f"{row['Condition']}", axis = 1)

    fastq = df[['simple_name', 'Reads']]

    # to join replicates
    replicates = df[['simple_name', 'condition']]

    # control_association
    fastq.to_csv("fastq.csv", sep = ',', header = True, index = False)
    replicates.to_csv("replicates.csv", sep = ',', header = True, index = False)

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-input", type = str, required = False)

    return parser.parse_args()

def main() : 

    args = get_args() 
    process_design(args.input)

if __name__ == "__main__" : 

    main()