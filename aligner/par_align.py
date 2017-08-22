#!/usr/bin/env python3.5
''' Simple parallelized aligner mockup.
# The basic workflow is structured as follows :
#           Fastq_file
#              |
#         fast_q_split()
#             / \
#         [splitfiles ...]
#           |   |       |
#         bwa() bwa()  bwa()
#           |   |       |
#          [bam_files ...]
#            \  |      /
#             merge()
#               |
#            check_merged()
'''
from os.path import abspath
import argparse
import parsl
import random
from parsl import *

workers = IPyParallelExecutor()
dfk = DataFlowKernel(workers)

#Here we will use the bash split utility to do the task of splitting the fast_q_file
#into n_chunks. We take the abs_path to the input fast_q_file and return a bunch of split
#files whose names are already in the [outputs] array.
@App('bash', dfk)
def fast_q_split(fast_q_file, n_chunks, outputs=[],
                 stdout=abspath('split.out'), stderr=abspath('split.err')):
    cmd_line = '''

    echo "Received file: {0}, splitting into :{1} pieces"
    lines=$(wc -l {0} | awk '{{print $1}}')
    split -l $(($lines / {1})) -d {0} {0}.split.

    '''

# Mocking up the bwa call. We cat the contents of the input file and pipe it through awk
# which just doubles the integers in each line into the output file
@App('bash', dfk)
def bwa(inputs=[], outputs=[], stdout=None, stderr=None):
    cmd_line = '''

    echo "Processing {inputs[0]} -> {outputs[0]}"
    cat {inputs[0]} | awk '{{print $1*2}}' &> {outputs[0]}

    '''

# The merge file simply concatenated the array of files passed in through the inputs list
# into the single file specified though the outputs keyword_arg.
# Note: On the commandline, {inputs} expands to the python list of files expressed as a
# string, eg "['.../bwa.01', '.../bwa.02', ...] which does not make sense in bash.
# So we translate this string with tr to ('.../bwa.01' '.../bwa.02' ...)
@App('bash', dfk)
def merge(inputs=[], outputs=[], stdout=abspath('merge.out'), stderr=abspath('merge.err')):
    cmd_line = '''

    echo "Merging {inputs} -> {outputs[0]}"
    input_array=($(echo {inputs} | tr -d '[],' ))
    cat ${{input_array[*]}} &> {outputs[0]}

    '''

# A simple python app function that takes the result from the merged data file
# and totals up the integers on each line.
@App('python', dfk)
def check_merged(inputs=[]):

    total = 0
    with open(inputs[0], 'r') as tfile :
        for line in tfile.readlines():
            total += int(line)
    return total


if __name__ == "__main__" :

    parser   = argparse.ArgumentParser()
    parser.add_argument("-s", "--split", default='10', help="Number of files to split the fastq file to")
    parser.add_argument("-f", "--fastq", default='mock.fastq', help="File path to the fastq file")
    parser.add_argument("-v", "--verbose", dest='verbose', action='store_true', help="Verbose output")
    args   = parser.parse_args()

    # Handle command-line options
    if args.verbose:
        parsl.set_stream_logger()

    n_chunks   = int(args.split)
    fastq_file = args.fastq

    # Call the first stage where we split the fastq file
    x, splitfiles = fast_q_split(abspath(fastq_file), n_chunks,
                                 outputs=['mock.fastq.split.{0:02d}'.format(i) for i in range(0, n_chunks)])

    # We pass each of the split fastq file to bwa, and these run in parallel
    bam_files = []
    for splitfile in splitfiles:
        bwa_rslt, bwa_out = bwa(inputs=[splitfile], outputs=[splitfile.filepath.replace('fastq.split', 'bam')],
                                stdout='{0}.out'.format(splitfile.filepath),
                                stderr='{0}.err'.format(splitfile.filepath))
        bam_files.extend(bwa_out)

    # The results from the bwa are in bam_files, and these are passed to the merge app
    status, merged = merge(inputs=bam_files, outputs=[abspath('merged.result')])

    # We do a quick check on the results.
    r = check_merged(inputs=merged)
    print("Total sum : ", r.result())

