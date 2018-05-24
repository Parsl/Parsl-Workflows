#!/usr/bin/env python3.5
''' Simple parallelized imputation mockup
# The basic workflow is structured as follows :
#      Ref Data per chromosome
#             / \
#         [per chrom ...]
#           |   |       |
#     shapeit2() ()  shapeit2()
#           |   |       |
#    [shapeit2 output_files ...]
#            \  |      /
#      impute2()  () impute2()
#            |   |  | 
#        [per 5GB region ]
#                |
#            check_merged()
'''
from os.path import abspath
import glob
import argparse
import parsl
import random
from parsl import *
import os

workers = ThreadPoolExecutor(max_workers=4)
dfk = DataFlowKernel(executors=[workers])

#workers = IPyParallelExecutor()
#dfk = DataFlowKernel(workers)

# submitting shpaeit tool to generate a hap and sample file output
@App('bash', dfk)
def shapeit(bed_inputs=[], bim_inputs=[], fam_inputs=[], map_inputs=[], 
            outputs=[]):
    cmd_line = '''

    echo "Received file: {bed_inputs[0]}, {bim_inputs[0]}, {fam_inputs[0]}, {gmap_inputs[0]}"
    shapeit --input-bed {bed_inputs[0]} {bim_inputs[0]} {fam_inputs[0]} --input-map {map_inputs[0]} --thread 8 --output-max {haps_outputs[0]} out.sample

    '''

    #print(cmd_line)
    return cmd_line

# Mocking up the bwa call. We cat the contents of the input file and pipe it through awk
# which just doubles the integers in each line into the output file
@App('bash', dfk)
def impute(hap_inputs=[], map_inputs=[], legend_inputs=[], haps_inputs=[], 
           impute_outputs=[], region=[]):
    cmd_line = '''

    echo "Processing {region[0]}: {region[1]} - {region[2]} -> {outputs[0]}"
    impute2 -m {map_inputs[0]} -h {hap_inputs[0]} -l {legend_inputs[0]} -known_haps_g {haps_inputs[0]} -use_prephased_g -int {region[1]} {region[2]} -Ne 20000 -buffer 250 -o {impute_outputs[0]}

    '''

    return cmd_line

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

    return cmd_line

def get_dir_files(directory, data_type):
    dir_files = []
    for i in sorted(glob.glob("%s/*.%s" % (directory, data_type))):
        if "chrALL_" in i or "chrMT_" in i or "chrX_" in i or "chrXY_" in i or "chrY_" in i:
            continue
        else:
            dir_files.append(i)
    return dir_files

def get_regions(regions_file, window_size):
    regions = []
    window_size = int(window_size)
    fh = open(regions_file, "r")
    for line in fh:
        values = line.split("\t")
        chr_name =  values[0]
        max_bp = int(values[1])
        start = 1
        end = window_size
        while end < max_bp:
            regions.append([chr_name, start, end])
            #print "%s\t%s\t%s" % (chr_name, start, end)
            start = end + 1
            end = end + window_size
        end = max_bp
        regions.append([chr_name, start, end])
    return regions

def get_chr_file(region_name, directory, separator, data_type):
    #print("%s/*_%s%s*%s" % (directory, region_name, separator, data_type))
    #directory = os.path.dirname(contents[0])
    results = glob.glob("%s/*_%s%s*%s" % (directory, region_name, separator, data_type))
    for i in results:
        return i
    return None

    

if __name__ == "__main__" :

    parser   = argparse.ArgumentParser()
    parser.add_argument("-b", "--bed-dir", dest='bed_dir', help="BED files directory")
    parser.add_argument("-m", "--bim-dir", dest='bim_dir', help="BIM files directory")
    parser.add_argument("-f", "--fam-dir", dest='fam_dir', help="FAM files directory")
    parser.add_argument("-g", "--map-dir", dest='map_dir', help="MAP files directory")
    parser.add_argument("-r", "--regions", dest="regions_file", help="genome dictionary file")
    parser.add_argument("-z", "--region-size", dest="region_size", help="region size to split chromosome")
    parser.add_argument("-a", "--hap-dir", dest="hap_dir", help="HAP files directory")
    parser.add_argument("-l", "--legend-dir", dest="legend_dir", help="Legend files directory")
    parser.add_argument("-v", "--verbose", dest='verbose', action='store_true', help="Verbose output")
    parser.add_argument("-o", "--output-dir", dest="output_dir", help="output directory")
    args   = parser.parse_args()

    # Handle command-line options
    if args.verbose:
        parsl.set_stream_logger()

    # Figure out the BED, BIM, FAM, MAP files from the given directories
    # for each chromosome. Assume all files are present
    bed_inputs = get_dir_files(args.bed_dir, "bed")
    bim_inputs = get_dir_files(args.bim_dir, "bim")
    fam_inputs = get_dir_files(args.fam_dir, "fam")
    map_inputs = get_dir_files(args.map_dir, "txt")

    # Call the shapeit tool on each chromosome
    shapeit_haps_outs = []
    shapeit_sample_outs = []
    for i in range(0,len(bed_inputs)):
        output_file = "%s/%s" % (args.output_dir, os.path.basename(bed_inputs[i]).replace('bed', 'haps'))
        x = shapeit(bed_inputs=[bed_inputs[i]], bim_inputs=[bim_inputs[i]], fam_inputs=[fam_inputs[i]], map_inputs=[map_inputs[i]], outputs=[output_file])
        shapeit_haps_outs.extend(x.outputs)

    #import sys
    #sys.exit()
    # call the impute2 on each 5MB region per chromosome using the shapeit_haps_outs as inputs
    impute2_outs = []
    regions = get_regions(args.regions_file, args.region_size)
    hap_inputs = get_dir_files(args.hap_dir, "hap")
    legend_inputs = get_dir_files(args.legend_dir, "legend")
    for region in regions:
        #print(region)
        hap_file = get_chr_file(region[0], os.path.dirname(hap_inputs[0]), ".", "hap")
        if hap_file is None:
            continue
        map_file = get_chr_file(region[0], os.path.dirname(map_inputs[0]), "_", "txt")
        legend_file = get_chr_file(region[0], os.path.dirname(legend_inputs[0]), ".", "legend")
        haps_shapeit_file = get_chr_file(region[0], os.path.dirname(shapeit_haps_outs[0].filepath), "_", "haps")
        output_file = "%s/%s" % (args.output_dir, os.path.basename(hap_file.replace('hap', '%s.%s.impute2_out.txt' % (region[1], region[2]))))
        x = impute(hap_inputs=[hap_file], map_inputs=[map_file], legend_inputs=[legend_file],
                               haps_inputs=[haps_shapeit_file], 
                               impute_outputs=[output_file], 
                               region=region)
        impute2_outs.extend(x.outputs)
    #print(impute2_outs)

    # The results from the bwa are in bam_files, and these are passed to the merge app
    x = merge(inputs=impute2_outs, outputs=[abspath('merged.result')])

