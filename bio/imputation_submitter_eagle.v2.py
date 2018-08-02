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
import os, sys

workers = ThreadPoolExecutor(max_workers=4)
dfk = DataFlowKernel(executors=[workers])

#workers = IPyParallelExecutor()
#dfk = DataFlowKernel(workers)

# submitting plink tool
# plink2 --threads 16 --chr $chr --bfile chr$chr --export vcf id-paste=iid bgz --out chr$chr.plink
@App('bash', dfk)
def plink2(b_inputs=[], chrom=[], plink_outputs=[]):
    out_prefix = plink_outputs[0].replace(".vcf.gz", "")
    b_prefix = b_inputs[0].replace(".bed", "")
    cmd_line = '''

    echo "Executing region: {chrom[0]}"
    plink2 --threads 16 --chr {chrom[0]} --bfile {b_prefix} --export vcf id-paste=iid bgz --out {out_prefix}

    '''
    #print(cmd_line)
    return cmd_line


# submitting shpaeit tool to generate a hap and sample file output
# eagle --vcfRef=vcfref_file.txt.gz --vcfTarget=chr$chr.plink.vcf.gz --vcfOutFormat=z --geneticMapFile=genetic_map_hg19_withX.txt.gz --chrom=$chrX_str --outPrefix=chr$chr --numThreads=16 2>&1 | tee chr$chr.eagle.log
@App('bash', dfk)
def eagle(plink_inputs=[], region=[], eagle_outputs=[], map_inputs=[], eagle_log_outputs=[], vcfRef_inputs=[]):
    out_prefix = eagle_outputs[0].replace(".vcf.gz", "")
    cmd_line = '''

    echo "Received file: {plink_inputs[0]}, {map_inputs[0]}"
    echo "Processing {region[0]}: {region[1]} - {region[2]} -> {outputs[0]}"
    eagle --vcfRef={vcfRef[0]} --vcfTarget={plink_inputs[0]} --vcfOutFormat=z --geneticMapFile={map_inputs[0]} --chrom={region[0]} --outPrefix={out_prefix} --numThreads=16 2>&1 | tee {eagle_log_outputs[0]}

    '''
    #print(cmd_line)
    return cmd_line

# submitting minimac
# minimac4-omp --refHaps $refhap --haps chr$chr.vcf.gz --format GT,DS,HDS,GP --passOnly --allTypedSites --window 1000000 --prefix chr$chr --log --cpu 16
@App('bash', dfk)
def minimac4(haps_inputs=[], refhaps=[], region=[], out_prefix=[]):
    cmd_line = '''

    echo "Received file: {haps_inputs[0]}"
    echo "Processing {region[0]}: {region[1]} - {region[2]} -> {outputs[0]}"
    minimac4-omp --refHaps {refhaps[0]} --haps {haps_inputs[0]} --format GT,DS,HDS,GP --passOnly --allTypedSites --window 1000000 --prefix {out_prefix[0]} --log --cpu 16

    '''

    #print(cmd_line)
    return cmd_line

# submitting eagle and minimac on each chunk
@App('bash', dfk)
def eagle_minimac(plink_inputs=[], hap_inputs=[], vcfRef_inputs=[], eagle_outputs=[], eagle_log_outputs=[], map_inputs=[], region=[], impute_outputs=[]):                          )
    out_prefix = eagle_outputs[0].replace(".vcf.gz", "")
    out_minimac_prefix = impute_outputs[0].replace(".vcf.gz", "")
    cmd_line = '''

    echo "Received file: {plink_inputs[0]}, {map_inputs[0]}"
    echo "Processing {region[0]}: {region[1]} - {region[2]} -> {outputs[0]}"
    eagle --vcfRef={vcfRef_inputs[0]} --vcfTarget={plink_inputs[0]} --vcfOutFormat=z --geneticMapFile={map_inputs[0]} --bpStart {region[1]} --bpEnd {region[2]} --chrom={region[0]} --outPrefix={out_prefix} --numThreads=36 2>&1 | tee {eagle_log_outputs[0]};
    minimac4-omp --refHaps {hap_inputs[0]} --haps {eagle_outputs[0]} --format GT,DS,HDS,GP --passOnly --allTypedSites --window 50000 --prefix {out_minimac_prefix[0]} --log --cpu 36
    '''
    #print(cmd_line)

    
    
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
    #print("%s/*%s%s*%s" % (directory, region_name, separator, data_type))
    results = glob.glob("%s/*%s%s*%s" % (directory, region_name, separator, data_type))
    for i in results:
        return i
    return None


if __name__ == "__main__" :

    parser   = argparse.ArgumentParser()
    parser.add_argument("-b", "--b-dir", dest='b_dir', help="BED files directory")
    parser.add_argument("-v", "--vcfref-dir", dest='vcfref_dir', help="VCF REF files directory")
    parser.add_argument("-r", "--regions", dest="regions_file", help="genome dictionary file")
    parser.add_argument("-z", "--region-size", dest="region_size", help="region size to split chromosome")
    parser.add_argument("-a", "--genetic-map-file", dest="genmap_file", help="Genetic MAP file reference")
    parser.add_argument("-m", "--hapmap-dir", dest='hapmap_dir', help="hapmap REF files directory")
    parser.add_argument("-v", "--verbose", dest='verbose', action='store_true', help="Verbose output")
    parser.add_argument("-o", "--output-dir", dest="output_dir", help="output directory")
    args   = parser.parse_args()

    # Handle command-line options
    if args.verbose:
        parsl.set_stream_logger()

    # Figure out the BED, BIM, FAM, MAP files from the given directories
    # for each chromosome. Assume all files are present
    b_inputs = get_dir_files(args.b_dir, "bed")
    #hapmap_inputs = get_dir_files(args.hap_dir, "gz")
    hapmap_dir = args.hapmap_dir
    genmap_file = args.genmap_file
    vcfref_inputs = get_dir_files(args.vcfref_dir, "gz")
    output_dir = args.output_dir
    region_size = 1000000

    # Call the plink tool on each chromosome
    plink_outs = []
    chrom_regions = [i for i in range(1,23)]
    for i in chrom_regions:
        output_file = "%s/%s" % (output_dir, os.path.basename(b_inputs[0]).replace('bed', 'plink.chr%s.vcf.gz' % i))
        print(output_file)
        x = plink2(b_inputs=[b_inputs[0]], chrom=[i], outputs=[output_file])
        plink_outs.extend(x.outputs)

    # call eagle tool on each chromosome and output from plink
    minimac_outs = []
    hap_inputs = get_dir_files(hapmap_dir, "m3vcf.gz")
    regions = get_regions(args.regions_file, region_size)

    for region in regions:
        plink_file = get_chr_file("plink.chr" + str(region[0]), os.path.dirname(plink_outs[0].filepath), ".", "gz")
        vcfref_file = get_chr_file("chr" + str(region[0]) + ".phase3_v5", os.path.dirname(vcfref_inputs[0]), "_", "gz")
        #output_file = "%s/%s" % (output_dir, os.path.basename(b_inputs[0]).replace('bed', 'eagle.chr%s.vcf.gz' % str(region[0])))
        #output_log_file = "%s/%s" % (output_dir, os.path.basename(b_inputs[0]).replace('bed', 'eagle.chr%s.log' % str(region[0])))
        hap_file = get_chr_file(region[0], os.path.dirname(hap_inputs[0]), ".", "m3vcf.gz")
        if hap_file is None:
            continue
        eagle_log_output = "%s/%s" % (args.output_dir, os.path.basename(hap_file.replace('m3vcf.gz', '%s.%s.%s.eagle_log_out.txt' % (region[0], region[1], region[2]))))
        eagle_output = "%s/%s" % (args.output_dir, os.path.basename(hap_file.replace('m3vcf.gz', '%s.%s.%s.eagle_out.txt' % (region[0], region[1], region[2]))))
        output_file = "%s/%s" % (args.output_dir, os.path.basename(hap_file.replace('m3vcf.gz', '%s.%s.%s.minimicac_out.txt' % (region[0], region[1], region[2]))))
        x = eagle_minimac(plink_inputs=[plink_file], 
                          hap_inputs=[hap_file],
                          vcfRef_inputs=[vcfref_file], 
                          map_inputs=[map_file],
                          region=region, 
                          impute_outputs=[output_file],
                          eagle_log_outputs=[eagle_log_output],
                          eagle_outputs=[eagle_output]
                          )
        #x = minimac(hap_inputs=[hap_file], map_inputs=[map_file], legend_inputs=[legend_file],
        #                       haps_inputs=[haps_shapeit_file],
        #                       impute_outputs=[output_file],
        #                       region=region)
        minimac_outs.extend(x.outputs)
 

    #print(minimac_outs)

    # The results from the bwa are in bam_files, and these are passed to the merge app
    #x = merge(inputs=minimac_outs, outputs=[abspath('merged.result')])
