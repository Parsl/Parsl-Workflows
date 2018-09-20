#!/usr/bin/env python3.5
from os.path import abspath
import glob
import argparse
import parsl
import random
from parsl import *
import os, sys

from parsl.config import Config
from parsl.executors.ipp import IPyParallelExecutor 
from libsubmit.providers import LocalProvider
from libsubmit.providers import TorqueProvider 
from libsubmit.channels import LocalChannel, SSHChannel
from libsubmit.launchers import SingleNodeLauncher
from libsubmit.launchers import MpiExecLauncher

'''
config = {
	"sites":[{
		"site": "MurphyTest",
		"auth":{"channel": "local","hostname":"va-murphy-login.kdi.local", "scriptDir":"/home/maddurir"},
		"execution": {	
			"executor":"ipp",
			"provider": "torque",
			"block": {
				"nodes": 2,
				"taskblocks": 20,
				"walltime":"00:01:00",
				"initBlocks": 1,
				"minBlocks": 1,
				"maxBlocks":1,
				"options": {
					"partition":"debug","overrides":"#PBS -l nodes=2:ppn=36"}
			}
		}
	}],
	"global": {"lazyErrors": True, "appCache": True, "retries":2}
}

'''

config = Config(
	executors = [
		IPyParallelExecutor(
		label='PLINKandEagle',
		provider=TorqueProvider(
			nodes_per_block = 1,
			tasks_per_node = 1,
			init_blocks = 1,
			#init_blocks = 1,
			#max_blocks=1,
			max_blocks=1,
			min_blocks=1,
			overrides='#PBS -l nodes=1;ppn=36',
			queue='batch',
			channel=LocalChannel(),
			launcher=SingleNodeLauncher(),
			walltime='4000:10:00'
		)
	)],

	retries=1
)

minimap_config = Config(
	executors = [
                IPyParallelExecutor(
		label='PLINKandEagle',
		provider=TorqueProvider(
			nodes_per_block = 1,
			tasks_per_node = 1,
			init_blocks = 0,
			max_blocks=25,
			min_blocks=0,
			overrides='module load apps/openmpi/gnu/3.0.0',
			queue='batch',
			channel=SSHChannel(hostname = "va-murphy-login.kdi.local",
                                           username = os.getenv("USER").split('@')[0],
                                           script_dir = os.getenv("HOME") + "/code-va/parsl-workflows/ssh_scripts"
                                          ),
			launcher=MpiExecLauncher(),
			walltime='240:00:00'
		  )
	        ),
		IPyParallelExecutor(
		label='minimac',
		provider=TorqueProvider(
			nodes_per_block = 1,
			tasks_per_node = 8,
			init_blocks = 0,
			max_blocks=30,
			min_blocks=0,
			overrides='module load apps/openmpi/gnu/3.0.0',
			queue='batch',
			channel=SSHChannel(hostname = "va-murphy-login.kdi.local",
                                           username = os.getenv("USER").split('@')[0],
                                           script_dir = os.getenv("HOME") + "/code-va/parsl-workflows/ssh_scripts"
                                          ),
			launcher=MpiExecLauncher(),
			walltime='4000:00:00'
		)
	)],

	retries=0
)

#local_config = Config(executors=[IPyParallelExecutor(label="local_ipp",provider=LocalProvider(channel=LocalChannel(),init_blocks=23,max_blocks=23))])

dfk = parsl.load(minimap_config)


#parsl.load(local_config)
#workers = ThreadPoolExecutor(max_workers=4)
#dfk = DataFlowKernel(executors=[workers])
#dfk = DataFlowKernel(config = config)

#workers = IPyParallelExecutor()
#dfk = DataFlowKernel(workers)

# submitting plink tool
# plink2 --threads 36 --chr $chr --bfile chr$chr --export vcf id-paste=iid bgz --out chr$chr.plink
@App('bash', executors=['PLINKandEagle'], cache=True)
def plink2(b_inputs=[], chrom=[], outputs=[],stdout=None, stderr=None):
    out_prefix = outputs[0].replace(".vcf.gz", "")
    b_prefix = b_inputs[0].replace(".bed", "")
    cmd_line = '''

    echo "Executing region: {chrom[0]}"
    module load apps/plink/2.00a
    plink2 --threads 36 --chr {chrom[0]} --bfile %s --export vcf id-paste=iid bgz --out %s

    ''' %(b_prefix,out_prefix)
    #print(cmd_line)
    return cmd_line


# submitting eagle to generate a hap and sample file output
# eagle --vcfRef=vcfref_file.txt.gz --vcfTarget=chr$chr.plink.vcf.gz --vcfOutFormat=z --geneticMapFile=genetic_map_hg19_withX.txt.gz --chrom=$chrX_str --outPrefix=chr$chr --numThreads=16 2>&1 | tee chr$chr.eagle.log
@App('bash', executors=['PLINKandEagle'], cache=True)
def eagle(plink_file, genmap_file, vcfRef, region, outputs=[], vcfRef_inputs=[],stdout=None, stderr=None):
    out_prefix = outputs[0].replace(".vcf.gz", "")
    cmd_line = '''

    echo "Received file: {plink_file}, {genmap_file}"
    module load apps/eagle/2.4
    echo "Processing {region}: -> {outputs[0]}"
    eagle --vcfRef={vcfRef} --vcfTarget={plink_file} --vcfOutFormat=z --geneticMapFile={genmap_file} --chrom={region} --outPrefix=%s --numThreads=36

    ''' %(out_prefix)
    print("EAGLE COMMAND LINE:",cmd_line)
    return cmd_line

@App('bash', executors=['PLINKandEagle'], cache=True)
def eagle2(plink_inputs=[], region=[], outputs=[],map_inputs=[],vcfRef_inputs=[],stdout=None, stderr=None):
    out_prefix = outputs[0].replace(".vcf.gz","")
    cmd_line = '''
    echo "Received file: {plink_inputs[0]}, {map_inputs[0]}"
    module load apps/eagle/2.4
    echo "Processing {region[0]} -> {outputs[0]}"
    eagle  --vcf={plink_inputs[0]} --vcfOutFormat=z --geneticMapFile={map_inputs[0]} --chrom={region[0]} --outPrefix=%s --numThreads=36
    ''' % (out_prefix)
    print("EAGLE COMMAND LINE:", cmd_line)
    return cmd_line

# submitting minimac
# minimac4-omp --refHaps $refhap --haps chr$chr.vcf.gz --format GT,DS,HDS,GP --passOnly --allTypedSites --window 1000000 --prefix chr$chr --log --cpu 16
@App('bash', executors=['minimac'])
def minimac4(haps_inputs=[], refhaps=[], region=[], prefix=[], outputs=[], stdout=None, stderr=None):
    cmd_line = '''
    echo "Received file: {haps_inputs[0]}"
    echo "Intermediate destination: {prefix[0]}"
    echo "Final destination: {outputs[0]}"
    touch {prefix[0]}
    module load apps/minimac4/f398ec6
    echo "Processing {region[0]}: {region[1]} - {region[2]} -> {prefix[0]}"
    minimac4 --noPhoneHome --refHaps {refhaps[0]} --haps {haps_inputs[0]} \
      --format GT,DS,HDS,GP --allTypedSites --chr {region[0]} \
      --start {region[1]} --end {region[2]} --window 1000000 \
      --topThreshold 0.01 --prefix {prefix[0]} --log --vcfBuffer 2000 --cpus 8 || \
      echo "Ignoring minimac exit code"
    if [ ! -f {prefix[0]}.dose.vcf.gz ] ; then
        echo "Minimac produced no results for {region[0]}, {region[1]}, {region[2]}. Check for errors."
    fi
    mv {prefix[0]}* $(dirname {outputs[0]})
    '''

    print(cmd_line)
    return cmd_line

# The merge file simply concatenated the array of files passed in through the inputs list
# into the single file specified though the outputs keyword_arg.
# Note: On the commandline, {inputs} expands to the python list of files expressed as a
# string, eg "['.../bwa.01', '.../bwa.02', ...] which does not make sense in bash.
# So we translate this string with tr to ('.../bwa.01' '.../bwa.02' ...)
@App('bash')
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
    print("RETURNING RESULTS", results)
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
    parser.add_argument("-p", "--p-dir", dest='plink_dir', help="PLINK RESULTS DIRECTORY")
    parser.add_argument("-e", "--e-dir", dest='eagle_dir', help="EAGLE RESULTS DIRECTORY")
    parser.add_argument("-vb", "--verbose", dest='verbose', action='store_true', help="Verbose output")
    parser.add_argument("-o", "--output-dir", dest="output_dir", help="Minimac output directory")
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
    region_size = args.region_size
    plink_outs = []
    plink_futures = []
    plink_dir = args.plink_dir
    eagle_dir = args.eagle_dir
    # Call the plink tool on each chromosome
    if plink_dir is not None and os.path.exists(plink_dir):
        plink_outs = glob.glob("%s/*vcf.gz" % plink_dir)
    # MVP Genotype data has chromosomes 24,25,26
    new_chrom_regions = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,25,26]
    chrom_regions = [i for i in range(1,23)]
    # Do smallest chromosomes first for quicker feedback
    chrom_regions.reverse()
    if len(plink_outs) < len(chrom_regions):
        for i in chrom_regions:
           output_file = "%s/%s" % (plink_dir, os.path.basename(b_inputs[0]).replace('bed', 'plink.chr%s.vcf.gz' % i))
           print(output_file)
           x = plink2(b_inputs=[b_inputs[0]], chrom=[i], outputs=[output_file], stdout='plink2.{}.out'.format(i),stderr='plink2.{}.err'.format(i))
           #x.result()
           plink_futures.append(x)
           plink_outs.append(output_file)
        for f in plink_futures:
           f.result()
    else:
        print("PLINK IS ALREADY RUN", plink_dir)
        print(plink_outs)

    # call eagle tool on each chromosome and output from plink
    #sys.exit()
    eagle_outs = []
    if eagle_dir is not None and os.path.exists(eagle_dir):
        eagle_outs = glob.glob("%s/*vcf.gz" % eagle_dir)

    eagle_futures = []
    if len(eagle_outs) < len(chrom_regions):
        for i in chrom_regions:
            plink_file = get_chr_file("plink.chr" + str(i), os.path.dirname(plink_outs[0]), ".", "gz")
            print ("PLINK_FILE:", plink_file)
            print("VCFREF_INPUTS[0]",vcfref_inputs[0])
            #vcfref_file = get_chr_file("ALL.chr" + str(i) + ".phase3_v5", os.path.dirname(vcfref_inputs[0]), "_", "gz")
            vcfref_file = "/lustre/valustre/mvp-core/world-shared/1000GENOMES/ALL.chr" +str(i)+".phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz"
            #VCFREF_INPUTS[0] /lustre/valustre/mvp-core/world-shared/1000GENOMES/1.1000g.Phase1.v3.No.Parameter.Estimates.m3vcf.gz

            #/lustre/valustre/mvp-core/world-shared/1000GENOMES/ALL.chr4.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz
            print ("VCFREF:", vcfref_file)
            print ("GENMAP_File:", genmap_file)
            output_file = "%s/%s" % (eagle_dir, os.path.basename(b_inputs[0]).replace('bed', 'eagle.chr%s.vcf.gz' % i))
            print ("Output File:", output_file)
            x = eagle2(plink_inputs=[plink_file], region=[i], outputs=[output_file], map_inputs=[genmap_file], vcfRef_inputs=[vcfref_file], stdout='eagle2.{}.out'.format(i),stderr='eagle2.{}.err'.format(i))
        #x = eagle(plink_file, genmap_file, vcfref_file, i, outputs=[output_file],stdout='eagle2.{}.out'.format(i),stderr='eagle2.{}.err'.format(i))
            eagle_futures.append(x)
            eagle_outs.append(output_file)
        for f in eagle_futures:
            f.result()
    else:
        print("EAGLE IS ALREADY RUN", eagle_dir)
        #sys.exit()

    # We're done with plink and eagle, so we can shut down that executor.
    #num_jobs = len(dfk.executors['PLINKandEagle'].provider.resources.keys())
    #dfk.executors['PLINKandEagle'].scale_in(num_jobs)
    #dfk.executors['PLINKandEagle'].shutdown()

    # call the minimac tool on each 1MB region per chromosome using the eagle_outs as inputs
    #sys.exit()
    minimac_outs = []
    minimac_futures = []
    print("HAP DIR:", hapmap_dir)
    regions = get_regions(args.regions_file, region_size)
    print (regions)
    #sys.exit()
    hap_inputs = get_dir_files(hapmap_dir, "vcf.gz")
    print("Running MINIMAC")
    print("HAP INPUTS:", hap_inputs)
    for region in regions:
        #print(region)
        try:
            region[0] = int(region[0][3:])
            if region[0] not in chrom_regions:
                continue
        except ValueError:
            continue

        hap_file = get_chr_file(region[0], os.path.dirname(hap_inputs[0]), ".", "vcf.gz")
        if hap_file is None:
            continue
        #haps_eagle_file = get_chr_file("eagle.chr" + region[0], os.path.dirname(eagle_outs[0].filepath), ".", "vcf.gz")
        #print ("HAP EAGLE FILE:", haps_eagle_file)
        #/lustre/valustre/mvp-core/world-shared/1000GENOMES/22.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz
        refhap_file = "/lustre/valustre/mvp-core/world-shared/1000GENOMES/" +str(region[0])+".1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz"
        print("REFHAP_FILE:", refhap_file)
        output_dir = "%s/%s" % (args.output_dir, os.path.basename(hap_file.replace('vcf.gz', '%s.%s.%s.minimimac_out' % (region[0], region[1], region[2]))))
        print("MINIMAC OUTPUT:", output_dir)
        #def minimac4(haps_inputs=[], refhaps=[], region=[], outputs=[], stdout=None, stderr=None):
        # Put the intermediate results in localscratch
        out_dest_dir = str(os.path.split(output_dir)[0])
        out_basename = str(os.path.split(output_dir)[1])
        prefix = "$localscratch" + '/{}'.format(out_basename)
        print("MINIMAC PREFIX:", prefix)
        #def minimac4(haps_inputs=[], refhaps=[], region=[], outputs=[], stdout=None, stderr=None):
        x = minimac4(haps_inputs=[hap_file], refhaps=[refhap_file],region=region, 
                prefix=[prefix], outputs=[output_dir],
                stdout='minimac4.{}_{}_{}.out'.format(region[0], region[1], region[2]),
                stderr='minimac4.{}_{}_{}.err'.format(region[0], region[1], region[2]))
        minimac_futures.append(x)
        minimac_outs.append(output_dir)
    #print(minimac_outs)
    for f in minimac_futures:
        f.result()
    # The results from the bwa are in bam_files, and these are passed to the merge app
    #x = merge(inputs=minimac_outs, outputs=[abspath('merged.result')])
