config = {
    "sites" : [{
        "site" : "Midway_SB",
        "auth" : {
            "channel" : "local"
        },
        "execution" : {
            "executor" : "ipp",
            "provider" : "slurm",
            "script_dir" : ".scripts",
            "block" : {
                "nodes" : 1,
                "taskBlocks" : 1,
                "walltime" : "00:30:00",
                "initBlocks" : 1,
                "minBlocks" : 0,
                "maxBlocks" : 1,
                "scriptDir" : ".",
                "options" : {
                    "partition" : "sandyb",
                    "overrides" : '''module unload python;
module load python/3.5.2+gcc-4.8;
source /scratch/midway2/yadunand/parsl_env_3.5.2_gcc/bin/activate'''
                  }
              }
          }
        },
        {
        "site" : "Midway_GPU",
        "auth" : {
            "channel" : "local"
        },
        "execution" : {
            "executor" : "ipp",
            "provider" : "slurm",
            "script_dir" : ".scripts",
            "block" : {
                "nodes" : 1,
                "taskBlocks" : 1,
                "walltime" : "01:00:00",
                "initBlocks" : 1,
                "minBlocks" : 0,
                "maxBlocks" : 1,
                "scriptDir" : ".",
                "options" : {
                    "partition" : "gpu2",
                    "overrides" : '''#SBATCH --gres=gpu:1
module unload python;
module load python/3.5.2+intel-16.0;
module load cuda/8.0;
module unload intel;
module load intel/17.0;
source /scratch/midway2/yadunand/parsl_env_3.5.2_gcc/bin/activate'''
                  }
              }
          }
        }
        ],
    "globals" : {   "lazyErrors" : True }
}
