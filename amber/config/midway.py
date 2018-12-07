from parsl.providers import SlurmProvider
from parsl.channels import LocalChannel
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SrunLauncher
from parsl.launchers import SingleNodeLauncher
from parsl.addresses import address_by_hostname


from parsl.data_provider.scheme import GlobusScheme


midway_htex = Config(
            executors=[
                HighThroughputExecutor(
                    label="midway_cpu",
                    worker_debug=True,
                    address=address_by_hostname(),
                    #address='172.25.180.72', # infiniband interface on midway
                    provider=SlurmProvider(
                        'sandyb',
                        launcher=SingleNodeLauncher(),
                        worker_init='source activate parsl_dev',
                        init_blocks=1,
                        max_blocks=1,
                        min_blocks=1,
                        nodes_per_block=1,
                        walltime='1:30:00'
                    ),
                ),
                HighThroughputExecutor(
                    label="midway_gpu",
                    worker_debug=True,
                    address=address_by_hostname(),
                    #address='172.25.180.72', # infiniband interface on midway
                    provider=SlurmProvider(
                        'gpu2',
                        launcher=SingleNodeLauncher(),
                        scheduler_options='''#SATCH --gres=gpu:1''',
                        worker_init='''source activate parsl_dev''',
                        init_blocks=1,
                        max_blocks=1,
                        min_blocks=1,
                        nodes_per_block=1,
                        walltime='1:30:00'
                    ),
                )
            ],
            strategy=None
        )
