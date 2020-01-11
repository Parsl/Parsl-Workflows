#!/usr/bin/env python3


"""
[Notes from Dan]
1 - a python code that runs and creates some directories
    called relax.xy, some directories called neb.xy, and a Makefile.
2 - We need to run VASP in each relax.xy directly - in our test case, there are 4 such
    directories.  Each VASP run will take 10-30 hours on O(100) cores.
3 - once these have finished, we run make, which uses the results in the relax.xy
    directories to build the inputs in the deb.xy directories - perhaps we could
    figure out what Make does and do it in python instead, but likely, we could
    just call Make from python...
4 - We can then run VASP in the deb.xy directories - in our test case, there are 17
    such directories, with similar VASP runtimes as before.
5 - Once these are done, we need to run some more python code that we don't actually
    have yet, but that a student here supposedly does have written and tested.

We will be working on Stampede 2, and everything we used can be installed via pip.

"""

from parsl import *
import os
import shutil
import random

workers = ThreadPoolExecutor(max_workers=8)
dfk = DataFlowKernel(workers)


def create_dirs(cwd):

    for dir in ['relax.01', 'relax.02', 'relax.03'] :
        rel_dir = '{0}/{1}'.format(cwd, dir)
        if os.path.exists(rel_dir):
            shutil.rmtree(rel_dir)
        os.makedirs(rel_dir)
        for i in range(0, random.randint(1,5)):
            rdir = '{0}/{1}'.format(rel_dir, i)
            os.makedirs(rdir)
            with open('{0}/results'.format(rdir, i), 'w') as f:
                f.write("{0} {1} - test data\n".format(i, dir))

    for dir in ['neb01', 'neb02', 'neb03', 'neb04'] :
        rel_dir = '{0}/{1}'.format(cwd, dir)
        if os.path.exists(rel_dir):
            shutil.rmtree(rel_dir)
        os.makedirs(rel_dir)
        with open('{0}/{1}.txt'.format(rel_dir, dir), 'w') as f:
            f.write("{0} test data\n".format(rel_dir))



@App('python', dfk)
def ls (pwd, outputs=[]):
    import os
    items = os.listdir(pwd)
    with open(outputs[0], 'w') as f:
        f.write(' '.join(items)) 
        f.write('\n')
    # Returning list of items in current dir as python object
    return items

@App('bash', dfk)
def catter (dir, outputs=[], stdout=None, stderr=None):
    cmd_line = 'cat {0}/*/results > {outputs[0]}'

if __name__ == "__main__" :

    pwd = os.getcwd()
    create_dirs(pwd)

    # Listing the cwd
    ls_fu, [res] = ls(pwd, outputs=['results'])

    dir_fus = {}
    for dir in ls_fu.result():
        if dir.startswith('relax') and os.path.isdir(dir):
            print("Launching {0}".format(dir))
            dir_fus[dir] = catter(dir, outputs=['{0}.allresults'.format(dir)],
                                  stderr='{0}.stderr'.format(dir))

    for dir in dir_fus:
        try :
            print(dir_fus[dir][0].result())
        except Exception as e :
            print ("Caught exception{0}  on {1}".format(e, dir))
