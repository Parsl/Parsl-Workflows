import parsl
from parsl.app.app import python_app, bash_app
from config.midway import midway_htex
from  jinja2 import Template, Environment
from jinja2.loaders import FileSystemLoader
import os

parsl.set_stream_logger()

parsl.load(midway_htex)
# Helper function to create an input file from a template
def create_template(template_name, output_name, contents):
    fs_loader = FileSystemLoader('config_files/templates')
    env = Environment(loader=fs_loader)
    template = env.get_template(template_name)
    t_path = os.path.join("config_files", output_name)
    t_file = open(t_path, 'w')
    t_file.write(template.render(contents))
    t_file.close()
    return t_path


#################   
# App definitions
#################
@bash_app(executors=['midway_cpu'])
def packmol(input_file=None, inputs=[], outputs=[], stdout=None, stderr=None):
    return "/scratch/midway2/chard/clean/packmol-17.163/packmol < %s" % input_file
    

@bash_app(executors=['midway_cpu'])
def tleap(input_file, inputs=[], outputs=[], stdout=None, stderr=None, mock=False ):
    return '''module load amber/16; tleap -f %s''' % input_file

## Minimization
@bash_app(executors=['midway_cpu'])
def minimization(input_file, inputs=[], outputs=[], stdout=None, stderr=None, mock=False ): 
    return '''module load amber; 
              mpirun -np 6 pmemd.MPI -O -i %s -o {outputs[0]} -c {inputs[0]} -p {inputs[1]} -r {outputs[1]} -ref {inputs[0]}''' % input_file

## Used for Heating, Equilibration, and Production
@bash_app(executors=['midway_gpu'])
def pmemd_cuda(input_file=None, inputs=[], outputs=[], stdout=None, stderr=None, ref=True ): 
    if ref:
        r = "-ref {inputs[1]}"
    else: 
        r = ""

    return '''pmemd.cuda -O -i %s -p {inputs[0]}  -c {inputs[1]} -o {outputs[0]} -r {outputs[1]} -x {outputs[2]} %s''' % (input_file, r)
       
# Helper functions
def round_temp(x, base=5):
    return "{0:.2f}".format(int(base * round(float(x)/base)))

def get_temp(heat_file):
    with open(heat_file, 'r') as f:
        for l in f:
            if l.startswith(" NSTEP"):
                t = l
            if "A V E R A G E" in l:
                break

    start_pos = t.find("TEMP(K)")
    end_pos = t.find("PRESS")

    temp = float(t[start_pos + 9 : end_pos])

    return round_temp(temp)


################
## Workflow code
################


print ("Running packmol")
input_path = create_template("packmol.jinja2", "packmol.inp", {"STRUCTURE_1": "il2y-f0-charge-balanced.pdb", "STRUCTURE_2": "DCA.pdb",
                    "STRUCTURE_3": "TMA.pdb", "STRUCTURE_4": "wat.pdb", "OUTPUT": "protein-w-IL.pdb"})

packmol_future = packmol(input_file=input_path, 
                inputs = ["il2y-f0-charge-balanced.pdb", "DCA.pdb", "TMA.pdb", "wat.pdb"], 
                outputs = [ "protein-w-IL.pdb"], 
                stdout="packmol.stdout", stderr="packmol.stderr",)


print ("Running tleap")
input_path = create_template("tleap.jinja2", "tleap.inp", {"DCA": "DCA.mol2", "TMA": "TMA.mol2",
            "DCA_FRCMOD": "DCA.frcmod", "TMA_FRCMOD": "TMA.frcmod",
            "DCA_LIB": "DCA.lib", "TMA_LIB": "TMA.lib",
            "PROT": packmol_future.outputs[0].result(),
            "OUTPUT1": "prot.parm7", "OUTPUT2": "prot.rst7"})

tleap_future = tleap(
                input_file=input_path,
                inputs=["DCA.mol2", "TMA.mol2", "DCA.frcmod", "TMA.frcmod", "DCA.lib", "TMA.lib", packmol_future.outputs[0]],
                outputs=["prot.parm7", "prot.rst7"],
                stdout="tleap.stdout", stderr="tleap.stderr")


print("Running minimization")
min_future = minimization(
                input_file="config_files/static/min.in",
                inputs= [tleap_future.outputs[1], tleap_future.outputs[0]],
                outputs = ['min.out', 'min.rst7'],
                stdout="min.stdout", stderr="min.stderr")

#min_future.result()

print("Constant volume heating")
heat_vol = pmemd_cuda(input_file='config_files/static/heat.1.in', 
        inputs=[tleap_future.outputs[0], min_future.outputs[1]], 
        outputs=['heat.1.out', 'heat.1.rst7', 'heat.1.nc'], 
        stdout='heat1.stdout', stderr='heat1.stderr')

heat_vol.result()

print("Constant pressure heat")
count = 2
heated = [None, heat_vol]

while True:
    # Note: hack here as it throws an error getting outputs[0].result
    if count == 2:
        temp = 100
    else:
        temp = get_temp(heated[count-1].outputs[0].filepath)
    
    print("Heat %s -- Temp %s" % (count, temp))
    input_file = create_template("heat.jinja2", "heat.%s.in" % count, { "VALUE_1" : temp})
    try:
        heated.append(pmemd_cuda(input_file=input_file, 
            inputs=[tleap_future.outputs[0], heated[count-1].outputs[1].filepath], 
            outputs=['heat.%s.out' % count, 'heat.%s.rst7' % count, 'heat.%s.nc' % count], 
            stdout='heat.%s.stdout' % count, stderr='heat.%s.stderr' % count))
        heated[count].result()
        print("No error: heating finished")
        break
    except:
        print ("Error heating, restarting from last temperature")
    count +=1


print("Running Equilibration")
equib = [None]# 
count = 1

for i in [10,5, 2.5, 1.25, 0.5, 0.1]:
    if count == 1:
        in_rst = heated[-1].outputs[1]
    else:
        in_rst = equib[count - 1].outputs[1]

    input_file = create_template("equilibration.jinja2", "equib.%s.in" % count, { "RESTRAINT_WT" : 10})
    
    equib.append(pmemd_cuda(input_file=input_file, 
        inputs=[tleap_future.outputs[0], in_rst], 
        outputs=['equib.%s.out' % count, 'equib.%s.rst7' % count, 'equib.%s.nc' % count], 
        stdout='equib.%s.stdout' % count, stderr='equib.%s.stderr' % count, 
        ref=True))
    equib[count].result()
    count +=1

print ("Running Production")
prod = [None]

for i in range (1,4):
    if i == 1:
        in_rst = equib[-1].outputs[1]
    else:
        in_rst = prod[i-1].outputs[1]   

    prod.append(pmemd_cuda(input_file="config_files/static/prod.in",
        inputs=[tleap_future.outputs[0], in_rst],
        outputs=['prod.%s.out' % (i), 'prod.%s.rst7' % (i), 'prod.%s.nc' % (i)],
        stdout='prod.%s.stdout' % (i), stderr='prod.%s.stderr' % (i),
        ref=False))
    prod[i].result()


print ("Finished")
