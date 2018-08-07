import parsl
import argparse

from parsl import App
from parsl.app.app import bash_app, python_app
from parsl.configs.local_threads import config

parsl.load(config)

#@App("bash") 
@bash_app
def setupVisit(DC2_CONFIGDIR, mock=True):
    """Setup phoSim inputs and prepare for parallelization
    """

    cmd = """export DC2_CONFIGDIR={0};
{0}/setupVisit.sh
sleep 1
""".format(DC2_CONFIGDIR)
    if mock:
        return "echo \"{}\"".format(cmd)
    else:
        return cmd


@bash_app
def trim(DC2_CONFIGDIR, trim_index=None, envs={}, mock=True):
    """Trim instance catalog for this visit
    """

    cmd = """export DC2_CONFIGDIR={0};
{0}/runTrim.sh {1}
sleep 1
""".format(DC2_CONFIGDIR, trim_index)
    if mock:
        return "echo \"{}\"".format(cmd)
    else:
        return cmd


@bash_app
def rayTrace(DC2_CONFIGDIR, trim, mock=True):
    """Trim instance catalog for this visit
    """

    cmd = """export DC2_CONFIGDIR={0};
{0}/runRaytrace.sh
sleep 0
""".format(DC2_CONFIGDIR)
    if mock:
        return "echo \"{}\"".format(cmd)
    else:
        return cmd

@bash_app
def e2adc(DC2_CONFIGDIR, trim, mock=True):
    """The e2adc processing step
    """

    cmd = """export DC2_CONFIGDIR={0};
{0}/runE2ADC.sh
sleep 0
""".format(DC2_CONFIGDIR)
    if mock:
        return "echo \"{}\"".format(cmd)
    else:
        return cmd

@python_app
def reg_files(e2adc, mock=True):

    if mock:
        print("Reg_files mocking")
        return True
    

def launchTrim(DC2_CONFIGDIR, setup_fu, mock=True):
    """Creates 22 trim sub-streams
    """

    trim_futs = []
    for i in range(22):
        trim_fu = trim(DC2_CONFIGDIR, trim_index=i, envs={}, mock=mock)
        trim_futs.append(trim_fu)

    return trim_futs

def launchSensors(DC2_CONFIGDIR, trim_futs, mock=True):
    """Launches raytrace on the list of trim_futs
    """
    raytrace_futs = {}
    for trim_fut in trim_futs:
        sensorlist = []
        for i in range(190):
            raytrace_fu = rayTrace(DC2_CONFIGDIR, trim_fut, mock=mock)
            e2adc_fu = e2adc(DC2_CONFIGDIR, raytrace_fu, mock=mock)
            reg_fu = reg_files(e2adc_fu, mock=mock)
            print(reg_fu.result())
            sensorlist.append(reg_fu)
        raytrace_futs[trim_fu] = sensorlist
    
    return raytrace_futs



if __name__ == "__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--width", default="10",
                        help="width of the pipeline")
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Count of apps to launch")
    parser.add_argument("-m", "--mock", action='store_true',
                        help="Enable mock run")

    args = parser.parse_args()

    DC2_FILTER = "r"
    DC2_FIELD = "WFD"
    WORKFLOW = "DC2-R2-0p"
    TASKNAME = "{}-test-{}".format(WORKFLOW, DC2_FILTER)

    DC2_PATH = "/global/projecta/projectdirs/lsst/production/DC2"
    DC2_WORKFLOW = "{}/{}".format(DC2_PATH, WORKFLOW)
    DC2_ROOT = "{}/{}".format(DC2_PATH, TASKNAME)
    DC2_CONFIGDIR = "{}/config".format(DC2_WORKFLOW)
    DC2_OUTPUT = "{}/output".format(DC2_ROOT)
    logRoot = "{}/logs".format(DC2_ROOT)

    
    setup_fu = setupVisit(DC2_CONFIGDIR, mock=args.mock)
    trim_futs = launchTrim(DC2_CONFIGDIR, setup_fu, mock=args.mock)
    sensor_futs = launchSensors(DC2_CONFIGDIR, trim_futs, mock=args.mock)
    print([i.result() for i in launch_trims])
    for trim in sensor_futs:
        print([i.result() for i in sensor_futs[trim]])
        
