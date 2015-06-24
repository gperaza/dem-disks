import numpy as np
import os
import subprocess
import csv
#import sys
import itertools
import multiprocessing

'''This script lunches a sweep of simulations exploring the behavior
of the 3 disk setup. The script relies on the C program disks to run
the simulations. The script contains some functions to analyze data
obtained from this simulations'''

def write_input_file(ac,tilt,dirname):
    f = open(dirname + "/input_file", "w")
    f.write("""#nParticles       1
#seed             123456   (seed for random generator)
#box_w            2        (in disk units)
#box_h            3        (in disk units)
#freq             80       (in Hz, for bottom movement)
#dimensionlessAc  """+ac+"""      (dimensionless acceleration of bottom)
#gravity          9.81     (m/s)
#gravityAngle     """+tilt+"""      (in fractions of PI)
#bGamma           0        (bulk dissipation)
#timestep         1e-6     (in s, timestep for integrator)
#relaxTime        5        (in s, time for relaxation)
#thermalTime      5
#runTime          5     (in s, time for simulation)
#timeForGraph     1e0      (in s, time between graphics)
#timeForWrite     1e-6     (in s, time between writes)
#meanR            0.02     (in m, mean disk radius)
#density          3.57     (density of the material)
#kn               4.5e6    (normal elastic constant)
#pr               0.37     (poisson ratio of the material)
#mu               0.1      (friction coefficient)
#vGamma           150      (viscoelastic dissipation)
#bCondType        1        (1: sinusoidal 2:random vib)
"""
)

def run_simulation((ac, tilt)):
    reRun = False
    #execPath = "/home/moukarzel/DEM-source/disks"
    execPath = "/home/gperaza/Documents/RotatingDisks/DEM-source/disks"
    dirname = "tilt."+tilt+".ac."+ac
    if os.path.isdir(dirname):
        if reRun:
            pass
        else:
            print "Skipping ac = ", ac, "and tilt = ", tilt, \
                  ". Already been simulated."
            return
    else:
        os.mkdir(dirname)

    write_input_file(ac,tilt,dirname)
    log = open(dirname+"/log","w")
    errorLog = open(dirname+"/error_log","w")
    print "running ac = ", ac, "tilt = ", tilt
    subprocess.check_call([execPath], cwd = dirname,
                          stdout = log, stderr = errorLog)
    log.close()
    errorLog.close()
    print "done"

def main(exType):
    if not os.path.isdir("3diskSimSweep"):
        os.mkdir("3diskSimSweep")
    os.chdir("3diskSimSweep")
    #acList = ["%06.4f"%x for x in np.arange(0, 1.25, 0.005)]
    acList = ["1.0000"]
    #tiltList = ["%05.3f"%x for x in np.arange(0, 0.166, 0.001)]
    #tiltList = ["%05.3f"%x for x in np.arange(0, 0.071, 0.001)]
    tiltList = ["0.000"]
    #vibList2 = [1e-6, 3e-6, 5e-6, 8e-6, 1e-5, 3e-5, 5e-5, 8e-5, 1e-4, 1.5e-4, 2e-4]
    #vibList = ["%09.7f"%x for x in vibList2]
    #reRun = False
    #if reRun:
    #    print "reRun is set to TRUE, old simulation will be lost.Continue?[y/n]"
    pool_size = 1
    pool = multiprocessing.Pool(processes=pool_size)
    if exType == 1: pool_outputs = pool.map(run_simulation, itertools.product(acList, tiltList))
    if exType == 2: pool_outputs = pool.map(run_simulation, itertools.product(vibList, tiltList))
    pool.close() # no more tasks
    pool.join()  # wrap up current tasks
    os.chdir("..")

if __name__ == "__main__":
    exType = 1 #Remember to change it in the input file
    main(exType)
