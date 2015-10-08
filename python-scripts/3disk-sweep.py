import os
import subprocess
import itertools
import multiprocessing
import numpy as np

'''This script lunches a sweep of simulations exploring the behavior
of the 3 disk setup.'''


def write_input_file(ac, tilt, dirname):
    f = open(dirname + "/input_file", "w")
    f.write("""#nParticles           1
#seed                 123456   (seed for random generator)
#box_w                2        (in disk units)
#box_h                3        (in disk units)
#freq                 80       (in Hz, for bottom movement)
#dimensionlessAc      """+ac+"""      (dimensionless acceleration of bottom)
#gravity              9.81     (m/s)
#gravityAngle         """+tilt+"""      (in fractions of PI)
#bGamma               0        (bulk dissipation)
#timestep             1e-6     (in s, timestep for integrator)
#relaxTime            10       (in s, time for relaxation)
#thermalTime          100
#runTime              1000     (in s, time for simulation)
#timeForGraph         1e0      (in s, time between graphics)
#timeForWriteRun      1e0      (in s, time between writes)
#timeForWriteThermal  1e0      (in s, time between writes)
#meanR                0.02     (in m, mean disk radius)
#density              3.57     (density of the material)
#kn                   4.5e6    (normal elastic constant)
#pr                   0.37     (poisson ratio of the material)
#mu                   0.1      (friction coefficient)
#vGamma               150      (viscoelastic dissipation)
#bCondType            1        (1: sinusoidal 2:random vib)
#relInitDisp          1        (between 1 and -1)
#wedge                1
""")


def run_simulation((ac, tilt)):
    reRun = False
    execPath = "/home/moukarzel/dem-code/disks"
    dirname = "tilt."+tilt+".ac."+ac
    if os.path.isdir(dirname):
        if not reRun:
            print "Skipping ac = ", ac, "and tilt = ", tilt, \
                ". Already been simulated."
            return
    else:
        os.mkdir(dirname)

    write_input_file(ac, tilt, dirname)
    log = open(dirname + "/log", "w")
    errorLog = open(dirname + "/error_log", "w")
    print "running ac = ", ac, "tilt = ", tilt
    subprocess.check_call([execPath], cwd=dirname, stdout=log, stderr=errorLog)
    log.close()
    errorLog.close()
    print "done ac = ", ac, "tilt = ", tilt


def main():
    if not os.path.isdir("3diskSimSweep"):
        os.mkdir("3diskSimSweep")
    os.chdir("3diskSimSweep")

    # List of acceleration amplitudes
    acList1 = ["1.0000"]
    acList2 = ["%06.4f"%x for x in np.arange(0, 1.25, 0.05)]
    acList3 = ["%06.4f"%x for x in np.arange(0, 1.25, 0.005)]

    # List of tilt angles.
    tiltList1 = ["0.000"]
    tiltList2 = ["%05.3f"%x for x in np.arange(0, 0.166, 0.01)]
    tiltList3 = ["%05.3f"%x for x in np.arange(0, 0.166, 0.001)]
    tiltList4 = ["%05.3f"%x for x in np.arange(0, 0.071, 0.001)]

    pool_size = 16
    pool = multiprocessing.Pool(processes=pool_size)
    pool.map(run_simulation, itertools.product(acList2, tiltList3))
    pool.close()  # no more tasks
    pool.join()  # wrap up current tasks
    os.chdir("..")

if __name__ == "__main__":
    main()
