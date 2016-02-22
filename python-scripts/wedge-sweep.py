import os
import subprocess
import itertools
import multiprocessing
import numpy as np


def write_input_file(ac, tilt, dirname):
    f = open(dirname + "/input_file", "w")
    f.write("""#nParticles           1
#seed                 100000   (seed for random generator)
#box_w                2        (in disk units)
#box_h                3        (in disk units)
#freq                 80       (in Hz, for bottom movement)
#dimensionlessAc      """+ac+"""      (dimensionless acceleration of bottom)
#gravity              9.81     (m/s)
#gravityAngle         """+tilt+"""      (in fractions of PI)
#bGamma               0        (bulk dissipation)
#timestep             1e-5     (in s, timestep for integrator)
#relaxTime            5       (in s, time for relaxation)
#thermalTime          10
#runTime              100     (in s, time for simulation)
#timeForGraph         1e9      (in s, time between graphics)
#timeForWriteRun      1e0      (in s, time between writes)
#timeForWriteThermal  1e9      (in s, time between writes)
#meanR                0.02     (in m, mean disk radius)
#density              3.57     (density of the material)
#kn                   4.5e4    (normal elastic constant)
#pr                   0.37     (poisson ratio of the material)
#mu                   0.1      (friction coefficient)
#vGamma               0.200    (viscoelastic dissipation)
#bCondType            1        (1: sinusoidal 2:random vib)
#relInitDisp          1        (between 1 and -1)
#wedge                1
""")


def run_simulation(ac, tilt):
    reRun = False
    execPath = "/home/gperaza/Documents/RotatingDisks/Code_and_scripts/DEM-source/disks"
    dirname = "tilt."+tilt+".ac."+ac
    if os.path.isdir(dirname):
        if not reRun:
            print("Skipping ac = ", ac, "and tilt = ", tilt,
                  ". Already been simulated.")
            return
    else:
        os.mkdir(dirname)
        os.mkdir(dirname+"/Collisions")

    write_input_file(ac, tilt, dirname)
    print("running ac = ", ac, "tilt = ", tilt)
    subprocess.check_call([execPath], cwd=dirname)
    print("done ac = ", ac, "tilt = ", tilt)


def main():
    # List of acceleration amplitudes
    acList = ["0.5", "0.6", "0.7", "0.8", "0.9", "1.0",
              "1.2", "1.4", "1.6", "1.8", "2.0"]
    # List of tilt angles.
    tiltList = ["%04.2f" % x for x in np.arange(0.01, 0.166, 0.02)]

    for tilt, ac in itertools.product(tiltList, acList):
        run_simulation(ac, tilt)

if __name__ == "__main__":
    main()
