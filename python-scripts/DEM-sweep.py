import os
import subprocess
import itertools
import multiprocessing
import numpy as np


def write_input_file(ac, seed, gamma, mu, dirname):
    #seed = "123456"
    f = open(dirname + "/input_file", "w")
    f.write("""#nParticles           60
#seed                 """+seed+"""   (seed for random generator)
#box_w                10        (in disk units)
#box_h                10        (in disk units)
#freq                 80       (in Hz, for bottom movement)
#dimensionlessAc      """+ac+"""      (dimensionless acceleration of bottom)
#gravity              9.81     (m/s)
#gravityAngle         0.0      (in fractions of PI)
#bGamma               0        (bulk dissipation)
#relaxTime            10       (in s, time for relaxation)
#thermalTime          50
#runTime              500      (in s, time for simulation)
#timeForGraph         1e0      (in s, time between graphics)
#timeForWriteRun      1e0      (in s, time between writes)
#timeForWriteThermal  1e9      (in s, time between writes)
#meanR                0.02     (in m, mean disk radius)
#density              3.57     (density of the material)
#kn                   1e4      (normal elastic constant)
#pr                   0.773    (elastic constant ratio)
#mu                   """+mu+"""      (friction coefficient)
#vGamma               """+gamma+"""      (viscoelastic dissipation)
#bCondType            1        (1: sinusoidal 2:random vib)
#relInitDisp          1        (between 1 and -1)
#wedge                0
""")


def run_simulation((ac, seed, gamma, mu)):
    reRun = False
    execPath = "/home/moukarzel/dem-code/disks"
    dirname = "seed."+seed+".ac."+ac
    if os.path.isdir(dirname):
        if not reRun:
            print("Skipping ac = ", ac, "with seed = ", seed,
                  ". Already been simulated.")
            return
    else:
        os.mkdir(dirname)
        os.mkdir(dirname+"/Collisions")

    write_input_file(ac, seed, gamma, mu, dirname)
    log = open(dirname + "/log", "w")
    errorLog = open(dirname + "/error_log", "w")
    error = False
    if subprocess.call([execPath], cwd=dirname, stdout=log, stderr=errorLog):
        error = True
        print("Error at seed = ", seed, "ac = ", ac)
    log.close()
    errorLog.close()
    if not error:
        print("done ac = ", ac, "seed = ", seed)


def main(gamma, mu):
    sweepdir = "DEM-Sweep-k-1e4-vg-"+gamma+"-mu-"+mu
    if not os.path.isdir(sweepdir):
        os.mkdir(sweepdir)
    os.chdir(sweepdir)

    # List of acceleration amplitudes
    acList = ["00.1", "00.3", "00.5", "00.7", "00.9", "01.0", "01.2", "01.5", 
              "02.0", "03.0", "04.0", "05.0", "06.0", "07.0", "08.0", "09.0"]

    # List of seeds.
    seedList = ["100000", "100001", "100002", "100003", "100004", "100005", "100006", "100007", "100008", "100009"]

    pool_size = 16
    pool = multiprocessing.Pool(processes=pool_size)
    pool.map(run_simulation, itertools.product(acList, seedList, [gamma], [mu]))
    pool.close()  # no more tasks
    pool.join()  # wrap up current tasks
    
    os.chdir("..")

if __name__ == "__main__":
    gammaList = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]
    muList = ["0.01", "0.03", "0.05", "0.07", "0.09", "0.10", "0.30", "0.50", "0.70", "0.90"]
    for mu in muList:
        main("0.1", mu)
