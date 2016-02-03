import os
import subprocess
import itertools
import multiprocessing
import numpy as np


def run_simulation((ac, seed)):
    reRun = True
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


def main():
    sweepdir = "DEM-Sweep"
    if not os.path.isdir(sweepdir):
        os.mkdir(sweepdir)
    os.chdir(sweepdir)

    # List of acceleration amplitudes
    acList = ["00.1", "00.5", "00.9", "01.0", "01.2", "01.5", "01.7", "02.0", 
              "03.0", "04.0", "05.0", "06.0", "07.0", "08.0", "09.0", "10.0"]

    # List of seeds.
    # random.seed(123456)
    # seedList = [str(random.randint(100000, 900000)) for x in range(10)]
    seedList = ["100000"]

    pool_size = 16
    pool = multiprocessing.Pool(processes=pool_size)
    pool.map(run_simulation, itertools.product(acList, seedList))
    pool.close()  # no more tasks
    pool.join()  # wrap up current tasks
    
    os.chdir("..")

if __name__ == "__main__":
    main()
