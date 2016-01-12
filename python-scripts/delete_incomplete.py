import os
import subprocess
import numpy as np
import shutil


def main():
    
    # List of acceleration amplitudes
    acList1 = ["1.0"]
    acList2 = ["%06.4f"%x for x in np.arange(0.8, 1.25, 0.05)]
    acList3 = ["%06.4f"%x for x in np.arange(0, 1.25, 0.005)]
    acList4 = ["0.8", "0.9", "1.0", "1.1", "1.2", "1.4", "1.6", "1.8", "2.0", "3.0", "4.0", "5.0"]

    # List of tilt angles.
    tiltList1 = ["0.000"]
    tiltList2 = ["%05.3f"%x for x in np.arange(0, 0.166, 0.01)]
    tiltList3 = ["%05.3f"%x for x in np.arange(0, 0.166, 0.001)]
    tiltList4 = ["%05.3f"%x for x in np.arange(0, 0.071, 0.001)]
    tiltList5 = ["%04.2f"%x for x in np.arange(0, 0.166, 0.01)]

    for ac in acList4:
        for tilt in tiltList5:
            dirname = "tilt."+tilt+".ac."+ac
            fname = dirname + "/ending_phase_space.out"
            if os.path.isdir(dirname):
                if not os.path.exists(fname):
                    print "directory empty"
                else:
                    if os.stat(fname).st_size == 0:
                        print "removing", fname
                        shutil.move(dirname, "Incomplete")

if __name__ == "__main__":
    main()
