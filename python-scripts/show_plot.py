import plots_3disks_lib as p3d
import sys
import matplotlib.pyplot as plt

def main(args):
    assert(len(args) == 3 or len(args) == 2)
    plotType = args[1]
    value = args[2]
    if plotType == "ac":
        print "plotting for ac", value
        fig, axes = plt.subplots(2)
        p3d.plot_vel_tilt_for_ac(value, ax=axes[0])
        p3d.plot_links_tilt_for_ac(value, ax=axes[1])
    elif plotType == "tilt":
        print "plotting for tilt", value
        fig, axes = plt.subplots(2)
        p3d.plot_vel_ac_for_tilt(value, ax=axes[0])
        p3d.plot_links_ac_for_tilt(value, ax=axes[1])
    elif plotType == "angle":
        print "plotting angle vs time"
        p3d.plot_angle_time(directory=value, show=True, plotRelax=True)
    else:
        print "Provide argumets as 'plotType value'"
    plt.show()
if __name__ == "__main__":
    main(sys.argv)
