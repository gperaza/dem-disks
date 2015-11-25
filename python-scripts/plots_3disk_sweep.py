import matplotlib.pyplot as plt
import os
import plots_3disks_lib as p3d


def get_tilts():
    tilts = []
    directoryList = [d for d in os.listdir(".")
                     if "tilt" in d and os.path.isdir(d)]
    for d in directoryList:
        dList = d.split(".")
        index = dList.index("tilt")
        tilt = dList[index+1]+"."+dList[index+2]
        if tilt not in tilts:
            tilts.append(tilt)
    tilts.sort()
    return tilts


def get_acs():
    acs = []
    directoryList = [d for d in os.listdir(".")
                     if "ac" in d and os.path.isdir(d)]
    for d in directoryList:
        dList = d.split(".")
        index = dList.index("ac")
        ac = dList[index+1]+"."+dList[index+2]
        if ac not in acs:
            acs.append(ac)
    acs.sort()
    return acs


def plot_angle_time_sweep():
    directoryList = [d for d in os.listdir(".")
                     if "ac" in d and os.path.isdir(d)]
    directoryList.sort()
    for d in directoryList:
        p3d.plot_angle_time(d, plotRelax=False)


def plot_vel_ac_for_tilts():
    tiltList = get_tilts()
    acList = get_acs()
    if len(acList) <= 1:
        return
    figAll = plt.figure(num=1, figsize=(12, 10))
    axAll = figAll.add_subplot(111)
    axAll.set_title("Rotational velocity vs acceleration")
    axAll.set_xlabel("Dimensionless acceleration")
    axAll.set_ylabel("Rotational velocity(rpm)")
    for tilt in tiltList:
        figOne = plt.figure(num=2, figsize=(12, 10))
        axVel = figOne.add_subplot(211)
        axVel.set_title("Rotational velocity vs acceleration")
        axVel.set_xlabel("Dimensionless acceleration")
        axVel.set_ylabel("Rotational velocity(rpm)")
        axLinks = figOne.add_subplot(212)
        axLinks.set_title("Probability of links states vs acceleration.")
        axLinks.set_xlabel("Dimensionless acceleration")
        axLinks.set_ylabel("State probability")
        p3d.plot_vel_ac_for_tilt(tilt, ax=axAll)
        p3d.plot_vel_ac_for_tilt(tilt, ax=axVel)
        p3d.plot_links_ac_for_tilt(tilt, ax=axLinks)
        axLinks.legend(loc="best")
        axVel.legend(loc="best")
        figOne.savefig("plot_vel_ac_for_tilt_"+tilt+".png", format="png")
        plt.close(figOne)
    # axAll.legend(loc="best")
    figAll.savefig("plot_vel_ac_all.png", format="png")
    plt.close(figAll)


def plot_vel_tilt_for_acs():
    acList = get_acs()
    tiltList = get_tilts()
    if len(tiltList) <= 1:
        return

    figAll = plt.figure(num=1, figsize=(12, 10))
    axAll = figAll.add_subplot(111)
    axAll.set_title("Rotational velocity vs tilt")
    axAll.set_xlabel("tilt (rho/pi)")
    axAll.set_ylabel("Rotational velocity(rpm)/ac^2")

    for ac in acList:
        if float(ac) == 0:
            continue
        p3d.plot_vel_tilt_for_ac(ac, ax=axAll, scale=1/float(ac)**2)

        figOne = plt.figure(num=2, figsize=(12, 10))

        axVel = figOne.add_subplot(211)
        axVel.set_title("Rotational velocity vs tilt")
        axVel.set_xlabel("tilt (rho/pi)")
        axVel.set_ylabel("Rotational velocity(rpm)")
        p3d.plot_vel_tilt_for_ac(ac, ax=axVel)
        axVel.legend(loc="best")

        axLinks = figOne.add_subplot(212)
        axLinks.set_title("Probability of link states vs acceleration.")
        axLinks.set_xlabel("tilt (rho/pi)")
        axLinks.set_ylabel("State probability")
        p3d.plot_links_tilt_for_ac(ac, ax=axLinks)
        axLinks.legend(loc="best")

        figOne.savefig("plot_vel_tilt_for_ac_"+ac+".png", format="png")
        plt.close(figOne)
    axAll.legend(loc="best")
    figAll.savefig("plot_vel_tilt_all.png", format="png")
    plt.close(figAll)


def main():
    # print("Generating all angular paths")
    # plot_angle_time_sweep()
    # print("Plotting vel vs ac")
    # plot_vel_ac_for_tilts()
    print("Plotting vel vs tilt")
    plot_vel_tilt_for_acs()
    print("Done")

if __name__ == "__main__":
    main()
