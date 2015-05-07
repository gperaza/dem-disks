import matplotlib.pyplot as plt
import csv, os, numpy as np
from matplotlib.collections import LineCollection
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec

def autolabel(rects):
    """Attach some text labels. To be used for histograms."""
    for rect in rects:
        height = rect.get_height()
        plt.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%.3f'%height,
                 ha='center', va='bottom')

def save_gnuplot(x, y, name):
    """Creates a gnuplot friendly file of y vs x."""
    dFile = open(name, "w")
    writer = csv.writer(dFile, delimiter = " ")
    for i in range(len(x)):
        writer.writerow([x[i], y[i]])
    dFile.close()

def disk_fall(directory="."):
    fall = False
    dataFile = directory + "/ending_phase_space.out"
    f = open(dataFile,"r")
    reader = csv.reader(f, delimiter = " ")
    for row in reader:
        if row[0] == "0":
            xi = float(row[26])
            xf = float(row[5])
            R = float(row[2])
    if abs(xf-xi) > R:
        fall = True
    f.close()
    return fall

def plot_coord_time(coord, directory=".", show=False, iTime=-float("inf"),
                    fTime = float("inf"), diskId="0", ax=None, color="b",
                    legend=True):
    """Plots any of the coordinates (x,y,w) of the upper disk as a
    function of time for the simulation inside <directory>. If show ==
    True then the plot is displayed instead of saved.

    """
    dataFile = directory + "/phase_space.out"
    x0 = []
    y0 = []
    w0 = []
    time = []
    readerFile = open(dataFile, "r")
    reader = csv.reader(readerFile, delimiter = " ")
    for row in reader:
        if row[0].startswith("#time:"):
            currTime = float(row[0].split(":")[1])
            if iTime <= currTime <= fTime:
                time.append(currTime)
        if row[0] == diskId:
            if iTime <= currTime <= fTime:
                x0.append(float(row[5]))
                y0.append(float(row[6]))
                w0.append(float(row[7]))
    readerFile.close()
    time = np.array(time)
    x0 = np.array(x0)
    deltaX = abs(max(x0) - min(x0))
    y0 = np.array(y0)
    deltaY = abs(max(y0) - min(y0))
    w0 = np.array(w0)
    if ax is not None:
        ax.set_xlabel("time(s)")
        ax.set_ylabel(coord + " of disk " + diskId)
        if coord == "x":
            ax.plot(time, x0, label="Delta=%e"%deltaX, color=color)
            if legend: ax.legend()
        elif coord == "y":
            ax.plot(time, y0, label="Delta=%e"%deltaY, color=color)
            if legend: ax.legend()
        elif coord == "w":
            ax.plot(time, w0, color=color)
        #plt.tight_layout()
        return
    plt.figure(num=1, figsize=(8,6))
    plt.xlabel("time(s)")
    plt.ylabel(coord)
    if coord == "x":
        plt.plot(time, x0, label="Delta=%lf" % (deltaX))
        plt.legend()
    elif coord == "y":
        plt.plot(time, y0, label="Delta=%lf" % (deltaY))
        plt.leyend()
    elif coord == "w":
        plt.plot(time, w0)
        plt.title(coord + " coordinate vs time")
    plt.tight_layout()
    if show:
        plt.show()
    else:
        plt.savefig(directory + "/plot_" + coord + "_time.png", format="png")
    plt.close()

def plot_angle_time(directory=".", show=False, plotRelax=False):
    """Plots the angle of the upper disk as a function of time for the
    simulation inside <directory>. If show == True then the plot is
    displayed instead of saved.

    """
    dataFile = directory + "/phase_space.out"
    angle = []
    time = []
    readerFile = open(dataFile, "r")
    reader = csv.reader(readerFile, delimiter = " ")
    for row in reader:
        if row[0].startswith("#time:"):
            currTime = float(row[0].split(":")[1])
            if plotRelax:
                time.append(currTime)
            elif currTime >= 0:
                time.append(currTime)
        if row[0] == "0": #0 is the id of the upper disk
            theta = float(row[7]) #this is the angle at this time
            theta0 = float(row[28])
            if plotRelax:
                angle.append(theta)
            elif currTime >= 0:
                angle.append(theta - theta0)
    readerFile.close()
    time = np.array(time)
    angle = np.array(angle)
    plt.figure(num=1, figsize=(8,6))
    plt.title("Rotated angle vs time")
    plt.xlabel("time(s)")
    plt.ylabel("angle(radians)")
    plt.plot(time, angle)
    if show:
        plt.show()
    else:
        plt.savefig(directory + "/plot_angle_time.png", format="png")
    plt.close()

def plot_center_path(directory=".", show=False, iTime=0, fTime=float("inf"),
                     ax=None):
    """Plots the path (y vs x) of the disk center of mass.  Color the
    point in the path reflects the disk's angular velocity.

    """
    dataFile = directory + "/phase_space.out"
    x0 = []
    y0 = []
    w1 = []
    readerFile = open(dataFile, "r")
    reader = csv.reader(readerFile, delimiter = " ")
    for row in reader:
        if row[0].startswith("#time:"):
            currTime = float(row[0].split(":")[1])
        if row[0] == "0" and iTime <= currTime <= fTime:
            x0.append(float(row[5]))
            y0.append(float(row[6]))
            w1.append(float(row[10]))
    readerFile.close()
    x0 = np.array(x0)
    y0 = np.array(y0)
    w1 = np.array(w1)
    points = np.array([x0, y0]).T.reshape(-1,1,2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=plt.get_cmap('coolwarm'))
    lc.set_array(w1)
    lc.set_linewidth(5)
    if ax is not None:
        ax.add_collection(lc)
    else:
        fig = plt.figure(figsize=(12,8))
        plt.gca().add_collection(lc)
        axcb = fig.colorbar(lc)
        axcb.set_label('angular velocity')
    plt.title("2D path of center of mass")
    plt.xlabel("x")
    plt.ylabel("y")
    #plt.axes().set_aspect('equal', 'datalim')
    plt.xlim(min(x0),max(x0))
    plt.plot(x0, y0)
    if ax is not None: return
    if show:
        plt.show()
    else:
        plt.savefig(directory + "/plot_path.png", format="png")
    plt.close(fig)

def plot_x_time(directory=".", show=False, plotRelax = False):
    """Plots the x coordinate of the upper disk as a function of time for
    the simulation inside <directory>. If show == True then the plot
    is displayed instead of saved.

    """
    dataFile = directory + "/phase_space.out"
    x0 = []
    time = []
    readerFile = open(dataFile, "r")
    reader = csv.reader(readerFile, delimiter = " ")
    for row in reader:
        if row[0].startswith("#time:"):
            currTime = float(row[0].split(":")[1])
            if plotRelax:
                time.append(currTime)
            elif currTime >= 0:
                time.append(currTime)
        if row[0] == "0": #0 is the id of the upper disk
            if plotRelax:
                x0.append(float(row[5]))
            elif currTime >= 0:
                x0.append(float(row[5]))
    readerFile.close()
    time = np.array(time)
    x0 = np.array(x0)
    plt.figure(num=1, figsize=(8,6))
    plt.title("x coordinate vs time")
    plt.xlabel("time(s)")
    plt.ylabel("x(m)")
    plt.plot(time, x0)
    if show:
        plt.show()
    else:
        plt.savefig(directory + "/plot_x_time.png", format="png")
    plt.close()

def plot_y_time(directory=".", show=False, plotRelax = False, diskId="0",
                ax=None):
    """Plots the y coordinate of the upper disk as a function of time for
    the simulation inside <directory>. If show == True then the plot
    is displayed instead of saved.

    """
    dataFile = directory + "/phase_space.out"
    y0 = []
    time = []
    readerFile = open(dataFile, "r")
    reader = csv.reader(readerFile, delimiter = " ")
    for row in reader:
        if row[0].startswith("#time:"):
            currTime = float(row[0].split(":")[1])
            if plotRelax:
                time.append(currTime)
            elif currTime >= 0:
                time.append(currTime)
        if row[0] == diskId:
            if plotRelax:
                y0.append(float(row[6]))
            elif currTime >= 0:
                y0.append(float(row[6]))
    readerFile.close()
    time = np.array(time)
    y0 = np.array(y0)
    if ax is not None:
        ax.set_title("y coordinate of disk " + diskId + " vs time")
        ax.set_xlabel("time(s)")
        ax.set_ylabel("y(m)")
        ax.plot(time, y0)
        return
    plt.figure(num=1, figsize=(8,6))
    plt.title("y coordinate vs time")
    plt.xlabel("time(s)")
    plt.ylabel("y(m)")
    plt.plot(time, y0)
    if show:
        plt.show()
    else:
        plt.savefig(directory + "/plot_y_time.png", format="png")
    plt.close()

def hist_links(directory='.', show=False):
    """Makes histograms from links statistics."""
    dataFile = directory + "/3disk_avglinkstat.out"
    readerFile = open(dataFile, "r")
    reader = csv.reader(readerFile, delimiter = " ")
    for row in reader:
        op_op = float(row[0])
        op_sl = float(row[1])
        op_cl = float(row[2])
        sl_op = float(row[3])
        sl_sl = float(row[4])
        sl_cl = float(row[5])
        cl_op = float(row[6])
        cl_sl = float(row[7])
        cl_cl = float(row[8])
    readerFile.close()
    bothOpen = op_op
    bothClosed = cl_cl
    bothSliding = sl_sl
    open1 = op_sl + op_cl
    closed1 = cl_op + cl_sl
    sliding1 = sl_op + sl_cl
    open2 = sl_op + cl_op
    closed2 = op_cl + sl_cl
    sliding2 = op_sl + cl_sl
    #First histogram, collected probabilities.
    link1 = [open1, closed1, sliding1]
    link2 = [open2, closed2, sliding2]
    both = [bothOpen, bothClosed, bothSliding]
    index = np.arange(3)
    bar_width = 0.3
    opacity = 0.4
    link1 = plt.bar(index, link1, bar_width,
                    alpha=opacity,
                    color='b',
                    label='Link 1')
    link2 = plt.bar(index + bar_width, link2, bar_width,
                    alpha=opacity,
                    color='r',
                    label='Link 2')
    both = plt.bar(index + 2*bar_width, both, bar_width,
                    alpha=opacity,
                    color='g',
                    label='Both')
    plt.ylabel('Time fraction')
    plt.title('Links status')
    plt.xticks(index + 1.5*bar_width,("Open", "Closed", "Sliding"))
    plt.legend()
    plt.tight_layout()
    autolabel(link1)
    autolabel(link2)
    autolabel(both)
    if show:
        plt.show()
    else:
        plt.savefig(directory + "/links1.png", format="png")
    plt.close()
    #Second histogram, individual probabilities.
    links = [op_op, op_sl, op_cl, sl_op, sl_sl, sl_cl, cl_op, cl_sl, cl_cl]
    index = np.arange(len(links))
    bar_width = 0.5
    opacity = 0.4
    links = plt.bar(index-bar_width/2, links, bar_width,
                    alpha=opacity,
                    color='b',
                    label='link1_link2')
    plt.ylabel('Time fraction')
    plt.title('Links status')
    plt.xticks(index, ("op_op", "op_sl", "op_cl",
                                      "sl_op", "sl_sl", "sl_cl",
                                      "cl_op", "cl_sl", "cl_cl"))
    plt.legend()
    plt.tight_layout()
    autolabel(links)
    if show:
        plt.show()
    else:
        plt.savefig(directory + "/links2.png", format="png")
    plt.close()

def plot_vel_ac_for_tilt(tilt, ax=None, show=False):
    """Plots the rotational velocity of the upper disk as a function of
    the input acceleration for a fixed tilt given by <tilt>.  Works
    from within a directory corresponding to a single sweep.  Plots on
    the given ax, if ax=None it generates a standalone plot.

    """
    dataFile = "ending_phase_space.out"
    ac = []
    vel = []
    directoryList = [d for d in os.listdir(".") if "tilt."+tilt in d]
    directoryList.sort()
    for directory in directoryList:
        if disk_fall(directory):
            break
        readerFile = open(directory + "/" + dataFile, "r")
        reader = csv.reader(readerFile, delimiter = " ")
        for row in reader:
            if row[0].startswith("#time:"):
                time = float(row[0].split(":")[1])
            if row[0].startswith("#dimen"):
                ac.append(float(row[0].split(":")[1]))
            if row[0] == "0": #0 is the id of the upper disk
                theta = float(row[7]) #this is the angle at this time
                theta0 = float(row[28])
        #Mean velocity obtained after traversing the whole file.
        vel.append((theta - theta0)/time)
        readerFile.close()
    vel = np.array(vel)
    ac = np.array(ac)
    vel = vel*(60/(2*np.pi)) #change units to rpm
    if ax is not None:
        ax.plot(ac, vel, label=tilt, marker = "o")
        return
    plt.figure(figsize=(8,6))
    plt.title("Rotational velocity vs acceleration for tilt="+tilt)
    plt.xlabel("Dimensionless acceleration")
    plt.ylabel("Rotational velocity(rpm)")
    plt.plot(ac, vel, marker = "o")
    if show:
        plt.show()
    else:
        plt.savefig("plot_vel_ac_for_tilt_"+tilt+".png", format="png")
    plt.close()

def plot_vel_tilt_for_ac(ac, ax=None, show=False):
    """Plots the rotational velocity of the upper disk as a function of
    the tilt for a fixed acceleration given by <ac>.  Works from
    within a directory corresponding to a single sweep.  Plots on the
    given ax, if ax=None it generates a standalone plot.

    """
    dataFile = "ending_phase_space.out"
    tilt = []
    vel = []
    directoryList = [d for d in os.listdir(".") if "ac."+ac in d]
    directoryList.sort()
    for directory in directoryList:
        if disk_fall(directory):
            break
        readerFile = open(directory + "/" + dataFile, "r")
        reader = csv.reader(readerFile, delimiter = " ")
        for row in reader:
            if row[0].startswith("#time:"):
                time = float(row[0].split(":")[1])
            if row[0].startswith("#gravityAngle"):
                tilt.append(float(row[0].split(":")[1]))
            if row[0] == "0": #0 is the id of the upper disk
                theta = float(row[7]) #this is the angle at this time
                theta0 = float(row[28])
        vel.append((theta - theta0)/time)
        readerFile.close()
    vel = np.array(vel)
    tilt = np.array(tilt)
    vel = vel*(60/(2*np.pi)) #change units to rpm
    tilt = tilt/np.pi #change units to rho/pi
    if ax is not None:
        ax.plot(tilt, vel, label=ac, marker = "o")
        return
    plt.figure(figsize=(8,6))
    plt.title("Rotational velocity vs tilt for acceleration="+ac)
    plt.xlabel("Tilt (rho/pi)")
    plt.ylabel("Rotational velocity(rpm)")
    plt.plot(tilt, vel, marker = "o")
    if show:
        plt.show()
    else:
        plt.savefig("plot_vel_tilt_for_ac_"+ac+".png", format="png")
    plt.close()

def plot_links_ac_for_tilt(tilt, ax=None, show=False):
    """Plots the fraction of time a link is open a function of the
    acceleration for a fixed tilt given by <tilt>.  Works from within
    a directory corresponding to a single sweep.  Plots on the given
    ax, if ax=None it generates a standalone plot.

    """
    ac = []
    open1 = []
    open2 = []
    bothOpen = []
    directoryList = [d for d in os.listdir(".") if "tilt."+tilt in d]
    directoryList.sort()
    for directory in directoryList:
        if disk_fall(directory):
            break
        dataFile = directory + "/3disk_avglinkstat.out"
        readerFile = open(dataFile, "r")
        reader = csv.reader(readerFile, delimiter = " ")
        for row in reader:
            op_op = float(row[0])
            op_sl = float(row[1])
            op_cl = float(row[2])
            sl_op = float(row[3])
            sl_sl = float(row[4])
            sl_cl = float(row[5])
            cl_op = float(row[6])
            cl_sl = float(row[7])
            cl_cl = float(row[8])
        readerFile.close()
        bothOpen.append(op_op)
        open1.append(op_sl + op_cl)
        open2.append(sl_op + cl_op)
        dataFile = directory + "/ending_phase_space.out"
        readerFile = open(dataFile, "r")
        reader = csv.reader(readerFile, delimiter = " ")
        for row in reader:
            if row[0].startswith("#dimen"):
                ac.append(float(row[0].split(":")[1]))
        readerFile.close()
    if ax is not None:
        ax.plot(ac, open1, label = "Link 1 opened, 2 closed.", marker = "o")
        ax.plot(ac, open2, label = "Link 2 opened, 2 closed.", marker = "o")
        ax.plot(ac, bothOpen, label = "Both links opened.", marker = "o")
        return
    plt.figure(figsize=(8,6))
    plt.title("Probability of links states vs acceleration for tilt="+tilt)
    plt.xlabel("Dimensionless acceleration")
    plt.ylabel("State probability")
    plt.plot(ac, open1, label = "Link 1 opened, 2 closed.", marker = "o")
    plt.plot(ac, open2, label = "Link 2 opened, 1 closed.", marker = "o")
    plt.plot(ac, bothOpen, label = "Both links opened.", marker = "o")
    if show:
        plt.show()
    else:
        plt.savefig("plot_links_ac_for_tilt_"+tilt+".png", format="png")
    plt.close()

def plot_links_tilt_for_ac(ac, ax=None, show=False):
    """Plots the fraction of time a link is open as a function of the tilt
    for a fixed acceleration given by <ac>. Works from within a
    directory corresponding to a single sweep. Plots on the given ax,
    if ax=None it generates a standalone plot.

    """
    tilt = []
    open1 = []
    open2 = []
    bothOpen = []
    directoryList = [d for d in os.listdir(".") if "ac."+ac in d]
    directoryList.sort()
    for directory in directoryList:
        if disk_fall(directory):
            break
        dataFile = directory + "/3disk_avglinkstat.out"
        readerFile = open(dataFile, "r")
        reader = csv.reader(readerFile, delimiter = " ")
        for row in reader:
            op_op = float(row[0])
            op_sl = float(row[1])
            op_cl = float(row[2])
            sl_op = float(row[3])
            sl_sl = float(row[4])
            sl_cl = float(row[5])
            cl_op = float(row[6])
            cl_sl = float(row[7])
            cl_cl = float(row[8])
        readerFile.close()
        bothOpen.append(op_op)
        open1.append(op_sl + op_cl)
        open2.append(sl_op + cl_op)
        dataFile = directory + "/ending_phase_space.out"
        readerFile = open(dataFile, "r")
        reader = csv.reader(readerFile, delimiter = " ")
        for row in reader:
            if row[0].startswith("#gravityAngle"):
                tilt.append(float(row[0].split(":")[1]))
        readerFile.close()
    tilt = np.array(tilt)
    tilt = tilt/np.pi #change units (rho/pi)
    if ax is not None:
        ax.plot(tilt, open1, label = "Link 1 opened, 2 closed.", marker = "o")
        ax.plot(tilt, open2, label = "Link 2 opened, 2 closed.", marker = "o")
        ax.plot(tilt, bothOpen, label = "Both links opened.", marker = "o")
        return
    plt.figure(figsize=(8,6))
    plt.title("Probability of links states vs tilt for ac="+ac)
    plt.xlabel("Tilt (rho/pi)")
    plt.ylabel("State probability")
    plt.plot(tilt, open1, label = "Link 1 opened, 2 closed.", marker = "o")
    plt.plot(tilt, open2, label = "Link 2 opened, 1 closed.", marker = "o")
    plt.plot(tilt, bothOpen, label = "Both links opened.", marker = "o")
    if show:
        plt.show()
    else:
        plt.savefig("plot_links_tilt_for_ac_"+ac+".png", format="png")
    plt.close()

def plot_compare_y(show=False, iTime=-float("inf"), fTime=float("inf")):
    fig = plt.figure(figsize=(24,16))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    plot_coord_time("y", iTime=iTime, fTime=fTime, diskId = "0", ax=ax1)
    plot_coord_time("y", iTime=iTime, fTime=fTime, diskId = "1", ax=ax3)
    plot_coord_time("x", iTime=iTime, fTime=fTime, diskId = "0", ax=ax4)
    plot_coord_time("w", iTime=iTime, fTime=fTime, diskId = "0", ax=ax2)
    plt.tight_layout()
    if show:
        plt.show()
    else:
        plt.savefig("compare.png", format="png")

def plot_supperpose(show=False, iTime=-float("inf"), fTime=float("inf"),
                    coord1="y", coord2="y", diskId1="0", diskId2="0"):
    fig = plt.figure(figsize=(12,8))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    plot_coord_time(coord1, diskId=diskId1, iTime=iTime, fTime=fTime, ax=ax1,
                    legend=False, color="b")
    plot_coord_time(coord2, diskId=diskId2, iTime=iTime, fTime=fTime, ax=ax2,
                    legend=False, color="r")
    if show:
        plt.show()
    else:
        plt.savefig("compare_" + coord1 + "_" + diskId1 + "_"
                    + coord2 + "_" + diskId2 + ".png", format="png")

def plot_compare_supperpose(show=False,
                            iTime=-float("inf"), fTime=float("inf")):
    fig = plt.figure(figsize=(24,16))
    ax1 = fig.add_subplot(221)
    ax1b = ax1.twinx()
    ax2 = fig.add_subplot(222)
    ax2b = ax2.twinx()
    ax3 = fig.add_subplot(223)
    ax3b = ax3.twinx()
    ax4 = fig.add_subplot(224)
    plot_coord_time(coord="y", diskId="1", iTime=iTime, fTime=fTime, ax=ax1,
                    legend=False, color="b")
    plot_coord_time(coord="y", diskId="0", iTime=iTime, fTime=fTime, ax=ax1b,
                    legend=False, color="r")
    plot_coord_time(coord="y", diskId="1", iTime=iTime, fTime=fTime, ax=ax2,
                    legend=False, color="b")
    plot_coord_time(coord="x", diskId="0", iTime=iTime, fTime=fTime, ax=ax2b,
                    legend=False, color="r")
    plot_coord_time(coord="y", diskId="1", iTime=iTime, fTime=fTime, ax=ax3,
                    legend=False, color="b")
    plot_coord_time(coord="w", diskId="0", iTime=iTime, fTime=fTime, ax=ax3b,
                    legend=False, color="r")
    plot_center_path(iTime=9.9, fTime=9.9125, ax=ax4)
    plt.tight_layout()
    if show:
        plt.show()
    else:
        plt.savefig("compare_sup.png", format="png")

def link_stats_coord(coord, directory=".", show=False, iTime=-float("inf"),
                     fTime = float("inf")):
    dataFile = directory + "/3disk_linkstat.out"
    time = []
    linkStatus1 = []
    linkStatus2 = []
    readerFile = open(dataFile, "r")
    reader = csv.reader(readerFile, delimiter = " ")
    for row in reader:
        currTime = float(row[0])
        if iTime <= currTime <= fTime:
            time.append(currTime)
            linkStatus1.append(int(row[1]) + int(row[2]))
            linkStatus2.append(int(row[7]) + int(row[8]))
    readerFile.close()
    time = np.array(time)
    linkStatus1 = np.array(linkStatus1)
    linkStatus2 = np.array(linkStatus2)
    colors1 = linkStatus1/2.0
    colors2 = linkStatus2/2.0
    dx = abs(time[1]-time[0])
    gs = gridspec.GridSpec(10,1)
    gs.update(hspace=0.05)
    ax1 = plt.subplot(gs[:-2,:])
    ax2 = plt.subplot(gs[8,:])
    ax3 = plt.subplot(gs[9,:])
    plot_coord_time(coord, iTime=iTime, fTime=fTime, ax=ax1)
    for t, color in zip(time, colors1):
        ax2.add_artist(Rectangle(xy=(t, 0), color=plt.cm.binary(color),
                                width=dx, height=1))
    for t, color in zip(time, colors2):
        ax3.add_artist(Rectangle(xy=(t, 0), color=plt.cm.binary(color),
                                width=dx, height=1))

    ax1.xaxis.set_major_locator(plt.NullLocator())
    ax1.set_xlim([iTime, fTime])
    ax2.set_ylabel("1")
    ax2.set_xlim([iTime, fTime])
    ax2.set_ylim([0,1])
    ax2.xaxis.set_major_locator(plt.NullLocator())
    ax2.yaxis.set_major_locator(plt.NullLocator())
    ax3.set_ylabel("2")
    ax3.set_xlim([iTime, fTime])
    ax3.set_ylim([0,1])
    ax3.yaxis.set_major_locator(plt.NullLocator())
    ax3.set_xlabel("time(s)")
    if show:
        plt.show()