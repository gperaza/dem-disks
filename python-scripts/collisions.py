# A script to identify and extract data for collisions.

import numpy as np
import csv

# We need to identify collisions first. The following code uses the
# link's status to search for collision and generate a list with
# collision times.
linksFileName = "3disk_linkstat.out"
linksFile = open(linksFileName, "r")
reader = csv.reader(linksFile, delimiter=" ")
linksFileNameW = "3disk_linkstat_col.out"
linksFileW = open(linksFileNameW, "w")
writer = csv.writer(linksFileW, delimiter=" ")
touching1 = 0
touching2 = 0
colliding = 0
collisionList = []
writeList = []
fn = []
ft = []
ftt = []
fnt = []
for row in reader:
    preTouching1 = touching1
    preTouching2 = touching2
    touching1 = int(row[1])
    touching2 = int(row[7])
    if not colliding:  # Search for a collision.
        if (touching1 + touching2) == 1 and (preTouching1 + preTouching2) == 0:
            if touching1:
                colliding = 1
            else:
                colliding = 2
            sTime = float(row[0])
            writeList.append(["#StartCollision"])
            writeList.append(row)
            fnt.append(float(row[3]))
            ftt.append(float(row[4]))
    else:  # Does collision continues?
        if touching1 == preTouching1 and touching2 == preTouching2:
            # Collision continues, do nothing.
            writeList.append(row)
            fnt.append(float(row[3]))
            ftt.append(float(row[4]))
        elif not touching1 and not touching2:
            # If disk is flying, collision ended successfully.
            eTime = float(row[0])
            collisionList.append((sTime, eTime, str(colliding)))
            colliding = False
            writeList.append(row)
            writeList.append(["#EndCollision"])
            fnt.append(float(row[3]))
            ftt.append(float(row[4]))
            for row2 in writeList:
                writer.writerow(row2)
            writeList = []
            for f in fnt:
                fn.append(f)
            for f in ftt:
                ft.append(f)
            ftt = []
            fnt = []
        else:
            # All other possibilities invalidate the collision.
            colliding = False
            writeList = []
            ftt = []
            fnt = []
linksFile.close()
linksFileW.close()

# We have the list with collision times, use this to extract relevant
# information from the phase space file.
phaseFile = open("phase_space.out", "r")
reader = csv.reader(phaseFile, delimiter=" ")
collisionFile = open("collisions.out", "w")
writer = csv.writer(collisionFile, delimiter=" ")
for row in reader:
    if row[0].startswith("#time:"):
        iTime = float(row[0].split(":")[1])
        break
while collisionList[0][0] < iTime:
    collisionList.pop(0)

phaseFile.seek(0)
sTime = collisionList[0][0]
fTime = collisionList[0][1]
collidingWith = collisionList[0][2]
collisionNo = 0
dt = 1e-8
starting = False
ending = False
for row in reader:
    # Check if collision is starting or ending.
    if row[0].startswith("#time:"):
        currTime = float(row[0].split(":")[1])
        if abs(currTime - sTime) < dt:
            colliding = True
            starting = True
            writer.writerow(["#StartCollision:"+str(collisionNo)])
            sTime = currTime
        elif abs(currTime - fTime) < dt:
            ending = True
            colliding = False
            eTime = currTime

    # Get the data.
    if row[0] == "0" and (colliding or ending):
        # get data
        rad0 = float(row[2])
        mass0 = float(row[3])
        imoment0 = float(row[4])
        x0 = float(row[5])
        y0 = float(row[6])
        w0 = float(row[7])
        vx0 = float(row[8])
        vy0 = float(row[9])
        vw0 = float(row[10])
    if row[0] == collidingWith and (colliding or ending):
        # get data
        rad1 = float(row[2])
        mass1 = float(row[3])
        imoment1 = float(row[4])
        x1 = float(row[5])
        y1 = float(row[6])
        w1 = float(row[7])
        vx1 = float(row[8])
        vy1 = float(row[9])
        vw1 = float(row[10])
        # calculate data
        ex = np.array([1, 0, 0])
        ey = np.array([0, 1, 0])
        ez = np.array([0, 0, 1])
        r0 = np.array([x0, y0, 0])
        r1 = np.array([x1, y1, 0])
        r12 = r0 - r1
        v0 = np.array([vx0, vy0, 0])
        v1 = np.array([vx1, vy1, 0])
        v12 = v0 - v1
        eNorm = r12/np.sqrt(np.dot(r12, r12))
        eTan = np.cross(eNorm, ez)
        v0N = np.dot(v0, eNorm)
        v0T = np.dot(v0, eTan)
        v0TT = v0T + vw0*rad0
        v12N = np.dot(v12, eNorm)
        v12T = np.dot(v12, eTan)
        v12TT = v12T + vw0*rad0
        # write data
        writer.writerow([v12N, v12T, v12TT, fn.pop(0), ft.pop(0)])
        if starting:
            v0NS = v0N
            v0TTS = v0TT
            v0TS = v0T
            vw0S = vw0
            starting = False
        # we need to get normal and tangential velocities, normal and
        # tangential forces, normal and tangential impulses, Cartesian
        # velocities, Cartesian impulses.
        if ending:
            deltaV0N = v0N - v0NS
            deltaV0T = v0T - v0TS
            deltaV0TT = v0TT - v0TTS
            deltaVW = vw0 - vw0S
            writer.writerow(["#EndCollision:"+str(collisionNo),
                             eTime - sTime, deltaV0N, deltaV0T,
                             deltaV0TT, deltaVW, imoment0, mass0, rad0])
            collisionNo += 1
            ending = False
            if collisionNo < len(collisionList):
                sTime = collisionList[collisionNo][0]
                fTime = collisionList[collisionNo][1]

phaseFile.close()
collisionFile.close()
assert(collisionNo == len(collisionList))
assert(len(fn) == 0)

collisionsFile = open("collisions.out", "r")
reader = csv.reader(collisionsFile, delimiter=" ")
for row in reader:
    pass
