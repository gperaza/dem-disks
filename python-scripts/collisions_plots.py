import numpy as np
import matplotlib.pyplot as plt
import csv

dataFileN = "collisions.out"
dataFile = open(dataFileN, "r")
reader = csv.reader(dataFile, delimiter=" ")
# count total collisions
for row in reader:
    if row[0].startswith("#EndCollision"):
        nCollisions = int(row[0].split(":")[1])
print "Total number of collisions:", nCollisions
dataFile.seek(0)  # reset input file

# For all the collisions.
# Plot scatter plots for the pre and post collision variables.
for row in reader:
    if row[0].startswith("#EndCollision"):
        tc = float(row[1])
        rad0 = float(row[8])
        mEff = float(row[7])  # since the other disk is fixed
        iMoment = float(row[6])
        alpha = mEff/3.0  # for collisions between disks
        deltaV0N = float(row[2])
        deltaV0T = float(row[3])
        deltaVW = float(row[5])*rad0/2.0
        print deltaV0T, deltaVW
# get N collisions at random, use a defined seed.
# Plot time evolution of N collisions.
