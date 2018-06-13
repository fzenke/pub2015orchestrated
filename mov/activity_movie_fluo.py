#!/usr/bin/python

from numpy import *
import sys
import os.path


N = 4096
aside = int(sqrt(N))

inputfile = sys.argv[1]
outputdir = sys.argv[2]
patternfile = "../src/data/rf1.pat"

filename = os.path.basename(inputfile)
dotpos = filename.rindex('.')
basename = filename[:dotpos]
extension = filename[dotpos+1:]

sliceSize = float(sys.argv[3])

interval_begin = 50
interval_end   = 360


dt = 1e-3
tau = 1.0
ctime = 0.
frameRate = 25
frameSize = 1./frameRate

decayMultiplier = exp(-dt/tau)


# read patterns
patterns = list()
pat = list()
f = open(patternfile, 'r')
while 1:
    lines = f.readlines(10000)
    if not lines:
        break
    for line in lines:
        a = line.split()
        if len(a):
            if a[0][0] == "#": # skip comments
                continue
            i = int(a[0])
            g = 1.0
            if len(a)>1:
                g = float(a[1])
            pat.append(i)
        elif len(pat): # new pattern if current not empty (ignore extra blank lines)
            patterns.append(pat)
            pat = list()

f.close()

densepat = zeros((len(patterns),aside*aside))
for i,pat in enumerate(patterns):
    for u in pat:
        densepat[i,u] = 1

activity = zeros(N,dtype=float)

sliceNumber = 0
subFrame = 0 

def write_file(a,fname):
    o = open(fname,'w')
    for i in range(aside):
        for j in range(aside):
            o.write("%d %d %f\n"%(i,j,a[i+j*aside]))

    for pat in patterns:
        o.write("\n\n")
        for i in pat:
            o.write("%d %d %f\n"%(i%aside,i/aside,a[i]))


    o.close()

f = open(inputfile, 'r')
ftime = 0
ctime = 0
frame = 0
while 1:
    lines = f.readlines(100000)
    if not lines:
        break
    for line in lines:
        a,b = line.split()
        t = float(a)
        i = int(b)
        if t<(interval_begin-10*tau):
            ftime = t
            frame = ftime*frameRate
            continue
        if t>interval_end+10*tau:
            break

        while (ctime < t):
            activity *= decayMultiplier
            ctime += dt

            if ftime < ctime and ctime>=interval_begin and ftime<=interval_end+dt:
                sliceNumber = frame/frameRate/sliceSize
                subFrame = frame%(sliceSize*frameRate)

                outputFileName = "%s/%s.%09d.%03d.%s"%(outputdir,basename,sliceNumber,subFrame,"act")
                print(outputFileName)
                write_file(activity,outputFileName)

                frame += 1
                ftime = frame*frameSize
        activity[i] += 1


