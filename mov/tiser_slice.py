#!/usr/bin/python3

from numpy import *
import sys
import os.path


inputfile = sys.argv[1]
outputdir = sys.argv[2]

filename = os.path.basename(inputfile)
dotpos = filename.rindex('.')
basename = filename[:dotpos]
extension = filename[dotpos+1:]

sliceSize = float(sys.argv[3])

sliceNumber = 0
outputFileName = "%s/%s.%09d.%s"%(outputdir,basename,sliceNumber,extension)
outputFile = open(outputFileName,'w')

f = open(inputfile, 'r')
while 1:
    lines = f.readlines(100000)
    if not lines:
        break
    for line in lines:
        a = line.split()
        t = float(a[0])
        currentTime = sliceNumber*sliceSize
        if t >= currentTime:
            outputFile.close()
            outputFileName = "%s/%s.%09d.%s"%(outputdir,basename,sliceNumber,extension)
            print(outputFileName)
            outputFile = open(outputFileName,'w')
            sliceNumber += 1
        outputFile.write(line)

f.close()
