#!/usr/bin/python

from numpy import *
import getopt
import sys

# some default parameters
dt = 10e-3
windowsize = 1
popsize = 4096
ntriggers = 20
inputfile = ""
triggerfile = ""
maxtime  = 3600
skiptime = 0
normalize = False
corr_mat = False
show_graph = False
verbose = False
patfile = ""
thr=-1 # -1 for auto

def usage():
    '''Display command line help message'''
    print '''
This program makes spike triggered averages from ras files

-h  this screen
-f  input ras file
-t  trigger file
-s  time to skip in beginning (default 0s)
-m  maximum time (default 3600s)
-w  window size (default 1s)
-r  temporal resolution (default 0.01s)
-n  normalize
-p  population size
-o  save patfile
-q  threshold for patfile -1 is auto
-g  graphical out
-c  plot correlation matrix
'''
    sys.exit(0)

try:
    opts, args = getopt.getopt(sys.argv[1:], "hf:t:s:m:w:r:np:q:o:gc", ["help","file=","trigger=","skip=","max=","windowsize=","resolution=","normalize=","popsize=","outfile=","qthr=","graph","corr"])
except getopt.GetoptError, err:
    print 'invalid command line'
    usage()
    sys.exit(2)
for opt, arg in opts:
    if opt in ("-t", "--trigger"):
        triggerfile = arg
    if opt in ("-f", "--file"):
        inputfile = arg
    if opt in ("-s", "--skip"):
        skiptime = float(arg)
    if opt in ("-m", "--max"):
        maxtime = float(arg)
    if opt in ("-w", "--windowsize"):
        windowsize = float(arg)
    if opt in ("-r", "--resolution"):
        dt = float(arg)
    if opt in ("-p", "--popsize"):
        popsize = int(arg)
    if opt in ("-n", "--normalize"):
        normalize = True
    if opt in ("-c", "--corr"):
        corr_mat = True
    if opt in ("-o", "--outfile"):
        patfile = arg
    if opt in ("-q", "--qthr"):
        thr = float(arg)
        print "Setting thr",thr
    if opt in ("-g", "--graph"):
        show_graph = True
    elif opt in ("-h", "--help"):
        usage()

if not show_graph and patfile=="":
    show_graph = True
    from pylab import *

print "windowsize %.2fs"%windowsize

if patfile!="":
    fvector = open(patfile, 'w')

ras = open(inputfile, 'r')
trigger = open(triggerfile, 'r')
hist = zeros((ntriggers,int(windowsize/dt),popsize))
tcounts = zeros(ntriggers)


print "Scanning ras file ..."

for triggerline in trigger:
    tokens = triggerline.split()
    triggert = float(tokens[0])
    triggerno = int(tokens[1])
    if len(tokens)>2: # interpret last entry as duration
        duration = min(float(tokens[2]),windowsize)
    else:
        duration = windowsize
    t = 0 
    if triggerno < ntriggers and triggert > skiptime:
        if verbose:
            print triggert,triggerno
        if duration < windowsize: # skip
            continue
        tcounts[triggerno] += 1
        while t<duration:
            rasline = ras.readline()
            a,b = rasline.split()
            t = float(a) - triggert
            i = int(b)
            if t > 0 and t < duration and i < popsize:
                hist[triggerno,int(t/dt),i] += 1

    if triggert > maxtime:
        break

print "%i triggers scanned"%tcounts.sum()
for i in range(ntriggers):
    print " - %i counts in %i"%(tcounts[i],i)
print "plotting..."

ras.close()
trigger.close()

countvectors = list()

for a in range(ntriggers):
    if tcounts[a]:
        h = hist[a].T/tcounts[a]/windowsize
    if normalize:
        h /= h.max()
    if tcounts[a]>0:
        countvectors.append(h.sum(1))
        # h = h[h.sum(1)!=0]
        if show_graph:
            figure()
            imshow(h,interpolation='nearest',origin='lower',aspect='auto')
            colorbar()
        if patfile!="":
            sigma = std(array(countvectors))
            if thr > -1:
                t = thr
            else: 
                t = mean(array(countvectors))
            for i,e in enumerate(h.sum(1)):
                if e > t:
                    fvector.write("%i %f\n"%(i,e))
            fvector.write("\n\n")

if corr_mat:
    d = array(countvectors)
    c = zeros((len(countvectors),len(countvectors)))
    for i in range(len(countvectors)):
        c[i,i] = dot(d[i],d[i])
        for j in range(0,i):
            c[i,j] = dot(d[i],d[j])
            c[j,i] = c[i,j]

    # c = corrcoef(d)
    print c
    if show_graph:
        figure()
        imshow(c,interpolation='nearest',origin='upper')
        colorbar()


if show_graph:
    show()
