#!/usr/bin/python

"""
Joe's mpi4py file example
"""

from mpi4py import MPI
import quiet
import sys
world=MPI.COMM_WORLD

rank=world.rank
size=world.size

runlist=quiet.L2Runlist("/data4/quiet/runlist_l2.txt")
segments=runlist.segments().values()
filenames=[segment.filename("/data4/quiet/level2") for segment in segments]
my_filenames=filenames[rank::size]   #This is a really handy construction that splits a list into parts for each process

for filename in my_filenames:
    print "Processor %d handling file %s" % (rank,filename)
    #Do the work here
    pass

sys.exit(0)
