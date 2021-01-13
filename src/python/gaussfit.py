#! /usr/bin/python

from pylab import *
from scipy import optimize
from copy import copy

def generate_hist():
	vals = []
	for i in range(100000):
		vals.append(normal())
	h = hist(vals, bins=100)
	clf()
	return h[0:2]

def gaussvalue(params, x):
	center, sigma, A = params
	return A * e**(-((x-center)**2) / (2 * sigma**2)) / (sigma * sqrt( 2 * pi))

def trygauss(params, args):
	bins, points = args
	gausspoints = gaussvalue(params, bins)
	return points - gausspoints
	

def fit_hist(bins, points):
	print bins
	print points
	plot(bins, points)
	show()

def get_fitparams(histogram,peakloc,peaksz,sigguess):    
 #  points, bins = generate_hist()
 #  bins = array(bins[:-1]) + 0.5 * (bins[1]-bins[0])
 #  points = array(points)

 #  center = 0.1
 #  sigma = 0.8
 #  A = 5900
 #  params = [center, sigma, A]
    
   bins = histogram[1][:] + 0.5 * (histogram[1][1]-histogram[1][0])
   points = histogram[0][:]
   params = [peakloc, sigguess, peaksz]
   results = optimize.leastsq(trygauss, copy(params), args = [bins, points])
   bestfitparams = results[0]
   return bestfitparams

#gausspoints = gaussvalue(bestfitparams, bins)
#plot(points, label="points")
#plot(gausspoints, label="gausspoints")
#legend()
#show()




