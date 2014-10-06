from matplotlib import pyplot as plt
import requests, json
import numpy as np
import os
import metatlas
from scipy.optimize import leastsq
from math import exp


def getEICForCompound(compound,myArray,runId,rtTol,client):
	mz = float(compound[u'mz'])
	mzTol = float(compound[u'mz_threshold'])
	mzMin = mz - mz*mzTol/1.0e6
	mzMax = mz + mz*mzTol/1.0e6
	rtMin = float(compound[u'rt_min'])
	rtMax = float(compound[u'rt_max'])
	rtPeak = float(compound[u'rt_peak'])
	payload = {'L':1,'P':1,'arrayname':myArray,'fileid':runId,
	'max_mz':mzMax,'min_mz':mzMin,
	'nsteps':10000,'queryType':'XICofFile'}
	url = 'https://metatlas.nersc.gov/api/run'
	r = client.get(url,params=payload)
	data = np.asarray(json.loads(r.content))
	xdata    = data[abs(data[:,0]-rtPeak)<rtTol,0]
	ydata    = data[abs(data[:,0]-rtPeak)<rtTol,1]
	peakArea = data[(data[:,0]>rtMin) & (data[:,0]<rtMax),1]
	if len(peakArea)>0:
		peakArea = sum(peakArea)
	else:
		peakArea = 0
	if len(xdata)>0:
		iMin = min(ydata)
		ydata = ydata - iMin
		iMax = max(ydata) + iMin
		ydata = ydata / iMax
	else:
		iMax = 0
	return {'eic':data,'xdata':xdata,'ydata':ydata,'name':compound[u'name'],'iMax':iMax,'peakArea':peakArea}

def createChromatogramPlots(data,compound,fitResult,ax):
	ax.plot(data['xdata'],data['ydata']*data['iMax'],'k-',data['xdata'], fitfunc(fitResult, data['xdata'])*data['iMax'],'r-',linewidth=2.0)
	ax.axvline(compound[u'rt_min'],linewidth=2, color='k') #original rtMin
	ax.axvline(compound[u'rt_max'],linewidth=2, color='k') #original rtMax
	#     ax.axvline(x=compound[u'rt_peak'],linewidth=2, color='b') #original rtPeak
	ax.axvline(x=fitResult[1],linewidth=2, color='r') #new rtPeak
	ax.axvspan(fitResult[1]-fitResult[3]*2, fitResult[1]+fitResult[2]*2, facecolor='c', alpha=0.5) #new rtBounds
	ax.set_xlabel('Time (min)')
	ax.set_ylabel('Intensity (au)')
	ax.set_title(compound[u'name'])

def fitACompound(compound,data):
	rtPeak = float(compound[u'rt_peak'])
	init  = [1.0, rtPeak, 0.1,0.1]
	out   = leastsq( errfunc, init, args=(data['xdata'], data['ydata'], rtPeak))
	fitResult = out[0]
	fitResult[2] = abs(fitResult[2])
	fitResult[3] = abs(fitResult[3])
	return fitResult


def fitfunc(p,x):
	a = np.zeros(x.shape)
	idx1 = x>=p[1]
	if len(idx1)>0:
		a[idx1]=p[0]*np.exp(-0.5*((x[idx1]-p[1])/p[2])**2)
	idx2 = x<p[1]
	if len(idx2)>0:
		a[idx2]=p[0]*np.exp(-0.5*((x[idx2]-p[1])/p[3])**2)
	return a

def errfunc(p,x,y,rtPeak):
	if (abs(p[2]) > 0.3) or (abs(p[3]) > 0.3) or (abs(p[2]) < 0.01) or (abs(p[3]) < 0.01):
		return 1e100
	else:
		return np.multiply((y-fitfunc(p,x)),np.exp(-0.5*((x-rtPeak)/0.1)**2))