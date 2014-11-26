from matplotlib import pyplot as plt
import requests, json
import numpy as np
import os
import metatlas
from scipy.optimize import leastsq
from math import exp

def shareExperiments(allUsers,allPerms,client, myExperimentID):
    for aUser in allUsers:
        payload = {"user":aUser,"perms":allPerms}
        sendData=json.dumps(payload)
        # print sendData
        url = 'https://metatlas.nersc.gov/api/experiment/%s/share/' % myExperimentID
        # print url
        r = client.post(url, data=sendData)
        return r
        # print r.content
    
def listMyExperiments(client):
    url = 'https://metatlas.nersc.gov/api/experiment'
    r = client.get(url)
    experiments = json.loads(r.content)
    return experiments

def authenticateUser(userFile):
    authURL = 'https://metatlas.nersc.gov/client/login/'
    file = open(userFile, 'r')
    userID = file.readline()[:-1]
    userPassword = file.readline()[:-1]
    file.close()

    client = requests.Session()
    # Retrieve the CSRF token first
    client.get(authURL)  # sets cookie
    csrftoken = client.cookies['csrftoken']
    login_data = dict(username=userID, password=userPassword, csrfmiddlewaretoken=csrftoken, next='/')
    r = client.post(authURL, data=login_data, headers=dict(Referer=authURL))
    return client

def getEICForCompounds(compound,myArray,files_I_want,rtTol,client,polarity):
	if isinstance(files_I_want,int):
		myList = str(files_I_want)
	else:
		myList = ','.join(map(str, files_I_want))
	mz = float(compound[u'mz'])
	mzTol = float(compound[u'mz_threshold'])
	mzMin = mz - mz*mzTol/1.0e6
	mzMax = mz + mz*mzTol/1.0e6
	rtMin = float(compound[u'rt_min'])-rtTol
	rtMax = float(compound[u'rt_max'])+rtTol
	rtPeak = float(compound[u'rt_peak'])

	payload = {'L':1,'P':polarity,'arrayname':myArray,'fileidlist':myList,
	          'max_mz':mzMax,'min_mz':mzMin,
	          'min_rt':rtMin,'max_rt':rtMax,
	          'nsteps':20000,'queryType':'XICofFile_mf'}
	url = 'https://metatlas.nersc.gov/api/run'
	r = client.get(url,params=payload)
	data = np.asarray(json.loads(r.content))
	return data
	# for myFile in files_I_want:
	#     x1 = data[:,0][(data[:,2]==myFile)]
	#     y1 = data[:,1][(data[:,2]==myFile)]
	#     idx = np.argsort(x1)
	#     plt.plot(x1[idx],y1[idx])
	# plt.xlabel('Time (min)')
	# plt.ylabel('TIC Intensity (au)')

def getEICForCompound(compound,myArray,runId,rtTol,client,polarity):
    #polarity is 1 for pos and 0 for neg
	mz = float(compound[u'mz'])
	mzTol = float(compound[u'mz_threshold'])
	mzMin = mz - mz*mzTol/1.0e6
	mzMax = mz + mz*mzTol/1.0e6
	rtMin = float(compound[u'rt_min'])
	rtMax = float(compound[u'rt_max'])
	rtPeak = float(compound[u'rt_peak'])
	payload = {'L':1,'P':polarity,'arrayname':myArray,'fileid':runId,
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
	ax.axvline(float(compound[u'rt_min']),linewidth=2, color='k') #original rtMin
	ax.axvline(float(compound[u'rt_max']),linewidth=2, color='k') #original rtMax
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
	if (abs(p[2]) > 0.5) or (abs(p[3]) > 0.5) or (abs(p[2]) < 0.01) or (abs(p[3]) < 0.01):
		return 1e100
	else:
		# return (y-fitfunc(p,x))**2
		return np.multiply((y-fitfunc(p,x))**2,np.exp(-0.5*((x-rtPeak)/0.1)**2))

#!/usr/bin/env python

'''Implementation of Sliding Window Minimum Algorithm. This module contains
one function `sliding_window_minimum`.
See http://people.cs.uct.ac.za/~ksmith/articles/sliding_window_minimum.html
for a longer explanation of the algorithm.'''

from collections import deque

def sliding_window_minimum(k, li):
    '''
    A iterator which takes the size of the window, `k`, and an iterable,
    `li`. Then returns an iterator such that the ith element yielded is equal
    to min(list(li)[max(i - k + 1, 0):i+1]).
    Each yield takes amortized O(1) time, and overall the generator takes O(k)
    space.
    __author__ = "Keegan Carruthers-Smith"
	__email__ = "keegan.csmith@gmail.com"
	__license__ = "MIT"
    '''

    window = deque()
    for i, x in enumerate(li):
        while window and window[-1][0] >= x:
            window.pop()
        window.append((x, i))
        while window[0][1] <= i - k:
            window.popleft()
        yield window[0][0]



 
def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
	
    import sys
    from numpy import NaN, Inf, arange, isscalar, asarray, array
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True
 
    return array(maxtab), array(mintab)

def splitFileListAcrossTwoDimensions(myFiles):
    import re
    fileInfo = {}
    fileInfo['date']=[]
    fileInfo['polarity']=[]
    fileInfo['conc']=[] #this is the concentration of material
    fileInfo['temp']=[] #this is the temperature of the inubation
    fileInfo['group']=[] #this is a unique class of time by group

    for file in myFiles:
        oldName = file
        file = file[:-5]
        file = re.sub('blank\d','blank_0',file)
        file = re.sub('RT','23',file) #substitute the room temperature with 23 for consistency with naming
        fileInfo['polarity'].append(file[-3:])
        file = file[:-4]
        fileInfo['date'].append(file[:6])
        file = file[7:]
        temperature = re.findall(('_[0-9]+$'),file)[0]
        fileInfo['temp'].append(temperature.replace('_',''))
        file = re.sub(temperature+'$','',file)
        file = file[:-2] # strip off the _replicate
        file = file.replace('_','.')
        file = file.replace('bla','blank')
        fileInfo['conc'].append(file)
        fileInfo['group'].append(fileInfo['conc'][-1]+'@'+fileInfo['temp'][-1])
    return fileInfo

def groupFilesAcrossTwoDimensions(fileInfo):
# This block is an example to teach people how to use regular expressions to put files into groups according to their filename
# by naming your files in a consistent manner, it can make analysis of N-way comparisons much quicker
	import re
	fileInfo['date']=[]
	fileInfo['polarity']=[]
	fileInfo['conc']=[] #this is the concentration of material
	fileInfo['temp']=[] #this is the temperature of the inubation
	fileInfo['group']=[] #this is a unique class of time by group

	for file in fileInfo['name']:
	    oldName = file
	    file = file[:-5]
	    file = re.sub('blank\d','blank_0',file)
	    file = re.sub('RT','23',file) #substitute the room temperature with 23 for consistency with naming
	    fileInfo['polarity'].append(file[-3:])
	    file = file[:-4]
	    fileInfo['date'].append(file[:6])
	    file = file[7:]
	    temperature = re.findall(('_[0-9]+$'),file)[0]
	    fileInfo['temp'].append(temperature.replace('_',''))
	    file = re.sub(temperature+'$','',file)
	    file = file[:-2] # strip off the _replicate
	    file = file.replace('_','.')
	    file = file.replace('bla','blank')
	    fileInfo['conc'].append(file)
	    fileInfo['group'].append(fileInfo['conc'][-1]+'@'+fileInfo['temp'][-1])

	fileInfo['pos_groups']={}
	myGroups = np.unique(fileInfo['group'])
	# # print fileInfo['group']
	for group in myGroups:
	    indices = [i for i, elem in enumerate(fileInfo['group']) if group == elem]
	    fileInfo['pos_groups'][group] = []
	    for i in indices:
	        if fileInfo['polarity'][i] == 'pos':
	            fileInfo['pos_groups'][group].append(fileInfo['fid'][i])
	    # print fileInfo['pos_groups'][group]
	# print fileInfo['pos_groups'].keys()

	fileInfo['neg_groups']={}
	myGroups = np.unique(fileInfo['group'])
	# # print fileInfo['group']
	for group in myGroups:
	    indices = [i for i, elem in enumerate(fileInfo['group']) if group == elem]
	    fileInfo['neg_groups'][group] = []
	    for i in indices:
	        if fileInfo['polarity'][i] == 'neg':
	            fileInfo['neg_groups'][group].append(fileInfo['fid'][i])
	    # print group
	    # print fileInfo['neg_groups'][group]
	# print fileInfo['neg_groups'].keys()
	return fileInfo
