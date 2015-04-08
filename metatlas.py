from matplotlib import pyplot as plt
import requests, json
import numpy as np
import os
# import metatlas
from scipy.optimize import leastsq
from math import exp
import getpass
from elements import ELEMENTS
import re

def authenticateUser(client,username):
    password = getpass.getpass()
    authURL = 'https://metatlas.nersc.gov/client/login/'
    # Retrieve the CSRF token first
    client.get(authURL)  # sets cookie
    csrftoken = client.cookies['csrftoken']
    login_data = dict(username=username, password=password, csrfmiddlewaretoken=csrftoken, next='/')
    r = client.post(authURL, data=login_data, headers=dict(Referer=authURL))
    return client
# def authenticateUser(userFile):
#     authURL = 'https://metatlas.nersc.gov/client/login/'
#     file = open(userFile, 'r')
#     userID = file.readline()[:-1]
#     userPassword = file.readline()[:-1]
#     file.close()

#     client = requests.Session()
#     # Retrieve the CSRF token first
#     client.get(authURL)  # sets cookie
#     csrftoken = client.cookies['csrftoken']
#     login_data = dict(username=userID, password=userPassword, csrfmiddlewaretoken=csrftoken, next='/')
#     r = client.post(authURL, data=login_data, headers=dict(Referer=authURL))
#     return client

def export_peakData_to_spreadsheet(filename,export_fileIds,fileInfo,data,dictData):
    import csv
    export_filenames = []
    for i,myFile in enumerate(export_fileIds):
        for j,fid in enumerate(fileInfo['fid']):
            if fid == myFile:
                export_filenames.append(fileInfo['name'][j])
    fid = open(filename,'wb')
    fid.write('%s\t' % 'compound')
    for filename in export_filenames:
        fid.write('%s\t' % filename)
    fid.write('\n')
    for i,datum in enumerate(data):
        fid.write('%s\t' % dictData[u'compounds'][i]['name'])
        mz = float(dictData[u'compounds'][i][u'mz'])
        mzTol = float(dictData[u'compounds'][i][u'mz_threshold'])
        mzMin = mz - mz*mzTol/1.0e6
        mzMax = mz + mz*mzTol/1.0e6
        rtMin = float(dictData[u'compounds'][i][u'rt_min'])
        rtMax = float(dictData[u'compounds'][i][u'rt_max'])
        for j,myFile in enumerate(export_fileIds):
            if datum.size>0:
                idx = np.logical_and( datum[:,2]==myFile, datum[:,0]>=rtMin, datum[:,0]<=rtMax )
                if np.sum(idx)>0:
                    x1 = datum[:,0][idx]
                    y1 = datum[:,1][idx]
                    # y1 = y1 - np.min(y1)
                    myname = dictData[u'compounds'][i]['name']
                    if myname.startswith('ist'):
                        y1 = y1[:]
                    else:    
                        y1 = y1[:] / fileInfo['normalization_factor'][j]
                    fid.write('%5.2f\t' % np.sum(y1))
                else:
                    fid.write('%5.2f\t' % 0)
            else:
                fid.write('%5.2f\t' % 0)
        fid.write('\n')
    fid.close()

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



def getEICForCompounds_oneLighter(compound,myArray,files_I_want,rtTol,client,polarity):
    if isinstance(files_I_want,int):
        myList = str(files_I_want)
    else:
        myList = ','.join(map(str, files_I_want))
    mz = float(compound[u'mz'])
    mzTol = float(compound[u'mz_threshold'])
    mzMin = mz - mz*mzTol/1.0e6 - 1.003355
    mzMax = mz + mz*mzTol/1.0e6 - 1.003355
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
    if r.content:
        data = np.asarray(json.loads(r.content))
        return data
    else:
        return []
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

def createChromatogramPlots_dataOnly(data,compound,ax):
    ax.plot(data['xdata'],data['ydata']*data['iMax'],'k-',linewidth=2.0)
    ax.axvline(float(compound[u'rt_min']),linewidth=2, color='r',alpha=0.75) #original rtMin
    ax.axvline(float(compound[u'rt_max']),linewidth=2, color='r',alpha=0.75) #original rtMax
    ax.axvline(float(compound[u'rt_peak']),linewidth=2, color='g',alpha=0.75) #original rtPeak
    #     ax.axvline(x=compound[u'rt_peak'],linewidth=2, color='b') #original rtPeak
    # ax.axvline(x=fitResult[1],linewidth=2, color='r') #new rtPeak
    # ax.axvspan(fitResult[1]-fitResult[3]*2, fitResult[1]+fitResult[2]*2, facecolor='c', alpha=0.5) #new rtBounds
    ax.set_xlabel('Time (min)',weight='bold')
    ax.set_ylabel('Intensity (au)',weight='bold')
    ax.set_title(compound[u'name'])

def fitACompound(compound,data):
    rtPeak = float(compound[u'rt_peak'])
    rtMin = float(compound[u'rt_min'])
    rtMax = float(compound[u'rt_max'])
    init  = [1.0, rtPeak, 0.1,0.1]
    out   = leastsq( errfunc, init, args=(data['xdata'], data['ydata'], rtPeak, rtMin, rtMax))
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

def errfunc(p,x,y,rtPeak, rtMin, rtMax):
    if (abs(p[2]) > 1.5) or (abs(p[3]) > 1.5) or (abs(p[2]) < 0.001) or (abs(p[3]) < 0.001) or (p[1] > rtMax) or (p[1] < rtMin):
        return 1e100
    else:
        # return (y-fitfunc(p,x))**2
        # idx = x > rtMin and x < rtMax
        return np.multiply((y-fitfunc(p,x))**2,np.exp(-0.5*((x-rtPeak)/0.075)**2))

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

def chemformula_struct(formula):
    matEle = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    for idx, row in enumerate(matEle):
        id, num = row
        if num is '':
            matEle[idx] = (id, '1')
        else:
            matEle[idx] = (id, num)
        return matEle

def monoisotopicmass(formula):
    f = chemformula_struct(formula)
    m = np.zeros((len(f)))
    for i in range(len(f)):
        e = ELEMENTS[f[i][0]]
        maxIso = 0
        for iso in e.isotopes:
            if e.isotopes[iso].abundance > maxIso:
                maxIso = e.isotopes[iso].abundance
                if f[i][1] == '':
                    m[i] = 0
                else:
                    m[i] = float(e.isotopes[iso].mass) * float(f[i][1])
    return np.sum(m)

def chemformula_list(formula):
    f = chemformula_struct(formula)
    str = 'HCNOSPDX'
    fList = np.zeros((len(str)))
    for i in range(len(f)):
        if f[i][1] == '':
            fList[str.index(f[i][0])] = 0
        else:
            fList[str.index(f[i][0])] = f[i][1]
    return fList

def isotopic_pattern(fList,DAbund, num_D_sites):
    # fList is a list of formula vectors

    import numpy as np
    #loop through all the formulae to get the unique mz list
    mzvec = []
    for iii,fVec in enumerate(fList):
        x1,y1,m1 = isotope(fVec,DAbund)
        mzvec = np.union1d(mzvec,x1)
        labVec = fVec[:]
        labVec[6] = num_D_sites[iii]
        labVec[0] = fVec[0] - num_D_sites[iii]
        x2,y2,m2 = isotope(labVec,DAbund)
        mzvec = np.union1d(mzvec,x2)
    mzvec = np.round(mzvec*200)/200 #hack to round for unique
    mzvec = np.unique(mzvec)
    for i in range(len(mzvec)-1):
        if abs(mzvec[i]-mzvec[i+1])<0.1: #hack to group less than 0.1. should be parameter
            mzvec[i+1] = mzvec[i]
                
    mzvec = np.unique(mzvec)

    #go back through the formulae to put the intensity values onto the mz list
    isoY = np.zeros((len(fList)*2,len(mzvec)))
    for iii,fVec in enumerate(fList):
        x1,y1,m1 = isotope(fVec,DAbund)
        for i,x in enumerate(x1):
            idx = np.argmin(abs(x-mzvec))
            isoY[iii*2+0,idx] = y1[i]

        labVec = fVec[:]
        labVec[6] = num_D_sites[iii]
        labVec[0] = fVec[0] - num_D_sites[iii]
        x2,y2,m2 = isotope(labVec,DAbund)
        for i,x in enumerate(x2):
            idx = np.argmin(abs(x-mzvec))
            isoY[iii*2+1,idx] = y2[i]

    return mzvec, isoY

def isotope(fVec,DAbund):
    '''
    %
    % Calculates isotopic distributions including isotopic fine structure
    % of molecules using FFT and various scaling 'tricks'. Easily adopted
    % to molecules of any elemental composition (by altering MAX_ELEMENTS
    % and the nuclide matrix A). To simulate spectra, convolute with peak
    % shape using FFT.
    %
    % (C) 1999 by Magnus Palmblad, Division of Ion Physics, Uppsala Univ.
    % Acknowledgements:
    % Lars Larsson-Cohn, Dept. of Mathematical Statistics, Uppsala Univ.,
    % for help on theory of convolutions and FFT.
    % Jan Axelsson, Div. of Ion Physics, Uppsala Univ. for comments and ideas
    %
    % Contact Magnus Palmblad at magnus.palmblad@gmail.com if you should
    % have any questions or comments.
    %

    Converted to Python 1/10/08 by
    Brian H. Clowers bhclowers@gmail.com

    October 31, 2014
    Added Phosphorous and chemical formula parsing
    Added conditional specification of stable isotope composition
    Ben Bowen, ben.bowen@gmail.com

    fVec is a vector representing the chemical formula including deuterium
    # [H, C, N, O, S, P, D]
    DAbund is the amount of deuterium [0-1], 0.05 is typical
    '''
    import numpy as np
    import numpy.fft.fftpack as F
    # import time
    # import pylab as P


    def next2pow(x):
        return 2**int(np.ceil(np.log(float(x))/np.log(2.0)))

    scaleFactor = 100000
    MAX_ELEMENTS=7+1  # add 1 due to mass correction 'element'
    MAX_ISOTOPES=4    # maxiumum # of isotopes for one element
    CUTOFF=1e-4     # relative intensity cutoff for plotting

    WINDOW_SIZE = 500
    #WINDOW_SIZE=input('Window size (in Da) ---> ');

    #RESOLUTION=input('Resolution (in Da) ----> ');  % mass unit used in vectors
    RESOLUTION = 0.5
    if RESOLUTION < 0.00001:#  % minimal mass step allowed
      RESOLUTION = 0.00001
    elif RESOLUTION > 0.5:  # maximal mass step allowed
      RESOLUTION = 0.5

    R=0.00001/RESOLUTION#  % R is used to scale nuclide masses (see below)

    WINDOW_SIZE=WINDOW_SIZE/RESOLUTION;   # convert window size to new mass units
    WINDOW_SIZE=next2pow(WINDOW_SIZE);  # fast radix-2 fast-Fourier transform algorithm

    if WINDOW_SIZE < np.round(496708*R)+1:
      WINDOW_SIZE = nextpow2(np.round(496708*R)+1)  # just to make sure window is big enough

    # print 'Vector size: 1x%d'%WINDOW_SIZE

    #H378 C254 N65 O75 S6
    
    # M=np.array([378,254,65,75,6,0]) #% empiric formula, e.g. bovine insulin
    M=np.array(fVec) #% empiric formula, e.g. bovine insulin

    # isotopic abundances stored in matrix A (one row for each element)
    A=np.zeros((MAX_ELEMENTS,MAX_ISOTOPES,2));

    A[0][0,:] = [100783,0.9998443]#                 % 1H
    A[0][1,:] = [201410,0.0001557]#                 % 2H
    A[1][0,:] = [100000,0.98889]#                   % 12C
    A[1][1,:] = [200336,0.01111]#                   % 13C
    A[2][0,:] = [100307,0.99634]#                   % 14N
    A[2][1,:] = [200011,0.00366]#                   % 15N
    A[3][0,:] = [99492,0.997628]#                  % 16O
    A[3][1,:] = [199913,0.000372]#                  % 17O
    A[3][2,:] = [299916,0.002000]#                  % 18O
    A[4][0,:] = [97207,0.95018]#                   % 32S
    A[4][1,:] = [197146,0.00750]#                   % 33S
    A[4][2,:] = [296787,0.04215]#                   % 34S
    A[4][3,:] = [496708,0.00017]#                   % 36S
    A[5][0,:] = [97376,1.0]# Phosphorous
    A[6][0,:] = [100783,1.0-DAbund]#                 % 1H
    A[6][1,:] = [201410,DAbund]#                 % 2H
    A[7][0,:] = [100000,1.00000]#                   % for shifting mass so that Mmi is
    #                                             % near left limit of window
    mass_removed_vec = [0,11,13,15,31,30,0,-1]
    monoisotopic = 0.0
    for i,e in enumerate(fVec):
        monoisotopic = monoisotopic + ( (mass_removed_vec[i]*scaleFactor+A[i][0,0])*e / scaleFactor)

    Mmi=np.array([np.round(100783*R), np.round(100000*R),\
                 np.round(100307*R), np.round(99492*R), np.round(97207*R), np.round(97376*R), np.round(100783*R), 0])*M#  % (Virtual) monoisotopic mass in new units
    Mmi = Mmi.sum()
    #% mass shift so Mmi is in left limit of window:
    #print "Mmi",Mmi
    #print "Window", WINDOW_SIZE
    FOLDED=np.floor(Mmi/(WINDOW_SIZE-1))+1#  % folded FOLDED times (always one folding due to shift below)

    #% shift distribution to 1 Da from lower window limit:
    M[MAX_ELEMENTS-1]=np.ceil(((WINDOW_SIZE-1)-np.mod(Mmi,WINDOW_SIZE-1)+np.round(100000*R))*RESOLUTION)
    MASS_REMOVED=np.array(mass_removed_vec)*M#';  % correction for 'virtual' elements and mass shift
    MASS_REMOVED = MASS_REMOVED.sum()

    ptA=np.ones(WINDOW_SIZE);
    t_fft=0
    t_mult=0

    for i in xrange(MAX_ELEMENTS):
        tA=np.zeros(WINDOW_SIZE)
        for j in xrange(MAX_ISOTOPES):
            if A[i][j,0] != 0:
                tA[np.round(A[i][j,0]*R)]=A[i][j,1]#;  % put isotopic distribution in tA

        tA=F.fft(tA) # FFT along elements isotopic distribution  O(nlogn)
        tA=tA**M[i]#  % O(n)
        ptA = ptA*tA#  % O(n)#this is where it is messing UP

    ptA=F.ifft(ptA).real#;  % O(nlogn)

    start = (FOLDED*(WINDOW_SIZE-1)+1)*RESOLUTION+MASS_REMOVED,(FOLDED+1)*(WINDOW_SIZE-1)*RESOLUTION+MASS_REMOVED
    stop = WINDOW_SIZE - 1

    MA=np.linspace((FOLDED*(WINDOW_SIZE-1)+1)*RESOLUTION+MASS_REMOVED,(FOLDED+1)*(WINDOW_SIZE-1)*RESOLUTION+MASS_REMOVED, WINDOW_SIZE-1)

    ind=np.where(ptA>CUTOFF)[0]

    x = MA[ind]
    y = ptA[ind]

    for i,xi in enumerate(x):
        x[i] = monoisotopic + (i*1.003355)


    return x,y,monoisotopic