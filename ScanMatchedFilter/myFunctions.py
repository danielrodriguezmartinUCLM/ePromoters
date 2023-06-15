# Fecha 01/07/2022
import math
import os
import sys
import subprocess
from subprocess import Popen, PIPE, STDOUT
import metaseq
from collections import OrderedDict
import numpy as np
import scipy
from scipy import interpolate
import pylab as P
import pybedtools as pbt
from scipy.optimize import curve_fit
#-------------------------------------------------------------------------
def deleteFile(fileName):
    try:
        os.remove(fileName)
        print ">>>> File : "+fileName + " erased correctly\n"
    except OSError as e:
        print "[ERROR] deleting file : "+fileName+"\n"
#-------------------------------------------------------------------------
def writeStep(step):
        fConfig=open('Config.txt','w')
        fConfig.seek(0)
        fConfig.write(str(step))
        fConfig.close()
        #print "**** STEP "+str(step)+ " done correctly\n"
#-------------------------------------------------------------------------
def readStep():
    try:
        fConfig=open('Config.txt','r')
        paso=fConfig.read()
        fConfig.close()
        
    except IOError:
        writeStep(0)
        paso=0
    return paso
#-------------------------------------------------------------------------
def askStep(stepToExecute): # It asks if we want to carry out the Step, returns YES or NO
    savedStep=readStep()
    
    answer="YES"
    
    if int(stepToExecute)>int(savedStep):
        answer="YES"
    
    if int(stepToExecute)==int(savedStep):
        answer=raw_input (">>>> This step has already been done. Do you want to run this step again? (Yes-No) [No] ")
        answer=answer.upper()
        if answer != "YES":
            answer="NO"
    
    if int(stepToExecute)<int(savedStep):
        answer="NO"
        
    if answer=="YES":
        print ">>>> The step "+ str(stepToExecute) + " will take place.\n"
    else:
        print ">>>> The step "+ str(stepToExecute) + " has not been performed.\n"
    return answer
#-------------------------------------------------------------------------
def ExecuteStep(stepToExecute): # Perform the step if you haven't already
    stepSaved=readStep()
    
    execute="NO"
    if int(stepToExecute)>int(stepSaved):
        execute="YES"
        print ">>>> The step "+ str(stepToExecute) + " will be carried out.\n"
    return execute
#-------------------------------------------------------------------------
def ConvertLog(bigWigFileList): # There is an error.....It transforms it to LOG fine, but the text file does not save the data correctly
    
    try:
        #Check that the DEC file exists
        fConfig=open(bigWigFileList+".DEC",'r')
        fConfig.close
        # If it continues, the .DEC file DOES exist (Originals)
        
        os.remove(bigWigFileList)   # The file that is in Logarithmic is deleted
        os.rename(bigWigFileList+".DEC", bigWigFileList) # The file that is in Decimal is renamed to be able to convert back to LOG
        
    except IOError:
        print ">>>> The file does not exist "+bigWigFileList+".DEC....Conversion continues"
    
    finally:
        fTemp=open("Temp.txt","w")
        with open(bigWigFileList,"r") as fChromatinFiles_Hela:
            for line in fChromatinFiles_Hela:
                length=len(line)
                indexExtension=line.find(".")
                indexPath=line.find("/")
                
                NameC=line[0:indexPath-1]
                fileName=line[indexPath:length-(length-indexExtension)]
                extensionFile=line[indexExtension+1:length-1]
                
                fileDEC=fileName+"."+extensionFile
                fileLOG=fileName+"_LOG."+extensionFile

                to_Log(fileDEC,fileLOG) # Convert to LOG without deleting any files
                print ">>>> File: " + fileName + " converted to LOG correctly\n"
                
                fTemp.write(NameC+"\t"+fileLOG+"\n") # Write in the txt file the name of the new converted LOG file
        fTemp.close()
        os.rename(bigWigFileList, bigWigFileList+".DEC")
        os.rename("Temp.txt",bigWigFileList)   
#-------------------------------------------------------------------------        
def to_Log(fileR, fileW):
    # Open fileR as read and fileW as write
    # Converts fileR (Decimal) to fileW (Logarithmic)
    
    os.system("bigWigToBedGraph " + fileR + " temp.bedgraph")
    
    fileR="temp.bedgraph"
    
    ip = open(fileR,"r")
    op = open(fileW,"w")

    for line in ip:
        #print "read line from file: " +line
        if "#" in line:
            continue
        
        fields = line.strip().split("\t")
        signal = float(fields[3])
        
        if (signal == 0.0):
            op.write(fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t"+str(0.0)+ "\n") 
        else:
            signal = math.log(signal+0.01,10)
            op.write(fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t"+str(signal)+ "\n")

    ip.close()
    op.close()
    deleteFile("temp.bedgraph")
#-------------------------------------------------------------------------
def my_ReadFiles(fileRead, dataName): # Reads all the files and returns 1 or 2 data depending on the file read
    
    try:
        ip = open(fileRead, "r")
        if dataName=="chromosomeFiles": #Reading chromosome names and lengths
            chrSize = OrderedDict()
            for line in ip:
                name = line.split("\t")[0] # creates a list of the elements that are in a line of the file separated by tab
                size = int(line.strip().split("\t")[1])
                chrSize[name] = size            
            return1=chrSize
            return2=""
            
        if dataName=="metaproFiles": #Reading metaprofiles  
            metaprofiles = OrderedDict()
            for line in ip:
                signalType = line.split("\t")[0]
                filename = line.strip().split("\t")[1]
                metaprofiles[signalType] = readMetaProfile(filename)
            return1=metaprofiles
            return2=""
        
        if dataName=="peakFiles":  #Reading peak file list       
            peakFiles = OrderedDict()
            for line in ip:
                peakFiles[line.split("\t")[0]] = line.strip().split("\t")[1]
            return1=peakFiles
            return2=""

        if dataName=="bigwigFiles": #Reading bigwig files
            bigWigList = OrderedDict()
            bigWigFiles = OrderedDict()
            for line in ip:
                signalType = line.split("\t")[0]
                filename = line.strip().split("\t")[1]
                bigWigList[signalType] =  metaseq.genomic_signal(filename, "bigWig")
                bigWigFiles[signalType] = filename            
            return1=bigWigList
            return2=bigWigFiles

        ip.close()
    except:
        sys.stderr.write("Cannot open " + fileRead + "\n")
        sys.exit(1)

    return return1, return2
#-------------------------------------------------------------------------

#******************** PROGRAM FUNCTIONS MOVED TO THE MYFUNCTIONS FILE
def readMetaProfile(metaProfileFile): # Necessary function. Do not comment or delete
    metaProfile = []
    try:
        ip = open(metaProfileFile, "r")
    except:
        sys.stderr.write("ERROR: could not open " + metaProfileFile +"\n")
        sys.exit(1)

    for line in ip:
        try:
            metaProfile.append(float(line.rstrip()))
        except:
            sys.stderr.write("ERROR: could not convert this line in " + metaProfileFile + "\n")
    ip.close()

    return metaProfile
#-------------------------------------------------------------------------
def sameMarks(bigWigList, metaprofiles):
    #Check if all marks with signals have associated metaprofiles. 
    #Note that all metaprofiles need not have inputs.
    for currMark in bigWigList:
        if currMark not in metaprofiles:
            return False

    return True

def getRelevantSignals(bigWigFiles, bigWigList, currRegions, currSignal):
    binSize = 25
    signalList = []
    
    for currFeature in currRegions:
        print currFeature
        numBP = currFeature.end - currFeature.start
        numBins = numBP/binSize
        if numBP % binSize != 0:
            numBins += 1
            numBP += 25 - (numBP % binSize)
        signal = bigWigList[currSignal].array([currFeature], bins=numBins)
        signalList.append(signal[0])

    del bigWigList[currSignal] 
    bigWigList[currSignal] =  metaseq.genomic_signal(bigWigFiles[currSignal], "bigWig")
    return signalList, currRegions

def scoreWithMatchedFilter(signalList, regionList, metaprofile, currWidth, currSignal, opPrefix):
    finalWidth = 10.0/(len(metaprofile) - 10) * currWidth + currWidth
    binSize = 25
    idx = 0
    numPoints = len(metaprofile)
    profile = np.array(metaprofile[::-1])
    op = open(opPrefix + "_" + currSignal + "_" + str(currWidth) + ".bed", "w")
    print regionList
    currChr = regionList[0].chrom

    for currSignal in signalList:
        #Checking if current region is larger than width of metafilter being tested
        if regionList[idx].end - regionList[idx].start < finalWidth:
            idx += 1
            continue

        #Create interpolation with right x and y ranges
        y = signalList[idx]
        end = np.size(y) * 25 - 12.5
        x = np.linspace(12.5, end, num=np.size(y))
        f = interpolate.interp1d(x, y)
        intervalStart = regionList[idx].start

        #Calculating matched filter scores for each possible interval of current size within the region
        currStart = 12.5
        currEnd = currStart + finalWidth
        while currEnd < end:
            #Interpolating the signal for current interval so that MF score can be calculated
            xnew = np.linspace(currStart, currEnd, num=numPoints)
            try:
                ynew = f(xnew)
            except:
                print "Issue", xnew, "not within", x
                sys.exit()
            ynew = ynew.T
            
            #If all values within the current interval is 0, then don't calculate MF score
            if abs(np.max(ynew)) < 0.001:
                currStart += binSize
                currEnd += binSize
                continue

            score = np.dot(profile, ynew)
            op.write(currChr + "\t" + str(int(round(currStart + intervalStart))) + "\t" + str(int(round(currEnd + intervalStart))) + "\t" + str(score) + "\n")
            #print otherScores, bestStart
            currStart += binSize
            currEnd += binSize
        idx += 1

    op.close()
    return

def getProbabilityDistribution(scores):
    y, binEdges = np.histogram(scores, bins=100)
    x = np.zeros(100)
    y = y/float(len(scores))
    for i in range(0, len(binEdges) - 1):
        x[i] = (binEdges[i] + binEdges[i + 1])/2
        x = np.array(x)
    return x, y

def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def removeExtremePoints(scores, norm, median, stddev):
    n = len(scores)
    n = 20
    i = 0
    while scipy.stats.norm.sf(i) > 1.0/len(scores):
        i += 0.1
    extreme = median + i * stddev
    reducedScores = [currScore for currScore in scores if currScore <= extreme]
    return reducedScores


def calculateNormalizedScores(peakFiles, currWidth, opPrefix, pp):
    
    normDist = OrderedDict()
    meanDist = OrderedDict()
    stdDist = OrderedDict()
    ipFiles = OrderedDict()
    header = "#chrom\tstart\tend"
    if "H3K27ac" in peakFiles:
        header += "\t" + "H3K27ac"
    for currSignal in peakFiles:
        if currSignal != "H3K27ac":
            header += "\t" + currSignal
        
        #Subtracting peak region 
        os.system("awk '{print $1 \"\\t\" $2 \"\\t\" $3}' " + peakFiles[currSignal] + \
			" | intersectBed -a " + opPrefix + "_" + currSignal + "_" + str(currWidth) + ".bed -b " + \
			"- -v > " + opPrefix + "_" + currSignal + "_" + str(currWidth) + "_2.bed")
        
        #Now finding distribution of MF scores
        ip = open(opPrefix + "_" + currSignal + "_" + str(currWidth) + "_2.bed", "r")
        scores = []
        line=ip.readline(); # <<<------- WE HAVE ADDED THIS BECAUSE YOU SHOULD READ THE FILE
        for line in ip:
            scores.append(float(line.strip().split("\t")[3]))
        ip.close()
        #print scores
        x, y = getProbabilityDistribution(scores) 

        #Initial curve fit with median and std deviation of all scores
        median = np.median(scores)
        stddev = np.std(scores)
        yest = gauss_function(x,np.max(y), median, stddev)
        P.figure()
        P.plot(x, y,'b+:',label='MF scores')
        P.plot(x, yest,'ro:',label='fit')
        print np.max(y), median, stddev

        #Removing extreme points and refit
        scores = removeExtremePoints(scores, np.max(y), median, stddev)
        x, y = getProbabilityDistribution(scores)
        median = np.median(scores)
        stddev = np.std(scores)
        try:
            popt,pcov = curve_fit(gauss_function, x, y, p0=[np.max(y), median, stddev])
        except:
            popt = [np.max(y), median, stddev]
        yest2 = gauss_function(x,*popt)
        P.plot(x, y,'b*:',label='filtered MF scores')
        P.plot(x, yest2,'g*:',label='fit2')
        P.title("width of MF in " + currSignal + " = " + str(currWidth) + "bp")
        P.legend()
        P.savefig(pp, format="pdf")
        normDist[currSignal] = popt[0]
        meanDist[currSignal] = popt[1]
        stdDist[currSignal] = popt[2]
        print popt
        ipFiles[currSignal] = open(opPrefix + "_" + currSignal + "_" + str(currWidth) + ".bed", "r")
    header += "\n"

    if "H3K27ac" in peakFiles:
        opFile1 = open(opPrefix + "_pValue_" + str(currWidth) + ".bed", "w")
        opFile2 = open(opPrefix + "_zscore_" + str(currWidth) + ".bed", "w")
        opFile1.write(header)
        opFile2.write(header)
        lines = OrderedDict()
        for line in ipFiles["H3K27ac"]:
            fields = line.strip().split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            z_score = (float(fields[3]) - meanDist["H3K27ac"])/stdDist["H3K27ac"]
            pValue = scipy.stats.norm.sf(z_score)
            opFile1.write(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + str(pValue))
            opFile2.write(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + str(z_score))
            
            for currSignal in peakFiles:
                if currSignal == "H3K27ac":
                    continue
                if currSignal not in lines:
                    lines[currSignal] = ipFiles[currSignal].readline()
                while True:
                    if lines[currSignal] == "":
                        opFile1.write("\t" + str(0.5))
                        opFile2.write("\t" + str(0.0))
                        break
                    fields = lines[currSignal].strip().split("\t")
                    newStart = int(fields[1])
                    if newStart < start:
                        lines[currSignal] = ipFiles[currSignal].readline()
                    elif newStart > start:
                        opFile1.write("\t" + str(0.5))
                        opFile2.write("\t" + str(0.0))
                        break
                    else:
                        currScore = float(fields[3])
                        #print currScore, float(meanDist[currSignal])
                        #print float(stdDist[currSignal])
                        z_score = (currScore - meanDist[currSignal])/stdDist[currSignal]
                        pValue = scipy.stats.norm.sf(z_score)
                        opFile1.write("\t" + str(pValue))
                        opFile2.write("\t" + str(z_score))
                        break
            opFile1.write("\n")
            opFile2.write("\n")
        opFile1.close()
        opFile2.close()

    return

def readHeader(ip):
    line = ip.readline()
    fields = line.strip().split("\t")
    return fields

def getPositivesMF(opPrefix, peakFiles, currWidth):
    ip = open(opPrefix + "_pValue_" + str(currWidth) + ".bed", "r")
    #Get header
    header = readHeader(ip)
    
    #Open output files and choose correct column
    op = OrderedDict()
    columnNumber = OrderedDict()
    for currSignal in peakFiles:
        columnNumber[currSignal] = header.index(currSignal)
        op[currSignal] = open(opPrefix + "_" + currSignal + "_MFpositives_" + str(currWidth) + ".bed", "w")
    
    #Check what passes the p-Value in the MF for corresponding signal
    for line in ip:
        fields = line.strip().split("\t")
        for currSignal in peakFiles:
            if float(fields[columnNumber[currSignal]]) <= 0.001:
                op[currSignal].write(fields[0] + "\t" + fields[1] + "\t" + fields[2] + "\t" + fields[columnNumber[currSignal]] + "\n")
    
    for currSignal in peakFiles:
        op[currSignal].close()
    
    ip.close()
    return

def smoothSignals(signal, win=10):
    length = len(signal)
    newSignal = [0] * length
    for idx in range(10, length-10):
        newSignal[idx] = np.mean(signal[idx-10:idx+11])
    return newSignal

def getMaxima(signal):
    length = len(signal)
    idx = 0

    maximaIndices = []
    for idx in range(10, length - 10):
        if signal[idx] > max(signal[idx-3:idx]) and signal[idx] >= max(signal[idx+1:idx+4]):
            maximaIndices.append(idx)

    return maximaIndices

def findPossiblePairings(maximaIndices, binSize, minWidth, maxWidth):
    minSep = minWidth/binSize
    maxSep = maxWidth/binSize
    allPairings = []
    for idx1 in range(0, len(maximaIndices)):
        for idx2 in range(idx1 + 1, len(maximaIndices)):
            if ((maximaIndices[idx2] - maximaIndices[idx1]) < minSep):
                continue
            if ((maximaIndices[idx2] - maximaIndices[idx1]) <= maxSep):
                currPairing = [maximaIndices[idx1], maximaIndices[idx2]]
                allPairings.append(currPairing)
            if ((maximaIndices[idx2] - maximaIndices[idx1]) > maxSep):
                break
    return allPairings

def filterPossibleMatches(bigWigFiles, bigWigList, metaprofile, opPrefix, currSignal, minWidth, maxWidth):
    
    currRegions = pbt.BedTool(opPrefix + "_" + currSignal + "_MFpositives.bed")
    
    signalList = getSignals(bigWigFiles, bigWigList, currRegions, currSignal)
    del bigWigList[currSignal] 
    bigWigList[currSignal] =  metaseq.genomic_signal(bigWigFiles[currSignal], "bigWig")
    print len(signalList)
    binSize = 25

    op = open(opPrefix + "_" + currSignal + "_tentPositives.bed", "w")
    newWidth = minWidth
    ipFiles = OrderedDict()
    while newWidth <= maxWidth:
        ipFiles[newWidth] = open(opPrefix + "_pValue_" + str(newWidth) + ".bed", "r")
        newWidth += binSize 
    
    for idx in range(0, len(signalList)):
        if idx % 1000 ==0:
            print idx
        region = currRegions[idx]
        smoothedSignal = smoothSignals(signalList[idx], win=10)
        maximaIndices = getMaxima(smoothedSignal)
        if len(maximaIndices) < 2:
            continue
        pairings = findPossiblePairings(maximaIndices, binSize, minWidth, maxWidth)
        for currPairing in pairings:
            newWidth = binSize * (currPairing[1] - currPairing[0])
            if newWidth < minWidth or newWidth > maxWidth:
                continue
            currStart = np.floor((region.start - 500) + (currPairing[0] * binSize) - 5.0/(len(metaprofile) - 10) * newWidth) + 25
            currEnd = np.floor(region.start - 500 + currPairing[1] * binSize + 5.0/(len(metaprofile) - 10) * newWidth)
            currChr = region.chrom
            
            for line in ipFiles[newWidth]:
                if "#" in line:
                    headerFields = line.strip().split("\t")
                    currSignalIdx = headerFields.index(currSignal)
                    continue
                fields = line.strip().split("\t")
                start = int(fields[1])
                chrom = fields[0]
                end = int(fields[2])
                if start > currStart + binSize:
                    break
                if start <= currStart and end >= currEnd:
                    #accepting this enhancer if it passes cutoff
                    if float(fields[currSignalIdx]) <= 0.001:
                        op.write(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + str(fields[currSignalIdx]) + "\t" + str(newWidth) + "\n")
                        break
    newWidth = minWidth
    while newWidth <= maxWidth:
        ipFiles[newWidth].close()
        newWidth += binSize 
    op.close()
    os.system("sortBed -i " + opPrefix + "_" + currSignal + "_tentPositives.bed > " + opPrefix + "_" + currSignal + "_tentPositives2.bed")
    ip = open(opPrefix + "_" + currSignal + "_tentPositives2.bed" ,"r")
    op = open(opPrefix + "_" + currSignal + "_finalPositives.bed", "w")
    currStart = currEnd = currPvalue = currWidth = 0
    currChrom = ""
    for line in ip:
        fields = line.strip().split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        pValue = float(fields[3])
        width = int(fields[-1])
        if chrom == currChrom and start <= currStart and end >= currEnd:
            currStart = start 
            currEnd = end
            currPvalue = pValue
            currWidth = width
            currChrom = chrom
        elif chrom == currChrom and currStart <= start and currEnd >= end:
            continue
        else:
            if currChrom != "":
                op.write(currChrom + "\t" + str(currStart) + "\t" + str(currEnd) + "\t" + str(currPvalue) + "\t" + str(currWidth) + "\n")
            currStart = start 
            currEnd = end
            currPvalue = pValue
            currWidth = width
            currChrom = chrom
    op.close()
    ip.close()
    os.system("rm " + opPrefix + "_" + currSignal + "_tentPositives2.bed " + opPrefix + "_" + currSignal + "_tentPositives.bed")
    return

def getGaussianFits(bigWigList, opPrefix):
    #Finding distribution of MF scores
    normDist = OrderedDict()
    meanDist = OrderedDict()
    stdDist = OrderedDict()
    for currSignal in bigWigList:
        print currSignal
        currWidth = 350
        normDist[currSignal] = {}
        meanDist[currSignal] = {}
        stdDist[currSignal] = {}
        while currWidth <= 1100:
            ip = open(opPrefix + "_" + currSignal + "_" + str(currWidth) + "_2.bed", "r")
            scores = []
            for line in ip:
                scores.append(float(line.strip().split("\t")[3]))
            ip.close()
            #print scores
            x, y = getProbabilityDistribution(scores) 

            #Initial curve fit with median and std deviation of all scores
            median = np.median(scores)
            stddev = np.std(scores)

            #Removing extreme points and refit
            scores = removeExtremePoints(scores, np.max(y), median, stddev)
            x, y = getProbabilityDistribution(scores)
            median = np.median(scores)
            stddev = np.std(scores)
            try:
                popt,pcov = curve_fit(gauss_function, x, y, p0=[np.max(y), median, stddev])
            except:
                popt = [np.max(y), median, stddev]
            yest2 = gauss_function(x,*popt)
            normDist[currSignal][currWidth] = popt[0]
            meanDist[currSignal][currWidth] = popt[1]
            stdDist[currSignal][currWidth] = popt[2]
            currWidth += 25

    return normDist, meanDist, stdDist

def getSignals(bigWigFiles, bigWigList, currRegions, currSignal):
    binSize = 25
    idx = 0 
    signalList = []
    for currFeature in currRegions:
        if currFeature.start  < 500:
            currFeature.start = 500
        newFeature = pbt.BedTool(currFeature.chrom + " " + str(currFeature.start - 500) + " " + str(currFeature.end + 500), from_string=True)[0]
        numBP = newFeature.end - newFeature.start
        numBins = numBP/binSize
        if numBP % binSize != 0:
            numBins += 1
            newFeature.end += 25 - (newFeature.end % binSize)
        signal = bigWigList[currSignal].array([newFeature], bins=numBins)
        signalList.append(signal[0])
        idx += 1
        if idx%1000 == 0:
            del bigWigList[currSignal] 
            bigWigList[currSignal] =  metaseq.genomic_signal(bigWigFiles[currSignal], "bigWig")

    return signalList    
# ***************************************************************************************************************
