#!/usr/bin/env python2.7
# Fecha 01/07/2022....

import matplotlib
matplotlib.use('Agg')
import sys
import os
from collections import OrderedDict
import operator 

import pybedtools as pbt
import metaseq
import numpy as np
from scipy import interpolate
from matplotlib.backends.backend_pdf import PdfPages
import pylab as P
import scipy
import training
import testData
import math

#***************  Code to eliminate Warnings ***************************
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
#****************************************************************************

from myFunctions import *

#-----------------------------------------------------------
def main(bigWigFileList, metaProfileList, chrNameFile, peakFileList, trainingPositives, trainingNegatives, opPrefix):
    os.system("clear")
    print "**************************************************************************************\n"
    print "*             GRADO DE INGENIERÍA ELECTRÓNICA INDUSTRIAL Y AUTOMÁTICA                *\n"
    print "*                                                                                    *\n"
    print "*       Desarrollo de un modelo de inteligencia artificial para el estudio de        *\n"
    print "*                               problemas epigenéticos                               *\n"
    print "*                                                                                    *\n"
    print "*  Fechas de desarrollo: 26/11/2021 - 03/07/2022          Curso académico 2021/2022  *\n"
    print "*                                                                                    *\n"
    print "*  Alumno:    Daniel Rodriguez Martín                                                *\n"
    print "*  Profesor:  Dr. José Mnauel Gilpérez Aguilar                                       *\n"
    print "**************************************************************************************\n"    
    print "\n"
    print "\n"
    print "--------------------------------------------------------------------------------------\n"
    print "STEP 1: *********** CONVERT bigWig FILES TO dB and READING FILES *********************\n"
    print "                 Tenga paciencia. Este paso dura 24h aproximadamente\n
    print "--------------------------------------------------------------------------------------\n"
    print "\n"
    doStep=AskStep(1) # We ask the first time, if we say YES... all the steps will be done without asking again
    # If we say NO, it is because he had already done it (he will skip the first step) and he will not ask the rest again
    if (doStep=="SI"):
        print ">>>> Log conversion begin !\n"
        ConvertLog(bigWigFileList)  # Added to convert bigWig files to dB  
        WriteStep(1) # writes that step 1 has been done
        
    # ***************************************  STARTING THE READING OF THE FILES  ************************************
    # This should always be done at the beginning
    chromosomes, myData=my_ReadFiles(chrNameFile,"chromosomeFiles") #Reading chromosome names and lengths
    chrName=  max(chromosomes.iteritems(), key=operator.itemgetter(1))[0]
    metaprofiles, myData=my_ReadFiles(metaProfileList,"metaproFiles") #Reading metaprofiles
    bigWigList, bigWigFiles = my_ReadFiles(bigWigFileList, "bigwigFiles")     #Reading bigwig files
    #Sanity check for signal from marks and metaprofiles
    if not(sameMarks(bigWigList, metaprofiles)):
        print ">>>> Keys do not match between input signals and metaprofiles\n"
        sys.exit(1)
    
    peakFiles, myData=my_ReadFiles(peakFileList, "peakFiles")
    
    if not(sameMarks(peakFiles, metaprofiles)):
        print ">>>> Keys do not match between input peaks and metaprofiles\n"
        sys.exit(1)

    print">>>> Reading files successful!\n"
    # ***************************************  END OF THE READING OF THE FILES  ************************************
    print "\n"
    print "--------------------------------------------------------------------------------------\n"
    print "STEP 2: ******************  CONVERTING BIGWIG TO BEDGRAPH  ***************************\n"
    print "--------------------------------------------------------------------------------------\n"
    print "\n"
    pp = PdfPages(opPrefix + "_figures.pdf") # 
    signalList = OrderedDict()
    regionList = OrderedDict()

    print ">>>> Identifying regions"
    bigWigList, bigWigFiles = my_ReadFiles(bigWigFileList+".DEC", "bigwigFiles")     #Reading bigwig files
    #*********************************************************************************************
    
    doStep=AskStep(2)
    if (doStep=="YES"):
        if "H3K27ac" in bigWigFiles:
            os.system("bigWigToBedGraph -chrom=" + chrName + " " + bigWigFiles["H3K27ac"] + " temp_H3K27ac.bedgraph")
            os.system("awk '{if ($4 != 0) print $0}' temp_H3K27ac.bedgraph | mergeBed -i - -d 1000 > temp2_H3K27ac.bedgraph")
            os.system("rm temp_H3K27ac.bedgraph")
            print ">>>> Los archivos BedGraph se han creado correctamente\n"
        WriteStep(2) # writes that step 2 has been done
    # **********************************************************************************************
    if "H3K27ac" in bigWigFiles:
        currRegions = pbt.BedTool("temp2_H3K27ac.bedgraph")
    print "\n"
    print "--------------------------------------------------------------------------------------\n"    
    print "STEP 3: ***********************       GET SIGNALS         ****************************\n"
    print "              Have patience. This step lasts approximately 10 days (240h)             \n"
    print "--------------------------------------------------------------------------------------\n"
    print "\n"
    doStep=ExecuteStep(3)
    if (doStep=="yes"):
        for currSignal in bigWigFiles:
            signalList[currSignal], regionList[currSignal] = getRelevantSignals(bigWigFiles, bigWigList, currRegions, currSignal)
        print "Get background MF scores"
        #Calculate matched filter scores for different signals
        for currSignal in bigWigFiles:
            currWidth = 350
            while currWidth <= 1100:
                scoreWithMatchedFilter(signalList[currSignal], regionList[currSignal], metaprofiles[currSignal], currWidth, currSignal, opPrefix)
                currWidth += 25
        del regionList, signalList
        WriteStep(3) 
    #--------------------------------------------------------------------------------
    print "\n"
    print "--------------------------------------------------------------------------------------\n"
    print "STEP 4: **************           NORMALIZING SCORES          *************************\n"
    print "                 Have patience. This step lasts approximately 48 hours.               \n"
    print "--------------------------------------------------------------------------------------\n"
    print "\n"
    doStep=ExecuteStep(4)
    if (doStep=="YES"):      
        currWidth = 350
        while currWidth <= 1100:
             calculateNormalizedScores(peakFiles, currWidth, opPrefix, pp)
             #getPositivesMF(opPrefix, peakFiles, currWidth) #********* COMMENTED ON THE ORIGINAL PROGRAM
             currWidth += 25
        WriteStep(4)
        
    #--------------------------------------------------------------------------------
    # THIS CODE SHOULD ALWAYS BE PERFORMED
    #Reading training data and then training model
    #Please ensure that same features are in positive score file and negative score file (and they are in same order)  
    print ">>>> Training\n"
    positiveFeatures, positiveScores = training.getScores(trainingPositives)
    negativeFeatures, negativeScores = training.getScores(trainingNegatives)
    positiveScores = training.chooseRelevantColumns(positiveScores, positiveFeatures, bigWigFiles)
    negativeScores = training.chooseRelevantColumns(negativeScores, negativeFeatures, bigWigFiles)

    Zpos, Zneg = training.calculateZscores(positiveScores, negativeScores, negativeScores)
    
    trainingScores = Zpos + Zneg
    trainingResults = ([1] * len(Zpos)) + ([0] * len(Zneg))
    SVM_model = training.performSVM(trainingScores, trainingResults)
    normDist, meanDist, stdDist = getGaussianFits(bigWigList, opPrefix) # Esto debemos hacerlo siempre
    #--------------------------------------------------------------------------------
    print "\n"
    print "--------------------------------------------------------------------------------------\n"
    print "STEP 5: **************   GET RELEVANT REGIONS TO PREDICT REGIONS    ******************\n"
    print "                Have patience. This step lasts approximately 48 hours.                \n"
    print "--------------------------------------------------------------------------------------\n"
    print "\n"
    # This point will ALWAYS be done
    op = open(opPrefix + "_testScores.dat" ,"w")
    signalList = OrderedDict()
    regionList = []
    #skip = True
    for currChr in chromosomes:
        print ">>>> performing conversion bigWigToBedGraph...\n"
        os.system("bigWigToBedGraph -chrom=" + currChr + " " + bigWigFiles["H3K27ac"] + " temp_H3K27ac.bedgraph")
        os.system("awk '{if ($4 != 0) print $0}' temp_H3K27ac.bedgraph | mergeBed -i - -d 1000 > temp2_H3K27ac.bedgraph")
        os.system("rm temp_H3K27ac.bedgraph")
        currRegions = pbt.BedTool("temp2_H3K27ac.bedgraph")
        print ">>>> Testing...\n"
        signalList = testData.getTestRegions(bigWigFiles, bigWigList, currRegions)    
        testData.getPossiblePairings(signalList, currRegions, metaprofiles, meanDist, stdDist, op)
    op.close()
    WriteStep(5)
    #--------------------------------------------------------------------------------
    testData.scoreTestData(opPrefix, SVM_model)
    pp.close()  # Esto va ha crear el archivo "MetaProfiles.txt_figures.pdf"
    print "\n"
    print "**************************************************************************************\n"
    print "*                               SUCCESSFULLY FINISHED                                *\n"
    print "*              Modified by:                 Daniel Rodriguez Martín                  *\n"
    print "**************************************************************************************\n"
    print "\n"
    return 0

if __name__ == "__main__":
    if len(sys.argv) != 8:
        sys.stderr.write("Usage: " + sys.argv[0] + " <fileList> <metaProfileList> <chrNameFile> <peakFileList> <positiveScores> <negativeScores> <opPrefix>\n")
        sys.stderr.write("where:\n")
        sys.stderr.write("\t<fileList> is a file with the list of chromatin signals in the format (2 column tab delimited)\n")
        sys.stderr.write("\t\tSignalType\tFilename\n")
        sys.stderr.write("\t<metaProfileList> is a file with the list of metaprofiles in the format (2 column tab delimited)\n")
        sys.stderr.write("\t\tSignalType\tFilename\n")
        sys.stderr.write("\t<chrNameFile> is the chromosome size file\n")
        sys.stderr.write("\t<peakFileList> is the name of the peak file list for each signal type\n")
        sys.stderr.write("\t<positiveScores> is the file containing the scores for all training positives\n")
        sys.stderr.write("\t<negativeScores> is the file containing the scores for all training negatives\n")
        sys.stderr.write("\t<outputPrefix> is the prefix for all chromosome names\n")
        sys.exit()

    bigWigFileList, metaProfileList, chrNameFile, peakFileList, trainingPositives, trainingNegatives, opPrefix = sys.argv[1:]
    main(bigWigFileList, metaProfileList, chrNameFile, peakFileList, trainingPositives, trainingNegatives, opPrefix)