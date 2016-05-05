#!/usr/bin/python
#bionetworks.py
#last update : 14 Aug 2014

__author__ = '''Hyunju Kim'''


import networkx as nx
import os
import sys
import random as ran
from math import log
from optparse import OptionParser, OptionGroup
from scipy import *
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
#import bionet



################# begin : build_historyList ############################
# *** generate decimal state's index for k-previous states with a given list, aList ***#
def build_historyList(aList, historyLength):
    historyList = []
    historyUnit = aList[:historyLength]
    #print "Unit", historyUnit
    aList[:historyLength] = []
    
    #print "history-Unit", aList
    
    historyState = 0
    for s in range(historyLength):
        historyState += historyUnit[s] * power(2, s)
    #print "unit", historyUnit[s]
    #print historyState
    historyList.append(historyState)

    for x in aList:
        historyState = historyState / 2 + x * power(2, historyLength - 1)
        historyList.append(historyState)
        #print x, historyState
    return historyList
################# end : build_historyList ########################


################## begin : comAI ########################
def comAI(timeSeries, historyLength):

    aList = list(timeSeries)
    historyList = build_historyList(aList, historyLength) # aList becomes aList[historyLength:] after historyList function
    
    #*** declare dic for distribution to compute Active Information ***#
    count_currState_hiState = defaultdict(int)
    count_hiState = defaultdict(int)
    count_currState = defaultdict(int)
    #print "test", len(aList)
    #*** obtain the distribution for each pattern to compute Active Information ***#
    for s in range(len(aList)):
        #print s
        count_currState_hiState[(aList[s], historyList[s])] += 1
        count_hiState[historyList[s]] += 1
        count_currState[aList[s]] += 1

    #*** obtain the distribution for each pattern to compute Active Information ***#
    AI = 0
    for s in range(len(aList)):
        prob_currState_hiState = float(count_currState_hiState[(aList[s], historyList[s])]) / float(len(aList))
        #print "joint", prob_currState_hiState
        prob_hiState = float(count_hiState[historyList[s]]) / float(len(aList))
        #print "his", prob_hiState, historyList[s]
        prob_currState = float(count_currState[aList[s]]) / float(len(aList))
        #print "curr", prob_currState
        AI = AI + log( prob_currState_hiState / ( prob_currState * prob_hiState)) / log(2.0) # since the summation is over not all possible pattern of currState_hiState but aList, there is no prob_currState_hiState multiplied by the log term.

    AI = AI / float(len(aList))
    return AI
################## end : comAI ########################


################## begin : comTE ########################
def comTE(timeSeriesA, timeSeriesB, historyLength):
    
    sourList = list(timeSeriesA)
    tarList = list(timeSeriesB)
    historyList = build_historyList(tarList, historyLength) # tarList becomes tarList[historyLength:] after historyList function
    sourList[:historyLength - 1] = [] # sourList becomes sourList[historyLength-1:]
    #print timeSeriesA
    #print timeSeriesB
    #print tarList
    #print sourList
    #*** declare dic for distribution to compute Transfer Entropy ***#
    count_tarCurrState_tarHiState_sourPrevState = defaultdict(int)
    count_tarCurrState_tarHiState = defaultdict(int)
    count_tarHiState_sourPrevState = defaultdict(int)
    count_tarHiState = defaultdict(int)
    
    #*** obtain the distribution for each pattern to compute Transfer Entropy ***#
    for s in range(len(tarList)):
        #        print "len s", len(sourList)
        #        print "len t", len(tarList)
#        print s
#        print "s", sourList[s]
#        print "h", historyList[s]
#        print "t", tarList[s]
        count_tarCurrState_tarHiState_sourPrevState[(tarList[s], historyList[s], sourList[s])] += 1
        count_tarCurrState_tarHiState[(tarList[s], historyList[s])] += 1
        count_tarHiState_sourPrevState[(historyList[s], sourList[s])] += 1
        count_tarHiState[historyList[s]] += 1
    
    #*** obtain the distribution for each pattern to compute Active Information ***#
    TE = 0
    for s in range(len(tarList)):
        prob_tarCurrState_tarHiState_sourPrevState = float(count_tarCurrState_tarHiState_sourPrevState[(tarList[s], historyList[s], sourList[s])]) / float(len(tarList))
        prob_tarCurrState_tarHiState = float(count_tarCurrState_tarHiState[(tarList[s], historyList[s])]) / float(len(tarList))
        prob_tarHiState_sourPrevState = float(count_tarHiState_sourPrevState[(historyList[s], sourList[s])]) / float(len(tarList))
        prob_tarHiState = float(count_tarHiState[historyList[s]]) / float(len(tarList))
        TE = TE + log( (prob_tarCurrState_tarHiState_sourPrevState * prob_tarHiState) / ( prob_tarHiState_sourPrevState * prob_tarCurrState_tarHiState)) / log(2.0) # since the summation is over not all possible pattern of tarCurrState_tarHiState_sourPrevState but tarList, there is no prob_tarCurrState_tarHiState_sourPrevState multiplied by the log term.
    TE = TE / float(len(tarList))
    

    return TE
################## end : comAI ########################

#aList = [0,0,1, 0, 0, 1, 1, 1]
#print comAI(aList, 3)
#tarList = [1,1,1,0,1,1,0,1,0,0,1,1,0,0,1,1,1,0,1,1,0,1]
#sourList = [1,0,1,0,1,0,1,0,1,1,0,1,1,0,1,1,0,1,1,0,0,1]
#print comTE(sourList, tarList, 1)

