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
def comAI(timeSeriesNode, historyLength, nodes_list):

    count_currState_hiState = defaultdict(int)
    count_hiState = defaultdict(int)
    count_currState = defaultdict(int)
    for si in range(0, pow(2,len(nodes_list))):
        aList = list(timeSeriesNode[si])
        historyList = build_historyList(aList, historyLength) # aList becomes aList[historyLength:] after historyList function
        #print "test", len(aList)
        #*** obtain the distribution for each pattern to compute Active Information ***#
        for s in range(len(aList)):
            #print s
            count_currState_hiState[(aList[s], historyList[s])] += 1
            count_hiState[historyList[s]] += 1
            count_currState[aList[s]] += 1

    #        print count_currState_hiState
    #        print count_currState
    #        print count_hiState
    AI = 0
    for si in range(0, pow(2,len(nodes_list))):
        aList = list(timeSeriesNode[si])
        historyList = build_historyList(aList, historyLength) # aList becomes aList[historyLength:] after historyList function
        sampleLength = len(aList) * pow(2,len(nodes_list))
        for s in range(len(aList)):
            prob_currState_hiState = float(count_currState_hiState[(aList[s], historyList[s])]) / float(sampleLength)
            #print "joint", prob_currState_hiState
            prob_hiState = float(count_hiState[historyList[s]]) / float(sampleLength)
            #print "his", prob_hiState, historyList[s]
            prob_currState = float(count_currState[aList[s]]) / float(sampleLength)
            #print "curr", prob_currState
            AI = AI + log( prob_currState_hiState / ( prob_currState * prob_hiState)) / log(2.0) # since the summation is over not all possible pattern of currState_hiState
    AI = AI / float(sampleLength)
    return AI
################## end : comAI ########################


################## begin : comTE ########################
def comTE(timeSeriesNodeA, timeSeriesNodeB, historyLength, nodes_list):
    
    #*** declare dic for distribution to compute Transfer Entropy ***#
    count_tarCurrState_tarHiState_sourPrevState = defaultdict(int)
    count_tarCurrState_tarHiState = defaultdict(int)
    count_tarHiState_sourPrevState = defaultdict(int)
    count_tarHiState = defaultdict(int)
    
    for si in range(0, pow(2,len(nodes_list))):
        sourList = list(timeSeriesNodeA[si])
        tarList = list(timeSeriesNodeB[si])
        historyList = build_historyList(tarList, historyLength) # tarList becomes tarList[historyLength:] after historyList function
        sourList[:historyLength - 1] = [] # sourList becomes sourList[historyLength-1:]
        #*** obtain the distribution for each pattern to compute Transfer Entropy ***#
        for s in range(len(tarList)):
            count_tarCurrState_tarHiState_sourPrevState[(tarList[s], historyList[s], sourList[s])] += 1
            count_tarCurrState_tarHiState[(tarList[s], historyList[s])] += 1
            count_tarHiState_sourPrevState[(historyList[s], sourList[s])] += 1
            count_tarHiState[historyList[s]] += 1

#    print count_tarCurrState_tarHiState_sourPrevState
#    print count_tarCurrState_tarHiState
#    print count_tarHiState_sourPrevState
#    print count_tarHiState
    #*** obtain the distribution for each pattern to compute Active Information ***#
    TE = 0
    for si in range(0, pow(2,len(nodes_list))):
        sourList = list(timeSeriesNodeA[si])
        tarList = list(timeSeriesNodeB[si])
        historyList = build_historyList(tarList, historyLength) # tarList becomes tarList[historyLength:] after historyList function
        sourList[:historyLength - 1] = [] # sourList becomes sourList[historyLength-1:]
        sampleLength = len(tarList) * pow(2,len(nodes_list))
        for s in range(len(tarList)):
            prob_tarCurrState_tarHiState_sourPrevState = float(count_tarCurrState_tarHiState_sourPrevState[(tarList[s], historyList[s], sourList[s])]) / float(sampleLength)
            prob_tarCurrState_tarHiState = float(count_tarCurrState_tarHiState[(tarList[s], historyList[s])]) / float(sampleLength)
            prob_tarHiState_sourPrevState = float(count_tarHiState_sourPrevState[(historyList[s], sourList[s])]) / float(sampleLength)
            prob_tarHiState = float(count_tarHiState[historyList[s]]) / float(sampleLength)
            TE = TE + log( (prob_tarCurrState_tarHiState_sourPrevState * prob_tarHiState) / ( prob_tarHiState_sourPrevState * prob_tarCurrState_tarHiState)) / log(2.0) # since the summation is over not all possible pattern of tarCurrState_tarHiState_sourPrevState but tarList, there is no prob_tarCurrState_tarHiState_sourPrevState multiplied by the log term.
    TE = TE / float(sampleLength)
    return TE
################## end : comAI ########################

#aList = [0,0,1, 0, 0, 1, 1, 1]
#print comAI(aList, 3)
#tarList = [1,1,1,0,1,1,0,1,0,0,1,1,0,0,1,1,1,0,1,1,0,1]
#sourList = [1,0,1,0,1,0,1,0,1,1,0,1,1,0,1,1,0,1,1,0,0,1]
#print comTE(sourList, tarList, 1)

