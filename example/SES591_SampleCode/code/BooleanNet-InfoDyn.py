#!/usr/bin/python
#bioinfo.py

__author__ = '''Hyunju Kim'''

import sys
import os
import random as ran
from math import log
from optparse import OptionParser, OptionGroup
from scipy import *
import matplotlib.pyplot as plt
import numpy as np
import itertools
from collections import defaultdict
import operator
import draw_plots

import input_net as inet
import updating_rule as ur
import time_evol as tev
import info_dyn as info


def main(args):


    ## to obtain biological sequence for the Fission Yeast Cell-Cycle Net starting from biological inital state
    EDGE_FILE = '../data/fission-net/fission-net-edges.txt'
    NODE_FILE = '../data/fission-net/fission-net-nodes.txt'
    BIO_INIT_FILE = '../data/fission-net/fission-net-bioSeq-initial.txt'
    
    net = inet.read_network_from_file(EDGE_FILE, NODE_FILE)
    nodes_list = inet.build_nodes_list(NODE_FILE)


    #input_file_name1 = 'time-series/%s-step%d-trans0.dat'%(network_index, maxStep)
    #input_file1 = open( input_file_name1, 'r')
    
    Nbr_Initial_States = np.power(2,len(nodes_list))
    maxStep = 20
    Nbr_States = 2
    historyLength = 5
    
    result_ai = open('../results/fission-net/ai-step%d-trans0-h%d.dat'%(maxStep, historyLength),'w')
    result_te = open('../results/fission-net/te-step%d-trans0-h%d.dat'%(maxStep, historyLength),'w')

    timeSeries = tev.time_series(net, nodes_list, Nbr_Initial_States, Nbr_States, MAX_TimeStep=20)


    print 'AI'
    AI = {}
    for n in nodes_list:
        AI[n] = info.compute_AI(timeSeries[n], historyLength, Nbr_Initial_States, Nbr_States)
        result_ai.write('%s\t%f\n'%(n, AI[n]))
        print n, AI[n]
    print 'done AI'


    print 'TE'
    TE =  defaultdict(float)
    for v in nodes_list:
        for n in nodes_list:
            TE[(v, n)] = info.compute_TE(timeSeries[v], timeSeries[n], historyLength, Nbr_Initial_States, Nbr_States)
            result_te.write('%s\t%s\t%f\n'%(v, n,TE[(v, n)] ))
            print v, n, TE[(v, n)]
    print 'done TE'



if __name__=='__main__':
    main(sys.argv)
