#!/usr/bin/python
#attractor-analysis.py
#last update: 17 DEC 2015

__author__ = '''Hyunju Kim'''


import os
import sys
import numpy as np
import networkx as nx
from collections import OrderedDict

import input_net_new as inet



################# begin: sigmoid_updating ######################
def sigmoid_updating(net, prevState):

    '''
        Arguments:
                1. net
                2. prevState
        Return:
               1. currState
    '''

    currState = {}
    #### compute the current states of nodes in the net ####
    for v in net.nodes():   #v is target node
        #### compute weighted sum for node v over its neighbors u ####
        eSum = 0
        for u in net.predecessors_iter(v): #applies to all nodes targeting v
            w_uv = 1.0*net[u][v]['weight']
            eSum += w_uv * prevState[u] #sums the weights of all edges targeting v for source nodes that were ON in prev state
        #### determine the current state for v as a function of eSum and threshold of v ####
        if eSum < net.node[v]['threshold']:
            currState[v] = 0
        if eSum > net.node[v]['threshold']:
            currState[v] = 1
        if eSum == net.node[v]['threshold']:
            currState[v] = prevState[v]

    return currState

################# end: sigmoid_updating ########################

def main():
    print "updating_rule module is the main code."
    EDGE_FILE = 'C:\Users\Kelle Dhein\C.-elegans\example\SES591_SampleCode\data\elegans\elegans-net-edges-new-names.dat'
    NODE_FILE = 'C:\Users\Kelle Dhein\C.-elegans\example\SES591_SampleCode\data\elegans\elegans-net-nodes-new-names.dat'

    net = inet.read_network_from_file(EDGE_FILE, NODE_FILE)

    #prevState = {'a':0.0, 'b':0.0, 'c':1.0}
    prevState = {'cdk-2/cyclinE':1.0, 'cki-1':1.0, 'cdc-14/fzy-1':1.0, 'fzr-1':1.0, 'cdk-1/cyclinB':0.0, 'lin-35/efl-1/dpl-1':1.0, 'cul-1/lin-23':0.0, 'cdc-25.1':0.0}
    # prevState['cdk-2/cyclinE'] = float(sys.argv[1])
    # prevState['cki-1'] = float(sys.argv[2])
    # prevState['cdc-14/fzy-1'] = float(sys.argv[3])
    # prevState['fzr-1'] = float(sys.argv[4])
    # prevState['cdk-1/cyclinB'] = float(sys.argv[5])
    # prevState['lin-35/efl-1/dpl-1'] = float(sys.argv[6])
    # prevState['cul-1/lin-23'] = float(sys.argv[7])
    # prevState['cdc-25.1'] = float(sys.argv[8])
    print "network state @ previous step", OrderedDict(sorted(prevState.items(), key=lambda t: t[0]))

    currState = sigmoid_updating(net, prevState)
    print "network state @ current step", OrderedDict(sorted(currState.items(), key=lambda t: t[0]))

if __name__=='__main__':
    main()
