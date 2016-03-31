#!/usr/bin/python
#attractor-analysis.py
#last update: 17 DEC 2015

__author__ = '''Hyunju Kim'''


import os
import sys
import numpy as np
import networkx as nx
import re

################# BEGIN: read_network_from_file(EDGE_FILE, NODE_FILE) ########################
def read_network_from_file(EDGE_FILE, NODE_FILE):

    '''
     Arguments:
               1. EDGE_FILE => source  target  weight
               2. NODE_FILE => node  threshold
     Return:
               1. net => a network
    '''

    #### read edge list with its weight from EDGE_FILE
    net = nx.read_edgelist(EDGE_FILE, create_using=nx.DiGraph(), nodetype=str, data=(('weight', float),))

    #### read node list with its threshold from NODE_FILE
    for line in open(NODE_FILE, 'r').readlines():
        items = [x.strip() for x in re.split('\s+', line.rstrip())]
        #items = [x.strip() for x in line.rstrip().split('\t')] #H's style with tabs
        if line[0] == '#' or line=='':
            continue
        net.add_node(items[0], threshold=float(items[1]))

    return net
################# END: read_network_from_file(EDGE_FILE, NODE_FILE) ########################

################# BEGIN: build_nodes_list(NODE_FILE) ########################
def build_nodes_list(NODE_FILE):

    '''
        Arguments:
                1. NODE_FILE => node  threshold
        Return:
                1. nodes_list => a list of nodes in the fixed order
    '''
    nodes_list = []
    for line in open(NODE_FILE, 'r').readlines():
        items = [x.strip() for x in re.split('\s+', line.rstrip())]
        #items = [x.strip() for x in line.rstrip().split('\t')] #H's style with tabs
        if line[0] == '#' or line=='':
            continue
        nodes_list.append(items[0])

    return nodes_list
################# END: build_nodes_list(NODE_FILE) ########################

################# BEGIN: read_init_from_file(BIO_INIT_FILE) ########################
def read_init_from_file(BIO_INIT_FILE):

    '''
        Arguments:
        1. BIO_INIT_FILE => node  initial-state-for-bio-seq #user must input data path to bio_init_file in time-evol_eleg.py
        Return:
        1. bio_initStates => dict { node: inital-binary-value, ... }
        '''
    bio_initStates = {}
    for line in open(BIO_INIT_FILE, 'r').readlines():
        items = [x.strip() for x in line.rstrip().split('\t')]
        if line[0] == '#' or line=='':
            continue
        bio_initStates[items[0]] = int(items[1])

    return bio_initStates
################# END: read_init_from_file(BIO_INIT_FILE) ########################



def main():
    print "input_net_eleg module is the main code."
    EDGE_FILE = 'C:\Users\Kelle Dhein\C.-elegans\example\SES591_SampleCode\data\elegans\elegans-net-edges-new-names.dat'
    NODE_FILE = 'C:\Users\Kelle Dhein\C.-elegans\example\SES591_SampleCode\data\elegans\elegans-net-nodes-new-names.dat'

    net = read_network_from_file(EDGE_FILE, NODE_FILE)
    nodes_list = build_nodes_list(NODE_FILE)

    for u, v in net.edges():
        print u, v, net[u][v]['weight']

    print nodes_list



if __name__=='__main__':
    main()
