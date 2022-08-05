"""
Functions for reading gene sets for performing various analysis for CoRe.
"""

import numpy as np
import pandas as pd
import json

def read_all_interactions(datafile):
    interactions = {}

    rf = open(datafile,'r')
    all_lines = rf.readlines()
    rf.close()

    for l in all_lines:
        ns = l.rstrip('\r\n').split(',')

        interactions[ns[0]] = []

        for nn in ns[1:]:
            interactions[ns[0]].append(nn)

    return interactions

def read_interactions(datafile):
    interactions = {}

    rf = open(datafile,'r')
    all_lines = rf.readlines()
    rf.close()

    for l in all_lines:
        ns = l.rstrip('\r\n').split(',')

        if len(ns)>1:
            interactions[ns[0]] = []

            for nn in ns[1:]:
                interactions[ns[0]].append(nn)

    return interactions

def read_interactions_for_GO(datafile):
    f = open(datafile)
    data_dict = json.load(f)
    f.close()

    data = {}

    for d in data_dict.keys():
        data[d] = pd.DataFrame(data_dict[d])

    return data
