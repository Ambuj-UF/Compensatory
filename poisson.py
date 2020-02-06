# -*- coding: utf-8 -*-
# Copyright 2019 by Ambuj Kumar, Iowa State University.
# All rights reserved.
#
# Bug reports welcome: ambuj@iastate.edu

import math
from Bio import SeqIO

def _poisson(distance, long):
    
    pdist = 0
    if 1 - float(distance)/long != 0:
        pdist = -math.log(1 - float(distance)/long)
    else:
        pdist = 0

    return pdist


def _poisson_dist(seq1, seq2):
    
    """Poisson distance"""
    
    dist = 0
    gap = 0
    for i, amino in enumerate(seq1):
        if seq1[i] == seq2[i] and seq1[i] != '-' and seq2[i] != '-':
            dist = dist + 1
        if seq1[i] == '-' and seq2[i] == '-':
            gap = gap + 1

    return _poisson(dist, len(seq1) - gap)
    
    
def calc_distance(id1, id2, seq1, seq2):
    with open("raw_distance.txt", "a+") as fp:
        dObj = _poisson_dist(seq1, seq2)
        fp.write("%s\t%s\t%s\n" %(id1, id2, dObj))
    return True









