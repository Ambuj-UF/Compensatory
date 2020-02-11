# -*- coding: utf-8 -*-
# Copyright 2019 by Ambuj Kumar, Iowa State University.
# All rights reserved.
#
# Bug reports welcome: ambuj@iastate.edu


"""
    
    
                                    Compensatory position calculation
    
    In Intra and Inter mode, program estimates site coevolution statistics using CAPS method.
    Program also allows user to choose between standard pearson correlation coefficient method
    implimented in CAPS program and Gobels site correlation test algorithm. 
    
    Program allows user to restrict the execution to the selected list of alignment positions
    to reduce the execution runtime in the cases where the user wants to test correlation
    between specific alignment sites.
    
    It can also detect the presence/absence of interatomic interactions, such as the hydrogen bonds,
    at the coevolving sites.
    
    
    
    Based upon 
    
    'CAPS: coevolution analysis using protein sequences'
    Fares MA, McNally D. Bioinformatics. 2006 Nov 15;22(22):2821-2.
    
    Gobel U, Sander C, Schneider R and Valencia A. (1994) Proteins, 18, 309-317.
    
    Manish C.Saraf, Gregory L.Moore and Costas D.Maranas.
    Using multiple sequence correlation analysis to characterize functionally important protein regions.
    Protein Engineering vol. 16 no. 6 pp. 397-406, 2003
    
    
    
    
                "Site Correlation" - Inter and Intra Mode
    
    BLOSUM, KMAT and other biophysical property values are normalized by
    the time of divergence between sequences. Values (V) are thus
    weighted for the transition between amino acids e and k using the time
    (t) since the divergence between sequences i and j:
    
                      Delta_ij = V*(t^-1)_ij
    
    The assumption made in equation 1 is that the different types of amino 
    acid transitions (slight to radical amino acid changes) in a particular 
    site follow a Poisson distribution along time. The greater the time since 
    the divergence between sequences i and j the greater the probability of 
    having a radical change. A linear relationship is thus assumed between the 
    property values and time.
    
    The next step is the estimation of the mean Delta parameter for each site (DC~)
    of the alignment, so that:
    
                   DC~ = (1/L)sigma[1->L]Delta_ij
    
    Here L stands for the total number of pairwise sequence comparisons, and thus:
    
                          L = N*(N-1)/2
    
    Where N is the total number of sequences in the alignment.
    
    The variability of each pairwise amino acid transition compared to that of 
    the site column is estimated as:
    
                    Var_site = (Delta_ij - DC~)^2
    
    The coevolution between amino acid sites (A and B) is estimated thereafter by 
    measuring the correlation in the pairwise amino acid variabilities of two amino
    acid positions (Var_site).
    
    Pair-wise comparisons including gaps in any or both sites at any sequence are 
    excluded from the analysis.
    
    #remove highly diverged sequences

    """


#from __future__ import print_function

import time
import sys
import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment

from scipy.stats import spearmanr
from itertools import combinations




class CoevolError(Exception):
    pass




def Blossum(a1, a2):
    
    """Blosum matrix data return function"""
    
    if a1 == '?' or a1 == 'X':
        a1 = '-'
    if a2 == '?' or a2 == 'X':
        a2 = '-'

    data_Matrix = {
                   'A': [4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2,  0],
                   'C': [0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2,  0],
                   'D': [-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3,  0],
                   'E': [-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2,  0],
                   'F': [-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3,  0],
                   'G': [0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3,  0],
                   'H': [-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2,  0],
                   'I': [-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1,  0],
                   'K': [-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2,  0],
                   'L': [-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1,  0],
                   'M': [-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1,  0],
                   'N': [-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2,  0],
                   'P': [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3,  0],
                   'Q': [-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1,  0],
                   'R': [-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2,  0],
                   'S': [1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2,  0],
                   'T': [0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2,  0],
                   'V': [0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1,  0],
                   'W': [-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2,  0],
                   'Y': [-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7,  0],
                   '-': [0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0]
    }

    tag = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-']

    return data_Matrix[a1][tag.index(a2)]
    
    
def KMAT(a1, a2):

    if a1 == '?':
        a1 = '-'
    if a2 == '?':
        a2 = '-'
        
    data_Matrix = { "A": [ 4, -2, -2, -2,  1, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1,  1,  0, -3, -3,  0, -2, -1,  0, -4, 0],
                    "R": [-2,  5,  0, -2, -4,  1,  0, -2,  1, -3, -2,  2, -1, -2, -2,  0, -1, -2, -2, -2, -1,  0, -1, -4, 0],
                    "N": [-2,  0,  6,  1, -3,  2, -1,  0,  1, -3, -4,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4, 0],
                    "D": [-2, -2,  1,  8, -3,  0,  3, -1,  1, -3, -4, -1, -3, -3, -1,  1, -2, -4, -3, -3,  4,  1, -1, -4, 0],
                    "C": [ 1, -4, -3, -3, 12, -2, -3, -3, -3, -2, -1, -3, -1, -2, -3, -1, -1, -1, -2, -1, -3, -3, -2, -4, 0],
                    "Q": [-1,  1,  2,  0, -2,  4,  1, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2,  0, -2,  0,  3, -1, -4, 0],
                    "E": [ 0,  0, -1,  3, -3,  1,  4, -3,  0, -4, -3,  0, -2, -4, -1, -1, -4, -3, -2, -2,  1,  4, -1, -4, 0],
                    "G": [ 0, -2,  0, -1, -3, -2, -3,  8, -3, -4, -4, -2, -3, -4, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4, 0],
                    "H": [-1,  1,  1,  1, -3,  0,  0, -3,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4, 0],
                    "I": [-1, -3, -3, -3, -2, -3, -4, -4, -3,  4,  3, -4,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4, 0],
                    "L": [-1, -2, -4, -4, -1, -2, -3, -4, -3,  3,  5, -3,  2,  0, -3, -3, -1, -2, -1,  3, -4, -3, -1, -4, 0],
                    "K": [-1,  2,  0, -1, -3,  1,  0, -2, -1, -4, -3,  5, -2, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4, 0],
                    "M": [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -2,  7,  1, -3, -1, -1, -3,  0,  1, -3, -1, -1, -4, 0],
                    "F": [-1, -2, -3, -3, -2, -3, -4, -4, -1,  0,  0, -3,  1,  5, -4, -2, -2,  1,  3,  0, -3, -3, -1, -4, 0],
                    "P": [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -3, -4,  7, -1, -1, -4, -3, -3, -2, -1, -2, -4, 0],
                    "S": [ 1,  0,  1,  1, -1,  0, -1,  0, -1, -2, -3,  0, -1, -2, -1,  4,  2, -3, -1, -2,  0,  0,  0, -4, 0],
                    "T": [ 0, -1,  0, -2, -1, -1, -4, -2, -2, -1, -1, -1, -1, -2, -1,  2,  6, -2, -2,  0, -1, -1,  0, -4, 0],
                    "W": [-3, -2, -4, -4, -1, -2, -3, -2, -2, -3, -2, -3, -3,  1, -4, -3, -2, 11,  3, -3, -4, -3, -2, -4, 0],
                    "Y": [-3, -2, -2, -3, -2,  0, -2, -3,  2, -1, -1, -2,  0,  3, -3, -1, -2,  3,  7, -1, -3, -2, -1, -4, 0],
                    "V": [ 0, -2, -3, -3, -1, -2, -2, -3, -3,  3,  3, -2,  1,  0, -3, -2,  0, -3, -1,  4, -3, -2, -1, -4, 0],
                    "B": [-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4, 0],
                    "Z": [-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4, 0],
                    "X": [ 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4, 0],
                    "*": [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1, 0],
                    "-": [0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0]
                    }
                    
    tag = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*', '-']

    return data_Matrix[a1][tag.index(a2)]


def biophys():
    
    """
        Contains Hydrophobicity and molecular weight data
        
        Column1 => Hydrophobicity data at pH 7
        Column2 => Hydrophobicity data at pH 2
        Column3 => Molecular weight
        """
    
    data_Matrix = {
                    'F': [100, 92, 165.2],
                    'I': [99, 100, 131.2],
                    'W': [97, 84, 204.2],
                    'L': [97, 100, 131.2],
                    'V': [76, 79, 117.2],
                    'M': [74, 74, 149.2],
                    'Y': [63, 49, 181.2],
                    'C': [49, 52, 121.2],
                    'A': [41, 47, 89.1],
                    'T': [13, 13, 119.1],
                    'H': [8, -42, 155.2],
                    'G': [0, 0, 75.1],
                    'S': [-5, -7, 105.1],
                    'Q': [-10, -18, 146.2],
                    'R': [-14, -26, 174.2],
                    'K': [-23, -37, 146.2],
                    'N': [-28, -41, 132.1],
                    'E': [-31, 8, 147.1],
                    'P': [-46, -46, 115.1],
                    'D': [-55, -18, 133.1],
                    '-': [0, 0, 0]

                    }

    return data_Matrix

def kidera():
    """
        Kidera factor from Dan
        """
    data_Matrix = { "A": [-1.56, -1.67, -0.97, -0.27, -0.93, -0.78, -0.2, -0.08, 0.21, -0.48],
                    "R": [0.22, 1.27, 1.37, 1.87, -1.7, 0.46, 0.92, -0.39, 0.23, 0.93],
                    "N": [1.14, -0.07, -0.12, 0.81, 0.18, 0.37, -0.09, 1.23, 1.1, -1.73],
                    "D": [0.58, -0.22, -1.58, 0.81, -0.92, 0.15, -1.52, 0.47, 0.76, 0.7],
                    "C": [0.12, -0.89, 0.45, -1.05, -0.71, 2.41, 1.52, -0.69, 1.13, 1.1],
                    "Q": [-0.47, 0.24, 0.07, 1.1, 1.1, 0.59, 0.84, -0.71, -0.03, -2.33],
                    "E": [-1.45, 0.19, -1.61, 1.17, -1.31, 0.4, 0.04, 0.38, -0.35, -0.12],
                    "G": [1.46, -1.96, -0.23, -0.16, 0.1, -0.11, 1.32, 2.36, -1.66, 0.46],
                    "H": [-0.41, 0.52, -0.28, 0.28, 1.61, 1.01, -1.85, 0.47, 1.13, 1.63],
                    "I": [-0.73, -0.16, 1.79, -0.77, -0.54, 0.03, -0.83, 0.51, 0.66, -1.78],
                    "L": [-1.04, 0, -0.24, -1.1, -0.55, -2.05, 0.96, -0.76, 0.45, 0.93],
                    "K": [-0.34, 0.82, -0.23, 1.7, 1.54, -1.62, 1.15, -0.08, -0.48, 0.6],
                    "M": [-1.4, 0.18, -0.42, -0.73, 2, 1.52, 0.26, 0.11, -1.27, 0.27],
                    "F": [-0.21, 0.98, -0.36, -1.43, 0.22, -0.81, 0.67, 1.1, 1.71, -0.44],
                    "P": [2.06, -0.33, -1.15, -0.75, 0.88, -0.45, 0.3, -2.3, 0.74, -0.28],
                    "S": [0.81, -1.08, 0.16, 0.42, -0.21, -0.43, -1.89, -1.15, -0.97, -0.23],
                    "T": [0.26, -0.7, 1.21, 0.63, -0.1, 0.21, 0.24, -1.15, -0.56, 0.19],
                    "W": [0.3, 2.1, -0.72, -1.57, -1.16, 0.57, -0.48, -0.4, -2.3, -0.6],
                    "Y": [1.38, 1.48, 0.8, -0.56, 0, -0.68, -0.31, 1.03, -0.05, 0.53],
                    "V": [-0.74, -0.71, 2.04, -0.4, 0.5, -0.81, -1.07, 0.06, -0.46, 0.65],
                    "-": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                }
                
    return data_Matrix
    
def calculate_tree_distance(treefile, treeformat):
    """
        Create distance matrix from phylogeny tree file
        """
        
    tree = dendropy.Tree.get_from_path(treefile, treeformat)
    pdm = tree.phylogenetic_distance_matrix()
    optimize_dist = dict()
    pairs = combinations(pdm, 2)
    for (t1, t2) in pairs:
        optimize_dist[str(t1)[1:-1] + "-" + str(t2)[1:-1]] = pdm(t1, t2)
        
    return optimize_dist
        
    
def _thetaEK(records, optimize_dist):
    
    """
    Divergence time based optimized blosum scores for amino acid alignment
    @records - Alignment record object
    @optimize_dist - divergence time
        """
    #posScore = dict()
    posScore = {"B62": dict(), "KMAT": dict(), "K1": dict(), "K2": dict(), "K3": dict(), "K4": dict(), "K5": dict(), "K6": dict(), "K7": dict(), "K8": dict(), "K9": dict(), "K10": dict()}
    msaObj = MultipleSeqAlignment(records)
    KIDERA = kidera()
    for i in range(len(msaObj[1])):
        posVector = msaObj[:, i:i+1]
        print(i)
        if len(set([str(x.seq) for x in posVector])) == 1:
            continue
        pairs = combinations(range(len(posVector)), 2)
        posScore["B62"][i] = list()
        posScore["KMAT"][i] = list()
        for key, val in posScore.items():
            posScore[key][i] = list()
        for (m, n) in pairs: # This for loop has N*(N-1) number of pairs cause N*(N-1) iterations where N is the number of sequences. This needs to be parallelized
            if str(posVector[m].seq) == "?":
                posVector[m].seq = "-"
            if str(posVector[n].seq) == "?":
                posVector[n].seq = "-"
            posScore["B62"][i].append(float(Blossum(str(posVector[m].seq).upper(), str(posVector[n].seq).upper()))/optimize_dist[str(posVector[m].id) + "-" + str(posVector[n].id)])
            posScore["KMAT"][i].append(float(KMAT(str(posVector[m].seq).upper(), str(posVector[n].seq).upper()))/optimize_dist[str(posVector[m].id) + "-" + str(posVector[n].id)])
            for j in range(10):
                posScore["K" + str(j+1)][i].append((KIDERA[str(posVector[m].seq).upper()][j]-KIDERA[str(posVector[n].seq).upper()][j])/optimize_dist[str(posVector[m].id) + "-" + str(posVector[n].id)])


    return posScore


def _meanTheta(posScore):
    meanDist = dict()
    for key, val in posScore.items():
        meanDist[key] = dict()
        for inkey, inval in val.items():
            meanDist[key][inkey] = _average(inval)

    return meanDist


def _variability(posScore):
    varDist = dict()
    for key, val in posScore.items():
        varDist[key] = dict()
        for inkey, inval in val.items():
            mVal = _average(inval)
            varDist[key][inkey] = list()
            for inval_in in inval:
                varDist[key][inkey].append((inval_in - mVal)**2)

    return varDist



############################################ Poisson distance calculation #############################################

def _poisson(distance, long):
    
    pdist = 0
    if 1 - float(distance)/long != 0:
        pdist = -math.log(1 - float(distance)/long)
    else:
        pdist = 0

    return pdist

from sys import getsizeof

def _poisson_dist(pairObj2):
    
    """Poisson distance"""
    
    dist = 0
    gap = 0
    seq1 = records[pairObj2[0] ].seq
    seq2 = records[pairObj2[1] ].seq
    for i, amino in enumerate(seq1):
        if seq1[i] == seq2[i] and seq1[i] != '-' and seq2[i] != '-':
            dist = dist + 1
        if seq1[i] == '-' and seq2[i] == '-':
            gap = gap + 1
    #print dist, len(seq1), gap, seq1, "\n\n", seq2
    return ((str(records[pairObj2[0] ].id) + "-" + str(records[pairObj2[1] ].id)), _poisson(dist, len(seq1) - gap))


def _relative_distance(distance, max_dist):
    for key, val in distance.items():
        if max_dist != 0:
            distance[key] = float(val)/max_dist

    return distance



############################################## Divergence time ###########################################


def find_max_dist(dictObj):
    maxVal = 0
    for key, val in dictObj.items():
        if val > maxVal:
            maxVal = val
        
    return maxVal
    
def _average(s): return sum(s) * 1.0 / len(s)

def _optimize(records):

    """Optimize distance as divergence time"""
    import multiprocessing
    from multiprocessing.dummy import Pool as ThreadPool
    from functools import partial

    n=multiprocessing.cpu_count()
    print(n)
    distance = dict()
    seqNameVec = list()
    testStore = dict()
    start = time.time()
    pairs = list(combinations(range(0,len(records),1), 2))


    
    distance = dict()

    
    
    pool = ThreadPool(n)
    mid = time.time()
    DistDict=pool.imap(_poisson_dist, pairs)
    distance = {key: val for key, val in DistDict}
    

    '''
    for numrec1,rec in enumerate(records): # This for loop and the inner for loop makes N*(N-1) iterations where N is the number of sequences. This needs to be parallelized
        testStore[rec.id] = "NA"
        for numrec1,inrec in enumerate(records):
            try:
                testStore[inrec.id]
                continue
            except KeyError:
                dObj = _poisson_dist([numrec1, numrec1])
                if dObj == 0:
                    dObj = 0.01
                distance[str(rec.id) + "-" + str(inrec.id)] = dObj
    '''

    finish = time.time()
    print(finish - start)
    max_dist = find_max_dist(distance)
    optimize_dist = _relative_distance(distance, max_dist)
    return optimize_dist
    


############################################## Correlation stuff ###########################################

def spearman_test(col1, col2): return spearmanr(col1,col2)

def _correlation(Propkey, varDistVal, storeCorrel):
    toolbar_width = len(varDistVal)*(len(varDistVal)-1)/2
    c = 0
    testStore = dict()
    #corScoreObj = dict()
    for key, val in varDistVal.items(): # This for loop and the inner for loop makes N*(N-1) iterations where N is the number of sequences. This needs to be parallelized
        
        testStore[key] = "NA"
        for inkey, inval in varDistVal.items():
            
            
            try:
                testStore[inkey]
                continue
            except:
                try:
                    storeCorrel[str(key)+"-"+str(inkey)]
                except:
                    storeCorrel[str(key)+"-"+str(inkey)] = dict()
                    
                c = c + 1
    
                try:
                    #print key, inkey, spearman_test(val, inval)
                    storeCorrel[str(key) + "-" + str(inkey)][Propkey] = spearman_test(val, inval)
                except ZeroDivisionError:
                    continue

                p = str((float(c)/toolbar_width)*100)[:5]
                sys.stdout.write("\r%s%%" %p)
                sys.stdout.flush()
            
    return storeCorrel
    
    
############################################################################################################

def execute_coevol(records, rCut=0.9):

    if records == None:
        raise CoevolError("Empty input alignment dataset")
        
    print("Calculating distance\n")
    optimize_dist = _optimize(records)
    posScore = _thetaEK(records, optimize_dist)
    meanDist = _meanTheta(posScore)
    varDist = _variability(posScore)
    storeCorrel = dict()
    for key, val in varDist.items():
        print("Analyzing", key, "\n")
        corScore = _correlation(key, val, storeCorrel)
        
    return corScore

records = list(SeqIO.parse(open("ns1.fas", "rU"), "fasta"))

output = execute_coevol(records)

#output = execute_coevol("ns1.fas")

with open("output.txt", "w") as fp:
    fp.write("Res1\tRes2\t" + "\t".join(output[output.keys()[0]].keys()) + "\n")
    for key, val in output.items():
        pos = key.split("-")
        pos1 = str(int(pos[0])+1)
        pos2 = str(int(pos[1])+1)
        str_write = pos1 + "\t" + pos2 + "\t" + "\t".join([str(round(inval[0], 2)) for inkey, inval in val.items()])
        fp.write("%s\n" %str_write)
        



if(__name__=='__main__'):
    execute_coevol()