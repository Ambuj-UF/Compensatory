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

    """




from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import os

from itertools import combinations



    
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
        with open("prop_vals.txt", "w") as fp:
            fp.write("Matrix\tDiff")
        posVector = msaObj[:, i:i+1]
        if len(set([str(x.seq) for x in posVector])) == 1:
            continue
        
        posScore["B62"][i] = list()
        posScore["KMAT"][i] = list()
        for key, val in posScore.items():
            posScore[key][i] = list()
            
        pairs = combinations(range(len(posVector)), 2)
        for (m, n) in pairs: # This needs paralelization for each pair <<<<<<<<<<<<<<<<<<<<<<<<
            # inputs to property_calc.py - str(m.seq), str(n.seq), optimize_dist <<<<<<<<<<<<<<<<<<<<<<<<
            print "Parellize this loop"

        with open("prop_vals.txt", "a+") as fp:
            for lines in fp:
                if "Matrix" in lines:
                    continue
                lobj = lines.strip().split("\t")
                posScore[lobj[0]][i].append(lobj[1])

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


def _relative_distance(distance, max_dist):
    for key, val in distance.items():
        if max_dist != 0:
            distance[key] = float(val)/max_dist

    return distance


def find_max_dist(dictObj):
    maxVal = 0
    for key, val in dictObj.items():
        if val > maxVal:
            maxVal = val
        
    return maxVal
    
def _average(s): return sum(s) * 1.0 / len(s)

def _optimize(records):

    """Optimize distance as divergence time"""
    with open("distance.txt", "w") as fp:
        fp.write("raw distance data\n")
    pairs = list(combinations(records, 2)) # If there are N sequences in an alignment, it will priduce N*(N-1)/2 pairs
    distance = dict()
    # send pairs as input to poisson.py < input - str(pairs[i][0].id), str(pairs[i][1].id), str(pairs[i][0].seq), str(pairs[i][1].seq) for each ith pair <<<<<<<<<<<<<<<<<<<<<<<<<<<
    with open("distance.txt", "r") as fp:
        for lines in fp:
            if "raw distance" in fp:
                continue
            lobj = lines.strip().split("\t")
            distance[str(lobj[0]) + "-" + str(lobj[1])] = float(lobj[2])

    max_dist = find_max_dist(distance)
    optimize_dist = _relative_distance(distance, max_dist)

    return optimize_dist


def _correlation(Propkey, varDistVal):
    toolbar_width = len(varDistVal)*(len(varDistVal)-1)/2
    c = 0
    storeCorrel = {Propkey: dict()}
    #corScoreObj = dict()
    pairs = combinations(varDistVal, 2)
    with open("correlations.txt", "w") as fp:
        fo.write("pos1\tpos2\tCorrelation Value")
    
    for (key1, key2) in pairs: # <<<<<<<<< parallelize this section for each pair of keys
        #run correly.py on each pair. Input - key1, key2, varDistVal[key1], varDistVal[key2] <<<<<<<<<<<<<
        print "Parrelize this"
        
    with open("correlations.txt", "r") as fp:
        for lines in fp:
            if "pos1" in lines:
                continue
            lobj = lines.strip().split("\t")
            storeCorrel[Propkey][lobj[0]+"-"+lobj[1]] = lobj[1]
            
    return storeCorrel
    
    

        
############################################################################################################

def execute_coevol(filename, rCut=0.9):

    records = list(SeqIO.parse(open(filename, "rU"), "fasta"))
    records = records[:1000]
    if records == None:
        raise CoevolError("Empty input alignment dataset")
    print "Calculating distance\n"
    optimize_dist = _optimize(records)
    posScore = _thetaEK(records, optimize_dist)
    meanDist = _meanTheta(posScore)
    varDist = _variability(posScore)
    storeCorrel = dict()
    for key, val in varDist.items():
        print "Analyzing", key, "correlations\n"
        storeCorrel = _correlation(key, val, storeCorrel)
        
    return storeCorrel

output = execute_coevol("2K1.fas")

with open("output.txt", "w") as fp:
    fp.write("Res1\tRes2\t" + "\t".join(output[output.keys()[0]].keys()) + "\n")
    for key, val in output.items():
        pos = key.split("-")
        pos1 = str(int(pos[0])+1)
        pos2 = str(int(pos[1])+1)
        str_write = pos1 + "\t" + pos2 + "\t" + "\t".join([str(round(inval[0], 2)) for inkey, inval in val.items()])
        fp.write("%s\n" %str_write)
        

#print "\n\n"
#print "Printing B62, KMAT and kidera factors correlations scores as dictionary below:\n\n"
##for key, val in x.items():
#print key, val, "\n"

