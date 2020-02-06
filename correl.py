# -*- coding: utf-8 -*-
# Copyright 2019 by Ambuj Kumar, Iowa State University.
# All rights reserved.
#
# Bug reports welcome: ambuj@iastate.edu

from scipy.stats import spearmanr

def spearman_test(key1, key2, list1, list2):
    with open("correlations.txt", "a+") as fp:
        out = spearman_test(key1, key2, list1, list2)
        correl = round(out[0], 2)
        fp.write("%s\t%s\%s" %(key1+1, key2+1, correl))
    
        
