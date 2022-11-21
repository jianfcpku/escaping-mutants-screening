#!/usr/bin/python

import pandas as pd
import numpy as np
import os
import sys
import re

var_f = sys.argv[1]
var_o = sys.argv[2]

Alleles_d =pd.read_table(var_f,header=None,index_col=None,sep ='\t')
Alleles_d.index = Alleles_d[1]
ACCG = ['A','T','C','G']

def Allele_line_process (A_line):
    AAs = [A_line[3]]
    AAs.extend(re.split(',',A_line[4]))
    rr = re.split(',',A_line[8])
    aa_dict = {'A':0, 'T':0,'C':0,'G':0}
    for mm in range(len(AAs)):
        aa_dict[AAs[mm]] = rr[mm]
    Final_line = [A_line[0],A_line[1], A_line[3], A_line[7]]
    for mm in ACCG:
        Final_line.append(aa_dict[mm])
    
    return Final_line

#Alleles_d_cut = Alleles_d.loc[(Alleles_d[4]!='<*>')&(Alleles_d[7]!='INDEL'),]
Alleles_d_cut = Alleles_d.loc[(Alleles_d[4]!='<*>'),]

Alleles_d_cut_d = pd.DataFrame(Alleles_d_cut.apply(Allele_line_process,axis=1).tolist())
Alleles_d_cut_d.columns = ["Chrom", "Pos", "Ref", "Depth", "A", "T", "C", "G"]

Alleles_d_cut_d.to_csv(var_o, header=True, index=None,sep='\t')




