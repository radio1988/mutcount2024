import itertools, re, gzip, os, sys, warnings
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt; 
from pathlib import Path

# Read CMD
# python 1.create_patterns.py WA_WE
FMT=sys.argv[1] # WA
print("CMD: ", sys.argv)

# Params
## C: Common Seq
## S: Site Seq
C0 = "CTCATGGTTCGGACTTACTTAAA"
C1 = "AAGATGAGGCATTCCGATGGCTTAGAGAAA"
C2 = "TCGCGGTTGATAAGC"
C3 = "AAGGACGGTAACTCGATTTTGAGGAAATGGCAG"
C4 = "TCACACCTTTTTGAAGATTTGTACTGT"
C5 = "CTATTTAGAGCTATAGAG"
C6 = "ATCAGGTATATCACG"
C7 = "GGGGGCACTTTGGAAACCCAAATT"
C8 = "AGAAAGTCCTCTGCACCC"
Commons = [C0, C1, C2, C3, C4, C5, C6, C7, C8]

name2Seq_Commons = {}
for i in range(0, 9):
    name2Seq_Commons["C{0}".format(i)] = Commons[i]
print(name2Seq_Commons)
# name2Seq_Commons['C0']: "CTCATGGTTCGGACTTACTTAAA" 

mutSeqDict = {}
# mutSeqDict["WWWWWWWW"] = "CTCATGGTTCGGACTTACTTAAAACACCCAAGATGAGGCATTCCGATGGCTTAGAGAAAACCCCATCGCGGTTGATAAGCACACCTAAGGACGGTAACTCGATTTTGAGGAAATGGCAGACTCCTTCACACCTTTTTGAAGATTTGTACTGTTCTCCGCTATTTAGAGCTATAGAGACTCCAATCAGGTATATCACGACGCCGGGGGGCACTTTGGAAACCCAAATTTCACCAAGAAAGTCCTCTGCACCC"

if FMT == "AE":
    # MC1, MC2 (Nov, 2019)
    S0 = ["GCTCCT", "GAGGAA"]
    S1 = ["GCTCCG", "GAGGAA"]
    S2 = ["GCGCCA", "GAAGAG"]
    S3 = ["GCACCA", "GAGGAG"]
    S4 = ["GCACCC", "GAGGAA"]
    S5 = ["GCTCCT", "GAAGAG"]
    S6 = ["GCCCCC", "GAGGAA"]
    S7 = ["GCTCCT", "GAGGAG"]
    num2Name = dict(zip([0,1], ["A", "E"]))
    mutSeqDict["WWWWWWWW"] = "CTCATGGTTCGGACTTACTTAAAACACCCAAGATGAGGCATTCCGATGGCTTAGAGAAAACCCCATCGCGGTTGATAAGCACACCTAAGGACGGTAACTCGATTTTGAGGAAATGGCAGACTCCTTCACACCTTTTTGAAGATTTGTACTGTTCTCCGCTATTTAGAGCTATAGAGACTCCAATCAGGTATATCACGACGCCGGGGGGCACTTTGGAAACCCAAATTTCACCAAGAAAGTCCTCTGCACCC"
    outname="MutSeqs257.AE.fasta"
elif FMT == "WA":
    # WT_A (Nov, 2020)
    S0 = ["ACACCC", "GCTCCT"]
    S1 = ["ACCCCA", "GCTCCG"]
    S2 = ["ACACCT", "GCGCCA"]
    S3 = ["ACTCCT", "GCACCA"]
    S4 = ["TCTCCG", "GCACCC"]
    S5 = ["ACTCCA", "GCTCCT"]
    S6 = ["ACGCCG", "GCCCCC"]
    S7 = ["TCACCA", "GCTCCT"]
    num2Name = dict(zip([0,1], ["W", "A"])) # {0: 'W', 1: 'A'}
    mutSeqDict["EEEEEEEE"] = "CTCATGGTTCGGACTTACTTAAAGAGGAAAAGATGAGGCATTCCGATGGCTTAGAGAAAGAGGAATCGCGGTTGATAAGCGAAGAGAAGGACGGTAACTCGATTTTGAGGAAATGGCAGGAGGAGTCACACCTTTTTGAAGATTTGTACTGTGAGGAACTATTTAGAGCTATAGAGGAAGAGATCAGGTATATCACGGAGGAAGGGGGCACTTTGGAAACCCAAATTGAGGAGAGAAAGTCCTCTGCACCC"
    outname="MutSeqs257.WA.fasta"
elif FMT == "WE":
    # WT_E (Nov, 2020)
    S0=['ACACCC', 'GAGGAA']
    S1=['ACCCCA', 'GAGGAA']
    S2=['ACACCT', 'GAAGAG']
    S3=['ACTCCT', 'GAGGAG']
    S4=['TCTCCG', 'GAGGAA']
    S5=['ACTCCA', 'GAAGAG']
    S6=['ACGCCG', 'GAGGAA']
    S7=['TCACCA', 'GAGGAG']
    num2Name = dict(zip([0,1], ["W", "E"]))
    mutSeqDict["AAAAAAAA"] = "CTCATGGTTCGGACTTACTTAAAGCTCCTAAGATGAGGCATTCCGATGGCTTAGAGAAAGCTCCGTCGCGGTTGATAAGCGCGCCAAAGGACGGTAACTCGATTTTGAGGAAATGGCAGGCACCATCACACCTTTTTGAAGATTTGTACTGTGCACCCCTATTTAGAGCTATAGAGGCTCCTATCAGGTATATCACGGCCCCCGGGGGCACTTTGGAAACCCAAATTGCTCCTAGAAAGTCCTCTGCACCC"
    outname="MutSeqs257.WE.fasta"

Sites = [S0, S1, S2, S3, S4, S5, S6, S7]
S = {}
# S[S0]: ["ACACCC", "GCTCCT"]
for i in range(0,8):
    S["S{0}".format(i)] = Sites[i]
print(S)
# {'S0': ['ACACCC', 'GCTCCT'],...}


# Get MutSeq
combinations = list(itertools.product([0,1], repeat=8)) # 0-255
# [(0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 1),...]
print(combinations[0], combinations[1], combinations[255])

for i in range(0, 256):
    #print("\n", i)
    
    combination = combinations[i] #  (0, 0, 0, 0, 0, 0, 0, 1)
    #print(combination)
    
    letter_lst = [num2Name[x] for x in combination] # W,W,...,A
    letter = "".join(letter_lst) # WWWWWWWA
    #print(">{0}".format(letter))
    
    site_seqs = []
    # ['ACACCC','ACCCCA',..., 'GCTCCT']
    for j in range(0, 8):
        Skey = "S{0}".format(j) # 'S0'
        seq = S[Skey] # ["ACACCC", "GCTCCT"]
        #print(Skey, seq)
        
        select = combination[j] # 0 or 1
        #print(Skey, seq[select])
        site_seqs.append(seq[select]) # one of the seqs
        
    common_seqs = Commons
    
    mut_seqs = list(zip(common_seqs, site_seqs)) # first 8
    mut_seqs.append(common_seqs[8]) # 9th common seq
    mut_seqs
    # [('CTCATGGTTCGGACTTACTTAAA', 'ACACCC'),..., 'AGAAAGTCCTCTGCACCC']
    join1 = [''.join(x) for x in mut_seqs]
    mutSeq = "".join(join1)
    # print(mutSeq)
    # 'CTCATGGTTCGGACTTACTTAAAACACCCAAGATGAGGCATTCCGATGGCTTAGAGAAAACCCCATCGCGGTTGATAAGCACACCTAAGGACGGTAACTCGATTTTGAGGAAATGGCAGACTCCTTCACACCTTTTTGAAGATTTGTACTGTTCTCCGCTATTTAGAGCTATAGAGACTCCAATCAGGTATATCACGACGCCGGGGGGCACTTTGGAAACCCAAATTGCTCCTAGAAAGTCCTCTGCACCC'
    
    mutSeqDict[letter] = mutSeq
    

print(len(mutSeqDict), " Mut seqs")


# Save Mut.Fasta
file = open(outname, 'w')
for key in mutSeqDict:
    seq = mutSeqDict[key]
    file.writelines(">{0}\n{1}\n".format(key, seq))    
file.close()

