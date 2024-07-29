import itertools, re, gzip, os, sys, warnings
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt; 
from pathlib import Path


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return ''.join([complement[base] for base in dna[::-1]])


# Read CMD
fname1 = sys.argv[1]
fname2 = sys.argv[2]
name = sys.argv[3]
print("CMD: ", sys.argv)

# Mkdir
Path("pie1").mkdir(parents=True, exist_ok=True)
Path("pie2").mkdir(parents=True, exist_ok=True)
Path("count").mkdir(parents=True, exist_ok=True)


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

C = {}
for i in range(0, 9):
    C["C{0}".format(i)] = Commons[i]
print(C)

S0 = ["GCTCCT", "GAGGAA"]
S1 = ["GCTCCG", "GAGGAA"]
S2 = ["GCGCCA", "GAAGAG"]
S3 = ["GCACCA", "GAGGAG"]
S4 = ["GCACCC", "GAGGAA"]
S5 = ["GCTCCT", "GAAGAG"]
S6 = ["GCCCCC", "GAGGAA"]
S7 = ["GCTCCT", "GAGGAG"]
Sites = [S0, S1, S2, S3, S4, S5, S6, S7]

S = {}
for i in range(0,8):
    S["S{0}".format(i)] = Sites[i]
print(S)


# Get MutSeq
combinations = list(itertools.product([0,1], repeat=8)) # 0-255
print(combinations[0],
combinations[1],
combinations[255])

dictionary = dict(zip([0,1], ["A", "E"]))
dictionary

mutSeqDict = {}
for i in range(0, 256):
    #print("\n", i)
    
    combination = combinations[i]
    #print(combination)
    
    letter_lst = [dictionary[x] for x in combination]
    letter = "".join(letter_lst)
    #print(">{0}".format(letter))
    
    site_seqs = []
    for j in range(0, 8):
        Skey = "S{0}".format(j)
        seq = S[Skey]
        #print(Skey, seq)
        
        select = combination[j]
        #print(Skey, seq[select])
        site_seqs.append(seq[select])
        
    common_seqs = Commons
    
    mut_seqs = list(zip(common_seqs, site_seqs))
    mut_seqs.append(common_seqs[8])
    mut_seqs
    join1 = [''.join(x) for x in mut_seqs]
    join2 = "".join(join1)
    mutSeq = join2
    # print(mutSeq)
    
    mutSeqDict[letter] = mutSeq
    
mutSeqDict["WWWWWWWW"] = "CTCATGGTTCGGACTTACTTAAAACACCCAAGATGAGGCATTCCGATGGCTTAGAGAAAACCCCATCGCGGTTGATAAGCACACCTAAGGACGGTAACTCGATTTTGAGGAAATGGCAGACTCCTTCACACCTTTTTGAAGATTTGTACTGTTCTCCGCTATTTAGAGCTATAGAGACTCCAATCAGGTATATCACGACGCCGGGGGGCACTTTGGAAACCCAAATTTCACCAAGAAAGTCCTCTGCACCC"

print(len(mutSeqDict), " Mut seqs")


# Save Mut.Fasta
file = open("MutSeqs257.fasta", 'w')
for key in mutSeqDict:
    seq = mutSeqDict[key]
    file.writelines(">{0}\n{1}\n".format(key, seq))    
file.close()


# Read Fastq and Exact Match
def exact_match(fname1, fname2, mutSeqDict=mutSeqDict, test = False):
    count_dict = {}
    n = 0
    with gzip.open(fname1, 'rt') as R1, gzip.open(fname2, 'rt') as R2: # read PE data
        for line1, line2 in zip(R1, R2):
            n += 1
            line1 = line1.strip()
            line2 = line2.strip()

            if n % 4 == 2: # for seq line
                #print("{0}\n{1}\n".format(line1, line2))
                read1 = line1#[5:145] # params no trimming
                read2 = line2#[5:145] # params no trimming
                #read1 = line1[24:86] # test
                #read2 = line2[17:60] # test
                read2_rc = reverse_complement(read2)
                #print("{0}\n{1}\n{2}\n".format(read1, read2, read2_rc))

                match_flag = 0
                mem_key = ""
                for key in mutSeqDict:
                    mutseq = mutSeqDict[key] # not simply mapping, length issue, 251 bp

                    if re.search(read1, mutseq): # R1 search
                        if re.search(read2_rc, mutseq): # R2 search
                            #print(n, "read1 both match", key, read1, mutseq)
                            #print(n, "read2 both match", key, read2, mutseq)
                            #print(key)
                            count_dict[key] = count_dict.get(key, 0) + 1 # .get allows you to specify a default value if the key does not exist.
                            if match_flag: # Bug, if match more than 2 patterns (because of testing shortcuts, should not have this in production code)
                                warnings.warn("More than two matches found for fastq line{0}: {1}, {2}".format(n, mem_key, key))
                            mem_key = key
                            match_flag += 1

            if test:
                if n > 120000:
                    break
    return [count_dict, n//4]


# Exact Match and Count
[count_dict, total] = exact_match(fname1, fname2, test = 0)

# Pie plot discarding 'Other' mutations
labels = list(count_dict.keys())
sizes = list(count_dict.values())

fig1, ax1 = plt.subplots()
ax1.pie(sizes, 
        labels=labels, 
        autopct='%1.1f%%',
        shadow=True, startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.savefig("pie1/" + name + ".all.pdf")

# Pie plot including 'Other' mutations
labels = list(count_dict.keys())
sizes = list(count_dict.values())
size_other = total - sum(sizes)
labels.append("other")
sizes.append(size_other)

fig1, ax1 = plt.subplots()
ax1.pie(sizes, 
        labels=labels, 
        autopct='%1.1f%%',
        shadow=True, startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.savefig("pie2/" + name + ".other.pdf")

# Save

file = open("count/"+name+".txt", 'w')
file.writelines("# ID,count\n")
for key in count_dict:
    count = count_dict[key]
    file.writelines("{0}, {1}\n".format(key, count))
file.writelines("Others,{0}\n".format(size_other))
file.close()
