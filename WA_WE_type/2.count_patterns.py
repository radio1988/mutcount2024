import itertools, re, gzip, os, sys, warnings
from Bio import SeqIO
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt; 
from pathlib import Path

# Dec, 2020: supports 2 additional formats: WE and WA, plus original AE

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return ''.join([complement[base] for base in dna[::-1]])


def exact_match(fname1, fname2, mutSeqDict, test = False):
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


# Read CMD
# python 2.count_patterns.py fname.R1.fastq.gz fname.R2.fastq.gz  MutSeqs511.WA_WE.fasta OUTNAME
fname1 = sys.argv[1]
fname2 = sys.argv[2]
patternFname=sys.argv[3]
name = sys.argv[4]
print("CMD: ", sys.argv)

# Mkdir
Path("pie1").mkdir(parents=True, exist_ok=True)
Path("pie2").mkdir(parents=True, exist_ok=True)
Path("count").mkdir(parents=True, exist_ok=True)



# Read Fastq and Exact Match
patternFname = 'MutSeqs511.WA_WE.fasta'
records = list(SeqIO.parse(patternFname, "fasta"))
mutSeqDict = {}
for r in records:
    mutSeqDict[r.id] = str(r.seq)




# Exact Match and Count
[count_dict, total] = exact_match(fname1, fname2, mutSeqDict=mutSeqDict, test = False)

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
for key, v in sorted(count_dict.items()):
    count = count_dict[key]
    file.writelines("{0}, {1}\n".format(key, count))
file.writelines("Others,{0}\n".format(size_other))
file.close()
