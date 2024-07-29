#  Steps to run test example

1. `python 1.create_patterns.py WE &> MutSeqs257.WE.fasta.log`

2. `python 1.create_patterns.py WA &> MutSeqs257.WA.fasta.log`

3. manually merge MutSeqs257.WE.fasta and MutSeqs257.WA.fasta into MutSeqs511.WA_WE.fasta

4. `mkdir test_output; mv MutSeqs* test_output; cd test_output`

5. `python ../2.count_patterns.py ../fastq/test.R1.fq.gz ../fastq/test.R2.fq.gz MutSeqs511.WA_WE.fasta test &> count.log`



