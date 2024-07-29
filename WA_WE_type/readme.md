#  Steps to run the test example

1. `cd example_run`

2. `python ../1.create_patterns.py WE &> MutSeqs257.WE.fasta.log`

3. `python ../1.create_patterns.py WA &> MutSeqs257.WA.fasta.log`

4. manually merge MutSeqs257.WE.fasta and MutSeqs257.WA.fasta into MutSeqs511.WA_WE.fasta
	-  `cat MutSeqs257.WE.fasta  MutSeqs257.WA.fasta  > MutSeqs514.WA_WE.fasta`
	- remove the duplicated AAAAAAAA, EEEEEEEE, WWWWWWWW from `MutSeqs514.WA_WE.fasta` and create `MutSeqs511.WA_WE.fasta`

5. `python ../2.count_patterns.py fastq/test.R1.fq.gz fastq/test.R2.fq.gz MutSeqs511.WA_WE.fasta test &> count.log`


