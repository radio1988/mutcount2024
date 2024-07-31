# mutcount2024

This repo hosts the code to quantify exact match counts for Phosphosite Scanning sequencing data of HCM1, developed in Jennifer Benanti's lab.

- The code in `AE_type/` can be used to reproduce the results from screening of the HCM1 AE library

- The code in `WA_WE_type/` folder can be used to reproduce the results for from screening the combined HCM1 WA and WE libraries

Each folder contains a readme file that details how to run the code on example data


## Details
- AE: library containing all A and E substitutions at each phosphosite, in all possible combinations, in addition to the completely wild type sequence, resulting in 257 unique sequences
-	WA: library containing the wild type sequence (W) or A substitutions at each phosphosite, in all possible combinations
-	WE: library containing the wild type sequence (W) or E substitutions at each phosphosite, in all possible combinations
- WA_WE: library that includes all sequences in the WA and WE libraries, resulting in 511 unique sequences  
-	Only exact matches in read fragments are counted
-	All read fragments without an exact match are counted as 'others'
