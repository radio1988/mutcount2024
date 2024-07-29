# mutcount2024

This repo host the code to quantify exact match count for Phosphosite Scanning sequencing data developed by Jennifer Benanti lab.

- The code in `AE_type/` folder is older, but it was the exact code used to reproduce the results for AE type

- The code in `WA_WE_type/` folder is newer, and it can be used to reproduce the results for WA_WE type

Under each folder, a readme file can be found that shows how to run the code on the example data


## Details
- AE: library containing all A and E substitutions at each phosphosite, in all possible combinations, in addition WWWWWWWW was added, creating 257 patterns

- WA: library containing the wild type sequence (W) or A substitutions at each phosphosite, in all possible combinations from WWWWWWWW to AAAAAAAA, in addition EEEEEEEE were added, making 257 combinations

- WE: library containing the wild type sequence (W) or E substitutions at each phosphosite, in all possible combinations from WWWWWWWW to EEEEEEEE, in addition AAAAAAAA were added, making 257 combinations

- WA_WE: All 257 combinations from WA and 257 combinations from WE were merged, then one copy of the 3 duplicated combinations (WWWWWWWW, AAAAAAAA, and EEEEEEEE) were removed, making 511 unique combinations

- Only exact matches in read fragments are counted

- All read fragments without an exact match are counted as 'others' 