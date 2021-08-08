# parallel_rTASSEL
Parallelized rTASSEL GLM and MLM using GNU parallel\
rTASSEL is available here: https://github.com/maize-genetics/rTASSEL/

## Required files 
All files were produced with TASSEL 5 GUI version\
Phenotype file was produced as the intersect of all traits with 5 PCs\
as described here: https://avikarn.com/2019-07-22-GWAS/
* HapMap file: genotype.hmp.txt
* Phenotype file: phenotype.txt (includes 5 PCs)
* Kinship file:     kinship.txt (for MLM)
* Empty folders for saving output: \
  results/ \
  plots/   \
  effects/
  
## Running parallel_rTASSEL
<u>num:</u> the number of individual phenotypes in the phenotype file\
parallel_rTASSEL will run each column of the phenotype (from 0 to num)\
on a separate thread and will output each result to the appropriate folder.\
The plots are saved as a PNG files and result stats are saved as TSV files.\
\
To run rTASSEL GLM with GNU parallel:
```
sudo parallel Rscript parallel_rTASSEL_GLM.R {} ::: {0..num}
```
To run rTASSEL MLM with GNU parallel:
```
sudo parallel Rscript parallel_rTASSEL_MLM.R {} ::: {0..num}
```

## Note
Running 36 threads with ~250k SNPs HapMap file requires ~50-60gb RAM.
