# NarrowSense_h2_in_URPs_GP2

### Step 1 - install all the required programs and packages through miniconda3
... (pending)





### Step 2 - get your genotype array data ready 
Before starting we are assuming a couple of things:

1. That your genotype array data has been already called and QCed. The basic parameters expected and more relevant information about how to perform QC in admixed populations is described elsewhere (https://github.com/MataLabCCF/GWASQC).
2. Your phenotype file is free of missing data.
3. Your phenotype file is composed of the individual IDs (IIDs) of your samples, alongside the phenotype (PD status, control = 1, case =2) and basic covariates like Sex (male = 1, female = 2) and Age (quantitative variable).

Here is an example on how the phenotype file should look like:
```
IID	SEX	STATUS	AGE
sample1ID	2	1	38
sample2ID	1	1	33
sample3ID	1	1	45
sample4ID	2	2	66
sample5ID	2	2	72
```

### Step 2.a
Run the following Rscript in order to compute the principal components through PC-AiR (DOI 10.1002/gepi.21896) and PC-Relate (http://dx.doi.org/10.1016/j.ajhg.2015.11.022.) methods (available through the R package GENESIS.

We choose this particular methods since underrepresented populations tend to be admixed, and classical methods to tease apart the related and unrelated datasets could be confused in the presence of genetic admixture. Since both relatedness and admixture are continuum of genetic distance. 
The yield of this approach compared to other more classical methods is not that much different in terms of adjustment for population structure (both can diminish the discovery of spurious association in single variant association testing, see: doi: https://doi.org/10.1101/2025.05.27.25328444, results section). However, since here we are interested in capturing as accurate as possible the h2 estimates, methods that account for population structure in an ancestry-aware approach, are highly encouraged. 

Example of command line:
```
pwd
/working_directory
mkdir pca_and_such
cd pca_and_such/
Rscript pcs_pipeline_h2_project.r \
  --input_plink_file genotyped_data_plink2_prefix \
  --input_format pfile \
  --covar_file covar.txt \
  --removeHighLDregions yes \
  --n_cores 4 2>&1 | tee pcs_pipeline_h2_project.log
```
The " 2>&1 | tee commandLineRun.log  " part is just to save all the print statements into a log file. You can skip it if you want.

The script by default will create an output folder named outFolder_pca_andSuch/ with all the downstream files that we will need.

the input for the script could be plink1 files, so you would run --input_plink_file genotuped_data_plink1_prefix --input_format bfile

you can also modify the default parameters of the script, for more information type: Rscript pcs_pipeline_h2_project.r --help

however, we advise you to keep the default parameters as they are and only change the flags in the example of the command line.
