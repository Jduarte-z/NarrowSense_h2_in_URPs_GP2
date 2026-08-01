# NarrowSense_h2_in_URPs_GP2

### Step 1 - get your genotype data ready 

Before starting we are assuming a couple of things:

1. That your genotype array data has been already called and QCed. The basic parameters expected and more relevant information about how to perform QC in admixed populations is described elsewhere (https://github.com/MataLabCCF/GWASQC).
2. Your phenotype file is free of missing data.
3. Your phenotype file is composed of the individual IDs (IIDs) of your samples, alongside the phenotype (PD status, control = 1, case =2) and basic covariates like Sex (male = 1, female = 2) and Age (quantitative variable).

Here is an example on how the phenotype file should look like:
```IID	SEX	STATUS	AGE
sample1ID	2	1	38
sample2ID	1	1	33
sample3ID	1	1	45
sample4ID	2	2	66
sample5ID	2	2	72
```
