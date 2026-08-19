# Tutorial for the estimation of Narrow-Sense SNP Heritability in Underrepresented Populations within the Global Parkinson's Genetics Program

*A hands-on tutorial: concepts, objectives, and methods.*

---

## 1. Background

### 1.1 The quantitative genetic model

Under the classical quantitative genetics framework, an individual's phenotype for a complex trait can be decomposed as:

$$
P = G + E = (A + D + I) + E
$$

where $$P$$ is the observed phenotype, $$G$$ genetics, and $$E$$ environment. Within $$"G"$$ $$A$$ is the additive genetic effect (the sum of average effects of alleles), $$D$$ is the dominance effect (within-locus allelic interaction), $$I$$ is the epistatic effect (between-locus interaction), and $$E$$ is the environmental deviation.

Because these components vary across individuals, the corresponding decomposition of the total phenotypic variance in a given population is:

$$
VP = VA + VD + VI + VE + 2·Cov(G,E) + VGxE
$$

The last two terms (the genotype–environment covariance and the genotype-by-environment interaction variance) are usually left out of the heritability analysis that we are about to.

Heritability estimates should not be interpreted as environment-independent estimates, but rather as population and its particular
distribution of environment specific.  

In this GitHub, we do not provide a comprehensive coverage of these topics. For useful references you could refer to books like: Falconer, D.S. & Mackay, T.F.C. (1996). Introduction to Quantitative Genetics, and Lynch, M. & Walsh, B. (1998). Genetics and Analysis of Quantitative Traits.


### 1.2 Broad and narrow-sense heritability

Heritability is a **ratio**: it is the proportion of phenotypic variance attributable to genetic variance.

- **Broad-sense heritability:** $$H² = (VA + VD + VI) / VP$$
- **Narrow-sense heritability:** $$h² = VA / VP$$

Two properties are worth to mention:

1. Heritability is a property of a **population in an environment**, not of a trait and certainly not of an individual. The same trait can have different heritabilities in multiple cohorts that have multiple environmental and/or genetic components. 

2. Heritability *within* groups carries no information about the causes of mean differences *between* groups, see [J.G. Schraiber & M.D. Edge (2024)](https://www.pnas.org/doi/10.1073/pnas.2319496121)


### 1.3 Twin studies, GWAS, and the "missing heritability" problem

Since the mid-20th century, heritability of complex traits was estimated primarily from twin and family designs. The classical ACE model partitions phenotypic variance into additive genetic ($$A$$), shared environmental ($$C$$), and unique environmental ($$E$$) components by exploiting the difference in genetic sharing between monozygotic (~1.0) and dizygotic (~0.5) twin pairs. Falconer's estimator, $$h² = 2(r_MZ − r_DZ)$$, is nominally a narrow-sense estimator, though in practice dominance, epistasis, shared environment, and assortative mating all load onto the $$A$$ term.

With the advent of genome-wide association studies, heritability could be estimated from genotype data in samples of "unrelated individuals". These methods estimate a distinct quantity:

**SNP heritability ($$h²_SNP$$, sometimes $$h²_g$$)**: the proportion of phenotypic variance explained by the additive effects of the genotyped (or imputed) variants included in the model (usually, Single Nucleotide Polymorphisms).

$$h²_SNP$$ is expected to be smaller than $$h²$$ because:

- Genotyping arrays interrogate common variation and tag rare causal variants poorly.
- Even where causal variants are present, the additive random-effect model assumes independent, exchangeable SNP effects, which is an approximation.
- Non-additive variance is excluded by construction (and dominant, epistatic effects, etc. are left out)

The gap between twin-based and GWAS-based estimates is the **"missing heritability" problem**. It reflects a combination of genuinely untagged variance, non-additive contributions, and the nature of family-based estimates (e.g. shared environment and assortative mating). A useful entry point to this literature is Manolio et al. (2009), [*Finding the missing heritability of complex diseases*](https://pmc.ncbi.nlm.nih.gov/articles/PMC2942068/), since we don't cover this aspect extensively in this repository.  

**Hence, with this tutorial, you will be able to estimate $$h²_SNP$$, not $$h²$$**

### 1.4 Binary traits and the liability scale

All of the previous discussions were in principle developed for quantitative traits, but Parkinson's disease is a dichotomous outcome, which introduces two complications.

**Scale.** Under the liability threshold model, each individual has an unobserved, normally distributed liability; individuals whose liability exceeds a threshold determined by the population prevalence $$K$$ are affected. Heritability estimated on the observed 0/1 scale depends on the case-control ratio in the sample and is not comparable across studies. Estimates must therefore be reported on the **liability scale**, which is a property of the population rather than of the study design.

**Ascertainment.** Case-control cohorts are, by design, enriched for cases relative to the population ($$P >> K$$). This oversampling biases the naive observed-to-liability transformation and, more seriously, biases REML-based estimators when covariates such as ancestry principal components are included in the model.

For PD we use $$K = 0.005$$ as the primary population prevalence, consistent with the values used in previous published literature, and report sensitivity across $$K ∈ {0.005, 0.01, 0.02}$$. Prevalence is not measured with precision in underrepresented populations, and this uncertainty propagates directly into the liability-scale estimate. 

### 1.5 Why estimate $$h²_SNP$$, and why in underrepresented populations?

The estimation of $$h2_SNP$$ is relevant by its own means regardless of the population. Consider that it is computed using the same resources allocated for conducting genome-wide association studies, and since it represents the phenotypic variance explained by the same genetic variation screened in a GWAS framework, its value helps us to get closer to a comprehensive mapping of the genetic architecture of a complex trait in a given population. For example, giving a potential upper limit estimate on the discovery power that GWAS could have, or the virtual upper limit of the performance of polygenic risk scores derived from the GWASs in question (since the accuracy of a PRS would approach heritability estimates as sample size and increase genetic variation is captured without bound). 

And beyond this point, exploring $$h2_SNP$$ estimation in underrepresented populations becomes meaningful in several ways. First, it will increase the methodological capacity so it does not depend on European reference data and pipelines. Nearly every heritability estimator in common use was developed and validated on European-ancestry samples, and their behavior in underrepresented, admixed cohorts is comparatively uncharacterized, especially for Parkinson's Disease. Additionally, it will help us to showcase the heterogeneity, challenges, limitations and future directions when computing and interpreting these type of estimates from diverse genomic data. 

---

## 2. Objectives

### Repository objective

Provide a reproducible, step-by-step workflow for estimating SNP heritability in underrepresented and admixed populations from raw genotype, phenotype, and GWAS summary data.

### Project aims

1. **Lay out the analytical framework** for estimating $$h²_SNP$$ in underrepresented populations (URPs), with explicit management of admixture, relatedness, PD case ascertainment bias, liability scale transformation and different heritability models. 
2. **Collaborate internationally** to build the capacity needed to estimate $$h2_SNP$$ as accurately as possible in an increasingly diverse genomic landscape. Powered by and for URPs. 
3. **Characterize heterogeneity, methodological uncertainty and challenges** when computing this estimates, in light of a field whose methods were developed and validated mainly on European-ancestry data.

---

## 3. Methods overview

A comprehensive review of heritability estimation methods is beyond the scope of this repository. A curated reference list is available [here](https://docs.google.com/document/d/1mcoMHNsUat0rzlDItxcBWFDSfRwRMm4rFGBU8VqcIG4/edit?usp=sharing). Disclaimer, the list is not intended to be comprehensive and its reading is not mandatory, but it definitely will help to satisfy your curiosity to a certain extent. 

On that note, we can make a broad classification of the current methods for estimating $$h2_SNP$$ based on the concept similarity that they exploit. Plus some caveats on their challenges when applied to URPs. 

1. **Methods that exploit linkage disequilibrium among genetic variants (LD scores)**: they require the presence of an LD reference panel appropriate to the population in question and tend to give more conservative (lower) heritability estimates. Classically exemplified in the Linkage Disequilibrium Score Regression (LDSC) software. Some notes about it:
   - **(a)** For European populations this is the go-to method when sample sizes are big, since the core logic is based on regressing the per-variant chi-square of GWAS summary statistics on the LD scores from the reference panel. So you will only need summary statistics and European reference panel-derived LD scores.
   - **(b)** When the mean chi-square of a given GWAS summary statistics is greater than 1, that means that mean association metrics for each variant across the genome are inflated. And this could be due to population stratification and/or the polygenicity of a given trait. Precisely regressing the per-variant chi-square from a GWAS on population-matched LD scores could help to differentiate these two, while at the same time computing the heritability of the complex trait: confounding raises the regression **intercept** (it inflates chi-squares regardless of how much LD a variant tags), whereas true polygenic signal raises the **slope** (variants tagging more of the genome carry more signal), and it is the slope that yields the heritability estimate. Just FYI, The authors of LDSC recommend a mean chi-square of **at least** 1.02 in order for a GWAS to be suitable for LDSC regression (in the pipeline we compute it). 
   - **(c)** For URPs, appropriate reference panels are usually lacking, and the admixture that characterizes these cohorts violates some of the mathematical assumptions of the method. Hence, individual data is required to compute an in-sample reference panel and adjust the LD scores for the long range LD that is present in admixed populations (using the same principal components that are included in the GWAS that generates the summary statistics to be used, see below). Two adjustments are involved: correcting the pairwise $$r2$$ for those principal components, and widening the window over which LD scores are computed (20 cM instead of the 1 cM default), because admixture LD decays far more slowly than LD within a single ancestral population. This is the main logic behind the adaptation of LDSC for admixed populations called covariate-LDSC (cov-LDSC). And in case you were wondering, when cov-LDSC is benchmarked in more "homogeneous populations" they yield very similar estimates to the classic LDSC (for example using an in-sample European ref panel instead of One Thousand Genomes Project that is the standard).
2. **Methods that exploit genotypic similarity among unrelated individuals**: they require the computation of Genomic Relatedness Matrices (GRMs), tend to give higher heritability estimates, and their computational cost grows steeply with sample size (the GRM is $$N x N$$, and GREML additionally inverts an $$N x N$$ matrix at every iteration). Some notes about them:
   - **(a)** The most famous example for these type of methods is embodied by the Genome-Wide Complex Trait Analysis Genomic relatedness matrix restricted maximum likelihood or GCTA-GREML toolkit. It fits a linear mixed model in which the genetic contribution of each individual enters as a random effect whose covariance is proportional to the GRM, and estimates the genetic and residual variance components by restricted maximum likelihood (REML).
   - **(b)** Alternative methods regress the phenotypic cross-product of each pair of individuals on their genotypic similarity, known as Haseman-Elston (HE) regression. The Phenotype Correlation-Genotype Correlation (PCGC) method extends HE regression by adding a correction for case ascertainment under a liability-threshold model, which is what makes it appropriate for ascertained binary traits. Because these are moment-based estimators, they do not require specifying a full likelihood for the observed data (PCGC still assumes the liability-threshold model), they are more robust under case ascertainment bias, and they are a bit more computationally efficient. However, they could be a bit more unstable when sample sizes are on the lower side (to keep in mind).
   - **Both (a) and (b)** estimate heritability from a GRM, a table of how genetically similar each pair of participants is. On the logic that if a trait is heritable, people who share more genome should also resemble each other more in that trait. Note that this signal comes from chance variation in genome-wide sharing among individuals who are **not** close relatives, which is why relatives are pruned out beforehand. What the matrix has to capture, then, is genetic similarity measured against an *ancestry-matched* expectation, but the standard calculation compares everyone to a single cohort-average allele frequency, which in URPs that tend to be admixed describes no one accurately. For example, two unrelated participants who both carry high proportions of a given ancestry will deviate from the average in the same direction at thousands of variants, and shared ancestry could be mistaken for excess genome sharing, breaking some of the assumptions made by the methods in question. Hence, for URPs, an ancestry-aware GRM is highly advisable to compute and benchmark alongside classical methods (see below).

Furthermore, there are a couple of different additional "axes" to classify $$h2_SNP$$ estimation methods, in particular:

1. **Based on the Heritability model**: every estimator carries a prior assumption about how heritability is distributed across the genome. Specifically, how variant's expected contribution depends on its allele frequency and how much LD it sits in. The general form for describing this is $$E[h²_j] ∝ w_j · [f_j(1−f_j)]^(1+α)$$, where $$f_j$$ is allele frequency and $$w_j$$ is an LD-based weight.
	- **(a)** All of the above methods tend to do $$w_j = 1$$ and, conventionally, $$α = −1$$ (the so called GCTA model). This makes $$E[h²_j]$$ constant: every SNP is expected to contribute equally, regardless of frequency or LD. LDSC assumes essentially the same thing.
	- **(b)** The human default model described by the creators of [LDAK](https://dougspeed.com/heritability-model/), is often described as an alternative to the classic GCTA model. It sets $$α = −0.25$$ and it is currently recommended for the GRM based methods. However, the GRM must be computed by the LDAK software itself, which limits their adaptation for an ancestry-aware GRM (a currently limitation of the field). 
2. **Based on genome-wide partition**: this alludes whether a single variance component for all genetic variants is used, or a split into bins based on different parameters is implemented. 
	- **(a)** All of the above methods (in the way we intend them to be used) use a single variance component, aka, they use genome-wide data all at once. However an alternative is to decompose genome-wide data in bins based on LD scores or minor allele frequency, since it has been shown that heritability tend to vary based on different thresholds for this data. However, the cost of including more parameters could be larger standard errors and the need for bigger sample sizes that underrepresented cohorts frequently don't have. However, they have not been benchmarked in URPs, and based on your input discussion, we can brainstorm if it worth to include in our analysis (since the original project is based on single variance component methods). 

---

## 4. Putting the benchmarking pipeline together 

Based on the brief overview of the main methods available, you can quickly realize that even for European populations there is no gold-standard method for estimating SNP heritability. And for URPs there are additional problems to solve. Hence, considering that each estimator differ in what they assume and their specific requirements and limitations in diverse populations, we aim to get a spread of estimates across a pre-specified set of methods and across a somewhat wide range of different plausible PD prevalences. 

Our analysis will follow two prongs, a) using methods that exploit  LD scores adapted for URPs (steps 1-7), and b) methods that exploit genotypic similarity among unrelated samples (steps 8-12). 

For a) we are going to focus on cov-LDSC that requires an in-sample reference panel, computation of adjusted LD scores by principal components and the summary statistics from a GWAS ran on the same samples. The caveat of this method in URPs is that the principal components included in the adjustment of the LD scores in the in-sample reference panel, must be the same that were included in the GWAS regression model. And it could not be from a generalized linear mixed model (GLM) (because the random effects of GLMs break the math assumed by LDSC) Hence, in order to enforce this requirement, we will run a GWAS with plink (a generalized linear model) with the same universe of samples and SNPs to be included in the reference panel. 

For b) we are going to use GCTA-GREML and PCGC, with both a regular GRM and an ancestry-aware GRM computed through the R package PC-Relate. Hence, a total of four different GRM based heritability results will be reported. Additionally, because PCGC is computationally cheap enough to be re-run many times, we will derive estimates under the null hypothesis of no genetic signal by shuffling the sample labels of the classic GRM kinships and of the ancestry-aware GRM ones.

Here is a summary of the analysis to be performed:

| # | Arm | GRM | Model | Role in the design |
|:---|:---|:---|:---|:---|
| 1 | **cov-LDSC**, in-sample LD scores | — | α =−1 | Summary-statistic based, with the caveat of in-sample reference panel, validated in admixed cohorts by Luo et al., and completely independent of GRM construction. |
| 2 | **GREML** + standard GCTA GRM | GCTA | α =−1 | One of the field's reference implementations. This is the number most directly comparable to older published PD heritability estimates, like Keller et al., |
| 3 | **GREML** + PC-Relate GRM | PC-Relate | α =−1 | Isolates GRM construction. Same estimator, same SNPs, same covariates as arm 2 — only ancestry residualisation differs. |
| 4 | **PCGC** + standard GCTA GRM | GCTA | α =−1 | Isolates ascertainment handling. Same GRM as arm 2 — only the liability-scale estimator differs. |
| 5 | **PCGC** + PC-Relate GRM | PC-Relate | α =−1 | Both corrections applied together. |
| 6 | **PCGC permutation** | PC-Relate kinships, shuffled | α =−1 | Null control. |
| 7 | **PCGC permutation** | GCTA kinships, shuffled | α =−1 | Null control. |

Here, now you may be noticing two things: 1) that all of these methods rely on a uniform heritability model (alpha = -1). And 2) none of them incorporate genome-wide partition adaptations. 

And the reason of not including these two other "axes" of heritability estimation methods rely on the conditions that they impose to the analysis. For example, if we wanted to incorporate a different heritability model, like the LDAK Human Default model (alpha = -0.25) we would have to restrict the benchmark to the GRM computed by this software, since the SNP weights and the thinning are encoded in its computation itself, so no ancestry-aware side by side benchmark for this model could be done (one of the limitations of current methods). On the other hand, with respect to genome-wide partitioning, they are not part of the main analysis since it has been shown to require higher sample sizes. 

However, they have not been benchmarked yet in URPs, and we are also open to discuss if include these type as part of sensitivity analysis in the future. Hence, your input as collaborator is crucial, and if you consider that these are worth exploring we can think about ways on how to do it. Just take into consideration that the initial analytical plan approved by the Project Proposal, Approval, and Execution Working Group within GP2 covers the seven arms listed above, so any addition would be reported as sensitivity analysis rather than as part of the primary approach. 
 
---

# Starting the analysis 

Without further ado, lets begin with the steps needed to be undertaken in the dry lab!

We will start with cov-LDSC and then move in the same order as described in the table of section 4. with the overview of the pipeline. 

An assumption that I'm making here is that you are working in a Linux machine and that you have already installed miniconda3. If you haven't yet, please take a look at this: https://www.anaconda.com/docs/getting-started/miniconda/install/linux-install.

I don't anticipate it will make that much of a difference working with slightly newer or older versions of conda, but just FYI, I was working with conda 26.5.3 version.

An **important disclaimer**: this tutorial was thought so you can run the analysis through the terminal of your institution computing cluster, and all the names of folders and files used here are generic, so relative paths could be use and provide a (hopefully) easy run. Down below I show some examples when working with SLURM and sending multiple jobs to your computing cluster administrator system so the analysis could speed up (with several jobs in parallel when feasible). However, if you don't have access to a computing cluster, and are working in cloud services like verily, I also provide script here in order to run it interactively with the terminal (just consider that this part will take longer).

I highly encourage you to follow the same naming and folder structure patterns described here so we can maximize the reproducibility of the pipeline. Additionally, here we use a bunch of different programs and environments that require specific versions. How to download them and use them will be covered as soon as we need them. Just consider that in order to avoid conflict between possible versions of the programs that you may have, stick with the naming convention and paths that we use here, in that way, everything could be contained in a specific folder dedicated to this analysis, and (hopefully) we don't mess up with your current program's set up. 

Another 

# First wave of analysis: Covariate-adjusted Linkage Disequilibrium score regression (cov-LDSC)

Standard LDSC estimates $$h²_SNP$$ by regressing GWAS chi-square statistics on LD scores, typically computed from an external reference panel such as 1000 Genomes. Two of its assumptions fail in admixed cohorts: no matched reference panel exists, and long-range admixture LD violates the assumption that LD is negligible beyond a short genomic window.

cov-LDSC [(Luo et al. 2021)](https://pubmed.ncbi.nlm.nih.gov/33987664/) addresses both. Principal components are projected out of the genotypes **before** LD is computed, so that the LD scores are adjusted for the same covariates included in the association model. LD scores are then computed in-sample rather than from an external panel.

Practical consequences for implementation:

- **Window size.** In European samples, mean LD scores plateau beyond ~1 cM. In admixed populations they continue to rise without covariate adjustment; after adjustment they plateau at ~20 cM. A 20-cM window is the recommended default.
- **Number of PCs.** Ten PCs is the published recommendation, though the appropriate number depends on the structure of the specific cohort and should be checked empirically.
- **Subsampling.** In-sample LD scores can be computed on a random subset of individuals. However in order to use the same genetic data for all the main analysis, we are planning on using the full cohort as its own reference for the LD scores. Additionally, filtering the genetic data for the whole cohort and use it in the reference panel will also serve as the basis for running a quick GWAS with plink in the same universe of samples and SNPs. Nonetheless, if it becomes computationally unscalable, we can revise this and look for an alternative. 
- **MAF filter.** Restrict to MAF > 0.01 for LD score computation.
- **Estimand.** Because LD is computed from array genotypes rather than refernce panel sequence data, cov-LDSC targets the same estimand as GCTA ($$h²_g$$) rather than LDSC's usual $$h²_common$$. This makes it directly more comparable to the GRM-based estimate below.


## Step 0. Project Structure 

I will start by creating our working directory

<details>
<summary>terminal example: creating the working directory and the genetic_data folder</summary>

```console
(base) [duarte@node1 ~]$ mkdir working_directory_h2
(base) [duarte@node1 ~]$ #then, get insde and create the folder to store genetic data 
(base) [duarte@node1 ~]$ cd working_directory_h2/
(base) [duarte@node1 working_directory_h2]$ mkdir genetic_data
(base) [duarte@node1 working_directory_h2]$ # should look like this:
(base) [duarte@node1 working_directory_h2]$ ls
genetic_data
(base) [duarte@node1 working_directory_h2]$ pwd
/home/duarte/working_directory_h2

```

</details>

Notice the (base) only means that my default conda environment is activated, we will get to create a couple ones for the analysis that are needed to run. But more on that when we get to those sections. 


## Step 1. Get your imputed data 

For our analysis we will need to access **imputed** and **only genotyped data**. The example here is working with .vcf files that were returned from TOPMed Imputation Server, that we used for imputing the data from [LARGE-PD](https://large-pd.org/), used in the papers from [LARGE-PD phase 1](https://pubmed.ncbi.nlm.nih.gov/34227697/) and [LARGE-PD phase 2] (https://pubmed.ncbi.nlm.nih.gov/40791673/).

For the purposes of this project and in order to maximize the sample size for heritability estimation, we merge both phases and were imputed together with [TOPMed] (https://imputation.biodatacatalyst.nhlbi.nih.gov/#!).

Here we are going to create a symlink [(like a type of shortcut in Linux)] (https://www.reddit.com/r/linux/comments/1a3xxq/someone_please_explain_symlinks_to_me_i_know/?rdt=51432) to where the actual imputed data sits, you can implement a similar approach to ensure that the downstream path architecture holds. 

<details>
<summary>terminal example: creating the symlink to the imputed data and checking that it resolves</summary>

```console
(base) [duarte@node1 working_directory_h2]$ pwd
/home/duarte/working_directory_h2
(base) [duarte@node1 working_directory_h2]$ ls
genetic_data
(base) [duarte@node1 genetic_data]$ # here is where we store our imputed data:
(base) [duarte@node1 genetic_data]$ ls /home/duarte/backups/ridiculously_large_path/imputed_data 
(base) [duarte@node1 genetic_data]$ ln -s /home/duarte/backups/ridiculously_large_path/imputed_data ./imputed
(base) [duarte@node1 genetic_data]$ # you will only have to change the actual path for your imputed data, keep the ./ imputed folder name as it is (this is important) 
(base) [duarte@node1 genetic_data]$ # then, you will see the symlink and try to ls it
(base) [duarte@node1 genetic_data]$ ls 
imputed
(base) [duarte@node1 genetic_data]$ ls imputed/
22test.imiss           chr11.dose.vcf.gz      chr_12.zip.md5         chr_14.zip             chr16.info.gz          chr18.dose.vcf.gz.csi  chr1.dose.vcf.gz       chr_20.zip.md5         chr_22.zip            chr3.info.gz          chr5.dose.vcf.gz.csi  chr7.dose.vcf.gz      chr_8.zip.md5            index_vcf.sh
22test.lmiss           chr11.dose.vcf.gz.csi  chr13.dose.vcf.gz      chr_14.zip.md5         chr_16.zip             chr18.info.gz          chr1.dose.vcf.gz.csi   chr21.dose.vcf.gz      chr_22.zip.md5        chr_3.zip             chr5.info.gz          chr7.dose.vcf.gz.csi  chr9.dose.vcf.gz         logs
22test.log             chr11.info.gz          chr13.dose.vcf.gz.csi  chr15.dose.vcf.gz      chr_16.zip.md5         chr_18.zip             chr1.info.gz           chr21.dose.vcf.gz.csi  chr2.dose.vcf.gz      chr_3.zip.md5         chr_5.zip             chr7.info.gz          chr9.dose.vcf.gz.csi     qc_report.txt
22test.nosex           chr_11.zip             chr13.info.gz          chr15.dose.vcf.gz.csi  chr17.dose.vcf.gz      chr_18.zip.md5         chr_1.zip              chr21.info.gz          chr2.dose.vcf.gz.csi  chr4.dose.vcf.gz      chr_5.zip.md5         chr_7.zip             chr9.info.gz             quality-control.html
chr10.dose.vcf.gz      chr_11.zip.md5         chr_13.zip             chr15.info.gz          chr17.dose.vcf.gz.csi  chr19.dose.vcf.gz      chr_1.zip.md5          chr_21.zip             chr2.info.gz          chr4.dose.vcf.gz.csi  chr6.dose.vcf.gz      chr_7.zip.md5         chr_9.zip                statistics
chr10.dose.vcf.gz.csi  chr12.dose.vcf.gz      chr_13.zip.md5         chr_15.zip             chr17.info.gz          chr19.dose.vcf.gz.csi  chr20.dose.vcf.gz      chr_21.zip.md5         chr_2.zip             chr4.info.gz          chr6.dose.vcf.gz.csi  chr8.dose.vcf.gz      chr_9.zip.md5            unzip.err
chr10.info.gz          chr12.dose.vcf.gz.csi  chr14.dose.vcf.gz      chr_15.zip.md5         chr_17.zip             chr19.info.gz          chr20.dose.vcf.gz.csi  chr22.dose.vcf.gz      chr_2.zip.md5         chr_4.zip             chr6.info.gz          chr8.dose.vcf.gz.csi  getTopmed_noAFcheck.err  unzip.out
chr_10.zip             chr12.info.gz          chr14.dose.vcf.gz.csi  chr16.dose.vcf.gz      chr_17.zip.md5         chr_19.zip             chr20.info.gz          chr22.dose.vcf.gz.csi  chr3.dose.vcf.gz      chr_4.zip.md5         chr_6.zip             chr8.info.gz          getTopmed_noAFcheck.out  unzip.sh
chr_10.zip.md5         chr_12.zip             chr14.info.gz          chr16.dose.vcf.gz.csi  chr18.dose.vcf.gz      chr_19.zip.md5         chr_20.zip             chr22.info.gz          chr3.dose.vcf.gz.csi  chr5.dose.vcf.gz      chr_6.zip.md5         chr_8.zip             getTopmed_noAFcheck.sh
```

</details>

Here we have confirmed that we have the shortcut active for further analysis.


## Step 2. Get your only genotyped data

Only genotyped data is important for the computation of principal components. 

For this step we are assuming that your genotype array data has been already called and QCed. The basic parameters expected and more relevant information about how to perform QC in admixed populations with your whole cohort together is described elsewhere (https://github.com/MataLabCCF/GWASQC). Since heritability estimations are more stable at bigger sample sizes, having the genetic data of your cohort altogether will be ideal. If your genotyped data was QCed with other frameworks that split individuals by ancestry (like [genotools] (https://github.com/dvitale199/GenoTools) and imputed by ancestry, we can discuss and brainstorm how to proceed. 

Your genotype data could be in plink or vcf format, here I'm giving an example on how extract the "only typed" variants form the imputed data (in this case from topmed), but if you already have the genotyped data in a different set of files it is also fine, just create the symlink in the "genotype" folder to the genetic data like we did with the imputed. 

Example for submitting 22 different jobs if you are working in a cluster. Here I'm using a python script to submit the 22 jobs for me, but also a slurm array would work, I just personally find it more useful to submit individual jobs for debugging purposes. 

Before running the script, make sure you have the latest version of plink2 installed. Here in my example I'm downloading plink2 and plink1 in bundle since we are going to need them both for further analysis. And in order to not mess up or interfere with other possible versions of your programs here I'm using relative paths. 

<details>
<summary>terminal example: downloading plink2 and plink1 into ./programs and checking both executables</summary>

```console
(base) [duarte@node1 working_directory_h2]$ pwd
/home/duarte/working_directory_h2
(base) [duarte@node1 working_directory_h2]$ ls 
genetic_data
(base) [duarte@node1 working_directory_h2]$ # create a folder to store the programs that we need 
(base) [duarte@node1 working_directory_h2]$ mkdir programs 
(base) [duarte@node1 working_directory_h2]$ cd programs/
(base) [duarte@node1 programs]$ # download plink2 and then plink1, unzip them and test that the executables are working
(base) [duarte@node1 programs]$ wget https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260808.zip 
Resolving s3.amazonaws.com (s3.amazonaws.com)... 52.216.56.240, 52.217.141.232, 52.217.192.120, ...
Connecting to s3.amazonaws.com (s3.amazonaws.com)|52.216.56.240|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 7593268 (7.2M) [application/zip]
Saving to: ‘plink2_linux_x86_64_20260808.zip’

plink2_linux_x86_64_20260808.zip                                                      100%[========================================================================================================================================================================================================================>]   7.24M  46.7MB/s    in 0.2s    

‘plink2_linux_x86_64_20260808.zip’ saved [7593268/7593268]

(base) [duarte@node1 programs]$ unzip plink2_linux_x86_64_20260808.zip 
Archive:  plink2_linux_x86_64_20260808.zip
  inflating: plink2                  
  inflating: vcf_subset              
  inflating: intel-simplified-software-license.txt  
(base) [duarte@node1 programs]$ ./plink2 --version 
PLINK v2.0.0-a.7.3LM 64-bit Intel (8 Aug 2026)
(base) [duarte@node1 programs]$ # plink2 is working, now download plink1 
(base) [duarte@node1 programs]$ wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip
Resolving s3.amazonaws.com (s3.amazonaws.com)... 16.15.213.131, 16.15.253.242, 52.216.24.22, ...
Connecting to s3.amazonaws.com (s3.amazonaws.com)|16.15.213.131|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 6321240 (6.0M) [application/zip]
Saving to: ‘plink_linux_x86_64_20250819.zip’

plink_linux_x86_64_20250819.zip                                                       100%[========================================================================================================================================================================================================================>]   6.03M  36.5MB/s    in 0.2s    

‘plink_linux_x86_64_20250819.zip’ saved [6321240/6321240]

(base) [duarte@node1 programs]$ unzip plink_linux_x86_64_20250819.zip 
Archive:  plink_linux_x86_64_20250819.zip
  inflating: plink                   
  inflating: LICENSE                 
  inflating: toy.ped                 
  inflating: toy.map                 
  inflating: prettify                
(base) [duarte@node1 programs]$ ./plink --version 
PLINK v1.9.0-b.7.11 64-bit (19 Aug 2025)
(base) [duarte@node1 programs]$ #plink 1 is working as well. when calling them in different scripts, we will use this relative path so it doesn't interfere with your other possible versions 
(base) [duarte@node1 programs]$ ls
intel-simplified-software-license.txt  LICENSE  plink  plink2  plink2_linux_x86_64_20260808.zip  plink_linux_x86_64_20250819.zip  prettify  toy.map  toy.ped  vcf_subset
(base) [duarte@node1 programs]$ # now, lets create the folder of the genotyped genetic data
(base) [duarte@node1 programs]$ cd ../genetic_data/
(base) [duarte@node1 genetic_data]$ mkdir genotyped
(base) [duarte@node1 genetic_data]$ ls
genotyped  imputed
(base) [duarte@node1 genetic_data]$ cd genotyped/

```

</details>

Here is the example to extract the genotyped data from the imputed files:

Just make sure that the only genotyped variants INFO field corresponds to the actual information in your files. In my case, TOPMed's vcfs carry the field "TYPED" for the only genotyped variants. 

And an important note, sometimes imputation servers miss some variant IDs we will enforce the flag --set-all-var-ids (paired with a safety value of new id max allele length). 

As a python script to submit 22 different jobs

<details>
<summary>script: <code>extract_genotyped_data.py</code> — writes one SLURM .pbs per chromosome and submits all 22</summary>

```python
import os

runFolder = f"./"
os.makedirs(f"{runFolder}/logs", exist_ok=True)

for i in range(1,23):
    mem = "30G" if i in {1,2,3,4,5} else "20G"
    vcfFileIn = f"../imputed/chr{i}.dose.vcf.gz" 
    plink2FileOut = f"chr{i}_typed" 
    threads = 8
    email = f"your_mail@email.org"
    partition = f"defq"

    fileSbatch = open(f"{runFolder}/logs/{i}extract_onlyTyped.pbs", "w")

    fileSbatch.write(f"#!/bin/sh\n")
    fileSbatch.write(f"#SBATCH --mail-type=END,FAIL\n")
    fileSbatch.write(f"#SBATCH --mail-user={email}\n")
    fileSbatch.write(f"#SBATCH --ntasks=1\n")
    fileSbatch.write(f"#SBATCH --cpus-per-task={threads}\n")
    fileSbatch.write(f"#SBATCH --mem={mem}\n")
    fileSbatch.write(f"#SBATCH --partition={partition}\n")
    fileSbatch.write(f"#SBATCH --job-name={i}extract_onlyTyped\n") 
    fileSbatch.write(f"#SBATCH -o {runFolder}/logs/{i}extract_onlyTyped.out\n")
    fileSbatch.write(f"#SBATCH -e {runFolder}/logs/{i}extract_onlyTyped.err\n\n")


    fileSbatch.write(f"../../programs/plink2 --threads {threads} --vcf {vcfFileIn} --require-info TYPED --make-pgen --out {plink2FileOut} --set-all-var-ids @:#:\\$r:\\$a --new-id-max-allele-len 1000 \n")


    fileSbatch.close()
    os.system(f"sbatch {runFolder}/logs/{i}extract_onlyTyped.pbs")

```

</details>

For running it do:

<details>
<summary>terminal example: saving the python script with nano and submitting the 22 jobs</summary>

```console
(base) [duarte@node1 genotyped]$ # do a nano with the name of the script, hit enter, then ctrl + X, and then hit "y" to save it
(base) [duarte@node1 genotyped]$ nano extract_genotyped_data.py
(base) [duarte@node1 genotyped]$ nano extract_genotyped_data.py
(base) [duarte@node1 genotyped]$ nano extract_genotyped_data.py #then the ctrl + X, then Y,
(base) [duarte@node1 genotyped]$ ls 
extract_genotyped_data.py
(base) [duarte@node1 genotyped]$ python extract_genotyped_data.py 
Submitted batch job 6109323
Submitted batch job 6109324
Submitted batch job 6109325
Submitted batch job 6109326
Submitted batch job 6109327
Submitted batch job 6109328
Submitted batch job 6109329
Submitted batch job 6109330
Submitted batch job 6109331
Submitted batch job 6109332
Submitted batch job 6109333
Submitted batch job 6109334
Submitted batch job 6109335
Submitted batch job 6109336
Submitted batch job 6109337
Submitted batch job 6109338
Submitted batch job 6109339
Submitted batch job 6109340
Submitted batch job 6109341
Submitted batch job 6109342
Submitted batch job 6109343
Submitted batch job 6109344

```

</details>

If you prefer to work in series through your terminal here is the equivalent for loop, it might take more time to run, but it will do the job as well.

<details>
<summary>script: <code>extract_genotyped_data.sh</code> — same extraction as a sequential for loop, no scheduler needed</summary>

```bash
#!/bin/bash
set -euo pipefail

threads=8 #depends on how many you have available in your machine 
mkdir -p ./logs

for i in $(seq 22 -1 1); do
    echo "=== chr${i} ==="
    ../../programs/plink2 \
        --threads "${threads}" \
        --vcf "../imputed/chr${i}.dose.vcf.gz" \
        --require-info TYPED \
        --make-pgen \
        --out "chr${i}_typed" \
        --set-all-var-ids '@:#:$r:$a' \
        --new-id-max-allele-len 1000 \
        > "./logs/${i}extract_onlyTyped.out" 2> "./logs/${i}extract_onlyTyped.err"
done

echo "done"
```

</details>

In the terminal do

<details>
<summary>terminal example: saving and running the bash version</summary>

```console
(base) [duarte@node1 genotyped]$ nano extract_genotyped_data.sh
(base) [duarte@node1 genotyped]$ bash extract_genotyped_data.sh 
(base) [duarte@node1 genotyped]$ # if it doesn't run, try to make sure it is an executable
(base) [duarte@node1 genotyped]$ chmod +x extract_genotyped_data.sh #and then try again to bash it
```

</details>

After running either of the two options, you will end up with one PLINK2 fileset per chromosome, 22 in total.

<details>
<summary>terminal example: what the output directory looks like</summary>

```console
(base) [duarte@node1 genotyped]$ ls chr*_typed.*
chr10_typed.log   chr15_typed.psam  chr20_typed.log   chr4_typed.psam
chr10_typed.pgen  chr15_typed.pvar  chr20_typed.pgen  chr4_typed.pvar
chr10_typed.psam  chr16_typed.log   chr20_typed.psam  chr5_typed.log
chr10_typed.pvar  chr16_typed.pgen  chr20_typed.pvar  chr5_typed.pgen
chr11_typed.log   chr16_typed.psam  chr21_typed.log   chr5_typed.psam
chr11_typed.pgen  chr16_typed.pvar  chr21_typed.pgen  chr5_typed.pvar
chr11_typed.psam  chr17_typed.log   chr21_typed.psam  chr6_typed.log
chr11_typed.pvar  chr17_typed.pgen  chr21_typed.pvar  chr6_typed.pgen
chr12_typed.log   chr17_typed.psam  chr22_typed.log   chr6_typed.psam
chr12_typed.pgen  chr17_typed.pvar  chr22_typed.pgen  chr6_typed.pvar
chr12_typed.psam  chr18_typed.log   chr22_typed.psam  chr7_typed.log
chr12_typed.pvar  chr18_typed.pgen  chr22_typed.pvar  chr7_typed.pgen
chr13_typed.log   chr18_typed.psam  chr2_typed.log    chr7_typed.psam
chr13_typed.pgen  chr18_typed.pvar  chr2_typed.pgen   chr7_typed.pvar
chr13_typed.psam  chr19_typed.log   chr2_typed.psam   chr8_typed.log
chr13_typed.pvar  chr19_typed.pgen  chr2_typed.pvar   chr8_typed.pgen
chr14_typed.log   chr19_typed.psam  chr3_typed.log    chr8_typed.psam
chr14_typed.pgen  chr19_typed.pvar  chr3_typed.pgen   chr8_typed.pvar
chr14_typed.psam  chr1_typed.log    chr3_typed.psam   chr9_typed.log
chr14_typed.pvar  chr1_typed.pgen   chr3_typed.pvar   chr9_typed.pgen
chr15_typed.log   chr1_typed.psam   chr4_typed.log    chr9_typed.psam
chr15_typed.pgen  chr1_typed.pvar   chr4_typed.pgen   chr9_typed.pvar
```

</details>

After this, you can proceed to merge all the 22 files into a single file that could be used in downstream analysis:

Here an example using this .sh script 

<details>
<summary>script: <code>merge_onlyTyped.sh</code> — merges the 22 per-chromosome filesets into one genome-wide pgen fileset</summary>

```bash
#!/bin/bash
set -euo pipefail

# Merges the 22 per-chromosome filesets produced by run_extract_onlyTyped.sh
# (chr1_typed .. chr22_typed) into a single genome-wide pgen fileset.
# Run from the same directory as run_extract_onlyTyped.sh.

threads=8 #depends on how many you have available in your machine
mkdir -p ./logs

out_prefix="allchr_typed"
merge_list="./merge_list_onlyTyped.txt"

# Build the fileset list in ascending chromosome order, one prefix per line,
# and fail early if any chromosome is missing or was left incomplete.
: > "${merge_list}"
for i in $(seq 1 22); do
    for ext in pgen pvar psam; do
        if [[ ! -s "chr${i}_typed.${ext}" ]]; then
            echo "ERROR: missing or empty chr${i}_typed.${ext}" >&2
            exit 1
        fi
    done
    echo "chr${i}_typed" >> "${merge_list}"
done

echo "=== merging chr1-22 into ${out_prefix} ==="
../../programs/plink2 \
    --threads "${threads}" \
    --pmerge-list "${merge_list}" pfile \
    --out "${out_prefix}" \
    > "./logs/merge_onlyTyped.out" 2> "./logs/merge_onlyTyped.err"

echo "=== merged fileset: ${out_prefix}.pgen / .pvar / .psam ==="
echo "variants: $(grep -cv '^#' "${out_prefix}.pvar")"
echo "samples:  $(grep -cv '^#' "${out_prefix}.psam")"

echo "done"
```

</details>

In the terminal do

<details>
<summary>terminal example: saving and running the merge script</summary>

```console
(base) [duarte@node1 genotyped]$ nano merge_onlyTyped.sh
(base) [duarte@node1 genotyped]$ bash merge_onlyTyped.sh 
(base) [duarte@node1 genotyped]$ ls *all*
allchr_typed.log allchr_typed.pgen allchr_typed.psam allchr_typed.pvar

```

</details>

## Step 3. Get your phenotype and covariate data ready 

For this step we are assuming that you have your phenotype and covariate data in a single file without missing data and the phenotype of PD coded as control = 1, case =2; and that you have the basis covariates of Sex (male = 1, female = 2) and Age (quantitative variable).

Here is an example on how the phenotype file should look like:

<details>
    <summary>covar.txt</summary>
IID	SEX	DISEASE	AGE <br>
sample1ID	2	1	38 <br>
sample2ID	1	1	33 <br>
sample3ID	1	1	45 <br>
sample4ID	2	2	66 <br>
sample5ID	2	2	72
</details>

Once you made sure that this is the shape of your covariate file, we can move onto the next step that involves PCA and such. Hence, we need to make sure our covariate file is at hand in the project structure.

Here an example:
```
(base) [duarte@node1 working_directory_h2]$ mkdir pca_and_such
(base) [duarte@node1 working_directory_h2]$ cd pca_and_such/
(base) [duarte@node1 pca_and_such]$ # copy your covariate file to this folder, or provide the absolute path to the next steps 
(base) [duarte@node1 pca_and_such]$ # here I put my covariate file directly
(base) [duarte@node1 pca_and_such]$ nano covar.txt
(base) [duarte@node1 pca_and_such]$ ls covar.txt 
covar.txt
```


## Step 4. Run the scripts for generating PCs, do relationship control and downstream covariate formating. 

For running computing the in-sample reference LD scores, we need to adjust them for the principal components, and we will use the only typed data for this purpose. 

Here we choose the R packages PC-Air and PC-Relate (available through GENESIS) to compute the PCs (and GRMs, see below), mainly because URPs can be very admixed cohorts and classical methods for computing PCs and GRMs tend to have problems when a dataset contains a diverse setting of genetic admixture. The yield of this approach compared to other more "classical methods" is not that much different in terms of adjustment for population structure (both can control spurious association in single variant association testing, see: https://doi.org/10.1101/2025.05.27.25328444, results section). However, since here we are interested in capturing as accurate as possible the h2 estimates, methods that account for population structure in an ancestry-aware approach, are highly encouraged. 

In order to compute the PCs with GENESIS we will need a specific conda environment, here an example (note, it will take a couple minutes to install)

```
(base) [duarte@node1 working_directory_h2]$ conda create -n genesis \
  -c conda-forge -c bioconda \
  r-base=4.5 \
  r-data.table \
  r-tidyverse \
  r-optparse \
  bioconductor-genesis \
  bioconductor-snprelate \
  bioconductor-gdsfmt \
  bioconductor-seqarray \
  bioconductor-seqvartools \
  bioconductor-gwastools \
  plink2 \
  bcftools \
  htslib \
  -y

#activate your environment
(base) [duarte@node1 working_directory_h2]$ conda activate genesis
(genesis) [duarte@node1 working_directory_h2]$ 


```
Once your environment is active, lets move to run PCA first 

### Step 4.a run the PCA pipeline 

<details>
<summary><strong>View script: <code>pcs_pipeline_h2_project.r</code></strong></summary>

<br>

```r
library(SNPRelate)
library(GENESIS)
library(GWASTools)
library(gdsfmt)
library(BiocParallel)
library(Matrix)

######################################
# COMMAND-LINE ARGUMENTS 
######################################

cli       <- commandArgs(trailingOnly = TRUE)
show_help <- any(cli %in% c("--help", "-h"))
cli       <- cli[!cli %in% c("--help", "-h")]

# Split c("--a", "1", "--b", "2") into list(a = "1", b = "2"):
# odd positions are the names, even positions the values.
cli_args <- list()
if (length(cli) > 0) {
  if (length(cli) %% 2 != 0 || !all(grepl("^--", cli[c(TRUE, FALSE)])))
    stop("Arguments must be given as --name value pairs. Run with --help.")
  cli_args        <- as.list(cli[c(FALSE, TRUE)])
  names(cli_args) <- sub("^--", "", cli[c(TRUE, FALSE)])
}

# Every default that passes through arg() is recorded here, so --help can list
# the options and unknown flags can be caught.
arg_defaults <- list()

# The --name given on the command line, or default when it was not given.
# A numeric default keeps values numeric. so maf 0.05 stays as a number
arg <- function(name, default) {
  arg_defaults[[name]] <<- default
  if (is.null(cli_args[[name]])) return(default)
  value <- if (is.numeric(default)) as.numeric(cli_args[[name]]) else cli_args[[name]]
  message("[arg] ", name, " = ", value, "   (default: ", default, ")")
  value
}

########
# CONFIG 
########

# Output folder — every file this script writes lands here 
# Created as soon as the arguments are resolved. Inputs are read from wherever they are; outputs never escape this folder.
out_dir <- arg("out_dir", "outFolder_pca_andSuch")


# Every output name below is a bare file name: out_path puts it in out_dir.
out_path <- function(name) file.path(out_dir, name)

# Log file — everything the script prints is copied here (see LOGGING below).
log_file <- out_path(arg("log_file", "pipeline_run.log"))

# Input: raw genotype data prefix — pipeline entry point
input_plink_file   <- arg("input_plink_file",  "inputPlinkFile")
input_format  <- arg("input_format", "pfile")

# The plink flag matching input_format; reused by every step below
input_flag <- if (input_format == "bfile") "--bfile" else "--pfile"

# Input: sample covariate table, merged with the PCs at the end
# Tab-separated with a header; its first column is the sample ID (renamed to IID on read just in case). Read-only, so it is NOT placed inside out_dir
covar_file <- arg("covar_file", "covar.txt")

# Variant ID re-parsing (plink2) — first thing done to the raw input
# The raw input is copied once with rebuilt variant IDs (chr:pos:ref:alt) and every downstream step reads that copy, never the raw input itself.
# The copy keeps the input's own format: pfile in -> pfile out, bfile -> bfile.
#this is important since this list is the one we are going to use when filtering imputed data (see the github for more details) 

var_id_template       <- arg("var_id_template", "@:#:$r:$a")
new_id_max_allele_len <- arg("new_id_max_allele_len", 1000)
prep_setids           <- out_path(arg("prep_setids", "input_setVarIDs"))

# Data preparation 
plink2_bin          <- arg("plink2_bin", "../programs/plink2")
highld_regions_file <- out_path(arg("highld_regions_file", "highLD_regions_grindeLab_hg38.tsv"))
maf_min             <- arg("maf_min",      0.01)
indep_window        <- arg("indep_window", 200)
indep_step          <- arg("indep_step",   50)
indep_r2            <- arg("indep_r2",     0.2)

# Optional: exclude high-LD / long-range-LD regions before PCA.
#   "yes" -> write the regions tsv and drop those regions (extra prep step)
#   "no"  -> keep all regions; the tsv and the exclusion step are skipped
removeHighLDregions <- arg("removeHighLDregions", "yes")

# Intermediate plink2 output prefixes. Each branch uses its own names so the
# high-LD-removed and high-LD-kept variants never overwrite one another.
prep_ld <- out_path(arg("prep_ld", "ld"))
if (removeHighLDregions == "yes") {
  prep_noHighLD <- out_path(arg("prep_noHighLD", "outBfile_noHighLD_regions"))
  prep_maf      <- out_path(arg("prep_maf",    "outBfile_noHighLD_regions_commonVars"))
  prep_pruned   <- out_path(arg("prep_pruned", "outBfile_noHighLD_regions_commonVars_indep-pairwise"))
} else {
  prep_maf      <- out_path(arg("prep_maf",    "outBfile_withHighLD_regions_commonVars"))
  prep_pruned   <- out_path(arg("prep_pruned", "outBfile_withHighLD_regions_commonVars_indep-pairwise"))
}
# KING uses the re-ID'd data as bed. If the input was already bed, the re-ID'd copy is bed too and is used directly; if it was a pfile, it is converted to this prefix in prep step 5
king_bed      <- if (input_format == "bfile") prep_setids else out_path(arg("king_bed", "inputPlinkFile_bed"))

# PLINK prefixes consumed by the R analysis (derived from prep outputs)
pca_plink_prefix  <- prep_pruned     
king_plink_prefix <- king_bed        

# Outputs: GDS files 
pca_gds  <- out_path(arg("pca_gds",  "onlyTyped4PCA.gds"))
king_gds <- out_path(arg("king_gds", "onlyTyped4King.gds"))

# Outputs: RDS files
king_mat_rds     <- out_path(arg("king_mat_rds",     "KINGmat_rds"))
pcair_r1_rds     <- out_path(arg("pcair_r1_rds",     "mypcair_1round_results_rds_1stRound"))
pcrel_r1_rds     <- out_path(arg("pcrel_r1_rds",     "mycprelate_1stround"))
pcrel_r1_mat_rds <- out_path(arg("pcrel_r1_mat_rds", "mypcrel_1round_mat_rds"))
pcair_r2_rds     <- out_path(arg("pcair_r2_rds",     "mypcair_r2_results_rds"))
#pcrel_r2_rds     <- out_path("mypcrel_r2_rds")
#pending to make the pcrel_r2_mat_rds that is in kinship scale 2, for the h2 project if other methods are used that require grms

# Outputs: TSV tables handed to downstream tools 
# Long-format pairwise kinship (ID1 ID2 kinship, no header) -> NAToRA input in downstream steps (outside this script)
kinship_pairs_tsv  <- out_path(arg("kinship_pairs_tsv",  "pcrelate_r1_kinship_pairs.tsv"))
# PC-AiR round-2 PCs as covariates (sample.id PC1 PC2 ...), with header
pcair_r2_covar_tsv <- out_path(arg("pcair_r2_covar_tsv", "pcair_r2_covariates.tsv"))
# The above left-joined onto covar_file, keyed on IID
pcair_r2_covar_merged_tsv <- out_path(arg("pcair_r2_covar_merged_tsv", "pcair_r2_covariates_merged.tsv"))

# Analysis parameters 
n_cores           <- arg("n_cores",           1)
n_pcs             <- arg("n_pcs",             10)
snp_block_size    <- arg("snp_block_size",    10000)
sample_block_size <- arg("sample_block_size", 7000)

# Output filename patterns (sprintf style)
# PC-pair scatter plots: %02d, %02d = the two PC numbers on the axes
pcpair_png_pattern_r1 <- out_path(arg("pcpair_png_pattern_r1", "plot_PC%02d_PC%02d.png"))
pcpair_png_pattern_r2 <- out_path(arg("pcpair_png_pattern_r2", "plot_PC%02d_PC%02d_2ndround.png"))
# SNP–PC correlation plots/tables: %02d = PC number
snpcorr_png_pattern_r1 <- out_path(arg("snpcorr_png_pattern_r1", "snpcorr_1stRound_PC%02d.png"))
snpcorr_tsv_pattern_r1 <- out_path(arg("snpcorr_tsv_pattern_r1", "snpcorr_1stRound_PC%02d.tsv"))
snpcorr_png_pattern_r2 <- out_path(arg("snpcorr_png_pattern_r2", "snpcorr_2ndRound_PC%02d.png"))
snpcorr_tsv_pattern_r2 <- out_path(arg("snpcorr_tsv_pattern_r2", "snpcorr_2ndRound_PC%02d.tsv"))

# Plot dimensions (pixels / dpi) 
pcpair_plot_width  <- arg("pcpair_plot_width",  1500)
pcpair_plot_height <- arg("pcpair_plot_height", 1500)
pcpair_plot_res    <- arg("pcpair_plot_res",    300)
snpcorr_plot_width  <- arg("snpcorr_plot_width",  1800)
snpcorr_plot_height <- arg("snpcorr_plot_height", 700)
snpcorr_plot_res    <- arg("snpcorr_plot_res",    150)

############################################################
# RESOLVE THE COMMAND LINE — runs once, before any real work
############################################################

# --help lists every option with the default it just registered above.
if (show_help) {
  cat("\nUsage: Rscript PCs_pipeline_h2_4project.r [--name value ...]\n\n")
  cat("Options, with their defaults. Output names are bare file names:\n")
  cat("they are always written inside --out_dir.\n\n")
  for (name in names(arg_defaults))
    cat(sprintf("  --%-27s %s\n", name, arg_defaults[[name]]))
  cat("\n")
  quit(save = "no")
}

# A misspelled flag would otherwise be ignored in silence, so refuse to start.
unknown <- setdiff(names(cli_args), names(arg_defaults))
if (length(unknown) > 0)
  stop("Unknown argument(s): --", paste(unknown, collapse = ", --"),
       "\nRun with --help to see the valid ones.")

# Output folder — every file this script writes lands here.
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#########################################################
# LOGGING — everything the script prints goes to log_file
#########################################################
# Started here, right after out_dir exists and the arguments are known, so the
# whole run is recorded. Two sinks are needed because R has two streams:
#   stdout -> print(), cat(), auto-printed results (split = TRUE keeps them on screen too)
#   stderr -> message(), warning(), and the error that aborts the run
# Closed at the very bottom of the script.
# Note: plink2 writes straight to the terminal, not through R, so its output is
# not in here — it is in the .log file plink2 writes next to each --out prefix.
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
message("[log] ", format(Sys.time()), " — writing this run to ", log_file)
message("[log] command line: Rscript PCs_pipeline_h2_4project.r ",
        paste(commandArgs(trailingOnly = TRUE), collapse = " "))

#########################################################
# HELPERS — generic names  for common files 
#########################################################

# Convert a PLINK bed/bim/fam triple (given by prefix) to a GDS file,
# then open it, print a summary, and close it again.
# Like the plink2 prep and the KING matrix, an existing output is reused: if the
# .gds file is already there the conversion is skipped and only the summary is
# printed, so re-runs are cheap. The check is on the file name alone, so delete
# the .gds (or use a different --pca_gds / --king_gds) whenever the plink
# outputs it was built from have changed, otherwise the stale one is reused.
bed2gds <- function(prefix, out_gds) {
  if (file.exists(out_gds)) {
    message("GDS already present — skipping conversion: ", out_gds)
  } else {
    snpgdsBED2GDS(
      bed.fn    = paste0(prefix, ".bed"),
      bim.fn    = paste0(prefix, ".bim"),
      fam.fn    = paste0(prefix, ".fam"),
      out.gdsfn = out_gds
    )
  }
  gds <- snpgdsOpen(out_gds)
  snpgdsSummary(out_gds)
  snpgdsClose(gds)
}

# Does a PLINK bed/bim/fam triple already exist for this prefix?
bed_exists <- function(prefix) all(file.exists(paste0(prefix, c(".bed", ".bim", ".fam"))))

# Run one plink2 command; stop the pipeline if it exits non-zero.
run_plink <- function(args) {
  args <- as.character(args)
  message("plink2 ", paste(args, collapse = " "))
  status <- system2(plink2_bin, args)
  if (status != 0L)
    stop("plink2 failed (exit ", status, "): ", paste(args, collapse = " "))
}


###########################
# DATA PREPARATION — plink2 
###########################
# Step 0 rewrites all variant IDs as chr:pos:ref:alt (prep_setids); 
#everything after it works from that copy. It then builds the two bed datasets the analysis consumes:
#   pca_plink_prefix  (optionally high-LD removed, maf-filtered, LD-pruned) -> PCA
#   king_plink_prefix (re-ID'd data converted to bed)                       -> KING
# High-LD removal is controlled by removeHighLDregions in CONFIG.
# Skipped automatically when both already exist, so re-runs are cheap.

# plink2 data prep
if (bed_exists(pca_plink_prefix) && bed_exists(king_plink_prefix)) {
  message("Prep outputs already present — skipping plink2 data preparation.")
} else {
  # 0) Rebuild every variant ID as chr:pos:ref:alt on the raw input, keeping the input's format (pfile -> pgen, bfile -> bed). 
  #    Every step after this one reads prep_setids, so the raw input is touched exactly once.
  #    shQuote protects '@:#:$r:$a' from the shell system2 runs the command in.
  run_plink(c(input_flag, input_plink_file,
              "--set-all-var-ids", shQuote(var_id_template),
              "--new-id-max-allele-len", new_id_max_allele_len,
              "--write-snplist", "allow-dups",
              if (input_format == "bfile") "--make-bed" else "--make-pgen",
              "--out", prep_setids))
#this snplist is the one that would be used later alongside the hapmap var list to subset the imputed data for running the h2 estimation

  # 1) Optionally drop high-LD / long-range-LD regions from the re-ID'd data -> bed.
  #    When removeHighLDregions = "no", this whole step (and the tsv) is skipped and the maf filter reads straight from the re-ID'd input instead.
  if (removeHighLDregions == "yes") {
    # High-LD regions to exclude (Grinde lab, hg38, https://github.com/GrindeLab/PCA).
    # Columns: chrom  st  end  label   (tab-separated, no header)
    highld_regions <- c(
      "1\t47761741\t51822307\tanderson1_price1_michigan1",
      "2\t129125957\t139525961\ttopmedLCT_michigan3_priveceliac1_price3_privepopres1_raskalct",
      "2\t182309767\t189427029\tmichigan4_anderson3_price4",
      "3\t47483506\t49987563\tanderson4_michigan5_price5",
      "3\t83368159\t86868160\tanderson5_michigan6_price6",
      "3\t161899518\t163699518\tpriveceliac4",
      "5\t98636396\t101136397\tmichigan9_price9",
      "5\t129636408\t132636409\tanderson7_michigan10_price10",
      "5\t136136412\t139136412\tmichigan11_price11",
      "6\t23691793\t38924246\tpriveceliac2_topmedMHC_raskahla_fellay2_anderson8_michigan12_price12_privepopres2",
      "6\t139637170\t142137170\tanderson10_michigan14_price14",
      "8\t6455071\t13598120\tpriveceliac3_privepopres3_topmedinversion_anderson12_fellay3_michigan16_price16_raskainv",
      "8\t110918595\t113918595\tanderson14_michigan18_price18",
      "11\t88127184\t91127184\tanderson16_michigan21_price21",
      "12\t110577812\t113099475\tprice23_michigan23",
      "14\t47061047\t47961047\tpriveceliac5",
      "17\t42394456\t46567318\ttopmedinversion",
      "20\t33948533\t36438183\tanderson18_michigan24_price24"
    )
    writeLines(highld_regions, highld_regions_file)

    run_plink(c(input_flag, prep_setids,
                "--exclude", "bed1", highld_regions_file,
                "--make-bed", "--out", prep_noHighLD))
    maf_input <- c("--bfile", prep_noHighLD)
  } else {
     # maf filter reads the re-ID'd input directly in case high-LD regions are kept
    message("removeHighLDregions = 'no' — keeping high-LD regions; tsv and exclusion step skipped.")
    maf_input <- c(input_flag, prep_setids)
  }

  # 2) Minor-allele-frequency filter
  run_plink(c(maf_input,
              "--maf", maf_min,
              "--make-bed", "--out", prep_maf))

  # 3) LD pruning: build the list of variants to prune
  run_plink(c("--bfile", prep_maf,
              "--indep-pairwise", indep_window, indep_step, indep_r2,
              "--out", prep_ld))

  # 4) Apply the prune list -> LD-pruned bed (pca_plink_prefix)
  run_plink(c("--bfile", prep_maf,
              "--exclude", paste0(prep_ld, ".prune.out"),
              "--make-bed", "--out", prep_pruned))

  # 5) Convert the re-ID'd pfile to bed for KING / SNP-correlations (king_plink_prefix).
  #    Skipped when the input was already bed — king_bed then points at the re-ID'd bed directly.
  if (input_format != "bfile") {
    run_plink(c("--pfile", prep_setids,
                "--make-bed", "--out", king_bed))
  }
}
#################################################################################
# DATA PREPARATION — convert plink files to gds files that R can read
#################################################################################
# Each conversion is skipped when its .gds already exists (see bed2gds above).
bed2gds(pca_plink_prefix,  pca_gds)
bed2gds(king_plink_prefix, king_gds)

#######################################################################################################################################
# Run KING — under their robust method, the kinship coefficient is computed to identify the set of unrelated individuals
#######################################################################################################################################
# KING is time-consuming: reuse a matrix from a previous run if one exists,
# otherwise compute it and cache it to king_mat_rds.
if (file.exists(king_mat_rds)) {
  message("Reusing existing KING matrix: ", king_mat_rds)
  KINGmat <- readRDS(king_mat_rds)
} else {
  message("Computing KING matrix -> ", king_mat_rds)
  gds_king <- snpgdsOpen(king_gds)

  king <- snpgdsIBDKING(gds_king, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,type=c("KING-robust"), family.id=NULL, verbose=TRUE)

  KINGmat <- king$kinship
  rownames(KINGmat) <- colnames(KINGmat) <- king$sample.id
  snpgdsClose(gds_king)

  saveRDS(KINGmat, king_mat_rds)
}

# KING matrix (square numeric matrix)
 # n_samples x n_samples
dim(KINGmat)  
# top-left corner                  
KINGmat[1:5, 1:5]               



#########################################################
# Run the first Round of PC-AiR using the input from KING
#########################################################
# Same reuse rule as the KING matrix above: if pcair_r1_rds is already there it
# is loaded and PC-AiR is not rerun. Delete it to force a fresh run — in
# particular after changing --n_pcs, which the file name does not record.
if (file.exists(pcair_r1_rds)) {
  message("Reusing existing PC-AiR round 1: ", pcair_r1_rds)
  mypcair_1round <- readRDS(pcair_r1_rds)
} else {
  message("Running PC-AiR round 1 -> ", pcair_r1_rds)
  #load the gds file
  geno_reader <- GdsGenotypeReader(filename = pca_gds)
  geno_data <- GenotypeData(geno_reader)

  #load the KING matrix
  KINGmat <- readRDS(king_mat_rds)

  #run pcair
  mypcair_1round <- pcair(geno_data, kinobj = KINGmat, divobj = KINGmat, num.cores=n_cores, eigen.cnt = n_pcs)
  saveRDS(mypcair_1round, pcair_r1_rds)
  close(geno_data)
}
summary(mypcair_1round)

#plot the pcs
for (i in seq(1, n_pcs - 1, by = 2)) {
  png(sprintf(pcpair_png_pattern_r1, i, i + 1),
      width = pcpair_plot_width, height = pcpair_plot_height, res = pcpair_plot_res)
  plot(mypcair_1round, vx = i, vy = i + 1)
  dev.off()
}


###########################################################
# Run the first Round of correlations of genotypes with PCs
###########################################################
#consider that is the full dataset, no ld exclusion, no ld pruninng, no maf filtering, nothing of that since we want to se the correlation genome-wide with all the genotyped variants


genofile <- snpgdsOpen(king_gds)

chr    <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
pos    <- read.gdsn(index.gdsn(genofile, "snp.position"))
snp_id <- read.gdsn(index.gdsn(genofile, "snp.id"))

cr1 <- snpgdsPCACorr(
  pcaobj     = mypcair_1round$vectors[, 1:n_pcs],
  gdsobj     = genofile,
  eig.which  = 1:n_pcs,
  num.thread = n_cores
)

snpgdsClose(genofile)

# align annotation to exactly the SNPs snpgdsPCACorr used
snp_idx_1round     <- match(cr1$snp.id, snp_id)
chr_used_1round    <- chr[snp_idx_1round]
pos_used_1round   <- pos[snp_idx_1round]
snp_id_used_1round <- snp_id[snp_idx_1round]

# chromosome x-axis: midpoint SNP index per chromosome, labeled chr1, chr2...
chr_levels_1round <- sort(unique(chr_used_1round))
chr_mids_1round   <- tapply(seq_along(chr_used_1round), chr_used_1round, median)
chr_labels_1round <- paste0("chr", names(chr_mids_1round))

# alternating colors by chromosome
# one color per chromosome
chr_colors        <- palette()[((as.integer(chr_levels_1round) - 1) %% length(palette())) + 1]
names(chr_colors) <- chr_levels_1round
col_vec           <- chr_colors[as.character(chr_used_1round)]
for (i in 1:n_pcs) {
  abs_corr_1round <- abs(cr1$snpcorr[i, ])

  # PNG with chromosome x-axis
  png(sprintf(snpcorr_png_pattern_r1, i),
      width = snpcorr_plot_width, height = snpcorr_plot_height, res = snpcorr_plot_res)
  plot(abs_corr_1round,
       col  = col_vec,
       pch  = 20,
       cex  = 0.3,
       ylim = c(0, 1),
       xaxt = "n",
       xlab = "Chromosome",
       ylab = paste0("PC", i, " Correlation"),
       main = paste0("PC-AiR 1st round | PC", i, " — SNP correlation"))
  axis(1, at = chr_mids_1round, labels = chr_labels_1round, las = 2, cex.axis = 0.7)
  dev.off()
  # TSV sorted by absolute correlation descending
  df_1round <- data.frame(
    snp_id     = snp_id_used_1round,
    chromosome = chr_used_1round,
    position   = pos_used_1round,
    abs_corr_1round   = abs_corr_1round
  )
  df_1round <- df_1round[order(df_1round$abs_corr_1round, decreasing = TRUE), ]
  write.table(df_1round, sprintf(snpcorr_tsv_pattern_r1, i),
              sep = "\t", row.names = FALSE, quote = FALSE)
}



##################################################################################
# Run the first Round of PC-Relate using the input from the first round of PC-Air
##################################################################################
# This step writes two files and both are reused when present:
#   pcrel_r1_mat_rds — the matrix, the only thing the rest of the script needs,
#                      so when it exists everything below is skipped outright
#   pcrel_r1_rds     — the raw pcrelate object; cached on its own so a run that
#                      died between the two saves resumes at the conversion
#                      instead of redoing pcrelate, the slowest step here
# Delete both to force a full rerun (e.g. after changing --n_pcs).
if (file.exists(pcrel_r1_mat_rds)) {
  message("Reusing existing PC-Relate round 1 matrix: ", pcrel_r1_mat_rds)
  mypcrel_1round_mat <- readRDS(pcrel_r1_mat_rds)
} else {
  if (file.exists(pcrel_r1_rds)) {
    message("Reusing existing PC-Relate round 1: ", pcrel_r1_rds)
    mypcrel_1round <- readRDS(pcrel_r1_rds)
  } else {
    message("Running PC-Relate round 1 -> ", pcrel_r1_rds)
    mypcair_1round <- readRDS(pcair_r1_rds)
    geno_reader <- GdsGenotypeReader(filename = pca_gds)
    geno_data   <- GenotypeData(geno_reader)
    geno_iter   <- GenotypeBlockIterator(geno_data, snpBlock = snp_block_size)

    mypcrel_1round <- pcrelate(geno_iter, pcs = mypcair_1round$vectors[, 1:n_pcs],
                                sample.block.size = sample_block_size,
                                training.set = mypcair_1round$unrels,
                                BPPARAM = BiocParallel::MulticoreParam(n_cores))
    saveRDS(mypcrel_1round, pcrel_r1_rds)
    close(geno_data)
  }

  # Convert for round 2 PC-AiR input: scaleKin=1
  # scaleKin=1 keeps the kinship-coefficient scale expected by pcair's kinobj
  mypcrel_1round_mat <- pcrelateToMatrix(mypcrel_1round,
                                         scaleKin = 1,
                                         thresh=NULL,
                                         verbose = TRUE)
  saveRDS(mypcrel_1round_mat, pcrel_r1_mat_rds)
}
#######################################################################################################################
# Run the second round of PC-AiR using the kinships from PC-Relate round 1 that are more accurate than those from KING. 
#however, the divergence signal is still required to come from KING
#######################################################################################################################
# kinobj  → PC-Relate r1 (accurate recent relatedness)
# divobj  → KING (retains divergence signal across ancestry groups)
# Reused when pcair_r2_rds exists, like every step above. Delete it to rerun.
if (file.exists(pcair_r2_rds)) {
  message("Reusing existing PC-AiR round 2: ", pcair_r2_rds)
  mypcair_r2 <- readRDS(pcair_r2_rds)
} else {
  message("Running PC-AiR round 2 -> ", pcair_r2_rds)
  KINGmat <- readRDS(king_mat_rds)
  geno_reader  <- GdsGenotypeReader(filename = pca_gds)
  geno_data    <- GenotypeData(geno_reader)

  mypcair_r2 <- pcair(
    geno_data,
    kinobj    = mypcrel_1round_mat,
    divobj    = KINGmat,
    num.cores = n_cores,
    eigen.cnt = n_pcs
  )
  saveRDS(mypcair_r2, pcair_r2_rds)
  close(geno_data)
}
summary(mypcair_r2)

##plot the new PCs
for (i in seq(1, n_pcs - 1, by = 2)) {
  png(sprintf(pcpair_png_pattern_r2, i, i + 1),
      width = pcpair_plot_width, height = pcpair_plot_height, res = pcpair_plot_res)
  plot(mypcair_r2, vx = i, vy = i + 1)
  dev.off()
}
#############################################################
# Run the secound Round of correlations of genotypes with PCs
#############################################################

genofile <- snpgdsOpen(king_gds)

chr    <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
pos    <- read.gdsn(index.gdsn(genofile, "snp.position"))
snp_id <- read.gdsn(index.gdsn(genofile, "snp.id"))

cr2 <- snpgdsPCACorr(
  pcaobj     = mypcair_r2$vectors[, 1:n_pcs],
  gdsobj     = genofile,
  eig.which  = 1:n_pcs,
  num.thread = n_cores
)

snpgdsClose(genofile)

# align annotation to exactly the SNPs snpgdsPCACorr used
snp_idx_2r     <- match(cr2$snp.id, snp_id)
chr_used_2r    <- chr[snp_idx_2r]
pos_used_2r   <- pos[snp_idx_2r]
snp_id_used_2r <- snp_id[snp_idx_2r]

# chromosome x-axis: midpoint SNP index per chromosome, labeled chr1, chr2...
chr_levels_2r <- sort(unique(chr_used_2r))
chr_mids_2r   <- tapply(seq_along(chr_used_2r), chr_used_2r, median)
chr_labels_2r <- paste0("chr", names(chr_mids_2r))

# alternating colors by chromosome
# one color per chromosome
chr_colors        <- palette()[((as.integer(chr_levels_2r) - 1) %% length(palette())) + 1]
names(chr_colors) <- chr_levels_2r
col_vec           <- chr_colors[as.character(chr_used_2r)]
for (i in 1:n_pcs) {
  abs_corr_2r <- abs(cr2$snpcorr[i, ])
  # PNG with chromosome x-axis
  png(sprintf(snpcorr_png_pattern_r2, i),
      width = snpcorr_plot_width, height = snpcorr_plot_height, res = snpcorr_plot_res)
  plot(abs_corr_2r,
       col  = col_vec,
       pch  = 20,
       cex  = 0.3,
       ylim = c(0, 1),
       xaxt = "n",
       xlab = "Chromosome",
       ylab = paste0("PC", i, " Correlation"),
       main = paste0("PC-AiR 2nd round | PC", i, " — SNP correlation"))
  axis(1, at = chr_mids_2r, labels = chr_labels_2r, las = 2, cex.axis = 0.7)
  dev.off()

  # TSV sorted by absolute correlation descending
  df_2r <- data.frame(
    snp_id     = snp_id_used_2r,
    chromosome = chr_used_2r,
    position   = pos_used_2r,
    abs_corr_2r   = abs_corr_2r
  )
  df_2r <- df_2r[order(df_2r$abs_corr_2r, decreasing = TRUE), ]
  write.table(df_2r, sprintf(snpcorr_tsv_pattern_r2, i),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# ###############################################################         #############################################
# # Run the second round of PC-Relate, this is for the h2 project         still under development, not implemented yet here 
# ###############################################################         #############################################
# geno_reader <- GdsGenotypeReader(filename = pca_gds)
# geno_data   <- GenotypeData(geno_reader)
# geno_iter   <- GenotypeBlockIterator(geno_data, snpBlock = snp_block_size)
# mypcrel_r2 <- pcrelate(
#   geno_iter,
#   pcs          = mypcair_r2$vectors[, 1:n_pcs],
#   training.set = mypcair_r2$unrels,
#   BPPARAM      = BiocParallel::MulticoreParam(n_cores)
# )
# saveRDS(mypcrel_r2, pcrel_r2_rds)
# close(geno_data)



# ##################################################################################
# # Construct the format to then input to natora to identify the related individuals 
# ##################################################################################
print(paste0("reading the pcrel_r1_mat_rds file to then create the long format pairwise kinship table"))
mypcrel_r1 <- readRDS(pcrel_r1_mat_rds)
# n_samples x n_samples
dim(mypcrel_r1)    
# top-left corner                 
mypcrel_r1[1:5, 1:5]               

# Long-format pairwise kinship: ID1  ID2  kinship ####────────────
# - upper triangle only  -> each unique pair once, self-kinship (diagonal) dropped
# - all pairs kept, negative values left as-is
# - tab-separated, no header
# Output name: kinship_pairs_tsv in CONFIG

ids <- rownames(mypcrel_r1)
# dsyMatrix -> base symmetric matrix
M   <- as.matrix(mypcrel_r1)              
  # row/col indices, i < j
ut  <- which(upper.tri(M), arr.ind = TRUE)

pairs <- data.frame(
  ID1     = ids[ut[, "row"]],
  ID2     = ids[ut[, "col"]],
  kinship = M[ut],
  stringsAsFactors = FALSE
)

cat("pairs to write:", nrow(pairs), "\n") 
head(pairs)

write.table(pairs, kinship_pairs_tsv,
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# ##################################################################################
# # construct the covariate table with the pcs and the other covariates 
# ##################################################################################
#here we are using the round 2 of pcair since those carry the cleaner kinship estimates from PCrelate
print(paste0("reading the pcair_r2_rds file to then create the covariate table with the pcs and the other covariates")) 
pcair_r2 <- readRDS(pcair_r2_rds)
pcs <- as.data.frame(pcair_r2$vectors)
str(pcair_r2)

colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
pcs <- cbind(sample.id = rownames(pcs), pcs)
write.table(pcs, pcair_r2_covar_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

pcs <- read.table(pcair_r2_covar_tsv, sep = "\t", header = TRUE)
head(pcs)

psam <- read.table(covar_file,
                   sep = "\t", header = TRUE, check.names = FALSE)
colnames(psam)[1] <- "IID"

head(psam)

merged <- merge(psam, pcs, by.x = "IID", by.y = "sample.id", all.x = TRUE)

print(paste0("header merged df"))
head(merged)
colSums(is.na(merged))
colMeans(is.na(merged))

write.table(merged, pcair_r2_covar_merged_tsv, sep ='\t', row.names = FALSE, quote = FALSE)

# Close the log (see LOGGING). Messages first, then stdout, in reverse order of
# how they were opened. If the run dies earlier, R closes them when it exits.
message("[log] ", format(Sys.time()), " — done.")
sink(type = "message")
sink()
close(log_con)

```
</details>


Here an example on how to run it
```
(genesis) [duarte@node1 working_directory_h2]$ cd pca_and_such/
(genesis) [duarte@node1 pca_and_such]$ pwd
/home/duarte/working_directory_h2/pca_and_such
(genesis) [duarte@node1 pca_and_such]$ nano pcs_pipeline_h2_project.r 
(genesis) [duarte@node1 pca_and_such]$ ls
covar.txt  pcs_pipeline_h2_project.r
(genesis) [duarte@node1 pca_and_such]$ # then run the script resolving to the relative paths of the project structure
(genesis) [duarte@node1 pca_and_such]$ Rscript pcs_pipeline_h2_project.r \
  --input_plink_file ../genetic_data/genotyped/allchr_typed \
  --input_format pfile \
  --covar_file covar.txt \
  --removeHighLDregions yes \
  --n_cores 4 
```
The script has some default settings, we encourage you to leave them as they are for the first runs. If during the assessment of your results we see advisable to change some of the default parameters we can see which options to tweak. 

Also make sure to input how many cores your machine has available for this analysis. For some context, I was working here with approximately 7.1k samples and a near to 1M variants from the genotyped data and the script took 1.1 hours to finalize 

After the run is completed, get inside the output folder and inspect some of the diagnostic plots. 

The script will generate PCA plots and genome-wide correlations with each PCA across two rounds of PC-Air. This particular R package recommends two rounds of PC-Air, the first one using a GRM derived from KING's robust method and the second one using the GRM derived from PC-Relate (but since PC-Relate requires the input from PC-Air, that is why a first round of KING is performed), more details on how the logic flows can be found [here] (https://bioconductor.org/packages/devel/bioc/vignettes/GENESIS/inst/doc/pcair.html).

For the purposes of this analysis, we mostly care about the results form the second round of PC-Air. So, it is highly advisable that you inspect the first 10 PCs visually and make sure that substantial genetic variation is captured across the axes of each PC, and inspect the genome-wide correlation plots to make sure that each PC has an overall even contribution across the genome. "Problematic PCs" could be the ones that are driven by specific regions of the genome or that potentially don't carry that much "variability" to split samples in your cohort. More details on how to interpret the genome-wide correlation plots could be found [here] (https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011242). Don’t hesitate to reach out if you have questions about this.

We will also appreciate if you could share the plots generated here with us so we can help to assess them as well. (Note, there is no need to collect them by hand: Step 13 at the end of the tutorial provides a single command that finds all the "shareables" across the whole project and puts them in one folder ready to zip and send). 


### Step 4.b remove related individuals 

This is important since for cov-LDSc, the reference panel and the summary statistics from your gwas must contain the same variant IDs, and in my personal experience this is an easy way to make sure about it (working on the assumption that your data was previously aligned to the human reference fasta file before imputation, a common step during pre-imputation QC). 



Using the GRM derived from the PC-AiR and PC-Relate runs, we can more accurately estimate the amount of related individuals. 

In order to remove the least amount of samples as possible, we are going to use the logic behind network-based relatedness-pruning, using the tool named NAToRA (https://spj.science.org/doi/10.1016/j.csbj.2022.04.009). 

Hence, now you can turn off your GENESIS package and run the tool:

```
(genesis) [duarte@node1 outFolder_pca_andSuch]$ conda deactivate 
(base) [duarte@node1 outFolder_pca_andSuch]$ git clone https://github.com/ldgh/NAToRA_Public
Cloning into 'NAToRA_Public'...
remote: Enumerating objects: 135, done.
remote: Counting objects: 100% (9/9), done.
remote: Compressing objects: 100% (7/7), done.
remote: Total 135 (delta 3), reused 6 (delta 2), pack-reused 126 (from 1)
Receiving objects: 100% (135/135), 126.52 MiB | 64.65 MiB/s, done.
Resolving deltas: 100% (47/47), done.
(base) [duarte@node1 outFolder_pca_andSuch]$ # then put to run the tool, sometimes it fails because you need to have installed networkx, you can solve it with pip install
(base) [duarte@node1 outFolder_pca_andSuch]$ python NAToRA_Public/NAToRA_Public.py -i pcrelate_r1_kinship_pairs.tsv -o NAToRA_output_pcrel -c 0.0884

```

### Step 4.c format the covariate file for downstream analysis

Different tools to benchmark have different requirements on the type of covariate file format that the need, run this R script to subset the covariate file to the unrelated pairs identified through NAToRA, and create the different variations needed for the covariate file

<details>
<summary><strong>View script: <code>make_covariate_files.r</code></strong></summary>

<br>

```r
#!/usr/bin/env Rscript

# Build the covariate/phenotype files for the unrelated subset.
#
# Usage:
#   Rscript make_covariate_files.R [covariates.tsv] [toRemove.txt] [outdir] [prefix]
#
# Defaults:
#   covariates.tsv = pcair_r2_covariates_merged.tsv
#   toRemove.txt   = NAToRA_output_pcrel_toRemove.txt
#   outdir         = .
#   prefix         = covariate_file_no_related_pairs
#
# Everything printed to the screen is also written to <outdir>/<prefix>.log

args <- commandArgs(trailingOnly = TRUE)

covar_file  <- if (length(args) >= 1) args[1] else "./pcair_r2_covariates_merged.tsv"
remove_file <- if (length(args) >= 2) args[2] else "./NAToRA_output_pcrel_toRemove.txt"
outdir      <- if (length(args) >= 3) args[3] else "."
prefix      <- if (length(args) >= 4) args[4] else "covariate_file_no_related_pairs"

n_pcs <- 10

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ------------------------------------------------------------------- logging
# say() prints to the screen and appends the same line to the log file.
log_file <- file.path(outdir, paste0(prefix, ".log"))
cat("", file = log_file)   # start a fresh log on every run

say <- function(...) {
  msg <- paste0(..., "\n")
  cat(msg)
  cat(msg, file = log_file, append = TRUE)
}

# Errors are logged too, then re-raised so the exit status stays non-zero.
die <- function(...) {
  say("ERROR: ", ...)
  stop(..., call. = FALSE)
}

say("=== make_covariate_files.R ===============================")
say("Run started : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("Working dir : ", getwd())

if (!file.exists(covar_file))  die("Covariate file not found: ", covar_file)
if (!file.exists(remove_file)) die("Removal list not found: ", remove_file)

# ---------------------------------------------------------------- read input
# Everything is read as character so IDs and numeric values are written back
# exactly as they came in (no rounding, no scientific notation).
covar <- read.table(covar_file, header = TRUE, sep = "\t",
                    colClasses = "character", check.names = FALSE,
                    comment.char = "", quote = "")

to_remove <- readLines(remove_file)
to_remove <- trimws(to_remove)
to_remove <- to_remove[nzchar(to_remove)]

# --------------------------------------------------------- check the columns
pc_names <- paste0("PC", seq_len(n_pcs))
needed   <- c("IID", "SEX", "AGE", "DISEASE", pc_names)
missing  <- setdiff(needed, names(covar))
if (length(missing) > 0) {
  die("Missing column(s) in ", covar_file, ": ", paste(missing, collapse = ", "))
}

# --------------------------------------------------------------- filter rows
n_before <- nrow(covar)
keep     <- !(covar$IID %in% to_remove)
covar    <- covar[keep, , drop = FALSE]
n_after  <- nrow(covar)

if (n_after == 0) die("No samples left after removing the related individuals.")

FID <- covar$IID
IID <- covar$IID

# ------------------------------------------------------------- write outputs
out <- function(suffix) file.path(outdir, paste0(prefix, suffix))

written <- character(0)   # collected for the log

write_tsv <- function(df, path, col.names = FALSE) {
  write.table(df, path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = col.names)
  written <<- c(written, path)
}

# 1. FID IID SEX  (categorical covariate, no header)
write_tsv(data.frame(FID, IID, covar$SEX), out("_covar.tsv"))

# 2. FID IID  (no header)
write_tsv(data.frame(FID, IID), out("_IIDs.txt"))

# 3. FID IID PC1..PC10  (no header)
write_tsv(data.frame(FID, IID, covar[, pc_names]), out("_PCs.tsv"))

# 4. FID IID DISEASE  (no header)
write_tsv(data.frame(FID, IID, covar$DISEASE), out("_pheno.tsv"))

# 5. FID IID AGE PC1..PC10  (quantitative covariates, no header)
write_tsv(data.frame(FID, IID, covar$AGE, covar[, pc_names]), out("_qcovar.tsv"))

# 6. Full covariate table with a proper header (IID column becomes FID + IID)
full <- data.frame(FID, covar, check.names = FALSE)
write_tsv(full, out(".tsv"), col.names = TRUE)

# 7. Single-column list of the unrelated IIDs
iid_file <- file.path(outdir, "unrelated_IIDs.txt")
writeLines(IID, iid_file)
written <- c(written, iid_file)

# -------------------------------------------------------------- sample counts
disease  <- covar$DISEASE
n_cases  <- sum(disease == "2", na.rm = TRUE)
n_ctrls  <- sum(disease == "1", na.rm = TRUE)
n_other  <- n_after - n_cases - n_ctrls
prop     <- n_cases / (n_cases + n_ctrls)

say("")
say("--- Summary ---------------------------------------------")
say("Covariate file          : ", covar_file)
say("Removal list            : ", remove_file,
    " (", length(unique(to_remove)), " unique IDs)")
say("Samples before removal  : ", n_before)
say("Samples removed         : ", n_before - n_after)
say("Total sample size       : ", n_after)
say("Cases    (DISEASE = 2)  : ", n_cases)
say("Controls (DISEASE = 1)  : ", n_ctrls)
if (n_other > 0) {
  say("Other/missing DISEASE   : ", n_other, " (excluded from the proportion)")
}
say("Disease proportion      : ", sprintf("%.6f", prop),
    sprintf(" (%.2f%%)", 100 * prop))
say("")
say("--- Files written ---------------------------------------")
for (f in written) say("  ", f, "  (", nrow(covar), " rows)")
say("")
say("Run finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("Log file    : ", log_file)
say("=========================================================")

```
</details>

Then hit run:

```
(base) [duarte@node1 outFolder_pca_andSuch]$ nano make_covariate_files.r
(base) [duarte@node1 outFolder_pca_andSuch]$ Rscript  make_covariate_files.r
```


## Step 5. construction of the reference panel and computation of LD scores
### Step 5.a. prepare the environment and required files
Before starting the construction of the reference panel, we must download the versions of the programs we need, specifically cov-LDSC and LDSC conda environment (since these programs were written in python2). 


So lets go to our programs folder again:


<details>
<summary>terminal example:</summary>

```console
cd programs/
(base) [duarte@node1 programs]$ pwd
/home/duarte/working_directory_h2/programs
(base) [duarte@node1 programs]$ git clone https://github.com/immunogenomics/cov-ldsc.git
Cloning into 'cov-ldsc'...
remote: Enumerating objects: 1190, done.
remote: Counting objects: 100% (78/78), done.
remote: Compressing objects: 100% (33/33), done.
remote: Total 1190 (delta 24), reused 73 (delta 23), pack-reused 1112 (from 1)
Receiving objects: 100% (1190/1190), 49.56 MiB | 55.77 MiB/s, done.
Resolving deltas: 100% (491/491), done.
Updating files: 100% (1088/1088), done.
(base) [duarte@node1 programs]$ #then get the github for the classic ldsc and create the conda environment
(base) [duarte@node1 programs]$ #the conda environments suits both, cov-ldsc to create the reference panel with the LD scores
(base) [duarte@node1 programs]$ #and then to run the ldsc regression to calculate h2
(base) [duarte@node1 programs]$ git clone https://github.com/bulik/ldsc.git
Cloning into 'ldsc'...
remote: Enumerating objects: 7658, done.
remote: Counting objects: 100% (310/310), done.
remote: Compressing objects: 100% (38/38), done.
remote: Total 7658 (delta 286), reused 273 (delta 272), pack-reused 7348 (from 1)
Receiving objects: 100% (7658/7658), 56.44 MiB | 1.32 MiB/s, done.
Resolving deltas: 100% (2724/2724), done.
Updating files: 100% (1093/1093), done.
(base) [duarte@node1 programs]$ ls
cov-ldsc  intel-simplified-software-license.txt  ldsc  LICENSE  plink  plink2  plink2_linux_x86_64_20260808.zip  plink_linux_x86_64_20250819.zip  prettify  toy.map  toy.ped  vcf_subset
(base) [duarte@node1 programs]$ cd ldsc/
(base) [duarte@node1 ldsc]$ #create the conda environment with the requirements the authors of LDSC recommend 
(base) [duarte@node1 ldsc]$ #creating the environment takes a couple minutes
(base) [duarte@node1 ldsc]$ conda env create --file environment.yml
(base) [duarte@node1 ldsc]$ conda activate ldsc 
(ldsc) [duarte@node1 ldsc]$ python2 ldsc.py 
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \

*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \

Error: no analysis selected.
ldsc.py -h describes options.
Total time elapsed: 0.0s

(ldsc) [duarte@node1 ldsc]$ #it worked! forget the error message since we are not running the tool itself right now, you can also try: python2 ldsc.py -h
(ldsc) [duarte@node1 ldsc]$ #to check the flags available, just another way of showing that the conda environment and python versions are working appropriately.


```
</details>

Reference panels for underrepresented populations are challenging to acquire (hence the name "underrepresented"). Hence, a convenient way to overcome this problem is use each cohort's own genetic data to construct an in house reference panel. 

However, this comes with its own challenges. Especially since the majority of our cohort's genetic data is genotype array data and the imputed data on top. So we must ensure the biggest coverage possible genome-wide, while preserving high quality variants. 

The solution we propose to this particular problem consists on generating the list of imputed variants that carry an INFO R2 bigger or equal than 0.8, complemented by the union of genotyped variants and high-confidence variants identified through the hap map project.
Underrepresented populations may suffer more frequently from low quality imputed variants, so in order to increase the coverage of the reference panel, we extracted the well imputed variants, and complemented them with the union of the genotyped variants and the hapmap3 variants, to then compute
the adjusted LD scores.

Additionally, we must also update the centimorgans map of our reference panel, since it is fundamental to compute the LD scores adjusted for global ancestry estimated derived from PCs (with a window of 20 cM as recommended by the authors of cov-LDSC).

Alongside this repository, we have attached the files named hpm3snplist.bed and genetic_maps.b38_shapeit4.tar.gz, make sure to download them and put them inside the directory intended to construct the reference panel.

The genetic map was obtained from Shapeit 4 documentation. Likewise, the list of snps from the hap map project was obtained from GWASLab documentation. 

Example of command line:
<details>
<summary>terminal example:</summary>

```console
(base) [duarte@node1 working_directory_h2]$ pwd
/home/duarte/working_directory_h2
(base) [duarte@node1 working_directory_h2]$ mkdir reference_panel
(base) [duarte@node1 working_directory_h2]$ cd reference_panel/
(base) [duarte@node1 reference_panel]$ wget https://github.com/Jduarte-z/NarrowSense_h2_in_URPs_GP2/raw/refs/heads/main/genetic_maps.b38_shapeit4.tar.gz
Resolving github.com (github.com)... 140.82.113.4
Connecting to github.com (github.com)|140.82.113.4|:443... connected.
HTTP request sent, awaiting response... 302 Found
Location: https://raw.githubusercontent.com/Jduarte-z/NarrowSense_h2_in_URPs_GP2/refs/heads/main/genetic_maps.b38_shapeit4.tar.gz [following]
Resolving raw.githubusercontent.com (raw.githubusercontent.com)... 185.199.108.133, 185.199.109.133, 185.199.110.133, ...
Connecting to raw.githubusercontent.com (raw.githubusercontent.com)|185.199.108.133|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 23440558 (22M) [application/octet-stream]
Saving to: ‘genetic_maps.b38_shapeit4.tar.gz’

genetic_maps.b38_shapeit4.tar.gz                                                100%[====================================================================================================================================================================================================>]  22.35M  --.-KB/s    in 0.1s    

‘genetic_maps.b38_shapeit4.tar.gz’ saved [23440558/23440558]

(base) [duarte@node1 reference_panel]$ tar -xvf genetic_maps.b38_shapeit4.tar.gz
chr10.b38.gmap.gz
chr11.b38.gmap.gz
chr12.b38.gmap.gz
chr13.b38.gmap.gz
chr14.b38.gmap.gz
chr15.b38.gmap.gz
chr16.b38.gmap.gz
chr17.b38.gmap.gz
chr18.b38.gmap.gz
chr19.b38.gmap.gz
chr1.b38.gmap.gz
chr20.b38.gmap.gz
chr21.b38.gmap.gz
chr22.b38.gmap.gz
chr2.b38.gmap.gz
chr3.b38.gmap.gz
chr4.b38.gmap.gz
chr5.b38.gmap.gz
chr6.b38.gmap.gz
chr7.b38.gmap.gz
chr8.b38.gmap.gz
chr9.b38.gmap.gz
chrX.b38.gmap.gz
chrX_par1.b38.gmap.gz
chrX_par2.b38.gmap.gz
(base) [duarte@node1 reference_panel]$ gunzip *
(base) [duarte@node1 reference_panel]$ wget https://github.com/Jduarte-z/NarrowSense_h2_in_URPs_GP2/raw/refs/heads/main/hpm3snplist.bed
Resolving github.com (github.com)... 140.82.113.4
Connecting to github.com (github.com)|140.82.113.4|:443... connected.
HTTP request sent, awaiting response... 302 Found
Location: https://raw.githubusercontent.com/Jduarte-z/NarrowSense_h2_in_URPs_GP2/refs/heads/main/hpm3snplist.bed [following]
Resolving raw.githubusercontent.com (raw.githubusercontent.com)... 185.199.109.133, 185.199.110.133, 185.199.111.133, ...
Connecting to raw.githubusercontent.com (raw.githubusercontent.com)|185.199.109.133|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 25204459 (24M) [text/plain]
Saving to: ‘hpm3snplist.bed’

hpm3snplist.bed                                                                 100%[====================================================================================================================================================================================================>]  24.04M  --.-KB/s    in 0.1s    

‘hpm3snplist.bed’ saved [25204459/25204459]

(base) [duarte@node1 reference_panel]$ ls
chr10.b38.gmap  chr12.b38.gmap  chr14.b38.gmap  chr16.b38.gmap  chr18.b38.gmap  chr1.b38.gmap   chr21.b38.gmap  chr2.b38.gmap  chr4.b38.gmap  chr6.b38.gmap  chr8.b38.gmap  chrX.b38.gmap       chrX_par2.b38.gmap             hpm3snplist.bed
chr11.b38.gmap  chr13.b38.gmap  chr15.b38.gmap  chr17.b38.gmap  chr19.b38.gmap  chr20.b38.gmap  chr22.b38.gmap  chr3.b38.gmap  chr5.b38.gmap  chr7.b38.gmap  chr9.b38.gmap  chrX_par1.b38.gmap  genetic_maps.b38_shapeit4.tar

```
</details>
 

### Step 5.b compute reference panel and adjusted LD scores 

Then, considering the availability of a cluster, here is the draft for submitting 22 jobs to process the imputed data in parallel. Also, below we provide a standalone bash script to do it interactively in the terminal as well. 

Whichever version you run, the logic per chromosome is the same four steps:

1. **Identify the well-imputed variants.** Read the imputed VCF, keep only variants with `INFO R2 >= 0.8`, drop the high-LD regions, drop duplicated positions, and write out the surviving variant IDs as a snplist.
2. **Build the combined variant list.** Convert that snplist and your genotyped variant list into position ranges, and concatenate them with the HapMap3 list. This is the union described above: well-imputed ∪ genotyped ∪ HapMap3.
3. **Extract the panel.** Go back to the imputed VCF and pull out those variants, for the unrelated samples only, into PLINK1 binary format (cov-LDSC reads `bed/bim/fam`, not `pgen`). Then apply the MAF > 0.01 filter and write the centimorgan positions into the `.bim` using the SHAPEIT4 genetic map — this is what makes the 20 cM window meaningful.
4. **Compute the covariate-adjusted LD scores** with cov-LDSC, passing the same PCs that will go into the GWAS.

Before running it, notice that every input this step needs was already produced earlier in the tutorial. It is worth checking that all of them are where the script expects them:

| input | where it comes from | path used below |
|---|---|---|
| imputed VCFs | Step 1 symlink | `../genetic_data/imputed/chr${i}.dose.vcf.gz` |
| unrelated sample list (`FID IID`, no header) | Step 4.c | `../pca_and_such/outFolder_pca_andSuch/covariate_file_no_related_pairs_IIDs.txt` |
| PC covariates (`FID IID PC1..PC10`, no header) | Step 4.c | `../pca_and_such/outFolder_pca_andSuch/covariate_file_no_related_pairs_PCs.tsv` |
| Population proportion of PD in your cohort | Step 4.c | `../pca_and_such/outFolder_pca_andSuch/covariate_file_no_related_pairs.log` |
| high-LD regions | Step 4.a (only written when `--removeHighLDregions yes`) | `../pca_and_such/outFolder_pca_andSuch/highLD_regions_grindeLab_hg38.tsv` |
| genotyped variant IDs | Step 4.a | `../pca_and_such/outFolder_pca_andSuch/input_setVarIDs.snplist` |
| HapMap3 positions | Step 5.a download | `./hpm3snplist.bed` |
| genetic maps | Step 5.a download | `./chr${i}.b38.gmap` |

**One thing to keep in mind:** the covariate file handed to `--cov` here must contain exactly the same principal components that you will later include in the GWAS regression model. This is the core requirement of cov-LDSC, and it is the reason we built a single covariate file back in Step 4.c instead of letting each tool make its own.

As a python script to submit 22 different jobs

<details>
<summary>script: <code>refPanel_covLDscores.py</code> — writes one SLURM .pbs per chromosome and submits all 22</summary>

```python
# Construction of the in-sample reference panel and computation of the
# covariate-adjusted LD scores.
# Run it inside /working_directory_h2/reference_panel

import os

# ---- cluster resources: adjust to what your scheduler offers -----------------
threads   = 4
mem       = "75G"      # the plink2 steps read whole imputed vcfs, they are memory hungry
mem_mb    = 75000      # the same number in MB, handed to plink's --memory
email     = "your_mail@email.org"
partition = "bigmem"

# ---- shared inputs: every one of these was produced in a previous step -------
# unrelated samples that survived NAToRA, in FID IID format
keepList_unrelated = "../pca_and_such/outFolder_pca_andSuch/covariate_file_no_related_pairs_IIDs.txt"
# FID IID PC1..PC10, no header. These MUST be the same PCs used in the GWAS
covariateFile      = "../pca_and_such/outFolder_pca_andSuch/covariate_file_no_related_pairs_PCs.tsv"
# the high LD regions list, written by the Rscript of the PCA pipeline
HLA_regions        = "../pca_and_such/outFolder_pca_andSuch/highLD_regions_grindeLab_hg38.tsv"
# the variant IDs of your genotyped data, also written by the Rscript of the PCA pipeline
genotyped_varsList = "../pca_and_such/outFolder_pca_andSuch/input_setVarIDs.snplist"
# downloaded in step 5.a. Originally from gwaslab:
# https://github.com/Cloufield/gwaslab/tree/main/src/gwaslab/data/hapmap3_SNPs
# here we use one that we pre-processed to match our pipeline requirements
hapMap3_bedList    = "./hpm3snplist.bed"

# ---- programs, as relative paths from ./reference_panel ---------------------
plink2  = "../programs/plink2"
plink1  = "../programs/plink"
covldsc = "../programs/cov-ldsc/ldsc.py"

for i in range(1, 23):
    runFolder               = f"./chr{i}"
    vcfFileIn               = f"../genetic_data/imputed/chr{i}.dose.vcf.gz"
    geneticMap              = f"./chr{i}.b38.gmap"
    plink1FileOut           = f"{runFolder}/chr{i}.plink1.extracted"
    plink1_cmUpdatedFileOut = f"{runFolder}/chr{i}.plink1.extracted_cmUpdated"

    os.makedirs(f"{runFolder}/logs", exist_ok=True)
    os.makedirs(f"{runFolder}/out",  exist_ok=True)

    fileSbatch = open(f"{runFolder}/logs/{i}refPanel_covLDscores.pbs", "w")

    fileSbatch.write(f"#!/bin/sh\n")
    fileSbatch.write(f"#SBATCH --mail-type=END,FAIL\n")
    fileSbatch.write(f"#SBATCH --mail-user={email}\n")
    fileSbatch.write(f"#SBATCH --ntasks=1\n")
    fileSbatch.write(f"#SBATCH --cpus-per-task={threads}\n")
    fileSbatch.write(f"#SBATCH --mem={mem}\n")
    fileSbatch.write(f"#SBATCH --partition={partition}\n")
    fileSbatch.write(f"#SBATCH --job-name={i}refPanel_covLDscores\n")
    fileSbatch.write(f"#SBATCH -o {runFolder}/logs/{i}refPanel_covLDscores.out\n")
    fileSbatch.write(f"#SBATCH -e {runFolder}/logs/{i}refPanel_covLDscores.err\n\n")

    fileSbatch.write(
        # 1) well-imputed variants: R2 >= 0.8, outside high-LD regions, no duplicates
        f"{plink2} --vcf {vcfFileIn} "
        f"--double-id "
        f"--threads {threads} "
        f"--exclude range {HLA_regions} "
        f"--extract-if-info \"R2>=0.8\" "
        f"--write-snplist "
        f"--rm-dup exclude-all list "
        f"--set-all-var-ids @:#:\\$r:\\$a "
        f"--new-id-max-allele-len 100 "
        f"--memory {mem_mb} "
        f"--out {runFolder}/chr{i}.step1 \n\n"

        # 2) combined variant list = hapmap3 + well imputed + genotyped
        #    the awk turns chr:pos:ref:alt IDs into 'chrom start end' ranges
        f"awk -F':' '{{print $1, $2, $2}}' {runFolder}/chr{i}.step1.snplist > {runFolder}/wellImputed_bed_varList.txt \n"
        f"awk -F':' '{{print $1, $2, $2}}' {genotyped_varsList} > {runFolder}/genotyped_bed_varList.txt \n"
        f"cat {hapMap3_bedList} {runFolder}/wellImputed_bed_varList.txt {runFolder}/genotyped_bed_varList.txt > {runFolder}/combined_bed_varList.txt \n\n"

        # 3a) extract those variants, unrelated samples only, into plink1 format
        f"{plink2} --vcf {vcfFileIn} "
        f"--double-id "
        f"--threads {threads} "
        f"--extract bed1 {runFolder}/combined_bed_varList.txt "
        f"--exclude range {HLA_regions} "
        f"--keep {keepList_unrelated} "
        f"--keep-allele-order "
        f"--set-all-var-ids @:#:\\$r:\\$a "
        f"--new-id-max-allele-len 1000 "
        f"--rm-dup exclude-all list "
        f"--memory {mem_mb} "
        f"--make-bed --out {plink1FileOut} \n\n"

        # 3b) MAF filter and write the centimorgan positions into the .bim,
        #     which is what makes the 20 cM window meaningful
        f"{plink1} --bfile {plink1FileOut} "
        f"--threads {threads} "
        f"--maf 0.01 "
        f"--cm-map {geneticMap} {i} "
        f"--keep-allele-order "
        f"--memory {mem_mb} "
        f"--make-bed --out {plink1_cmUpdatedFileOut} \n\n"

        # 4) covariate-adjusted LD scores, 20 cM window
        #    cov-LDSC is python2, so it runs inside the ldsc conda environment
        #    if 'conda' is not on the PATH of your compute nodes, replace the line
        #    below with the explicit path to your miniconda3/etc/profile.d/conda.sh
        f"source \"$(conda info --base)/etc/profile.d/conda.sh\"\n"
        f"conda activate ldsc \n"
        f"python2 {covldsc} "
        f"--ld-wind-cm 20.0 "
        f"--cov {covariateFile} "
        f"--bfile {plink1_cmUpdatedFileOut} "
        f"--out {runFolder}/out/covldsc_chr{i} \n"
    )

    fileSbatch.close()
    os.system(f"sbatch {runFolder}/logs/{i}refPanel_covLDscores.pbs")
```

</details>

For running it do:

<details>
<summary>terminal example: saving the python script with nano and submitting the 22 jobs</summary>

```console
(base) [duarte@node1 reference_panel]$ pwd
/home/duarte/working_directory_h2/reference_panel
(base) [duarte@node1 reference_panel]$ # do a nano with the name of the script, hit enter, then ctrl + X, and then hit "y" to save it
(base) [duarte@node1 reference_panel]$ nano refPanel_covLDscores.py
(base) [duarte@node1 reference_panel]$ python refPanel_covLDscores.py
Submitted batch job 6112401
Submitted batch job 6112402
Submitted batch job 6112403
Submitted batch job 6112404
Submitted batch job 6112405
Submitted batch job 6112406
Submitted batch job 6112407
Submitted batch job 6112408
Submitted batch job 6112409
Submitted batch job 6112410
Submitted batch job 6112411
Submitted batch job 6112412
Submitted batch job 6112413
Submitted batch job 6112414
Submitted batch job 6112415
Submitted batch job 6112416
Submitted batch job 6112417
Submitted batch job 6112418
Submitted batch job 6112419
Submitted batch job 6112420
Submitted batch job 6112421
Submitted batch job 6112422

```

</details>

If you prefer to work in series through your terminal here is the equivalent for loop, it might take considerably more time to run, but it will do the job as well.

<details>
<summary>script: <code>refPanel_covLDscores.sh</code> — same reference panel and LD scores as a sequential for loop, no scheduler needed</summary>

```bash
#!/bin/bash
# Sequential version of refPanel_covLDscores.py
# Runs chr1 .. chr22 one after the other in the terminal (no scheduler).
# Run it inside /working_directory_h2/reference_panel

# note: -u is deliberately left out, conda's activation script trips on unset variables
set -eo pipefail

# ---- resources: adjust to what your machine has ------------------------------
threads=4
mem_mb=75000

# ---- shared inputs: every one of these was produced in a previous step -------
keepList_unrelated="../pca_and_such/outFolder_pca_andSuch/covariate_file_no_related_pairs_IIDs.txt"   # FID IID
covariateFile="../pca_and_such/outFolder_pca_andSuch/covariate_file_no_related_pairs_PCs.tsv"         # FID IID PC1..PC10
HLA_regions="../pca_and_such/outFolder_pca_andSuch/highLD_regions_grindeLab_hg38.tsv"
genotyped_varsList="../pca_and_such/outFolder_pca_andSuch/input_setVarIDs.snplist"
hapMap3_bedList="./hpm3snplist.bed"

# ---- programs, as relative paths from ./reference_panel ---------------------
plink2="../programs/plink2"
plink1="../programs/plink"
covldsc="../programs/cov-ldsc/ldsc.py"

# ---- fail early if anything is missing, before burning hours on chr1 ---------
for f in "${keepList_unrelated}" "${covariateFile}" "${HLA_regions}" \
         "${genotyped_varsList}" "${hapMap3_bedList}" "${plink2}" "${plink1}" "${covldsc}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing or empty: ${f}" >&2
        exit 1
    fi
done

# ---- conda env for cov-LDSC, activated once ---------------------------------
# if 'conda' is not on your PATH, replace this with the explicit path to your
# miniconda3/etc/profile.d/conda.sh
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ldsc

for i in $(seq 1 22); do

    echo "=================== chr${i} ==================="

    runFolder="./chr${i}"
    vcfFileIn="../genetic_data/imputed/chr${i}.dose.vcf.gz"
    geneticMap="./chr${i}.b38.gmap"
    plink1FileOut="${runFolder}/chr${i}.plink1.extracted"
    plink1_cmUpdatedFileOut="${runFolder}/chr${i}.plink1.extracted_cmUpdated"

    mkdir -p "${runFolder}/logs" "${runFolder}/out"

    # 1) well-imputed variants: R2 >= 0.8, outside high-LD regions, no duplicates
    "${plink2}" --vcf "${vcfFileIn}" \
        --double-id \
        --threads "${threads}" \
        --exclude range "${HLA_regions}" \
        --extract-if-info "R2>=0.8" \
        --write-snplist \
        --rm-dup exclude-all list \
        --set-all-var-ids '@:#:$r:$a' \
        --new-id-max-allele-len 1000 \
        --memory "${mem_mb}" \
        --out "${runFolder}/chr${i}.step1"

    # 2) combined variant list = hapmap3 + well imputed + genotyped
    #    the awk turns chr:pos:ref:alt IDs into 'chrom start end' ranges
    awk -F':' '{print $1, $2, $2}' "${runFolder}/chr${i}.step1.snplist" > "${runFolder}/wellImputed_bed_varList.txt"
    awk -F':' '{print $1, $2, $2}' "${genotyped_varsList}"              > "${runFolder}/genotyped_bed_varList.txt"
    cat "${hapMap3_bedList}" \
        "${runFolder}/wellImputed_bed_varList.txt" \
        "${runFolder}/genotyped_bed_varList.txt"  > "${runFolder}/combined_bed_varList.txt"

    # 3a) extract those variants, unrelated samples only, into plink1 format
    "${plink2}" --vcf "${vcfFileIn}" \
        --double-id \
        --threads "${threads}" \
        --extract bed1 "${runFolder}/combined_bed_varList.txt" \
        --exclude range "${HLA_regions}" \
        --keep "${keepList_unrelated}" \
        --keep-allele-order \
        --set-all-var-ids '@:#:$r:$a' \
        --new-id-max-allele-len 1000 \
        --rm-dup exclude-all list \
        --memory "${mem_mb}" \
        --make-bed --out "${plink1FileOut}"

    # 3b) MAF filter and write the centimorgan positions into the .bim,
    #     which is what makes the 20 cM window meaningful
    "${plink1}" --bfile "${plink1FileOut}" \
        --threads "${threads}" \
        --maf 0.01 \
        --cm-map "${geneticMap}" "${i}" \
        --keep-allele-order \
        --memory "${mem_mb}" \
        --make-bed --out "${plink1_cmUpdatedFileOut}"

    # 4) covariate-adjusted LD scores, 20 cM window
    python2 "${covldsc}" \
        --ld-wind-cm 20.0 \
        --cov "${covariateFile}" \
        --bfile "${plink1_cmUpdatedFileOut}" \
        --out "${runFolder}/out/covldsc_chr${i}"

done

echo "All chromosomes done."
```

</details>

In the terminal do

<details>
<summary>terminal example: saving and running the bash version</summary>

```console
(base) [duarte@node1 reference_panel]$ nano refPanel_covLDscores.sh
(base) [duarte@node1 reference_panel]$ bash refPanel_covLDscores.sh
=================== chr1 ===================
(base) [duarte@node1 reference_panel]$ # if it doesn't run, try to make sure it is an executable
(base) [duarte@node1 reference_panel]$ chmod +x refPanel_covLDscores.sh #and then try again to bash it
```

</details>

After running either of the two options, you will end up with one folder per chromosome inside `reference_panel`, each one holding its own intermediate files and its LD scores.

<details>
<summary>terminal example: what the reference_panel directory looks like</summary>

```console
(base) [duarte@node1 reference_panel]$ ls
chr1   chr12  chr15  chr18  chr20  chr4  chr7  chr1.b38.gmap   chr12.b38.gmap  chr15.b38.gmap  chr18.b38.gmap  chr20.b38.gmap  chr4.b38.gmap  chr7.b38.gmap  hpm3snplist.bed
chr10  chr13  chr16  chr19  chr21  chr5  chr8  chr10.b38.gmap  chr13.b38.gmap  chr16.b38.gmap  chr19.b38.gmap  chr21.b38.gmap  chr5.b38.gmap  chr8.b38.gmap  refPanel_covLDscores.py
chr11  chr14  chr17  chr2   chr22  chr6  chr9  chr11.b38.gmap  chr14.b38.gmap  chr17.b38.gmap  chr2.b38.gmap   chr22.b38.gmap  chr6.b38.gmap  chr9.b38.gmap  genetic_maps.b38_shapeit4.tar
chr3                                                           chr3.b38.gmap

(base) [duarte@node1 reference_panel]$ # a closer look at one chromosome
(base) [duarte@node1 reference_panel]$ ls chr22/
chr22.plink1.extracted.bed             chr22.plink1.extracted_cmUpdated.bed   chr22.step1.log        genotyped_bed_varList.txt
chr22.plink1.extracted.bim             chr22.plink1.extracted_cmUpdated.bim   chr22.step1.rmdup.list  logs
chr22.plink1.extracted.fam             chr22.plink1.extracted_cmUpdated.fam   chr22.step1.snplist     out
chr22.plink1.extracted.log             chr22.plink1.extracted_cmUpdated.log   combined_bed_varList.txt
chr22.plink1.extracted.rmdup.list                                             wellImputed_bed_varList.txt

(base) [duarte@node1 reference_panel]$ # and the LD scores themselves
(base) [duarte@node1 reference_panel]$ ls chr22/out/
covldsc_chr22.l2.ldscore.gz  covldsc_chr22.l2.M  covldsc_chr22.l2.M_5_50  covldsc_chr22.log

(base) [duarte@node1 reference_panel]$ zcat chr22/out/covldsc_chr22.l2.ldscore.gz | head -3
CHR	SNP	BP	L2
22	22:16392862:G:T	16392862	19.482
22	22:16393163:T:C	16393163	19.563
```

</details>

Two quick sanity checks before moving on. First, confirm that all 22 chromosomes actually produced LD scores, since a single chromosome that ran out of memory is easy to miss:

<details>
<summary>terminal example: checking that the 22 chromosomes finished</summary>

```console
(base) [duarte@node1 reference_panel]$ ls chr*/out/*.l2.ldscore.gz | wc -l
22
(base) [duarte@node1 reference_panel]$ # total number of variants in the reference panel
(base) [duarte@node1 reference_panel]$ cat chr*/chr*.plink1.extracted_cmUpdated.bim | wc -l
8102072
```

</details>

Second, remember the point from the cov-LDSC section above: the whole reason for the 20 cM window is that in admixed cohorts the mean LD score keeps rising with window size until the PC adjustment is applied. If you want to verify that your cohort's LD scores have actually plateaued, you can re-run a single chromosome across a few window sizes (for example 1, 5, 10, 20 and 50 cM) and plot the mean $$L2$$ against the window. I don't believe this would be needed for our project but it is good to know just in case. 

And FYI, considering the test dataset that we used in LARGE-PD with 7.1K samples and more than 8M variants (post filtering). The longest chromosomes (1-5) took a maximum of 8 hours each, with the configurations for the computing cluster specified. Hence, this is one of the most computationally extensive steps


## Step 6 run a GWAS with the same universe of SNPs and samples used for the reference panel 

Recall the caveat from section 4: cov-LDSC regresses the chi-square statistics of a GWAS on LD scores that were adjusted for principal components, and the method only holds if **those are the same principal components that were included in the GWAS regression model**. 

The good news is that we have already done the work that makes this automatic. The `bed/bim/fam` filesets we built in Step 5.b are already restricted to the unrelated samples, already MAF-filtered, and already carry the final variant IDs. So instead of running the GWAS on the imputed data and then trying to reconcile the two variant sets afterwards, **we run the GWAS directly on the reference panel filesets themselves**. The universes then cannot drift apart, because they are literally the same files. 

Note: the authors of cov-LDSC argue that you can use a subset of your samples to compute the reference panel (between 1-5K).However, since we like to benchmark several other methods using this same universe of genetic variants, we are proposing on using the whole cohort in both the GWAS and reference panel. Because this full parsed genetic dataset will be feed into GRM based methods to compute them either way. 

Let us create a folder for this step, at the same level as `reference_panel`:

<details>
<summary>terminal example: creating the gwas folder</summary>

```console
(base) [duarte@node1 reference_panel]$ cd ..
(base) [duarte@node1 working_directory_h2]$ mkdir gwas
(base) [duarte@node1 working_directory_h2]$ cd gwas/
(base) [duarte@node1 gwas]$ pwd
/home/duarte/working_directory_h2/gwas
(base) [duarte@node1 working_directory_h2]$ ls ..
genetic_data  gwas  pca_and_such  programs  reference_panel
```

</details>

**A note on the covariates.** The model below uses `SEX`, `AGE` and `PC1..PC10`, all read from the full covariate table written in Step 4.c. In my own run I also carry a `PHASE` covariate, because LARGE-PD merges two recruitment phases that were genotyped separately, and batch effects of that kind must be adjusted for. Most cohorts will not have that column, so in the scripts below it is exposed as an `extra_covars` variable that is **empty by default**. If your cohort uses additional covariates remember that whatever you add here must also be present when the same covariate file is used in the GRM-based arms later on, otherwise the arms stop being comparable.

As a python script to submit 22 different jobs

<details>
<summary>script: <code>plinkGWAS.py</code> — writes one SLURM .pbs per chromosome and submits all 22</summary>

```python
# Quick GWAS on the same samples and variants that make up the reference panel.
# Run it inside /working_directory_h2/gwas

import os

# ---- cluster resources: adjust to what your scheduler offers -----------------
threads   = 8
email     = "your_mail@email.org"
partition = "defq"

# ---- covariates --------------------------------------------------------------
# The full covariate table written in step 4.c (has a header: FID IID SEX DISEASE AGE PC1..PC10)
covar = "../pca_and_such/outFolder_pca_andSuch/covariate_file_no_related_pairs.tsv"

# Any extra cohort-specific covariate, comma separated and WITHOUT a trailing comma.
# Leave it as "" if you have none. Example: "PHASE" for a cohort that merges two
# recruitment phases genotyped separately, or "BATCH,ARRAY".
extra_covars = ""

pcs         = ",".join(f"PC{n}" for n in range(1, 11))
covar_names = ",".join(x for x in ["SEX", "AGE", extra_covars, pcs] if x)

# ---- programs, as relative paths from ./gwas --------------------------------
plink2 = "../programs/plink2"

os.makedirs("./logs", exist_ok=True)
os.makedirs("./out",  exist_ok=True)

for i in range(1, 23):
    mem = "75G" if i in {1, 2, 3, 4, 5} else "50G"

    # the reference panel filesets from step 5.b, used here as-is
    plink_files = f"../reference_panel/chr{i}/chr{i}.plink1.extracted_cmUpdated"
    out         = f"./out/plinkGWAS_chr{i}.result"

    fileSbatch = open(f"./logs/{i}plinkGWAS.pbs", "w")

    fileSbatch.write(f"#!/bin/sh\n")
    fileSbatch.write(f"#SBATCH --mail-type=END,FAIL\n")
    fileSbatch.write(f"#SBATCH --mail-user={email}\n")
    fileSbatch.write(f"#SBATCH --ntasks=1\n")
    fileSbatch.write(f"#SBATCH --cpus-per-task={threads}\n")
    fileSbatch.write(f"#SBATCH --mem={mem}\n")
    fileSbatch.write(f"#SBATCH --partition={partition}\n")
    fileSbatch.write(f"#SBATCH --job-name={i}plinkGWAS\n")
    fileSbatch.write(f"#SBATCH -o ./logs/{i}plinkGWAS.out\n")
    fileSbatch.write(f"#SBATCH -e ./logs/{i}plinkGWAS.err\n\n")

    fileSbatch.write(
        f"{plink2} --bfile {plink_files} "
        # hide-covar keeps only the ADD line per variant, which is all we need
        # cols=+a1freq adds A1_FREQ, needed downstream by gwaslab
        f"--glm hide-covar cols=+a1freq omit-ref "
        f"--pheno {covar} --pheno-name DISEASE "
        f"--covar {covar} --covar-name {covar_names} "
        # puts every covariate on the same scale, which helps the model converge
        f"--covar-variance-standardize "
        f"--threads {threads} "
        f"--out {out} \n"
    )

    fileSbatch.close()
    os.system(f"sbatch ./logs/{i}plinkGWAS.pbs")
```

</details>

For running it do:

<details>
<summary>terminal example: saving the python script with nano and submitting the 22 jobs</summary>

```console
(base) [duarte@node1 gwas]$ nano plinkGWAS.py
(base) [duarte@node1 gwas]$ python plinkGWAS.py
Submitted batch job 6114512
Submitted batch job 6114513
Submitted batch job 6114514
Submitted batch job 6114515
Submitted batch job 6114516
Submitted batch job 6114517
Submitted batch job 6114518
Submitted batch job 6114519
Submitted batch job 6114520
Submitted batch job 6114521
Submitted batch job 6114522
Submitted batch job 6114523
Submitted batch job 6114524
Submitted batch job 6114525
Submitted batch job 6114526
Submitted batch job 6114527
Submitted batch job 6114528
Submitted batch job 6114529
Submitted batch job 6114530
Submitted batch job 6114531
Submitted batch job 6114532
Submitted batch job 6114533

```

</details>

If you prefer to work in series through your terminal here is the equivalent for loop. This step is much lighter than Step 5.b, so running it sequentially is perfectly reasonable.

<details>
<summary>script: <code>plinkGWAS.sh</code> — same GWAS as a sequential for loop, no scheduler needed</summary>

```bash
#!/bin/bash
# Sequential version of plinkGWAS.py
# Runs chr1 .. chr22 one after the other in the terminal (no scheduler).
# Run it inside /working_directory_h2/gwas

set -euo pipefail

threads=8

# ---- covariates -------------------------------------------------------------
covar="../pca_and_such/outFolder_pca_andSuch/covariate_file_no_related_pairs.tsv"

# Any extra cohort-specific covariate, comma separated and WITHOUT a trailing comma.
# Leave empty if you have none. Example: extra_covars="PHASE"
extra_covars=""

pcs="PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
if [[ -n "${extra_covars}" ]]; then
    covar_names="SEX,AGE,${extra_covars},${pcs}"
else
    covar_names="SEX,AGE,${pcs}"
fi

# ---- programs, as relative paths from ./gwas --------------------------------
plink2="../programs/plink2"

mkdir -p ./logs ./out

# ---- fail early if anything is missing --------------------------------------
if [[ ! -s "${covar}" ]]; then
    echo "ERROR: missing or empty covariate file: ${covar}" >&2
    exit 1
fi

for i in $(seq 1 22); do

    echo "=================== chr${i} ==================="

    plink_files="../reference_panel/chr${i}/chr${i}.plink1.extracted_cmUpdated"
    out="./out/plinkGWAS_chr${i}.result"

    if [[ ! -s "${plink_files}.bed" ]]; then
        echo "ERROR: reference panel fileset not found: ${plink_files}.bed" >&2
        exit 1
    fi

    "${plink2}" --bfile "${plink_files}" \
        --glm hide-covar cols=+a1freq omit-ref \
        --pheno "${covar}" --pheno-name DISEASE \
        --covar "${covar}" --covar-name "${covar_names}" \
        --covar-variance-standardize \
        --threads "${threads}" \
        --out "${out}" \
        > "./logs/${i}plinkGWAS.out" 2> "./logs/${i}plinkGWAS.err"

done

echo "All chromosomes done."
```

</details>

In the terminal do

<details>
<summary>terminal example: saving and running the bash version</summary>

```console
(base) [duarte@node1 gwas]$ nano plinkGWAS.sh
(base) [duarte@node1 gwas]$ bash plinkGWAS.sh
=================== chr1 ===================
(base) [duarte@node1 gwas]$ # if it doesn't run, try to make sure it is an executable
(base) [duarte@node1 gwas]$ chmod +x plinkGWAS.sh #and then try again to bash it
```

</details>

Either way you end up with one association file per chromosome inside `./out`. The `.glm.logistic.hybrid` extension tells you that PLINK2 used its Firth-fallback logistic regression, which is the sensible default for a case-control trait.

Note: something to consider here is that we are using plink1 files to run the regressions. Hence, we are working with hardcalls rather than dosages, that in the context of GWAS is not ideal. However, since association testing is not the main goal of this project (and we will not put spotlight on the GWAS results, we just need to generate the sumtats in the same conditions as the reference panel), we move forward with this caveat. 

<details>
<summary>terminal example: what the GWAS output looks like</summary>

```console
(base) [duarte@node1 gwas]$ ls out/ | head -6
plinkGWAS_chr10.result.DISEASE.glm.logistic.hybrid
plinkGWAS_chr10.result.log
plinkGWAS_chr11.result.DISEASE.glm.logistic.hybrid
plinkGWAS_chr11.result.log
plinkGWAS_chr12.result.DISEASE.glm.logistic.hybrid
plinkGWAS_chr12.result.log

(base) [duarte@node1 gwas]$ # sanity check: 22 association files
(base) [duarte@node1 gwas]$ ls out/*.glm.logistic.hybrid | wc -l
22

(base) [duarte@node1 gwas]$ head -3 out/plinkGWAS_chr22.result.DISEASE.glm.logistic.hybrid
#CHROM	POS	ID	REF	ALT	A1	OMITTED	A1_FREQ	FIRTH?	TEST	OBS_CT	OR	LOG(OR)_SE	Z_STAT	P	ERRCODE
22	10516675	chr22:10516675:G:A	G	A	A	G	0.231447	N	ADD	6412	1.02214	0.0487731	0.449276	0.653259	.
22	10516726	chr22:10516726:C:T	C	T	T	C	0.190228	N	ADD	6412	0.981135	0.0521884	-0.365184	0.714993	.
```

</details>

Take a moment to confirm that `OBS_CT` is close to the number of unrelated samples reported by `make_covariate_files.r` in Step 4.c. If it is much smaller, something went wrong with the ID matching between the covariate file and the `.fam` of the reference panel.

### Step 6.a parse the summary statistics with GWASLab

The 22 association files now have to become a single set of summary statistics in the format that LDSC reads. We use [GWASLab](https://cloufield.github.io/gwaslab/) for this, because it harmonises the alleles, rebuilds the variant IDs, sorts by coordinate and exports straight to the LDSC format in one pass, while also giving us the Manhattan and QQ plots we want to inspect anyway.

GWASLab is python3, whereas LDSC is python2, so it needs its own conda environment.

<details>
<summary>terminal example: creating the gwaslab environment</summary>

```console
(base) [duarte@node1 gwas]$ conda create -n gwaslab -c conda-forge python=3.10 -y
(base) [duarte@node1 gwas]$ conda activate gwaslab
(gwaslab) [duarte@node1 gwas]$ # we pin the version so that everybody in the project parses the sumstats identically
(gwaslab) [duarte@node1 gwas]$ pip install gwaslab==4.1.6

```

</details>

Then the parsing script itself. Note the `@` in the input path: GWASLab expands it to the chromosome number, so the 22 files are read and concatenated in a single call.

<details>
<summary>script: <code>parse_plink_gwas.py</code> — harmonises the 22 association files and exports them in LDSC format</summary>

```python
# Parse the per-chromosome PLINK2 association files into a single set of
# summary statistics in LDSC format.
# Run it inside /working_directory_h2/gwas, with the gwaslab environment active.

import os
import numpy as np
import gwaslab as gl

os.makedirs("./out/gwaslab_output", exist_ok=True)

threads = 4

# the @ is expanded by gwaslab to the chromosome number, so all 22 files are read at once
input_file  = "./out/plinkGWAS_chr@.result.DISEASE.glm.logistic.hybrid"
# rename this with your own cohort so the shared files are self-describing
output_file = "./out/gwaslab_output/gwas_sumstats_allchr"

ss = gl.Sumstats(
    input_file,
    snpid="ID",
    build="38",
    chrom="#CHROM",
    pos="POS",
    ea="A1",
    eaf="A1_FREQ",
    nea="OMITTED",
    p="P",
    n="OBS_CT",
    OR="OR",
    se="LOG(OR)_SE",
    sep="\t",
    z="Z_STAT",
    na_values=["NA", "."],
    verbose=True,
)

# ---- strict data clean -------------------------------------------------------
# every fix_* call with remove=True drops the records it cannot repair, which is
# what we want here: LDSC is unforgiving about malformed rows
ss.fix_chr(remove=True)
ss.fix_pos(remove=True)
ss.fix_allele(remove=True)
ss.fix_id(fixchrpos=False, fixid=True, fixsep=False, forcefixid=True, overwrite=True)
ss.normalize_allele(threads=threads)
ss.sort_coordinate()

# ---- diagnostic plots --------------------------------------------------------
ss.plot_mqq(
    mode='mqq',
    sig_level=1e-6,
    anno_sig_level=1e-6,
    anno="GENENAME",
    build="38",
    sig_line=True,
    additional_line=[5e-8],
    additional_line_color=['black'],
    font_family="DejaVu Sans",
    fontsize=11,
    anno_fontsize=12,
    colors=["#000000", "#ABABAB"],
    save="./out/gwaslab_output/gwas_sumstats_allchr.png",
    save_kwargs={"dpi": 300, "facecolor": "white"},
)

# ---- tidy output in LDSC format ---------------------------------------------
ss.sort_coordinate()
ss.sort_column()
ss.to_format(path=output_file, fmt="ldsc")

# ---- is this GWAS suitable for LD score regression? -------------------------
# the authors of LDSC recommend a mean chi-square of at least 1.02, see section 3
chi2 = np.asarray(ss.data["Z"], dtype=float) ** 2
print("variants kept :", len(chi2))
print("mean chi2     :", round(float(np.nanmean(chi2)), 4))
print("lambda GC     :", round(float(np.nanmedian(chi2) / 0.4549), 4))
```

</details>

And to run it:

<details>
<summary>terminal example: running the parser and checking the output</summary>

```console
(gwaslab) [duarte@node1 gwas]$ nano parse_plink_gwas.py
(gwaslab) [duarte@node1 gwas]$ python parse_plink_gwas.py
(gwaslab) [duarte@node1 gwas]$ ls out/gwaslab_output/
gwas_sumstats_allchr.ldsc.tsv.gz  gwas_sumstats_allchr.log  gwas_sumstats_allchr.png

(gwaslab) [duarte@node1 gwas]$ zcat out/gwaslab_output/gwas_sumstats_allchr.ldsc.tsv.gz | head -3
SNP	A1	A2	N	Z
chr1:11063:T:G	G	T	6412	0.4128
chr1:13259:G:A	A	G	6412	-1.2074
```

</details>

Two things to check here before moving on to the regression itself:

- **The mean chi-square.** If it comes out below ~1.02, LD score regression could give unstable estimates, since the GWAS could be underpowered to take out the most of LD score regression. For sure something to come back to in case we stumble upon this scenario in one of the cohorts, specially with samples sizes on the lower side... If you hit this, get in touch and we will discuss how to present it.
- **The Manhattan and QQ plots.** Inspect `gwas_sumstats_allchr.png`. What you want to see is a QQ plot that lifts gently along its whole length (polygenic signal) rather than one that is inflated from the very first quantile (uncorrected stratification). It would be nice if you could share this plot with us as well. 


## Step 7. run LDSC to obtain the heritability estimates

Everything is finally in place. We have the two ingredients the regression needs, and both of them were built in-sample and share the same principal components:

- the **summary statistics** from Step 6, in LDSC format,
- the **covariate-adjusted LD scores** from Step 5.b, one set per chromosome.

Note that from here on we call the *classic* `ldsc.py` from the Bulik-Sullivan repository, not the cov-LDSC one. The cov-LDSC job was to compute LD scores that are adjusted for the PCs and use a 20 cM window, and that job is done. The regression of chi-square on LD scores is the same operation in both, so we hand our adjusted scores to the standard implementation. This is also why we cloned both repositories back in Step 5.a, and why one conda environment serves both.

This step is cheap, so there is no scheduler version here, just run it in the terminal.

### Step 7.a prepare the folder and find your sample prevalence

Two prevalences go into the liability-scale transformation:

- `--samp-prev` is the **proportion of cases in your own sample**. This is not something you should estimate or guess: `make_covariate_files.r` already computed it in Step 4.c, on exactly the unrelated sample set that went into the GWAS. It is printed in that script's log as `Disease proportion`.
- `--pop-prev` is the **population prevalence of PD**, `K`. As explained in section 1.4, this is not measured with precision in underrepresented populations, so we do not pick one number. We run the whole thing three times, across `K ∈ {0.005, 0.01, 0.02}`, and report the spread. `K = 0.005` is our primary value.

<details>
<summary>terminal example: creating the folder and reading off the sample prevalence</summary>

```console
(base) [duarte@node1 gwas]$ cd ..
(base) [duarte@node1 working_directory_h2]$ mkdir covldsc_h2
(base) [duarte@node1 working_directory_h2]$ cd covldsc_h2/
(base) [duarte@node1 covldsc_h2]$ pwd
/home/duarte/working_directory_h2/covldsc_h2
(base) [duarte@node1 covldsc_h2]$ # the sample prevalence was already computed for you in step 4.c
(base) [duarte@node1 covldsc_h2]$ grep -A6 "Summary" ../pca_and_such/outFolder_pca_andSuch/covariate_file_no_related_pairs.log
--- Summary ---------------------------------------------
Covariate file          : ./pcair_r2_covariates_merged.tsv
Removal list            : ./NAToRA_output_pcrel_toRemove.txt (689 unique IDs)
Samples before removal  : 7101
Samples removed         : 689
Total sample size       : 6412
Cases    (DISEASE = 2)  : 3788
Controls (DISEASE = 1)  : 2624
Disease proportion      : 0.590768 (59.08%)

#the numbers are just examples...
```

</details>

### Step 7.b run the regression across the three prevalences

The script below reads that proportion straight out of the log 

<details>
<summary>script: <code>run_covldsc_h2.sh</code> — runs the LD score regression across the three population prevalences</summary>

```bash
#!/bin/bash
# Covariate-adjusted LD score regression for h2.
# Run it inside /working_directory_h2/covldsc_h2

# note: -u is deliberately left out, conda's activation script trips on unset variables
set -eo pipefail

# ---- inputs, all produced in previous steps ---------------------------------
sumstats="../gwas/out/gwaslab_output/gwas_sumstats_allchr.ldsc.tsv.gz"
ldscores="../reference_panel/chr@/out/covldsc_chr@"    # the @ is expanded by ldsc to 1..22
covar_log="../pca_and_such/outFolder_pca_andSuch/covariate_file_no_related_pairs.log"

# ---- programs ---------------------------------------------------------------
# the classic LDSC repository, cloned in step 5.a
ldsc="../programs/ldsc/ldsc.py"

# ---- population prevalences to scan (section 1.4) ---------------------------
# 0.005 is our primary value, the other two are the sensitivity analysis
pop_prevs="0.005 0.01 0.02"

mkdir -p ./logs

# ---- sample prevalence: read it from the step 4.c log ------------------------
# the line looks like: "Disease proportion      : 0.590768 (59.08%)"
samp_prev=$(awk '/Disease proportion/ {print $4}' "${covar_log}")

if [[ -z "${samp_prev}" ]]; then
    echo "ERROR: could not read the disease proportion from ${covar_log}" >&2
    echo "       set samp_prev by hand if you built your covariate files another way" >&2
    exit 1
fi
echo "sample prevalence read from step 4.c : ${samp_prev}"

# ---- fail early if anything is missing --------------------------------------
for f in "${sumstats}" "${ldsc}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing or empty: ${f}" >&2
        exit 1
    fi
done
if [[ ! -s "../reference_panel/chr22/out/covldsc_chr22.l2.ldscore.gz" ]]; then
    echo "ERROR: LD scores not found, did step 5.b finish for all 22 chromosomes?" >&2
    exit 1
fi

# ---- conda env, the same one used for cov-LDSC in step 5 --------------------
# if 'conda' is not on your PATH, replace this with the explicit path to your
# miniconda3/etc/profile.d/conda.sh
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ldsc

# ---- one regression per population prevalence -------------------------------
for k in ${pop_prevs}; do

    echo "=================== population prevalence K = ${k} ==================="

    python2 "${ldsc}" \
        --h2 "${sumstats}" \
        --ref-ld-chr "${ldscores}" \
        --w-ld-chr   "${ldscores}" \
        --samp-prev  "${samp_prev}" \
        --pop-prev   "${k}" \
        --out "./logs/cov_ldsc_popPrev_${k}"

done

echo "done. The estimates are in ./logs/cov_ldsc_popPrev_*.log"
```

</details>

In the terminal do

<details>
<summary>terminal example: saving and running the regression</summary>

```console
(base) [duarte@node1 covldsc_h2]$ nano run_covldsc_h2.sh
(base) [duarte@node1 covldsc_h2]$ bash run_covldsc_h2.sh
sample prevalence read from step 4.c : 0.590768
=================== population prevalence K = 0.005 ===================
=================== population prevalence K = 0.01 ===================
=================== population prevalence K = 0.02 ===================
done. The estimates are in ./logs/cov_ldsc_popPrev_*.log
(base) [duarte@node1 covldsc_h2]$ # if it doesn't run, try to make sure it is an executable
(base) [duarte@node1 covldsc_h2]$ chmod +x run_covldsc_h2.sh #and then try again to bash it
```

</details>

### Step 7.c read the output

LDSC writes its results into the `.log` file rather than a table, so this is where the actual heritability estimate lives.

<details>
<summary>terminal example: the LDSC log for the primary prevalence</summary>

```console
(base) [duarte@node1 covldsc_h2]$ cat logs/cov_ldsc_popPrev_0.005.log
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
*********************************************************************
Call:
./ldsc.py \
--h2 ../gwas/out/gwaslab_output/gwas_sumstats_allchr.ldsc.tsv.gz \
--ref-ld-chr ../reference_panel/chr@/out/covldsc_chr@ \
--out ./logs/cov_ldsc_popPrev_0.005 \
--samp-prev 0.590768 \
--w-ld-chr ../reference_panel/chr@/out/covldsc_chr@ \
--pop-prev 0.005

Reading summary statistics from ../gwas/out/gwaslab_output/gwas_sumstats_allchr.ldsc.tsv.gz ...
Read summary statistics for 8087341 SNPs.
Reading reference panel LD Score from ../reference_panel/chr@/out/covldsc_chr@ ... (ldscore_fromlist)
Read reference panel LD Scores for 8102072 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from ../reference_panel/chr@/out/covldsc_chr@ ... (ldscore_fromlist)
Read regression weight LD Scores for 8102072 SNPs.
After merging with reference panel LD, 8074915 SNPs remain.
After merging with regression SNP LD, 8074915 SNPs remain.
Total Observed scale h2: 0.1783 (0.0241)
Lambda GC: 1.0219
Mean Chi^2: 1.0431
Intercept: 0.9982 (0.0079)
Ratio < 0 (usually indicates GC correction).
Total Liability scale h2: 0.2416 (0.0327)
Total time elapsed: 36.0s
```

</details>

And to put the three runs side by side:

<details>
<summary>terminal example: collecting the three estimates into one table</summary>

```console
(base) [duarte@node1 covldsc_h2]$ for k in 0.005 0.01 0.02; do
>   printf "K=%-6s %s\n" "${k}" "$(grep 'Total Liability scale h2' logs/cov_ldsc_popPrev_${k}.log)"
> done
K=0.005  Total Liability scale h2: 0.2416 (0.0327)
K=0.01   Total Liability scale h2: 0.2698 (0.0365)
K=0.02   Total Liability scale h2: 0.3016 (0.0408)
```

</details>

Four numbers in that log deserve your attention:

| field | what to look for |
|---|---|
| `After merging with reference panel LD, N SNPs remain` | This should be close to the size of your reference panel. If it collapses to a few thousand, your summary statistics and your LD scores are not using the same variant IDs, and everything below it is meaningless. |
| `Intercept` | Should sit near 1. Substantially above 1 points to residual confounding that the ten PCs did not absorb, which in an admixed cohort is exactly the failure mode we are watching for. |
| `Total Observed scale h2` | The estimate on the 0/1 scale. Not comparable across studies, since it depends on your case-control ratio — do not report this one on its own. |
| `Total Liability scale h2` | **This is the number we want**, together with its standard error in parentheses. Reported per prevalence. |

Finally, please share all three `.log` files with us. Together with the diagnostics from Steps 4 and 6, they are what lets us tell apart a genuinely low heritability from an underpowered GWAS.


# Second wave of analysis: GRM based estimates 

With the summary-statistic and LD scores arm behind us, we now move to arms 2 to 5 of the design table in section 4, the ones that estimate heritability from a Genomic Relatedness Matrix. This part of the tutorial covers **arm 2: GREML on a standard GCTA GRM**, which is one of the classical field's reference implementation and serves as a complementary source to validate LDSC estimates. 

Two things are worth stating before we start:

- **Same samples, same variants.** We build the GRM from the very same `bed/bim/fam` filesets produced in Step 5.b and used for the GWAS in Step 6. They are already restricted to the unrelated individuals, already MAF-filtered at 0.01, and already carry the final variant IDs. So arm 2 and arm 1 differ in the estimator, not in the data.
- **Same covariates.** The PCs and the other covariates handed to GREML come from the same files written in Step 4.c that fed the GWAS. If you added a cohort-specific covariate back in Step 6 (a batch or phase variable), it has to travel with you here too.

One consequence of the first point is worth anticipating: the GRM is computed over the full reference panel variant set, imputed variants included, which is a lot more variants than a genotyped-only GRM would use. That is deliberate, but it does make this the second most computationally demanding step of the tutorial after Step 5.b.

## Step 8. download GCTA and compute a classic GRM

### Step 8.a get the program

As with plink, we keep the executable inside `./programs` so it does not interfere with any GCTA you may already have installed.

<details>
<summary>terminal example: downloading GCTA and checking the executable</summary>

```console
(base) [duarte@node1 working_directory_h2]$ cd programs/
(base) [duarte@node1 programs]$ wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.95.3-linux-x86_64.zip 
(base) [duarte@node1 programs]$ unzip gcta-1.95.3-linux-x86_64.zip 
(base) [duarte@node1 programs]$ ./gcta-1.95.3-linux-x86_64/gcta 
*******************************************************************
* Genome-wide Complex Trait Analysis (GCTA)
* version v1.95.3 Linux
* Built at Jul 21 2026 11:23:29, by GCC 8.4
* (C) 2010-present, Yang Lab, Westlake University
* Please report bugs to Jian Yang <jian.yang@westlake.edu.cn>
*******************************************************************
Hostname: r01-ccf.org

Error: no analysis has been launched by the option(s)
Please see online documentation at https://yanglab.westlake.edu.cn/software/gcta/
```

</details>

Same story as with LDSC back in Step 5.a: the error message just means we launched GCTA without asking it to do anything, which is exactly the check we wanted. The executable works.

### Step 8.b create the folder for this analysis

<details>
<summary>terminal example: creating the greml folder</summary>

```console
(base) [duarte@node1 programs]$ cd ..
(base) [duarte@node1 working_directory_h2]$ mkdir greml
(base) [duarte@node1 working_directory_h2]$ cd greml/
(base) [duarte@node1 greml]$ pwd
/home/duarte/working_directory_h2/greml
(base) [duarte@node1 greml]$ ls ..
covldsc_h2  genetic_data  greml  gwas  pca_and_such  programs  reference_panel
```

</details>

### Step 8.c compute one GRM per chromosome

A genome-wide GRM could in principle be built in a single GCTA call, but on a cohort of this size that means holding the whole variant set in memory at once, and we will do it only when it is necesary (see GRM with PC-Relate below). It is both faster and far more forgiving to build one GRM per chromosome and merge them afterwards, which is what GCTA's `--mgrm` is designed for. Splitting also means that a chromosome that runs out of memory can be rerun on its own instead of restarting everything.

Note that we do not pass `--autosome` or `--maf` here. The per-chromosome filesets are already autosome-only by construction (we only ever looped over 1 to 22) and were already filtered at MAF > 0.01 in Step 5.b. Adding the flags again would be harmless but redundant, and keeping them out makes it explicit that the variant set is inherited from the reference panel.

As a python script to submit 22 different jobs

<details>
<summary>script: <code>make_grms_per_chr.py</code> — writes one SLURM .pbs per chromosome and submits all 22</summary>

```python
# Build one GCTA GRM per chromosome, from the same filesets used for the
# reference panel and the GWAS.
# Run it inside /working_directory_h2/greml

import os

# ---- cluster resources: adjust to what your scheduler offers -----------------
threads   = 4
mem       = "25G"
email     = "your_mail@email.org"
partition = "bigmem"

# ---- programs, as relative paths from ./greml -------------------------------
gcta = "../programs/gcta-1.95.3-linux-x86_64/gcta"

for i in range(1, 23):
    runFolder = f"./chr{i}"
    # the reference panel filesets from step 5.b, used here as-is
    plinkFile = f"../reference_panel/chr{i}/chr{i}.plink1.extracted_cmUpdated"
    grmOut    = f"{runFolder}/out/chr{i}"

    os.makedirs(f"{runFolder}/logs", exist_ok=True)
    os.makedirs(f"{runFolder}/out",  exist_ok=True)

    fileSbatch = open(f"{runFolder}/logs/{i}make_grms_per_chr.pbs", "w")

    fileSbatch.write(f"#!/bin/bash\n")
    fileSbatch.write(f"#SBATCH --mail-type=END,FAIL\n")
    fileSbatch.write(f"#SBATCH --mail-user={email}\n")
    fileSbatch.write(f"#SBATCH --ntasks=1\n")
    fileSbatch.write(f"#SBATCH --cpus-per-task={threads}\n")
    fileSbatch.write(f"#SBATCH --mem={mem}\n")
    fileSbatch.write(f"#SBATCH --partition={partition}\n")
    fileSbatch.write(f"#SBATCH --job-name={i}make_grms_per_chr\n")
    fileSbatch.write(f"#SBATCH -o {runFolder}/logs/{i}make_grms_per_chr.out\n")
    fileSbatch.write(f"#SBATCH -e {runFolder}/logs/{i}make_grms_per_chr.err\n\n")

    fileSbatch.write(f"""set -uo pipefail

BFILE="{plinkFile}"
GRM_OUT="{grmOut}"

# Confirm that the chromosome-specific PLINK files exist.
for extension in bed bim fam; do
    if [[ ! -s "${{BFILE}}.${{extension}}" ]]; then
        echo "ERROR: Missing or empty file:" >&2
        echo "  ${{BFILE}}.${{extension}}" >&2
        exit 1
    fi
done

echo "Number of samples:"
wc -l "${{BFILE}}.fam"

echo "Number of variants:"
wc -l "${{BFILE}}.bim"

echo
echo "Constructing chromosome {i} GRM..."

{gcta} --bfile "${{BFILE}}" \\
    --make-grm \\
    --thread-num {threads} \\
    --out "${{GRM_OUT}}"

# Verify the required GCTA output files.
for extension in grm.bin grm.N.bin grm.id; do
    if [[ ! -s "${{GRM_OUT}}.${{extension}}" ]]; then
        echo "ERROR: GCTA did not generate:" >&2
        echo "  ${{GRM_OUT}}.${{extension}}" >&2
        exit 1
    fi
done

echo
echo "Successfully generated chromosome {i} GRM."
echo "Finished: $(date)"
""")

    fileSbatch.close()
    os.system(f"sbatch {runFolder}/logs/{i}make_grms_per_chr.pbs")
```

</details>

For running it do:

<details>
<summary>terminal example: saving the python script with nano and submitting the 22 jobs</summary>

```console
(base) [duarte@node1 greml]$ nano make_grms_per_chr.py
(base) [duarte@node1 greml]$ python make_grms_per_chr.py
Submitted batch job 6116740
Submitted batch job 6116741
Submitted batch job 6116742
Submitted batch job 6116743
Submitted batch job 6116744
Submitted batch job 6116745
Submitted batch job 6116746
Submitted batch job 6116747
Submitted batch job 6116748
Submitted batch job 6116749
Submitted batch job 6116750
Submitted batch job 6116751
Submitted batch job 6116752
Submitted batch job 6116753
Submitted batch job 6116754
Submitted batch job 6116755
Submitted batch job 6116756
Submitted batch job 6116757
Submitted batch job 6116758
Submitted batch job 6116759
Submitted batch job 6116760
Submitted batch job 6116761

```

</details>

If you prefer to work in series through your terminal here is the equivalent for loop.

<details>
<summary>script: <code>make_grms_per_chr.sh</code> — same per-chromosome GRMs as a sequential for loop, no scheduler needed</summary>

```bash
#!/bin/bash
# Sequential version of make_grms_per_chr.py
# Runs chr1 .. chr22 one after the other in the terminal (no scheduler).
# Run it inside /working_directory_h2/greml

set -uo pipefail

threads=4

# ---- programs, as relative paths from ./greml -------------------------------
gcta="../programs/gcta-1.95.3-linux-x86_64/gcta"

if [[ ! -x "${gcta}" ]]; then
    echo "ERROR: GCTA binary not found or not executable: ${gcta}" >&2
    exit 1
fi

for i in $(seq 1 22); do

    BFILE="../reference_panel/chr${i}/chr${i}.plink1.extracted_cmUpdated"
    GRM_OUT="./chr${i}/out/chr${i}"
    mkdir -p "./chr${i}/out" "./chr${i}/logs"

    for ext in bed bim fam; do
        if [[ ! -s "${BFILE}.${ext}" ]]; then
            echo "SKIP chr${i}: missing or empty ${BFILE}.${ext}" >&2
            continue 2
        fi
    done

    echo "=== chr${i}: $(wc -l < "${BFILE}.fam") samples, $(wc -l < "${BFILE}.bim") variants ==="
    echo "chr${i}: started $(date)"

    if ! "${gcta}" --bfile "${BFILE}" \
              --make-grm \
              --thread-num "${threads}" \
              --out "${GRM_OUT}" \
              > "./chr${i}/logs/chr${i}.grm.log" 2>&1; then
        echo "ERROR chr${i}: gcta exited non-zero, see ./chr${i}/logs/chr${i}.grm.log" >&2
        continue
    fi

    for ext in grm.bin grm.N.bin grm.id; do
        if [[ ! -s "${GRM_OUT}.${ext}" ]]; then
            echo "ERROR chr${i}: gcta did not generate ${GRM_OUT}.${ext}" >&2
            continue 2
        fi
    done

    echo "chr${i}: OK, finished $(date)"
    echo
done

echo "All chromosomes done."
```

</details>

In the terminal do

<details>
<summary>terminal example: saving and running the bash version</summary>

```console
(base) [duarte@node1 greml]$ nano make_grms_per_chr.sh
(base) [duarte@node1 greml]$ bash make_grms_per_chr.sh
```

</details>

Either way you end up with one folder per chromosome, each holding a GRM in GCTA's three-file format.

<details>
<summary>terminal example: what the greml directory looks like</summary>

```console
(base) [duarte@node1 greml]$ ls
chr1   chr11  chr13  chr15  chr17  chr19  chr20  chr22  chr4  chr6  chr8  make_grms_per_chr.py
chr10  chr12  chr14  chr16  chr18  chr2   chr21  chr3   chr5  chr7  chr9  make_grms_per_chr.sh
(base) [duarte@node1 greml]$ # sanity check: 22 GRMs, and all of them with the same people
(base) [duarte@node1 greml]$ ls chr*/out/*.grm.bin | wc -l
22
```

</details>

The three files are the GRM itself (`.grm.bin`, the lower triangle of the matrix in binary), the number of variants used for each pair (`.grm.N.bin`), and the sample IDs in matrix order (`.grm.id`). All three are needed, and the merge in the next step will refuse to run unless the `.grm.id` files agree across all 22 chromosomes.

## Step 9. merge the 22 GRMs and run GREML across the prevalence grid

Now we hand the 22 chromosome GRMs to GCTA so it can add them into a single genome-wide GRM, and then fit the REML model on it. As in Step 7, the model itself is fit once on the observed scale, and only the liability-scale transformation changes across the prevalence grid — so the three `K` values are cheap.

The script needs three files from Step 4.c, and it is worth knowing what each one is because GCTA is strict about their shape and gives unhelpful errors when they are wrong:

| flag | file | shape |
|---|---|---|
| `--pheno` | `covariate_file_no_related_pairs_pheno.tsv` | `FID IID DISEASE`, no header. Must be coded 1 = control, 2 = case, otherwise `--prevalence` is silently ignored and you get no liability-scale row |
| `--qcovar` | `covariate_file_no_related_pairs_qcovar.tsv` | `FID IID AGE PC1..PC10`, no header. The quantitative covariates |
| `--covar` | `covariate_file_no_related_pairs_covar.tsv` | `FID IID SEX`, no header. The categorical covariates |

<details>
<summary>script: <code>merge_grms_runGREML.sh</code> — merges the 22 GRMs and runs single-GRM GREML across the prevalence grid</summary>

```bash
#!/bin/bash
#
# Merge 22 per-chromosome GCTA GRMs into a single genome-wide GRM and run
# single-GRM GREML across a prevalence grid.
#
# Run it inside /working_directory_h2/greml:
#     bash merge_grms_runGREML.sh
#
# All paths default to the directory containing this script, so it can be
# invoked from anywhere. Any of WORKDIR, GRM_DIR, PCA_DIR, OUT_ROOT,
# GCTA or THREADS can be overridden from the environment, e.g.
#     THREADS=12 GRM_DIR=/some/other/path bash merge_grms_runGREML.sh

set -euo pipefail


# ---------------------------------------------------------------------
# Paths and runtime settings
# ---------------------------------------------------------------------
# Anchor to the script's own location rather than the caller's CWD.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${WORKDIR:-${SCRIPT_DIR}}"

THREADS="${THREADS:-4}"

# The GCTA executable downloaded in step 8.a, kept in ./programs like plink
GCTA="${GCTA:-${WORKDIR}/../programs/gcta-1.95.3-linux-x86_64/gcta}"

# Directory holding the per-chromosome GRMs written by
# make_grms_per_chr.py as chr<N>/out/chr<N>.*
GRM_DIR="${GRM_DIR:-${WORKDIR}}"

# Root for everything this script writes.
OUT_ROOT="${OUT_ROOT:-${WORKDIR}/gcta_singleGRM}"
MERGED_DIR="${OUT_ROOT}/genomewide"
GREML_DIR="${OUT_ROOT}/greml"
LOG_DIR="${OUT_ROOT}/logs"

mkdir -p "${MERGED_DIR}" "${GREML_DIR}" "${LOG_DIR}"


# ---------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------
# Everything written to stdout/stderr from this point on goes to BOTH the
# terminal and a persistent log. Fixed filename: each run overwrites the
# previous log.

SCRIPT_LOG="${LOG_DIR}/merge_GREML.log"

exec > >(tee "${SCRIPT_LOG}") 2>&1

echo "=============================================================="
echo "Script      : ${SCRIPT_DIR}/$(basename "${BASH_SOURCE[0]}")"
echo "Started     : $(date)"
echo "Host        : $(hostname)"
echo "Working dir : ${WORKDIR}"
echo "GRM dir     : ${GRM_DIR}"
echo "Output root : ${OUT_ROOT}"
echo "Threads     : ${THREADS}"
echo "Log file    : ${SCRIPT_LOG}"
echo "=============================================================="

if [[ ! -x "${GCTA}" ]]; then
    echo "ERROR: GCTA binary not found or not executable: ${GCTA}" >&2
    exit 1
fi

GRM_LIST="${MERGED_DIR}/chromosome_grms.txt"
MERGED_GRM="${MERGED_DIR}/genomewide_single_grm"


# ---------------------------------------------------------------------
# GCTA input files, all written by make_covariate_files.r in step 4.c
# ---------------------------------------------------------------------
PCA_DIR="${PCA_DIR:-${WORKDIR}/../pca_and_such/outFolder_pca_andSuch}"

# No header: FID IID phenotype
# Binary traits must be coded 1=control, 2=case (missing as NA or -9),
# otherwise --prevalence is ignored and no liability-scale row is written.
PHENO="${PCA_DIR}/covariate_file_no_related_pairs_pheno.tsv"

# Quantitative covariates, no header: FID IID AGE PC1 PC2 ... PC10
QCOVAR="${PCA_DIR}/covariate_file_no_related_pairs_qcovar.tsv"

# Optional categorical covariates, no header: FID IID SEX
COVAR="${PCA_DIR}/covariate_file_no_related_pairs_covar.tsv"

# Prevalence grid. The first value is the primary (field-standard for PD);
# the remaining values are declared sensitivity analyses. REML itself is
# fit on the observed scale, so only the liability-scale transformation
# (V(G)/Vp_L) differs between these runs.
PREVALENCES=(0.005 0.01 0.02)

# Output prefix stem; one prefix per prevalence is derived from this.
GREML_STEM="${GREML_DIR}/PD_singleGRM_GREML"
SUMMARY_TSV="${GREML_DIR}/PD_singleGRM_GREML_prevalence_grid.tsv"


# ---------------------------------------------------------------------
# Validate chromosome GRMs and create GRM list
# ---------------------------------------------------------------------
: > "${GRM_LIST}"

REFERENCE_ID="${GRM_DIR}/chr1/out/chr1.grm.id"

if [[ ! -s "${REFERENCE_ID}" ]]; then
    echo "ERROR: Reference ID file is missing: ${REFERENCE_ID}" >&2
    echo "       Set GRM_DIR if the per-chromosome GRMs live elsewhere." >&2
    exit 1
fi

for CHR in {1..22}; do
    GRM_PREFIX="${GRM_DIR}/chr${CHR}/out/chr${CHR}"

    for extension in grm.bin grm.N.bin grm.id; do
        if [[ ! -s "${GRM_PREFIX}.${extension}" ]]; then
            echo "ERROR: Missing chromosome ${CHR} output:" >&2
            echo "       ${GRM_PREFIX}.${extension}" >&2
            exit 1
        fi
    done

    # All chromosome GRMs should contain exactly the same people
    # in exactly the same order.
    if ! cmp -s "${REFERENCE_ID}" "${GRM_PREFIX}.grm.id"; then
        echo "ERROR: Sample IDs/order differ for chromosome ${CHR}." >&2
        echo "Compare:" >&2
        echo "  ${REFERENCE_ID}" >&2
        echo "  ${GRM_PREFIX}.grm.id" >&2
        exit 1
    fi

    echo "${GRM_PREFIX}" >> "${GRM_LIST}"
done

echo
echo "Chromosome GRMs:"
cat "${GRM_LIST}"

echo
echo "Number of samples per chromosome GRM:"
wc -l "${REFERENCE_ID}"


# ---------------------------------------------------------------------
# Merge 22 chromosome GRMs into ONE genome-wide GRM
# ---------------------------------------------------------------------
echo
echo "########## Merging chromosome GRMs ##########"

"${GCTA}" \
    --mgrm "${GRM_LIST}" \
    --make-grm \
    --thread-num "${THREADS}" \
    --out "${MERGED_GRM}"

for extension in grm.bin grm.N.bin grm.id; do
    if [[ ! -s "${MERGED_GRM}.${extension}" ]]; then
        echo "ERROR: Merged GRM file was not created:" >&2
        echo "       ${MERGED_GRM}.${extension}" >&2
        exit 1
    fi
done

echo "Merged GRM contains:"
echo "  Samples: $(wc -l < "${MERGED_GRM}.grm.id")"


# ---------------------------------------------------------------------
# Validate phenotype/covariate files
# ---------------------------------------------------------------------
for required_file in "${PHENO}" "${QCOVAR}"; do
    if [[ ! -s "${required_file}" ]]; then
        echo "ERROR: Missing required file: ${required_file}" >&2
        echo "       Set PCA_DIR if the covariate files live elsewhere." >&2
        exit 1
    fi
done

USE_COVAR=1
if [[ ! -s "${COVAR}" ]]; then
    USE_COVAR=0
    echo "WARNING: categorical covariate file is missing or empty:" >&2
    echo "         ${COVAR}" >&2
    echo "WARNING: running GREML WITHOUT --covar (no sex adjustment)." >&2
fi


# ---------------------------------------------------------------------
# Run classical single-GRM GREML across the prevalence grid
# ---------------------------------------------------------------------
FAILED_K=()

for PREVALENCE in "${PREVALENCES[@]}"; do

    GREML_OUT="${GREML_STEM}_K${PREVALENCE}"

    echo
    echo "=============================================================="
    echo "GREML | prevalence K = ${PREVALENCE}"
    echo "Output prefix: ${GREML_OUT}"
    echo "=============================================================="

    GREML_COMMAND=(
        "${GCTA}"
        --reml
        --grm "${MERGED_GRM}"
        --pheno "${PHENO}"
        --qcovar "${QCOVAR}"
        --prevalence "${PREVALENCE}"
        --thread-num "${THREADS}"
        --out "${GREML_OUT}"
    )

    # Add categorical covariates only when the file exists and is nonempty.
    if [[ "${USE_COVAR}" -eq 1 ]]; then
        GREML_COMMAND+=(--covar "${COVAR}")
    fi

    echo "Running:"
    printf '%q ' "${GREML_COMMAND[@]}"
    echo

    # Do not let one non-converging prevalence abort the whole grid.
    if ! "${GREML_COMMAND[@]}"; then
        echo "ERROR: GCTA returned a nonzero exit status for K=${PREVALENCE}." >&2
        FAILED_K+=("${PREVALENCE}")
        continue
    fi

    if [[ ! -s "${GREML_OUT}.hsq" ]]; then
        echo "ERROR: GREML did not generate ${GREML_OUT}.hsq" >&2
        FAILED_K+=("${PREVALENCE}")
        continue
    fi

    echo
    echo "########## GREML results | K = ${PREVALENCE} ##########"
    cat "${GREML_OUT}.hsq"

done


# ---------------------------------------------------------------------
# Collate the grid into a single summary table
# ---------------------------------------------------------------------
echo
echo "########## Prevalence grid summary ##########"

printf 'K\tsource\testimate\tSE\n' > "${SUMMARY_TSV}"

for PREVALENCE in "${PREVALENCES[@]}"; do
    HSQ="${GREML_STEM}_K${PREVALENCE}.hsq"
    [[ -s "${HSQ}" ]] || continue

    awk -v k="${PREVALENCE}" '
        $1 == "V(G)/Vp"   { printf "%s\t%s\t%s\t%s\n", k, $1, $2, $3 }
        $1 == "V(G)/Vp_L" { printf "%s\t%s\t%s\t%s\n", k, $1, $2, $3 }
        $1 == "Pval"      { printf "%s\t%s\t%s\tNA\n", k, $1, $2 }
        $1 == "n"         { printf "%s\t%s\t%s\tNA\n", k, $1, $2 }
    ' "${HSQ}" >> "${SUMMARY_TSV}"
done

column -t "${SUMMARY_TSV}" || cat "${SUMMARY_TSV}"

echo
echo "Summary table written to: ${SUMMARY_TSV}"

if [[ "${#FAILED_K[@]}" -gt 0 ]]; then
    echo
    echo "ERROR: GREML failed for prevalence value(s): ${FAILED_K[*]}" >&2
    echo "Finished with errors: $(date)"
    exit 1
fi

echo
echo "Finished: $(date)"
echo "Full log: ${SCRIPT_LOG}"

# Give the tee subprocess a moment to flush before the shell exits.
sleep 1
```

</details>

In the terminal do

<details>
<summary>terminal example: saving and running the merge plus GREML script</summary>

```console
(base) [duarte@node1 greml]$ nano merge_grms_runGREML.sh
(base) [duarte@node1 greml]$ bash merge_grms_runGREML.sh
==============================================================
Script      : /home/duarte/working_directory_h2/greml/merge_grms_runGREML.sh
Host        : node1.lerner.ccf.org
Working dir : /home/duarte/working_directory_h2/greml
GRM dir     : /home/duarte/working_directory_h2/greml
Output root : /home/duarte/working_directory_h2/greml/gcta_singleGRM
Threads     : 4
Log file    : /home/duarte/working_directory_h2/greml/gcta_singleGRM/logs/merge_GREML.log
==============================================================

Chromosome GRMs:
/home/duarte/working_directory_h2/greml/chr1/out/chr1
/home/duarte/working_directory_h2/greml/chr2/out/chr2
...
/home/duarte/working_directory_h2/greml/chr22/out/chr22

Number of samples per chromosome GRM:
6412 /home/duarte/working_directory_h2/greml/chr1/out/chr1.grm.id

########## Merging chromosome GRMs ##########
Merged GRM contains:
  Samples: 6412
```

</details>

Once it finishes, the results live in `gcta_singleGRM/greml/`, one `.hsq` per prevalence plus the collated table.

Reading the `.hsq`, four things matter:

| field | what to look for |
|---|---|
| `n` | The number of individuals that actually entered the model. It should match the unrelated sample count from Step 4.c. If it is smaller, some IDs failed to match between the GRM and the phenotype file |
| `V(G)/Vp` | The observed-scale estimate. Note it is identical across the three runs, as it must be: the prevalence only enters the transformation, not the REML fit |
| `V(G)/Vp_L` | **This is the number we report**, together with its SE. One per prevalence |
| `Pval` | The likelihood ratio test against the null of no genetic variance |

Two things are worth noticing when you compare this with the cov-LDSC result from Step 7. First,the GREML liability estimate will very likely come out higher than the cov-LDSC one. And second, that gap is not an error either, it is the expected behaviour described in section 3, and characterising its size in underrepresented cohorts is one of the aims of this project.

As before, please share the three `.hsq` files and the summary table with us.

## Step 10. build the ancestry-aware GRM with PC-Relate

This is **arm 3** of the design table: the same estimator, the same samples, the same variants and the same covariates as arm 2, with exactly one thing changed, namely how the GRM is built. That is the whole point of running both.

Recall the argument from there: a standard GCTA GRM compares every individual against a single cohort-wide allele frequency. In an admixed cohort that reference describes nobody, and two unrelated people who happen to share a high proportion of the same ancestry will look genetically similar for reasons that have nothing to do with recent genealogy. PC-Relate instead estimates **individual-specific allele frequencies** from the principal components (comming from PC-Air), and measures each pair's similarity against their own expected frequencies. The ancestry component is removed rather than absorbed.

There is one practical difference from arm 2 that shapes this step. GCTA lets us build 22 chromosome GRMs and add them together afterwards with `--mgrm`. PC-Relate has no such option: it estimates individual-specific allele frequencies across the whole genome at once, so it needs a single genome-wide fileset. Hence Step 10.a.

### Step 10.a merge the 22 reference panel filesets into one

We merge the same per-chromosome `bed/bim/fam` filesets from Step 5.b that fed the GWAS and the GCTA GRMs. Because the 22 filesets cover disjoint chromosomes and contain identical samples, this is a straight concatenation. 

Two things to be aware of before you launch it. The merged `.bed` is large: roughly 15 GB for 7.1k samples and 8M variants. Make sure you have the disk space. 

<details>
<summary>script: <code>merge_reference_panel.sh</code> — merges the 22 per-chromosome filesets into one genome-wide fileset</summary>

```bash
#!/bin/bash
# Merge the 22 per-chromosome reference panel filesets from step 5.b into a
# single genome-wide fileset, which is what PC-Relate needs.
# Run it inside /working_directory_h2/reference_panel

set -euo pipefail

threads=4
mem_mb=50000

# ---- programs, as relative paths from ./reference_panel ---------------------
plink1="../programs/plink"

out_prefix="merged_chr1_22"
merge_list="./merge_list.txt"

# ---- fail early if any chromosome is missing or incomplete ------------------
for i in $(seq 1 22); do
    for ext in bed bim fam; do
        f="chr${i}/chr${i}.plink1.extracted_cmUpdated.${ext}"
        if [[ ! -s "${f}" ]]; then
            echo "ERROR: missing or empty ${f}" >&2
            echo "       did step 5.b finish for all 22 chromosomes?" >&2
            exit 1
        fi
    done
done

# ---- build the merge list --------------------------------------------------
# chr1 goes in --bfile, so the list holds chr2..chr22, one bed/bim/fam triple
# per line
: > "${merge_list}"
for i in $(seq 2 22); do
    stem="chr${i}/chr${i}.plink1.extracted_cmUpdated"
    echo "${stem}.bed ${stem}.bim ${stem}.fam" >> "${merge_list}"
done

echo "=== merging chr1-22 into ${out_prefix} ==="
"${plink1}" \
    --bfile chr1/chr1.plink1.extracted_cmUpdated \
    --merge-list "${merge_list}" \
    --make-bed \
    --out "${out_prefix}" \
    --memory "${mem_mb}" \
    --threads "${threads}"

echo "=== merged fileset: ${out_prefix}.bed / .bim / .fam ==="
echo "variants: $(wc -l < "${out_prefix}.bim")"
echo "samples:  $(wc -l < "${out_prefix}.fam")"

echo "done"
```

</details>

In the terminal do

<details>
<summary>terminal example: saving and running the merge</summary>

```console
(base) [duarte@node1 greml]$ cd ../reference_panel/
(base) [duarte@node1 reference_panel]$ nano merge_reference_panel.sh
(base) [duarte@node1 reference_panel]$ bash merge_reference_panel.sh
=== merging chr1-22 into merged_chr1_22 ===
=== merged fileset: merged_chr1_22.bed / .bim / .fam ===
variants: 8102072
samples:  6412
done
(base) [duarte@node1 reference_panel]$ du -h merged_chr1_22.*
13G	merged_chr1_22.bed
289M	merged_chr1_22.bim
162K	merged_chr1_22.fam
3.1K	merged_chr1_22.log
```

</details>

Check that the variant count matches the total you got at the end of Step 5.b. If it is lower, a chromosome silently failed to merge and the `.log` will say which.

### Step 10.b compute the PC-Relate GRM

This step goes back to the `genesis` conda environment we created in Step 4, since PC-Relate lives in the same GENESIS package that gave us PC-AiR.

<details>
<summary>terminal example: creating the folder and activating the environment</summary>

```console
(base) [duarte@node1 reference_panel]$ cd ..
(base) [duarte@node1 working_directory_h2]$ mkdir greml_pcrelate
(base) [duarte@node1 working_directory_h2]$ cd greml_pcrelate/
(base) [duarte@node1 greml_pcrelate]$ conda activate genesis
(genesis) [duarte@node1 greml_pcrelate]$ pwd
/home/duarte/working_directory_h2/greml_pcrelate
```

</details>

The script reuses several things we already computed, so nothing here is recalculated from scratch:

| input | where it comes from |
|---|---|
| `../reference_panel/merged_chr1_22` | Step 10.a, just now |
| `mypcair_r2_results_rds` | Step 4.a, the second round of PC-AiR |
| `unrelated_IIDs.txt` | Step 4.c, the samples that survived NAToRA |
| `pcair_r2_covariates_merged.tsv` | Step 4.a, the merged covariate table |

**A word on covariates before you run it.** PC-Relate itself takes no covariates at all: it takes PCs and genotypes and returns kinship coefficients. So if your cohort carries extra covariates beyond `SEX`, `AGE` and the PCs, **the GRM does not change and you do not need to touch that part of the script**. Covariates enter later, when GCTA fits the model, through the files that `make_covariate_files.r` already wrote in Step 4.c. The one place where extra covariates do matter is the block at the very end of the script, which writes `pd.ldak.covar`, an all-numeric covariate file with categorical variables expressed as 0/1 dummies. **This file is required later**: arms 4 to 7 run PCGC through the LDAK program, and LDAK does not read the GCTA-style covariate files that `make_covariate_files.r` wrote in Step 4.c. Hence declare the extra covariates that you want in the `extra_covars` variable at the top of that block rather than editing the code.

A point of vocabulary that trips people up here: **section 4 rules out the LDAK heritability *model*, not the LDAK *software***. What we excluded from the benchmark is the LDAK-Thin model with `α = −0.25`, because its SNP weights are baked into its own GRM construction and so cannot be benchmarked side by side with an ancestry-aware GRM. We still use the LDAK program itself, since it carries the PCGC implementation we need, run on our own GRMs under the standard `α = −1` model.

Unlike the PCA pipeline in Step 4, this script has no command-line interface. Open the config block at the top and edit it before running, in particular `n_cores`.

For 7.1K samples and ~8M variants, the total computational time took me around ~4 hours at 12 threads and 94G of memory. 

<details>
<summary><strong>View script: <code>grm_with_pcrelate.r</code></strong></summary>

<br>

```r
# Ancestry-aware GRM for the h2 project, built with PC-Relate.
# Run it inside /working_directory_h2/greml_pcrelate, with the genesis
# conda environment active.

library(SNPRelate); library(GENESIS); library(GWASTools)
library(gdsfmt); library(BiocParallel); library(Matrix); library(data.table)

# --- config -------------------------------------------------------------
pca_out_dir  <- "../pca_and_such/outFolder_pca_andSuch/"   # output of the PCA pipeline, step 4
h2_out_dir   <- "./outFolder_h2"
dir.create(h2_out_dir, showWarnings = FALSE, recursive = TRUE)

# the cov-LDSC variant set: imputed r2>0.8 U genotyped U HapMap3, merged in step 10.a
h2_input     <- "../reference_panel/merged_chr1_22"        # plink prefix
h2_input_fmt <- "bfile"

pcair_r2_rds <- file.path(pca_out_dir, "mypcair_r2_results_rds")
unrel_file   <- file.path(pca_out_dir, "unrelated_IIDs.txt")   # one IID per line
covar_file   <- file.path(pca_out_dir, "pcair_r2_covariates_merged.tsv")
force_pcrel  <- FALSE                     # TRUE to recompute even if the .rds exists

n_pcs        <- 10
n_cores      <- 12                        # set this to what your machine actually has
maf_min      <- 0.01                      # matches the cov-LDSC LD scores
pcrel_scale  <- "variant"                 # 'overall' for sensitivity

op <- function(...) file.path(h2_out_dir, ...)


# --- logging ------------------------------------------------------------
# On by default: everything the script prints (cat/print) plus the messages and
# warnings raised by GENESIS et al. is mirrored to a timestamped log under
# h2_out_dir, while still going to the console. Override the path with
# H2_LOG=/path/to/file, or set H2_LOG=none to turn logging off.

log_path  <- Sys.getenv("H2_LOG", unset = "")
log_close <- function() invisible(NULL)          # no-op unless logging is on

if (!identical(log_path, "none")) {
  if (!nzchar(log_path))
    log_path <- op(format(Sys.time(), "grm_with_pcrel_%Y%m%d-%H%M%S.log"))

  log_con <- file(log_path, open = "wt")
  sink(log_con, split = TRUE)                    # cat/print -> console AND log

  # Writing straight to log_con while it is a split sink target would echo the
  # text to the console as well, so suspend the split around each write.
  log_note <- function(txt) {
    sink()
    cat(txt, file = log_con, sep = ""); flush(log_con)
    sink(log_con, split = TRUE)
  }

  # message() and warning() write to stderr, which sink() cannot split, so
  # mirror them by hand and let them through to the console untouched. This is
  # what captures pcrelate(verbose = TRUE) progress.
  globalCallingHandlers(
    message = function(m) log_note(conditionMessage(m)),
    warning = function(w) log_note(paste0("Warning: ", conditionMessage(w), "\n"))
  )

  log_close <- function() {
    while (sink.number() > 0) sink()
    if (isOpen(log_con)) close(log_con)
  }

  # An error unwinds the sinks so the console is left in a sane state, and the
  # message lands in the log rather than only on the terminal.
  options(error = function() {
    log_note(paste0("\nERROR: ", geterrmessage()))
    log_close()
    if (!interactive()) quit(save = "no", status = 1)
  })

  cat("logging to:", log_path, "\n")
}

cat("started:", format(Sys.time()), "\n")


# --- load the PC-AiR round 2 results, no need to recompute --------------

pcair_r2  <- readRDS(pcair_r2_rds)
unrel_ids <- readLines(unrel_file)

pcs_all <- pcair_r2$vectors[, 1:n_pcs, drop = FALSE]
stopifnot(all(unrel_ids %in% rownames(pcs_all)))
pcs_sub <- pcs_all[unrel_ids, , drop = FALSE]   # rownames = sample IDs; GENESIS matches on these

cat("unrelated samples:", length(unrel_ids), " PCs:", ncol(pcs_sub), "\n")


# --- convert to GDS -----------------------------------------------------
# No MAF filter or --autosome here: the merged fileset already inherits both
# from step 5.b, exactly as in the GCTA arm.

h2_gds <- op("h2_common.gds")

if (!file.exists(h2_gds)) {
  snpgdsBED2GDS(bed.fn = paste0(h2_input, ".bed"),
                bim.fn = paste0(h2_input, ".bim"),
                fam.fn = paste0(h2_input, ".fam"),
                out.gdsfn = h2_gds)
}
snpgdsSummary(h2_gds)

M <- nrow(fread(paste0(h2_input, ".bim"), select = 1))
cat("SNPs in GRM:", M, "\n")


# --- PC-Relate (cached: skipped if the .rds is already on disk) ----------
pcrel_rds <- op(sprintf("pcrel_h2_%s.rds", pcrel_scale))

if (!force_pcrel && file.exists(pcrel_rds)) {
  cat("loading cached PC-Relate:", pcrel_rds, "\n")
  pcrel_h2 <- readRDS(pcrel_rds)
} else {
  geno_reader <- GdsGenotypeReader(filename = h2_gds)
  geno_data   <- GenotypeData(geno_reader)
  geno_iter   <- GenotypeBlockIterator(geno_data, snpBlock = 10000)

  pcrel_h2 <- pcrelate(
    geno_iter,
    pcs                = pcs_sub,
    scale              = pcrel_scale,   # 'variant' -> per-SNP HWE standardization
    ibd.probs          = FALSE,         # not needed for a GRM; large speedup
    sample.include     = unrel_ids,
    training.set       = unrel_ids,     # all included samples are mutually unrelated
    maf.thresh         = maf_min,       # on individual-specific AF, not cohort MAF
    maf.bound.method   = "filter",
    sample.block.size  = 5000,
    small.samp.correct = TRUE,
    BPPARAM            = BiocParallel::MulticoreParam(n_cores),
    verbose            = TRUE
  )
  saveRDS(pcrel_h2, pcrel_rds)
  close(geno_data)
}

# the cached object must still cover the samples the rest of the script assumes
stopifnot(setequal(union(pcrel_h2$kinBtwn$ID1, pcrel_h2$kinSelf$ID),
                   unrel_ids))


# --- how much data did the ISAF filter cost, and to whom? ----------------
# maf.thresh operates on individual-specific allele frequencies, so it is doing
# something different from the cohort-level --maf. In a multi-way admixed cohort
# a variant at cohort MAF 0.08 can have an ISAF near zero for individuals with a
# high proportion of ancestry X, and those cells get dropped. If the pairs
# contributing the fewest SNPs are systematically the most ancestrally distinct
# ones, that is worth reporting, and a reason to also run maf.thresh = 0.05 as a
# sensitivity analysis.

kb <- as.data.table(pcrel_h2$kinBtwn)
s  <- kb[sample(.N, min(1e6, .N))]

pc1 <- setNames(pcs_sub[, 1], rownames(pcs_sub))
s[, `:=`(d1 = abs(pc1[ID1] - pc1[ID2]),          # ancestry distance within the pair
         m1 = abs((pc1[ID1] + pc1[ID2]) / 2))]   # how extreme the pair is on PC1

cat("--- SNPs contributing per pair ---\n")
print(summary(s$nsnp))
cat("as a fraction of M:", round(min(s$nsnp) / M, 3), "at the minimum\n")
print(c(cor_with_ancestry_distance = cor(s$nsnp, s$d1),
        cor_with_ancestry_extremity = cor(s$nsnp, s$m1)))


# --- build the GRM and check it -----------------------------------------
grm_rds <- op(sprintf("grm_pcrelate_%s.rds", pcrel_scale))

if (!force_pcrel && file.exists(grm_rds)) {
  cat("loading cached GRM:", grm_rds, "\n")
  grm <- readRDS(grm_rds)
} else {
  # scaleKin = 2 puts the matrix on the GRM scale GCTA expects (diagonal ~ 1),
  # NOT the kinship scale (diagonal ~ 0.5) used internally by PC-AiR
  grm <- pcrelateToMatrix(pcrel_h2, scaleKin = 2, thresh = NULL, verbose = TRUE)
  grm <- as.matrix(grm)                 # dsyMatrix -> base symmetric
  saveRDS(grm, grm_rds)
}
stopifnot(!any(is.na(grm)), isSymmetric(unname(grm)))

# The diagonal is the check that matters here: scale='variant' plus ISAF
# instability will produce wild diagonals far more readily than negative
# eigenvalues.
cat("--- diagonal ---\n"); print(summary(diag(grm)))
cat("outside [0.8, 1.5]:", sum(diag(grm) < 0.8 | diag(grm) > 1.5), "\n")

off <- grm[upper.tri(grm)]
cat("--- off-diagonal ---\n")
print(c(mean = mean(off), sd = sd(off),
        q999 = unname(quantile(off, 0.999)), max = max(off)))

# Off-diagonal SD should be roughly 1/sqrt(Meff). A max above ~0.1 means NAToRA
# left a related pair in: go and check it rather than proceeding.

# PSD check. Cholesky rather than eigen, because it is far cheaper and answers
# the same question.
psd_ok <- !inherits(try(chol(grm), silent = TRUE), "try-error")
cat("PSD check: positive definite:", psd_ok, "\n")
if (!psd_ok) {
  lmin <- min(eigen(grm, symmetric = TRUE, only.values = TRUE)$values)
  cat("lambda_min:", lmin, " | REML breaks if |lambda_min| > sigma_e^2/sigma_g^2 (~4)\n")
}

# If this passes, stop thinking about PSD. If it fails with lambda_min well
# under 1, still proceed: gcta --reml-bendV handles it in place if the
# likelihood misbehaves, and PCGC never inverts anything.


# --- write the matrix in GCTA's gzipped text format ----------------------
# This produces <stem>.grm.gz and <stem>.grm.id, which GCTA reads with
# --grm-gz <stem>. It is a big file: n*(n+1)/2 rows.

n        <- nrow(grm)
grm_stem <- op(sprintf("grm_pcrelate_%s", pcrel_scale))
txt      <- paste0(grm_stem, ".grm")
if (file.exists(txt))                file.remove(txt)
if (file.exists(paste0(txt, ".gz"))) file.remove(paste0(txt, ".gz"))

for (start in seq(1, n, by = 500)) {
  end <- min(start + 499, n)
  ii  <- rep(start:end, times = start:end)
  jj  <- unlist(lapply(start:end, seq_len))
  fwrite(data.table(i = ii, j = jj, nsnp = M, x = grm[cbind(ii, jj)]),
         txt, sep = "\t", col.names = FALSE, append = TRUE)
}
system2("gzip", c("-f", txt))

# The .grm.id must be in the same row order as the matrix, and carry the
# FID/IID pair from the fam file. GCTA matches on the pair, and writing the IID
# twice against a fam that has real FIDs silently intersects to zero samples.
fam <- fread(paste0(h2_input, ".fam"), header = FALSE,
             col.names = c("FID","IID","PAT","MAT","SEX","PHENO"))
k <- match(rownames(grm), fam$IID)
stopifnot(!any(is.na(k)))

fwrite(data.table(FID = fam$IID[k], IID = fam$IID[k]),
       paste0(grm_stem, ".grm.id"), sep = "\t", col.names = FALSE)

cat("GRM written to:", paste0(grm_stem, ".grm.gz"), "\n")



# --- phenotype sanity check ---------------------------------------------
cv <- fread(covar_file)
cv <- cv[match(rownames(grm), cv$IID), ]
stopifnot(!any(is.na(cv$IID)))
cv[, `:=`(FID = fam$IID[k], IID = fam$IID[k])]

# 1 = control, 2 = case (PLINK/GCTA convention)
cv[, pheno := fifelse(DISEASE == 2, 2L, 1L)]
cat("cases:", sum(cv$pheno == 2), " controls:", sum(cv$pheno == 1), "\n")


# --- covariate file in LDAK format, REQUIRED for arms 4 to 7 -------------
# Arms 4 to 7 run PCGC through the LDAK program, and LDAK does not read the
# GCTA-style covariate files that make_covariate_files.r wrote in step 4.c.
# It wants a single all-numeric covariate file, with categorical variables
# expressed as 0/1 dummies, which is what this block builds.
#
# Note this is the LDAK *software*, not the LDAK heritability model: section 4
# rules out the LDAK-Thin model (alpha = -0.25), not the program.
#
# Declare any cohort-specific categorical covariate here rather than editing
# the code below. Leave it empty if you have none.
extra_covars <- character(0)        # e.g. c("PHASE") or c("PHASE", "BATCH")

cv[, sex_d := as.integer(SEX == 2)]
cv[, age_n := as.numeric(AGE)]
cat("AGE missing:", sum(is.na(cv$age_n)), "of", nrow(cv), "\n")

dummy_cols <- character(0)
for (v in extra_covars) {
  if (!v %in% names(cv))
    stop("extra covariate not found in ", covar_file, ": ", v)
  lv <- sort(unique(stats::na.omit(cv[[v]])))
  if (length(lv) != 2)
    stop(v, " has ", length(lv), " levels. This shortcut only builds 0/1 dummies ",
         "for binary variables; for more levels, build one dummy column per level yourself.")
  newcol <- paste0(tolower(v), "_d")
  set(cv, j = newcol, value = as.integer(cv[[v]] == lv[2]))
  dummy_cols <- c(dummy_cols, newcol)
}

fwrite(cv[, c("FID", "IID", "sex_d", "age_n", dummy_cols, paste0("PC", 1:n_pcs)),
          with = FALSE],
       op("pd.ldak.covar"), sep = "\t", col.names = FALSE)


# --- close the log ------------------------------------------------------
cat("finished:", format(Sys.time()), "\n")
log_close()
```

</details>

Here an example on how to run it interactively in the terminal

<details>
<summary>terminal example: running the PC-Relate GRM script</summary>

```console
(genesis) [duarte@node1 greml_pcrelate]$ nano grm_with_pcrelate.r
(genesis) [duarte@node1 greml_pcrelate]$ Rscript grm_with_pcrelate.r 
```

</details>

Given the runtime, on a cluster you will almost certainly want to submit this rather than watch it in a terminal that could drop the connection. Unlike the earlier steps there is nothing to parallelise here, since PC-Relate needs the whole genome at once, so this is a single job rather than a python script writing 22 of them.

<details>
<summary>script: <code>get_grm_pcrel.pbs</code> — SLURM job for the PC-Relate GRM</summary>

```bash
#!/bin/sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_mail@email.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=94G
#SBATCH --partition=defq
#SBATCH --job-name=get_grm_pcrel
#SBATCH -o ./logs/get_grm_pcrel.out
#SBATCH -e ./logs/get_grm_pcrel.err

# the genesis environment created in step 4. If 'conda' is not on the PATH of
# your compute nodes, replace this with the explicit path to your
# miniconda3/etc/profile.d/conda.sh
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate genesis

Rscript grm_with_pcrelate.r
```

</details>

<details>
<summary>terminal example: submitting the job</summary>

```console
(base) [duarte@node1 greml_pcrelate]$ mkdir -p ./logs
(base) [duarte@node1 greml_pcrelate]$ nano get_grm_pcrel.pbs
(base) [duarte@node1 greml_pcrelate]$ sbatch get_grm_pcrel.pbs
Submitted batch job 6118904
(base) [duarte@node1 greml_pcrelate]$ # check on it later with
(base) [duarte@node1 greml_pcrelate]$ squeue -u duarte
(base) [duarte@node1 greml_pcrelate]$ tail -20 logs/get_grm_pcrel.out
```

</details>

Two things to keep in sync when you edit the resources. First, **`--cpus-per-task` in the job header and `n_cores` in the R config block must match.** They are set independently, and asking BiocParallel for 12 workers on a 4-CPU allocation will oversubscribe the node and slow everything down rather than speed it up. Second, the memory is driven by the size of the merged fileset from Step 10.a, so if your cohort is larger than the numbers quoted above, scale `--mem` accordingly rather than assuming 94G will do.

Note that the job submits from `base`, not from an activated `genesis`: the environment is activated inside the job script, on the compute node where it actually matters.

The script prints a set of diagnostics as it goes, and they are worth reading rather than scrolling past. They are important to share with the team. 

| check | what it means | what to do if it fails |
|---|---|---|
| **diagonal spread** | Each diagonal element is an individual's self-similarity and should sit near 1. Values far outside `[0.8, 1.5]` mean the individual-specific allele frequencies are unstable for that person, usually because their ancestry is poorly captured by the ten PCs | Report the count. If it is more than a handful of people, more PCs may be needed |
| **off-diagonal max** | The most similar pair in the matrix. Above roughly 0.1 means a related pair survived NAToRA | Go back and check that pair before proceeding, do not just carry on |
| **SNPs per pair vs ancestry** | The correlation between how many SNPs a pair contributed and how ancestrally distinct that pair is. Strongly negative means the ISAF filter is systematically discarding data from the most admixed individuals | Report it, and consider rerunning with `maf.thresh = 0.05` as a sensitivity analysis |
| **PSD check** | Whether the matrix is positive definite. REML assumes it is | If it fails with `lambda_min` well under 1, proceed anyway: GCTA's `--reml-bendV` handles it, and PCGC never inverts the matrix at all (aka, doesn't really care) |

Please send us these diagnostics along with the estimates. Comparing them across cohorts is informative, since the behaviour of ancestry-aware GRMs in admixed PD cohorts is essentially uncharacterised.

The GRM itself comes out in GCTA's gzipped text format, as `outFolder_h2/grm_pcrelate_variant.grm.gz` plus its `.grm.id`. That is what arm 3 feeds to GCTA in the following steps. 

## Step 11. run GREML on the PC-Relate GRM

This completes arm 3. It is the same GCTA REML fit, the same phenotype and covariate files and the same prevalence grid as Step 9 — the only thing that changes is which GRM goes into `--grm`. Everything else being identical.

We keep this arm inside `greml_pcrelate`, next to the GRM it consumes, so that the two arms never accidentally read each other's files.

**One conversion first.** The R script in Step 10.b wrote the GRM in GCTA's gzipped *text* format (`grm_pcrelate_variant.grm.gz` plus `.grm.id`), because that is the format an R matrix can be streamed into. GCTA reads that with `--grm-gz`, but it parses the whole text file on every run, and we are about to run it three times. So the script below converts it once to GCTA's binary format and then works from that. The conversion is skipped if the binary GRM is already there, so re-runs are cheap.

<details>
<summary>script: <code>run_greml_pcrelate.sh</code> — converts the GRM to binary and runs GREML across the prevalence grid</summary>

```bash
#!/bin/bash
#
# Single-GRM GREML on the PC-Relate variant GRM, across a prevalence grid.
#
# Run it inside /working_directory_h2/greml_pcrelate:
#     bash run_greml_pcrelate.sh
#
# Paths are anchored to this script's own directory, so it works from any
# CWD. Override GRM_TEXT, PCA_DIR, GCTA, THREADS or OUTDIR from the
# environment:
#     THREADS=12 bash run_greml_pcrelate.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The GCTA executable downloaded in step 8.a, kept in ./programs like plink
GCTA="${GCTA:-${SCRIPT_DIR}/../programs/gcta-1.95.3-linux-x86_64/gcta}"
THREADS="${THREADS:-4}"

# The gzipped text GRM written by grm_with_pcrelate.r in step 10.b, and the
# binary version of it that this script derives once and then reuses.
GRM_TEXT="${GRM_TEXT:-${SCRIPT_DIR}/outFolder_h2/grm_pcrelate_variant}"
GRM="${GRM_TEXT}_bin"

# The covariate files written by make_covariate_files.r in step 4.c
PCA_DIR="${PCA_DIR:-${SCRIPT_DIR}/../pca_and_such/outFolder_pca_andSuch}"

OUTDIR="${OUTDIR:-${SCRIPT_DIR}/greml_out}"

PHENO="${PCA_DIR}/covariate_file_no_related_pairs_pheno.tsv"
QCOVAR="${PCA_DIR}/covariate_file_no_related_pairs_qcovar.tsv"
COVAR="${PCA_DIR}/covariate_file_no_related_pairs_covar.tsv"

PREVALENCES=(0.005 0.01 0.02)

STEM="${OUTDIR}/greml_pcrel_prev"
LOG="${OUTDIR}/greml_pcrel.log"
SUMMARY="${OUTDIR}/greml_pcrel_prevalence_grid.tsv"

mkdir -p "${OUTDIR}"

# Duplicate all stdout/stderr into the log. Overwrites on each run.
exec > >(tee "${LOG}") 2>&1

echo "=============================================================="
echo "Started : $(date)"
echo "Host    : $(hostname)"
echo "GRM     : ${GRM}"
echo "Threads : ${THREADS}"
echo "Log     : ${LOG}"
echo "=============================================================="

if [[ ! -x "${GCTA}" ]]; then
    echo "ERROR: GCTA binary not found or not executable: ${GCTA}" >&2
    exit 1
fi

# ---------------------------------------------------------------------
# Convert the text GRM to binary, once
# ---------------------------------------------------------------------
if [[ -s "${GRM}.grm.bin" && -s "${GRM}.grm.id" ]]; then
    echo "Binary GRM already present, skipping conversion: ${GRM}"
else
    for f in "${GRM_TEXT}.grm.gz" "${GRM_TEXT}.grm.id"; do
        if [[ ! -s "${f}" ]]; then
            echo "ERROR: missing or empty input: ${f}" >&2
            echo "       did grm_with_pcrelate.r finish in step 10.b?" >&2
            exit 1
        fi
    done

    echo "Converting the text GRM to binary -> ${GRM}"
    "${GCTA}" \
        --grm-gz "${GRM_TEXT}" \
        --make-grm \
        --thread-num "${THREADS}" \
        --out "${GRM}"
fi

# ---------------------------------------------------------------------
# Validate every input before spending time on REML
# ---------------------------------------------------------------------
for f in "${GRM}.grm.bin" "${GRM}.grm.id" "${PHENO}" "${QCOVAR}" "${COVAR}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing or empty input: ${f}" >&2
        exit 1
    fi
done

echo "Samples in GRM: $(wc -l < "${GRM}.grm.id")"

# Arms 2 and 3 must run on exactly the same people, otherwise the comparison
# between them is not telling us about GRM construction any more. Warn rather
# than fail, since you may legitimately be running arm 3 before arm 2.
ARM2_ID="${SCRIPT_DIR}/../greml/gcta_singleGRM/genomewide/genomewide_single_grm.grm.id"
if [[ -s "${ARM2_ID}" ]]; then
    if cmp -s <(sort "${ARM2_ID}") <(sort "${GRM}.grm.id"); then
        echo "Sample set matches the standard GCTA GRM from arm 2."
    else
        echo "WARNING: the sample set differs from the arm 2 GRM:" >&2
        echo "         arm 2: $(wc -l < "${ARM2_ID}") samples, ${ARM2_ID}" >&2
        echo "         arm 3: $(wc -l < "${GRM}.grm.id") samples, ${GRM}.grm.id" >&2
        echo "         the two arms are no longer directly comparable." >&2
    fi
else
    echo "NOTE: arm 2 GRM not found, skipping the sample-set comparison."
fi

# ---------------------------------------------------------------------
# Run GREML across the prevalence grid
# ---------------------------------------------------------------------
FAILED=()

for K in "${PREVALENCES[@]}"; do
    OUT="${STEM}_${K}"

    echo
    echo "########## GREML | K = ${K} ##########"

    # Do not let one non-converging K abort the rest of the grid.
    if ! "${GCTA}" --reml \
        --grm "${GRM}" \
        --pheno "${PHENO}" \
        --qcovar "${QCOVAR}" \
        --covar "${COVAR}" \
        --prevalence "${K}" \
        --thread-num "${THREADS}" \
        --out "${OUT}"; then
        echo "ERROR: GCTA exited nonzero for K=${K}" >&2
        FAILED+=("${K}")
        continue
    fi

    if [[ ! -s "${OUT}.hsq" ]]; then
        echo "ERROR: no .hsq written for K=${K}" >&2
        FAILED+=("${K}")
        continue
    fi

    echo
    cat "${OUT}.hsq"
done

# ---------------------------------------------------------------------
# Collate the grid into a single summary table
# ---------------------------------------------------------------------
echo
echo "########## Prevalence grid summary ##########"

printf 'K\tsource\testimate\tSE\n' > "${SUMMARY}"

for K in "${PREVALENCES[@]}"; do
    HSQ="${STEM}_${K}.hsq"
    [[ -s "${HSQ}" ]] || continue
    awk -v k="${K}" '
        $1 == "V(G)/Vp"   { printf "%s\t%s\t%s\t%s\n", k, $1, $2, $3 }
        $1 == "V(G)/Vp_L" { printf "%s\t%s\t%s\t%s\n", k, $1, $2, $3 }
        $1 == "Pval"      { printf "%s\t%s\t%s\tNA\n",  k, $1, $2 }
        $1 == "n"         { printf "%s\t%s\t%s\tNA\n",  k, $1, $2 }
    ' "${HSQ}" >> "${SUMMARY}"
done

column -t "${SUMMARY}" || cat "${SUMMARY}"

echo
echo "Summary: ${SUMMARY}"

if [[ "${#FAILED[@]}" -gt 0 ]]; then
    echo "ERROR: GREML failed for K = ${FAILED[*]}" >&2
    exit 1
fi

echo "Finished: $(date)"

# Give the tee subprocess a moment to flush before the shell exits.
sleep 1
```

</details>

In the terminal do

<details>
<summary>terminal example: saving and running GREML on the PC-Relate GRM</summary>

```console
(genesis) [duarte@node1 greml_pcrelate]$ conda deactivate
(base) [duarte@node1 greml_pcrelate]$ nano run_greml_pcrelate.sh
(base) [duarte@node1 greml_pcrelate]$ bash run_greml_pcrelate.sh
==============================================================
Host    : node1.lerner.ccf.org
GRM     : /home/duarte/working_directory_h2/greml_pcrelate/outFolder_h2/grm_pcrelate_variant_bin
Threads : 4
Log     : /home/duarte/working_directory_h2/greml_pcrelate/greml_out/greml_pcrel.log
==============================================================
Converting the text GRM to binary -> /home/duarte/working_directory_h2/greml_pcrelate/outFolder_h2/grm_pcrelate_variant_bin
Samples in GRM: 6412
Sample set matches the standard GCTA GRM from arm 2.

########## GREML | K = 0.005 ##########
...
```

</details>

Note that the `genesis` environment is not needed here, only GCTA, so you can deactivate it. And on a second run the conversion step is skipped:

<details>
<summary>terminal example: the output folder and the collated grid</summary>

```console
(base) [duarte@node1 greml_pcrelate]$ ls greml_out/
greml_pcrel.log                       greml_pcrel_prev_0.01.hsq  greml_pcrel_prev_0.02.log
greml_pcrel_prevalence_grid.tsv       greml_pcrel_prev_0.01.log  greml_pcrel_prev_0.005.hsq
greml_pcrel_prev_0.02.hsq             greml_pcrel_prev_0.005.log
```

</details>

Read the `.hsq` exactly as described at the end of Step 9: `n` should match your unrelated sample count, $$V(G)/Vp$$ is identical across the three runs by construction, and $$V(G)/Vp_L$$ is the number we report.

What is new here is the comparison. Put the two grids side by side at the primary prevalence:

<details>
<summary>terminal example: arm 2 against arm 3 at K = 0.005</summary>

```console
(base) [duarte@node1 greml_pcrelate]$ grep 'V(G)/Vp_L' \
>   ../greml/gcta_singleGRM/greml/PD_singleGRM_GREML_K0.005.hsq \
>   greml_out/greml_pcrel_prev_0.005.hsq
../greml/gcta_singleGRM/greml/PD_singleGRM_GREML_K0.005.hsq:V(G)/Vp_L	0.286547	0.029183
greml_out/greml_pcrel_prev_0.005.hsq:V(G)/Vp_L	0.261584	0.027979
```

</details>

Please share both prevalence grids and both sets of `.hsq` files.

## Step 12. run PCGC on both GRMs, with a permutation null

This one step closes out four of the seven arms at once, because PCGC is cheap enough to run many times:

| arm | GRM | what it is |
|---|---|---|
| 4 | standard GCTA | PCGC observed estimate. Isolates ascertainment handling: same GRM as arm 2, only the estimator differs |
| 5 | PC-Relate | PCGC observed estimate. Both corrections applied together |
| 7 | standard GCTA, shuffled | Null control for arm 4 |
| 6 | PC-Relate, shuffled | Null control for arm 5 |

The script below does the observed run and the permutation replicates for **one** GRM. You run it twice, once per GRM, and the arm labels above fall out.

Why bother with the permutations at all? GREML gives us a likelihood ratio test and gives us a p-value for the observed estimate compared to a null, but PCGC is a moment-based estimator, and through LDAK we can shuffle the correspondence between phenotypes and kinships destroying any real genetic signal while leaving the sample size, the ascertainment and the ancestry structure exactly as they are. The spread of the resulting estimates *is* the noise floor. 

### Step 12.a get LDAK

PCGC is implemented in LDAK, so we need that program. As with plink and GCTA, it goes in `./programs`.

Remember the distinction made in Step 10.b: we are using the LDAK *software* for its PCGC implementation, not the LDAK heritability *model*. Every run below stays on the standard `α = −1` model, on our own GRMs.

<details>
<summary>terminal example: downloading LDAK and checking the executable</summary>

```console
(base) [duarte@node1 working_directory_h2]$ cd programs/
(base) [duarte@node1 programs]$ wget https://github.com/dougspeed/LDAK/raw/main/ldak6.3.linux  
(base) [duarte@node1 programs]$ chmod a+x ldak6.3.linux 
```

</details>

### Step 12.b create the folder

PCGC consumes both GRMs, one from `greml` and one from `greml_pcrelate`, so it gets its own folder at the same level rather than living inside either one.

<details>
<summary>terminal example: creating the pcgc folder</summary>

```console
(base) [duarte@node1 programs]$ cd ..
(base) [duarte@node1 working_directory_h2]$ mkdir pcgc
(base) [duarte@node1 working_directory_h2]$ cd pcgc/
(base) [duarte@node1 pcgc]$ ls ..
covldsc_h2  genetic_data  greml  greml_pcrelate  gwas  pca_and_such  pcgc  programs  reference_panel
```

</details>

### Step 12.c run PCGC and the permutation null

The script takes one argument, `gcta` or `pcrelate`, which selects the GRM and labels every output accordingly. Everything else is shared.

Three things are worth knowing before you launch it.

**The GRM has to be covariate-adjusted first.** LDAK's `--adjust-grm` projects the covariates out of the kinship matrix, and PCGC needs that done when ancestry PCs are among the covariates. It depends only on the covariates, not on the prevalence, so it runs once and is reused across the whole grid.

**Unlike GREML, the prevalence enters the estimation itself.** In Step 9 you saw $$V(G)/Vp$$ come out identical across the three $$K$$ values, because GREML fits on the observed scale and only transforms afterwards. PCGC does not work that way: the ascertainment correction uses $$K$$ inside the moment equations, so all three grids differ from scratch.

**Start small.** The defaults are 100 replicates at each of three prevalences, so 303 LDAK runs per GRM, 606 in total. Do a smoke test first with `NPERM=10` to confirm everything parses if you want, then launch the full thing.

<details>
<summary>script: <code>run_pcgc.sh</code> — PCGC observed estimate plus a permutation null, for one GRM</summary>

```bash
#!/bin/bash
#
# LDAK PCGC on a covariate-adjusted GRM, across a prevalence grid, with a
# permutation null at each K.
#
# Run it inside /working_directory_h2/pcgc, once per GRM:
#     bash run_pcgc.sh gcta        # arms 4 and 7
#     bash run_pcgc.sh pcrelate    # arms 5 and 6
#
# Structure:
#   Step 1  --adjust-grm, once (depends on covariates only, not on K)
#   Step 2  per K: one observed run (rep 0), then NPERM runs with
#           --permute YES, which breaks the correspondence between the
#           phenotypes and the kinship matrix and so estimates h2 under
#           the null of no genetic signal
#   Step 3  per K: mean and SD of the permuted estimates, the mean of
#           LDAK's reported SEs, and the calibration ratio
#
# Paths anchor to this script's directory, so it works from any CWD.
# Override LDAK, GRM_IN, ADJ_GRM, COVAR, PCA_DIR, OUTDIR, NPERM from the
# environment, e.g.
#     NPERM=10 bash run_pcgc.sh gcta              # quick smoke test
#     NPERM=100 OUTDIR=/scratch/pcgc bash run_pcgc.sh pcrelate

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------
# Which GRM: 'gcta' (arms 4 and 7) or 'pcrelate' (arms 5 and 6)
# ---------------------------------------------------------------------
WHICH_GRM="${1:-}"

case "${WHICH_GRM}" in
    gcta)
        DEFAULT_GRM_IN="${SCRIPT_DIR}/../greml/gcta_singleGRM/genomewide/genomewide_single_grm"
        ;;
    pcrelate)
        DEFAULT_GRM_IN="${SCRIPT_DIR}/../greml_pcrelate/outFolder_h2/grm_pcrelate_variant_bin"
        ;;
    *)
        echo "Usage: bash run_pcgc.sh {gcta|pcrelate}" >&2
        echo "  gcta     -> arms 4 and 7, the standard GCTA GRM from step 9" >&2
        echo "  pcrelate -> arms 5 and 6, the PC-Relate GRM from step 10.b" >&2
        exit 1
        ;;
esac

# The LDAK executable downloaded in step 12.a, kept in ./programs
LDAK="${LDAK:-${SCRIPT_DIR}/../programs/ldak6.3.linux}"

# Input GRM and the covariate-adjusted GRM it produces. --adjust-grm is
# mandatory before PCGC when ancestry PCs are in the covariate file.
GRM_IN="${GRM_IN:-${DEFAULT_GRM_IN}}"

OUTDIR="${OUTDIR:-${SCRIPT_DIR}/out_${WHICH_GRM}}"
KLOG_DIR="${OUTDIR}/pcgc_screen_logs"

ADJ_GRM="${ADJ_GRM:-${OUTDIR}/grm_${WHICH_GRM}_adj}"

# The all-numeric covariate file written by grm_with_pcrelate.r in step 10.b.
# Both arms use the same covariates: that is the point of the comparison.
COVAR="${COVAR:-${SCRIPT_DIR}/../greml_pcrelate/outFolder_h2/pd.ldak.covar}"

# The phenotype file written by make_covariate_files.r in step 4.c
PCA_DIR="${PCA_DIR:-${SCRIPT_DIR}/../pca_and_such/outFolder_pca_andSuch}"
PHENO="${PHENO:-${PCA_DIR}/covariate_file_no_related_pairs_pheno.tsv}"

PREVALENCES=(0.005 0.01 0.02)
NPERM="${NPERM:-100}"

STEM="${OUTDIR}/pcgc_grm_${WHICH_GRM}_adj"
LOG="${OUTDIR}/pcgc_ldak_${WHICH_GRM}.log"
RESULTS="${OUTDIR}/pcgc_all_estimates_${WHICH_GRM}.tsv"   # long format, one row per run
SUMMARY="${OUTDIR}/pcgc_null_summary_${WHICH_GRM}.tsv"    # one row per K

mkdir -p "${OUTDIR}" "${KLOG_DIR}"

# Duplicate all stdout/stderr into the master log. Overwrites on each run.
exec > >(tee "${LOG}") 2>&1

echo "=============================================================="
echo "Started    : $(date)"
echo "Host       : $(hostname)"
echo "GRM choice : ${WHICH_GRM}"
echo "Input GRM  : ${GRM_IN}"
echo "Adj GRM    : ${ADJ_GRM}"
echo "Covar      : ${COVAR}"
echo "Pheno      : ${PHENO}"
echo "Prevalences: ${PREVALENCES[*]}"
echo "Permutations per K: ${NPERM}"
echo "Log        : ${LOG}"
echo "=============================================================="

if [[ ! -x "${LDAK}" ]]; then
    echo "ERROR: LDAK binary not found or not executable: ${LDAK}" >&2
    exit 1
fi

for f in "${GRM_IN}.grm.bin" "${GRM_IN}.grm.id" "${COVAR}" "${PHENO}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing or empty input: ${f}" >&2
        exit 1
    fi
done

# Arms 4 and 5 must run on exactly the same people as each other and as the
# GREML arms, otherwise the comparison stops being about the estimator.
echo "Samples in input GRM: $(wc -l < "${GRM_IN}.grm.id")"
echo "Samples in phenotype file: $(wc -l < "${PHENO}")"


# ---------------------------------------------------------------------
# Step 1: covariate-adjust the GRM (once)
# ---------------------------------------------------------------------
echo
echo "########## --adjust-grm ##########"

"${LDAK}" --adjust-grm "${ADJ_GRM}" \
    --grm "${GRM_IN}" \
    --covar "${COVAR}" \
    --kinship-details NO

if [[ ! -s "${ADJ_GRM}.grm.bin" ]]; then
    echo "ERROR: adjusted GRM was not written: ${ADJ_GRM}.grm.bin" >&2
    exit 1
fi

echo "Adjusted GRM samples: $(wc -l < "${ADJ_GRM}.grm.id")"


# ---------------------------------------------------------------------
# Helper: parse "Total heritability 0.1651 (SE 0.0417)" from a screen log
# Prints "<estimate>\t<SE>" or nothing if the line is absent.
# ---------------------------------------------------------------------
parse_h2 () {
    grep -m1 'Total heritability' "$1" | awk '
        {
            se = $5
            gsub(/[()]/, "", se)
            printf "%s\t%s\n", $3, se
        }
    '
}


# ---------------------------------------------------------------------
# Step 2: observed run + permutation replicates, per K
# ---------------------------------------------------------------------
printf 'K\trep\tmode\th2_liability\tSE_analytic\n' > "${RESULTS}"

TMP_SCREEN="$(mktemp)"
trap 'rm -f "${TMP_SCREEN}"' EXIT

FAILED=()

for K in "${PREVALENCES[@]}"; do

    echo
    echo "=============================================================="
    echo "K = ${K}"
    echo "=============================================================="

    # ---- observed (rep 0) -------------------------------------------
    OBS_OUT="${STEM}_prev_${K}"
    OBS_SCREEN="${KLOG_DIR}/pcgc_K${K}_observed.screen"

    echo "--- observed run ---"

    # set -o pipefail is what makes this catch a failing LDAK rather than
    # a successful tee
    if ! "${LDAK}" --pcgc "${OBS_OUT}" \
        --grm "${ADJ_GRM}" \
        --pheno "${PHENO}" \
        --covar "${COVAR}" \
        --prevalence "${K}" \
        --kinship-details NO 2>&1 | tee "${OBS_SCREEN}"; then
        echo "ERROR: LDAK exited nonzero for K=${K} (observed)" >&2
        FAILED+=("${K}:observed")
        continue
    fi

    OBS_PARSED="$(parse_h2 "${OBS_SCREEN}" || true)"
    if [[ -z "${OBS_PARSED}" ]]; then
        echo "ERROR: no 'Total heritability' line for K=${K} (observed)" >&2
        FAILED+=("${K}:observed")
        continue
    fi
    printf '%s\t0\tobserved\t%s\n' "${K}" "${OBS_PARSED}" >> "${RESULTS}"

    # ---- permutation replicates -------------------------------------
    # All replicates share one LDAK output prefix (the .pcgc files are
    # overwritten each time; only the parsed estimate is retained). Full
    # screen output is concatenated into one per-K file for auditing,
    # rather than flooding the master log with NPERM x LDAK output.
    PERM_OUT="${STEM}_prev_${K}_perm"
    PERM_SCREEN="${KLOG_DIR}/pcgc_K${K}_perm_all.screen"
    : > "${PERM_SCREEN}"

    echo "--- ${NPERM} permutation replicates ---"

    PERM_FAILS=0

    for REP in $(seq 1 "${NPERM}"); do

        if ! "${LDAK}" --pcgc "${PERM_OUT}" \
            --grm "${ADJ_GRM}" \
            --pheno "${PHENO}" \
            --covar "${COVAR}" \
            --prevalence "${K}" \
            --permute YES \
            --kinship-details NO > "${TMP_SCREEN}" 2>&1; then
            echo "WARNING: LDAK exited nonzero for K=${K} rep=${REP}" >&2
            PERM_FAILS=$((PERM_FAILS + 1))
            {
                echo "===== K=${K} rep=${REP} (NONZERO EXIT) ====="
                cat "${TMP_SCREEN}"
            } >> "${PERM_SCREEN}"
            continue
        fi

        {
            echo "===== K=${K} rep=${REP} ====="
            cat "${TMP_SCREEN}"
        } >> "${PERM_SCREEN}"

        PARSED="$(parse_h2 "${TMP_SCREEN}" || true)"
        if [[ -z "${PARSED}" ]]; then
            echo "WARNING: no 'Total heritability' line for K=${K} rep=${REP}" >&2
            PERM_FAILS=$((PERM_FAILS + 1))
            continue
        fi

        printf '%s\t%s\tpermuted\t%s\n' "${K}" "${REP}" "${PARSED}" >> "${RESULTS}"

        if (( REP % 10 == 0 )); then
            echo "  ... ${REP}/${NPERM} replicates done"
        fi
    done

    if (( PERM_FAILS > 0 )); then
        echo "WARNING: ${PERM_FAILS}/${NPERM} permutation replicates failed at K=${K}" >&2
    fi
done


# ---------------------------------------------------------------------
# Step 3: per-K null summary
# ---------------------------------------------------------------------
echo
echo "########## Per-K null summary ##########"

awk -F'\t' '
    NR == 1 { next }

    $3 == "observed" {
        obs[$1]    = $4
        obs_se[$1] = $5
        next
    }

    $3 == "permuted" {
        k = $1; h = $4 + 0; s = $5 + 0
        n[k]++
        sum[k]    += h
        sumsq[k]  += h * h
        sumse[k]  += s
        if (h < 0) neg[k]++
        if (k in obs && h >= obs[k] + 0) ge[k]++
        if (!(k in seen)) { seen[k] = 1; order[++nk] = k }
    }

    END {
        printf "K\th2_obs\tSE_analytic_obs\tn_perm\tmean_perm\tSD_perm\tmean_SE_analytic\tratio_SD_over_SE\tfrac_negative\tperm_p\n"

        for (i = 1; i <= nk; i++) {
            k = order[i]
            N = n[k]
            if (N < 2) continue

            mean   = sum[k] / N
            var    = (sumsq[k] - N * mean * mean) / (N - 1)
            sd     = (var > 0) ? sqrt(var) : 0
            meanse = sumse[k] / N
            ratio  = (meanse > 0) ? sprintf("%.4f", sd / meanse) : "NA"
            fneg   = (k in neg) ? neg[k] / N : 0
            g      = (k in ge) ? ge[k] : 0
            p      = (k in obs) ? sprintf("%.4f", (1 + g) / (N + 1)) : "NA"
            o      = (k in obs) ? obs[k] : "NA"
            ose    = (k in obs) ? obs_se[k] : "NA"

            printf "%s\t%s\t%s\t%d\t%.6f\t%.6f\t%.6f\t%s\t%.2f\t%s\n",
                   k, o, ose, N, mean, sd, meanse, ratio, fneg, p
        }
    }
' "${RESULTS}" > "${SUMMARY}"

column -t "${SUMMARY}" || cat "${SUMMARY}"

echo
echo "All estimates (long format): ${RESULTS}"
echo "Null summary              : ${SUMMARY}"
echo "Per-K screen logs         : ${KLOG_DIR}"
echo
echo "Reading the summary:"
echo "  mean_perm        should sit near 0; a systematic offset is estimator bias"
echo "  SD_perm          the empirical null SE -- the noise floor for h2_obs"
echo "  ratio_SD_over_SE ~1 means LDAK's analytic SEs are calibrated here;"
echo "                   >1 means they are anticonservative"
echo "  frac_negative    PCGC is unconstrained, so expect ~0.5 under the null."
echo "                   Near 0 suggests estimates are being floored, which"
echo "                   invalidates mean_perm and SD_perm."
echo "  perm_p           (1 + #{h2_perm >= h2_obs}) / (n_perm + 1)"

if [[ "${#FAILED[@]}" -gt 0 ]]; then
    echo
    echo "ERROR: observed run failed for: ${FAILED[*]}" >&2
    exit 1
fi

echo
echo "Finished: $(date)"

# Give the tee subprocess a moment to flush before the shell exits.
sleep 1
```

</details>

In the terminal do

<details>
<summary>terminal example: smoke test, then the two full runs</summary>

```console
(base) [duarte@node1 pcgc]$ nano run_pcgc.sh
(base) [duarte@node1 pcgc]$ # you can smoke test first: 10 replicates instead of 100
(base) [duarte@node1 pcgc]$ #NPERM=10 bash run_pcgc.sh gcta
(base) [duarte@node1 pcgc]$ # if looks right, now the real thing for both GRMs
(base) [duarte@node1 pcgc]$ bash run_pcgc.sh gcta
(base) [duarte@node1 pcgc]$ bash run_pcgc.sh pcrelate
```

</details>

Each run writes into its own `out_gcta` or `out_pcrelate` folder, so the two never collide.

The summary table is the deliverable here, and each column answers a specific question:

| column | what it tells you | what "good" looks like |
|---|---|---|
| `mean_perm` | Is PCGC biased in this cohort? Shuffling removes all genetic signal, so the average estimate should be zero | Near 0. A systematic offset is estimator bias, and it should be subtracted from your reading of `h2_obs` |
| `SD_perm` | The empirical noise floor. This is the honest standard error, measured rather than derived | Compare `h2_obs` against this, not against the analytic SE |
| `ratio_SD_over_SE` | Are LDAK's reported SEs trustworthy here? | Near 1. Above 1 means the analytic SEs are anticonservative, and confidence intervals built from them are too narrow |
| `frac_negative` | A sanity check on the null itself. PCGC is unconstrained, so under the null half the estimates should fall below zero | Near 0.5. Near 0 means estimates are being floored somewhere, which invalidates both `mean_perm` and `SD_perm` |
| `perm_p` | A non-parametric test of the observed estimate against its own null | With 100 replicates the smallest attainable value is 1/101 ≈ 0.0099 |

Note that last point: `perm_p` cannot go below `1/(NPERM+1)`, so a value of 0.0099 with 100 replicates means "no permutation exceeded the observed estimate", not "p = 0.0099 exactly". If you need to resolve smaller p-values, raise `NPERM`.

Finally, comparing arms 4 and 5 against arms 2 and 3 is where this gets interesting. Arm 4 against arm 2 holds the GRM fixed and swaps GREML for PCGC, so any difference is about how the two estimators handle case ascertainment. Arm 5 against arm 3 asks the same question on the ancestry-aware GRM. If PCGC comes out systematically lower than GREML in your cohort, that is the ascertainment correction doing its job, and it is one of the headline comparisons this project set out to make.

Please share both null summary tables, both long-format estimate files, and the master logs. Step 13 collects them for you along with everything else, so there is no need to go hunting for them by hand.

## Step 13. collect the shareables into one folder

Throughout the tutorial we have asked you to send us particular files. Rather than making you collect all of that by hand across nine folders and 22 chromosome subfolders, this last script walks the project tree and gathers everything into a single folder.

It stops there, on purpose. Nothing is compressed and nothing is sent. You get a plain folder you can open and read through, and when you are satisfied with what is in it, you zip it yourself and send it over.

**It is built around one rule: nothing that identifies a participant leaves your institution.** The bundle carries logs, heritability estimates, variance components, prevalence grids, diagnostics and plots. It does not carry covariate or phenotype files, GRMs or their `.grm.id` files, kinship tables, sample ID lists, or plink filesets. After copying, the script re-scans the finished bundle and deletes anything matching an individual-level pattern that slipped through, so the exclusion does not depend on the copy rules alone being right.

Here is what it collects and where each piece comes from:

| What | Where it comes from | Why we need it |
|---|---|---|
| Sample counts, cases, controls, disease proportion | Step 4.c log | The prevalence every liability-scale number in the bundle was computed with |
| PC scatter plots, both rounds | Step 4.a | To assess ancestry structure and whether the PCs behave |
| SNP–PC correlation plots, plus the top 1000 rows per PC | Step 4.a | To spot a PC driven by a single locus rather than by ancestry |
| plink and cov-LDSC logs, per chromosome | Steps 2, 5.b | Variant counts at each filtering stage |
| GWAS logs, Manhattan and QQ plots | Steps 6, 6.a | Telling an underpowered GWAS apart from uncorrected stratification |
| cov-LDSC estimates, three prevalences | Step 7 | **Arm 1** |
| GREML `.hsq` files and prevalence grid | Step 9 | **Arm 2** |
| PC-Relate GRM diagnostics log | Step 10.b | Diagonal, off-diagonal, other diagnostics |
| GREML `.hsq` files and prevalence grid | Step 11 | **Arm 3** |
| PCGC null summaries, long-format estimates, logs | Step 12 | **Arms 4 to 7** |

Two things are deliberately left out for size rather than privacy: the GWAS summary statistics and the 22 LD score files. Both are aggregate and perfectly shareable, they are just very large. If we need them for your cohort we will ask.

We are also not collecting log files that contain individual-level data. The one that matters here is the PCA pipeline log, `pipeline_run.log`, which prints a corner of the KING and PC-Relate matrices and the head of the kinship pair table, all of it labelled with sample IDs. It is skipped, and there is nothing for you to redact. Every other log in the pipeline reports counts, variances and estimates only, and those all travel.

<details>
<summary>script: <code>collect_shareables.sh</code> — gathers every shareable output into one folder</summary>

```bash
#!/bin/bash
# Collect every shareable result of the h2 tutorial into a single folder,
# ready to be zipped and sent back to the coordinating team.
#
# Run it inside /working_directory_h2:
#     bash collect_shareables.sh <COHORT_NAME>
#
# e.g.
#     bash collect_shareables.sh LARGE-PD
#
# WHAT GOES IN
#   Logs, heritability estimates (.hsq, LDSC .log, PCGC summaries), prevalence
#   grid tables, GRM diagnostics, and the diagnostic plots. All of it is either
#   aggregate (counts, variances, estimates) or variant-level.
#
# WHAT STAYS OUT, deliberately
#   Nothing that carries a sample identifier or a per-sample value:
#     - covariate and phenotype files (*_covar.tsv, *_qcovar.tsv, *_pheno.tsv,
#       *_PCs.tsv, *_IIDs.txt, unrelated_IIDs.txt, pcair_r2_covariates*.tsv,
#       pd.ldak.covar)
#     - GRMs and their id files (*.grm.bin, *.grm.N.bin, *.grm.gz, *.grm.id)
#     - pairwise kinship tables and the NAToRA removal lists
#     - plink filesets (*.bed/*.bim/*.fam/*.pgen/*.psam/*.pvar) and *.gds/*_rds
#     - log files that carry individual-level data. The PCA pipeline log
#       (pca_and_such/outFolder_pca_andSuch/pipeline_run.log) prints a corner
#       of the KING and PC-Relate matrices and the head of the kinship pair
#       table, all of it labelled with sample IDs, so it is not collected.
#   Also left out for size rather than privacy: the GWAS summary statistics and
#   the 22 cov-LDSC LD score sets. They are aggregate and can be sent on
#   request:
#     gwas/out/gwaslab_output/gwas_sumstats_allchr.ldsc.tsv.gz
#     reference_panel/chr*/out/covldsc_chr*.l2.ldscore.gz  (plus .l2.M*)
#
# OUTPUT
#   A plain folder, not an archive. Look through it, and once you are happy
#   with what is in there, zip it yourself and send it over.

set -euo pipefail
shopt -s nullglob          # a glob that matches nothing disappears instead of
                           # being passed through literally

# ---------------------------------------------------------------------
# Arguments and paths
# ---------------------------------------------------------------------

if [[ $# -lt 1 ]]; then
    echo "usage: bash collect_shareables.sh <COHORT_NAME>" >&2
    echo "   e.g. bash collect_shareables.sh LARGE-PD" >&2
    exit 1
fi

COHORT="$1"
# Keep the label filesystem-safe: anything that is not a letter, digit, dot,
# underscore or dash becomes a dash.
COHORT_SAFE="$(echo "${COHORT}" | tr -c 'A-Za-z0-9._-' '-' | sed 's/-\{2,\}/-/g; s/^-//; s/-$//')"

# Anchor to the script's own location so it works from anywhere.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${WORKDIR:-${SCRIPT_DIR}}"
cd "${WORKDIR}"

STAMP="$(date +%Y%m%d)"
BUNDLE_NAME="h2_shareables_${COHORT_SAFE}_${STAMP}"
BUNDLE="${WORKDIR}/${BUNDLE_NAME}"

if [[ -e "${BUNDLE}" ]]; then
    echo "ERROR: ${BUNDLE} already exists." >&2
    echo "       Remove it or rename it before running again, so nothing" >&2
    echo "       from an older run is silently mixed into the new bundle." >&2
    exit 1
fi

mkdir -p "${BUNDLE}"

MANIFEST="${BUNDLE}/00_MANIFEST.tsv"
MISSING="${BUNDLE}/00_MISSING.txt"
README="${BUNDLE}/00_README.txt"

: > "${MISSING}"

echo "=============================================================="
echo "Cohort   : ${COHORT}"
echo "Workdir  : ${WORKDIR}"
echo "Bundle   : ${BUNDLE}"
echo "Started  : $(date)"
echo "=============================================================="
echo

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

# grab <critical: yes|no> <description> <file...>
#
# Copies each file into the bundle, preserving its path relative to the working
# directory, so nothing collides across the 22 chromosome folders and the
# result mirrors the project tree the tutorial describes.
#
# When a glob matches nothing and the item was marked critical, the description
# is recorded in 00_MISSING.txt rather than failing the run: a cohort that
# stopped after arm 2 should still be able to send what it has.
grab() {
    local critical="$1" desc="$2"
    shift 2
    local n=0 f rel

    for f in "$@"; do
        [[ -f "${f}" ]] || continue
        rel="${f#./}"
        mkdir -p "${BUNDLE}/$(dirname "${rel}")"
        cp -p "${f}" "${BUNDLE}/${rel}"
        n=$((n + 1))
    done

    if [[ ${n} -eq 0 ]]; then
        if [[ "${critical}" == "yes" ]]; then
            printf 'MISSING   %s\n' "${desc}" >> "${MISSING}"
            printf '  !! %-55s no files found\n' "${desc}"
        else
            printf '  -- %-55s none\n' "${desc}"
        fi
    else
        printf '  ok %-55s %3d file(s)\n' "${desc}" "${n}"
    fi
}

# ---------------------------------------------------------------------
# Step 4 — PCA, relatedness, and the sample summary
# ---------------------------------------------------------------------
echo "[Step 4] PCA and sample summary"

PCA_DIR="pca_and_such/outFolder_pca_andSuch"

# The sample-size summary: total N, cases, controls, disease proportion.
# Counts only, no identifiers. This is the file the prevalence used by LDSC and
# every liability-scale transformation is read from.
grab yes "sample count / prevalence summary" \
    "${PCA_DIR}"/covariate_file_no_related_pairs.log

# PC-pair scatter plots, both rounds of PC-AiR.
grab yes "PC-pair scatter plots (round 1)" \
    "${PCA_DIR}"/plot_PC[0-9][0-9]_PC[0-9][0-9].png
grab yes "PC-pair scatter plots (round 2)" \
    "${PCA_DIR}"/plot_PC[0-9][0-9]_PC[0-9][0-9]_2ndround.png

# SNP-PC correlation plots. These are how a PC driven by a single locus
# (an inversion, an unremoved high-LD region) is spotted.
grab yes "SNP-PC correlation plots (round 1)" \
    "${PCA_DIR}"/snpcorr_1stRound_PC[0-9][0-9].png
grab yes "SNP-PC correlation plots (round 2)" \
    "${PCA_DIR}"/snpcorr_2ndRound_PC[0-9][0-9].png

# The full snpcorr tables are one row per variant and run to hundreds of MB
# across 20 files, so only the strongest 1000 correlations per PC travel. That
# is enough to identify a locus driving a PC. Variant-level, no sample data.
echo "  .. trimming SNP-PC correlation tables to the top 1000 rows"
n_corr=0
for f in "${PCA_DIR}"/snpcorr_[12]*Round_PC[0-9][0-9].tsv; do
    [[ -f "${f}" ]] || continue
    mkdir -p "${BUNDLE}/${PCA_DIR}"
    base="$(basename "${f}" .tsv)"
    head -n 1001 "${f}" > "${BUNDLE}/${PCA_DIR}/${base}_top1000.tsv"
    n_corr=$((n_corr + 1))
done
printf '  ok %-55s %3d file(s)\n' "SNP-PC correlation tables (top 1000)" "${n_corr}"
[[ ${n_corr} -eq 0 ]] && echo "MISSING   SNP-PC correlation tables" >> "${MISSING}"

echo

# ---------------------------------------------------------------------
# Step 2 — genotyped extraction logs
# ---------------------------------------------------------------------
echo "[Step 2] genotyped data extraction"

grab no "genotype extraction logs" \
    genetic_data/genotyped/logs/*.log \
    genetic_data/genotyped/logs/*.out \
    genetic_data/genotyped/logs/*.err
grab no "genotype extraction plink logs" \
    genetic_data/genotyped/*.log

echo

# ---------------------------------------------------------------------
# Step 5 — reference panel and LD scores (logs only)
# ---------------------------------------------------------------------
echo "[Step 5] reference panel and cov-LDSC LD scores"

# plink logs per chromosome: variant counts at every filtering stage.
grab yes "reference panel plink logs (per chromosome)" \
    reference_panel/chr[0-9]/*.log \
    reference_panel/chr[0-9][0-9]/*.log

# The cov-LDSC run log per chromosome: window size, variants scored.
grab yes "cov-LDSC LD score logs (per chromosome)" \
    reference_panel/chr[0-9]/out/*.log \
    reference_panel/chr[0-9][0-9]/out/*.log

grab no "reference panel scheduler logs" \
    reference_panel/chr[0-9]/logs/*.log \
    reference_panel/chr[0-9]/logs/*.out \
    reference_panel/chr[0-9]/logs/*.err \
    reference_panel/chr[0-9][0-9]/logs/*.log \
    reference_panel/chr[0-9][0-9]/logs/*.out \
    reference_panel/chr[0-9][0-9]/logs/*.err

# The .l2.M / .l2.M_5_50 files are a single integer each: the number of
# variants behind the LD scores. Tiny, and needed to sanity-check the
# regression, so they travel even though the scores themselves do not.
grab no "LD score variant counts (.l2.M)" \
    reference_panel/chr[0-9]/out/*.l2.M \
    reference_panel/chr[0-9]/out/*.l2.M_5_50 \
    reference_panel/chr[0-9][0-9]/out/*.l2.M \
    reference_panel/chr[0-9][0-9]/out/*.l2.M_5_50

# The merge into a single genome-wide fileset (Step 10.a).
grab no "genome-wide merge log" \
    reference_panel/merged_chr1_22.log

echo

# ---------------------------------------------------------------------
# Step 6 — GWAS
# ---------------------------------------------------------------------
echo "[Step 6] GWAS and summary statistic parsing"

grab yes "GWAS plink logs" \
    gwas/out/*.log
grab no "GWAS scheduler logs" \
    gwas/logs/*.log gwas/logs/*.out gwas/logs/*.err

# Manhattan and QQ. The QQ is the single most informative plot for telling an
# underpowered GWAS apart from an uncorrected-stratification one.
grab yes "GWAS Manhattan / QQ plot" \
    gwas/out/gwaslab_output/*.png
grab yes "GWASLab parsing log" \
    gwas/out/gwaslab_output/*.log

echo

# ---------------------------------------------------------------------
# Step 7 — arm 1, cov-LDSC
# ---------------------------------------------------------------------
echo "[Step 7] arm 1 - cov-LDSC heritability"

# LDSC writes the estimate into the .log itself, so these three files ARE the
# arm 1 result.
grab yes "cov-LDSC estimates, 3 prevalences (arm 1)" \
    covldsc_h2/logs/cov_ldsc_popPrev_*.log

echo

# ---------------------------------------------------------------------
# Steps 8-9 — arm 2, GREML on the standard GCTA GRM
# ---------------------------------------------------------------------
echo "[Steps 8-9] arm 2 - GREML on the GCTA GRM"

grab no "per-chromosome GRM logs" \
    greml/chr[0-9]/out/*.log \
    greml/chr[0-9][0-9]/out/*.log
grab no "per-chromosome GRM scheduler logs" \
    greml/chr[0-9]/logs/*.log \
    greml/chr[0-9]/logs/*.out \
    greml/chr[0-9]/logs/*.err \
    greml/chr[0-9][0-9]/logs/*.log \
    greml/chr[0-9][0-9]/logs/*.out \
    greml/chr[0-9][0-9]/logs/*.err

grab no "genome-wide GRM merge log" \
    greml/gcta_singleGRM/genomewide/*.log

grab yes "GREML estimates .hsq, 3 prevalences (arm 2)" \
    greml/gcta_singleGRM/greml/*.hsq
grab yes "GREML prevalence grid table (arm 2)" \
    greml/gcta_singleGRM/greml/*prevalence_grid.tsv
grab no "GREML run logs (arm 2)" \
    greml/gcta_singleGRM/greml/*.log \
    greml/gcta_singleGRM/logs/*.log

echo

# ---------------------------------------------------------------------
# Steps 10-11 — arm 3, PC-Relate GRM and GREML on it
# ---------------------------------------------------------------------
echo "[Steps 10-11] arm 3 - PC-Relate GRM and GREML"

# The GRM diagnostics: diagonal summary, off-diagonal mean/SD/quantiles, PSD
# check, SNPs per pair. All aggregate, and the tutorial asks for them by name.
grab yes "PC-Relate GRM diagnostics log (arm 3)" \
    greml_pcrelate/outFolder_h2/grm_with_pcrel_*.log

grab yes "GREML estimates .hsq, 3 prevalences (arm 3)" \
    greml_pcrelate/greml_out/*.hsq
grab yes "GREML prevalence grid table (arm 3)" \
    greml_pcrelate/greml_out/*prevalence_grid.tsv
grab no "GREML run logs (arm 3)" \
    greml_pcrelate/greml_out/*.log
grab no "PC-Relate scheduler logs" \
    greml_pcrelate/logs/*.log \
    greml_pcrelate/logs/*.out \
    greml_pcrelate/logs/*.err

echo

# ---------------------------------------------------------------------
# Step 12 — arms 4 to 7, PCGC and the permutation null
# ---------------------------------------------------------------------
echo "[Step 12] arms 4-7 - PCGC and permutation null"

grab yes "PCGC null summary tables (arms 4-7)" \
    pcgc/out_gcta/pcgc_null_summary_*.tsv \
    pcgc/out_pcrelate/pcgc_null_summary_*.tsv
grab yes "PCGC all estimates, long format (arms 4-7)" \
    pcgc/out_gcta/pcgc_all_estimates_*.tsv \
    pcgc/out_pcrelate/pcgc_all_estimates_*.tsv
grab yes "PCGC master logs" \
    pcgc/out_gcta/pcgc_ldak_*.log \
    pcgc/out_pcrelate/pcgc_ldak_*.log
grab no "PCGC per-prevalence screen logs" \
    pcgc/out_gcta/pcgc_screen_logs/* \
    pcgc/out_pcrelate/pcgc_screen_logs/*

echo

# ---------------------------------------------------------------------
# Safety net: refuse to ship anything that looks like individual-level data
# ---------------------------------------------------------------------
echo "[safety] scanning the bundle for individual-level files"

# These patterns are deliberately extension-anchored. A bare '*covar*' would
# also match covariate_file_no_related_pairs.log, which is the sample-count
# summary we most want to keep - the data files are the .tsv/.txt siblings.
LEAKS="$(find "${BUNDLE}" -type f \( \
        -name '*.grm.bin'    -o -name '*.grm.N.bin' -o -name '*.grm.id'   -o \
        -name '*.grm.gz'     -o -name '*.grm'       -o -name '*.bed'      -o \
        -name '*.bim'        -o -name '*.fam'       -o -name '*.pgen'     -o \
        -name '*.psam'       -o -name '*.pvar'      -o -name '*.gds'      -o \
        -name '*_rds'        -o -name '*.rds'                             -o \
        -name '*covar*.tsv'  -o -name '*covar*.txt' -o -name '*.covar'    -o \
        -name '*pheno*.tsv'  -o -name '*_PCs.tsv'   -o -name '*_IIDs.txt' -o \
        -name 'unrelated_IIDs.txt' -o -name '*kinship_pairs*'             -o \
        -name 'NAToRA_output*'     -o -name 'pipeline_run.log' \
    \) 2>/dev/null || true)"

if [[ -n "${LEAKS}" ]]; then
    echo "  !! individual-level files reached the bundle. Removing them:" >&2
    while IFS= read -r f; do
        [[ -n "${f}" ]] || continue
        echo "     removed: ${f#${BUNDLE}/}" >&2
        rm -f "${f}"
    done <<< "${LEAKS}"
    echo "     (this safety net is doing its job)" >&2
else
    echo "  ok nothing individual-level found"
fi

# The safety net above deletes on a name match, so confirm it did not take the
# sample-count summary with it. That file is the one the whole liability-scale
# comparison is anchored on, and losing it silently would be expensive.
if [[ -f "${PCA_DIR}/covariate_file_no_related_pairs.log" \
   && ! -f "${BUNDLE}/${PCA_DIR}/covariate_file_no_related_pairs.log" ]]; then
    echo "  !! the sample-count summary was removed by the safety net." >&2
    echo "     Restoring it: it holds counts only, no identifiers." >&2
    mkdir -p "${BUNDLE}/${PCA_DIR}"
    cp -p "${PCA_DIR}/covariate_file_no_related_pairs.log" \
          "${BUNDLE}/${PCA_DIR}/"
fi

# Drop any empty directories the filtering left behind.
find "${BUNDLE}" -type d -empty -delete 2>/dev/null || true
mkdir -p "${BUNDLE}"

echo

# ---------------------------------------------------------------------
# Manifest and README
# ---------------------------------------------------------------------
echo "[manifest] listing and checksumming"

if command -v sha256sum >/dev/null 2>&1; then
    HASHER="sha256sum"
elif command -v shasum >/dev/null 2>&1; then
    HASHER="shasum -a 256"
else
    HASHER=""
fi

{
    printf 'file\tbytes\tsha256\n'
    while IFS= read -r f; do
        rel="${f#${BUNDLE}/}"
        [[ "${rel}" == "00_MANIFEST.tsv" ]] && continue
        size="$(wc -c < "${f}" | tr -d ' ')"
        if [[ -n "${HASHER}" ]]; then
            hash="$(${HASHER} "${f}" | awk '{print $1}')"
        else
            hash="NA"
        fi
        printf '%s\t%s\t%s\n' "${rel}" "${size}" "${hash}"
    done < <(find "${BUNDLE}" -type f | sort)
} > "${MANIFEST}"

N_FILES="$(( $(wc -l < "${MANIFEST}") - 1 ))"

if [[ -s "${MISSING}" ]]; then
    N_MISSING="$(wc -l < "${MISSING}" | tr -d ' ')"
else
    N_MISSING=0
    echo "Nothing expected was missing. All arms produced their outputs." > "${MISSING}"
fi

cat > "${README}" <<EOF
h2 in underrepresented populations - results bundle
===================================================

Cohort        : ${COHORT}
Generated     : $(date)
Host          : $(hostname)
Files         : ${N_FILES}
Missing items : ${N_MISSING}   (see 00_MISSING.txt)

The folder layout mirrors the project tree from the tutorial, so every file
sits where you would expect to find it in working_directory_h2.

WHERE THE HEADLINE NUMBERS ARE
  arm 1  cov-LDSC          covldsc_h2/logs/cov_ldsc_popPrev_*.log
  arm 2  GREML + GCTA GRM  greml/gcta_singleGRM/greml/*.hsq
                           greml/gcta_singleGRM/greml/*prevalence_grid.tsv
  arm 3  GREML + PC-Relate greml_pcrelate/greml_out/*.hsq
                           greml_pcrelate/greml_out/*prevalence_grid.tsv
  arms 4-7  PCGC + null    pcgc/out_gcta/pcgc_null_summary_gcta.tsv
                           pcgc/out_pcrelate/pcgc_null_summary_pcrelate.tsv
                           pcgc/out_*/pcgc_all_estimates_*.tsv

SAMPLE COUNTS AND PREVALENCE
  pca_and_such/outFolder_pca_andSuch/covariate_file_no_related_pairs.log
  Total N, cases, controls, and the sample disease proportion that every
  liability-scale transformation in this bundle was computed with.

DIAGNOSTICS WE WILL READ
  GWAS QQ / Manhattan      gwas/out/gwaslab_output/*.png
  PC scatters, 2 rounds    pca_and_such/outFolder_pca_andSuch/plot_PC*.png
  SNP-PC correlations      pca_and_such/outFolder_pca_andSuch/snpcorr_*.png
  PC-Relate GRM behaviour  greml_pcrelate/outFolder_h2/grm_with_pcrel_*.log

WHAT IS NOT HERE
  No individual-level data of any kind: no covariate or phenotype files, no
  GRMs or their .grm.id files, no kinship pair tables, no plink filesets, no
  sample ID lists. The bundle is aggregate and variant-level only.

  Log files that carry individual-level data are not collected either. The
  PCA pipeline log (outFolder_pca_andSuch/pipeline_run.log) is the one that
  matters: it prints a corner of the KING and PC-Relate matrices and the head
  of the kinship pair table, all labelled with sample IDs. Every other log in
  the pipeline reports counts and variances only, and those are all here.

  Left out for size, not privacy - ask us and we will tell you whether we need
  them:
    - GWAS summary statistics (gwas/out/gwaslab_output/*.ldsc.tsv.gz)
    - cov-LDSC LD scores (reference_panel/chr*/out/*.l2.ldscore.gz)
    - the full SNP-PC correlation tables; only the top 1000 rows per PC are
      included, as *_top1000.tsv

VERIFYING
  00_MANIFEST.tsv lists every file with its size and sha256, as the script
  left the folder. If you add or remove anything before zipping, it will no
  longer match - which is fine, just tell us what you changed.
EOF

echo "  ok ${N_FILES} files, manifest written"
echo

# ---------------------------------------------------------------------
# Final report
# ---------------------------------------------------------------------
echo "=============================================================="
echo "Done: $(date)"
echo
echo "  folder : ${BUNDLE_NAME}/"
echo "  files  : ${N_FILES}"
echo "  size   : $(du -sh "${BUNDLE}" | cut -f1)"
echo

if [[ "${N_MISSING}" -gt 0 ]]; then
    echo "  ${N_MISSING} expected item(s) were NOT found:"
    echo
    sed 's/^/    /' "${MISSING}"
    echo
    echo "  If you have not run those steps yet, that is expected - finish"
    echo "  them and run this script again. If you have, something did not"
    echo "  write its output and it is worth checking before you send."
else
    echo "  Every expected output was found."
fi

echo
echo "  Nothing has been compressed. Have a look through the folder, and"
echo "  when you are happy with what is in it:"
echo
echo "      zip -r ${BUNDLE_NAME}.zip ${BUNDLE_NAME}"
echo
echo "  then send us the zip."
echo "=============================================================="
```

</details>

Save it at the top level of `working_directory_h2`, next to the nine analysis folders, and run it with your cohort name:

<details>
<summary>terminal example: collecting the shareables</summary>

```console
(base) [duarte@node1 working_directory_h2]$ nano collect_shareables.sh
(base) [duarte@node1 working_directory_h2]$ bash collect_shareables.sh LARGE-PD
==============================================================
Cohort   : LARGE-PD
Workdir  : /home/duarte/working_directory_h2
Bundle   : /home/duarte/working_directory_h2/h2_shareables_LARGE-PD
==============================================================

[Step 4] PCA and sample summary
  ok sample count / prevalence summary                         1 file(s)
  ok PC-pair scatter plots (round 1)                           5 file(s)
  ok PC-pair scatter plots (round 2)                           5 file(s)
  ok SNP-PC correlation plots (round 1)                       10 file(s)
  ok SNP-PC correlation plots (round 2)                       10 file(s)
  .. trimming SNP-PC correlation tables to the top 1000 rows
  ok SNP-PC correlation tables (top 1000)                     20 file(s)

[Step 2] genotyped data extraction
  ok genotype extraction logs                                 46 file(s)
  ok genotype extraction plink logs                           23 file(s)

[Step 5] reference panel and cov-LDSC LD scores
  ok reference panel plink logs (per chromosome)              66 file(s)
  ok cov-LDSC LD score logs (per chromosome)                  22 file(s)
  ok reference panel scheduler logs                           44 file(s)
  ok LD score variant counts (.l2.M)                          44 file(s)
  ok genome-wide merge log                                     1 file(s)

[Step 6] GWAS and summary statistic parsing
  ok GWAS plink logs                                          22 file(s)
  ok GWAS scheduler logs                                      44 file(s)
  ok GWAS Manhattan / QQ plot                                  1 file(s)
  ok GWASLab parsing log                                       1 file(s)

[Step 7] arm 1 - cov-LDSC heritability
  ok cov-LDSC estimates, 3 prevalences (arm 1)                 3 file(s)

[Steps 8-9] arm 2 - GREML on the GCTA GRM
  ok per-chromosome GRM logs                                  22 file(s)
  ok per-chromosome GRM scheduler logs                        44 file(s)
  ok genome-wide GRM merge log                                 1 file(s)
  ok GREML estimates .hsq, 3 prevalences (arm 2)               3 file(s)
  ok GREML prevalence grid table (arm 2)                       1 file(s)
  ok GREML run logs (arm 2)                                    4 file(s)

[Steps 10-11] arm 3 - PC-Relate GRM and GREML
  ok PC-Relate GRM diagnostics log (arm 3)                     3 file(s)
  ok GREML estimates .hsq, 3 prevalences (arm 3)               3 file(s)
  ok GREML prevalence grid table (arm 3)                       1 file(s)
  ok GREML run logs (arm 3)                                    4 file(s)
  ok PC-Relate scheduler logs                                  2 file(s)

[Step 12] arms 4-7 - PCGC and permutation null
  ok PCGC null summary tables (arms 4-7)                       2 file(s)
  ok PCGC all estimates, long format (arms 4-7)                2 file(s)
  ok PCGC master logs                                          2 file(s)
  ok PCGC per-prevalence screen logs                          12 file(s)

[safety] scanning the bundle for individual-level files
  ok nothing individual-level found

[manifest] listing and checksumming
  ok 475 files, manifest written

==============================================================

  folder : h2_shareables_LARGE-PD/
  files  : 475
  size   : 38M

  Every expected output was found.

  Nothing has been compressed. Have a look through the folder, and
  when you are happy with what is in it:

      zip -r h2_shareables_LARGE-PD.zip h2_shareables_LARGE-PD

  then send us the zip.
==============================================================
```

</details>

Three files inside the bundle are worth opening before you send it:

- **`00_README.txt`** tells us and you which file holds which arm's result.
- **`00_MISSING.txt`** lists anything the script expected and did not find. If you have run every step, this should say nothing was missing. If it names an arm you did run, something did not write its output and it is worth checking before sending rather than after.
- **`00_MANIFEST.tsv`** lists every file with its size and `sha256`, so we can confirm the transfer arrived intact.

Everything in the folder is plain text or png, so `less` and any image viewer are all you need to go through it. Once it looks right to you, zip it and send it over. If `00_MISSING.txt` flagged something you could not resolve, just tell us what and why.



---
# Step 14: Wrapping up

# Notes after you finish the tutorial 

Congrats and thank you so much for your collaboration!

Once you have reached this point I hope you had a smooth run and had fun while doing it! 

You should have ended with the following project structure. Only the directories are shown here, not the files inside them, so do not worry if `ls` gives you more than this:


<details>
<summary>the full project tree (directories only)</summary>

```text
working_directory_h2
├── covldsc_h2
│   └── logs
├── genetic_data
│   ├── genotyped
│   │   └── logs
│   └── imputed -> /path/to/your/imputed_data
├── greml
│   ├── chr1
│   │   ├── logs
│   │   └── out
│   ├── ...(21 additional chromosome folders)
│   └── gcta_singleGRM
│       ├── genomewide
│       ├── greml
│       └── logs
├── greml_pcrelate
│   ├── greml_out
│   ├── logs
│   └── outFolder_h2
├── gwas
│   ├── logs
│   └── out
│       └── gwaslab_output
├── pca_and_such
│   └── outFolder_pca_andSuch
│       └── NAToRA_Public
│           └── ...
├── pcgc
│   ├── out_gcta
│   │   └── pcgc_screen_logs
│   └── out_pcrelate
│       └── pcgc_screen_logs
├── programs
│   ├── cov-ldsc
│   │   └── ...
│   ├── gcta-1.95.3-linux-x86_64
│   │   └── ...
│   └── ldsc
│       └── ...
└── reference_panel
    ├── chr1
    │   ├── logs
    │   └── out
    └── ...(21 additional chromosome folders)
```

</details>

Two small things you may see that are not in the tree above. Unzipping the GCTA archive can leave a `programs/__MACOSX/` folder behind, which is a harmless artifact of how the archive was packaged and can be deleted. And `programs/` also holds the loose executables (`plink2`, `plink`, `ldak6.3.linux`), which do not show up here because only folders are listed.


# Plotting

If you are interested in generating a visual representation of the h2 estimates of your cohort across the different methods and prevalences, you can edit the following script:

The script has empty values where each estimator should go, and you can hard code the specific values and edit some legends inside the script (I apologize for the hard coding). Nonetheless, I can provide a script that hunts for the specific values in the text files and does it for you. I will be working on that for the immediate future. But in the meantime, this is the first draft. 


<details>
<summary><strong>View script: <code>plotting_h2_estimates.r</code></strong></summary>

<br>

```r

library(ggplot2)
e <- function(short, label, family, h2, se, panel)
  data.frame(short, label, family, h2, se, panel, stringsAsFactors = FALSE)
# Panel order, left to right. Edit here to reorder or add a prevalence.
prev <- c("K = 0.5%", "K = 1%", "K = 2%")
est <- rbind(
  # ---- K = 0.5% -------------------------------------------------------
  e("cov-LDSC",          "cov-LDSC - in sample ref. panel", "cov-LDSC", 0, 0, prev[1]),
  e("GREML\nGCTA",       "GREML + standard GCTA GRM",       "GREML",    0, 0, prev[1]),
  e("GREML\nPC-Rel",     "GREML + PC-Relate GRM",           "GREML",    0, 0, prev[1]),
  e("PCGC\nGCTA",        "PCGC + standard GCTA GRM",        "PCGC",     0, 0, prev[1]),
  e("PCGC\nPC-Rel",      "PCGC + PC-Relate GRM",            "PCGC",     0, 0, prev[1]),
  e("PCGC\nGCTA GRM Permutation", "PCGC shuffle kinship matrix (GCTA)",     "PCGC",     0, 0, prev[1]),
  e("PCGC\nPC-Rel GRM Permutation", "PCGC shuffle kinship matrix (PC-Rel)",     "PCGC",     0, 0, prev[1]),
  # ---- K = 1% ---------------------------------------------------------
  # TODO: paste the K = 0.01 re-run here, replacing NA, NA.
  e("cov-LDSC",          "cov-LDSC - in sample ref. panel", "cov-LDSC", 0, 0, prev[2]),
  e("GREML\nGCTA",       "GREML + standard GCTA GRM",       "GREML",    0, 0, prev[2]),
  e("GREML\nPC-Rel",     "GREML + PC-Relate GRM",           "GREML",    0, 0, prev[2]),
  e("PCGC\nGCTA",        "PCGC + standard GCTA GRM",        "PCGC",     0, 0, prev[2]),
  e("PCGC\nPC-Rel",      "PCGC + PC-Relate GRM",            "PCGC",     0, 0, prev[2]),
  e("PCGC\nGCTA GRM Permutation", "PCGC shuffle kinship matrix (GCTA)",     "PCGC",     0, 0, prev[2]),
  e("PCGC\nPC-Rel GRM Permutation", "PCGC shuffle kinship matrix (PC-Rel)",     "PCGC",     0, 0, prev[2]),
  # ---- K = 2% ---------------------------------------------------------
  # TODO: paste the K = 0.02 re-run here, replacing NA, NA.
  e("cov-LDSC",          "cov-LDSC - in sample ref. panel", "cov-LDSC", 0, 0, prev[3]),
  e("GREML\nGCTA",       "GREML + standard GCTA GRM",       "GREML",    0, 0, prev[3]),
  e("GREML\nPC-Rel",     "GREML + PC-Relate GRM",           "GREML",    0, 0, prev[3]),
  e("PCGC\nGCTA",        "PCGC + standard GCTA GRM",        "PCGC",     0, 0, prev[3]),
  e("PCGC\nPC-Rel",      "PCGC + PC-Relate GRM",            "PCGC",     0, 0, prev[3]),
  e("PCGC\nGCTA GRM Permutation", "PCGC shuffle kinship matrix (GCTA)",     "PCGC",     0, 0, prev[3])
  e("PCGC\nPC-Rel GRM Permutation", "PCGC shuffle kinship matrix (PC-Rel)",     "PCGC",     0, 0, prev[3]),
)

# 2. PALETTE
pal <- c(
  "cov-LDSC - in sample ref. panel"  = "#35618F",
  "GREML + standard GCTA GRM"       = "#B26A22",
  "GREML + PC-Relate GRM"           = "#DDA96A",
  "PCGC + standard GCTA GRM" = "#7E4260",
  "PCGC + PC-Relate GRM"     = "#B98BA3",
  "PCGC shuffle kinship matrix (GCTA)"     = "#2d0b1c",
  "PCGC shuffle kinship matrix (PC-Rel)"     = "#30051b"
)

# 3. TEXT

plot_title    <- "SNP heritability of PD in LARGE-PD" # Make sure to edit the name of your cohort
plot_subtitle <- "n=7,178, ~8M snps (maf 0.01, imputed r2>0.8 + hap map3 + only genotyped)" #and here the right numbers for sample size and number of snps
y_lab         <- expression(paste("SNP ", italic(h)^2, " (liability scale)"))
legend_title  <- "Estimator"
plot_caption  <- paste(
  "Estimates across different prevalences (K) of PD",
  "Box spans the point estimate +/- 1 SE; whiskers +/- 1.96 SE.",
  "All models use alpha = -1, except from LDAK GRM that uses the Human Default Model",
  "Kinship matrix shuffling: 100 permutations per K, mean h2 estimate and SD",
  sep = "\n"
)
out_file <- "h2_estimates.png"
out_w    <- 16
out_h    <- 6.5
out_dpi  <- 320

# 4. DERIVED QUANTITIES
est <- transform(
  est,
  lower = h2 - se,
  upper = h2 + se,
  ymin  = h2 - 1.96 * se,
  ymax  = h2 + 1.96 * se
)
# Factor levels follow row order.
# unique() each estimator appears once per panel:
# factor() rejects duplicated levels.
est$short <- factor(est$short, levels = unique(est$short))
est$label <- factor(est$label, levels = unique(est$label))
est$panel <- factor(est$panel, levels = prev)
stopifnot(all(levels(est$label) %in% names(pal)))
# catches a typo in `pal`
# Each estimator must appear exactly once per panel, or the boxes silently
# shift sideways and the panels stop being comparable.
stopifnot(all(table(est$short, est$panel) == 1))

# 5. PLOT

p <- ggplot(est, aes(x = short, fill = label)) +
  # h2 = 0 reference.
  geom_hline(yintercept = 0, linetype = "22", colour = "grey45", linewidth = 0.4) +
  geom_boxplot(
    aes(ymin = ymin, lower = lower, middle = h2, upper = upper, ymax = ymax),
    stat      = "identity",
    width     = 0.55,
    colour    = "grey20",
    linewidth = 0.45,
    alpha     = 0.9
  ) +
  geom_text(
    aes(y = ymax, label = sprintf("%.3f", h2)),
    vjust  = -0.9,
    size   = 2.7,
    colour = "grey25"
  ) +
  scale_fill_manual(values = pal, name = legend_title) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.16))) +
  # Shared y-axis across panels
  facet_wrap(~ panel, nrow = 1) +
  labs(title = plot_title, subtitle = plot_subtitle,
       x = NULL, y = y_lab, caption = plot_caption) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.35),
    axis.text.x        = element_text(size = 7.5, lineheight = 0.95,
                                      colour = "grey20"),
    strip.text         = element_text(size = 11, colour = "grey30",
                                      face = "plain",
                                      margin = margin(t = 2, b = 8)),
    panel.spacing.x    = unit(1.4, "lines"),
    axis.title.y       = element_text(margin = margin(r = 9)),
    axis.line.x        = element_line(colour = "grey70", linewidth = 0.4),
    plot.title         = element_text(face = "bold", size = 14.5, hjust = 0.5,
                                      margin = margin(b = 3)),
    plot.subtitle      = element_text(size = 10, colour = "grey35", hjust = 0.5,
                                      margin = margin(b = 14)),
    plot.caption       = element_text(size = 8, colour = "grey45", hjust = 0,
                                      margin = margin(t = 20)),
    plot.caption.position = "plot",
    plot.title.position   = "plot",
    legend.position    = "bottom",
    legend.title       = element_text(size = 9.5),
    legend.text        = element_text(size = 8.8),
    legend.key.size    = unit(0.85, "lines"),
    legend.box.margin  = margin(t = 8, b = 4),
    plot.margin        = margin(14, 18, 20, 14)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
print(p)
ggsave(out_file, p, width = out_w, height = out_h, dpi = out_dpi, bg = "white")
cat("wrote", out_file, "\n")

```
</details>
