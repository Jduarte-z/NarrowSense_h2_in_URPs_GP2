# Estimation of Narrow-Sense SNP Heritability in Underrepresented Populations within the Global Parkinson's Genetics Program

*A hands-on tutorial: concepts, objectives, and methods.*

---

## 1. Background

### 1.1 The quantitative genetic model

Under the classical quantitative genetics framework, an individual's phenotype for a complex trait can be decomposed as:

```
P = G + E = (A + D + I) + E
```

where `P` is the observed phenotype, `G` genetics, and `E` environment. Within "`G1`"`A` is the additive genetic effect (the sum of average effects of alleles), `D` is the dominance effect (within-locus allelic interaction), `I` is the epistatic effect (between-locus interaction), and `E` is the environmental deviation.

Because these components vary across individuals, the corresponding decomposition of the total phenotypic variance in a given population is:

```
VP = VA + VD + VI + VE + 2·Cov(G,E) + VGxE
```

The last two terms (the genotype–environment covariance and the genotype-by-environment interaction variance) are usually left out of the heritability analysis that we are about to.

Hheritability estimates should not be interpreted as environment-independent estimates, but rather as population and its particular
distribution of environment specific.  

In this GitHub, we do not provide a comprehensive coverage of these topics. For useful references you coud refer to books like: Falconer, D.S. & Mackay, T.F.C. (1996). Introduction to Quantitative Genetics, and Lynch, M. & Walsh, B. (1998). Genetics and Analysis of Quantitative Traits.


### 1.2 Broad and narrow-sense heritability

Heritability is a **ratio**: it is the proportion of phenotypic variance attributable to genetic variance.

- **Broad-sense heritability:** `H² = (VA + VD + VI) / VP`
- **Narrow-sense heritability:** `h² = VA / VP`

Two properties are worth to mention:

1. Heritability is a property of a **population in an environment**, not of a trait and certainly not of an individual. The same trait can have different heritabilities in multiple cohorts that have multiple environmental and/or genetic components. 

2. Heritability *within* groups carries no information about the causes of mean differences *between* groups, see [J.G. Schraiber & M.D. Edge (2024)](https://www.pnas.org/doi/10.1073/pnas.2319496121)


### 1.3 Twin studies, GWAS, and the "missing heritability" problem

Since the mid-20th century, heritability of complex traits was estimated primarily from twin and family designs. The classical ACE model partitions phenotypic variance into additive genetic (A), shared environmental (C), and unique environmental (E) components by exploiting the difference in genetic sharing between monozygotic (~1.0) and dizygotic (~0.5) twin pairs. Falconer's estimator, `h² = 2(r_MZ − r_DZ)`, is nominally a narrow-sense estimator, though in practice dominance, epistasis, shared environment, and assortative mating all load onto the A term.

With the advent of genome-wide association studies, heritability could be estimated from genotype data in samples of "unrelated individuals". These methods estimate a distinct quantity:

**SNP heritability (`h²_SNP`, sometimes `h²_g`)**: the proportion of phenotypic variance explained by the additive effects of the genotyped (or imputed) variants included in the model (usually, Single Nucleotide Polymorphisms).

`h²_SNP` is generally smaller than `h²` because:

- Genotyping arrays interrogate common variation and tag rare causal variants poorly.
- Even where causal variants are present, the additive random-effect model assumes independent, exchangeable SNP effects, which is an approximation.
- Non-additive variance is excluded by construction (and dominant, epistatic effects, etc. are left out)

The gap between twin-based and GWAS-based estimates is the **"missing heritability" problem**. It reflects a combination of genuinely untagged variance, non-additive contributions, and the nature of family-based estimates (e.g. shared environment and assortative mating). A useful entry point to this literature is Manolio et al. (2009), [*Finding the missing heritability of complex diseases*](https://pmc.ncbi.nlm.nih.gov/articles/PMC2942068/), since we don't cover this aspect extensively in this repository.  

**Hence, with this tutorial, you will be able to estimate `h²_SNP`, not `h²`**

### 1.4 Binary traits and the liability scale

Parkinson's disease is a dichotomous outcome, which introduces two complications absent from the continuous-trait framework above.

**Scale.** Under the liability threshold model, each individual has an unobserved, normally distributed liability; individuals whose liability exceeds a threshold determined by the population prevalence `K` are affected. Heritability estimated on the observed 0/1 scale depends on the case-control ratio in the sample and is not comparable across studies. Estimates must therefore be reported on the **liability scale**, which is a property of the population rather than of the study design.

**Ascertainment.** Case-control cohorts are, by design, enriched for cases relative to the population (`P >> K`). This oversampling biases the naive observed-to-liability transformation and, more seriously, biases REML-based estimators when covariates such as ancestry principal components are included in the model.

For PD we use `K = 0.005` as the primary population prevalence, consistent with the values used in previous published literature, and report sensitivity across `K ∈ {0.005, 0.01, 0.02}`. Prevalence is not measured with precision in underrepresented populations, and this uncertainty propagates directly into the liability-scale estimate. 

### 1.5 Why estimate h²_SNP, and why in underrepresented populations?

The estimation of `h2_SNP` is relevant by its own means regardless of the population. Consider that it is computed using the same resources allocated for conducting genome-wide association studies, and since it represents the phenotypic variance explained by the same genetic variation screened in a GWAS framework, its value helps us to get closer to a comprehensive mapping of the genetic architecture of a complex trait in a given population. For example, giving a potential upper limit estimate on the discovery power that GWAS could have, or the virtual upper limit of the performance of polygenic risk scores derived from the GWASs in question (since the accuracy of a PRS would approach heritability estimates as sample size and increase genetic variation is captured without bound). 

And beyond this point, exploring `h2_SNP` estimation in underrepresent populations becomes meaningful in several ways. First, it will increase the methodological capacity so it does not depend on European reference data and pipelines. Nearly every heritability estimator in common use was developed and validated on European-ancestry samples, and their behavior in underrepresented, admixed cohorts is comparatively uncharacterized, especially for Parkinson's Disease. Additionally, it will help us to showcase the heterogeneity, challenges, limitations and future directions when computing and interpreting these type of estimates from diverse genomic data. 

---

## 2. Objectives

### Repository objective

Provide a reproducible, step-by-step workflow for estimating SNP heritability in underrepresented and admixed populations from raw genotype, phenotype, and GWAS summary data.

### Project aims

1. **Lay out the analytical framework** for estimating `h²_SNP` in underrepresented populations, with explicit management of admixture, relatedness, PD case ascertainment bias, liability scale transformation and different heritability models. 
2. **Collaborate internationally** to build the capacity needed to estimate `h2_SNP` as accurately as possible in an increasingly diverse genomic landscape. Powered by and for URPs. 
3. **Characterize heterogeneity, methodological uncertainty and challenges** when computing this estimates, in light of a field whose methods were developed and validated mainly on European-ancestry data.

---

## 3. Methods overview

A comprehensive review of heritability estimation methods is beyond the scope of this repository. A curated reference list is available [here](https://docs.google.com/document/d/1mcoMHNsUat0rzlDItxcBWFDSfRwRMm4rFGBU8VqcIG4/edit?usp=sharing). Disclaimer, the list is not intended to be comprehensive and its reading is not mandatory, but it definitely will help to satisfy your curiosity to a certain extent. 

On that note, we can make a broad classification of the current methods for estimating `h2_SNP` based on the concept similarity that they exploit. Plus some caveats on their challenges when applied to URPs. 

1. **Methods that exploit linkage disequilibrium among genetic variants (LD scores)**: they require the presence of an LD reference panel appropriate to the population in question and tend to give more conservative (lower) heritability estimates. Classically exemplified in the Linkage Disequilibrium Score Regression (LDSC) software. Some notes about it:
   - **(a)** For European populations this is the go-to method when sample sizes are big, since the core logic is based on regressing the chi-square of GWAS summary statistics on the LD scores from the reference panel. So you will only need summary statistics and European reference panel-derived LD scores.
   - **(b)** When the mean chi-square of a given GWAS summary statistics is greater than 1, that means that mean association metrics for each variant across the genome are inflated. And this could be due to population stratification and/or the polygenicity of a given trait. Precisely regressing the mean chi-square from a GWAS on a population-matched LD scores could help to differentiate these two, while at the same time computing the heritability of the complex trait. Just FYI, The authors of LDSC recommend a mean chi-square of 1.02 in order to be suitable for LDSC regression. 
   - **(c)** For URPs, appropriate reference panel are usually lacking, and the admixture that characterizes these cohorts violates some of the mathematical assumptions of the method. Hence, Individual data is required to compute an in-sample reference panel and adjust the LD scores for the long range LD that is present in admixed populations (using the same principal components that are included in the GWAS that generates the summary statistics to be used, see below). This is the main logic behind the adaptation of LDSC for admixed population called covariate-LDSC (cov-LDSC)
2. **Methods that exploit genotypic similarity among unrelated individuals**: they require the computation of Genomic Relatedness Matrices (GRMs), tend to give higher heritability estimates and scale up quickly as sample sizes increase. Some notes about them:
   - **(a)** The most famous example for these type of methods is embodied by the Genome-Wide Complex Trait Analysis Genomic relatedness matrix restricted maximum likelihood or GCTA-GREML toolkit. Which uses a linear mixed model framework to regress phenotypic similarity between individuals on their genotypic similarity.
   - **(b)** Alternative methods, use a regression of phenotype correlation on genotype correlations, known as Haseman-Elston regressions for quantitative traits or the Phenotype Correlations Genotype Correlations (PGCG) method for binary traits. These methods require less assumptions that the GCTA-GREML (do not require assuming an entire probabilistic model), are more robust under case ascertainment bias and a bit more computationally efficient.
   - **Both (a) and (b)** estimate heritability from a GRM, a table of how genetically similar each pair of participants is. On the logic that if a trait is heritable, people who share more genome should also resemble each other more in that trait. This requires the matrix to capture *recent shared genealogy*, but the standard calculation compares everyone to a single cohort-average allele frequency, which in URPs that tend to be admixed, describes no one accurately. For example, two unrelated participants who both carry high proportions of a given ancestry will deviate from the average in the same direction at thousands of variants, and shared ancestry could be mistaken as *recent shared genealogy*, breaking some of the assumptions made by the methods in question. Hence, for URPs, an ancestry-aware GRM is highly advisable to compute and benchmarked alongside classical methods (see below).

Furthermore, there are a couple of different additional "axes" to classify `h2_SNP` estimation methods, in particular:

1. **Based on the Heritability model**: every estimator carries a prior assumption about how heritability is distributed across the genome. Specifically, how variant's expected contribution depends on its allele frequency and how much LD it sits in. The general form for describing this is `E[h²_j] ∝ w_j · [f_j(1−f_j)]^(1+α)`, where `f_j` is allele frequency and `w_j` is an LD-based weight.
	- **(a)** All of the above methods tend to do `w_j = 1` and, conventionally, `α = −1` (the so called GCTA model) This makes `E[h²_j]` constant: every SNP is expected to contribute equally, regardless of frequency or LD. LDSC assumes essentially the same thing.
	- **(b)** The human default model described by the creators of [LDAK](https://dougspeed.com/heritability-model/), is often described as an alternative to the classic GCTA model. 
2. **Based on genome-wide partition**: this alludes whether a single variance component for all genetic variants is used, or a split into bins based on different parameters is implemented. 
	- **(a)** All of the above methods (in the way we intend them to be used) use a single variance component, aka, they use genome-wide data all at once. However an alternative is to decompose genome-wide data in bins based on LD scores or minor allele frequency, since it has been shown that heritability tend to vary based on different thresholds for this data. However, the cost of including more parameters tend to be larger standard errors and the need for bigger sample sizes that underrepresented cohorts frequently don't have. 


## 4. Putting the benchmarking pipeline together 

Based on the brief overview of the main methods available, you can quickly realize that even for European populations there is no gold-standard method for estimating SNP heritability. And for URPs there are additional problems to solve. Hence, considering that each estimator differ in what they assume and their specific requirements and limitations in diverse populations, we aim to get a spread of estimates across a pre-specified set of methods and across a wide range of different plausible PD prevalences. 

Our analysis will follow two prongs, a) using methods that exploit  LD scores adapted for URPs, and b) methods that exploit genotypic similarity among unrelated samples. 

For a) we are going to focus on cov-LDSC that requires an in-sample reference panel, computation of adjusted LD scores by principal components and the summary statistics from a GWAS ran on the same samples.  

For b) we are going on GCTA-GREML and PCGC, with both a regular GRM and an ancestry-aware GRM computed through the R package PC-Relate. Hence, a total of four different GRM based heritability results will be reported. Additionally, given that the PCGC method allows for the permutation of the GRM kinships, we will derive estimates under the null hypothesis of no genetic signal when shuffling the classic GRM kinships and the ancestry-aware GRM ones.

Here is a more detailed description on the analysis to be performed:

| # | Arm | GRM | Model | Role in the design |
|---|---|---|---|---|
| 1 | **cov-LDSC**, in-sample LD scores | — | α =−1 | Summary-statistic based, with the caveat of in-sample reference panel, validated in admixed cohorts by Luo et al., and completely independent of GRM construction. |
| 2 | **GREML** + standard GCTA GRM | GCTA | α =−1 | One of the field's reference implementation. This is the number most directly comparable to older published PD heritability estimates, like Keller et al., |
| 3 | **GREML** + PC-Relate GRM | PC-Relate | α =−1 | Isolates GRM construction. Same estimator, same SNPs, same covariates as arm 2 — only ancestry residualisation differs. |
| 4 | **PCGC** + standard GCTA GRM | GCTA | α =−1 | Isolates ascertainment handling. Same GRM as arm 2 — only the liability-scale estimator differs. |
| 5 | **PCGC** + PC-Relate GRM | PC-Relate | α =−1 | Both corrections applied together. |
| 6 | **PCGC permutation** | PC-Relate kinships, shuffled | α =−1 | Null control. |
| 7 | **PCGC permutation** | GCTA kinships, shuffled | α =−1 | Null control. |

### 3.1 Covariate-adjusted LD score regression (cov-LDSC)

Standard LDSC estimates `h²_SNP` by regressing GWAS chi-square statistics on LD scores, typically computed from an external reference panel such as 1000 Genomes. Two of its assumptions fail in admixed cohorts: no matched reference panel exists, and long-range admixture LD violates the assumption that LD is negligible beyond a short genomic window.

cov-LDSC [(Luo et al. 2021)](https://pubmed.ncbi.nlm.nih.gov/33987664/) addresses both. Principal components are projected out of the genotypes **before** LD is computed, so that the LD scores are adjusted for the same covariates included in the association model. LD scores are then computed in-sample rather than from an external panel.

Practical consequences for implementation:

- **Window size.** In European samples, mean LD scores plateau beyond ~1 cM. In admixed American samples they continue to rise without covariate adjustment; after adjustment they plateau at ~20 cM. A 20-cM window is the recommended default.
- **Number of PCs.** Ten PCs is the published recommendation, though the appropriate number depends on the structure of the specific cohort and should be checked empirically.
- **Subsampling.** In-sample LD scores can be computed on a random subset of individuals; unbiased estimates were obtained with as few as 1,000 samples. However in order to use the same genetic data for all the main and sensitivity analysis we are planning on using the full cohort as its own reference for the LD scores. Nonetheless, if it becomes computationally unscalable, we can revise this.
- **MAF filter.** Restrict to MAF > 0.01 for LD score computation.
- **Estimand.** Because LD is computed from array genotypes rather than refernce panel sequence data, cov-LDSC targets the same estimand as GCTA (`h²_g`) rather than LDSC's usual `h²_common`. This makes it directly comparable to the GRM-based estimate below.

### 3.2 GRM-based estimation: GREML and PCGC

GCTA-GREML is the complementary, individual-level estimator to benchmark against cov-LDSC, and even the authors of cov-LDSC benchmarked their real-world data and simulations with this two methods. 
However, 
The original proposal specified GCTA-GREML as the complementary, individual-level estimator. GREML is well suited to the modest sample sizes typical of underrepresented cohorts and uses individual-level genotypes directly.

However, for ascertained case-control data, REML-based estimators give inconsistent liability-scale estimates, and the bias is aggravated — not mitigated — by the inclusion of fixed-effect covariates such as ancestry PCs (Golan et al. 2014; Weissbrod et al. 2018). 

We therefore use **PCGC regression** (Phenotype-Correlation–Genotype-Correlation), implemented in LDAK, as an additional tool to benchmark alongside cov-LDSC and GREML. PCGC estimates liability-scale heritability directly and remains unbiased under ascertainment and in the presence of covariates.

One honest caveat: PCGC's variance behavior depends on the relationship between sample and population prevalence. It is essentially unbiased when sample prevalence exceeds population prevalence, but MSE increases substantially when population prevalence is small and sample prevalence does not exceed it (Ojavee et al. 2022). A case-enriched PD cohort with `K = 0.005` sits in the favorable side, which is the justification for this choice.

### 3.3 Sensitivity analyses

The estimate is a point on a surface defined by several modelling choices, and the honest deliverable is the surface, not the point. Axes to vary:

| Axis | Range |
|---|---|
| Population prevalence `K` | 0.005 (primary), 0.01, 0.02 |
| Genetic architecture model | GCTA-uniform weighting (primary), LDAK-parametric weighting |
| cov-LDSC window size | ≥ 20 cM; confirm plateau criterion (<1% increase per cM) |
| cov-LDSC PC count | Vary; confirm LD score stability |
| GRM construction | Relatedness threshold, MAF filter, LD-pruning scope |

**Note on architecture:** GCTA-uniform weighting (α = −1) is used as the primary model for comparability with the IPDGC literature and to match the cov-LDSC estimand. The LDAK-weighted model is retained as a sensitivity dimension rather than dismissed; the MAF-axis distinction between the two largely washes out when analysis is restricted to common variants.

---

## 4. Repository structure

<!-- TODO: fill in once the scripts are organized -->

```
.
├── 00_qc/
├── 01_ancestry_relatedness/
├── 02_gwas/
├── 03_covldsc/
├── 04_pcgc/
└── docs/
```

## 5. Prerequisites

<!-- TODO: pin versions -->

| Tool | Version | Purpose |
|---|---|---|
| PLINK | 1.9 / 2.0 | Genotype QC and manipulation |
| R | ≥ 4.x | GENESIS, SNPRelate, GWASTools |
| GENESIS | | PC-AiR, PC-Relate |
| cov-LDSC | immunogenomics fork | LD score computation and regression |
| LDAK | | PCGC regression |

Compute environment: SLURM-based HPC cluster. Several steps (notably `pcrelate` at N > 5,000) are not feasible in an interactive session.

## 6. Data access

<!-- TODO: GP2 / LARGE-PD data access statement and tier description -->

## 7. Citation

<!-- TODO -->

## 8. Contact

<!-- TODO -->

---

## Key references

- Bulik-Sullivan et al. (2015) — LD score regression
- Conomos et al. (2016) — PC-Relate, model-free relatedness estimation
- Golan et al. (2014) — PCGC regression
- Huang et al. (2025) — Interpreting SNP heritability in admixed populations
- Lee et al. (2011) — Observed-to-liability scale transformation
- Luo et al. (2021) — cov-LDSC
- Ojavee et al. (2022) — Liability-scale heritability for low-prevalence disease
- Schraiber & Edge (2024) — Within- vs. between-group heritability
- Speed et al. (2017) — Re-evaluation of SNP heritability and architecture models
- Weissbrod et al. (2018) — SNP heritability and genetic correlation in case-control studies
- Yang et al. (2010, 2011) — GCTA and genome partitioning
