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

The last two terms — the genotype–environment covariance and the genotype-by-environment interaction variance — are usually left out of the heritability analysis that we are about to.

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

> **SNP heritability (`h²_SNP`, sometimes `h²_g`)** — the proportion of phenotypic variance explained by the additive effects of the genotyped (or imputed) variants included in the model (usually, Single Nucleotide Polymorphisms).

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
2. **Characterize heterogeneity, methodological uncertainty and challenges** when computing this estimates, in light of a field whose methods were developed and validated mainly on European-ancestry data.

---

## 3. Methods overview

A comprehensive review of heritability estimation methods is beyond the scope of this repository. A curated reference list is available [here](https://docs.google.com/document/d/1mcoMHNsUat0rzlDItxcBWFDSfRwRMm4rFGBU8VqcIG4/edit?usp=sharing). Disclaimer, the list is not intended to be comprehensive and its reading is not mandatory, but it definitely will help to satisfy your curiosity to a certain extent. 

On that note, we can make a broad classification of the current methods for estimating `h2_SNP` based on the concept similarity that they exploit. And some notes on their challenges when applied to URPs. 

	1. Methods that exploit genotypic similarity among unrelated individuals: require the computation of Genomic Relatedness Matrices (GRMs), tend to give higher heritability estimates and scale up quickly as sample sizes increase. 
		a. The most famous example for these type of methods is embodied by the Genome-Wide Complex Trait Analysis Genomic relatedness matrix restricted maximum likelihood or GCTA-GREML toolkit. Which uses a linear mixed model framework to regress phenotypic similarity between individuals on their genotypic similarity. 
		b. And the alternative that uses a regression of phenotype correlation on genotype correlations, known as Haseman-Elston regressions for quantitative traits or the Phenotype Correlations Genotype Correlations (PGCG) method for binary traits. These methods require less assumptions that the GCTA-GREML (do not require assuming an entire probabilistic model), are more robust under case ascertainment bias and a bit more computationally efficient. 
		c. Both a) and b) methods estimate heritability from a GRM, a table of how genetically similar each pair of participants is. On the logic that if a trait is heritable, people who share more genome should also resemble ach other more in that trait. This requires the matrix to capture *recent shared genealogy*, but the standard calculation compares everyone to a single cohort-average allele frequency, which in URPs that tend to be admixed describes no one accurately. For example, two unrelated participants who both carry high proportions of a given ancestry will deviate from the average in the same direction at thousands of variants, and shared ancestry could be mistaken as shared recent genealogy, breaking some of the assumptions made by the methods in question. Hence, for URPs, an ancestry-aware GRM is highly advisable to compute and benchmarked alongside classical methods (see below). 
	
	2. Methods that that exploit linkage disequilibrium among genetic variants (LD scores): require the presence of an LD reference panel appropriate to the population in question, tend to give more conservative (lower) heritability estimates. Classically exemplified in the Linkage Disequilibrium Score Regression method (LDSC)
		a. For European populations this is the go-to method when sample sizes are big, since the core of the method is based on regressing the chi-square of GWAS summary statistics on the LD scores from the reference panel. So you will only need summary statistics and One-Thousand Genomes Project derived LD scores.
		b. For URPs, appropriate reference panel are usually lacking, and the admixture that characterizes several of the URPs violates some of the mathematical assumptions these methods require. Hence, Individual data is required to compute an in-sample reference panel and adjust the LD scores for the long range LD that is present in admixed populations (using the same principal components that are included in the GWAS that generates the summary statistics to be used, see below). 

