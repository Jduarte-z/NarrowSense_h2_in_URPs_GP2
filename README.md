# Estimation of narrow sense heritability in Underrepresented Populations within the Global Parkinson's Genetics Program 

## Brief introduction of important concepts, objectives and methods

## Intro
A complex trait, under the lens of quantitative population genetics, could be explained in the following way:
Y ~ A + D + I + E

Where Y is the phenotype of interest and is a function of: "A", the additive genetic effects that contribute to the development of the trait; "D", the dominant genetic effects; "I", the epistatic interactions; and "E", the enviromental component. 

Since these parameters tend to vary across individuals and populations. At any given cohort, and across individuals, the variance (V) of the parameters could be writen like this:

VY ~ VA + VD + VI + VE

Where VY is the total phenotypic variance, VA is the additive genetic variance, VD is the dominant genetic variance, VI is variance in epistatic interactions, and VE the environmnental variance. 

Hence, the heritability of a phenotype or trait could be defined as the total phenotypic variance that could be explain by the variance in genetic effects (since variance is a squared metric, it is written like H2 and h2, see below discussion to know the difference). And is specific of the population and environment at which the individulas of interest belong.

Heritability could be broken down into broad-sense (H2) and narrow-sense heritability (h2). Braod-sence heritability is essentially VA + VD + VI. And narrow-sense heritability referss exclusivelly to the additive genetic effects (VA). 
Classically, since the mid-20th century, heritability of complex traits has been studied through twin studies, in which researchers gathered multiple monzygotic (MZ) and dizygotic (DZ) twins, and modeled how much extra phenotypic similarity in MZ pairs must be due to genetic variance (given their higher genetic sharing) compared to DZ twins, and partitioning the total phenotypic variance into genetic and environmental components. Furthermore, with the advent of genome-wide association studies, the possibility of estimating the heritability of complex traits through genotyping or sequencing data became available. However, GWAS heritability estimates have been classically focused on additive genetic effects of common single nucleotide polymorphisms (SNPs) (an intrinsic limitation of the most commonly used genotyping tools). Hence, for many traits and phenotypes, the heritability estimated through twin studies (closer to H2) versus the one estimated through genome-wide association studies (closer to h2) are discordant. Being twin data the one giving much higher estimates. A problem known as the "missing heritability of complex diseases" (the "problem" is basically that GWAS-based heritability estimates tend to be lower than the ones from twin data).  

The discussion about the reasons behind the "missing heritability" problem are outside of the scope of this tutorial, but useful information could be found elsewhere (https://pmc.ncbi.nlm.nih.gov/articles/PMC2942068/). 

## Objective of this tutorial

This repository is intended to showcase the steps needed to undertake the estimation of narrow sense heritability (h2) in underrepresented populations using raw genetic, phenotypic and genome-wide association data. 

## Objective of the project 

Our main aims are:
1) layout the analytical framework to estimate h2 explained by single nucleotide polymorphisms in underrepresented populations. 
2) estimate the h2 corresponding to genome-wide significant loci discovered so far.
3) provide a landscape of the potential heterogeneity and challenges regarding this estimations in underrepresented populations in light of current European-dominant field. 


## Methods overview 

Nowadays, multiple tools stand out to infer h2 in the field of quantitative population genomics. The BIG majority of these have being developed using data derived from European populations. However, recent adaptations for underrepresented and admixed populations have being designed. In this repository, an explicit step by step, hands-on adaptation is described. Specifically in the context of underrepresented cohorts within the field of Parkinson's Disease (PD). 

A comprehensive review of the methods available for h2 estimation is out of the scope of this repository. Relevant references to this topics have been gathered in this document: https://docs.google.com/document/d/1mcoMHNsUat0rzlDItxcBWFDSfRwRMm4rFGBU8VqcIG4/edit?usp=sharing (not mandatory to read but useful to satisfy curiosity, to a certain extent). However, this list is not intended to cover all the literature available and should NOT be considered comprehensive. 

The initial framework that was approved by GP2's Project Proposal, Approval, and Execution Working Group consisted on benchmarking two methods:

1) Covariate Linkage Disequilibrium Score regression (cov-LDSC), that corresponds to the adaptation of the classic LDSC method but for admixed/underrepresented populations. Classic LDSC is designed to work only with GWAS summary statistics and a reference panel LD scores. However, since a reference panel is usually not available for URPs and some of the mathematical assumptions are violated in the presence of admixture, cov-LDSC solves these issues by computing an in-sample reference panel and adjusting the LD scores by genetic principal component analysis. 

2) Genome-wide complex trait analysis (GCTA) Genomic relatedness matrix (GRM) restricted maximum likelihood (GREML), that would be used as a complementary estimate to the ones provided by cov-LDSC. Since at lower sample sizes (an usual scenario for URPs), would be computationally feasible and potentially provide a useful alternative estimate to analyze, because it includes raw individual data in the estimation of h2 itself, contrary to cov-LDSC [note: cov-LDSC also uses raw individual level data, but in the construction of the in-sample reference panel, the calculations themselves only use the LD scores and summary statistics, aka, de-identified data). 

These two methods could be considered amongst the most popular ones in the literature to estimate narrow sense heritability. However, there are plenty of newer versions and different approaches to answer this question (see the document with relevant references above). So as the project advances, and based on your insight as well, we can discuss and to include additional complementary tools that have different assumptions, advantages, and disadvantages as part of sensitivity analysis.




