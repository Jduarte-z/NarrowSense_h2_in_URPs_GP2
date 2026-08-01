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
conda activate h2_project

pwd
/working_directory

mkdir pca_and_such

cd pca_and_such/

ls
covar.txt genotyped_data_plink2_prefix.pgen genotyped_data_plink2_prefix.psam genotyped_data_plink2_prefix.pvar pcs_pipeline_h2_project.r

Rscript pcs_pipeline_h2_project.r \
  --input_plink_file genotyped_data_plink2_prefix \
  --input_format pfile \
  --covar_file covar.txt \
  --removeHighLDregions yes \
  --n_cores 4 2>&1 | tee pcs_pipeline_h2_project.log
```
The script in question:
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
plink2_bin          <- arg("plink2_bin", "plink2")
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
# HELPERS — generic names  for common files 
#########################################################

# Convert a PLINK bed/bim/fam triple (given by prefix) to a GDS file,
# then open it, print a summary, and close it again.
bed2gds <- function(prefix, out_gds) {
  snpgdsBED2GDS(
    bed.fn    = paste0(prefix, ".bed"),
    bim.fn    = paste0(prefix, ".bim"),
    fam.fn    = paste0(prefix, ".fam"),
    out.gdsfn = out_gds
  )
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
    # Columns: chrom  start  end  label   (tab-separated, no header)
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
#load the gds file
geno_reader <- GdsGenotypeReader(filename = pca_gds)
geno_data <- GenotypeData(geno_reader)

#load the KING matrix
KINGmat <- readRDS(king_mat_rds)

#run pcair
mypcair_1round <- pcair(geno_data, kinobj = KINGmat, divobj = KINGmat, num.cores=n_cores, eigen.cnt = n_pcs)
saveRDS(mypcair_1round, pcair_r1_rds)
close(geno_data)
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

# Convert for round 2 PC-AiR input: scaleKin=1
# scaleKin=1 keeps the kinship-coefficient scale expected by pcair's kinobj
mypcrel_1round_mat <- pcrelateToMatrix(mypcrel_1round,
                                       scaleKin = 1,
                                       thresh=NULL,
                                       verbose = TRUE)
saveRDS(mypcrel_1round_mat, pcrel_r1_mat_rds)
#######################################################################################################################
# Run the second round of PC-AiR using the kinships from PC-Relate round 1 that are more accurate than those from KING. 
#however, the divergence signal is still required to come from KING
#######################################################################################################################
# kinobj  → PC-Relate r1 (accurate recent relatedness)
# divobj  → KING (retains divergence signal across ancestry groups)
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
# # Run the second round of PC-Relate, this is for the h2 project         still under development, not implemented yet 
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

```
</details>
The " 2>&1 | tee commandLineRun.log  " part is just to save all the print statements into a log file. You can skip it if you want.

The script by default will create an output folder named outFolder_pca_andSuch/ with all the downstream files that we will need.

the input for the script could be plink1 files, so you would run --input_plink_file genotuped_data_plink1_prefix --input_format bfile

you can also modify the default parameters of the script, for more information type: Rscript pcs_pipeline_h2_project.r --help

however, we advise you to keep the default parameters as they are and only change the flags in the example of the command line.

