# rae3: recurrent copy-number alterations in tumour cohorts

A reimplementation of **RAE** (Taylor BS et al. 2008, *Functional copy-number alterations in cancer*,
PLoS ONE 3(9):e3179), which finds the DNA copy-number alterations that recur across a cohort of tumours.
For each tumour, RAE fits soft detectors for gain, amplification, hemizygous loss and homozygous
deletion. It then builds a unified breakpoint profile across the cohort and tests each region against a
permutation background model. 


The repository holds two front ends to the same C++17 core:

|         |                                                                 |
|---------|-----------------------------------------------------------------|
| `rae3/` | the R package; `rae3/inst/include/rae/` is the header-only core |
| `cli/`  | `rae`, a standalone command-line tool built from the same core  |

For the same input and seed, R and the command-line tool give byte-identical output.

## Install

R (≥ 4.0, a C++17 compiler):

```r
# install.packages("remotes")
remotes::install_github("RodrigoGM/rae3", subdir = "rae3")
```

Command line (CMake ≥ 3.16):

```
git clone https://github.com/RodrigoGM/rae3.git && cd rae3
cmake -S cli -B build && cmake --build build
./build/rae --help
```

## Quick start

```
./build/rae cohort.seg                       # hg38, sequencing preset
./build/rae cohort.seg --build hg19 --preset cgh --out results/cohort.rae   # array CGH
```

```r
library(rae3)
res <- rae(seg_data_frame, build = "hg19")                  # or a segmented QDNAseq object
res <- rae(qdnaseq_object, bin_noise = TRUE)                # noise (MAPD) from the unsegmented bins
res <- rae(seg_data_frame, purity_ploidy = pp)              # optional purity/ploidy: ID, purity, ploidy
res$roi; res$samples
rae_plot_genome(res); rae_plot_frequency(res)
res <- rae_read("results/cohort.rae")                       # read a command-line run into R
```

**Input:** a `.seg` table of log2 ratios (`ID chrom loc.start loc.end num.mark seg.mean`, optional
`seg.sd`), autosomes only. Optional: bins or per-sample noise statistics, purity/ploidy, matched normals,
a population-CNV catalogue.

## Presets

- **`sequencing`** (default): places each tumour's detectors on its copy-number ladder. The one-copy step
  comes from supplied purity and ploidy, or is estimated from the spacing of the segment means. On
  low-noise sequencing data the noise model RAE was designed with often finds no single-copy peak. On
  BriTROC-1 (low-pass whole-genome sequencing, 265 tumours), agreement on amplifications with the
  published calls was κ 0.65 with the estimated ladder and 0.95 with published purity/ploidy, against
  0.52 for `cgh`.
- **`cgh`**: for array CGH. It behaves like the released RAE code (beta 0.91), which reproduces the
  published results.
- **`taylor2008`**: for array CGH. It follows the method as the 2008 text describes it.

`?rae_params` lists every parameter. The tags `G-xx` and `D-xx` in the documentation refer to the
development specification (gap register and code comparison), which is not part of this repository.

## Output

Every run writes the same eight files: `<prefix>.regions.tsv`, `.roi.tsv` (regions of interest, with
q-values and overlap with common population CNVs), `.samples.tsv`, `.segments.tsv`, `.calls.tsv`,
`.null.tsv`, `.run.json` (all parameters, the seed and which fallbacks fired) and `.report.html`.

**Fraction of the genome altered.** CNA burden is the concept; the fraction of the genome altered (FGA) is
how it is measured. `.samples.tsv` gives two FGAs per tumour:

- **`rae_fga`**: counts the segments the RAE calls gained or lost. This is the CNA burden of Taylor et al.
  2010 and Hieronymus et al. 2014, and it reproduces their published values (median ratio 0.998) and their
  5.41 % prognostic cut-off. `n_altered` is the number of those segments.
- **`fga`**: counts the segments whose log2 ratio is above `--fga-threshold` (0.2) or below minus
  `--fga-threshold-loss` (by default the same value). With 0.2 this is the cBioPortal / MSK-IMPACT FGA.
  For single cells, set both higher, just short of the one-copy steps (log2 3/2 = 0.58 up, log2 1/2 = −1
  down).

Both divide by the genome each tumour's segments cover; `--fga-assembly` divides by the whole assembly.
See *Output* in `?rae`.

## Validation

- **Unit tests:** about 400 test expectations, including R/CLI identity.
- **Simulated cohort:** 23 of 26 pre-registered targets recovered with no false positives, and p-values
  calibrated under a global null.
- **Two published prostate cohorts (`cgh` preset):**
  - CNA burden (`rae_fga`): Pearson 0.9999 (Hieronymus 2014, n = 104) and 0.997 (Taylor 2010, n = 147).
  - Taylor 2010 regions: 58 of 66 focal events (Table S3) at q ≤ 0.10, and all 7 arm-level events
    (Table S4).

## Genome builds

hg38 (default), hg19, T2T-CHM13, `auto` (the best-fitting bundled build) and `custom` (no centromere or
hotspot table).

## Citation

Cite the method as Taylor BS, Barretina J, Socci ND, DeCarolis P, Ladanyi M, Meyerson M, Singer S, Sander
C. *Functional copy-number alterations in cancer.* PLoS ONE 2008;3(9):e3179, and describe this software as
a reimplementation of RAE (rae3). Cite Taylor et al. 2010 (Cancer Cell 18:11–22) and Hieronymus et al. 2014
(PNAS 111:11139–11144) when using the discrete calls or the CNA burden.

## Bundled data

The core includes filtered or derived copies of public reference data. Cite them when you use the
corresponding output:

- **Recombination hotspots:** HapMap phase II release 21 (McVean et al. 2004; Winckler et al. 2005; Frazer
  et al. 2007), lifted to each build; used by the null model.
- **Chromosome sizes, centromeres, cytobands:** UCSC Genome Browser.
- **Common population CNVs** (overlap columns of `roi.tsv`):
  - DGV Gold Standard 2016, hg19 (MacDonald et al. 2014);
  - gnomAD-SV v4.1, hg38 (CC0);
  - for T2T-CHM13: HPRC release 2 long-read deletions, plus gnomAD-SV gains lifted to T2T-CHM13.

Cancer-gene labels for the plots come from OncoKB. They are not bundled: `rae_fetch_cancer_genes()`
downloads them. OncoKB is free for academic research use; see its terms.

## Licence

BSD 3-Clause (`LICENSE`).
