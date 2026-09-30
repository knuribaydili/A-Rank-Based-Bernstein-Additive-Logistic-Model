# B-RCM++: A Rank-Based Bernstein Additive Logistic Model with Spearman-Informed Regularization for Interpretable Binary Classification

Code and output files accompanying

> Baydili, K.N. A Rank-Based Bernstein Additive Logistic Model with Spearman-Informed Regularization for Interpretable Binary Classification. Submitted to *Mathematics* (MDPI).

Every table and figure in the paper can be regenerated from the three R scripts in this repository, and each one is traced below to the output file behind it.

---

## 1. Contents

| File | Purpose | Produces |
|---|---|---|
| `BRCM_simulation.R` | The estimator, the 18 simulation comparators, the factorial Monte Carlo study (216 scenarios), the screening analysis and the sensitivity analysis | Tables 4 to 14; Figures 1, 2, 4, 5, 6; Supplementary Tables S1 to S24, S30, S31 |
| `BRCM_realdata.R` | The benchmark study on five public datasets with 19 methods, including the explainable boosting machine (EBM) | Tables 3, 15 to 20; Figures 7 to 9; Supplementary Tables S25 to S29 and Figure S1 |
| `KOSUM_HIGHDIM.R` | A stand-alone driver that runs only the dedicated high-dimensional screening analysis and checks its output (*kosum* is Turkish for "run") | Table 8, Figure 3 |
| `BRCMsim/` | Output of the reported simulation run; written when the scripts are run; not included in this repository | see Section 6.1 |
| `highdim/` | Output of the reported screening run; written when the scripts are run; not included in this repository | see Section 6.2 |
| `BRCMreal/` | Output of the reported benchmark run; written when the scripts are run; not included in this repository | see Section 6.3 |

The implementation of B-RCM++ is identical in `BRCM_simulation.R` and `BRCM_realdata.R`; each script is self-contained so that either can be run on its own.

## 2. Requirements

The reported results were produced with **R 4.4.2** on **Windows 11 x64**.

R packages, with the versions used:

| Package | Version | Used for |
|---|---|---|
| mgcv | 1.9-1 | generalized additive model |
| rpart | 4.1.24 | decision tree |
| randomForest | 4.7-1.2 | random forest |
| e1071 | 1.7-16 | support vector machine |
| xgboost | 1.7.11.1 | XGBoost |
| lightgbm | 4.6.0 | LightGBM |
| catboost | 1.2.8 | CatBoost |
| mlbench | 2.1-9 | the five benchmark datasets |
| reticulate | 1.47.0 | bridge to Python, for EBM |
| openxlsx | 4.2.8.1 | Excel copies of the output tables |

```r
install.packages(c("mgcv", "rpart", "randomForest", "e1071", "xgboost",
                   "lightgbm", "mlbench", "reticulate", "openxlsx", "remotes"))
```

`catboost` is not distributed through CRAN. On Windows it can be installed from the binary release:

```r
remotes::install_url(
  "https://github.com/catboost/catboost/releases/download/v1.2.8/catboost-R-windows-x86_64-1.2.8.tgz",
  INSTALL_opts = c("--no-multiarch", "--no-test-load"))
```

For other systems see https://catboost.ai/docs/en/installation/r-installation-binary-installation. The estimator itself, ordinary and ridge-penalized logistic regression, the four ablations and all evaluation code use base R only.

**Python, for EBM only.** The explainable boosting machine is called from R through `reticulate`. The reported run used a conda environment named `brcm-ebm` with the `interpret` package installed from PyPI; the Python and `interpret` versions were not recorded in the session information of that run.

```r
reticulate::conda_create("brcm-ebm")
reticulate::py_install("interpret", envname = "brcm-ebm", pip = TRUE)
```

`BRCM_realdata.R` finds this environment by itself. If Python is installed elsewhere, point to the interpreter before sourcing the script:

```r
Sys.setenv(BRCM_PYTHON = "C:/path/to/python.exe")
```

Neither `BRCM_simulation.R` nor `KOSUM_HIGHDIM.R` needs Python.

## 3. Quick start

Each script runs a built-in self-test of 14 checks when it is sourced, covering among other things the exact reconstruction of predictions from the plotted components, the orthogonality of the interaction surfaces and the determinism of the fitting code. A run stops if any check fails.

To confirm that the environment works before committing to a long run, use a small replication count in a separate output folder:

```r
Sys.setenv(BRCM_NSIM = "10", BRCM_CORES = "4", BRCM_OUT = "BRCMtest")
source("BRCM_simulation.R")
```

To check the EBM installation without running the benchmark:

```r
Sys.setenv(BRCM_AUTORUN = "FALSE")
source("BRCM_realdata.R")
brcm_ebm_check()   # reports the interpreter found and the result of a trial fit
```

## 4. Reproducing the results

Each script writes every file it produces inside a single output folder in the working directory and nothing outside it. Settings are read from environment variables, which must be set before the script is sourced; start each run from a fresh R session.

### 4.1 Simulation study (Sections 3.1 to 3.6)

```r
Sys.setenv(BRCM_NSIM = "1000", BRCM_CORES = "12", BRCM_OUT = "BRCMsim")
source("BRCM_simulation.R")
```

Sourcing the file runs `brcm_run_all()`, which performs, in order: the factorial study of 216 scenarios with 1000 replications each; the screening analysis at p = 10, 20 and 50 with 500 replications per cell; the sensitivity analysis with 200 replications per setting; and finally writes all tables and figures. A run from scratch therefore reproduces Table 8 and Figure 3 as well.

### 4.2 Screening analysis only (Table 8, Figure 3)

```r
source("KOSUM_HIGHDIM.R")
```

This runs only the 18 cells of the screening analysis (p = 10, 20, 50; n = 500, 1000; three processes; 500 replications each) without touching the factorial study, and then checks the output for internal consistency: among other things, that the proportion of fits admitting at least one false surface never exceeds the mean number of false surfaces, and that under the null it equals the proportion admitting any surface. It writes `02_highdim_final.rds`, `T_highdim_screening_final.csv`, `F3_screening_by_dimension.png` and `sessionInfo_highdim_final.txt`.

The reported screening results come from this script rather than from the factorial run; the reason is given in Section 9. To redraw Figure 3 from the saved table without repeating the run:

```r
Sys.setenv(BRCM_AUTORUN = "FALSE")
source("BRCM_simulation.R")
draw_screening_figure(read.csv("T_highdim_screening_final.csv"), "F3_screening_by_dimension.png")
```

### 4.3 Benchmark study (Section 3.7)

```r
Sys.setenv(BRCM_CORES = "10", BRCM_OUT = "BRCMreal")
source("BRCM_realdata.R")
```

The five datasets are loaded from the `mlbench` package; no download is needed. Each dataset is evaluated by ten repeats of five-fold stratified cross-validation. At start-up the script prints `EBM reachable in all N workers`. If any worker cannot reach the Python interpreter the run stops, because an EBM row averaged over an unknown subset of repeats would not be comparable with the others.

### 4.4 Environment variables

| Variable | Default | Scripts | Meaning |
|---|---|---|---|
| `BRCM_NSIM` | 1000 | simulation | replications per factorial scenario |
| `BRCM_REPEATS` | 10 | benchmark | cross-validation repeats per dataset |
| `BRCM_FOLDS` | 5 | benchmark | outer folds |
| `BRCM_CORES` | all cores minus one | both | parallel workers |
| `BRCM_OUT` | `BRCMsim` / `BRCMreal` | both | output folder |
| `BRCM_SEED` | 1000 / 2000 | both | base random seed |
| `BRCM_AUTORUN` | TRUE | both | FALSE loads the functions without running anything |
| `BRCM_STRICT` | TRUE | both | stop if a required package, or EBM, is unavailable; FALSE allows an exploratory run without them |
| `BRCM_PROGRESS` | 20 | both | progress lines per scenario or dataset |
| `BRCM_OOFK` | 3 | both | inner folds for the out-of-fold interaction evidence |
| `BRCM_R2THRESH` | 0.01 | both | evidence threshold for admitting an interaction surface |
| `BRCM_LLTHRESH` | 0.002 | simulation | threshold for the alternative log-loss selection criterion |
| `BRCM_PYTHON` | none | benchmark | path to the Python interpreter for EBM |

The reported results use the defaults, with the numbers of workers given in Section 5.

## 5. Running time

All times refer to the machine used for the reported runs.

| Stage | Workers | Time |
|---|---|---|
| Factorial study | 12 for the first 183 scenarios, 28 for the remaining 33 | about 56 h for the first 183 scenarios; the total was not recorded |
| Screening analysis (`KOSUM_HIGHDIM.R`) | sequential | about 32 h |
| Sensitivity analysis | sequential | not recorded separately |
| Benchmark study | 10 | about 20 to 40 min |

Cost in the factorial study is dominated by the p = 20 scenarios: with 12 workers a scenario took about 7 min at p = 4 but roughly 30 to 60 min at p = 20. The screening and sensitivity analyses run their replications one after another, so adding workers does not shorten them; the p = 50 cells account for about 80% of the running time of the screening analysis.

**Interruptions.** The factorial and benchmark studies write a checkpoint after every scenario and every dataset. If a run stops, sourcing the script again with the same `BRCM_OUT` resumes from the last completed unit and loses at most the one in progress. A manifest stored with the checkpoints records the settings of the run, and a resumed run stops if they differ, so results computed under different settings cannot be mixed. The screening analysis has no checkpoints and must run to completion.

For runs of this length, disable sleep mode and keep the output folder on a local disk rather than a cloud-synchronized one. One of the runs behind this repository was interrupted while its output folder was on a synchronized drive, and was resumed from its checkpoints.

## 6. Output files and where they appear in the paper
The output files are not included in this repository; the table below records which file each table and figure of the paper comes from, once the scripts have been run.
The file names inside the output folders predate the final numbering of the paper, which is why, for example, Figure 4 of the paper is the file `F3_HL_forest_AUC.png`.

Axis labels follow the typographic convention of the journal: negative numbers are set with a minus sign (U+2212) rather than a hyphen, through the helper `axis_minus()` defined in both scripts. The figures published with the paper are exactly the files written by `write_figures()` and `write_figures_rd()`; no separate figure script is used.

**Redrawing the figures without repeating the analyses.** Both figure writers work from the saved results, refitting only the single models behind the component-function figures, so the complete set can be redrawn in a few minutes:

```r
Sys.setenv(BRCM_AUTORUN = "FALSE")
source("BRCM_simulation.R")
write_figures(readRDS("BRCMsim/sim.rds"), "BRCMsim")      # Figures 1, 2, 4, 5, 6 and the screening figure
source("BRCM_realdata.R")
write_figures_rd(readRDS("BRCMreal/realdata.rds"), "BRCMreal")   # Figures 7, 8, 9 and Figure S1
```

### 6.1 Simulation (`BRCMsim/`)

| Paper | File |
|---|---|
| Table 4 | `tables/T1_overall.csv` |
| Tables 5, 6, 12 | `tables/DATA_summary.csv` (per-scenario means, averaged as described in each table note) |
| Tables 7, 11 | `tables/T_diagnostics.csv` |
| Table 9 | `tables/DATA_paired_scenario.csv` |
| Table 10 | `tables/T_paired_summary.csv` |
| Table 13 | `tables/T_runtime.csv` |
| Table 14 | `tables/T_interaction_diagnostics.csv` (full grid; the paper shows selected rows) and `03_sensitivity.rds` |
| Figure 1 | `figures/F1_AUC_by_n.png` |
| Figure 2 | `figures/F2_LogLoss_by_DGP.png` |
| Figure 4 | `figures/F3_HL_forest_AUC.png` |
| Figure 5 | `figures/F4_component_functions.png` |
| Figure 6 | `figures/F5_interaction_surface.png` |
| Tables S1 to S24 | `tables/T02_detail_*.csv` to `tables/T25_detail_*.csv` (the 72 scenarios at p = 4) |
| Table S30 | `tables/DATA_paired_scenario.csv` |
| Table S31 | `tables/T_diagnostics.csv` |
| not in the paper | `tables/T26_detail_*.csv` to `tables/T73_detail_*.csv`, the same detail for p = 10 and p = 20 |
| everything above | `01_factorial.rds`, the replication-level results (see Section 7) |
| — | `Paper_Tables.xlsx`, an Excel copy of the tables; `sessionInfo.txt` |

### 6.2 Screening analysis (`highdim/`)

| Paper | File |
|---|---|
| Table 8 | `T_highdim_screening_final.csv` |
| Figure 3 | `F3_screening_by_dimension.png` |
| — | `02_highdim_final.rds`, `sessionInfo_highdim_final.txt` |

### 6.3 Benchmark (`BRCMreal/`)

| Paper | File |
|---|---|
| Table 3 | `tables/RT0_datasets.csv` |
| Table 15 | `tables/RT_overall.csv` |
| Table 16 | computed from `realdata.rds`; see below |
| Tables 17, 18 | `tables/RD_paired.csv` |
| Table 19 | `tables/RD_diagnostics.csv` |
| Table 20 | `tables/RD_interaction_summary.csv` (per fold: `tables/RD_interaction_by_fold.csv`) |
| Figure 7 | `figures/RF_AUC_by_dataset.png` |
| Figure 8 | `figures/RF_LogLoss_by_dataset.png` |
| Figure 9 | `figures/RF_pima_components.png` |
| Figure S1 | `figures/RF_HL_forest_AUC.png` |
| Tables S25 to S29 | `tables/RT1_PimaIndiansDiabetes.csv` to `tables/RT5_HouseVotes84.csv` |
| — | `realdata.rds` (all repeat-level results), one `.rds` per dataset, `RealData_Tables.xlsx` |

Table 16 restricts the comparison to the four datasets completed by every method shown, excluding Ionosphere, and omits ordinary logistic regression, which fails on two of those four. It is computed from the saved results:

```r
z  <- readRDS("BRCMreal/realdata.rds")
ds <- c("PimaIndiansDiabetes", "BreastCancer", "Sonar", "HouseVotes84")
ms <- setdiff(names(z$raw[[ds[1]]]), "LR")
m  <- function(k) sapply(ms, function(mm)
        mean(sapply(ds, function(d) mean(z$raw[[d]][[mm]][, k], na.rm = TRUE))))
tab <- data.frame(AUC = 100 * m("Test_AUC"), Brier = m("Test_Brier"), LogLoss = m("Test_LogLoss"))
round(tab[order(-tab$AUC), ], 4)
```

## 7. Large files

The replication-level results are not distributed with this repository: `01_factorial.rds` is 435 MB and exceeds GitHub's file-size limit. Running `BRCM_simulation.R` reproduces it, and the file is available from the author on request.

## 8. Known issues

### 8.1 A worker process can terminate on one configuration

In the factorial study, one configuration repeatedly terminated a parallel worker process outright rather than raising an R error: twenty predictors, n = 100 and prevalence 0.50. After the 70/30 split this leaves about 70 training observations for 20 predictors. In the order in which the script visits the scenarios, the cases were:

| Scenario | Name | What happened |
|---|---|---|
| 147 of 216 | `logit_mixed_p20_n100_p50` | terminated a worker in an early run, which then stopped |
| 183 of 216 | `logit_interact_p20_n100_p50` | terminated a worker in an early run, which then stopped |
| 192 of 216 | `logit_pure_int_p20_n100_p50` | in the reported run, 2 of 1000 replications could not be completed |

When this happens the master process reports `error reading from connection` or `error writing to connection`. Early versions of the script stopped the whole run at that point. The released version contains the failure instead: it rebuilds the cluster, reruns the affected chunk one replication at a time, and skips any replication that still terminates a worker when run on its own. A skipped replication is recorded rather than silently dropped. In the reported run this happened only in scenario 192, which therefore rests on **998 replications for every method**: `tables/T_diagnostics.csv` shows `Attempted = 998` for that scenario, and the paper states the exception in Sections 2 and 3.4. Every other scenario rests on 1000 replications.

The failure was not traced to a particular comparator. It is most likely a crash in the compiled code of one of the comparator libraries on a very small training partition. `brcm_diagnose()` in `BRCM_simulation.R` fits each method in a separate worker and can be used to isolate it.

Because every replication draws its data from its own fixed seed, rebuilding the cluster and rerunning a chunk reproduces exactly the values an uninterrupted run would have produced; only the skipped replications are absent.

At the end of a run in which the cluster was rebuilt, R prints warnings of the form `closing unused connection` (in a Turkish locale, `kullanılmayan bağlantı kapatılıyor`). They come from the discarded cluster and do not affect the results.

### 8.2 Other messages

`wilcox.test` may warn `requested conf.level not achievable`. This happens when nearly all paired differences between two methods are zero, typically when B-RCM++ and its interaction-free ablation make identical predictions; the estimate and the p-value are still returned, and the interval is slightly narrower than nominal.

Ordinary logistic regression fails in many replications at n = 100, particularly with twenty predictors, because the data are separable. These fits are recorded as failures in the diagnostics, not dropped.

CatBoost creates a `catboost_info/` folder in the working directory. It is not part of the output and can be deleted.

## 9. A correction made during review

The dedicated screening analysis was first run with a count of family-wise false selection that treated a fit as a false selection whenever any interaction surface was admitted. Under the null process that is the correct quantity, because every admitted surface is false; on the two interaction processes it is not, because recovering the generating pair also counts. The count was corrected to record only fits that admit at least one surface other than the generating pair, the mean number of false surfaces is now counted directly, and the screening analysis was rerun in full with `KOSUM_HIGHDIM.R`, extended at the same time from p = 10 and 20 to p = 10, 20 and 50. Table 8 and Figure 3 report the corrected run. The factorial study was not affected by the error and was not rerun.

## 10. Citation, licence and contact

If you use this code, please cite the paper:

```bibtex
@unpublished{Baydili2026BRCM,
  author = {Baydili, K{\"u}r{\c{s}}ad Nuri},
  title  = {A Rank-Based Bernstein Additive Logistic Model with Spearman-Informed
            Regularization for Interpretable Binary Classification},
  note   = {Manuscript submitted to Mathematics (MDPI)},
  year   = {2026}
}
```

**Licence.** The code and the output files in this repository are released under the MIT licence; see `LICENSE`.

**Use of generative AI.** Generative AI assistants were used to revise and debug the code, which was originally written by the author; see Section 2.7 of the paper.

**Contact.** Kürşad Nuri Baydili, Department of Biostatistics and Medical Informatics, Hamidiye Faculty of Medicine, University of Health Sciences Turkey, Istanbul, Türkiye. kursadnuri.baydili@sbu.edu.tr · ORCID [0000-0002-2785-0406](https://orcid.org/0000-0002-2785-0406)

