B-RCM++: A Rank-Based Bernstein Additive Logistic Model with Spearman-Informed Regularization
Replication code for the manuscript
> Baydili, K.N. \*A Rank-Based Bernstein Additive Logistic Model with Spearman-Informed Regularization for Interpretable Binary Classification.\* Manuscript submitted for publication, 2026.
This repository accompanies a manuscript that is currently under peer review. Full bibliographic details will be added here once the article is published.

---
Overview
B-RCM++ is a single penalized logistic model designed so that the object that is interpreted and the object that predicts are the same object. Each predictor enters through a smooth Bernstein expansion of its mid-rank empirical distribution; the ridge penalty on each predictor block is weighted by that predictor's absolute Spearman correlation with the outcome; and a small number of pure pairwise interaction surfaces are admitted only when out-of-fold evidence supports them.
Two identifiability constraints make the fitted object exactly decomposable:
every main-effect block is centred to sum to zero on the training sample, and
every interaction surface is orthogonalized against the intercept and the main-effect blocks of its own two predictors.
As a consequence, reconstructing the linear predictor from the displayed components returns the model's own predicted probabilities. This is not asserted but verified: each script contains an executable self-test (see the Self-test section below) that checks this to numerical tolerance.
This repository contains the two self-contained R scripts that produce every number, table and figure reported in the manuscript.
---
Repository contents
File	Description
`SUP-BRCM\_simulation.R`	Monte Carlo simulation study: 7 data-generating processes x 3 sample sizes x 3 prevalences = 63 scenarios, 1000 replications each, 18 methods.
`SUP-BRCM\_realdata.R`	Benchmark study on five public datasets under 10 repeats of stratified 5-fold cross-validation, same 18 methods.
`README.md`	This file.
`LICENSE`	MIT licence.
`sessionInfo.txt`	Environment record: R and package versions used for the reported results.
Each script is standalone: it defines the estimator, all comparator methods, the evaluation protocol, the tables and the figures, and requires no other file in this repository.
---
Requirements
R version. 4.5.1
Required packages (both scripts fail immediately if any is missing under the default strict mode):
Package	Used for
`mgcv`	generalized additive model comparator
`rpart`	classification tree comparator
`randomForest`	random forest comparator
`e1071`	support vector machine comparator
`xgboost`	XGBoost comparator
`lightgbm`	LightGBM comparator
`catboost`	CatBoost comparator
`mlbench`	benchmark datasets (`SUP-BRCM\_realdata.R` only)
Optional packages: `foreach`, `doParallel` (parallel execution), `openxlsx` (writes the combined `.xlsx` workbook; the individual `.csv` files are written regardless).
`parallel`, `stats` and `utils` are part of base R and need no installation. The estimator itself uses no external package: the rank transformation, Bernstein basis, penalized IRLS and interaction search are implemented in base R only.
Installing
```r
install.packages(c("mgcv", "rpart", "randomForest", "e1071",
                   "xgboost", "lightgbm", "mlbench",
                   "foreach", "doParallel", "openxlsx"))
```
`catboost` is not on CRAN. Its R package is distributed as a platform-specific binary from the project's GitHub releases page, so it must be installed separately from the command above.
Step 1. Open https://github.com/catboost/catboost/releases and find the release you want. Under that release's Assets list, locate the R binary for your system. The file names follow this pattern:
```
catboost-R-windows-x86\_64-{VERSION}.tgz     Windows
catboost-R-linux-x86\_64-{VERSION}.tgz       Linux
catboost-R-darwin-universal2-{VERSION}.tgz  macOS (Intel and Apple silicon)
```
Step 2. Install it, substituting the file name you found:
```r
install.packages("remotes")

remotes::install\_url(
  "https://github.com/catboost/catboost/releases/download/v1.2.7/catboost-R-windows-x86\_64-1.2.7.tgz",
  INSTALL\_opts = c("--no-multiarch", "--no-test-load"))
```
The example above installs version 1.2.7 for Windows; change the version number and the file name to match your platform and the release you chose. Note that the version appears twice in the URL: once after `download/v` and once inside the file name.
Step 3. Confirm the installation and note the version:
```r
packageVersion("catboost")
```
Any recent release works: the scripts call only `catboost.load\_pool`, `catboost.train` and `catboost.predict`, and impose no version requirement. The exact version used for the reported results is recorded in `sessionInfo.txt`.
If the installation stops with an error about hard-coded temporary paths, add `"--no-staged-install"` to `INSTALL\_opts`.
---
Reproducing the reported results
Both scripts run end to end from a single entry point and write their outputs to a directory created in the current working directory.
```bash
Rscript SUP-BRCM\_simulation.R
Rscript SUP-BRCM\_realdata.R
```
Equivalently, from an interactive R session:
```r
source("SUP-BRCM\_simulation.R")   # runs automatically and returns `sim`
source("SUP-BRCM\_realdata.R")     # runs automatically and returns `res`
```
Each script ends with `if (AUTORUN \&\& sys.nframe() == 0L) ...`, so sourcing it from the top level launches the full run, while sourcing it from inside another function does not.
Seeds. All randomness is seeded deterministically (`BASE\_SEED` = 1000 for the simulation, 2000 for the benchmark study), and fold assignment is deterministic and class-stratified. Re-running either script on the same machine with the same package versions reproduces the reported values exactly.
Quick smoke test (a few minutes)
```bash
BRCM\_NSIM=5 BRCM\_OUT=smoke\_test Rscript SUP-BRCM\_simulation.R
```
This runs the entire pipeline with 5 replications per scenario instead of 1000. It produces the same file structure with unstable numbers, and confirms that all packages, the parallel backend and the output paths are working.
---
Configuration
Both scripts read their settings from environment variables, so no line of code needs to be edited to change a run.
Variable	Default (simulation)	Default (benchmark)	Meaning
`BRCM\_NSIM`	`1000`	—	replications per scenario
`BRCM\_REPEATS`	—	`10`	repeats of cross-validation
`BRCM\_FOLDS`	—	`5`	outer folds per repeat
`BRCM\_OOFK`	`3`	`3`	inner folds for penalty selection and interaction scoring
`BRCM\_SEED`	`1000`	`2000`	base seed
`BRCM\_CORES`	all but one	all but one	parallel workers
`BRCM\_OUT`	`BRCM\_ciktilar`	`BRCM\_realdata\_ciktilar`	output directory
`BRCM\_R2THRESH`	`0.01`	`0.01`	out-of-fold R-squared threshold for admitting an interaction
`BRCM\_AUTORUN`	`TRUE`	`TRUE`	run automatically when sourced at top level
`BRCM\_STRICT`	`TRUE`	`TRUE`	require all comparator packages
Example — a shorter run on 4 cores writing to a different directory:
```bash
BRCM\_NSIM=100 BRCM\_CORES=4 BRCM\_OUT=my\_run Rscript SUP-BRCM\_simulation.R
```
Development mode
To load the functions without running anything (useful for inspecting the estimator or running the self-test):
```r
Sys.setenv(BRCM\_AUTORUN = "FALSE")
source("SUP-BRCM\_simulation.R")
brcm\_self\_test()
```
To run with only the base-R comparators, skipping any boosting library you have not installed:
```bash
BRCM\_STRICT=FALSE Rscript SUP-BRCM\_simulation.R
```
Note that a non-strict run does not reproduce the reported results, because the skipped methods are absent from the comparison.
---
Self-test
Both scripts define `brcm\_self\_test()`, which verifies eleven properties of the estimator and returns `TRUE` only if all pass:
rank map consistency: `u\_train(x)` equals `u\_test(x, x)`
an interaction is selected when a true interaction is present
each pure interaction surface is orthogonal to the intercept and to the main-effect blocks of its own two predictors
prediction equals the sum of the plotted components
the plotted interaction surface equals the model's own contribution
the plotted categorical bars equal the model's own contribution
continuous components sum to zero on the training data
the fitting procedure is deterministic and consumes no random numbers
a single-level categorical predictor is handled without error
an unseen test category is handled without error
the null option is genuine: `r2\_thresh = Inf` yields zero interactions
Properties 3 to 7 are the executable form of the identifiability and decomposability claims made in Section 2.1 of the manuscript. To run it:
```r
Sys.setenv(BRCM\_AUTORUN = "FALSE")
source("SUP-BRCM\_simulation.R")
brcm\_self\_test(verbose = TRUE)
```
Expected output: eleven `PASS` lines and `TRUE`.
---
Outputs
`SUP-BRCM\_simulation.R` (default directory `BRCM\_ciktilar/`)
Main tables
File	Content
`T1\_overall.csv`	overall performance of all 18 methods across the 63 scenarios
`T\_paired\_summary.csv`	paired comparison of B-RCM++ against each of the other 17 methods
`T\_diagnostics.csv`	per-scenario convergence rates and interaction-selection counts
`T\_interaction\_diagnostics.csv`	sensitivity to the interaction threshold and the Spearman offset
`T\_runtime.csv`	median fitting time per evaluation, by method
`T01\_detail\_\*.csv` … `T21\_detail\_\*.csv`	full median (IQR) results for each process x prevalence combination
Replication-level data
File	Content
`DATA\_summary.csv`	scenario x method x metric summaries
`DATA\_paired\_scenario.csv`	scenario-level paired differences with confidence intervals
`sim.rds`	complete replication-level object (large)
Figures
File	Content
`F1\_AUC\_by\_n.png`	mean AUC as a function of sample size
`F2\_LogLoss\_by\_DGP.png`	mean log-loss by data-generating process
`F3\_HL\_forest\_AUC.png`	scenario-level paired differences with 95% confidence intervals
`F4\_component\_functions.png`	fitted component functions on the nonlinear process
`F5\_interaction\_surface.png`	recovered pure interaction surface
`Paper\_Tables.xlsx` collects the tables in a single workbook if `openxlsx` is available.
`SUP-BRCM\_realdata.R` (default directory `BRCM\_realdata\_ciktilar/`)
File	Content
`RT0\_datasets.csv`	dataset dimensions, event rates and preprocessing rules
`RT\_overall.csv`	overall performance across the five datasets
`RT1\_\*.csv` … `RT5\_\*.csv`	per-dataset results for all 18 methods
`RD\_summary.csv`	dataset x method x metric summaries
`RD\_paired.csv`	repeat-level paired differences
`RD\_diagnostics.csv`	successful fits out of ten attempts, by method and dataset
`RF\_AUC\_by\_dataset.png`, `RF\_LogLoss\_by\_dataset.png`	performance by dataset
`RF\_HL\_forest\_AUC.png`	paired differences against every comparator
`RF\_pima\_components.png`, `RF\_pima\_interaction.png`	fitted components on Pima Indians Diabetes
`realdata.rds`	complete repeat-level object
`RealData\_Tables.xlsx`	combined workbook (requires `openxlsx`)
---
What the scripts implement
The estimator
`brcm\_fit(X, y, cat\_mask, ...)` returns a fitted model; `brcm\_predict(fit, Xnew)` returns predicted probabilities. Key arguments:
Argument	Default	Meaning
`m`	`10`	Bernstein degree for main effects
`m2`	`3`	Bernstein degree for interaction surfaces
`Kmax`	`4`	maximum number of admitted interaction surfaces
`lam\_grid`	`c(0.25, 0.5, 1, 2, 4)`	candidate global penalty values
`lam\_int`	`1.0`	fixed ridge parameter for interaction blocks
`r2\_thresh`	`0.01`	out-of-fold R-squared required to admit a surface
`min\_n\_inter`	`300`	minimum training size below which the interaction search is disabled
`c\_spear`	`0.05`	offset in the Spearman penalty weight
`top\_feat`	`8`	screening size for interaction candidates
Fold-local estimation
Every quantity that depends on the data is recomputed inside each training fold and applied unchanged to the held-out part. This applies to the mid-rank empirical distribution maps, the block centring constants, the Spearman penalty weights, the additive fit, the residuals used to score candidate interactions, and the orthogonalization operator. No quantity used to construct or score a candidate interaction is estimated on the observations against which it is scored, so the reported out-of-fold R-squared is a genuine held-out quantity.
Compared methods (18)
B-RCM++; logistic regression; ridge-penalized logistic regression; classification tree; generalized additive model; random forest; support vector machine; XGBoost; LightGBM; CatBoost; Platt-calibrated variants of the random forest, XGBoost, LightGBM and CatBoost; and four nested ablations of the proposed model (A1-LinearRidge, A2-RankLinear, A3-RankBernstein, A4-PlusSpearman) in which consecutive models differ by exactly one component.
Data-generating processes (7)
`logit\_mixed`, `logit\_nl` (nonlinear), `logit\_corr` (correlated predictors), `lpm\_mixed` (linear probability), `logit\_interact`, `logit\_pure\_int` (interaction with no main effects), `null\_mixed` (no signal).
Benchmark datasets (5)
Pima Indians Diabetes, Breast Cancer Wisconsin, Ionosphere, Sonar and Congressional Voting Records, all obtained from the `mlbench` package. Missing-data rules and outcome encodings are applied inside the training folds only.
---
Known behaviour worth noting
These are documented so that a reader who re-runs the code is not surprised by them; all are discussed in the manuscript.
The interaction search does not activate on every dataset. It requires at least `min\_n\_inter = 300` training observations and at least two continuous predictors. Under 5-fold cross-validation the training partitions of Ionosphere (280) and Sonar (166) fall below that threshold, and all predictors in Congressional Voting Records are categorical. On these three datasets B-RCM++ therefore reduces to its additive form by construction and is identical to the A4-PlusSpearman ablation. This is a structural consequence of the protocol, not evidence that those datasets lack interaction structure.
LightGBM fails on Ionosphere. That dataset contains a predictor with a single observed level, and the design-matrix expansion used to pass the data to the boosting library applies treatment contrasts, which are undefined for a factor with fewer than two levels. The call raises an error before any tree is grown, on every repeat.
Ordinary logistic regression fails on three datasets. Ionosphere, Sonar and Congressional Voting Records admit complete or quasi-complete separation at the fold sizes used here, so the maximum-likelihood estimate diverges.
The generalized additive model succeeds on only 2 of 10 repeats of Ionosphere, the smoothing basis being unidentifiable on folds in which several predictors are near-collinear.
Failure criterion. A repeat counts as successful only if the method returns finite predicted probabilities together with a finite operating threshold for every observation in the dataset, and both outcome classes are present in the resulting prediction vector. Anything else — an explicit error, a non-convergence flag or a non-finite value — is recorded as a failure. The same rule is applied identically to all 18 methods.
---
Creating the environment record
Before archiving, record the exact environment in which the reported results were produced:
```r
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
```
Run this in the same R session used for the final runs and add the resulting file to the repository. The Data Availability Statement of the manuscript refers to it.
---
License
Released under the MIT License. See the `LICENSE` file for the full text.
The benchmark datasets are redistributed by the `mlbench` package and originate from the UCI Machine Learning Repository; they are not included in this repository and remain subject to their own terms.
---
Citation
This archive has a permanent DOI and can be cited directly:
```bibtex
@software{brcmplusplus\_code,
  author    = {Baydili, K{\\"u}r{\\c s}ad Nuri},
  title     = {Replication code for B-RCM++: A Rank-Based Bernstein Additive
               Logistic Model with Spearman-Informed Regularization},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {\[ZENODO DOI]}
}
```
The accompanying manuscript is under review and not yet published. Until it appears, it may be referred to as:
```bibtex
@unpublished{brcmplusplus\_manuscript,
  author = {Baydili, K{\\"u}r{\\c s}ad Nuri},
  title  = {A Rank-Based Bernstein Additive Logistic Model with
            Spearman-Informed Regularization for Interpretable
            Binary Classification},
  year   = {2026},
  note   = {Manuscript submitted for publication}
}
```
---

Contact
Kürşad Nuri Baydili — kursadnuri.baydili@sbu.edu.tr or knuribaydili@hotmail.com
Department of Biostatistics and Medical Informatics, Hamidiye Faculty of Medicine, University of Health Sciences Turkey, Istanbul, Türkiye
https://orcid.org/0000-0002-2785-0406
