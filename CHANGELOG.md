# Changelog

All notable changes to this project are recorded here, grouped by date.

---

## 2026-04-25

- [`32b3d45`](https://github.com/noahlorinczcomi/bpact/commit/32b3d450429bf0c7f399fdd921f18e8a3ab0000b) Removed the GitHub Actions R CMD check workflow that was added on 2026-04-24.
- [`936bfb5`](https://github.com/noahlorinczcomi/bpact/commit/936bfb5a185c1890f9d654f3df8eeb3367992dde) CI/CD configuration changes.
- [`03aec16`](https://github.com/noahlorinczcomi/bpact/commit/03aec16caa3784d1a2fc45b1dcccb356f0d0e700) Added a CODEOWNERS file.
- [`da39360`](https://github.com/noahlorinczcomi/bpact/commit/da3936010de0863236bb05b482973f1e9c954f46) CI/CD configuration changes.

## 2026-04-24

- [`6ce4a16`](https://github.com/noahlorinczcomi/bpact/commit/6ce4a165708bb46074269e16b15125610278f053) CI/CD configuration changes.
- [`450a116`](https://github.com/noahlorinczcomi/bpact/commit/450a1160a00fcc9e438ab88bc7888b89735c37c2) CI/CD configuration changes.
- [`97cf815`](https://github.com/noahlorinczcomi/bpact/commit/97cf815752e7d8c50ad36dc8a4bfa64e5ec12b3f) Resolved all R CMD check errors. Added `R/imports.R` with explicit `@importFrom` directives for every function used from external packages (dplyr, rlang, stats, graphics, methods, utils) and declared global variables with `utils::globalVariables()` to eliminate undefined-variable NOTEs. Also fixed several `matrix()` calls that used the non-standard `nr=`/`nc=` shorthand instead of `nrow=`/`ncol=`.
- [`2433a87`](https://github.com/noahlorinczcomi/bpact/commit/2433a878fed5b494de8c7d26276151dd2a22092e) Added GitHub Actions R CMD check workflow. Fixed a bug in `ptrap()` where the CDF integral over the left ramp was algebraically incorrect; replaced with `y0/(2*(x1-a))*(x-a)^2`. Fixed the method-of-moments fit for the trapezoidal distribution in `mom_trap_distribution()`, which was mistakenly using log-transformed values and computing a sum-of-errors rather than a sum-of-squared-errors. Removed an unused `mu1` parameter from the internal likelihood function in `mh()`. Cleaned up example calls in documentation to pass `R CMD check`.

## 2025-09-12

- [`167838c`](https://github.com/noahlorinczcomi/bpact/commit/167838c8c330a4cae36865125fffedcf50a40351) Merge pull request [#2](https://github.com/noahlorinczcomi/bpact/pull/2) from noah-DG/main.
- [`de83ec0`](https://github.com/noahlorinczcomi/bpact/commit/de83ec0803982a9b5bcff63076b2068f7db76831) Fixed the same `pvalue` → `pval` typo in `posterior_gene()`'s `dplyr::rename()` call that had been fixed in `compositemh()` two commits earlier.
- [`e21d157`](https://github.com/noahlorinczcomi/bpact/commit/e21d157d29a2f0a187c729fa3e77692f95410826) Made `ld.df` and `gent.Rho` explicit required arguments in `compositemh()`, and `ld.df` a required argument in `posterior_gene()`, removing the previous approach of defaulting them to `NULL` and loading them with `data()` inside the function body, which did not work correctly.
- [`e0c9875`](https://github.com/noahlorinczcomi/bpact/commit/e0c9875acd7e2c3d485a3788e1134a83314c970e) Changed `ld.df` and `gent.Rho` defaults in `compositemh()` and `ld.df` in `posterior_gene()` from inline `data(...)` calls (invalid as default argument values) to `NULL`, with explicit `data()` loading at the top of the function body when `NULL` is passed.
- [`e7b35d0`](https://github.com/noahlorinczcomi/bpact/commit/e7b35d06b922f69e0514b92c60dcd1c26f6fec09) Merge pull request [#1](https://github.com/noahlorinczcomi/bpact/pull/1) from noah-DG/main.
- [`1d97edd`](https://github.com/noahlorinczcomi/bpact/commit/1d97edd763606dd6ea08916ad91b973d84f7c95c) Fixed remaining typo `pvalu` → `pval` in `compositemh()`'s `dplyr::rename()` call.
- [`5ec9eeb`](https://github.com/noahlorinczcomi/bpact/commit/5ec9eebc7298b305550046ced36d79f5682485a5) Partially fixed a typo in `compositemh()`'s `dplyr::rename()` call: `pvalue` → `pvalu` (still incorrect; corrected fully in the next commit).
- [`05d26c9`](https://github.com/noahlorinczcomi/bpact/commit/05d26c913471e3ec3f5220f7520fcf795eb3a9aa) Removed a `cat()` progress-printing statement from `functions.R`.

## 2025-09-11

- [`a66edb6`](https://github.com/noahlorinczcomi/bpact/commit/a66edb6e2f1e144dce2fa0626bf3881b362c0cae) Updated README.
- [`18e2a4e`](https://github.com/noahlorinczcomi/bpact/commit/18e2a4e0102f7ed8af2e16f49913014d87a5eb75) Added named column-mapping parameters (`gene_id`, `chromosome`, `position`, `pval`, `null_mean`, `null_variance`, `gamma_shape`, `gamma_rate`) to both `compositemh()` and `posterior_gene()`. Users can now pass data frames whose columns have non-default names; the functions internally rename and select the relevant columns using `dplyr::rename()` and `dplyr::select()`. Also improved parameter documentation for both functions.

## 2025-06-23

- [`3e6d712`](https://github.com/noahlorinczcomi/bpact/commit/3e6d712750b6757c4015ea1f512995f7c4ddd21a) Fixed a documentation typo in `compositemh()`: "refines the prior distributions" → "defines the prior distributions".
- [`cc8fce9`](https://github.com/noahlorinczcomi/bpact/commit/cc8fce9d7c8f77db86a0928d3d94d5eac0895722) Updated README.

## 2025-05-12

- [`c561430`](https://github.com/noahlorinczcomi/bpact/commit/c5614306a1ae2e632a31e8691683d8ccda5e58b3) Added `Matrix` to the `Imports` field in DESCRIPTION.
- [`970ce2b`](https://github.com/noahlorinczcomi/bpact/commit/970ce2b2de00bb6737a5de27fbcd0aa608af74fc) Updated README.

## 2025-04-30

- [`5df13ec`](https://github.com/noahlorinczcomi/bpact/commit/5df13ec52134d5a16aad8c1868213fa4c26586bf) Updated README.

## 2025-04-22

- [`8b966cd`](https://github.com/noahlorinczcomi/bpact/commit/8b966cd8b19cd78de0ddb7c63408d8115ea2b0b5) Updated README.
- [`c2d9832`](https://github.com/noahlorinczcomi/bpact/commit/c2d98322175754da8eab140321331d37ff3a1188) Added `flowchart.png` (package methodology diagram).
- [`f047704`](https://github.com/noahlorinczcomi/bpact/commit/f04770abbb8573b9c4a00ad41e13afcda998412f) Removed `flowchart.pdf`.
- [`570d607`](https://github.com/noahlorinczcomi/bpact/commit/570d607714e375a861d7063f4a9d510d3c3218ad) Added `flowchart.pdf` (package methodology diagram).

## 2025-04-15

- [`ce40faf`](https://github.com/noahlorinczcomi/bpact/commit/ce40fafcf23594f0d8cdaf2725e0753ec6042f7b) Updated README.

## 2025-04-02

- [`7633a64`](https://github.com/noahlorinczcomi/bpact/commit/7633a64573c65822c5914c1eee446ec9175e5187) Fixed a crash in `compositemh()` and `propshared()` when `max_size` exceeded the number of genes on a chromosome. `bigsnpr::snp_ldsplit()` now receives `max_size=min(c(nrow(Rho), max_size))` in both functions.

## 2025-03-18

- [`0f4fd08`](https://github.com/noahlorinczcomi/bpact/commit/0f4fd08818a5c7fcb860ceb013868bf48d791ff1) Updated README.

## 2025-03-17

- [`cf3d962`](https://github.com/noahlorinczcomi/bpact/commit/cf3d962ac0e81a330d05da688b8f01d43aa3ece3) Updated README.
- [`b8c385b`](https://github.com/noahlorinczcomi/bpact/commit/b8c385b92c41281be2003f2312af129bcf6d9a5a) Updated README.
- [`0078a15`](https://github.com/noahlorinczcomi/bpact/commit/0078a1515347f21255b9593b8c48495b7d37cefb) Updated README.
- [`a65c04d`](https://github.com/noahlorinczcomi/bpact/commit/a65c04db8475506702b9923c922e673c45160df6) Updated README.
- [`7475e75`](https://github.com/noahlorinczcomi/bpact/commit/7475e75020f55d1d9949e4bdea6a6d9483065b05) Removed `tutorial.R`.
- [`cdcb47a`](https://github.com/noahlorinczcomi/bpact/commit/cdcb47a23ded1b697b4f6df1ce0cec6f8d3552e7) Added example datasets `ad.gent` (Alzheimer's disease) and `lbd.gent` (Lewy body dementia) gene-based association test results, with documentation for both.
- [`1c7dcee`](https://github.com/noahlorinczcomi/bpact/commit/1c7dceeb3c68249f9801825d339a29593265bed2) Updated README.
- [`e03ef95`](https://github.com/noahlorinczcomi/bpact/commit/e03ef958de3a97e96b1c626468669c8b58891107) Added reference datasets `gent.Rho` (European 1000 Genomes gene-level test statistic correlation matrices) and `ld.df` (gene-level MAF-weighted LD scores), with documentation and DESCRIPTION dependency entries. Also added `tutorial.R` and the RStudio project file.

## 2025-03-07

- [`d73df47`](https://github.com/noahlorinczcomi/bpact/commit/d73df4731a54ea7c0af97033a0407fe77a90c92b) Added all core package source files: `R/functions.R` (871 lines) implementing `dtrap`, `ptrap`, `qtrap`, `rtrap`, `numerical_trap_distribution`, `mom_trap_distribution`, `mom_beta_distribution`, `mh`, `compositemh`, `posterior_gene`, `propshared`, and `shared_count_simex`; DESCRIPTION; NAMESPACE; and all man pages.
- [`dcb6fba`](https://github.com/noahlorinczcomi/bpact/commit/dcb6fbab1da660d295327598713ab83e30c82718) Initial commit (README placeholder).
