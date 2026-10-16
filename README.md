## Purpose

This is the code and data repository for the manuscript:

RM Flight, PT Thompson, HNB Moseley; "Evaluation of the quality of deposited datasets in the Metabolomics Workbench data repository." *Draft*.

## Results

Current results are available on the Zenodo repository: [10.5281/zenodo.23269018](https://doi.org/10.5281/zenodo.23269018)

The data available there include:

* A CSV file with all of the metadata extracted from each analysis, with the overall color grade, as well as individual metrics and color grades, with the analysis ID, study ID, and project ID for each one.
* A zip file containing quality-control and quality-assurance reports for each analysis.
* A zip file containing R serialized rds files for each analysis.

## Presentation

The slides used for the presentation at the International Electronic Conference on Metabolomics 2026 (IECM2026) are available [here](outputs/mw_qcq_iecm2026_presentation_2026-10-16.pdf).

## Overview

Here is an overview figure of the workflow that we are running.

Currently, we process 7240 datasets using 7000 lines of R code organized in 230 functions, all orchestrated using the `targets` package.

![Pipeline overview.](resources/qcqa_pipeline_lucid_alt2.png)