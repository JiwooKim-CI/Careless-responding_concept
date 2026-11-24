# A Causal Approach to Careless Responding

## Supplementary Material: Reproduction Code

This supplementary material accompanies the manuscript *A Causal Approach to Careless Responding*, submitted to *Methodology*.

The R script reproduces the two numerical illustrations presented in the manuscript.

## Numerical Illustrations

**Illustration 1: Response-Generating Mechanisms**

Examines mean bias under CRCR and CRNAR when all observed responses are retained or careless responses are deleted based on the true response state.

**Illustration 2: Detection Structures**

Examines selection bias and residual contamination following bogus-item, time-based, and pattern-based detection under CRCR. Each detection procedure is compared with random deletion at the same realized sensitivity and specificity.

## Requirements

The following R packages are required:

- dplyr
- tidyr
- knitr
- ggplot2
- ragg
- tikzDevice

A working LaTeX installation is required for TikZ figure generation and compilation.

## Reproducing the Results

Install the required packages and run `numerical_illustrations.R` from the working directory.

All data are generated within the script. No external datasets are required. The random seed is fixed at 123.

The script automatically generates replication-level results, summary tables, numerical verification checks, and the figures presented in the manuscript.

## Output Files

Results are saved in the `numerical_illustrations_results` directory.

| File | Description |
|---|---|
| stage1_generation.csv | Results for Illustration 1 |
| stage2_detection.csv | Bias decomposition for Illustration 2 |
| stage2_matched_structure.csv | Accuracy-matched random-deletion results |
| detector_accuracy.csv | Classification accuracy |
| generation_replications.csv | Replication-level generation results |
| detection_replications.csv | Replication-level detection results |
| all_results.rds | Complete numerical results |
| sessionInfo.txt | R and package version information |
| figures/ | Figure files and underlying data |

Additional summary tables and numerical verification results are also generated.

## Additional Information

Detailed simulation settings, parameter values, and numerical results are provided in the accompanying supplementary document.

The numerical illustrations demonstrate the theoretical implications of the proposed framework rather than establish the general superiority of particular detection procedures.

## Anonymization

Identifying information has been omitted for double-blind peer review.
