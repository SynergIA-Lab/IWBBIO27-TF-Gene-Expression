# Regulatory Biomarker Transfer Summary

## Data

| cohort | platform_class | endpoint | n_features | n_patients | n_events | censoring_fraction | analysis_status |
|---|---|---|---|---|---|---|---|
| TCGA_LUAD | rnaseq | PFI | 12650 | 504 | 207 | 0.589 | confirmatory_eligible |
| GSE31210 | microarray | RFS | 12650 | 204 | 54 | 0.735 | confirmatory_eligible |
| TCGA_LUAD | rnaseq | OS | 12650 | 504 | 182 | 0.639 | confirmatory_eligible |
| GSE72094 | microarray | OS | 12650 | 398 | 113 | 0.716 | confirmatory_eligible |
| TCGA_BRCA | rnaseq | PFI | 13100 | 1080 | 145 | 0.866 | confirmatory_eligible |
| GSE20685 | microarray | DMFS | 13100 | 319 | 75 | 0.765 | confirmatory_eligible |
| TCGA_BRCA | rnaseq | OS | 13100 | 1080 | 151 | 0.860 | confirmatory_eligible |
| GSE96058 | rnaseq | OS | 13100 | 3273 | 336 | 0.897 | confirmatory_eligible |
| TCGA_COAD | rnaseq | PFI | 12736 | 438 | 117 | 0.733 | confirmatory_eligible |
| GSE39582 | microarray | RFS | 12736 | 485 | 137 | 0.718 | confirmatory_eligible |
| TCGA_COAD | rnaseq | OS | 12736 | 439 | 98 | 0.777 | confirmatory_eligible |
| GSE17536 | microarray | OS | 12736 | 177 | 73 | 0.588 | confirmatory_eligible |

## Inference

| analysis | validation_cohort | n_patients | n_events | tau_days | model | primary_ipcw_cindex | lasso_sensitivity_ipcw_cindex | rsf_ipcw_cindex | primary_ci_lower | primary_ci_upper | primary_mean_dynamic_auc | primary_alpha_selection |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | gene_expression | 0.686 | 0.654 | 0.712 | 0.616 | 0.748 | 0.714 | cv |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | tf_expression | 0.653 | 0.648 | 0.667 | 0.570 | 0.734 | 0.674 | cv |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | collectri_ulm | 0.689 | 0.666 | 0.685 | 0.622 | 0.756 | 0.707 | cv |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | dorothea_viper | 0.701 | 0.691 | 0.689 | 0.630 | 0.764 | 0.705 | cv |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | clinical | 0.648 | nan | nan | 0.570 | 0.728 | nan | nan |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | clinical+gene_expression | 0.713 | nan | nan | 0.652 | 0.769 | nan | nan |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | clinical+tf_expression | 0.673 | nan | nan | 0.594 | 0.747 | nan | nan |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | clinical+collectri_ulm | 0.716 | nan | nan | 0.649 | 0.782 | nan | nan |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | clinical+dorothea_viper | 0.736 | nan | nan | 0.672 | 0.795 | nan | nan |
| LUAD_OS | GSE72094 | 398 | 113 | 1825.000 | gene_expression | 0.637 | 0.625 | 0.668 | 0.573 | 0.699 | 0.681 | cv |
| LUAD_OS | GSE72094 | 398 | 113 | 1825.000 | tf_expression | 0.608 | 0.596 | 0.633 | 0.536 | 0.680 | 0.675 | cv |
| LUAD_OS | GSE72094 | 398 | 113 | 1825.000 | collectri_ulm | 0.681 | 0.675 | 0.691 | 0.614 | 0.736 | 0.713 | cv |
| LUAD_OS | GSE72094 | 398 | 113 | 1825.000 | dorothea_viper | 0.671 | 0.661 | 0.656 | 0.611 | 0.728 | 0.703 | cv |
| LUAD_OS | GSE72094 | 393 | 111 | 1825.000 | clinical | 0.603 | nan | nan | 0.531 | 0.673 | nan | nan |
| LUAD_OS | GSE72094 | 393 | 111 | 1825.000 | clinical+gene_expression | 0.673 | nan | nan | 0.605 | 0.739 | nan | nan |
| LUAD_OS | GSE72094 | 393 | 111 | 1825.000 | clinical+tf_expression | 0.655 | nan | nan | 0.585 | 0.723 | nan | nan |
| LUAD_OS | GSE72094 | 393 | 111 | 1825.000 | clinical+collectri_ulm | 0.699 | nan | nan | 0.634 | 0.760 | nan | nan |
| LUAD_OS | GSE72094 | 393 | 111 | 1825.000 | clinical+dorothea_viper | 0.697 | nan | nan | 0.630 | 0.754 | nan | nan |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | gene_expression | 0.677 | 0.663 | 0.657 | 0.605 | 0.751 | 0.711 | cv |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | tf_expression | 0.617 | 0.579 | 0.634 | 0.545 | 0.687 | 0.625 | cv |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | collectri_ulm | 0.692 | 0.671 | 0.682 | 0.622 | 0.760 | 0.720 | cv |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | dorothea_viper | 0.685 | 0.679 | 0.654 | 0.619 | 0.756 | 0.716 | cv |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | clinical | 0.664 | nan | nan | 0.587 | 0.735 | nan | nan |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | clinical+gene_expression | 0.745 | nan | nan | 0.683 | 0.803 | nan | nan |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | clinical+tf_expression | 0.735 | nan | nan | 0.674 | 0.789 | nan | nan |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | clinical+collectri_ulm | 0.749 | nan | nan | 0.685 | 0.808 | nan | nan |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | clinical+dorothea_viper | 0.744 | nan | nan | 0.684 | 0.804 | nan | nan |
| BRCA_OS | GSE96058 | 3273 | 336 | 1825.000 | gene_expression | 0.622 | 0.616 | 0.627 | 0.587 | 0.655 | 0.642 | cv |
| BRCA_OS | GSE96058 | 3273 | 336 | 1825.000 | tf_expression | 0.604 | 0.594 | 0.602 | 0.570 | 0.635 | 0.616 | cv |
| BRCA_OS | GSE96058 | 3273 | 336 | 1825.000 | collectri_ulm | 0.622 | 0.609 | 0.610 | 0.589 | 0.655 | 0.655 | cv |
| BRCA_OS | GSE96058 | 3273 | 336 | 1825.000 | dorothea_viper | 0.612 | 0.611 | 0.573 | 0.580 | 0.642 | 0.634 | cv |
| BRCA_OS | GSE96058 | 3149 | 316 | 1825.000 | clinical | 0.729 | nan | nan | 0.697 | 0.761 | nan | nan |
| BRCA_OS | GSE96058 | 3149 | 316 | 1825.000 | clinical+gene_expression | 0.725 | nan | nan | 0.694 | 0.753 | nan | nan |
| BRCA_OS | GSE96058 | 3149 | 316 | 1825.000 | clinical+tf_expression | 0.719 | nan | nan | 0.688 | 0.747 | nan | nan |
| BRCA_OS | GSE96058 | 3149 | 316 | 1825.000 | clinical+collectri_ulm | 0.744 | nan | nan | 0.714 | 0.772 | nan | nan |
| BRCA_OS | GSE96058 | 3149 | 316 | 1825.000 | clinical+dorothea_viper | 0.731 | nan | nan | 0.699 | 0.760 | nan | nan |
| COAD_PFI | GSE39582 | 485 | 137 | 1825.000 | gene_expression | 0.592 | 0.581 | 0.578 | 0.544 | 0.639 | 0.612 | cv |
| COAD_PFI | GSE39582 | 485 | 137 | 1825.000 | tf_expression | 0.596 | 0.597 | 0.578 | 0.552 | 0.637 | 0.624 | cv |
| COAD_PFI | GSE39582 | 485 | 137 | 1825.000 | collectri_ulm | 0.596 | 0.591 | 0.601 | 0.543 | 0.640 | 0.615 | cv |
| COAD_PFI | GSE39582 | 485 | 137 | 1825.000 | dorothea_viper | 0.546 | 0.522 | 0.572 | 0.494 | 0.593 | 0.548 | cv |
| COAD_PFI | GSE39582 | 484 | 137 | 1825.000 | clinical | 0.623 | nan | nan | 0.581 | 0.666 | nan | nan |
| COAD_PFI | GSE39582 | 484 | 137 | 1825.000 | clinical+gene_expression | 0.646 | nan | nan | 0.605 | 0.694 | nan | nan |
| COAD_PFI | GSE39582 | 484 | 137 | 1825.000 | clinical+tf_expression | 0.646 | nan | nan | 0.607 | 0.691 | nan | nan |
| COAD_PFI | GSE39582 | 484 | 137 | 1825.000 | clinical+collectri_ulm | 0.640 | nan | nan | 0.595 | 0.688 | nan | nan |
| COAD_PFI | GSE39582 | 484 | 137 | 1825.000 | clinical+dorothea_viper | 0.630 | nan | nan | 0.586 | 0.674 | nan | nan |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | gene_expression | 0.510 | 0.516 | 0.497 | 0.439 | 0.585 | 0.518 | cv |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | tf_expression | 0.559 | 0.547 | 0.579 | 0.487 | 0.628 | 0.564 | cv |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | collectri_ulm | 0.551 | 0.541 | 0.534 | 0.488 | 0.624 | 0.584 | cv |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | dorothea_viper | 0.576 | 0.575 | 0.555 | 0.513 | 0.647 | 0.584 | cv |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | clinical | 0.761 | nan | nan | 0.707 | 0.818 | nan | nan |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | clinical+gene_expression | 0.711 | nan | nan | 0.644 | 0.774 | nan | nan |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | clinical+tf_expression | 0.725 | nan | nan | 0.660 | 0.788 | nan | nan |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | clinical+collectri_ulm | 0.748 | nan | nan | 0.685 | 0.803 | nan | nan |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | clinical+dorothea_viper | 0.745 | nan | nan | 0.687 | 0.804 | nan | nan |

Panel stability (mean pairwise Kuncheva over TCGA bootstrap resamples):

| analysis | representation | k | universe_size | n_resamples | kuncheva_mean |
|---|---|---|---|---|---|
| LUAD_PFI | gene_expression | 25 | 12650 | 50 | 0.065 |
| LUAD_PFI | tf_expression | 25 | 556 | 50 | 0.209 |
| LUAD_PFI | collectri_ulm | 25 | 556 | 50 | 0.249 |
| LUAD_PFI | dorothea_viper | 25 | 293 | 50 | 0.269 |
| LUAD_OS | gene_expression | 25 | 12650 | 50 | 0.164 |
| LUAD_OS | tf_expression | 25 | 556 | 50 | 0.321 |
| LUAD_OS | collectri_ulm | 25 | 556 | 50 | 0.323 |
| LUAD_OS | dorothea_viper | 25 | 293 | 50 | 0.400 |
| BRCA_PFI | gene_expression | 25 | 13100 | 50 | 0.083 |
| BRCA_PFI | tf_expression | 25 | 570 | 50 | 0.157 |
| BRCA_PFI | collectri_ulm | 25 | 570 | 50 | 0.208 |
| BRCA_PFI | dorothea_viper | 25 | 293 | 50 | 0.253 |
| BRCA_OS | gene_expression | 25 | 13100 | 50 | 0.099 |
| BRCA_OS | tf_expression | 25 | 570 | 50 | 0.159 |
| BRCA_OS | collectri_ulm | 25 | 570 | 50 | 0.228 |
| BRCA_OS | dorothea_viper | 25 | 293 | 50 | 0.274 |
| COAD_PFI | gene_expression | 25 | 12736 | 50 | 0.117 |
| COAD_PFI | tf_expression | 25 | 553 | 50 | 0.209 |
| COAD_PFI | collectri_ulm | 25 | 553 | 50 | 0.205 |
| COAD_PFI | dorothea_viper | 25 | 294 | 50 | 0.215 |
| COAD_OS | gene_expression | 25 | 12736 | 50 | 0.091 |
| COAD_OS | tf_expression | 25 | 553 | 50 | 0.202 |
| COAD_OS | collectri_ulm | 25 | 553 | 50 | 0.257 |
| COAD_OS | dorothea_viper | 25 | 294 | 50 | 0.202 |

## Primary Analysis (pooled, Holm)

Fixed-effect inverse-variance pooling of paired C-index differences; Holm across the pre-specified comparisons (role = primary). Per-cohort rows below are descriptive.

| comparison | scope | role | k_cohorts | fe_estimate | fe_ci_lower | fe_ci_upper | fe_p | fe_p_holm | re_estimate | re_p | i2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| collectri_ulm - gene_expression | all_cohorts | primary | 6 | 0.011 | -0.010 | 0.032 | 0.294 | 0.294 | 0.011 | 0.294 | 0.000 |
| collectri_ulm - gene_expression | endpoint_OS | sensitivity | 3 | 0.014 | -0.013 | 0.041 | 0.302 | nan | 0.016 | 0.290 | 0.108 |
| collectri_ulm - gene_expression | endpoint_progression | sensitivity | 3 | 0.007 | -0.026 | 0.039 | 0.689 | nan | 0.007 | 0.689 | 0.000 |
| collectri_ulm - tf_expression | all_cohorts | primary | 6 | 0.022 | 0.001 | 0.042 | 0.041 | 0.081 | 0.024 | 0.057 | 0.215 |
| collectri_ulm - tf_expression | endpoint_OS | sensitivity | 3 | 0.022 | -0.006 | 0.049 | 0.129 | nan | 0.023 | 0.211 | 0.264 |
| collectri_ulm - tf_expression | endpoint_progression | sensitivity | 3 | 0.022 | -0.009 | 0.052 | 0.170 | nan | 0.030 | 0.201 | 0.452 |

## Confirmatory Comparisons

Per-cohort results for the pre-specified comparisons (descriptive; inference is the pooled primary analysis).

| analysis | validation_cohort | n_patients | n_events | tau_days | comparison | comparison_status | endpoint_alignment | cindex_a | cindex_b | difference | ci_lower | ci_upper | bootstrap_se | p_two_sided | n_boot_valid | tau |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | collectri_ulm - gene_expression | confirmatory | progression-related endpoints paired | 0.689 | 0.686 | 0.004 | -0.057 | 0.061 | 0.031 | 0.982 | 1000 | 1825.000 |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | collectri_ulm - tf_expression | confirmatory | progression-related endpoints paired | 0.689 | 0.653 | 0.037 | -0.036 | 0.108 | 0.037 | 0.348 | 1000 | 1825.000 |
| LUAD_OS | GSE72094 | 398 | 113 | 1825.000 | collectri_ulm - gene_expression | confirmatory | OS paired with OS | 0.681 | 0.637 | 0.044 | -0.013 | 0.101 | 0.030 | 0.156 | 1000 | 1825.000 |
| LUAD_OS | GSE72094 | 398 | 113 | 1825.000 | collectri_ulm - tf_expression | confirmatory | OS paired with OS | 0.681 | 0.608 | 0.073 | 0.010 | 0.149 | 0.037 | 0.034 | 1000 | 1825.000 |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | collectri_ulm - gene_expression | confirmatory | progression-related endpoints paired | 0.692 | 0.677 | 0.015 | -0.050 | 0.076 | 0.033 | 0.676 | 1000 | 1825.000 |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | collectri_ulm - tf_expression | confirmatory | progression-related endpoints paired | 0.692 | 0.617 | 0.075 | 0.007 | 0.145 | 0.035 | 0.030 | 1000 | 1825.000 |
| BRCA_OS | GSE96058 | 3273 | 336 | 1825.000 | collectri_ulm - gene_expression | confirmatory | OS paired with OS | 0.622 | 0.622 | -0.000 | -0.031 | 0.033 | 0.017 | 0.988 | 1000 | 1825.000 |
| BRCA_OS | GSE96058 | 3273 | 336 | 1825.000 | collectri_ulm - tf_expression | confirmatory | OS paired with OS | 0.622 | 0.604 | 0.018 | -0.015 | 0.051 | 0.017 | 0.296 | 1000 | 1825.000 |
| COAD_PFI | GSE39582 | 485 | 137 | 1825.000 | collectri_ulm - gene_expression | confirmatory | progression-related endpoints paired | 0.596 | 0.592 | 0.004 | -0.042 | 0.053 | 0.024 | 0.872 | 1000 | 1825.000 |
| COAD_PFI | GSE39582 | 485 | 137 | 1825.000 | collectri_ulm - tf_expression | confirmatory | progression-related endpoints paired | 0.596 | 0.596 | -0.000 | -0.041 | 0.038 | 0.020 | 0.972 | 1000 | 1825.000 |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | collectri_ulm - gene_expression | confirmatory | OS paired with OS | 0.551 | 0.510 | 0.041 | -0.035 | 0.116 | 0.039 | 0.302 | 1000 | 1825.000 |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | collectri_ulm - tf_expression | confirmatory | OS paired with OS | 0.551 | 0.559 | -0.009 | -0.078 | 0.057 | 0.035 | 0.844 | 1000 | 1825.000 |

## Exploratory Analyses

| analysis | validation_cohort | n_patients | n_events | tau_days | comparison | comparison_status | endpoint_alignment | cindex_a | cindex_b | difference | ci_lower | ci_upper | bootstrap_se | p_two_sided | n_boot_valid | tau |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | gene_expression - tf_expression | exploratory | progression-related endpoints paired | 0.686 | 0.653 | 0.033 | -0.034 | 0.101 | 0.035 | 0.328 | 1000 | 1825.000 |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | gene_expression - dorothea_viper | exploratory | progression-related endpoints paired | 0.686 | 0.701 | -0.016 | -0.074 | 0.045 | 0.031 | 0.612 | 1000 | 1825.000 |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | tf_expression - dorothea_viper | exploratory | progression-related endpoints paired | 0.653 | 0.701 | -0.049 | -0.120 | 0.026 | 0.038 | 0.170 | 1000 | 1825.000 |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | collectri_ulm - dorothea_viper | exploratory | progression-related endpoints paired | 0.689 | 0.701 | -0.012 | -0.064 | 0.044 | 0.028 | 0.612 | 1000 | 1825.000 |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | clinical+gene_expression - clinical | exploratory | progression-related endpoints paired | 0.713 | 0.648 | 0.065 | -0.021 | 0.145 | 0.043 | 0.152 | 1000 | 1825.000 |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | clinical+tf_expression - clinical | exploratory | progression-related endpoints paired | 0.673 | 0.648 | 0.025 | -0.060 | 0.115 | 0.045 | 0.604 | 1000 | 1825.000 |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | clinical+collectri_ulm - clinical | exploratory | progression-related endpoints paired | 0.716 | 0.648 | 0.068 | -0.019 | 0.153 | 0.043 | 0.130 | 1000 | 1825.000 |
| LUAD_PFI | GSE31210 | 204 | 54 | 1825.000 | clinical+dorothea_viper - clinical | exploratory | progression-related endpoints paired | 0.736 | 0.648 | 0.088 | 0.001 | 0.169 | 0.042 | 0.046 | 1000 | 1825.000 |
| LUAD_OS | GSE72094 | 398 | 113 | 1825.000 | gene_expression - tf_expression | exploratory | OS paired with OS | 0.637 | 0.608 | 0.029 | -0.025 | 0.094 | 0.029 | 0.332 | 1000 | 1825.000 |
| LUAD_OS | GSE72094 | 398 | 113 | 1825.000 | gene_expression - dorothea_viper | exploratory | OS paired with OS | 0.637 | 0.671 | -0.034 | -0.094 | 0.025 | 0.030 | 0.270 | 1000 | 1825.000 |
| LUAD_OS | GSE72094 | 398 | 113 | 1825.000 | tf_expression - dorothea_viper | exploratory | OS paired with OS | 0.608 | 0.671 | -0.063 | -0.139 | -0.000 | 0.036 | 0.050 | 1000 | 1825.000 |
| LUAD_OS | GSE72094 | 398 | 113 | 1825.000 | collectri_ulm - dorothea_viper | exploratory | OS paired with OS | 0.681 | 0.671 | 0.010 | -0.028 | 0.048 | 0.020 | 0.636 | 1000 | 1825.000 |
| LUAD_OS | GSE72094 | 393 | 111 | 1825.000 | clinical+gene_expression - clinical | exploratory | OS paired with OS | 0.673 | 0.603 | 0.070 | 0.015 | 0.135 | 0.030 | 0.008 | 1000 | 1825.000 |
| LUAD_OS | GSE72094 | 393 | 111 | 1825.000 | clinical+tf_expression - clinical | exploratory | OS paired with OS | 0.655 | 0.603 | 0.052 | 0.011 | 0.099 | 0.023 | 0.012 | 1000 | 1825.000 |
| LUAD_OS | GSE72094 | 393 | 111 | 1825.000 | clinical+collectri_ulm - clinical | exploratory | OS paired with OS | 0.699 | 0.603 | 0.096 | 0.042 | 0.157 | 0.029 | 0.000 | 1000 | 1825.000 |
| LUAD_OS | GSE72094 | 393 | 111 | 1825.000 | clinical+dorothea_viper - clinical | exploratory | OS paired with OS | 0.697 | 0.603 | 0.094 | 0.044 | 0.154 | 0.028 | 0.000 | 1000 | 1825.000 |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | gene_expression - tf_expression | exploratory | progression-related endpoints paired | 0.677 | 0.617 | 0.060 | -0.018 | 0.133 | 0.039 | 0.128 | 1000 | 1825.000 |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | gene_expression - dorothea_viper | exploratory | progression-related endpoints paired | 0.677 | 0.685 | -0.007 | -0.068 | 0.052 | 0.032 | 0.798 | 1000 | 1825.000 |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | tf_expression - dorothea_viper | exploratory | progression-related endpoints paired | 0.617 | 0.685 | -0.067 | -0.142 | 0.002 | 0.037 | 0.056 | 1000 | 1825.000 |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | collectri_ulm - dorothea_viper | exploratory | progression-related endpoints paired | 0.692 | 0.685 | 0.007 | -0.051 | 0.057 | 0.028 | 0.806 | 1000 | 1825.000 |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | clinical+gene_expression - clinical | exploratory | progression-related endpoints paired | 0.745 | 0.664 | 0.081 | 0.029 | 0.136 | 0.028 | 0.002 | 1000 | 1825.000 |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | clinical+tf_expression - clinical | exploratory | progression-related endpoints paired | 0.735 | 0.664 | 0.071 | 0.023 | 0.120 | 0.025 | 0.000 | 1000 | 1825.000 |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | clinical+collectri_ulm - clinical | exploratory | progression-related endpoints paired | 0.749 | 0.664 | 0.085 | 0.041 | 0.129 | 0.022 | 0.000 | 1000 | 1825.000 |
| BRCA_PFI | GSE20685 | 319 | 75 | 1825.000 | clinical+dorothea_viper - clinical | exploratory | progression-related endpoints paired | 0.744 | 0.664 | 0.080 | 0.042 | 0.124 | 0.021 | 0.000 | 1000 | 1825.000 |
| BRCA_OS | GSE96058 | 3273 | 336 | 1825.000 | gene_expression - tf_expression | exploratory | OS paired with OS | 0.622 | 0.604 | 0.018 | -0.013 | 0.050 | 0.016 | 0.264 | 1000 | 1825.000 |
| BRCA_OS | GSE96058 | 3273 | 336 | 1825.000 | gene_expression - dorothea_viper | exploratory | OS paired with OS | 0.622 | 0.612 | 0.010 | -0.022 | 0.042 | 0.017 | 0.550 | 1000 | 1825.000 |
| BRCA_OS | GSE96058 | 3273 | 336 | 1825.000 | tf_expression - dorothea_viper | exploratory | OS paired with OS | 0.604 | 0.612 | -0.009 | -0.041 | 0.021 | 0.016 | 0.636 | 1000 | 1825.000 |
| BRCA_OS | GSE96058 | 3273 | 336 | 1825.000 | collectri_ulm - dorothea_viper | exploratory | OS paired with OS | 0.622 | 0.612 | 0.009 | -0.022 | 0.041 | 0.016 | 0.554 | 1000 | 1825.000 |
| BRCA_OS | GSE96058 | 3149 | 316 | 1825.000 | clinical+gene_expression - clinical | exploratory | OS paired with OS | 0.725 | 0.729 | -0.004 | -0.021 | 0.013 | 0.009 | 0.602 | 1000 | 1825.000 |
| BRCA_OS | GSE96058 | 3149 | 316 | 1825.000 | clinical+tf_expression - clinical | exploratory | OS paired with OS | 0.719 | 0.729 | -0.010 | -0.032 | 0.012 | 0.011 | 0.352 | 1000 | 1825.000 |
| BRCA_OS | GSE96058 | 3149 | 316 | 1825.000 | clinical+collectri_ulm - clinical | exploratory | OS paired with OS | 0.744 | 0.729 | 0.014 | 0.004 | 0.026 | 0.006 | 0.016 | 1000 | 1825.000 |
| BRCA_OS | GSE96058 | 3149 | 316 | 1825.000 | clinical+dorothea_viper - clinical | exploratory | OS paired with OS | 0.731 | 0.729 | 0.002 | -0.009 | 0.013 | 0.006 | 0.748 | 1000 | 1825.000 |
| COAD_PFI | GSE39582 | 485 | 137 | 1825.000 | gene_expression - tf_expression | exploratory | progression-related endpoints paired | 0.592 | 0.596 | -0.004 | -0.047 | 0.035 | 0.021 | 0.856 | 1000 | 1825.000 |
| COAD_PFI | GSE39582 | 485 | 137 | 1825.000 | gene_expression - dorothea_viper | exploratory | progression-related endpoints paired | 0.592 | 0.546 | 0.046 | 0.003 | 0.093 | 0.022 | 0.038 | 1000 | 1825.000 |
| COAD_PFI | GSE39582 | 485 | 137 | 1825.000 | tf_expression - dorothea_viper | exploratory | progression-related endpoints paired | 0.596 | 0.546 | 0.050 | 0.006 | 0.103 | 0.024 | 0.032 | 1000 | 1825.000 |
| COAD_PFI | GSE39582 | 485 | 137 | 1825.000 | collectri_ulm - dorothea_viper | exploratory | progression-related endpoints paired | 0.596 | 0.546 | 0.050 | 0.001 | 0.099 | 0.026 | 0.050 | 1000 | 1825.000 |
| COAD_PFI | GSE39582 | 484 | 137 | 1825.000 | clinical+gene_expression - clinical | exploratory | progression-related endpoints paired | 0.646 | 0.623 | 0.023 | -0.011 | 0.059 | 0.018 | 0.174 | 1000 | 1825.000 |
| COAD_PFI | GSE39582 | 484 | 137 | 1825.000 | clinical+tf_expression - clinical | exploratory | progression-related endpoints paired | 0.646 | 0.623 | 0.023 | -0.008 | 0.055 | 0.016 | 0.138 | 1000 | 1825.000 |
| COAD_PFI | GSE39582 | 484 | 137 | 1825.000 | clinical+collectri_ulm - clinical | exploratory | progression-related endpoints paired | 0.640 | 0.623 | 0.017 | -0.014 | 0.049 | 0.016 | 0.276 | 1000 | 1825.000 |
| COAD_PFI | GSE39582 | 484 | 137 | 1825.000 | clinical+dorothea_viper - clinical | exploratory | progression-related endpoints paired | 0.630 | 0.623 | 0.007 | -0.021 | 0.036 | 0.015 | 0.644 | 1000 | 1825.000 |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | gene_expression - tf_expression | exploratory | OS paired with OS | 0.510 | 0.559 | -0.050 | -0.125 | 0.025 | 0.038 | 0.184 | 1000 | 1825.000 |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | gene_expression - dorothea_viper | exploratory | OS paired with OS | 0.510 | 0.576 | -0.066 | -0.144 | 0.011 | 0.039 | 0.084 | 1000 | 1825.000 |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | tf_expression - dorothea_viper | exploratory | OS paired with OS | 0.559 | 0.576 | -0.016 | -0.097 | 0.060 | 0.040 | 0.648 | 1000 | 1825.000 |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | collectri_ulm - dorothea_viper | exploratory | OS paired with OS | 0.551 | 0.576 | -0.025 | -0.093 | 0.041 | 0.035 | 0.424 | 1000 | 1825.000 |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | clinical+gene_expression - clinical | exploratory | OS paired with OS | 0.711 | 0.761 | -0.050 | -0.085 | -0.014 | 0.018 | 0.012 | 1000 | 1825.000 |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | clinical+tf_expression - clinical | exploratory | OS paired with OS | 0.725 | 0.761 | -0.035 | -0.070 | -0.002 | 0.018 | 0.038 | 1000 | 1825.000 |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | clinical+collectri_ulm - clinical | exploratory | OS paired with OS | 0.748 | 0.761 | -0.013 | -0.033 | 0.009 | 0.011 | 0.246 | 1000 | 1825.000 |
| COAD_OS | GSE17536 | 177 | 73 | 1825.000 | clinical+dorothea_viper - clinical | exploratory | OS paired with OS | 0.745 | 0.761 | -0.015 | -0.043 | 0.015 | 0.015 | 0.354 | 1000 | 1825.000 |

## Interpretation

No results for this section.

## Limitations

Cohort-wise z-scoring, endpoint mismatches, low event counts, ARACNe leakage flags and any 'fallback' penalty selection must be reviewed before making claims.
