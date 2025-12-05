# Performance-of-ThyroScan-in-the-Diagnosis-of-Thyroid-Cancer-in-a-Chinese-Population

This code repository contains statistical analysis code for the paper "Diagnostic Performance of a Multigene Genomic Classifier for Thyroid Cancer in a Chinese Population: A Prospective, Blinded, Multicenter Study". It is primarily used to evaluate the diagnostic performance of the ThyroScan multigene genomic classifier.

## File Description

- `code_for_result_and_plot.r` - Main R code file containing statistical analysis, ROC curve plotting, and performance metrics calculation

## Overview

### Figure S1: ROC Curve Analysis
- Plot ROC curves using the `pROC` package
- Calculate AUC (Area Under the Curve) with 95% confidence intervals
- Generate ROC curve plots

### Diagnostic Performance Metrics Calculation
- Calculate accuracy, sensitivity, specificity, positive predictive value (PPV), and negative predictive value (NPV)
- Compute 95% confidence intervals for all metrics using Wilson's method
- Output formatted results

### Figure 2: PPV/NPV Analysis Across Prevalence Rates
- Predict the performance of the multigene genomic classifier across varying cancer prevalence rates based on given sensitivity and specificity
- Calculate 95% confidence intervals using the Clopper-Pearson exact method
- Plot curves showing PPV and NPV changes with prevalence rates

## Requirements

### R Package Dependencies
- `pROC` - ROC curve analysis
- `PropCIs` - Confidence interval calculation

Install required packages:
```r
install.packages(c("pROC", "PropCIs"))
```

### Data Requirements
The code requires the following variables (which are provided and replacable):
- `data$True_Label` - True labels (0/1)
- `data$Prediction_Score` - Prediction scores
- `tp`, `fp`, `fn`, `tn` - Confusion matrix elements
- `Se`, `Sp` - Sensitivity and specificity
