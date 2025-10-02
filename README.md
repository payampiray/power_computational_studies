# power_computational_studies

Analysis code for the following study:  
Piray, Payam, "Addressing low statistical power in computational modeling studies in psychology and neuroscience", *Nature Human Behaviour*, 2025.

## Overview

This repository contains the MATLAB scripts used to perform the simulations and generate all figures presented in the paper.  
The code implements simulation-based power analyses tailored to computational modeling studies in psychology and neuroscience.

## Usage

To reproduce all figures from the paper, run the following command in MATLAB:

```matlab
main_figures()
```

This master script will sequentially call the corresponding scripts for each figure in the paper.

## Python Library

A companion Python package implementing the method introduced in the paper is available here:  
[cbm_power](https://github.com/payampiray/cbm_power)

The Python library provides:
- Tools for estimating statistical power given sample size and model space size.
- Automated sample-size optimization procedures.
