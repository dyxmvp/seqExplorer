# Installation guide on HPC

This is an installation guide of seqExplorer for [HPC at Yale University](https://research.computing.yale.edu/services/high-performance-computing), but it should also work for other HPC. This guide shows how to set up seqExplorer on [Ruddle](https://research.computing.yale.edu/support/hpc/clusters/ruddle), [Grace](https://research.computing.yale.edu/support/hpc/clusters/grace), and [Farnam](https://research.computing.yale.edu/support/hpc/clusters/farnam). The installation takes about **30 ~ 40 minutes**.

### 1. Installation
Run the following commands:
```sh
srun --x11 --pty -p interactive -c 1 --mem=10g bash

module load R/3.5.0-foss-2016b-avx2
```
If you are using **Grace cluster**, add the following command:
```sh
module load Libs/HDF5/1.8.13-intel15
```
Then, run the following commands:
```sh
R

if (!require("devtools"))
  install.packages("devtools")

devtools::install_github('dyxmvp/seqExplorer')
```
Type **"yes"** and enter, if the system asks:

(1)
```sh
Would you like to use a personal library instead? (yes/No/cancel)
```
(2)
```sh
Would you like to create a personal library? (yes/No/cancel)
```
Select a **CRAN mirror** (for example: USA(CA 1)[https]), if the system asks:
```
--- Please select a CRAN mirror for the use in this session ---
```
After finished, **exit** R program using:
```R
q()
```

### 2. Getting Started
Run the following commands:
```sh
srun --x11 --pty -p interactive -c 5 --mem=40g bash

module load R/3.5.0-foss-2016b-avx2
```
If you are using **Grace cluster**, add the following command:
```sh
module load Libs/HDF5/1.8.13-intel15
```
Then, run the following commands:
```sh
R

library(seqExplorer)

seqExplorer()
```
Finally, **seqExplorer will launch with a browser.**

### 3. Exit
To close seqExplorer, use **ctrl+c** to return to R or use **ctrl+z** to close all R program. 
