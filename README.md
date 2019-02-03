# seqExplorer
Integrated and interactive software for single-cell RNA-seq data analysis.
![](https://github.com/dyxmvp/Demos/blob/master/seq-Explorer/demo1.gif)   |![](https://github.com/dyxmvp/Demos/blob/master/seq-Explorer/demo2.gif)
-------------------------|-------------------------
![](https://github.com/dyxmvp/Demos/blob/master/seq-Explorer/demo3.gif)   |![](https://github.com/dyxmvp/Demos/blob/master/seq-Explorer/demo4.gif)


## Installation

These instructions will get you a copy of the software up and running on HPC or local machine.

### Prerequisites

```
R (>= 3.4)
Drop-seq_tools (>= 2.0)
```

### Installing

-> To install from **R** or **Rstudio**, run the following code:

(1) First, install package “devtools”:
```R
if (!require("devtools"))
  install.packages("devtools")
```
(2) Then, install seqExplorer:
```R
devtools::install_github('dyxmvp/seqExplorer')
```
-> To install from **HPC**, see this [guide](HPC_installation_guide.md) for more details.

## Getting Started

Enter commands in R or R studio:

```R
library(seqExplorer)
seqExplorer()
```

## License

This software is licensed under the GPL-3. See the [LICENSE](LICENSE) file for more details.
