# seq-Explorer
Integrated and interactive software for single-cell RNA-seq data analysis.

![](https://github.com/dyxmvp/Demos/blob/master/seq-Explorer/demo1.gif)   |![](https://github.com/dyxmvp/Demos/blob/master/seq-Explorer/demo2.gif)
:-------------------------:|:-------------------------:
![](https://github.com/dyxmvp/Demos/blob/master/seq-Explorer/demo3.gif)   |![](https://github.com/dyxmvp/Demos/blob/master/seq-Explorer/demo4.gif)


## 1. Installation

These instructions will get you a copy of the software up and running on HPC or local machine.

### Prerequisites

```
R (>= 3.4)
Drop-seq_tools (>= 2.0)
Seurat (>= 2.3)
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

## 2. Getting Started

Enter commands in R or R studio:

```R
library(seqExplorer)
seqExplorer()
```

## 3. License

seq-Explorer is developed and maintained by Yanxiang Deng@[Rong Fan Lab](https://www.eng.yale.edu/fanlab/Rong_Fan_Group/Welcome.html), BME, Yale University. The software is licensed under the GPL-3.0. See the [LICENSE](LICENSE) file for more details.

## 4. Contact

seq-Explorer is hosted at [Github](https://github.com/dyxmvp/seqExplorer). Please e-mail [yanxiang.deng@yale.edu](mailto:yanxiang.deng@yale.edu) or submit pull requests from the [Github](https://github.com/dyxmvp/seqExplorer) with any questions about the software.
