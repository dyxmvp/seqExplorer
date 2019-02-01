# seqExplorer
Integrated and interactive software for single-cell RNA-seq data analysis

## Installation

These instructions will get you a copy of the software up and running on HPC or local machine.

### Prerequisites

```
R (>= 3.4)
Drop-seq_tools (>= 2.0)
Seurat (>= 2.3)
```

### Installing

-> To install from **R** or **Rstudio**, run the following code:

```R
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github('dyxmvp/seqExplorer')
```
-> To install from **HPC**, see this **[guide](HPC installing guide)** for more details.

## Getting Started

Enter commands in R or R studio:

```R
library(seqExplorer)
seqExplorer()
```

## License

This software is licensed under the GPL-3. See the [LICENSE](LICENSE) file for more details.
