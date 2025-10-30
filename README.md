# FlexibleCountModel
This repository provides the implementation of a flexible statistical model designed to handle count data with varying levels of zero inflation and outliers.

**_To fit discrete extended generalized Pareto distribution (degpd) and zero-inflated discrete extended generalized Pareto distribution (zidegpd), we developed a code for new families and run it using evgam package fucntion_**

**_Intallation_**
1. Download the code and set the working directory. e.g., setwd("C:\Users\atouqeer\Downloads\degpd-and-zidegpd-main\degpd-and-zidegpd-main")
2. Call all the C++ and R functions using the following code

```markdown
dest <- "./R/"      # this function all R function 
files = list.files(dest, full.names = T)
for (i in 1:length(files)) {
  source(files[i])
}


dest <- "./src/"  
files = list.files(dest, full.names = T)
for (i in 1:length(files)) {
  Rcpp::sourceCpp(files[i])
}
```


