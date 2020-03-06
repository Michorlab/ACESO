# ACESO

## To install ACESO:
```{r}
install.packages("devtools")  
library(devtools)  
install_github("Michorlab/ACESO") 
```
## mrgsolve package
We recommend the installation of the development version of mrgsolve package from the metrumresearchgroup github repository:
```{r}
devtools::install_github("metrumresearchgroup/mrgsolve")
```
Additionally, a C/C++ compiler is needed as mrgsolves uses code in C to speed up the simulations.
### Windows:
For Windows users, Rtools is needed to be able to use mrgsolve package. It can be downloaded from https://cran.r-project.org/bin/windows/Rtools/. Select the Rtools download compatible with your R version and check your PATH environment variable in R during the installation process.

### Mac OS X:

Install C++ and gfortran compilers found here: https://cran.r-project.org/bin/macosx/tools/

### Linux: 
For Debian/Ubuntu, you can install the core software development utilities required for R package development by executing:  

sudo apt-get install r-base-dev texlive-full  

Some packages may require installation of additional R build dependencies. To provide all components needed to build R itself from source you can execute:  

sudo apt-get build-dep r-base-core  

## Final Test
A good test to check that everything is set up correctly is to use the devtools package. In R:
```{r}
install.packages("devtools")
devtools::has_devel()
```
If returns TRUE, you are good to go with ACESO.

# Vignettes
In order to learn how to use ACESO, 4 different vignettes has been created. To check their names:
```{r}
vignettes(package="ACESO")
```
To select one of them:
```{r}
vignettes(package="ACESO","Two_type_branching_process")
```
