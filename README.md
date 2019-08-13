## Umediation
The Umediation R package enables the user to simulate unmeasured confounding in mediation analysis in order to see how the results of the mediation analysis would change in the presence of unmeasured confounding.

## Installation
Requirements:
* R v3.4 or higher
* You will need the proper compiling tools for your platform.
  * For Windows (Rtools installer): https://cran.r-project.org/bin/windows/Rtools/
  * For MacOSX (clang and gfortran): https://cran.r-project.org/bin/macosx/tools/
```
# The devtools package must be installed first
install.packages("devtools") 

# these will fail to install when already loaded, and install_github will sometimes 
# load these as part of its activity, and will then try to install them if they need 
# an update for one of the package dependencies
install.packages(c("Rcpp","RcppEigen", "curl"), quiet=T) 
# if there is an error you might need to set quiet=F to see more information about it

#finally install UMediationThread and any remaining dependencies
devtools::install_github("SharonLutz/UmediationThread",quiet=T)
```

The install process will involve compiling source code. If you are on MacOSX, and do not quiet the build process with quiet=T, this may involve the clang compiler issuing warnings about unknown pragmas similar to the text below. Do not be alarmed if you see these. If there is actually an error, it will be present among the last several messages issued by the compiler.

```
warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
```
## Example
Below, we simulate 4 unmeasured confounders U (2 normally distributed and 2 Bernouilli distributed random variables) on the binary exposure, A, normally distributed mediator, M, and normally distributed outcome Y adjusted for one normally distributed covariate and 2 binary distributed covariates.
```


library(UmediationThread)
?Umediation # For details on this function and how to choose input variables

testM<- Umediation(n=1000,Atype="D",Mtype="C",Ytype="C",Ctype=c("C","D","D"),Utype=c("C","D","D","C"),
interact=TRUE,muC=c(0.1,0.3,0.2),varC=c(1,1,1),muU=c(.1,0.3,0.2,.1),varU=c(1,1,1,1),gamma0=0,
gammaC=c(1,0.3,0.2),gammaU=c(1,0.3,0.2,0.4),varA=1,alpha0=0,alphaA=1,alphaC=c(0.3,0.2,0.2),
alphaU=c(0.3,0.2,0.3,0.2),varM=1,beta0=0,betaA=-1,betaM=1,betaI=1,betaC=c(0.3,0.2,0.1),
betaU=c(0.3,0.2,-1.3,0.2),varY=1,alpha=0.05,nSim=100,nBoot=400)

testM

```
## Speed up mediation with Threading and Eigen via C++
the Umediation command accepts the following parameters:
* use_cpp, a boolean(T, F, True, or False), which activates the use of Rcpp RcppEigen, and threading.
* num_cores, an integer specifying the number of cpus/threads you wish to use.
```
Example Using Rcpp with Eigen and 5 threads:

testM<- Umediation(n=1000,Atype="D",Mtype="C",Ytype="C",Ctype=c("C","D","D"),Utype=c("C","D","D","C"),
interact=TRUE,muC=c(0.1,0.3,0.2),varC=c(1,1,1),muU=c(.1,0.3,0.2,.1),varU=c(1,1,1,1),gamma0=0,
gammaC=c(1,0.3,0.2),gammaU=c(1,0.3,0.2,0.4),varA=1,alpha0=0,alphaA=1,alphaC=c(0.3,0.2,0.2),
alphaU=c(0.3,0.2,0.3,0.2),varM=1,beta0=0,betaA=-1,betaM=1,betaI=1,betaC=c(0.3,0.2,0.1),
betaU=c(0.3,0.2,-1.3,0.2),varY=1,alpha=0.05,nSim=100,nBoot=400, use_cpp = T, num_cores = 5)

```
### Important Considerations for Threading:

#### Determining Number of Cores and amount of RAM available on your system:
For Windows:
Opening the start menu and searching for System Information, or running %windir%\system32\msinfo32.exe at the command prompt should open a window giving a system summary. The Row beginning with "Processor" will list your CPU's name and describe the number of physical and logical cores. For our purposes, the number of logical cores is the value for consideration.

For Mac OSX:

From a terminal: 
```
sysctl hw.physicalcpu hw.logicalcpu #should show the number of physical and logical cores.
```

For linux:
There is a system file '/proc/cpuinfo' which contains information about the processors available on the system. There should be one entry per logical core.

#### What to do with this information?
It is advisable in a standard use case that you tailor your num_jobs variable to be the # of your logical CPU cores - 1, at the maximum, to leave 1 core free to handle the original calling R process and any background OS processes. If you use more procesess or threads than this, your system may become slow and relatively unresponsive and may lock up until the processing completes, and it will likely not benefit from the additional threads in terms of speed improvement. It may actually run slower due to overconsumption of system resources.

## Output
For this analysis, we can see that there is not a significant difference in the proportion of simulations for the mediated effect if the unmeasured confounders are included, but there is a large diffence in the inference for the direct effect if these unmeasured confounders are not included in the analysis.


```
$Results
                                                                               [,1]
Prop. of simulations w/ significant ACME excluding U                    1.000000000
Prop. of simulations w/ significant ACME including U                    1.000000000
Prop. of simulations where conclusions based on ACME match              1.000000000
Average ACME excluding U                                                2.061403587
Average ACME including U                                                1.494588641
Average absolute difference of ACME including U minus ACME excluding U  0.566814946
Prop. of simulations w/ significant ADE excluding U                     0.030000000
Prop. of simulations w/ significant ADE including U                     0.450000000
Prop. of simulations where conclusions based on ADE match               0.520000000
Average ADE excluding U                                                 0.006946501
Average ADE including U                                                -0.185434050
Average absolute difference of ADE including U minus ADE excluding U    0.192380551

$Correlations_Between_Variables
      A    M     Y    C1    C2    C3    U1    U2    U3    U4
A  1.00 0.56  0.49  0.31  0.04  0.07  0.37  0.05  0.03  0.16
M  0.56 1.00  0.87  0.33  0.08  0.09  0.35  0.11  0.15  0.26
Y  0.49 0.87  1.00  0.37  0.09  0.08  0.38  0.11 -0.08  0.27
C1 0.31 0.33  0.37  1.00  0.04  0.01 -0.03 -0.01  0.00 -0.01
C2 0.04 0.08  0.09  0.04  1.00  0.03 -0.03  0.02  0.00  0.02
C3 0.07 0.09  0.08  0.01  0.03  1.00  0.01  0.01  0.02 -0.01
U1 0.37 0.35  0.38 -0.03 -0.03  0.01  1.00 -0.01  0.00 -0.04
U2 0.05 0.11  0.11 -0.01  0.02  0.01 -0.01  1.00  0.01 -0.02
U3 0.03 0.15 -0.08  0.00  0.00  0.02  0.00  0.01  1.00  0.04
U4 0.16 0.26  0.27 -0.01  0.02 -0.01 -0.04 -0.02  0.04  1.00

$Warning
[1] "Warning: correlations are only valid if at least one of the variables is normally distributed."

```

## Warning: Do not try to access package internals directly or do so at your own risk!
If you try to run methods/functions that are not exported and intended for end users, and feed these functions environments, parameters, or values that are not correctly formed, it could result in an uncaught or uncatchable C++ exception or segmentation fault. If this occurs, it will kill your R session/terminal and if you were working within RStudio it will probably crash too.

## Reference
**Lutz SM**, Thwing A, Schmiege S, Kroehl M, Baker CD, Starling AP, Hokanson JE, Ghosh D. (2017) Examining the Role of Unmeasured Confounding in Mediation Analysis with Genetic and Genomic Applications. *BMC Bioinformatics.* 18(1):344.

