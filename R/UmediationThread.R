#' @import car
#' @import mediation
#' @import pbapply
#' @include generate_data_matrix.R
#' @include mediate_parallel.R
#' @include perform_mediation.R

#' @export
#' @title UmediationThread
#' @name UmediationThread
#' @aliases UmediationThread
#' @description The Umediation function examines the role of unmeasured confounding on the estimates for the average causal mediated effect (ACME) and average direct effect (ADE) in mediation analysis. User input specifies the relationship between the exposure A, the mediator M, the outcome Y, the measured confounder C, and the unmeasured confounder U. The function runs mediation analysis including and excluding the unmeasured confounder U in the model. Umediation allows the user to examine how the results of the mediation analysis would change if the unmeasured confounder U was included or not included in the model. The function allows for continuous or dichotomous exposure A, mediator M, outcome Y, measured confounder C, and unmeasured confounder U. Umediation allows for multiple measured confounders C and unmeasured confounders U. In addition, Umediation allows for an interaction between the exposure A and the mediator M on the outcome Y.
#' @author Sharon Lutz, Michael Gooch
#' @param n is the sample size of the population that is being simulated.
#' @param Atype is either "C" for continuous, normally distributed exposure A or "D" for dichotomous, binary exposure A (i.e. Atype="D").
#' @param Mtype is either "C" continuous, normally distributed mediator M or "D" for dichotomous, binary mediator M (i.e. Mtype="C").
#' @param Ytype is either "C" for continuous, normally distributed outcome Y or "D" for dichotomous, binary outcome Y (i.e. Ytype="D").
#' @param Ctype is either "C" for continuous, normally distributed measured confounder C or "D" for dichotomous, binary measured confounder C. Ctype can be a single value (i.e. Ctype="C") or a vector for multiple measured confounders (i.e. Ctype=c("C","C","D") ).
#' @param Utype is either "C" for continuous, normally distributed unmeasured confounder U or "D" for dichotomous, binary unmeasured confounder U. Utype can be a single value (i.e. Utype="C") or a vector for multiple measured confounders (i.e. Utype=c("C","C","D") ).
#' @param interact Using  the flag interact=TRUE allows for an interaction between the exposure A and the mediator M on the outcome Y (i.e. E[Y]=beta0+betaA*A+betaM*M+betaC*C+betaU*U+betaI*A*M). By default, interact=FALSE (i.e. E[Y]=beta0+betaA*A+betaM*M+betaC*C+betaU*U).
#' @param muC is the mean vector for the measured confounder C. For continuous measured confounder (i.e. Ctype="C"), muC is the mean of C. For dichotomous measured confounder C (i.e. Ctype="D"), muC is the probability C=1.
#' @param varC is the variance of the measured confounder C when Ctype="C". For multiple measured confounders, the length of varC must match muC and Ctype (i.e. Ctype=c("C","C","D") and muC=c(-0.1,0.2,0.3) and varC=c(1,1,1))
#' @param muU is the mean vector for the unmeasured confounder U. For continuous unmeasured confounder (i.e. Utype="C"), muU is the mean of U. For dichotomous unmeasured confounder U (i.e. Utype="D"), muU is the probability U=1.
#' @param varU is the variance of the unmeasured confounder U when Utype="C". For multiple unmeasured confounders, the length of varU must match muU and Utype (i.e. Utype=c("C","C","D") and muU=c(-0.1,0.2,0.3) and varU=c(1,1,1))
#' @param gamma0 specifies the intercept for the exposure A (i.e. logit(P(A=1)) or E[A]=gamma0+gammaU*U+gammaC*C).
#' @param gammaC specifies the relationship between the measured confounder C and the exposure A (i.e. logit(P(A=1)) or E[A]=gamma0+gammaU*U+gammaC*C).
#' @param gammaU specifies the relationship between the measured confounder U and the exposure A (i.e. logit(P(A=1)) or E[A]=gamma0+gammaU*U+gammaC*C).
#' @param varA is the variance of the exposure A when Atype="C". Default is varA=1.
#' @param alpha0 specifies the intercept for the mediator M (i.e. Logit(P(M=1)) or E[M]=alpha0+alphaA*A+alphaC*C+alphaU*U).
#' @param alphaA specifies the relationship between the exposure A and the mediator M (i.e  logit(P(M=1)) or E[M]=alpha0+alphaA*A+alphaC*C+alphaU*U).
#' @param alphaC specifies the relationship between the measured confounder C and the mediator M (i.e. E[M] or logit(P(M=1))=alpha0+alphaA*A+alphaC*C+alphaU*U).
#' @param alphaU specifies the relationship between the unmeasured confounder U and the mediator M (i.e. E[M] or logit(P(M=1))=alpha0+alphaA*A+alphaC*C+alphaU*U).
#' @param varM is the variance of the mediator M when Mtype="C". Default is varM=1.
#' @param beta0 specifies the intercept for the outcome Y (i.e. logit(P(Y=1)) or E[Y]=beta0+betaA*A+betaM*M+betaI*A*M+betaC*C+betaU*U).
#' @param betaA specifies the relationship between the exposure A and the outcome Y (i.e. E[Y] or logit(P(Y=1))=beta0+betaA*A+betaM*M+betaI*A*M+betaC*C+betaU*U).
#' @param betaM specifies the relationship between the mediator M and the outcome Y (i.e. E[Y] or logit(P(Y=1))=beta0+betaA*A+betaM*M+betaI*A*M+betaC*C+betaU*U).
#' @param betaI specifies the interaction between the mediator M and the exposure A on the outcome Y (i.e. E[Y] or logit(P(Y=1))= beta0+betaA*A+betaM*M+betaI*A*M+betaC*C+betaU*U) when the flag interact=TRUE.
#' @param betaC specifies the relationship between the measured confounder C and the outcome Y (i.e. E[Y] or logit(P(Y=1))=beta0+betaA*A+betaM*M+betaI*A*M+betaC*C+betaU*U).
#' @param betaU specifies the relationship between the unmeasured confounder U and the outcome Y (i.e. E[Y] or logit(P(Y=1))=beta0+betaA*A+betaM*M+betaI*A*M+betaC*C+betaU*U).
#' @param varY is the variance of the outcome Y when Ytype="C". Default is varY=1.
#' @param alpha is the significance level. Default value is alpha=0.05.
#' @param nSim is the number of simulations run for the function. The more simulations run, the more accurate the results, but this will make the function slower.
#' @param nBoot is the number of Monte Carlo draws for nonparametric bootstrap or quasi-Bayesian approximation for the mediate function.
#' @param seed sets the seed used for the random generator.
#' @param atreat sets the treatment group for the exposure A.
#' @param acontrol sets the control group for the exposure A.
#' @param use_cpp use RcppEigen (will enable threading if multi-processing is not activated)
#' @param use_multi_processing use multiple processes
#' @param num_jobs number of tasks or cores to use
#' @return The function outputs (1) the proportion of simulations where the average causal mediation effect (ACME) is significant when the model does NOT include U, (2) the proportion of simulations where the ACME is significant when the model includes U, and (3) the proportion of simulations where conclusions based on the ACME match (i.e. the ACME is significant when U is excluded from the model and included in the model or the ACME is not significant when U is excluded from the model and included in the model). The function also outputs (1) the average estimate of the average ACME when U is NOT included in the model, (2) the average ACME when U is included in the model, and (3) the average absolute difference for the ACME when U is included in the model and the ACME when U is excluded from the model. This is given for both the ACME and the average direct effect (ADE). The correlation between variables is also given to show how the change in betas, alphas, and gammas effect the relationship between these variables. Note: correlation is valid if at least on of the variables is normally distributed.
#' @examples
#' testM<- UmediationThread(n=1000,Atype="D",Mtype="C",Ytype="C",Ctype=c("C","D","D"),Utype=c("C","D","D","C"),
#' interact=TRUE,muC=c(0.1,0.3,0.2),varC=c(1,1,1),muU=c(.1,0.3,0.2,.1),varU=c(1,1,1,1),gamma0=0,
#' gammaC=c(1,0.3,0.2),gammaU=c(1,0.3,0.2,0.4),varA=1,alpha0=0,alphaA=1,alphaC=c(0.3,0.2,0.2),
#' alphaU=c(0.3,0.2,0.3,0.2),varM=1,beta0=0,betaA=-1,betaM=1,betaI=1,betaC=c(0.3,0.2,0.1),
#' betaU=c(0.3,0.2,-1.3,0.2),varY=1,alpha=0.05,nSim=100,nBoot=400)
#' 
#' testM
#' 
#' testM_cpp<- UmediationThread(n=1000,Atype="D",Mtype="C",Ytype="C",Ctype=c("C","D","D"),Utype=c("C","D","D","C"),
#' interact=TRUE,muC=c(0.1,0.3,0.2),varC=c(1,1,1),muU=c(.1,0.3,0.2,.1),varU=c(1,1,1,1),gamma0=0,
#' gammaC=c(1,0.3,0.2),gammaU=c(1,0.3,0.2,0.4),varA=1,alpha0=0,alphaA=1,alphaC=c(0.3,0.2,0.2),
#' alphaU=c(0.3,0.2,0.3,0.2),varM=1,beta0=0,betaA=-1,betaM=1,betaI=1,betaC=c(0.3,0.2,0.1),
#' betaU=c(0.3,0.2,-1.3,0.2),varY=1,alpha=0.05,nSim=100,nBoot=400 use_cpp=T, num_jobs=7)
#' 
#' testM_cpp
#' @keywords function mediation unmeasured confounding
#' @section Warning: 
#' library(mediation), library(car), and library(pbapply) are needed to run this function. 
UmediationThread <- function(
  n=100,Atype="D",Mtype="C",Ytype="C",Ctype="C",Utype="C",
  interact=FALSE,muC=0,varC=1,muU=0,varU=1,gamma0=0,gammaC=0,gammaU=0,varA=1,
  alpha0=0,alphaA=0,alphaC=0,alphaU=0,varM=1,
  beta0=0,betaA=0,betaM=0,betaI=0,betaC=0,betaU=0,varY=1,
  alpha=0.05,nSim=250,nBoot=400,seed=1,atreat=1,acontrol=0, 
  use_cpp=F, use_multi_processing=F, num_jobs=1){
  
  if(num_jobs < 1){
    stop("invalid parameter: num_jobs < 1")
  }
  if(num_jobs > 1 && !(use_cpp || use_multi_processing)){
    stop("incompatible parameters: use_cpp and use_multi_processing, are both false, and num_jobs != 1")
  }
  
  #######################################
  # Check Input for Errors
  #######################################
  ErrorCheck(n,Atype,Mtype,Ytype,Ctype,Utype,interact,muC,varC,muU,varU,gamma0,gammaC,gammaU,varA,alpha0,alphaA,alphaC,alphaU,varM,beta0,betaA,betaM,betaI,betaC,betaU,varY,alpha,nSim,nBoot,seed,atreat,acontrol)
  
  #######################################
  # Create Results Matrix
  #######################################
  
  Results<-matrix(0,nrow=12,ncol=1)
  rownames(Results)<-c("Prop. of simulations w/ significant ACME excluding U","Prop. of simulations w/ significant ACME including U","Prop. of simulations where conclusions based on ACME match","Average ACME excluding U","Average ACME including U","Average absolute difference of ACME including U minus ACME excluding U","Prop. of simulations w/ significant ADE excluding U","Prop. of simulations w/ significant ADE including U","Prop. of simulations where conclusions based on ADE match","Average ADE excluding U","Average ADE including U","Average absolute difference of ADE including U minus ADE excluding U")
  
  #######################################
  # Create Input Data Matrix
  #######################################
  
  data_matrix_gen = generate_data_matrix(
    n=n, interact=interact,
    Atype=Atype,Mtype=Mtype,Ytype=Ytype,Ctype=Ctype,Utype=Utype,
    muC=muC,muU=muU,VarA=VarA,varM=varM,varY=varY,varC=varC,varU=varU,
    alpha0=alpha0,alphaA=alphaA,alphaC=alphaC,alphaU=alphaU,
    beta0=beta0,betaA=betaA,betaM=betaM,betaI=betaI,betaC=betaC,betaU=betaU,
    gamma0=gamma0,gammaC=gammaC,gammaU=gammaU,
    nSim=nSim,seed=seed,nBoot=nBoot,
    use_cpp=use_cpp, use_multi_processing=use_multi_processing, num_jobs=num_jobs
  )
  data_matrix = data_matrix_gen$data
  matC = data_matrix_gen$matC
  #######################################
  # Run the mediation and collect the results
  #######################################
  
  if(use_multi_processing){
    options(mediate.jobs = num_jobs)
    if(parallel::detectCores() == 1){
      warning("your machine may not be suitable for multiprocessing, only 1 core was detected")
    }
    if(num_jobs < 2){
      stop("There is no point in using MultiProcessing with less than 2 jobs")
    }
    if((nSim / num_jobs) < 1.0){
      warning(paste("you don't have enough Simulations in nSim:", nSim, " to fully benefit from num_jobs:", num_jobs, sep=""))
    }
    result.matrix = mediate_parallel(data_matrix)
  } else {
    result.matrix = pbapply::pblapply(data_matrix, perform_mediation)
    dim(result.matrix) = dim(data_matrix)
  }
  rm(data_matrix)
  
  #######################################
  # Loop through the mediation result matrix
  #######################################
  
  for(si in 1:nSim){ #loop through all simulations
    
    # TODO access the med.out/ med.outU values for this simulation
    result_element = result.matrix[[si]]
    
    med.out = result_element[["med.out"]]
    med.outU = result_element[["med.outU"]]
    
    med.summary = mediation::summary.mediate(med.out)
    medU.summary = mediation::summary.mediate(med.outU)
    
    med.d.avg = med.summary[["d.avg"]]
    med.z.avg = med.summary[["z.avg"]]
    med.d.avg.p = med.summary[["d.avg.p"]]
    med.z.avg.p = med.summary[["z.avg.p"]]
    
    medU.d.avg = medU.summary[["d.avg"]]
    medU.z.avg = medU.summary[["z.avg"]]
    medU.d.avg.p = medU.summary[["d.avg.p"]]
    medU.z.avg.p = medU.summary[["z.avg.p"]]
    
    #Results
    Results["Average ACME excluding U",1]<-Results["Average ACME excluding U",1]+med.d.avg
    Results["Average ACME including U",1]<-Results["Average ACME including U",1]+medU.d.avg
    Results["Average absolute difference of ACME including U minus ACME excluding U",1]<-Results["Average absolute difference of ACME including U minus ACME excluding U",1]+abs(med.d.avg-medU.d.avg)
    
    Results["Average ADE excluding U",1]<-Results["Average ADE excluding U",1]+med.z.avg
    Results["Average ADE including U",1]<-Results["Average ADE including U",1]+medU.z.avg
    Results["Average absolute difference of ADE including U minus ADE excluding U",1]<-Results["Average absolute difference of ADE including U minus ADE excluding U",1]+abs(med.z.avg-medU.z.avg)
    
    if(med.d.avg.p<alpha){Results["Prop. of simulations w/ significant ACME excluding U",1]<- Results["Prop. of simulations w/ significant ACME excluding U",1]+1}
    if(medU.d.avg.p<alpha){Results["Prop. of simulations w/ significant ACME including U",1]<- Results["Prop. of simulations w/ significant ACME including U",1]+1}
    if(((med.d.avg.p<alpha)&(medU.d.avg.p<alpha))|((med.d.avg.p>alpha)&(medU.d.avg.p>alpha))){Results["Prop. of simulations where conclusions based on ACME match",1]<- Results["Prop. of simulations where conclusions based on ACME match",1]+1}
    
    if(med.z.avg.p<alpha){Results["Prop. of simulations w/ significant ADE excluding U",1]<- Results["Prop. of simulations w/ significant ADE excluding U",1]+1}
    if(medU.z.avg.p<alpha){Results["Prop. of simulations w/ significant ADE including U",1]<- Results["Prop. of simulations w/ significant ADE including U",1]+1}
    if(((med.z.avg.p<alpha)&(medU.z.avg.p<alpha))|((med.z.avg.p>alpha)&(medU.z.avg.p>alpha))){Results["Prop. of simulations where conclusions based on ADE match",1]<- Results["Prop. of simulations where conclusions based on ADE match",1]+1}
  }
  
  #######################################
  # Results
  #######################################
  
  matR<-Results/nSim
  listA<-list(matR,matC,paste("Warning: correlations are only valid if at least one of the variables is normally distributed."))
  names(listA)<-c("Results","Correlations_Between_Variables","Warning")
  listA
  
} #end of function
