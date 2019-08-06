

#' @export
#' @title ErrorCheck
#' @name ErrorCheck
#' @aliases ErrorCheck
#' @description Performs error checks on the input for the Umediation function.
#' @author Sharon Lutz
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
#' @return The function exits with an error message if the error checks are not met.
ErrorCheck<-function(n=1000,Atype="D",Mtype="C",Ytype="C",Ctype="C",Utype="C",
                     interact=FALSE,muC=0,varC=1,muU=0,varU=1,gamma0=0,gammaC=0,
                     gammaU=0,varA=1,alpha0=0,alphaA=0,alphaC=0,alphaU=0,varM=1,
                     beta0=0,betaA=0,betaM=0,betaI=0,betaC=0,betaU=0,varY=1,alpha=0.05,
                     nSim=300,nBoot=500,seed=1,atreat=1,acontrol=0){
    
    ##################
    # Models of exposure A,
    # mediator M, outcome Y
    ##################
    # logit(P(A=1)) or E[A]=gamma0+gammaC*C+gammaU*U
    # Logit(P(M=1)) or E[M]=alpha0+alphaA*A+alphaC*C+alphaU*U)
    # logit(P(Y=1)) or E[Y]=beta0+betaA*A+betaM*M+betaI*A*M+betaC*C+betaU*U

    ##################
    # check A,M,Y,C,U type vectors are only C or D
    ##################
    if(Atype!="C"& Atype!="D"){stop("Error: Atype must be D for dichtomous or C for continuous for exposure A.")}
    if(Mtype!="C"& Mtype!="D"){stop("Error: Mtype must be D for dichtomous or C for continuous for mediator M.")}
    if(Ytype!="C"& Ytype!="D"){stop("Error: Ytype must be D for dichtomous or C for continuous for outcome Y.")}
    for(jj in 1:length(Ctype)){if(Ctype[jj]!="C" & Ctype[jj]!="D"){stop("Error: Ctype must be D for dichtomous or C for continuous measured confounders C.")}}
    for(jj in 1:length(Utype)){if(Utype[jj]!="C" & Utype[jj]!="D"){stop("Error: Utype must be D for dichtomous or C for continuous unmeasured confounders U.")}}
 
     ##################
     # check A variables
     ##################
     if((length(Atype)!=length(varA))| (length(Atype)!=length(gamma0))| (length(Atype)!=length(alphaA))|(length(Atype)!=length(betaA)) |(length(Atype)!=1 ) ){stop("Error: The length of Atype does not equal the length of varA, gamma0, alphaA, and or betaA, which does not equal one. This function does not accomodate multiple exposures A.")}
     
     ##################
     # check M variables
     ##################
     if( (length(Mtype)!=length(varM))|(length(Mtype)!=length(alpha0)) |(length(Mtype)!=length(betaM)) |(length(Mtype)!=1 ) ){stop("Error: The length of Mtype does not equal the length of varM, alpha0, and or betaM, which does not equal one. This function does not accomodate multiple mediators M.")}
    
    ##################
     # check Y variables
     ##################
     if((length(Ytype)!=length(varY))|  (length(Ytype)!=length(beta0)) |(length(Ytype)!=1)){stop("Error: The length of Ytype does not equal the length of varY, beta0, and or betaI, which does not equal one. This function does not accomodate multiple outcomes Y.")}
   if(interact==TRUE & (length(betaI)!=1)){stop("Error: length of betaI must be 1 when interact flag is TRUE.")}
   
     ##################
     # check C variables
     ##################
     if((length(Ctype)!=length(muC)) |(length(Ctype)!=length(varC)) |(length(Ctype)!=length(gammaC)) |(length(Ctype)!=length(alphaC))|(length(Ctype)!=length(betaC))){stop(paste("Error: Vectors for Ctype, muC, varC, gammaC, alphaC, and betaC must all be the same length."))}
     for(j1 in 1:length(Ctype)){
         if(Ctype[j1]=="D"){
             if(muC[j1]<0|muC[j1]==0|muC[j1]==1|muC[j1]>1){stop(paste("Error: Prob(C",j1,"=1) is ",muC[j1],". Prob(C",j1,"=1) must be between 0 and 1 for Ctype",j1,". ",sep=""))}
         }}
         
     ##################
     # check U variables
     ##################
     if((length(Utype)!=length(muU)) |(length(Utype)!=length(varU)) |(length(Utype)!=length(gammaU))|(length(Utype)!=length(alphaU))|(length(Utype)!=length(betaU))){stop(paste("Error: Vectors for Utype, muU, varU, gammaU, alphaU, and betaU must all be the same length."))}
     for(j2 in 1:length(Utype)){
         if(Utype[j2]=="D"){
             if(muU[j2]<0|muU[j2]==0|muU[j2]==1|muU[j2]>1){stop(paste("Error: Prob(U",j2,"=1) is ",muU[j2],". Prob(U",j2,"=1) must be between 0 and 1 for Utype",j2,". ",sep=""))}
         }}
     
     ##################
     #check interact flag
     ##################
     if(interact!=FALSE & interact!=TRUE){stop("Error: interact must be TRUE or FALSE. This flag is case sensitive.")}
     if(length(interact)!=1){stop("Error: length of the vector for interact must be 1.")}
 
     ##################
     #check alpha
     ##################
     if(alpha<0 | alpha>1|alpha==0|alpha==1){stop("Error: alpha must be greater than 0 and less than 1")}
     if(length(alpha)!=1){stop("Error: length of the vector for alpha must be 1")}
     
     ##################
     #check n, nSim, nBoot, seed
     ##################
     if(length(n)!=1){stop("Error: length of the vector for n must be 1")}
     if(length(nSim)!=1){stop("Error: length of the vector for nSim must be 1")}
     if(length(nBoot)!=1){stop("Error: length of the vector for nBoot must be 1")}
     if(length(seed)!=1){stop("Error: length of the vector for seed must be 1")}
     if((floor(n)!=ceiling(n))|(n<0)|(n==0)){stop(paste("Error: the sample size n must be an integer greater than 0"))}
     if((floor(nSim)!=ceiling(nSim))|(nSim<0)|(nSim==0)){stop(paste("Error: the number of simulations nSim must be an integer greater than 0"))}
     if((floor(nBoot)!=ceiling(nBoot))|(nBoot<0)|(nBoot==0)){stop(paste("Error: the number of bootstrap samples nBoot must be an integer greater than 0"))}
      if((floor(seed)!=ceiling(seed))|(seed<0)|(seed==0)){stop(paste("Error: the seed must be an integer greater than 0"))}
      
      ##################
      #atreat and acontrol
      ##################
      if(length(atreat)!=1 |length(acontrol)!=1 ){stop("Error: length of the vector for atreat and acontrol must be 1")}
      
     ###################
}#end of function
