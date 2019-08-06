#' @import car
#' @import mediation

#' @export
#' @title Umediation
#' @name Umediation
#' @aliases Umediation
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
#' @return The function outputs (1) the proportion of simulations where the average causal mediation effect (ACME) is significant when the model does NOT include U, (2) the proportion of simulations where the ACME is significant when the model includes U, and (3) the proportion of simulations where conclusions based on the ACME match (i.e. the ACME is significant when U is excluded from the model and included in the model or the ACME is not significant when U is excluded from the model and included in the model). The function also outputs (1) the average estimate of the average ACME when U is NOT included in the model, (2) the average ACME when U is included in the model, and (3) the average absolute difference for the ACME when U is included in the model and the ACME when U is excluded from the model. This is given for both the ACME and the average direct effect (ADE). The correlation between variables is also given to show how the change in betas, alphas, and gammas effect the relationship between these variables. Note: correlation is valid if at least on of the variables is normally distributed.
#' @examples
#' testM<- Umediation(n=1000,Atype="D",Mtype="C",Ytype="C",Ctype=c("C","D","D"),Utype=c("C","D","D","C"),
#' interact=TRUE,muC=c(0.1,0.3,0.2),varC=c(1,1,1),muU=c(.1,0.3,0.2,.1),varU=c(1,1,1,1),gamma0=0,
#' gammaC=c(1,0.3,0.2),gammaU=c(1,0.3,0.2,0.4),varA=1,alpha0=0,alphaA=1,alphaC=c(0.3,0.2,0.2),
#' alphaU=c(0.3,0.2,0.3,0.2),varM=1,beta0=0,betaA=-1,betaM=1,betaI=1,betaC=c(0.3,0.2,0.1),
#' betaU=c(0.3,0.2,-1.3,0.2),varY=1,alpha=0.05,nSim=100,nBoot=400)
#' 
#' testM
#' 
#' testM_cpp<- Umediation(n=1000,Atype="D",Mtype="C",Ytype="C",Ctype=c("C","D","D"),Utype=c("C","D","D","C"),
#' interact=TRUE,muC=c(0.1,0.3,0.2),varC=c(1,1,1),muU=c(.1,0.3,0.2,.1),varU=c(1,1,1,1),gamma0=0,
#' gammaC=c(1,0.3,0.2),gammaU=c(1,0.3,0.2,0.4),varA=1,alpha0=0,alphaA=1,alphaC=c(0.3,0.2,0.2),
#' alphaU=c(0.3,0.2,0.3,0.2),varM=1,beta0=0,betaA=-1,betaM=1,betaI=1,betaC=c(0.3,0.2,0.1),
#' betaU=c(0.3,0.2,-1.3,0.2),varY=1,alpha=0.05,nSim=100,nBoot=400 use_cpp=T, num_cores=7)
#' 
#' testM_cpp
#' @keywords function mediation unmeasured confounding
#' @section Warning: 
#' library(mediation) and library(car) are needed to run this function. 
Umediation <- function(
  n=100,Atype="D",Mtype="C",Ytype="C",Ctype="C",Utype="C",interact=FALSE,muC=0,varC=1,muU=0,varU=1,gamma0=0,gammaC=0,gammaU=0,
  varA=1,alpha0=0,alphaA=0,alphaC=0,alphaU=0,varM=1,beta0=0,betaA=0,betaM=0,betaI=0,betaC=0,betaU=0,varY=1,alpha=0.05,nSim=250,
  nBoot=400,seed=1,atreat=1,acontrol=0, use_cpp=F, num_cores=getOption("mediate.threads", default = 1)){
  if(num_cores > 1 && !use_cpp){
    stop("incompatible parameters: use_cpp is false, and num_cores > 1")
  }
  options(mediate.threads = num_cores)
  
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
  # Loop through all simulations
  #######################################
  
  for(si in 1:nSim){ #loop through all simulations
    set.seed(seed+(si-1)) #set the seed to ensure reproducibility
    
    simStep<-1
    if(floor(si/simStep)==ceiling(si/simStep)){print(paste("simulation",si, "of",nSim))}
    
    #######################################
    # Generate the unmeasured confounder,
    # measured confounder C, exposure A,
    # mediator M, and outcome Y
    #######################################
    
    # generate unmeasured confounder U
    U<-matrix(0,nrow=n,ncol=length(Utype))
    for(iu in 1:length(Utype)){
      if(Utype[iu]=="C"){
        if(varU[iu]==0|varU[iu]<0){stop(paste("Error: Variance of U",iu," (varU) is 0 or negative.",sep=""))}
        U[,iu]<-rnorm(n,muU[iu],sqrt(varU[iu]))
      }
      if(Utype[iu]=="D"){
        if(muU[iu]<0|muU[iu]==0|muU[iu]==1|muU[iu]>1){stop(paste("Error: Prob(U",iu,"=1) (i.e. muU",iu,") is less than or equal to zero or greater than or equal to 1.",sep=""))}
        U[,iu]<-rbinom(n,1,muU[iu])
      }
      colnames(U)<-paste("U",length(Utype),Utype,sep="")
    }
    
    # generate measured confounder C
    CC<-matrix(0,nrow=n,ncol=length(Ctype))
    for(ic in 1:length(Ctype)){
      if(Ctype[ic]=="C"){
        if(varC[ic]==0|varC[ic]<0){stop(paste("Error: Variance of C",ic," (varC) is 0 or negative.",sep=""))}
        CC[,ic]<-rnorm(n,muC[ic],sqrt(varC[ic]))
      }
      if(Ctype[ic]=="D"){
        if(muC[ic]<0|muC[ic]==0|muC[ic]==1|muC[ic]>1){stop(paste("Error: Prob(C",ic,"=1) is less than or equal to zero or greater than or equal to 1.",sep=""))}
        CC[,ic]<-rbinom(n,1,muC[ic])
      }
      colnames(CC)<-paste("C",length(Ctype),Ctype,sep="")
    }
    
    # generate exposure A as a function of C and U
    # logit(P(A=1)) or E[A]=gamma0+gammaU*U+gammaC*C
    muA<-gamma0+CC%*%gammaC+U%*%gammaU
    if(Atype=="C"){
      A<-rnorm(n,muA,sqrt(varA)) # normally distributed mediator
      if(varA==0|varA<0){stop("Error: Variance of M (varM) is 0 or negative.")}
    }
    if(Atype=="D"){
      Ap<-exp(muA)/(1+exp(muA))
      if(max(Ap)=="NaN"|is.na(max(Ap))){stop("Error: Prob(A=1) is too close to 0 or 1. This can occur if the absolute value of gamma0 is too large.")}
      A<-rbinom(n,1,Ap)# binary mediator
      if(length(A[A==1])/length(A)<(2/n)|length(A[A==1])/length(A)>(1-2/n)){
        stop("Error: there was no enough varaibality in A (i.e. Prob(A=1)=1 or Prob(A=1)=0). Considering increasing the sample size n or using mean centered confounders (i.e. muC=0 and muU=0).")
      }
    }
    
    # generate mediator M
    # Logit(P(M=1)) or E[M]=alpha0+alphaA*A+alphaC*C+alphaU*U)
    muM<-alpha0+alphaA*A+CC%*%alphaC+U%*%alphaU
    if(Mtype=="C"){
      M<-rnorm(n,muM,sqrt(varM)) # normally distributed mediator
      if(varM==0|varM<0){stop("Error: Variance of M (varM) is 0 or negative.")}
    }
    if(Mtype=="D"){
      Mp<-exp(muM)/(1+exp(muM))
      if(max(Mp)=="NaN"|is.na(max(Mp))){stop("Error: Prob(M=1) is too close to 0 or 1. This can occur if the absolute value of alpha0 is too large.")}
      M<-rbinom(n,1,Mp)# binary mediator
      if(length(M[M==1])/length(M)<(5/n)|length(M[M==1])/length(M)>(1-5/n)){
        stop("Error: there was no enough varaibality in M (i.e. Prob(M=1)=1 or Prob(M=1)=0). Considering increasing the sample size n or using mean centered confounders (i.e. muC=0 and muU=0).")
      }
    }
    
    
    # generate outcome Y
    # logit(P(Y=1)) or E[Y]=beta0+betaA*A+betaM*M+betaC*C+betaU*U
    if(interact==TRUE){ muY<-beta0+betaA*A+betaM*M+betaI*A*M+CC%*%betaC+U%*%betaU}
    if(interact==FALSE){ muY<-beta0+betaA*A+betaM*M+CC%*%betaC+U%*%betaU}
    if(Ytype=="C"){
      if(varY==0|varY<0){stop("Error: Variance of Y (varY) is 0 or negative.")}
      Y<-rnorm(n,muY,sqrt(varY)) # normally distributed outcome
    }
    if(Ytype=="D"){
      Yp<-exp(muY)/(1+exp(muY))
      if(max(Yp)=="NaN"|is.na(max(Yp))){stop("Error: Prob(Y=1) is too close to 0 or 1. This can occur if the absolute value of beta0 is too large.")}
      Y<-rbinom(n,1,Yp) # binary outcome
      if(length(Y[Y==1])/length(Y)<(5/n)|length(Y[Y==1])/length(Y)>(1-5/n)){
        stop("Error: there was no enough varaibality in Y (i.e. Prob(Y=1)=1 or Prob(Y=1)=0). Considering increasing the sample size n or using mean centered confounders (i.e. muC=0 and muU=0).")
      }
    }
    
    #######################################
    # Check for collinearity
    #######################################
    
    if(length(Utype)>1|length(Ctype)>1){
      if(max(vif(lm(Y~A+M+CC+U))[,3])>10 & interact==FALSE ){stop("Error: There is strong evidence of collinearity (i.e. VIF>10) when A,M,C,U are regressed on Y. Consider reducing the absolute value of large gammas, alphas, or betas.")}
      
      if(max(vif(lm(Y~A+M+A*M+CC+U))[,3])>10 & interact==TRUE ){stop("Error: There is strong evidence of collinearity (i.e. VIF>10) when A,M,A*M,C,U are regressed on Y. Consider reducing the absolute value of large gammas, alphas, or betas.")}
      
      if(max(vif(lm(M~A+CC+U))[,3])>10){stop("Error: There is strong evidence of collinearity (i.e. VIF>10) when A,C,U are regressed on M. Consider reducing the absolute value of large gammas or alphas.")}
    }
    
    if(length(Utype)==1&length(Ctype)==1){
      if(max(vif(lm(Y~A+M+CC+U)))>10 & interact==FALSE ){stop("Error: There is strong evidence of collinearity (i.e. VIF>10) when A,M,C,U are regressed on Y. Consider reducing the absolute value of large gammas, alphas, or betas.")}
      
      if(max(vif(lm(Y~A+M+A*M+CC+U)))>10 & interact==TRUE ){stop("Error: There is strong evidence of collinearity (i.e. VIF>10) when A,M,A*M,C,U are regressed on Y. Consider reducing the absolute value of large gammas, alphas, or betas.")}
      
      if(max(vif(lm(M~A+CC+U)))>10){stop("Error: There is strong evidence of collinearity (i.e. VIF>10) when A,C,U are regressed on M. Consider reducing the absolute value of large gammas or alphas.")}
    }
    
    #######################################
    # Correlation for U,C,A,M,Y
    #######################################
    
    matA<-cbind(A,M,Y,CC,U)
    matC<-round(cor(matA),digits=2)
    colnames(matC)<-c("A","M","Y",paste("C",c(1:length(Ctype)),sep=""),paste("U",c(1:length(Utype)),sep=""))
    rownames(matC)<-c("A","M","Y",paste("C",c(1:length(Ctype)),sep=""),paste("U",c(1:length(Utype)),sep=""))
    
    #######################################
    # Mediation analysis w/ and w/out U
    #######################################
    
    # models without U
    if(Mtype=="C"){med.fit <-lm(M~A+CC)} #fit the model for the mediator
    if(Mtype=="D"){med.fit <-glm(M~A+CC,family=binomial(link = "logit"))}
    if(interact==FALSE){
      if(Ytype=="C"){out.fit<-lm(Y~M+A+CC)} #fit the model for the outcome
      if(Ytype=="D"){out.fit<-glm(Y~M+A+CC,family=binomial(link = "logit"))}
    }
    if(interact==TRUE){
      if(Ytype=="C"){out.fit<-lm(Y~M+A+CC+A*M)} #fit the model for the outcome
      if(Ytype=="D"){out.fit<-glm(Y~M+A+CC+A*M,family=binomial(link = "logit"))}
    }
    g_env = globalenv()
    g_env[["med.fit"]] = med.fit
    g_env[["out.fit"]] = out.fit
    if(use_cpp){
      med.out <- mediate_with_rcpp(med.fit, out.fit, treat = "A", mediator = "M",sim=nBoot)
    } else {
      med.out <- mediate(med.fit, out.fit, treat = "A", mediator = "M",sim=nBoot)
    }
    
    
    # models with U
    if(Mtype=="C"){med.fitU <-lm(M~A+CC+U)} #fit the model for the mediator
    if(Mtype=="D"){med.fitU <-glm(M~A+CC+U,family=binomial(link = "logit"))}
    if(interact==FALSE){
      if(Ytype=="C"){out.fitU<-lm(Y~M+A+CC+U)} #fit the model for the outcome
      if(Ytype=="D"){out.fitU<-glm(Y~M+A+CC+U,family=binomial(link = "logit"))}
    }
    if(interact==TRUE){
      if(Ytype=="C"){out.fitU<-lm(Y~M+A+CC+U+A*M)} #fit the model for the outcome
      if(Ytype=="D"){out.fitU<-glm(Y~M+A+CC+U+A*M,family=binomial(link = "logit"))}
    }
    
    g_env[["med.fitU"]] = med.fitU
    g_env[["out.fitU"]] = out.fitU
    
    if(use_cpp){
      med.outU <- mediate_with_rcpp(med.fitU, out.fitU, treat = "A", mediator = "M",sim=nBoot)
    } else {
      med.outU <- mediate(med.fitU, out.fitU, treat = "A", mediator = "M",sim=nBoot)
    }
    
    
    #Results
    Results["Average ACME excluding U",1]<-Results["Average ACME excluding U",1]+summary(med.out)$d.avg
    Results["Average ACME including U",1]<-Results["Average ACME including U",1]+summary(med.outU)$d.avg
    Results["Average absolute difference of ACME including U minus ACME excluding U",1]<-Results["Average absolute difference of ACME including U minus ACME excluding U",1]+abs(summary(med.out)$d.avg-summary(med.outU)$d.avg)
    
    Results["Average ADE excluding U",1]<-Results["Average ADE excluding U",1]+summary(med.out)$z.avg
    Results["Average ADE including U",1]<-Results["Average ADE including U",1]+summary(med.outU)$z.avg
    Results["Average absolute difference of ADE including U minus ADE excluding U",1]<-Results["Average absolute difference of ADE including U minus ADE excluding U",1]+abs(summary(med.out)$z.avg-summary(med.outU)$z.avg)
    
    if(summary(med.out)$d.avg.p<alpha){Results["Prop. of simulations w/ significant ACME excluding U",1]<- Results["Prop. of simulations w/ significant ACME excluding U",1]+1}
    if(summary(med.outU)$d.avg.p<alpha){Results["Prop. of simulations w/ significant ACME including U",1]<- Results["Prop. of simulations w/ significant ACME including U",1]+1}
    if(((summary(med.out)$d.avg.p<alpha)&(summary(med.outU)$d.avg.p<alpha))|((summary(med.out)$d.avg.p>alpha)&(summary(med.outU)$d.avg.p>alpha))){Results["Prop. of simulations where conclusions based on ACME match",1]<- Results["Prop. of simulations where conclusions based on ACME match",1]+1}
    
    if(summary(med.out)$z.avg.p<alpha){Results["Prop. of simulations w/ significant ADE excluding U",1]<- Results["Prop. of simulations w/ significant ADE excluding U",1]+1}
    if(summary(med.outU)$z.avg.p<alpha){Results["Prop. of simulations w/ significant ADE including U",1]<- Results["Prop. of simulations w/ significant ADE including U",1]+1}
    if(((summary(med.out)$z.avg.p<alpha)&(summary(med.outU)$z.avg.p<alpha))|((summary(med.out)$z.avg.p>alpha)&(summary(med.outU)$z.avg.p>alpha))){Results["Prop. of simulations where conclusions based on ADE match",1]<- Results["Prop. of simulations where conclusions based on ADE match",1]+1}
    
    
  } #end of simulation loop
  
  #######################################
  # Results
  #######################################
  
  matR<-Results/nSim
  listA<-list(matR,matC,paste("Warning: correlations are only valid if at least one of the variables is normally distributed."))
  names(listA)<-c("Results","Correlations_Between_Variables","Warning")
  listA
  
} #end of function
