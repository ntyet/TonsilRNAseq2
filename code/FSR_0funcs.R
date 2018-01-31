library(ggplot2) # gor ggplot
# library(ggplot2, lib.loc = "/home/ntyet/R/x86_64-pc-linux-gnu-library/3.2")# For condo server
library(sva)
# install.packages("ruv")
library(ruv)
# library(dSVA)
# library(svapls)
# library(dSVA, lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.3")# For condo2017 server
library(plyr) # for llply
library(reshape2) # for  melt
library(limma) # for is.fullrank
library(parallel) # for mclapply
library(MASS) # for mvrnorm function
source("../code/multtest.txt")
library(AUC)
library(data.table)
# counts <- counts[1:100, ]
set.seed(20170204)
## Forward Selection With Voom--------
# Input: FixCov, VarCov: Data frame with colnames
#        FixCov: is the fixed covariates that are not subjected to variable selection
#        VarCov: is the flexible covariates that are subjected to variable selection
#        Count data: Data frame with colnames
#        print.progress: logical factor T or F to indicate if we want to see what happens at each iteration of forward selection
fs <- function(counts, FixCov, VarCov, print.progress = F){
  if(!is.data.frame(FixCov)|!is.data.frame(VarCov)) stop("FixCov and VarCov have to be data frame with different column names")
  SelCov <- FixCov # Initiaize Set of Selected Covariates
  ConCov <- VarCov # Initialize Set of Considered Covariates
  PvList <- list() # List of P-values for the covariate selected at each iteration
  BestP5 <- list() # Number of p-value less than 0.05 for the covariate selected at each iteration
  if(ncol(ConCov)!=0){
  Iter <- 1
  repeat{
    if(print.progress){
      cat("----------------------------------------\n")
      cat("Iteration = ", Iter, "\n")
      cat("SelCov = ", names(SelCov), "\n")
    }
    BestP5[[Iter]] <- 0
    Ind <- 0       # Count the number of time that we cannot add new considered covariate to forward selection because of singularity of design matrix
    for(i in 1:ncol(ConCov)){
      # cat("i = ", i, "\n")
      AllCov <- cbind(SelCov, ConCov[i])
      names(AllCov)
      dm <- model.matrix(formula(paste0("~", paste0(names(AllCov), collapse = "+"))), data = AllCov)
      colnames(dm)[1] <- "Intercept"
      if(!is.fullrank(dm)|ncol(dm)==nrow(dm)){Ind <- Ind+1; next}
      
      vout <- voom(counts = counts, design = dm, lib.size = apply(counts, 2, quantile, .75), plot = F)
      suppressWarnings(fit <- lmFit(vout))
      if(is.factor(ConCov[,i]) | is.character(ConCov[,i])) {
        ct <- paste0(grep(paste0(names(ConCov)[i]),  x = colnames(dm), value = T), collapse = ",  ")
      }else{
        ct <- paste0(grep(paste0(names(ConCov)[i], "$"),  x = colnames(dm), value = T), collapse = ",  ") 
        }
      C.matrix <- eval(parse(text=paste0("makeContrasts(",  ct, ",levels = dm)")))
      fit <- contrasts.fit(fit, contrasts =C.matrix)
      fit <- eBayes(fit)
      tt <- topTableF(fit, sort ="none", n = Inf)
      pv <- tt$P.Value
      P5 <- max(sum(pv<=.05, na.rm = T), 1)/max(sum(pv>.75, na.rm = T)/5, 1) # sometime pv has missing value NA make P5 cannot calculable
      if(print.progress)cat("ConCov = ", names(ConCov)[i], ", P5 = ", P5,  "\n")
      if(BestP5[[Iter]] <=P5) {
        BestCov <- i
        BestP5[[Iter]] <- P5
        PvList[[Iter]] <- pv
        names(PvList)[Iter] <- names(BestP5)[Iter] <- names(ConCov)[i] 
        }
    }
    if(Ind == ncol(ConCov)) {
      BestP5[[Iter]] <- NULL
      if(print.progress){
        cat("Current Iteration, Iter = ", Iter,  "Cannot Add more Covariates in Forward Selection\n")
        cat("dim of design matrix now is ", dim(dm), "\n")
      }
      break
    }
    Iter <- Iter +1
    SelCov <- cbind(SelCov, ConCov[BestCov])
    ConCov <- ConCov[-BestCov]
    if(ncol(ConCov)==0)break
  }
  BestP5 <- do.call("c", BestP5)
  }else{BestP5 <- numeric(0)}
  # P <- data.frame(t(do.call("rbind", PvList)))
  BestP5
}

#Example-----------
# pm <- proc.time()
# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order", "neut", "lymp", "mono", "eosi", "baso")]
# out <- fs(counts = counts[1:50,], FixCov = FixCov, VarCov = VarCov, print.progress = T)
#----------------------------------------------



# Forward Selection for Real Data -----------
# This function running forward selection for real data, the output will be the list
# of  2 components whose first component is the covariates entering the model in the sequence together
# with its number of p-values ratio,  and the second component is a vector of number of covariates 
# entering the model when taking IV threshold alpha

# Input: counts: counts dataset
#        FixCov: a dataframe of fix covariates that are not subjected to variable selection
#        VarCov: a dataframe of flexible covariates that are subjected to variable selection
#        alpha:  a numeric vector of threshold values
#        print.progress: a logical variable either T or F to indicate that if we want to print the running process of fs function or not
# Output: is a list of 2 components: 
#        fsr: a vector of variable importance measurement in the order of the covariate entering the model in fs process. 
#             length of fsr equals number of column of VarCov
#        S: a vector of the number of selected variates at different threshold alpha. 
DataFS <- function(counts, FixCov, VarCov, alpha, print.progress = F){
  
  FS <- fs(counts = counts, FixCov = FixCov, VarCov = VarCov, print.progress = print.progress)
  S <- vapply(1:length(alpha), function(i)sum(cummin(FS) >alpha[i]), FUN.VALUE = 1.0)
  return(list(FS = FS, S = S, alpha = alpha))
}

#Example--------
# amax <-  15
# alpha <- c(seq(1, amax, length = 800))
# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# out <- DataFS(counts = counts[1:50,], FixCov = FixCov, VarCov = VarCov, alpha, print.progress = T)
# out$S
#out$fsr


#Forward Selection when Adding phony variables----------------------------------------------
# create projection matrix -----------
# Input: dm0 is the design matrix containing all considered covariates
# This function calculate projection matrix of the design matrix created by data frame AllCov0
ProjMat <- function(dm0){
  Hx <- dm0%*%chol2inv(chol(crossprod(dm0, dm0)))%*%t(dm0)
  Ix <- diag(nrow(dm0))
  Ix - Hx
  }


#Example-------
# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# AllCov0 <- cbind(FixCov, VarCov)
# dm0 <- model.matrix(formula(paste0("~", paste0(names(AllCov0), collapse = "+"))), data = AllCov0)
# colnames(dm0)[1] <- "Intercept"
# Px <- ProjMat(dm0)

### Phony variable generation-----------
# This function generates phony variables in 4 different cases: white noise, permutation of design row, orthogonalize of the design matrixs
# Input: dm0: is the design matrix including all variables considered
#        this is valid only if the number of variables considered is less than the number of samples
#        m: number of phony variables generated
#        option: either WN, OWN, RX, ORX
PhoVar <- function(dm0,  m, option){
  Px <- ProjMat(dm0)
  if(option =="WN"){
    Nois <- t(mvrnorm(n = m, mu = rep(0, nrow(dm0)), Sigma = diag(nrow(dm0))))
    if(m==1) Nois <- data.matrix(mvrnorm(n = m, mu = rep(0, nrow(dm0)), Sigma = diag(nrow(dm0))))
  }
  else if(option == "OWN"){
    Nois0 <- t(mvrnorm(n = m, mu = rep(0, nrow(dm0)), Sigma = diag(nrow(dm0))))
    if(m==1) Nois0 <- data.matrix(mvrnorm(n = m, mu = rep(0, nrow(dm0)), Sigma = diag(nrow(dm0))))
    Nois <- Px%*%Nois0
  } else if (option == "RX"){
    Nois <- data.matrix((data.frame(dm0)[-1])[sample(nrow(dm0)), sample(ncol((data.frame(dm0))[-1]), m)])
  } else {
    Nois0 <- data.matrix((data.frame(dm0)[-1])[sample(nrow(dm0)), sample(ncol((data.frame(dm0))[-1]), m)])
    Nois <- Px%*%Nois0
  }
  colnames(Nois) <- NULL
  Nois <- data.frame(Nois)
  Nois
}

#Example---------------------------------

# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# AllCov0 <- cbind(FixCov, VarCov)
# dm0 <- model.matrix(formula(paste0("~", paste0(names(AllCov0), collapse = "+"))), data = AllCov0)
# colnames(dm0)[1] <- "Intercept"
# PhoVar(dm0 = dm0, m = 2, option = "RX")


### Bootstrap function to generate data and calculate U, SP for the simulatied data-------
## Input: nrep: seed of the bootstrap
##        option: either white noise (WN), permutation of row of X (RX), or orthogonal version (OWN, ORX)
##        FixCov: Fix covariates that are not subjected to variable selection
##        VarCov:  Flexible covariates that are subjected to variable selection
##        m: number of phony variables need to generated
##        counts: counts data with G genes and s columns
##        alpha: a vector of threshold of VI measurement
# cor(covset[c("RFI", "RINa", "RINb", "Conca", "Concb", "lymp", 
#              "neut", "mono", "baso", "eosi")])
BootFS <- function(counts, FixCov, VarCov, m=ncol(VarCov), nrep, alpha, option = "OWN", print.progress = F){
  set.seed(nrep)
   # cat("nrep = ", nrep, "\n")
  AllCov0 <- cbind(FixCov, VarCov)
  dm0 <- model.matrix(formula(paste0("~", paste0(names(AllCov0), collapse = "+"))), data = AllCov0)
  colnames(dm0)[1] <- "Intercept"
  # Nois <- PhoVar(dm0 = dm0, m = m, option = option)
  # VarCov2 <- cbind(VarCov, Nois)
  # 
  # 
  repeat{
    Nois <- PhoVar(dm0 = dm0, m = m, option = option)
    VarCov2 <- cbind(VarCov, Nois)
    AllCov1 <- cbind(FixCov, VarCov2)
    dm1 <- model.matrix(formula(paste0("~", paste0(names(AllCov1), collapse = "+"))), data = AllCov1)
    if(is.fullrank(dm1)) break
  }
  
  pv <- fs(counts = counts, FixCov = FixCov, VarCov = VarCov2, print.progress = print.progress)
  pvmax <- cummin(pv)
  vor <- names(pv)
  U <- vapply(1:length(alpha), function(i)sum(vor%in%names(Nois) & pvmax>alpha[i]), FUN.VALUE = 1.0)
  SP <- vapply(1:length(alpha), function(i)sum( pvmax>alpha[i]), FUN.VALUE = 1.0)
  res <- c(U = U, SP = SP)
   cat("nrep = ", nrep, "\n")
  res
}


#Example---------------------
# amax <-  15
# alpha <- c(seq(1, amax, length = 800))
# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# option <- "ORX"
# out <- BootFS(counts = counts[1:50,], FixCov = FixCov, VarCov = VarCov, m = 20, nrep = 1, alpha = alpha ,option = option, print.progress = F)

###Output of Running B bootstrap samples------------------------
# This function will obtain number of phony variables selected (U) and total number of variables selected (SP) when 
# considering m phony variable for each threshold level alpha. The results will be a matrix of B x (2*length(alpha)) value
# Input: The same as input for bootstrap function above and
#        B: the number of bootstrap sample
#        ncores: number of cores for embrassingly parallel running

BBootFS <- function(B, ncores, counts, FixCov, VarCov, m, alpha, option = "OWN", print.progress = F){
  USP <- mclapply(1:B, function(nrep)BootFS(counts, FixCov, 
                                             VarCov, m=m, nrep, 
                                             alpha, option = option, 
                                             print.progress), 
                mc.cores = ncores)
  USP <- do.call("rbind", USP)
  list(USP=USP, mPhoCov = m)
}
#Example--------------------------------
# amax <-  15
# alpha <- c(seq(1, amax, length = 800))
# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# option <- "ORX"
# B <- 2
# ncores <- 1
# out <- BBootFS(B,ncores = ncores,  counts = counts[1:10,], FixCov = FixCov, VarCov = VarCov, m = 20, alpha = alpha ,option = option, print.progress = T)

# 
# 

# DataFSOut  <- res$res$WN$DataFSOut
# #  BBootFSOut <- res$res$WN$BBootFSOut
FSRP <- function(DataFSOut, BBootFSOut){
  S <- DataFSOut$S
  alpha <- DataFSOut$alpha
  nAlpha <- length(alpha)
  nVarCov <- length(DataFSOut$FS)
  nPhoCov <- BBootFSOut$mPhoCov
  B <- nrow(BBootFSOut$USP)
  USP <- apply(BBootFSOut$USP, 2, sum)/B
  U <- USP[1:nAlpha]
  SP <- USP[-c(1:nAlpha)]
  ghat.P.RE <- U/(1+SP)
  ghat.P.ER <- U/(1+S)
  out <- list(B = B, nPhoCov = nPhoCov, FS = DataFSOut$FS, alpha = alpha, ghat.P.RE = ghat.P.RE, ghat.P.ER = ghat.P.ER)
  out
  }
# 
# # estimate.kU <- function(FSRPOut){
# #   
# # }
# 
# 
# ## Optimal threshold alpha ---------------
BestAlpha <- function(FSRPOut, gam0=0.05){
  alpha <- FSRPOut$alpha
  FS <- FSRPOut$FS
  kT <- length(FS)
  kP <- FSRPOut$nPhoCov
  ghat.P.RE <- FSRPOut$ghat.P.RE
  alphamin.RE<-min(alpha[ghat.P.RE==max(ghat.P.RE)])  # alpha with largest ghat.RE
  kU.RE <- kT
  gam0.P <- kP*gam0/(kU.RE +kP*gam0)
  repeat{
    ind.RE <- (cummin(ghat.P.RE)<=gam0.P& alpha>=alphamin.RE)*1
    alphahat.RE <- alpha[(which(ind.RE==1))[1]]       # RE est. of alpha
    if(is.na(alphahat.RE)){
      warning(paste0("All Estimated FSRs are greater than ", gam0, ". This happens in BestAlpha function."))
      alphahat.RE <- 10000
    }
    temp <- kT - length(FS[cummin(FS)>=alphahat.RE])
    if(temp==kU.RE) {break}
    kU.RE <- temp
    gam0.P <- kP*gam0/(kU.RE +kP*gam0)
  }

  ghat.P.ER <- FSRPOut$ghat.P.ER
  alphamin.ER<-min(alpha[ghat.P.ER==max(ghat.P.ER)])  # alpha with largest ghat.RE
  kU.ER <- kT
  gam0.P <- kP*gam0/kU.ER
  repeat{
    ind.ER <- (cummin(ghat.P.ER)<=gam0.P& alpha>=alphamin.ER)*1
    alphahat.ER <- alpha[(which(ind.ER==1))[1]]       # ER est. of alpha
    if(is.na(alphahat.ER)){
      warning(paste0("All Estimated FSRs are greater than ", gam0, ". This happens in BestAlpha function."))
      alphahat.ER <- 10000
    }
    temp <- kT - length(FS[cummin(FS)>=alphahat.ER])
    if(temp==kU.ER) {break}
    kU.ER <- temp
    gam0.P <- kP*gam0/(kU.ER)
  }

  c(alphahat.RE  = alphahat.RE , alphahat.ER  = alphahat.ER, kU.RE = kU.RE, kU.ER = kU.ER )
}

## Input will be alpha, output from bootstrap
FSR <- function(FSRPOut, BestAlphaOut, gam0 = .05){
  # FSRPOut <- FSRP(DataFSOut, BBootFSOut)
  # BestAlphaOut <- BestAlpha(FSRPOut, gam0 = gam0)
  ghat.RE<-pmin(1, FSRPOut$ghat.P.RE/(1-FSRPOut$ghat.P.RE)*BestAlphaOut["kU.RE"]/FSRPOut$nPhoCov)      # gammahat_RE, SP-U plays role of S
  ghat<- pmin(1, FSRPOut$ghat.P.ER*BestAlphaOut["kU.ER"]/FSRPOut$nPhoCov)             # gammahat_ER
  out <- cbind(alpha = FSRPOut$alpha, ghat.RE = ghat.RE, ghat = ghat)
  rownames(out) <- NULL
  data.frame(out)
}


#### Function calculate FSR both RE and ER ------------
## Input will be alpha, output from bootstrap
# FSR <- function(DataFSOut, BBootFSOut){
#   S <- DataFSOut$S                
#   alpha <- DataFSOut$alpha
#   nAlpha <- length(alpha)
#   nVarCov <- length(DataFSOut$FS)
#   nPhoCov <- BBootFSOut$mPhoCov
#   B <- nrow(BBootFSOut$USP)
#   USP <- apply(BBootFSOut$USP, 2, sum)/B
#   U <- USP[1:nAlpha]
#   SP <- USP[-c(1:nAlpha)]
#   ghat.RE<-(nVarCov - S)*(U/nPhoCov)/(1+SP-U)        # gammahat_RE, SP-U plays role of S
#   ghat<-(nVarCov - S)*(U/nPhoCov)/(1+S)              # gammahat_ER
#   out <- cbind(alpha = alpha, ghat.RE = ghat.RE, ghat = ghat)
#   rownames(out) <- NULL
#   data.frame(out)
# }
# 

#Example----------------
# amax <-  15
# alpha <- c(seq(1, amax, length = 800))
# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# DataOut <- DataFS(counts = counts[1:10,], FixCov = FixCov, VarCov = VarCov, alpha, print.progress = T)
# option <- "ORX"
# B <- 2
# ncores <- 1
# BootOut <- BBoot(B,ncores = ncores,  counts = counts[1:10,], FixCov = FixCov, VarCov = VarCov, m = ncol(VarCov), alpha = alpha ,option = option, print.progress = T)
# FSROut <- FSR(DataOut, BootOut)

# Plot 2 FSR Estimates----------------
# PlotFSR <- function(FSROut, gam0, FileName){
#   pdf(paste0(FileName, ".pdf"))
#   plot(FSROut$alpha, FSROut$ghat, type="l",lty=2,ylab="Est.FSR",xlab="alpha")
#   lines(FSROut$alpha,FSROut$ghat.RE,col="red")
#   abline(h = gam0)
#   title("Est.FSR.RE (red) and FSR.ER (dash) Curves")
#   dev.off()
# }


PlotFSR <- function(FSROut, gam0, option, B, m){
  FSROut.melt <- melt(FSROut, id.vars ="alpha", measure.vars =  c("ghat", "ghat.RE"))
  names(FSROut.melt) <- c("alpha", "Formula", "Est.FSR")
  levels(FSROut.melt$Formula) <- c("ER", "RE")
  p <- ggplot(FSROut.melt, aes(x=alpha, y=Est.FSR, group = Formula)) +
    geom_line(aes(linetype=Formula, color = Formula), size = .5) +   
    scale_linetype_manual(values=c("dotdash", "solid"))+ 
    scale_color_manual(values=c("red", "blue"))+
    theme(legend.title=element_blank(), axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"))+
    geom_hline(yintercept=gam0) +
    labs(title = paste(paste0(option, ", B=", B, ", m=", m, ", gam0=", gam0)))+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab( expression(alpha))+
    ylab(expression(hat(gamma)))
  p
  }
# FSROut <- FSR(DataFSOut,BBootFSOut)
# PlotFSR(FSROut, gam0 = .05, option = "OWN", B = 240, m = 20)
#Example----------
#PlotFSR(FSROut = FSROut, gam0= .05, FileName = "Test")

## Optimal threshold alpha ---------------
# BestAlpha <- function(FSROut, gam0){
#   alpha <- FSROut$alpha
#   ghat.RE <- FSROut$ghat.RE
#   alphamin.RE<-min(alpha[ghat.RE==max(ghat.RE)])  # alpha with largest ghat.RE
#   ind.RE <- (cummin(ghat.RE)<=gam0& alpha>=alphamin.RE)*1
#   alphahat.RE <- alpha[(which(ind.RE==1))[1]]       # RE est. of alpha
#   
#   ghat <- FSROut$ghat
#   
#   alphamin<-min(alpha[ghat==max(ghat)])     # alpha with largest ghat
#   ind <- (cummin(ghat)<=gam0 & alpha>=alphamin)*1
#   
#   
#   alphahat.ER <- alpha[(which(ind==1))[1]]          # ER est. of alpha
#   if(is.na(alphahat.RE)){
#     warning(paste0("All Estimated FSRs are greater than ", gam0, ". This happens in BestAlpha function."))
#     alphahat.RE <- 10000
#     }
#   if(is.na(alphahat.ER)) {
#     warning(paste0("All Estimated FSRs are greater than ", gam0, ". This happens in BestAlpha function."))
#     alphahat.ER <- 10000
#     }
#   c(alphahat.RE  = alphahat.RE , alphahat.ER  = alphahat.ER )
# }

#Example ---------
# BestAlphaOut <- BestAlpha(FSROut, gam0 = .05)

# Function to calculate the selected covariates using FSR method
BestCov <- function(DataFSOut, BestAlphaOut){
  list(BestER = DataFSOut$FS[cummin(DataFSOut$FS)>=BestAlphaOut["alphahat.ER"]],
       BestRE = DataFSOut$FS[cummin(DataFSOut$FS)>=BestAlphaOut["alphahat.RE"]])
  }

#Example----------
# BestCovOut <- BestCov(DataFSOut, BestAlphaOut)

# Pvalue of all covariates in voom method -----------------
# This function calculate all pvalue of covariates in voom model with 
# Input: counts: counts data
#        AllCov: a data frame containing all covariates considered
# Output:  a dataframe of all p-value of all covariates in AllCov
VoomPv <- function(counts, AllCov){
  dm <- model.matrix(formula(paste0("~", paste0(names(AllCov), collapse = "+"))), data = AllCov)
  colnames(dm)[1] <- "Intercept"
  vout <- voom(counts = counts, design = dm, lib.size = apply(counts, 2, quantile, .75), plot = F)
  fit <- lmFit(vout)
  fit2 <- eBayes(fit)
  pvs <- ldply(1:ncol(AllCov), function(i){
    if(is.factor(AllCov[,i]) | is.character(AllCov[,i])) {
      ct <- paste0(grep(paste0(names(AllCov)[i]),  x = colnames(dm), value = T), collapse = ",  ")
    }else{
      ct <- paste0(grep(paste0(names(AllCov)[i], "$"),  x = colnames(dm), value = T), collapse = ",  ") 
    }
    C.matrix <- eval(parse(text=paste0("makeContrasts(",  ct, ",levels = dm)")))
    fit1 <- contrasts.fit(fit, contrasts =C.matrix)
    fit1 <- eBayes(fit1)
    tt <- topTableF(fit1, sort ="none", n = Inf)
    pv <- tt$P.Value
  })
  pvs <- t(pvs)
  qvs <- apply(pvs, 2, function(x)jabes.q(x))
  colnames(pvs) <- colnames(qvs) <- colnames(AllCov)
  rownames(pvs) <- rownames(qvs) <- rownames(fit$coefficients)
  
  ic <- limmaIC(vout, fit)
  res <- list(y = vout$E, 
              pvs = pvs, #VoomOut  = vout, 
              qvs = qvs,
              Beta = fit$coef, 
              Yhat = fit$coef%*%t(vout$design), 
              lib.size = vout$targets$lib.size, 
              weights = vout$weights, 
              sigma = fit$sigma, 
              s2.post = fit2$s2.post,
              df.prior = fit2$df.prior,
              s2.prior = fit2$s2.prior,
              df.total = fit2$df.total,
              design = vout$design, 
              aic = ic["aic"], bic = ic["bic"])
  res
}


## Function to calculate AIC BIC of output from limma voom-----

limmaIC <- function(vout, fit){
  y <- vout$E
  Yhat <- fit$coef%*%t(vout$design)
  weights <- vout$weights
  sigma <- fit$sigma
  npar <- nrow(vout$design) - ncol(vout$design) + 1
  
  minus2.loglik <- sum(laply(1:nrow(y), function(g){
    sum(log(2*pi) -log(weights[g,])+log(sigma[g]) + weights[g,]^2/sigma[g]^2*(y[g,] - Yhat[g,])^2)
  }))/nrow(y)
  aic <- minus2.loglik +2*npar
  bic <- minus2.loglik +log(nrow(vout$design))*npar
  c(aic = aic, bic = bic)
}


##Example----------------

# AllCov <- cbind(FixCov, VarCov[names(BestCovOut$BestER)])
# VoomPvOut <- VoomPv(counts, AllCov)
# 
# AllCov <- cbind(FixCov, VarCov[names(BestCovOut$BestRE)])
# VoomPvOut <- VoomPv(counts[1:10, ], AllCov)

#Plot histogram of p-values of all covariates in the model
# PlotVoomPv <- function(VoomPvOut, FileName){
#   pdf(paste0(FileName, ".pdf"), width = 20, height = 20)
#   par(mfrow = c(4,4 ))
#   for(i in colnames(VoomPvOut)){
#     hist(VoomPvOut[,i], nclass = 100, main = paste0(i, ", DEG=", sum(jabes.q(VoomPvOut[,i]) <= .05)), 
#          xlab = "pvalue")
#   }
#   dev.off()
# }

PlotVoomPv <- function(VoomPvOut, option= "OWN", B = 100, m = 3, ErrType= "ER", gam0 = 0.05, FDR.level = .05, bin.width = 0.05, title.name){
  VoomPvOut <- VoomPvOut$pvs
  VoomQvOut <- apply(VoomPvOut, 2, function(x)jabes.q(x))
  DEGs <- apply(VoomQvOut <= FDR.level, 2, sum)
  VoomPvOut.melt <- melt(VoomPvOut, measure.vars = .)
  names(VoomPvOut.melt)[2:3] <- c("Covariate", "pvalue")
  levels(VoomPvOut.melt$Covariate)<- paste(levels(VoomPvOut.melt$Covariate),  DEGs[levels(VoomPvOut.melt$Covariate)], sep = ", q.05 = ")
  p <- qplot(pvalue, data = VoomPvOut.melt, geom = "histogram", breaks=seq(0,1,by=bin.width))+ 
    facet_wrap(~ Covariate, scales = "free_y") +labs(title = title.name)+
    theme(plot.title = element_text(hjust = 0.5))
  p
}
#Example--------
# PlotVoomPv(VoomPvOut, option, B, m, ErrType)


# The function shows FSR algorithm for data analysis----------
FSRAnalysis <- function(counts, FixCov, VarCov, option = "OWN", B= 48, m = ncol(VarCov), amax= 15, gam0 = .05, ncores, print.progress){
  cat("Option=", option, ", B=", B, ", m=",m, "\n")
  pm <- proc.time()
  alpha <- c(seq(1, amax, length = 800))
  DataFSOut <- DataFS(counts = counts, FixCov = FixCov, VarCov=VarCov, alpha = alpha, print.progress = print.progress)
  BBootFSOut <- BBootFS(B= B, ncores = ncores, counts = counts, FixCov, VarCov = VarCov[names(DataFSOut$FS)], m= m, alpha = alpha, option = option, print.progress = print.progress)
  # saveRDS(BBootFSOut, file = "BBootFSOut.rds")
  # saveRDS(DataFSOut, file = "DataFSOut.rds")
  FSRPOut <- FSRP(DataFSOut=DataFSOut, BBootFSOut = BBootFSOut)
  BestAlphaOut <- BestAlpha(FSRPOut, gam0)
  FSROut <- FSR(FSRPOut, BestAlphaOut)
  PlotFSROut <- PlotFSR(FSROut = FSROut, gam0 = gam0, option = option, B = B, m = m)
  BestCovOut <- BestCov(DataFSOut = DataFSOut, BestAlphaOut = BestAlphaOut)
  VoomPvOutER <- VoomPv(counts = counts, 
                      AllCov = cbind(FixCov, VarCov[names(BestCovOut$BestER)]))
  PlotVoomPvOutER <-  PlotVoomPv(VoomPvOut = VoomPvOutER, option = option, B = B, m = m, ErrType = "ER", gam0 = gam0)
  VoomPvOutRE <- VoomPv(counts = counts, 
                        AllCov = cbind(FixCov, VarCov[names(BestCovOut$BestRE)]))
  PlotVoomPvOutRE <-  PlotVoomPv(VoomPvOut = VoomPvOutRE, option = option, B = B, m = m, ErrType = "RE", gam0 = gam0)
  res <- list(DataFSOut = DataFSOut, BBootFSOut = BBootFSOut, 
              FSROut = FSROut, PlotFSROut = PlotFSROut, 
              BestAlphaOut = BestAlphaOut, 
              BestCovOut = BestCovOut, 
              VoomPvOutER = VoomPvOutER, PlotVoomPvOutER = PlotVoomPvOutER, 
              VoomPvOutRE = VoomPvOutRE, PlotVoomPvOutRE = PlotVoomPvOutRE, 
              option = option, B = B, m = m)
  pm <- proc.time()-pm
  cat("Option=", option, ", B=", B, ", m=",m,  " runs in ", pm[3], " s\n")
  res
}

#Example------------

# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# option <- c("WN", "OWN", "RX", "ORX")
# B <- 2
# m <- ncol(VarCov)
# amax <- 15
# gam0 <- .05
# ncores <- 2
# option <- "OWN"
# print.progress <- F
# FSRAnalysisOut <- FSRAnalysis(counts[1:10,], FixCov, VarCov,
#                                             option = option, B= B, m = m,
#                                             amax= amax, gam0 = gam0, ncores=ncores,
#                                             print.progress=print.progress)

#Function to run FSR Analysis for all 4 options
FSRAnalysisAll <- function(counts, FixCov, VarCov, B, m , amax, gam0, ncores, print.progress){
  option <- c("RX", "ORX", "WN", "OWN")
  
  res <- llply(option, function(i)FSRAnalysis(counts, FixCov, VarCov, 
                                              option = i, B= B, m = m, 
                                              amax= amax, gam0 = gam0, ncores=ncores, 
                                              print.progress=print.progress))
  names(res) <- option
  out <- list(res = res, FixCov = FixCov, VarCov = VarCov)
  out
}
#Example---------------------------
# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# B <- 2
# m <- ncol(VarCov)
# amax <- 15
# gam0 <- .05
# ncores <- 2
# print.progress <- F
# FSRAnalysisAllOut <- FSRAnalysisAll(counts[1:10,], FixCov, VarCov,B= B, m = m,
#                                             amax= amax, gam0 = gam0, ncores=ncores,
#                                             print.progress=print.progress)


# A function wrap all things for FSRAnalysisAll for real data ----------------
FSRAnalysisAllSave <- function(counts, FixCov, VarCov, B, m, amax, gam0, ncores, print.progress){
  RealDataOutPath <- paste0("RealDataOut")
  dir.create(path = RealDataOutPath, showWarnings = F)
  RealDataOut <- FSRAnalysisAll(counts, FixCov, VarCov, B, m, amax, gam0, ncores, print.progress)
  saveRDS(RealDataOut, file = paste0(RealDataOutPath, "/",RealDataOutPath,  ".BData", B, ".mData", m, ".rds"))
  RealDataOut
}



# A Function to obtain names of selected covariates in each option. ---------

#Input: FSRAnalysisAllOut: is the output from FSRAnalysisAll function
#       option: The option to generate phony variables (WN, OWN, RX, ORX) 
#       ErrType: Type of FSR formula  used (ER, or RE)  
#Output: a list of 4 components. Each component is a list of 2 other components whose
#        one is a vector of name of selected covariates with FSR.ER, 
#        the other is a vector of name of selected covariates with FSR.RE

FSRSelCov <- function(FSRAnalysisAllOut, option, ErrType ){
  SelCovName <- names(((FSRAnalysisAllOut[[1]])[[option]]$BestCovOut)[[paste0("Best", ErrType)]])
  SelCovName
}


#Example---------------
#Example---------------------------
# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# B <- 2
# m <- ncol(VarCov)
# amax <- 15
# gam0 <- .05
# ncores <- 2
# print.progress <- F
# FSRAnalysisAllOut <- FSRAnalysisAll(counts[1:10,], FixCov, VarCov,B= B, m = m,
#                                     amax= amax, gam0 = gam0, ncores=ncores,
#                                     print.progress=print.progress)
# 
# FSRSelCovOut <- FSRSelCov(FSRAnalysisAllOut, option = "WN", ErrType = "ER")
# FSRSelCovOut



#Function to simulate count data from Real Data Analysis----------------
# This function create sim count data with nGene from the selected model
# of real data analysis
# Input: RealDataOut: is the list of final result from real data analysis
#        nGene: number of genes in the simulated data
#        option: the option that is used to generate phony variables in real data analysis 
#                option is one of WN, OWN, RX, ORX
#        ErrType: either "ER" or "RE" which is character showing formula of calculating estimated FSR
# SimCounts <- function(RealDataOut, nGene, option, ErrType="ER", nrep){
#   set.seed(2017+nrep)
#   res <- ((RealDataOut[[1]])[[option]])[[paste0("VoomPvOut",ErrType)]]
#   lib.size <- res$lib.size; design <- res$design
#   sigma <- res$sigma; Yhat0 <- res$Yhat; weights <- res$weights
#   VarMat0 <- sigma^2*1/weights
#   IndSample <- sample(1:nrow(Yhat0), size = nGene)
#   Yhat <- Yhat0[IndSample,]; VarMat <- VarMat0[IndSample,]
#   sim.counts <- vapply(1:nrow(Yhat), function(i){
#     y <- Yhat0[i,]+ mvrnorm(n = 1, mu = rep(0, ncol(Yhat0)), Sigma = diag(VarMat[i,])) # library(MASS)
#     unname(round(2^y*(lib.size+1)/10^6 )[ , drop = T])
#   }, FUN.VALUE = rep(1.0, ncol(Yhat0)))
#   t(sim.counts)
# }


## Another simulation setting--------------

## Previously, in the simulation setting, it looks like all covairates are significant
# for every genes. In realty, this may not be true, i.e, some covariates just
# significant for a particular subset of genes, while the other is not. 
# In this simulation setting, we will construct a simulation design
# that each covariates just significant with respect to a proportion of genes
# The number of DEGs with respect to each covariates is estimated by the method 
# estimate.m0 by Nettleton. Therefore, in the simulation, we will simulate 
# count data as follow. 
# First we will modifiy the regression coefficient from Real Data. 
# (m-m0_covariate) The regression coefficient for each covariate will be the same
# m0_covariate regression coefficient for each covariate will set to be equal to zero
# The assignment are based on ranking of the p-value for each covariates
# After this, we fix new regression coefficient. From this 
# regression coefficient, we calculate Yhat as before. 
# Then the simulation process keep the same by generating randomly 5000 genes!!!
#not use full here --------------------------
# YhatFunc <- function(RealDataOut, SimType = 0){
#   EEGene <- list()
#   res <- ((RealDataOut[[1]])[[option]])[[paste0("VoomPvOut",ErrType)]]
#   design <- res$design
#   G0 <- apply(res$pvs,2, function(x)estimate.m0(x));G <- nrow(res$pvs)
#   Beta <- res$Beta
#   AllCov <- cbind(RealDataOut$FixCov, RealDataOut$VarCov)
#   ConCov <- AllCov[,colnames(res$pvs)]
#   m0 <- apply(res$pvs, 2, function(x)estimate.m0(x))
#   for(i in 1:ncol(ConCov)){
#     if(is.factor(ConCov[,i]) | is.character(ConCov[,i])) {
#       ct <- grep(paste0(names(ConCov)[i]),  x = colnames(design), value = T)
#     }else{
#       ct <- grep(paste0(names(ConCov)[i], "$"),  x = colnames(design), value = T)
#     }
#     EEGene[[i]] <- (sort(res$pvs[,i], decreasing = T,  index.return = T))$ix[1:m0[i]]
#     Beta[EEGene[[i]], ct] <- 0
#   }
#   names(EEGene) <- colnames(res$pvs)
#   Yhat0 <- Beta%*%t(design)
#   if(SimType==0){Yhat0 <- res$Yhat; EEGene <- NULL}
#   return(list(Yhat0 = Yhat0, EEGene = EEGene))
# }


# SimCounts <- function(RealDataOut, nGene, option, ErrType="ER", nrep, SimType = 0){
#   set.seed(2017+nrep)
#   res <- ((RealDataOut[[1]])[[option]])[[paste0("VoomPvOut",ErrType)]]
#   lib.size <- res$lib.size; design <- res$design
#   sigma <- res$sigma;  weights <- res$weights; VarMat0 <- sigma^2*1/weights
#   YhatFuncOut <- YhatFunc(RealDataOut=RealDataOut,SimType =SimType)
#   Yhat0 <- YhatFuncOut$Yhat0
#   IndSample <- sample(1:nrow(Yhat0), size = nGene)
#   Yhat <- Yhat0[IndSample,]; VarMat <- VarMat0[IndSample,]
#   sim.counts <- vapply(1:nrow(Yhat), function(i){
#     y <- Yhat0[i,]+ mvrnorm(n = 1, mu = rep(0, ncol(Yhat0)), Sigma = diag(VarMat[i,])) # library(MASS)
#     unname(round(2^y*(lib.size+1)/10^6 )[ , drop = T])
#   }, FUN.VALUE = rep(1.0, ncol(Yhat0)))
#   SimCnt <- t(sim.counts)
#   list(SimCnt = SimCnt, YhatFuncOut = YhatFuncOut, IndSample = IndSample)
# }

##end not use full----------------
# SimCounts <- function(RealDataOut, nGene, option, ErrType="ER", nrep, SimType = 0){
#   set.seed(2017+nrep)
#   res <- ((RealDataOut[[1]])[[option]])[[paste0("VoomPvOut",ErrType)]]
#   lib.size <- res$lib.size; design <- res$design
#   sigma <- res$sigma;  weights <- res$weights; VarMat0 <- sigma^2*1/weights
SimCounts <- function(RealDataOut, nGene, option, ErrType="ER", nrep){
  set.seed(2017+nrep)
  res <- ((RealDataOut[[1]])[[option]])[[paste0("VoomPvOut",ErrType)]]
  EEGene <- list()
  AllCov <- cbind(RealDataOut$FixCov, RealDataOut$VarCov)
  ConCov <- AllCov[colnames(res$pvs)]
  lib.size <- res$lib.size; design <- res$design; Beta0 <- res$Beta; 
  sigma <- res$sigma;  weights <- res$weights; VarMat0 <- sigma^2*1/weights
  m0 <- apply(res$pvs, 2, function(x)estimate.m0(x))
  for(i in 1:ncol(ConCov)){
    if(is.factor(ConCov[,i]) | is.character(ConCov[,i])) {
      ct <- grep(paste0(names(ConCov)[i]),  x = colnames(design), value = T)
    }else{
      ct <- grep(paste0(names(ConCov)[i], "$"),  x = colnames(design), value = T)
    }
    EEGene[[i]] <- (sort(res$pvs[,i], decreasing = T,  index.return = T))$ix[1:m0[i]]
    Beta0[EEGene[[i]], ct] <- 0
  }
  names(EEGene) <- colnames(res$pvs)
  IndSample <- sample(1:nrow(Beta0), size = nGene)
  Beta <- Beta0[IndSample,]; VarMat <- VarMat0[IndSample,]
  sim.counts <- vapply(1:nrow(Beta), function(i){
    y <- design%*%Beta[i,]+ mvrnorm(n = 1, mu = rep(0, nrow(design)), Sigma = diag(VarMat[i,])) # library(MASS)
    unname(round(2^y*(lib.size+1)/10^6 )[ , drop = T])
  }, FUN.VALUE = rep(1.0, nrow(design)))
  SimCnt <- t(sim.counts)
  list(SimCnt = SimCnt, IndSample = IndSample, EEGene = EEGene)
}


#Example------------------
# RealDataOut <- readRDS("RFIAnalysis_Result.rds")
# nGene <- 100
# option <- "OWN"
# ErrType <- "ER"
# nrep <- 1
# Sc <- SimCounts(RealDataOut, nGene, option, ErrType, nrep)

# Function to run 1 simulation  for 1 option--------------------
FSR1Sim <- function(RealDataOut, nGene, option, ErrType, nrep, B, m, amax, gam0, ncores, print.progress, SimPath){
  cat("Sim = ", nrep, "\n")
  FixCov <- RealDataOut$FixCov
  VarCov <- RealDataOut$VarCov
  SimCnt <- SimCounts(RealDataOut, nGene, option,ErrType,nrep)
  # saveRDS(SimCnt, file = "SimCnt.rds")
  SimCntOut <- FSRAnalysisAll(counts= SimCnt$SimCnt, FixCov, VarCov, B, m, amax, gam0, ncores, print.progress)
  SimSelCov.ER <- llply(c("WN", "OWN", "RX", "ORX"),
                     function(option)FSRSelCov(FSRAnalysisAllOut = SimCntOut, option = option, ErrType = "ER"))
  names(SimSelCov.ER) <- c("WN.ER", "OWN.ER", "RX.ER", "ORX.ER")
  SimSelCov.RE <- llply(c("WN", "OWN", "RX", "ORX"),
                        function(option)FSRSelCov(FSRAnalysisAllOut = SimCntOut, option = option, ErrType = "RE"))
  names(SimSelCov.RE) <- c("WN.RE", "OWN.RE", "RX.RE", "ORX.RE")
  SimOut <- list(SimCnt = SimCnt, SimCntOut = SimCntOut, SimSelCov.ER = SimSelCov.ER, SimSelCov.RE = SimSelCov.RE)
  saveRDS(SimOut, file = paste0(SimPath, "/",SimPath, ".nrep", nrep, ".rds"))
  SimSelCov <- list(SimSelCov.ER = SimSelCov.ER, SimSelCov.RE = SimSelCov.RE)
  SimSelCov
}

# Function to run length(nSim) simulation from 1 option---------------------
# Input: nSim: is vector  like 1:10 indicating the where max(nSim) is the number of simulated data
#        RealDataOut: is the output from real data analysis from FSR method
#        nGene: number of genes in simulated dataset
#        option: the option to generate phony variables: "WN", "OWN", "RX", "ORX"
#        ErrType: Type of FSR estimator formula, either "ER" or "RE"
#        B: number of bootstrap samples during running FS for model including phony variables
#        m: number of phony variables added
#        amax: alpha max =15
#        gam0: = 0.05
#        ncores: number of cores used
#        print.progress: printing the progress of FS

# Output: a list of 2 components. First component is the vector of True Selected Covariates from RFI data
#         second component is a list of length(nSim) components, each of which corresponding to the set of 
#         selected covariate for the corresponding simulated data.

FSRnSim <- function(RealDataOut, nGene, option, ErrType, nSim, B, m, amax, gam0, ncores, print.progress, SimPath){
  TrueSelCov <- FSRSelCov(FSRAnalysisAllOut=RealDataOut, option=option, ErrType = ErrType)
  nSimSelCov <- ldply(nSim, function(nrep){
    SimSelCov <- FSR1Sim(RealDataOut, nGene, option, ErrType, nrep, B, m, amax, gam0, ncores, print.progress, SimPath)
    
    S.ER <-  laply(SimSelCov$SimSelCov.ER, function(x) length(intersect(TrueSelCov, x)))
    R.ER <-  laply(SimSelCov$SimSelCov.ER, function(x)length(x))
    U.ER <- R.ER - S.ER
    # FSP <- U/pmax(R, 1)
    FSP.ER <- U.ER/(R.ER+1) # i
    res1 <- c(S.ER, R.ER, U.ER, FSP.ER)
    names(res1) <- paste(rep(c("S", "R", "U", "FSP"), each = 4), rep(names(SimSelCov$SimSelCov.ER), 4), sep = ".")
    S.RE <-  laply(SimSelCov$SimSelCov.RE, function(x) length(intersect(TrueSelCov, x)))
    R.RE <-  laply(SimSelCov$SimSelCov.RE, function(x)length(x))
    U.RE <- R.RE - S.RE
    # FSP <- U/pmax(R, 1)
    FSP.RE <- U.RE/(R.RE+1) # i
    res2 <- c(S.RE, R.RE, U.RE, FSP.RE)
    names(res2) <- paste(rep(c("S", "R", "U", "FSP"), each = 4), rep(names(SimSelCov$SimSelCov.RE), 4), sep = ".")
    res <- c(res1, res2)
    res
    })
  res <- list(TrueSelCov = TrueSelCov, nSimSelCov = nSimSelCov)
  res
  }
#Example------------
# RealDataOut <- readRDS("RFIAnalysis_Result.rds")
# nGene <- 10
# option <- "WN"
# ErrType <- "ER"
# nSim <- 1:2
# B <- 1
# m <- 13
# amax <- 15
# gam0 <- 0.05
# ncores <- 1
# print.progress <- F
# FSRSimOut <- FSRnSim(RealDataOut, nGene, option, ErrType,nSim,B,m,amax,gam0,ncores,print.progress)


#Function wrapping up all simumation result ------

FSRnSimAll <- function(BData, mData, nGene, option, ErrType, nSim, BSim, mSim, amax, gam0, ncores, print.progress){
  RealDataOut <- readRDS(file = paste0("RealDataOut/RealDataOut", ".BData", BData, ".mData", mData, ".rds"))
  SimPath <- paste0("Sim_", "BData", BData, "_mData", mData, "_",option, "_gamma", ErrType, "_gam0_", gam0, "_nGene", nGene, "_Bsim", BSim, "_mSim", mSim)
  dir.create(path = SimPath, showWarnings = F)
  FSRSimOut <- FSRnSim(RealDataOut, nGene, option, ErrType,nSim,B=BSim,m=mSim,amax,gam0,ncores,print.progress, SimPath)
  saveRDS(FSRSimOut, file = paste0(SimPath, "/",SimPath, "_nSim", min(nSim), "_", max(nSim), ".rds"))
  FSRSimOut
}


# ##BackWard Selection Using Adding PseudoVariables--------------
# bs <- function(counts, FixCov, VarCov, print.progress = F){
#   if(!is.data.frame(FixCov)|!is.data.frame(VarCov)) stop("FixCov and VarCov have to be data frame with different column names")
#   DelCov <- VarCov[,0] # Initiaize Set of Deleted Covariates
#   ConCov <- VarCov # Initialize Set of Considered Covariates
#   PvList <- list() # List of P-values for the covariate selected at each iteration
#   WorstP5 <- list() # Number of p-value less than 0.05 for the covariate selected at each iteration
#   if(ncol(ConCov)!=0){
#     Iter <- 1
#     repeat{
#       if(print.progress){
#         cat("----------------------------------------\n")
#         cat("Iteration = ", Iter, "\n")
#         cat("DelCov = ", names(DelCov), "\n")
#       }
#       WorstP5[[Iter]] <- Inf
#       AllCov <- cbind(FixCov, ConCov)
#       dm <- model.matrix(formula(paste0("~", paste0(names(AllCov), collapse = "+"))), data = AllCov)
#       colnames(dm)[1] <- "Intercept"
#       # if(!is.fullrank(dm)|ncol(dm)==nrow(dm)){ Ind <- T;
#       # warning("Sigularity design matrix! Check generated pseudo-variables!"); break}
#       vout <- voom(counts = counts, design = dm, lib.size = apply(counts, 2, quantile, .75), plot = F)
#       suppressWarnings(fit <- lmFit(vout))
#       for(i in 1:ncol(ConCov)){
#         if(is.factor(ConCov[,i]) | is.character(ConCov[,i])) {
#           ct <- paste0(grep(paste0(names(ConCov)[i]),  x = colnames(dm), value = T), collapse = ",  ")
#         }else{
#           ct <- paste0(grep(paste0(names(ConCov)[i], "$"),  x = colnames(dm), value = T), collapse = ",  ") 
#         }
#         C.matrix <- eval(parse(text=paste0("makeContrasts(",  ct, ",levels = dm)")))
#         fit1 <- contrasts.fit(fit, contrasts =C.matrix)
#         fit1 <- eBayes(fit1)
#         tt <- topTableF(fit1, sort ="none", n = Inf)
#         pv <- tt$P.Value
#         P5 <- max(sum(pv<=.05, na.rm = T), 1)/max(sum(pv>.75, na.rm = T)/5, 1) # sometime pv has missing value NA make P5 cannot calculable
#         if(print.progress)cat("ConCov = ", names(ConCov)[i], ", P5 = ", P5,  "\n")
#         if(WorstP5[[Iter]] > P5) {
#           WorstCov <- i
#           WorstP5[[Iter]] <- P5
#           PvList[[Iter]] <- pv
#           names(PvList)[Iter] <- names(WorstP5)[Iter] <- names(ConCov)[i] 
#         }
#       }
#       Iter <- Iter +1
#       DelCov <- cbind(DelCov, ConCov[WorstCov])
#       ConCov <- ConCov[-WorstCov]
#       if(ncol(ConCov)==0)break
#     }
#     WorstP5 <- do.call("c", WorstP5)
#   }else{WorstP5 <- numeric(0)}
#   # P <- data.frame(t(do.call("rbind", PvList)))
#   WorstP5
# }
# 
# 
# # Backward Selection for Real Data -----------
# # This function running forward selection for real data, the output will be the list
# # of  2 components whose first component is the covariates entering the model in the sequence together
# # with its number of p-values ratio,  and the second component is a vector of number of covariates 
# # entering the model when taking IV threshold alpha
# 
# # Input: counts: counts dataset
# #        FixCov: a dataframe of fix covariates that are not subjected to variable selection
# #        VarCov: a dataframe of flexible covariates that are subjected to variable selection
# #        alpha:  a numeric vector of threshold values
# #        print.progress: a logical variable either T or F to indicate that if we want to print the running process of fs function or not
# # Output: is a list of 2 components: 
# #        fsr: a vector of variable importance measurement in the order of the covariate entering the model in fs process. 
# #             length of fsr equals number of column of VarCov
# #        S: a vector of the number of selected variates at different threshold alpha. 
# DataBS_sva <- function(counts, FixCov, VarCov, alpha, print.progress = F){
#   
#   BS <- bs(counts = counts, FixCov = FixCov, VarCov = VarCov, print.progress = print.progress)
#   S <- vapply(1:length(alpha), function(i)sum(cummax(BS) >alpha[i]), FUN.VALUE = 1.0)
#   return(list(BS = BS, S = S, alpha = alpha))
# }
# 
# #Example--------
# # pm <- proc.time()
# # amax <-  10
# # alpha <- c(seq(1, amax, length = 800))
# # FixCov <- covset["Line"]
# # VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# # DataBSOut <- DataBS(counts = counts, FixCov = FixCov, VarCov = VarCov, alpha, print.progress = T)
# # proc.time() - pm
# # out$S
# 
# 
# 
# ### Bootstrap function to generate data and calculate U, SP for the simulatied data-------
# ## Input: nrep: seed of the bootstrap
# ##        option: either white noise (WN), permutation of row of X (RX), or orthogonal version (OWN, ORX)
# ##        FixCov: Fix covariates that are not subjected to variable selection
# ##        VarCov:  Flexible covariates that are subjected to variable selection
# ##        m: number of phony variables need to generated
# ##        counts: counts data with G genes and s columns
# ##        alpha: a vector of threshold of VI measurement
# # cor(covset[c("RFI", "RINa", "RINb", "Conca", "Concb", "lymp", 
# #              "neut", "mono", "baso", "eosi")])
# BootBS <- function(counts, FixCov, VarCov, m=3, nrep, alpha, option = "OWN", print.progress = T){
#   set.seed(nrep)
#   # cat("nrep = ", nrep, "\n")
#   AllCov0 <- cbind(FixCov, VarCov)
#   dm0 <- model.matrix(formula(paste0("~", paste0(names(AllCov0), collapse = "+"))), data = AllCov0)
#   colnames(dm0)[1] <- "Intercept"
#  repeat{
#     Nois <- PhoVar(dm0 = dm0, m = m, option = option)
#     VarCov2 <- cbind(VarCov, Nois)
#     AllCov1 <- cbind(FixCov, VarCov2)
#     dm1 <- model.matrix(formula(paste0("~", paste0(names(AllCov1), collapse = "+"))), data = AllCov1)
#     if(is.fullrank(dm1)) break
#       }
#   pv <- bs(counts = counts, FixCov = FixCov, VarCov = VarCov2, print.progress = print.progress)
#   pvmax <- cummax(pv)
#   vor <- names(pv)
#   U <- vapply(1:length(alpha), function(i)sum(vor%in%names(Nois) & pvmax>alpha[i]), FUN.VALUE = 1.0)
#   SP <- vapply(1:length(alpha), function(i)sum( pvmax>alpha[i]), FUN.VALUE = 1.0)
#   res <- c(U = U, SP = SP)
#   cat("nrep = ", nrep, "\n")
#   res
# }
# 
# 
# #Example---------------------
# # amax <-  10
# # alpha <- c(seq(1, amax, length = 800))
# # FixCov <- covset["Line"]
# # VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# # option <- "ORX"
# # out <- BootFS(counts = counts[1:50,], FixCov = FixCov, VarCov = VarCov, m = 20, nrep = 1, alpha = alpha ,option = option, print.progress = F)
# 
# ###Output of Running B bootstrap samples------------------------
# # This function will obtain number of phony variables selected (U) and total number of variables selected (SP) when 
# # considering m phony variable for each threshold level alpha. The results will be a matrix of B x (2*length(alpha)) value
# # Input: The same as input for bootstrap function above and
# #        B: the number of bootstrap sample
# #        ncores: number of cores for embrassingly parallel running
# 
# BBootBS <- function(B, ncores, counts, FixCov, VarCov, m, alpha, option = "OWN", print.progress = F){
#   USP <- mclapply(1:B, function(nrep)BootBS(counts, FixCov, 
#                                             VarCov, m=m, nrep, 
#                                             alpha, option = option, 
#                                             print.progress), 
#                   mc.cores = ncores)
#   USP <- do.call("rbind", USP)
#   list(USP=USP, mPhoCov = m)
# }
# #Example--------------------------------
# # pm <- proc.time()
# # amax <-  10
# # alpha <- c(seq(1, amax, length = 800))
# # FixCov <- covset["Line"]
# # VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# # option <- "OWN"
# # B <- 40
# # ncores <- 4
# # m <- 3
# # BBootBSOut <- BBootBS(B,ncores = ncores,  counts = counts, FixCov = FixCov, VarCov = VarCov, m = 3, alpha = alpha ,option = option, print.progress = T)
# # proc.time() - pm
# # 
# # 
# 
# # DataFSOut  <- res$res$WN$DataFSOut
# # #  BBootFSOut <- res$res$WN$BBootFSOut
# FSRPBS <- function(DataBSOut, BBootBSOut){
#   S <- DataBSOut$S
#   alpha <- DataBSOut$alpha
#   nAlpha <- length(alpha)
#   nVarCov <- length(DataBSOut$BS)
#   nPhoCov <- BBootBSOut$mPhoCov
#   B <- nrow(BBootBSOut$USP)
#   USP <- apply(BBootBSOut$USP, 2, sum)/B
#   U <- USP[1:nAlpha]
#   SP <- USP[-c(1:nAlpha)]
#   ghat.P.RE <- U/(1+SP)
#   ghat.P.ER <- U/(1+S)
#   out <- list(B = B, nPhoCov = nPhoCov, BS = DataBSOut$BS, alpha = alpha, ghat.P.RE = ghat.P.RE, ghat.P.ER = ghat.P.ER)
#   out
# }
# # 
# # # estimate.kU <- function(FSRPOut){
# # #   
# # # }
# # 
# # 
# # FSRPBSOut <- FSRPBS(DataBSOut,BBootBSOut)
# 
# # ## Optimal threshold alpha ---------------
# BestAlphaBS <- function(FSRPBSOut, gam0=0.05){
#   alpha <- FSRPBSOut$alpha
#   BS <- FSRPBSOut$BS
#   kT <- length(BS)
#   kP <- FSRPBSOut$nPhoCov
#   ghat.P.RE <- FSRPBSOut$ghat.P.RE
#   alphamin.RE<-min(alpha[ghat.P.RE==max(ghat.P.RE)])  # alpha with largest ghat.RE
#   kU.RE <- kT
#   gam0.P <- kP*gam0/(kU.RE +kP*gam0)
#   repeat{
#     ind.RE <- (cummin(ghat.P.RE)<=gam0.P& alpha>=alphamin.RE)*1
#     alphahat.RE <- alpha[(which(ind.RE==1))[1]]       # RE est. of alpha
#     if(is.na(alphahat.RE)){
#       warning(paste0("All Estimated FSRs are greater than ", gam0, ". This happens in BestAlpha function."))
#       alphahat.RE <- 10000
#     }
#     temp <- kT - length(BS[cummax(BS)>=alphahat.RE])
#     if(temp==kU.RE) {break}
#     kU.RE <- temp
#     gam0.P <- kP*gam0/(kU.RE +kP*gam0)
#   }
#   
#   ghat.P.ER <- FSRPBSOut$ghat.P.ER
#   alphamin.ER<-min(alpha[ghat.P.ER==max(ghat.P.ER)])  # alpha with largest ghat.RE
#   kU.ER <- kT
#   gam0.P <- kP*gam0/kU.ER
#   repeat{
#     ind.ER <- (cummin(ghat.P.ER)<=gam0.P& alpha>=alphamin.ER)*1
#     alphahat.ER <- alpha[(which(ind.ER==1))[1]]       # ER est. of alpha
#     if(is.na(alphahat.ER)){
#       warning(paste0("All Estimated FSRs are greater than ", gam0, ". This happens in BestAlpha function."))
#       alphahat.ER <- 10000
#     }
#     temp <- kT - length(BS[cummax(BS)>=alphahat.ER])
#     if(temp==kU.ER) {break}
#     kU.ER <- temp
#     gam0.P <- kP*gam0/(kU.ER)
#   }
#   
#   c(alphahat.RE  = alphahat.RE , alphahat.ER  = alphahat.ER, kU.RE = kU.RE, kU.ER = kU.ER )
# }
# 
# # BestAlphaOut <- BestAlphaBS(FSRPBSOut)
# 
# 
# ## Input will be alpha, output from bootstrap
# FSR <- function(FSRPOut, BestAlphaOut, gam0 = .05){
#   # FSRPOut <- FSRP(DataFSOut, BBootFSOut)
#   # BestAlphaOut <- BestAlpha(FSRPOut, gam0 = gam0)
#   ghat.RE<-pmin(1, FSRPOut$ghat.P.RE/(1-FSRPOut$ghat.P.RE)*BestAlphaOut["kU.RE"]/FSRPOut$nPhoCov)      # gammahat_RE, SP-U plays role of S
#   ghat<- pmin(1, FSRPOut$ghat.P.ER*BestAlphaOut["kU.ER"]/FSRPOut$nPhoCov)             # gammahat_ER
#   out <- cbind(alpha = FSRPOut$alpha, ghat.RE = ghat.RE, ghat = ghat)
#   rownames(out) <- NULL
#   data.frame(out)
# }
# 
# 
# #### Function calculate FSR both RE and ER ------------
# ## Input will be alpha, output from bootstrap
# # FSR <- function(DataFSOut, BBootFSOut){
# #   S <- DataFSOut$S                
# #   alpha <- DataFSOut$alpha
# #   nAlpha <- length(alpha)
# #   nVarCov <- length(DataFSOut$FS)
# #   nPhoCov <- BBootFSOut$mPhoCov
# #   B <- nrow(BBootFSOut$USP)
# #   USP <- apply(BBootFSOut$USP, 2, sum)/B
# #   U <- USP[1:nAlpha]
# #   SP <- USP[-c(1:nAlpha)]
# #   ghat.RE<-(nVarCov - S)*(U/nPhoCov)/(1+SP-U)        # gammahat_RE, SP-U plays role of S
# #   ghat<-(nVarCov - S)*(U/nPhoCov)/(1+S)              # gammahat_ER
# #   out <- cbind(alpha = alpha, ghat.RE = ghat.RE, ghat = ghat)
# #   rownames(out) <- NULL
# #   data.frame(out)
# # }
# # 
# 
# #Example----------------
# # amax <-  15
# # alpha <- c(seq(1, amax, length = 800))
# # FixCov <- covset["Line"]
# # VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# # DataOut <- DataFS(counts = counts[1:10,], FixCov = FixCov, VarCov = VarCov, alpha, print.progress = T)
# # option <- "ORX"
# # B <- 2
# # ncores <- 1
# # BootOut <- BBoot(B,ncores = ncores,  counts = counts[1:10,], FixCov = FixCov, VarCov = VarCov, m = ncol(VarCov), alpha = alpha ,option = option, print.progress = T)
# # FSROut <- FSR(DataOut, BootOut)
# 
# # Plot 2 FSR Estimates----------------
# # PlotFSR <- function(FSROut, gam0, FileName){
# #   pdf(paste0(FileName, ".pdf"))
# #   plot(FSROut$alpha, FSROut$ghat, type="l",lty=2,ylab="Est.FSR",xlab="alpha")
# #   lines(FSROut$alpha,FSROut$ghat.RE,col="red")
# #   abline(h = gam0)
# #   title("Est.FSR.RE (red) and FSR.ER (dash) Curves")
# #   dev.off()
# # }
# 
# # FSROut <- FSR(FSRPOut = FSRPBSOut, BestAlphaOut)
# 
# # PlotFSR(FSROut = FSROut, gam0 = .05, option = option, B = B, m=m)
# # FSROut <- FSR(DataFSOut,BBootFSOut)
# # PlotFSR(FSROut, gam0 = .05, option = "OWN", B = 240, m = 20)
# #Example----------
# #PlotFSR(FSROut = FSROut, gam0= .05, FileName = "Test")
# 
# ## Optimal threshold alpha ---------------
# # BestAlpha <- function(FSROut, gam0){
# #   alpha <- FSROut$alpha
# #   ghat.RE <- FSROut$ghat.RE
# #   alphamin.RE<-min(alpha[ghat.RE==max(ghat.RE)])  # alpha with largest ghat.RE
# #   ind.RE <- (cummin(ghat.RE)<=gam0& alpha>=alphamin.RE)*1
# #   alphahat.RE <- alpha[(which(ind.RE==1))[1]]       # RE est. of alpha
# #   
# #   ghat <- FSROut$ghat
# #   
# #   alphamin<-min(alpha[ghat==max(ghat)])     # alpha with largest ghat
# #   ind <- (cummin(ghat)<=gam0 & alpha>=alphamin)*1
# #   
# #   
# #   alphahat.ER <- alpha[(which(ind==1))[1]]          # ER est. of alpha
# #   if(is.na(alphahat.RE)){
# #     warning(paste0("All Estimated FSRs are greater than ", gam0, ". This happens in BestAlpha function."))
# #     alphahat.RE <- 10000
# #     }
# #   if(is.na(alphahat.ER)) {
# #     warning(paste0("All Estimated FSRs are greater than ", gam0, ". This happens in BestAlpha function."))
# #     alphahat.ER <- 10000
# #     }
# #   c(alphahat.RE  = alphahat.RE , alphahat.ER  = alphahat.ER )
# # }
# 
# #Example ---------
# # BestAlphaOut <- BestAlpha(FSROut, gam0 = .05)
# 
# # Function to calculate the selected covariates using FSR method
# BestCovBS <- function(DataBSOut, BestAlphaOut){
#   list(BestER = DataBSOut$BS[cummax(DataBSOut$BS)>=BestAlphaOut["alphahat.ER"]],
#        BestRE = DataBSOut$BS[cummax(DataBSOut$BS)>=BestAlphaOut["alphahat.RE"]])
# }
# 
# # BestCovBSOut <- BestCovBS(DataBSOut, BestAlphaOut)
# # 
# # BestCovOut
# 
# #Example--------
# # PlotVoomPv(VoomPvOut, option, B, m, ErrType)
# 
# 
# # The function shows FSR algorithm for data analysis----------
# FSRAnalysisBS <- function(counts, FixCov, VarCov, option = "OWN", B= 100, m = 3, amax= 10, gam0 = .05, ncores, print.progress){
#   cat("Option=", option, ", B=", B, ", m=",m, ", gam0 = ", gam0, "\n")
#   pm <- proc.time()
#   alpha <- c(seq(1, amax, length = 800))
#   DataBSOut <- DataBS(counts = counts, FixCov = FixCov, VarCov=VarCov, alpha = alpha, print.progress = print.progress)
#   BBootBSOut <- BBootBS(B= B, ncores = ncores, counts = counts, FixCov, VarCov, m= m, alpha = alpha, option = option, print.progress = print.progress)
#   # saveRDS(BBootBSOut, file = "BBootBSOut.rds")
#   # saveRDS(DataBSOut, file = "DataBSOut.rds")
#   FSRPBSOut <- FSRPBS(DataBSOut=DataBSOut, BBootBSOut = BBootBSOut)
#   BestAlphaOut <- BestAlphaBS(FSRPBSOut, gam0)
#   FSROut <- FSR(FSRPBSOut, BestAlphaOut)
#   PlotFSROut <- PlotFSR(FSROut = FSROut, gam0 = gam0, option = option, B = B, m = m)
#   BestCovOut <- BestCovBS(DataBSOut = DataBSOut, BestAlphaOut = BestAlphaOut)
#   VoomPvOutER <- VoomPv(counts = counts, 
#                         AllCov = cbind(FixCov, VarCov[names(BestCovOut$BestER)]))
#   PlotVoomPvOutER <-  PlotVoomPv(VoomPvOut = VoomPvOutER, option = option, B = B, m = m, ErrType = "ER", gam0 = gam0)
#   VoomPvOutRE <- VoomPv(counts = counts, 
#                         AllCov = cbind(FixCov, VarCov[names(BestCovOut$BestRE)]))
#   PlotVoomPvOutRE <-  PlotVoomPv(VoomPvOut = VoomPvOutRE, option = option, B = B, m = m, ErrType = "RE", gam0 = gam0)
#   res <- list(DataBSOut = DataBSOut, BBootBSOut = BBootBSOut, 
#               FSROut = FSROut, PlotFSROut = PlotFSROut, 
#               BestAlphaOut = BestAlphaOut, 
#               BestCovOut = BestCovOut, 
#               VoomPvOutER = VoomPvOutER, PlotVoomPvOutER = PlotVoomPvOutER, 
#               VoomPvOutRE = VoomPvOutRE, PlotVoomPvOutRE = PlotVoomPvOutRE, 
#               option = option, B = B, m = m)
#   pm <- proc.time()-pm
#   cat("Option=", option, ", B=", B, ", m=",m,  " runs in ", pm[3], " s\n")
#   res
# }
# 
# #Example------------
# 
# # FixCov <- covset["Line"]
# # VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# # # option <- c("WN", "OWN", "RX", "ORX")
# # B <- 40
# # m <- 3
# # amax <- 10
# # gam0 <- .05
# # ncores <- 4
# # option <- "OWN"
# # print.progress <- T
# # FSRAnalysisOut <- FSRAnalysisBS(counts[sample(nrow(counts), 3000),], FixCov, VarCov,
# #                                             option = option, B= B, m = m,
# #                                             amax= amax, gam0 = gam0, ncores=ncores,
# #                                             print.progress=print.progress)
# # FSRAnalysisOut$BestAlphaOut
# # FSRAnalysisOut$PlotVoomPvOutRE
# # FSRAnalysisOut$DataBSOut$BS
# #Function to run FSR Analysis for all 4 options
# FSRAnalysisBSAll <- function(counts, FixCov, VarCov, B, m , amax, gam0, ncores, print.progress){
#   option <- c("RX", "ORX", "WN", "OWN")
#   
#   res <- llply(option, function(i)FSRAnalysisBS(counts, FixCov, VarCov, 
#                                               option = i, B= B, m = m, 
#                                               amax= amax, gam0 = gam0, ncores=ncores, 
#                                               print.progress=print.progress))
#   names(res) <- option
#   out <- list(res = res, FixCov = FixCov, VarCov = VarCov)
#   out
# }
# #Example---------------------------
# # FixCov <- covset["Line"]
# # VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# # B <- 2
# # m <- ncol(VarCov)
# # amax <- 15
# # gam0 <- .05
# # ncores <- 2
# # print.progress <- F
# # FSRAnalysisAllOut <- FSRAnalysisAll(counts[1:10,], FixCov, VarCov,B= B, m = m,
# #                                             amax= amax, gam0 = gam0, ncores=ncores,
# #                                             print.progress=print.progress)
# 
# 
# # A function wrap all things for FSRAnalysisAll for real data ----------------
# FSRAnalysisBSAllSave <- function(counts, FixCov, VarCov, B, m, amax, gam0, ncores, print.progress, RealDataOutPath="RealDataOutBS"){
#   # RealDataOutPath <- paste0("RealDataOutBS")
#   dir.create(path = RealDataOutPath, showWarnings = F)
#   RealDataOut <- FSRAnalysisBSAll(counts, FixCov, VarCov, B, m, amax, gam0, ncores, print.progress)
#   saveRDS(RealDataOut, file = paste0(RealDataOutPath, "/",RealDataOutPath,  "_BData", B, "_mData", m, "_gam0_", gam0, ".rds"))
#   RealDataOut
# }
# 
# 
# ## Simulation using Backward Selection-----------
# 
# # SimCounts <- function(RealDataOut, nGene, option, ErrType="RE", nrep){
# #   set.seed(2017+nrep)
# #   res <- ((RealDataOut[[1]])[[option]])[[paste0("VoomPvOut",ErrType)]]
# #   EEGene <- list()
# #   AllCov <- cbind(RealDataOut$FixCov, RealDataOut$VarCov)
# #   ConCov <- AllCov[colnames(res$pvs)]
# #   lib.size <- res$lib.size; design <- res$design; Beta0 <- res$Beta; 
# #   sigma <- res$sigma;  weights <- res$weights; VarMat0 <- sigma^2*1/weights
# #   m0 <- apply(res$pvs, 2, function(x)estimate.m0(x))
# #   for(i in 1:ncol(ConCov)){
# #     if(is.factor(ConCov[,i]) | is.character(ConCov[,i])) {
# #       ct <- grep(paste0(names(ConCov)[i]),  x = colnames(design), value = T)
# #     }else{
# #       ct <- grep(paste0(names(ConCov)[i], "$"),  x = colnames(design), value = T)
# #     }
# #     EEGene[[i]] <- (sort(res$pvs[,i], decreasing = T,  index.return = T))$ix[1:m0[i]]
# #     Beta0[EEGene[[i]], ct] <- 0
# #   }
# #   names(EEGene) <- colnames(res$pvs)
# #   IndSample <- sample(1:nrow(Beta0), size = nGene)
# #   Beta <- Beta0[IndSample,]; VarMat <- VarMat0[IndSample,]
# #   sim.counts <- vapply(1:nrow(Beta), function(i){
# #     y <- design%*%Beta[i,]+ mvrnorm(n = 1, mu = rep(0, nrow(design)), Sigma = diag(VarMat[i,])) # library(MASS)
# #     unname(round(2^y*(lib.size+1)/10^6 )[ , drop = T])
# #   }, FUN.VALUE = rep(1.0, nrow(design)))
# #   SimCnt <- t(sim.counts)
# #   list(SimCnt = SimCnt, IndSample = IndSample, EEGene = EEGene)
# # }
# 
# 
# #Example------------------
# # RealDataOut <- readRDS("RFIAnalysis_Result.rds")
# # nGene <- 100
# # option <- "OWN"
# # ErrType <- "ER"
# # nrep <- 1
# # Sc <- SimCounts(RealDataOut, nGene, option, ErrType, nrep)
# 
# # Function to run 1 simulation  for 1 option--------------------
# 
# FSR1Sim <- function(RealDataOut, nGene, option, ErrType, nrep, B, m, amax, gam0, ncores, print.progress, SimPath){
#   cat("Sim = ", nrep, "\n")
#   FixCov <- RealDataOut$FixCov
#   VarCov <- RealDataOut$VarCov
#   SimCnt <- SimCounts(RealDataOut, nGene, option,ErrType,nrep)
#   # saveRDS(SimCnt, file = "SimCnt.rds")
#   SimCntOut <- FSRAnalysisAll(counts= SimCnt$SimCnt, FixCov, VarCov, B, m, amax, gam0, ncores, print.progress)
#   SimSelCov.ER <- llply(c("WN", "OWN", "RX", "ORX"),
#                         function(option)FSRSelCov(FSRAnalysisAllOut = SimCntOut, option = option, ErrType = "ER"))
#   names(SimSelCov.ER) <- c("WN.ER", "OWN.ER", "RX.ER", "ORX.ER")
#   SimSelCov.RE <- llply(c("WN", "OWN", "RX", "ORX"),
#                         function(option)FSRSelCov(FSRAnalysisAllOut = SimCntOut, option = option, ErrType = "RE"))
#   names(SimSelCov.RE) <- c("WN.RE", "OWN.RE", "RX.RE", "ORX.RE")
#   SimOut <- list(SimCnt = SimCnt, SimCntOut = SimCntOut, SimSelCov.ER = SimSelCov.ER, SimSelCov.RE = SimSelCov.RE)
#   saveRDS(SimOut, file = paste0(SimPath, "/",SimPath, ".nrep", nrep, ".rds"))
#   SimSelCov <- list(SimSelCov.ER = SimSelCov.ER, SimSelCov.RE = SimSelCov.RE)
#   SimSelCov
# }
# 
# # Function to run length(nSim) simulation from 1 option---------------------
# # Input: nSim: is vector  like 1:10 indicating the where max(nSim) is the number of simulated data
# #        RealDataOut: is the output from real data analysis from FSR method
# #        nGene: number of genes in simulated dataset
# #        option: the option to generate phony variables: "WN", "OWN", "RX", "ORX"
# #        ErrType: Type of FSR estimator formula, either "ER" or "RE"
# #        B: number of bootstrap samples during running FS for model including phony variables
# #        m: number of phony variables added
# #        amax: alpha max =15
# #        gam0: = 0.05
# #        ncores: number of cores used
# #        print.progress: printing the progress of FS
# 
# # Output: a list of 2 components. First component is the vector of True Selected Covariates from RFI data
# #         second component is a list of length(nSim) components, each of which corresponding to the set of 
# #         selected covariate for the corresponding simulated data.
# 
# FSRnSim <- function(RealDataOut, nGene, option, ErrType, nSim, B, m, amax, gam0, ncores, print.progress, SimPath){
#   TrueSelCov <- FSRSelCov(FSRAnalysisAllOut=RealDataOut, option=option, ErrType = ErrType)
#   nSimSelCov <- ldply(nSim, function(nrep){
#     SimSelCov <- FSR1Sim(RealDataOut, nGene, option, ErrType, nrep, B, m, amax, gam0, ncores, print.progress, SimPath)
#     
#     S.ER <-  laply(SimSelCov$SimSelCov.ER, function(x) length(intersect(TrueSelCov, x)))
#     R.ER <-  laply(SimSelCov$SimSelCov.ER, function(x)length(x))
#     U.ER <- R.ER - S.ER
#     # FSP <- U/pmax(R, 1)
#     FSP.ER <- U.ER/(R.ER+1) # i
#     res1 <- c(S.ER, R.ER, U.ER, FSP.ER)
#     names(res1) <- paste(rep(c("S", "R", "U", "FSP"), each = 4), rep(names(SimSelCov$SimSelCov.ER), 4), sep = ".")
#     S.RE <-  laply(SimSelCov$SimSelCov.RE, function(x) length(intersect(TrueSelCov, x)))
#     R.RE <-  laply(SimSelCov$SimSelCov.RE, function(x)length(x))
#     U.RE <- R.RE - S.RE
#     # FSP <- U/pmax(R, 1)
#     FSP.RE <- U.RE/(R.RE+1) # i
#     res2 <- c(S.RE, R.RE, U.RE, FSP.RE)
#     names(res2) <- paste(rep(c("S", "R", "U", "FSP"), each = 4), rep(names(SimSelCov$SimSelCov.RE), 4), sep = ".")
#     res <- c(res1, res2)
#     res
#   })
#   res <- list(TrueSelCov = TrueSelCov, nSimSelCov = nSimSelCov)
#   res
# }
# 
# 
# FSR1SimBS <- function(RealDataOut, nGene, option, ErrType, nrep, B, m, amax, gam0, ncores, print.progress, SimPath){
#   cat("Sim = ", nrep, "\n")
#   FixCov <- RealDataOut$FixCov
#   VarCov <- RealDataOut$VarCov
#   SimCnt <- SimCounts(RealDataOut, nGene, option,ErrType,nrep)
#   saveRDS(SimCnt, file = "SimCnt.rds")
#   SimCntOut <- FSRAnalysisBSAll(counts= SimCnt$SimCnt, FixCov, VarCov, B, m, amax, gam0, ncores, print.progress)
#   SimSelCov.ER <- llply(c("WN", "OWN", "RX", "ORX"),
#                         function(option)FSRSelCov(FSRAnalysisAllOut = SimCntOut, option = option, ErrType = "ER"))
#   names(SimSelCov.ER) <- c("WN.ER", "OWN.ER", "RX.ER", "ORX.ER")
#   SimSelCov.RE <- llply(c("WN", "OWN", "RX", "ORX"),
#                         function(option)FSRSelCov(FSRAnalysisAllOut = SimCntOut, option = option, ErrType = "RE"))
#   names(SimSelCov.RE) <- c("WN.RE", "OWN.RE", "RX.RE", "ORX.RE")
#   SimOut <- list(SimCnt = SimCnt, SimCntOut = SimCntOut, SimSelCov.ER = SimSelCov.ER, SimSelCov.RE = SimSelCov.RE)
#   saveRDS(SimOut, file = paste0(SimPath, "/",SimPath, ".nrep", nrep, ".rds"))
#   SimSelCov <- list(SimSelCov.ER = SimSelCov.ER, SimSelCov.RE = SimSelCov.RE)
#   SimSelCov
#   }
# 
# # Function to run length(nSim) simulation from 1 option---------------------
# # Input: nSim: is vector  like 1:10 indicating the where max(nSim) is the number of simulated data
# #        RealDataOut: is the output from real data analysis from FSR method
# #        nGene: number of genes in simulated dataset
# #        option: the option to generate phony variables: "WN", "OWN", "RX", "ORX"
# #        ErrType: Type of FSR estimator formula, either "ER" or "RE"
# #        B: number of bootstrap samples during running FS for model including phony variables
# #        m: number of phony variables added
# #        amax: alpha max =15
# #        gam0: = 0.05
# #        ncores: number of cores used
# #        print.progress: printing the progress of FS
# 
# # Output: a list of 2 components. First component is the vector of True Selected Covariates from RFI data
# #         second component is a list of length(nSim) components, each of which corresponding to the set of 
# #         selected covariate for the corresponding simulated data.
# 
# FSRnSimBS <- function(RealDataOut, nGene, option, ErrType, nSim, B, m, amax, gam0, ncores, print.progress, SimPath){
#   TrueSelCov <- FSRSelCov(FSRAnalysisAllOut=RealDataOut, option=option, ErrType = ErrType)
#   nSimSelCov <- ldply(nSim, function(nrep){
#     SimSelCov <- FSR1SimBS(RealDataOut, nGene, option, ErrType, nrep, B, m, amax, gam0, ncores, print.progress, SimPath)
# 
#     S.ER <-  laply(SimSelCov$SimSelCov.ER, function(x) length(intersect(TrueSelCov, x)))
#     R.ER <-  laply(SimSelCov$SimSelCov.ER, function(x)length(x))
#     U.ER <- R.ER - S.ER
#     # FSP <- U/pmax(R, 1)
#     FSP.ER <- U.ER/(R.ER+1) # i
#     res1 <- c(S.ER, R.ER, U.ER, FSP.ER)
#     names(res1) <- paste(rep(c("S", "R", "U", "FSP"), each = 4), rep(names(SimSelCov$SimSelCov.ER), 4), sep = ".")
#     S.RE <-  laply(SimSelCov$SimSelCov.RE, function(x) length(intersect(TrueSelCov, x)))
#     R.RE <-  laply(SimSelCov$SimSelCov.RE, function(x)length(x))
#     U.RE <- R.RE - S.RE
#     # FSP <- U/pmax(R, 1)
#     FSP.RE <- U.RE/(R.RE+1) # i
#     res2 <- c(S.RE, R.RE, U.RE, FSP.RE)
#     names(res2) <- paste(rep(c("S", "R", "U", "FSP"), each = 4), rep(names(SimSelCov$SimSelCov.RE), 4), sep = ".")
#     res <- c(res1, res2)
#     res
#   })
#   res <- list(TrueSelCov = TrueSelCov, nSimSelCov = nSimSelCov)
#   res
# }
# #Example------------
# # RealDataOut <- readRDS("RFIAnalysis_Result.rds")
# # nGene <- 10
# # option <- "WN"
# # ErrType <- "ER"
# # nSim <- 1:2
# # B <- 1
# # m <- 13
# # amax <- 15
# # gam0 <- 0.05
# # ncores <- 1
# # print.progress <- F
# # FSRSimOut <- FSRnSim(RealDataOut, nGene, option, ErrType,nSim,B,m,amax,gam0,ncores,print.progress)
# 
# 
# #Function wrapping up all simumation result ------
# 
# FSRnSimAllBS <- function(BData, mData, nGene, option, ErrType, nSim, BSim, mSim, amax, gam0Data, gam0, ncores, print.progress, RealDataOutPath = "RealDataOutBS"){
#   RealDataOut <- readRDS(file = paste0(RealDataOutPath, "/", RealDataOutPath,   "_BData", BData, "_mData", mData, "_gam0_", gam0Data, ".rds"))
#   if(RealDataOutPath == "RealDataOutBS"){
#     SimPath <- paste0("SimBS_", "BData", BData, "_mData", mData, "_",option, "_gamma", ErrType, "_gam0_", gam0, "_nGene", nGene, "_Bsim", BSim, "_mSim", mSim)
#   }else{
#     SimPath <- paste0("sva_SimBS_", "BData", BData, "_mData", mData, "_",option, "_gamma", ErrType, "_gam0_", gam0, "_nGene", nGene, "_Bsim", BSim, "_mSim", mSim)  
#     }
#   dir.create(path = SimPath, showWarnings = F)
#   FSRSimOut <- FSRnSimBS(RealDataOut, nGene, option, ErrType,nSim,B=BSim,m=mSim,amax,gam0,ncores,print.progress, SimPath)
#   saveRDS(FSRSimOut, file = paste0(SimPath, "/",SimPath, "_nSim", min(nSim), "_", max(nSim), ".rds"))
#   FSRSimOut
# }



### sva covariates added 03/06/2017-------------

# 
# set.seed(1)
# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# AllCov <- cbind(FixCov, VarCov)
# dm <- model.matrix(formula(paste0("~", paste0(names(AllCov), collapse = "+"))), data = AllCov)
# offset <- apply(counts, 2, quantile, .75)
# offset <- exp(log(offset) - mean(log(offset)))
# counts1 <- sweep(x = counts, MARGIN = 2, STATS = offset, FUN = "/")
# svaout <- svaseq(dat = counts1, mod = dm)
# sv  <- svaout$sv
# sv <- data.frame(sv)
# colnames(sv) <- paste0("sva", 1:ncol(sv))
# VarCov_sva <- cbind(VarCov, sv)
# saveRDS(VarCov_sva, file = "VarCov_sva.rds")
# rm(scount, FixCov, VarCov, AllCov , dm, offset, counts1, sv)

# sim <- readRDS("SimBS_BData6_mData1_OWN_gammaER_gam0_0.01_nGene2000_Bsim100_mSim1/SimBS_BData6_mData1_OWN_gammaER_gam0_0.01_nGene2000_Bsim100_mSim1.nrep1.rds")

svacov <- function(counts, FixCov, VarCov){
  set.seed(1)
  AllCov <- cbind(FixCov, VarCov)
  dm <- model.matrix(formula(paste0("~", paste0(names(AllCov), collapse = "+"))), data = AllCov)
  vout <- voom(counts = counts, lib.size = apply(counts, 2, quantile, .75), design = dm)
  n.sv <- num.sv(dat = vout$E,mod  = dm)
  if(n.sv ==0){
    sv <- data.frame(dm[,0])
  }else{
    if(n.sv < 3){
      svaout <- sva(dat = vout$E, mod = dm, n.sv = n.sv)
    }else{
      svaout <- sva(dat =  vout$E, mod = dm, n.sv = 2)
    }
    sv  <- svaout$sv
    sv <- data.frame(sv)
    colnames(sv) <- paste0("sva", 1:ncol(sv))
  }
  VarCov_sva <- cbind(VarCov, sv)
  VarCov_sva
}



##BackWard Selection Using Adding PseudoVariables--------------
bs <- function(counts, FixCov, VarCov, print.progress = F, svamethod = F){
  if(!is.data.frame(FixCov)|!is.data.frame(VarCov)) stop("FixCov and VarCov have to be data frame with different column names")
   if(svamethod)VarCov <- svacov(counts = counts, FixCov = FixCov, VarCov = VarCov)
  DelCov <- VarCov[,0] # Initiaize Set of Deleted Covariates
  ConCov <- VarCov # Initialize Set of Considered Covariates
  PvList <- list() # List of P-values for the covariate selected at each iteration
  WorstP5 <- list() # Number of p-value less than 0.05 for the covariate selected at each iteration
  if(ncol(ConCov)!=0){
    Iter <- 1
    repeat{
      if(print.progress){
        cat("----------------------------------------\n")
        cat("Iteration = ", Iter, "\n")
        cat("DelCov = ", names(DelCov), "\n")
      }
      WorstP5[[Iter]] <- Inf
      AllCov <- cbind(FixCov, ConCov)
      dm <- model.matrix(formula(paste0("~", paste0(names(AllCov), collapse = "+"))), data = AllCov)
      colnames(dm)[1] <- "Intercept"
      # if(!is.fullrank(dm)|ncol(dm)==nrow(dm)){ Ind <- T;
      # warning("Sigularity design matrix! Check generated pseudo-variables!"); break}
      vout <- voom(counts = counts, design = dm, lib.size = apply(counts, 2, quantile, .75), plot = F)
      suppressWarnings(fit <- lmFit(vout))
      for(i in 1:ncol(ConCov)){
        if(is.factor(ConCov[,i]) | is.character(ConCov[,i])) {
          ct <- paste0(grep(paste0(names(ConCov)[i]),  x = colnames(dm), value = T), collapse = ",  ")
        }else{
          ct <- paste0(grep(paste0(names(ConCov)[i], "$"),  x = colnames(dm), value = T), collapse = ",  ") 
        }
        C.matrix <- eval(parse(text=paste0("makeContrasts(",  ct, ",levels = dm)")))
        fit1 <- contrasts.fit(fit, contrasts =C.matrix)
        fit1 <- eBayes(fit1)
        tt <- topTableF(fit1, sort ="none", n = Inf)
        pv <- tt$P.Value
        P5 <- max(sum(pv<=.05, na.rm = T), 1)/max(sum(pv>.75, na.rm = T)/5, 1) # sometime pv has missing value NA make P5 cannot calculable
        if(print.progress)cat("ConCov = ", names(ConCov)[i], ", P5 = ", P5,  "\n")
        if(WorstP5[[Iter]] > P5) {
          WorstCov <- i
          WorstP5[[Iter]] <- P5
          PvList[[Iter]] <- pv
          names(PvList)[Iter] <- names(WorstP5)[Iter] <- names(ConCov)[i] 
        }
      }
      Iter <- Iter +1
      DelCov <- cbind(DelCov, ConCov[WorstCov])
      ConCov <- ConCov[-WorstCov]
      if(ncol(ConCov)==0)break
    }
    WorstP5 <- do.call("c", WorstP5)
  }else{WorstP5 <- numeric(0)}
  # P <- data.frame(t(do.call("rbind", PvList)))
  res <- list(WorstP5 = WorstP5, VarCov = VarCov, FixCov = FixCov, svamethod = svamethod)
  res
}



## Old Backward selection algorithm using voom----
##BackWard Selection Using Adding PseudoVariables--------------
jabes.bs <- function(counts, FixCov, VarCov, print.progress = F){
  if(!is.data.frame(FixCov)|!is.data.frame(VarCov)) stop("FixCov and VarCov have to be data frame with different column names")
  DelCov <- VarCov[,0] # Initiaize Set of Deleted Covariates
  ConCov <- VarCov # Initialize Set of Considered Covariates
  PvList <- list() # List of P-values for the covariate selected at each iteration
  WorstP5 <- list() # Number of p-value less than 0.05 for the covariate selected at each iteration
  Q5 <- list()
  if(ncol(ConCov)!=0){
    Iter <- 1
    repeat{
      if(print.progress){
        cat("----------------------------------------\n")
        cat("Iteration = ", Iter, "\n")
        cat("DelCov = ", names(DelCov), "\n")
      }
      WorstP5[[Iter]] <- Inf
      AllCov <- cbind(FixCov, ConCov)
      dm <- model.matrix(formula(paste0("~", paste0(names(AllCov), collapse = "+"))), data = AllCov)
      colnames(dm)[1] <- "Intercept"
      # if(!is.fullrank(dm)|ncol(dm)==nrow(dm)){ Ind <- T;
      # warning("Sigularity design matrix! Check generated pseudo-variables!"); break}
      vout <- voom(counts = counts, design = dm, lib.size = apply(counts, 2, quantile, .75), plot = F)
      suppressWarnings(fit <- lmFit(vout))
      for(i in 1:ncol(AllCov)){
        if(is.factor(AllCov[,i]) | is.character(AllCov[,i])) {
          ct <- paste0(grep(paste0(names(AllCov)[i]),  x = colnames(dm), value = T), collapse = ",  ")
        }else{
          ct <- paste0(grep(paste0(names(AllCov)[i], "$"),  x = colnames(dm), value = T), collapse = ",  ") 
        }
        C.matrix <- eval(parse(text=paste0("makeContrasts(",  ct, ",levels = dm)")))
        fit1 <- contrasts.fit(fit, contrasts =C.matrix)
        fit1 <- eBayes(fit1)
        tt <- topTableF(fit1, sort ="none", n = Inf)
        pv <- tt$P.Value
        P5 <- max(sum(pv<=.05, na.rm = T), 1) # sometime pv has missing value NA make P5 cannot calculable
        if(i ==1) Q5[[Iter]] <- sum(jabes.q(pv)<=.05)
        if(print.progress)cat("ConCov = ", names(AllCov)[i], ", P5 = ", P5,  "\n")
        if(WorstP5[[Iter]] > P5) {
          WorstCov <- i
          WorstP5[[Iter]] <- P5
          PvList[[Iter]] <- pv
          names(PvList)[Iter] <- names(WorstP5)[Iter] <- names(AllCov)[i] 
        }
      }
      Iter <- Iter +1
      DelCov <- cbind(DelCov, AllCov[WorstCov])
      ConCov <- ConCov[-(WorstCov-1)]
      if(WorstCov==1)break
    }
    WorstP5 <- do.call("c", WorstP5)
    Q5 <- do.call("c", Q5)
  }else{WorstP5 <- numeric(0)}
  # P <- data.frame(t(do.call("rbind", PvList)))
  
  BestCov <- VarCov[setdiff(names(VarCov),names(WorstP5[1:(which(Q5==max(Q5))[1]-1)]))]
  res <- list(WorstP5 = WorstP5, Q5 = Q5,  VarCov = VarCov, FixCov = FixCov, BestCov = BestCov)
  res
}

# Backward Selection for Real Data -----------
# This function running forward selection for real data, the output will be the list
# of  2 components whose first component is the covariates entering the model in the sequence together
# with its number of p-values ratio,  and the second component is a vector of number of covariates 
# entering the model when taking IV threshold alpha

# Input: counts: counts dataset
#        FixCov: a dataframe of fix covariates that are not subjected to variable selection
#        VarCov: a dataframe of flexible covariates that are subjected to variable selection
#        alpha:  a numeric vector of threshold values
#        print.progress: a logical variable either T or F to indicate that if we want to print the running process of fs function or not
# Output: is a list of 2 components: 
#        fsr: a vector of variable importance measurement in the order of the covariate entering the model in fs process. 
#             length of fsr equals number of column of VarCov
#        S: a vector of the number of selected variates at different threshold alpha. 
DataBS <- function(counts, FixCov, VarCov, alpha, print.progress = F, svamethod = F){
  # if(svamethod)VarCov <- svacov(counts = counts, FixCov = FixCov, VarCov = VarCov)
  BSout <- bs(counts = counts, FixCov = FixCov, VarCov = VarCov, print.progress = print.progress, svamethod = svamethod)
  BS <- BSout$WorstP5
  S <- vapply(1:length(alpha), function(i)sum(cummax(BS) >alpha[i]), FUN.VALUE = 1.0)
  return(list(BS = BS, S = S, alpha = alpha, FixCov = BSout$FixCov, VarCov = BSout$VarCov, svamethod = BSout$svamethod))
}

#Example--------
# pm <- proc.time()
# amax <-  10
# alpha <- c(seq(1, amax, length = 800))
# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# DataBSOut <- DataBS(counts = counts, FixCov = FixCov, VarCov = VarCov, alpha, print.progress = T)
# proc.time() - pm
# out$S



### Bootstrap function to generate data and calculate U, SP for the simulatied data-------
## Input: nrep: seed of the bootstrap
##        option: either white noise (WN), permutation of row of X (RX), or orthogonal version (OWN, ORX)
##        FixCov: Fix covariates that are not subjected to variable selection
##        VarCov:  Flexible covariates that are subjected to variable selection
##        m: number of phony variables need to generated
##        counts: counts data with G genes and s columns
##        alpha: a vector of threshold of VI measurement
# cor(covset[c("RFI", "RINa", "RINb", "Conca", "Concb", "lymp", 
#              "neut", "mono", "baso", "eosi")])
BootBS <- function(counts, FixCov, VarCov, m=3, nrep, alpha, option = "OWN", print.progress = T){
  set.seed(nrep)
  # cat("nrep = ", nrep, "\n")
  # 
  # # if(svamethod)VarCov <- svacov(counts = counts, FixCov = FixCov, VarCov = VarCov)
  
  AllCov0 <- cbind(FixCov, VarCov)
  dm0 <- model.matrix(formula(paste0("~", paste0(names(AllCov0), collapse = "+"))), data = AllCov0)
  colnames(dm0)[1] <- "Intercept"
  repeat{
    Nois <- PhoVar(dm0 = dm0, m = m, option = option)
    VarCov2 <- cbind(VarCov, Nois)
    AllCov1 <- cbind(FixCov, VarCov2)
    dm1 <- model.matrix(formula(paste0("~", paste0(names(AllCov1), collapse = "+"))), data = AllCov1)
    if(is.fullrank(dm1)) break
  }
  bsout <- bs(counts = counts, FixCov = FixCov, VarCov = VarCov2, print.progress = print.progress, svamethod = F)
  pv <- bsout$WorstP5
  pvmax <- cummax(pv)
  vor <- names(pv)
  U <- vapply(1:length(alpha), function(i)sum(vor%in%names(Nois) & pvmax>alpha[i]), FUN.VALUE = 1.0)
  SP <- vapply(1:length(alpha), function(i)sum( pvmax>alpha[i]), FUN.VALUE = 1.0)
  res <- c(U = U, SP = SP)
  # cat("nrep = ", nrep, "\n")
  res
}


#Example---------------------
# amax <-  10
# alpha <- c(seq(1, amax, length = 800))
# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# option <- "ORX"
# out <- BootFS(counts = counts[1:50,], FixCov = FixCov, VarCov = VarCov, m = 20, nrep = 1, alpha = alpha ,option = option, print.progress = F)

###Output of Running B bootstrap samples------------------------
# This function will obtain number of phony variables selected (U) and total number of variables selected (SP) when 
# considering m phony variable for each threshold level alpha. The results will be a matrix of B x (2*length(alpha)) value
# Input: The same as input for bootstrap function above and
#        B: the number of bootstrap sample
#        ncores: number of cores for embrassingly parallel running

BBootBS <- function(B, ncores, counts, FixCov, VarCov, m, alpha, option = "OWN", print.progress = F){
  USP <- mclapply(1:B, function(nrep)BootBS(counts, FixCov, 
                                            VarCov, m=m, nrep, 
                                            alpha, option = option, 
                                            print.progress), 
                  mc.cores = ncores)
  USP <- do.call("rbind", USP)
  list(USP=USP, mPhoCov = m)
}
#Example--------------------------------
# pm <- proc.time()
# amax <-  10
# alpha <- c(seq(1, amax, length = 800))
# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# option <- "OWN"
# B <- 40
# ncores <- 4
# m <- 3
# BBootBSOut <- BBootBS(B,ncores = ncores,  counts = counts, FixCov = FixCov, VarCov = VarCov, m = 3, alpha = alpha ,option = option, print.progress = T)
# proc.time() - pm
# 
# 

# DataFSOut  <- res$res$WN$DataFSOut
# #  BBootFSOut <- res$res$WN$BBootFSOut
FSRPBS <- function(DataBSOut, BBootBSOut){
  S <- DataBSOut$S
  alpha <- DataBSOut$alpha
  nAlpha <- length(alpha)
  nVarCov <- length(DataBSOut$BS)
  nPhoCov <- BBootBSOut$mPhoCov
  B <- nrow(BBootBSOut$USP)
  USP <- apply(BBootBSOut$USP, 2, sum)/B
  U <- USP[1:nAlpha]
  SP <- USP[-c(1:nAlpha)]
  ghat.P.RE <- U/(1+SP)
  ghat.P.ER <- U/(1+S)
  out <- list(B = B, nPhoCov = nPhoCov, BS = DataBSOut$BS, alpha = alpha, ghat.P.RE = ghat.P.RE, ghat.P.ER = ghat.P.ER)
  out
}
# 
# # estimate.kU <- function(FSRPOut){
# #   
# # }
# 
# 
# FSRPBSOut <- FSRPBS(DataBSOut,BBootBSOut)

# ## Optimal threshold alpha ---------------
BestAlphaBS <- function(FSRPBSOut, gam0=0.05){
  alpha <- FSRPBSOut$alpha
  BS <- FSRPBSOut$BS
  kT <- length(BS)
  kP <- FSRPBSOut$nPhoCov
  ghat.P.RE <- FSRPBSOut$ghat.P.RE
  alphamin.RE<-min(alpha[ghat.P.RE==max(ghat.P.RE)])  # alpha with largest ghat.RE
  kU.RE <- kT
  gam0.P <- kP*gam0/(kU.RE +kP*gam0)
  repeat{
    ind.RE <- (cummin(ghat.P.RE)<=gam0.P& alpha>=alphamin.RE)*1
    alphahat.RE <- alpha[(which(ind.RE==1))[1]]       # RE est. of alpha
    if(is.na(alphahat.RE)){
      warning(paste0("All Estimated FSRs are greater than ", gam0, ". This happens in BestAlpha function."))
      alphahat.RE <- 10000
    }
    temp <- kT - length(BS[cummax(BS)>=alphahat.RE])
    if(temp==kU.RE) {break}
    kU.RE <- temp
    gam0.P <- kP*gam0/(kU.RE +kP*gam0)
  }
  
  ghat.P.ER <- FSRPBSOut$ghat.P.ER
  alphamin.ER<-min(alpha[ghat.P.ER==max(ghat.P.ER)])  # alpha with largest ghat.RE
  kU.ER <- kT
  gam0.P <- kP*gam0/kU.ER
  repeat{
    ind.ER <- (cummin(ghat.P.ER)<=gam0.P& alpha>=alphamin.ER)*1
    alphahat.ER <- alpha[(which(ind.ER==1))[1]]       # ER est. of alpha
    if(is.na(alphahat.ER)){
      warning(paste0("All Estimated FSRs are greater than ", gam0, ". This happens in BestAlpha function."))
      alphahat.ER <- 10000
    }
    temp <- kT - length(BS[cummax(BS)>=alphahat.ER])
    if(temp==kU.ER) {break}
    kU.ER <- temp
    gam0.P <- kP*gam0/(kU.ER)
  }
  
  c(alphahat.RE  = alphahat.RE , alphahat.ER  = alphahat.ER, kU.RE = kU.RE, kU.ER = kU.ER )
}

# BestAlphaOut <- BestAlphaBS(FSRPBSOut)


## Input will be alpha, output from bootstrap
FSR <- function(FSRPOut, BestAlphaOut, gam0 = .05){
  # FSRPOut <- FSRP(DataFSOut, BBootFSOut)
  # BestAlphaOut <- BestAlpha(FSRPOut, gam0 = gam0)
  ghat.RE<-pmin(1, FSRPOut$ghat.P.RE/(1-FSRPOut$ghat.P.RE)*BestAlphaOut["kU.RE"]/FSRPOut$nPhoCov)      # gammahat_RE, SP-U plays role of S
  ghat<- pmin(1, FSRPOut$ghat.P.ER*BestAlphaOut["kU.ER"]/FSRPOut$nPhoCov)             # gammahat_ER
  out <- cbind(alpha = FSRPOut$alpha, ghat.RE = ghat.RE, ghat = ghat)
  rownames(out) <- NULL
  data.frame(out)
}


#### Function calculate FSR both RE and ER ------------
## Input will be alpha, output from bootstrap
# FSR <- function(DataFSOut, BBootFSOut){
#   S <- DataFSOut$S                
#   alpha <- DataFSOut$alpha
#   nAlpha <- length(alpha)
#   nVarCov <- length(DataFSOut$FS)
#   nPhoCov <- BBootFSOut$mPhoCov
#   B <- nrow(BBootFSOut$USP)
#   USP <- apply(BBootFSOut$USP, 2, sum)/B
#   U <- USP[1:nAlpha]
#   SP <- USP[-c(1:nAlpha)]
#   ghat.RE<-(nVarCov - S)*(U/nPhoCov)/(1+SP-U)        # gammahat_RE, SP-U plays role of S
#   ghat<-(nVarCov - S)*(U/nPhoCov)/(1+S)              # gammahat_ER
#   out <- cbind(alpha = alpha, ghat.RE = ghat.RE, ghat = ghat)
#   rownames(out) <- NULL
#   data.frame(out)
# }
# 

#Example----------------
# amax <-  15
# alpha <- c(seq(1, amax, length = 800))
# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# DataOut <- DataFS(counts = counts[1:10,], FixCov = FixCov, VarCov = VarCov, alpha, print.progress = T)
# option <- "ORX"
# B <- 2
# ncores <- 1
# BootOut <- BBoot(B,ncores = ncores,  counts = counts[1:10,], FixCov = FixCov, VarCov = VarCov, m = ncol(VarCov), alpha = alpha ,option = option, print.progress = T)
# FSROut <- FSR(DataOut, BootOut)

# Plot 2 FSR Estimates----------------
# PlotFSR <- function(FSROut, gam0, FileName){
#   pdf(paste0(FileName, ".pdf"))
#   plot(FSROut$alpha, FSROut$ghat, type="l",lty=2,ylab="Est.FSR",xlab="alpha")
#   lines(FSROut$alpha,FSROut$ghat.RE,col="red")
#   abline(h = gam0)
#   title("Est.FSR.RE (red) and FSR.ER (dash) Curves")
#   dev.off()
# }

# FSROut <- FSR(FSRPOut = FSRPBSOut, BestAlphaOut)

# PlotFSR(FSROut = FSROut, gam0 = .05, option = option, B = B, m=m)
# FSROut <- FSR(DataFSOut,BBootFSOut)
# PlotFSR(FSROut, gam0 = .05, option = "OWN", B = 240, m = 20)
#Example----------
#PlotFSR(FSROut = FSROut, gam0= .05, FileName = "Test")

## Optimal threshold alpha ---------------
# BestAlpha <- function(FSROut, gam0){
#   alpha <- FSROut$alpha
#   ghat.RE <- FSROut$ghat.RE
#   alphamin.RE<-min(alpha[ghat.RE==max(ghat.RE)])  # alpha with largest ghat.RE
#   ind.RE <- (cummin(ghat.RE)<=gam0& alpha>=alphamin.RE)*1
#   alphahat.RE <- alpha[(which(ind.RE==1))[1]]       # RE est. of alpha
#   
#   ghat <- FSROut$ghat
#   
#   alphamin<-min(alpha[ghat==max(ghat)])     # alpha with largest ghat
#   ind <- (cummin(ghat)<=gam0 & alpha>=alphamin)*1
#   
#   
#   alphahat.ER <- alpha[(which(ind==1))[1]]          # ER est. of alpha
#   if(is.na(alphahat.RE)){
#     warning(paste0("All Estimated FSRs are greater than ", gam0, ". This happens in BestAlpha function."))
#     alphahat.RE <- 10000
#     }
#   if(is.na(alphahat.ER)) {
#     warning(paste0("All Estimated FSRs are greater than ", gam0, ". This happens in BestAlpha function."))
#     alphahat.ER <- 10000
#     }
#   c(alphahat.RE  = alphahat.RE , alphahat.ER  = alphahat.ER )
# }

#Example ---------
# BestAlphaOut <- BestAlpha(FSROut, gam0 = .05)

# Function to calculate the selected covariates using FSR method
BestCovBS <- function(DataBSOut, BestAlphaOut){
  list(BestER = DataBSOut$BS[cummax(DataBSOut$BS)>=BestAlphaOut["alphahat.ER"]],
       BestRE = DataBSOut$BS[cummax(DataBSOut$BS)>=BestAlphaOut["alphahat.RE"]])
}

# BestCovBSOut <- BestCovBS(DataBSOut, BestAlphaOut)
# 
# BestCovOut

#Example--------
# PlotVoomPv(VoomPvOut, option, B, m, ErrType)


# The function shows FSR algorithm for data analysis----------
FSRAnalysisBS <- function(counts, FixCov, VarCov, option = "OWN", B= 100, m = 3, amax= 10, gam0 = .05, ncores, print.progress, svamethod= F){
  cat("Option=", option, ", B=", B, ", m=",m, ", gam0 = ", gam0, "\n")
  pm <- proc.time()
  alpha <- c(seq(1, amax, length = 800))
  DataBSOut <- DataBS(counts = counts, FixCov = FixCov, VarCov=VarCov, alpha = alpha, print.progress = print.progress, svamethod = svamethod)
  BBootBSOut <- BBootBS(B= B, ncores = ncores, counts = counts, FixCov, VarCov = DataBSOut$VarCov, m= m, alpha = alpha, option = option, print.progress = print.progress)
  # saveRDS(BBootBSOut, file = "BBootBSOut.rds")
  # saveRDS(DataBSOut, file = "DataBSOut.rds")
  FSRPBSOut <- FSRPBS(DataBSOut=DataBSOut, BBootBSOut = BBootBSOut)
  BestAlphaOut <- BestAlphaBS(FSRPBSOut, gam0)
  FSROut <- FSR(FSRPBSOut, BestAlphaOut)
  PlotFSROut <- PlotFSR(FSROut = FSROut, gam0 = gam0, option = option, B = B, m = m)
  BestCovOut <- BestCovBS(DataBSOut = DataBSOut, BestAlphaOut = BestAlphaOut)
  VoomPvOutER <- VoomPv(counts = counts, 
                        AllCov = cbind(FixCov, (DataBSOut$VarCov)[names(BestCovOut$BestER)]))
  PlotVoomPvOutER <-  PlotVoomPv(VoomPvOut = VoomPvOutER, option = option, B = B, m = m, ErrType = "ER", gam0 = gam0)
  VoomPvOutRE <- VoomPv(counts = counts, 
                        AllCov = cbind(FixCov, (DataBSOut$VarCov)[names(BestCovOut$BestRE)]))
  PlotVoomPvOutRE <-  PlotVoomPv(VoomPvOut = VoomPvOutRE, option = option, B = B, m = m, ErrType = "RE", gam0 = gam0)
  res <- list(DataBSOut = DataBSOut, BBootBSOut = BBootBSOut, 
              FSROut = FSROut, PlotFSROut = PlotFSROut, 
              BestAlphaOut = BestAlphaOut, 
              BestCovOut = BestCovOut, 
              VoomPvOutER = VoomPvOutER, PlotVoomPvOutER = PlotVoomPvOutER, 
              VoomPvOutRE = VoomPvOutRE, PlotVoomPvOutRE = PlotVoomPvOutRE, 
              option = option, B = B, m = m)
  pm <- proc.time()-pm
  cat("Option=", option, ", B=", B, ", m=",m,  " runs in ", pm[3], " s\n")
  res
}

#Example------------

# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# # option <- c("WN", "OWN", "RX", "ORX")
# B <- 40
# m <- 3
# amax <- 10
# gam0 <- .05
# ncores <- 4
# option <- "OWN"
# print.progress <- T
# FSRAnalysisOut <- FSRAnalysisBS(counts[sample(nrow(counts), 3000),], FixCov, VarCov,
#                                             option = option, B= B, m = m,
#                                             amax= amax, gam0 = gam0, ncores=ncores,
#                                             print.progress=print.progress)
# FSRAnalysisOut$BestAlphaOut
# FSRAnalysisOut$PlotVoomPvOutRE
# FSRAnalysisOut$DataBSOut$BS
#Function to run FSR Analysis for all 4 options
FSRAnalysisBSAll <- function(counts, FixCov, VarCov, B, m , amax, gam0, ncores, print.progress, svamethod = F){
  # option <- c("RX", "ORX", "WN", "OWN")
  option <- c( "ORX", "OWN")
  
  res <- llply(option, function(i)FSRAnalysisBS(counts, FixCov, VarCov, 
                                                option = i, B= B, m = m, 
                                                amax= amax, gam0 = gam0, ncores=ncores, 
                                                print.progress=print.progress, svamethod = svamethod))
  names(res) <- option
  out <- list(res = res, FixCov = FixCov, VarCov = VarCov)
  out
}
#Example---------------------------
# FixCov <- covset["Line"]
# VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
# B <- 2
# m <- ncol(VarCov)
# amax <- 15
# gam0 <- .05
# ncores <- 2
# print.progress <- F
# FSRAnalysisAllOut <- FSRAnalysisAll(counts[1:10,], FixCov, VarCov,B= B, m = m,
#                                             amax= amax, gam0 = gam0, ncores=ncores,
#                                             print.progress=print.progress)


# A function wrap all things for FSRAnalysisAll for real data ----------------
FSRAnalysisBSAllSave <- function(counts, FixCov, VarCov, B, m, amax, gam0, ncores, print.progress, RealDataOutPath="RealDataOutBS"){
  # RealDataOutPath <- paste0("RealDataOutBS")
  dir.create(path = RealDataOutPath, showWarnings = F)
  svamethod <- F
  if(RealDataOutPath == "RealDataOutBS_sva")svamethod <- T
  RealDataOut <- FSRAnalysisBSAll(counts, FixCov, VarCov, B, m, amax, gam0, ncores, print.progress, svamethod = svamethod)
  saveRDS(RealDataOut, file = paste0(RealDataOutPath, "/",RealDataOutPath,  "_BData", B, "_mData", m, "_gam0_", gam0, ".rds"))
  RealDataOut
}


FSRAnalysisBSAllSave2 <- function(counts, FixCov, VarCov, B, m, amax, gam0, ncores, print.progress, RealDataOutPath="RealDataOutBS", svamethod = F){
  # RealDataOutPath <- paste0("RealDataOutBS")
  # dir.create(path = RealDataOutPath, showWarnings = F)
  RealDataOut <- FSRAnalysisBSAll(counts, FixCov, VarCov, B, m, amax, gam0, ncores, print.progress, svamethod = svamethod)
  saveRDS(RealDataOut, file = paste0("../output/",RealDataOutPath,  "_BData", B, "_mData", m, "_gam0_", gam0, "_sva", svamethod, ".rds"))
  RealDataOut
}

## Simulation using Backward Selection-----------

# SimCounts <- function(RealDataOut, nGene, option, ErrType="RE", nrep){
#   set.seed(2017+nrep)
#   res <- ((RealDataOut[[1]])[[option]])[[paste0("VoomPvOut",ErrType)]]
#   EEGene <- list()
#   AllCov <- cbind(RealDataOut$FixCov, RealDataOut$VarCov)
#   ConCov <- AllCov[colnames(res$pvs)]
#   lib.size <- res$lib.size; design <- res$design; Beta0 <- res$Beta; 
#   sigma <- res$sigma;  weights <- res$weights; VarMat0 <- sigma^2*1/weights
#   m0 <- apply(res$pvs, 2, function(x)estimate.m0(x))
#   for(i in 1:ncol(ConCov)){
#     if(is.factor(ConCov[,i]) | is.character(ConCov[,i])) {
#       ct <- grep(paste0(names(ConCov)[i]),  x = colnames(design), value = T)
#     }else{
#       ct <- grep(paste0(names(ConCov)[i], "$"),  x = colnames(design), value = T)
#     }
#     EEGene[[i]] <- (sort(res$pvs[,i], decreasing = T,  index.return = T))$ix[1:m0[i]]
#     Beta0[EEGene[[i]], ct] <- 0
#   }
#   names(EEGene) <- colnames(res$pvs)
#   IndSample <- sample(1:nrow(Beta0), size = nGene)
#   Beta <- Beta0[IndSample,]; VarMat <- VarMat0[IndSample,]
#   sim.counts <- vapply(1:nrow(Beta), function(i){
#     y <- design%*%Beta[i,]+ mvrnorm(n = 1, mu = rep(0, nrow(design)), Sigma = diag(VarMat[i,])) # library(MASS)
#     unname(round(2^y*(lib.size+1)/10^6 )[ , drop = T])
#   }, FUN.VALUE = rep(1.0, nrow(design)))
#   SimCnt <- t(sim.counts)
#   list(SimCnt = SimCnt, IndSample = IndSample, EEGene = EEGene)
# }


#Example------------------
# RealDataOut <- readRDS("RFIAnalysis_Result.rds")
# nGene <- 100
# option <- "OWN"
# ErrType <- "ER"
# nrep <- 1
# Sc <- SimCounts(RealDataOut, nGene, option, ErrType, nrep)

# Function to run 1 simulation  for 1 option--------------------

FSR1Sim <- function(RealDataOut, nGene, option, ErrType, nrep, B, m, amax, gam0, ncores, print.progress, SimPath, svamethod =F){
  cat("Sim = ", nrep, "\n")
  FixCov <- RealDataOut$FixCov
  VarCov <- RealDataOut$VarCov
  SimCnt <- SimCounts(RealDataOut, nGene, option,ErrType,nrep)
  # saveRDS(SimCnt, file = "SimCnt.rds")
  SimCntOut <- FSRAnalysisAll(counts= SimCnt$SimCnt, FixCov, VarCov, B, m, amax, gam0, ncores, print.progress, svamethod = svamethod)
  SimSelCov.ER <- llply(c("WN", "OWN", "RX", "ORX"),
                        function(option)FSRSelCov(FSRAnalysisAllOut = SimCntOut, option = option, ErrType = "ER"))
  names(SimSelCov.ER) <- c("WN.ER", "OWN.ER", "RX.ER", "ORX.ER")
  SimSelCov.RE <- llply(c("WN", "OWN", "RX", "ORX"),
                        function(option)FSRSelCov(FSRAnalysisAllOut = SimCntOut, option = option, ErrType = "RE"))
  names(SimSelCov.RE) <- c("WN.RE", "OWN.RE", "RX.RE", "ORX.RE")
  SimOut <- list(SimCnt = SimCnt, SimCntOut = SimCntOut, SimSelCov.ER = SimSelCov.ER, SimSelCov.RE = SimSelCov.RE)
  saveRDS(SimOut, file = paste0(SimPath, "/",SimPath, ".nrep", nrep, ".rds"))
  SimSelCov <- list(SimSelCov.ER = SimSelCov.ER, SimSelCov.RE = SimSelCov.RE)
  SimSelCov
}

# Function to run length(nSim) simulation from 1 option---------------------
# Input: nSim: is vector  like 1:10 indicating the where max(nSim) is the number of simulated data
#        RealDataOut: is the output from real data analysis from FSR method
#        nGene: number of genes in simulated dataset
#        option: the option to generate phony variables: "WN", "OWN", "RX", "ORX"
#        ErrType: Type of FSR estimator formula, either "ER" or "RE"
#        B: number of bootstrap samples during running FS for model including phony variables
#        m: number of phony variables added
#        amax: alpha max =15
#        gam0: = 0.05
#        ncores: number of cores used
#        print.progress: printing the progress of FS

# Output: a list of 2 components. First component is the vector of True Selected Covariates from RFI data
#         second component is a list of length(nSim) components, each of which corresponding to the set of 
#         selected covariate for the corresponding simulated data.

FSRnSim <- function(RealDataOut, nGene, option, ErrType, nSim, B, m, amax, gam0, ncores, print.progress, SimPath, svamethod = F){
  TrueSelCov <- FSRSelCov(FSRAnalysisAllOut=RealDataOut, option=option, ErrType = ErrType)
  nSimSelCov <- ldply(nSim, function(nrep){
    SimSelCov <- FSR1Sim(RealDataOut, nGene, option, ErrType, nrep, B, m, amax, gam0, ncores, print.progress, SimPath, svamethod = svamethod)
    
    S.ER <-  laply(SimSelCov$SimSelCov.ER, function(x) length(intersect(TrueSelCov, x)))
    R.ER <-  laply(SimSelCov$SimSelCov.ER, function(x)length(x))
    U.ER <- R.ER - S.ER
    # FSP <- U/pmax(R, 1)
    FSP.ER <- U.ER/(R.ER+1) # i
    res1 <- c(S.ER, R.ER, U.ER, FSP.ER)
    names(res1) <- paste(rep(c("S", "R", "U", "FSP"), each = 4), rep(names(SimSelCov$SimSelCov.ER), 4), sep = ".")
    S.RE <-  laply(SimSelCov$SimSelCov.RE, function(x) length(intersect(TrueSelCov, x)))
    R.RE <-  laply(SimSelCov$SimSelCov.RE, function(x)length(x))
    U.RE <- R.RE - S.RE
    # FSP <- U/pmax(R, 1)
    FSP.RE <- U.RE/(R.RE+1) # i
    res2 <- c(S.RE, R.RE, U.RE, FSP.RE)
    names(res2) <- paste(rep(c("S", "R", "U", "FSP"), each = 4), rep(names(SimSelCov$SimSelCov.RE), 4), sep = ".")
    res <- c(res1, res2)
    res
  })
  res <- list(TrueSelCov = TrueSelCov, nSimSelCov = nSimSelCov)
  res
}


FSR1SimBS <- function(RealDataOut, nGene, option, ErrType, nrep, B, m, amax, gam0, ncores, print.progress, SimPath, svamethod = F){
  # cat("Sim = ", nrep, "\n")
  FixCov <- RealDataOut$FixCov
  VarCov <- RealDataOut$VarCov
  svaCov <- grep(pattern = "sva", x = names(VarCov), value = T)
  VarCov <- VarCov[setdiff(names(VarCov), svaCov)]
  SimCnt <- SimCounts(RealDataOut, nGene, option,ErrType,nrep)
  # saveRDS(SimCnt, file = "SimCnt.rds")
  SimCntOut <- FSRAnalysisBSAll(counts= SimCnt$SimCnt, FixCov, VarCov, B, m, amax, gam0, ncores, print.progress, svamethod = svamethod)
  SimSelCov.ER <- llply(c("WN", "OWN", "RX", "ORX"),
                        function(option)FSRSelCov(FSRAnalysisAllOut = SimCntOut, option = option, ErrType = "ER"))
  names(SimSelCov.ER) <- c("WN.ER", "OWN.ER", "RX.ER", "ORX.ER")
  SimSelCov.RE <- llply(c("WN", "OWN", "RX", "ORX"),
                        function(option)FSRSelCov(FSRAnalysisAllOut = SimCntOut, option = option, ErrType = "RE"))
  names(SimSelCov.RE) <- c("WN.RE", "OWN.RE", "RX.RE", "ORX.RE")
  SimOut <- list(SimCnt = SimCnt, SimCntOut = SimCntOut, SimSelCov.ER = SimSelCov.ER, SimSelCov.RE = SimSelCov.RE)
  saveRDS(SimOut, file = paste0(SimPath, "/",SimPath, ".nrep", nrep, ".rds"))
  SimSelCov <- list(SimSelCov.ER = SimSelCov.ER, SimSelCov.RE = SimSelCov.RE)
  SimSelCov
}

# Function to run length(nSim) simulation from 1 option---------------------
# Input: nSim: is vector  like 1:10 indicating the where max(nSim) is the number of simulated data
#        RealDataOut: is the output from real data analysis from FSR method
#        nGene: number of genes in simulated dataset
#        option: the option to generate phony variables: "WN", "OWN", "RX", "ORX"
#        ErrType: Type of FSR estimator formula, either "ER" or "RE"
#        B: number of bootstrap samples during running FS for model including phony variables
#        m: number of phony variables added
#        amax: alpha max =15
#        gam0: = 0.05
#        ncores: number of cores used
#        print.progress: printing the progress of FS

# Output: a list of 2 components. First component is the vector of True Selected Covariates from RFI data
#         second component is a list of length(nSim) components, each of which corresponding to the set of 
#         selected covariate for the corresponding simulated data.

FSRnSimBS <- function(RealDataOut, nGene, option, ErrType, nSim, B, m, amax, gam0, ncores, print.progress, SimPath, svamethod = F){
  TrueSelCov <- FSRSelCov(FSRAnalysisAllOut=RealDataOut, option=option, ErrType = ErrType)
  nSimSelCov <- ldply(nSim, function(nrep){
    SimSelCov <- FSR1SimBS(RealDataOut, nGene, option, ErrType, nrep, B, m, amax, gam0, ncores, print.progress, SimPath, svamethod = svamethod)
    
    S.ER <-  laply(SimSelCov$SimSelCov.ER, function(x) length(intersect(TrueSelCov, x)))
    R.ER <-  laply(SimSelCov$SimSelCov.ER, function(x)length(x))
    U.ER <- R.ER - S.ER
    # FSP <- U/pmax(R, 1)
    FSP.ER <- U.ER/(R.ER+1) # i
    res1 <- c(S.ER, R.ER, U.ER, FSP.ER)
    names(res1) <- paste(rep(c("S", "R", "U", "FSP"), each = 4), rep(names(SimSelCov$SimSelCov.ER), 4), sep = ".")
    S.RE <-  laply(SimSelCov$SimSelCov.RE, function(x) length(intersect(TrueSelCov, x)))
    R.RE <-  laply(SimSelCov$SimSelCov.RE, function(x)length(x))
    U.RE <- R.RE - S.RE
    # FSP <- U/pmax(R, 1)
    FSP.RE <- U.RE/(R.RE+1) # i
    res2 <- c(S.RE, R.RE, U.RE, FSP.RE)
    names(res2) <- paste(rep(c("S", "R", "U", "FSP"), each = 4), rep(names(SimSelCov$SimSelCov.RE), 4), sep = ".")
    res <- c(res1, res2)
    res
  })
  res <- list(TrueSelCov = TrueSelCov, nSimSelCov = nSimSelCov)
  res
}
#Example------------
# RealDataOut <- readRDS("RFIAnalysis_Result.rds")
# nGene <- 10
# option <- "WN"
# ErrType <- "ER"
# nSim <- 1:2
# B <- 1
# m <- 13
# amax <- 15
# gam0 <- 0.05
# ncores <- 1
# print.progress <- F
# FSRSimOut <- FSRnSim(RealDataOut, nGene, option, ErrType,nSim,B,m,amax,gam0,ncores,print.progress)


#Function wrapping up all simumation result ------

FSRnSimAllBS <- function(BData, mData, nGene, option, ErrType, nSim, BSim, mSim, amax, gam0Data, gam0, ncores, print.progress, RealDataOutPath = "RealDataOutBS", svamethod = F){
  RealDataOut <- readRDS(file = paste0(RealDataOutPath, "/", RealDataOutPath,   "_BData", BData, "_mData", mData, "_gam0_", gam0Data, ".rds"))
  if(svamethod == F){
    SimPath <- paste0("SimBS_", "BData", BData, "_mData", mData, "_",option, "_gamma", ErrType, "_gam0_", gam0, "_nGene", nGene, "_Bsim", BSim, "_mSim", mSim)
  }else{
    SimPath <- paste0("sva_SimBS_", "BData", BData, "_mData", mData, "_",option, "_gamma", ErrType, "_gam0_", gam0, "_nGene", nGene, "_Bsim", BSim, "_mSim", mSim)  
  }
  dir.create(path = SimPath, showWarnings = F)
  FSRSimOut <- FSRnSimBS(RealDataOut, nGene, option, ErrType,nSim,B=BSim,m=mSim,amax,gam0,ncores,print.progress, SimPath, svamethod = svamethod)
  saveRDS(FSRSimOut, file = paste0(SimPath, "/",SimPath, "_nSim", min(nSim), "_", max(nSim), ".rds"))
  FSRSimOut
}



## Ignore all covariates except main factor of interest, using sva covariates instead----------

# This function create sva covariates from sva  or dSVA method
# input: counts: counts data 
# FixCov: main factors of interest, not subject to variable selection
# VarCov: adjusted covariates included to get the full model, then to calculate sva covariates
# method: either "sva" or "dSVA" 

cov_sva_only <- function(counts, FixCov, VarCov){
  set.seed(1)
  AllCov <- cbind(FixCov, VarCov)
  dm <- model.matrix(formula(paste0("~", paste0(names(AllCov), collapse = "+"))), data = AllCov)
  vout <- voom(counts = counts, lib.size = apply(counts, 2, quantile, .75), design = dm)
  sink("aux"); 
  n.sv <- num.sv(dat = vout$E,mod  = dm)
  sink(NULL); 
  if(n.sv ==0){
    sv1 <- sv2 <- data.frame(dm[,0])
  }else{
      sink("aux"); 
      svaout1 <- sva(dat =  vout$E, mod = dm, n.sv = n.sv)
      sink(NULL); 
      sv1  <- svaout1$sv
      sv1 <- data.frame(sv1)
      colnames(sv1) <- paste0("sva", 1:ncol(sv1))
      sink("aux"); 
      svaout2 <- dSVA(Y =  t(vout$E), X =dm[,-1], ncomp =  n.sv)
      sink(NULL); 
      sv2  <- svaout2$Z
      sv2 <- data.frame(sv2)
      colnames(sv2) <- paste0("sva", 1:ncol(sv2))
      
  }
  svacov <- list(sva = sv1, dSVA = sv2)
  svacov
}



##BackWard Selection Using Adding PseudoVariables--------------
bs_sva_only <- function(counts, FixCov, VarCov, print.progress = F, svamethod = F){
  if(!is.data.frame(FixCov)|!is.data.frame(VarCov)) stop("FixCov and VarCov have to be data frame with different column names")
  if(svamethod)VarCov <- cov_sva_only(counts = counts, FixCov = FixCov, VarCov = VarCov)
  DelCov <- VarCov[,0] # Initiaize Set of Deleted Covariates
  ConCov <- VarCov # Initialize Set of Considered Covariates
  PvList <- list() # List of P-values for the covariate selected at each iteration
  WorstP5 <- list() # Number of p-value less than 0.05 for the covariate selected at each iteration
  if(ncol(ConCov)!=0){
    Iter <- 1
    repeat{
      if(print.progress){
        cat("----------------------------------------\n")
        cat("Iteration = ", Iter, "\n")
        cat("DelCov = ", names(DelCov), "\n")
      }
      WorstP5[[Iter]] <- Inf
      AllCov <- cbind(FixCov, ConCov)
      dm <- model.matrix(formula(paste0("~", paste0(names(AllCov), collapse = "+"))), data = AllCov)
      colnames(dm)[1] <- "Intercept"
      # if(!is.fullrank(dm)|ncol(dm)==nrow(dm)){ Ind <- T;
      # warning("Sigularity design matrix! Check generated pseudo-variables!"); break}
      vout <- voom(counts = counts, design = dm, lib.size = apply(counts, 2, quantile, .75), plot = F)
      suppressWarnings(fit <- lmFit(vout))
      for(i in 1:ncol(ConCov)){
        if(is.factor(ConCov[,i]) | is.character(ConCov[,i])) {
          ct <- paste0(grep(paste0(names(ConCov)[i]),  x = colnames(dm), value = T), collapse = ",  ")
        }else{
          ct <- paste0(grep(paste0(names(ConCov)[i], "$"),  x = colnames(dm), value = T), collapse = ",  ") 
        }
        C.matrix <- eval(parse(text=paste0("makeContrasts(",  ct, ",levels = dm)")))
        fit1 <- contrasts.fit(fit, contrasts =C.matrix)
        fit1 <- eBayes(fit1)
        tt <- topTableF(fit1, sort ="none", n = Inf)
        pv <- tt$P.Value
        P5 <- max(sum(pv<=.05, na.rm = T), 1)/max(sum(pv>.75, na.rm = T)/5, 1) # sometime pv has missing value NA make P5 cannot calculable
        if(print.progress)cat("ConCov = ", names(ConCov)[i], ", P5 = ", P5,  "\n")
        if(WorstP5[[Iter]] > P5) {
          WorstCov <- i
          WorstP5[[Iter]] <- P5
          PvList[[Iter]] <- pv
          names(PvList)[Iter] <- names(WorstP5)[Iter] <- names(ConCov)[i] 
        }
      }
      Iter <- Iter +1
      DelCov <- cbind(DelCov, ConCov[WorstCov])
      ConCov <- ConCov[-WorstCov]
      if(ncol(ConCov)==0)break
    }
    WorstP5 <- do.call("c", WorstP5)
  }else{WorstP5 <- numeric(0)}
  # P <- data.frame(t(do.call("rbind", PvList)))
  res <- list(WorstP5 = WorstP5, VarCov = VarCov, FixCov = FixCov, svamethod = svamethod)
  res
}

DataBS_sva_only <- function(counts, FixCov, VarCov, alpha, print.progress = F, svamethod = F){
  # if(svamethod)VarCov <- svacov(counts = counts, FixCov = FixCov, VarCov = VarCov)
  BSout <- bs_sva_only(counts = counts, FixCov = FixCov, VarCov = VarCov, print.progress = print.progress, svamethod = svamethod)
  BS <- BSout$WorstP5
  S <- vapply(1:length(alpha), function(i)sum(cummax(BS) >alpha[i]), FUN.VALUE = 1.0)
  return(list(BS = BS, S = S, alpha = alpha, FixCov = BSout$FixCov, VarCov = BSout$VarCov, svamethod = BSout$svamethod))
}

# this function is to calculate the correlation of main
# factor of interest and the predicted values obtaining 
# from the linear regression model of the main factor of interest on 
# the residual of TrueCov on Projection matrix surrogate variates
# i.e., the part of TrueCov that cannot be explained by surrogate variables

cor_FixCov_resid_sv <- function(FixCov, sv, TrueCov0){
  if(ncol(sv)==0|ncol(TrueCov0)==0){
    out <- 0
  }else{
    sv <- data.frame(sv)
    sv <- model.matrix(as.formula(paste0("~", paste0(colnames(sv), collapse = "+"))), data = sv)
    P <- sv%*%solve(t(sv)%*%sv)%*%t(sv)
    Id <- diag(nrow(sv))
    XTrueCov <- model.matrix(as.formula(paste0("~", paste0(colnames(TrueCov0), collapse = "+"))), data = TrueCov0)
    lmout <- lm(data.matrix(FixCov) ~(Id-P)%*%XTrueCov)
    out <- as.numeric(unname(cor(data.matrix(FixCov), lmout$fitted.values)))  
    }
    out
}

## function calculate pauc, Rt, Vt, St, and correlation of sva covariates and the true covariate used to simulate the data from the case only use sva without additional covariates ----------
# just for the Line main effect ------
# sva_only_compare <- function(sim, TrueSelCov, UnknownCov = NULL){
#   FixCov <- sim$SimCntOut$FixCov
#   TrueCov <- sim$SimCntOut$VarCov[TrueSelCov]
#   AllCov <- cbind(sim$SimCntOut$FixCov, sim$SimCntOut$VarCov)
#   # DataBSOut <- DataBS_sva_only(counts = sim$SimCnt$SimCnt, FixCov = FixCov, VarCov = FixCov[, 0], alpha = seq(1, 3, length = 10), print.progress = F, svamethod = T)
#   # VoomPvOut <- VoomPv(counts = sim$SimCnt$SimCnt, AllCov = cbind(DataBSOut$FixCov, DataBSOut$VarCov))
#   # PlotVoomPvOut <- PlotVoomPv(VoomPvOut)
#   svacov0 <- cov_sva_only(counts = sim$SimCnt$SimCnt,FixCov = FixCov, VarCov = FixCov[,0])
#   names(svacov0) <- paste0(names(svacov0), "0")
#   svacov1 <- cov_sva_only(counts = sim$SimCnt$SimCnt,FixCov = FixCov, VarCov = TrueCov)
#   names(svacov1) <- paste0(names(svacov1), "1")
#   VoomPvOut0 <- llply(svacov0, function(i)VoomPv(counts = sim$SimCnt$SimCnt, AllCov = cbind(FixCov, i)))
#   VoomPvOut1 <- llply(svacov1, function(i)VoomPv(counts = sim$SimCnt$SimCnt, AllCov = cbind(FixCov, TrueCov, i)))
#   numsva0 <- llply(svacov0, function(i)ncol(i))
#   names(numsva0) <- paste0(names(numsva0), ".numsv")
#   numsva0 <- do.call("c", numsva0)
#   corFixCov0 <- llply(svacov0, function(i)cor_FixCov_resid_sv(FixCov = FixCov, sv = i, TrueCov = sim$SimCntOut$VarCov[c(TrueSelCov, UnknownCov)]))
#   names(corFixCov0) <- paste0(names(corFixCov0), ".corFixCov")
#   corFixCov0 <- do.call("c", corFixCov0)
#   
#   corFixCov1 <- llply(svacov1, function(i)cor_FixCov_resid_sv(FixCov = FixCov, sv = cbind(sim$SimCntOut$VarCov[TrueSelCov],i), TrueCov = sim$SimCntOut$VarCov[c(TrueSelCov, UnknownCov)]))
#   names(corFixCov1) <- paste0(names(corFixCov1), ".corFixCov")
#   corFixCov1 <- do.call("c", corFixCov1)
#   
#   
#   numsva1 <- llply(svacov1, function(i)ncol(i))
#   names(numsva1) <- paste0(names(numsva1), ".numsv")
#   numsva1 <- do.call("c", numsva1)
#   
#   pauc_out_general <- function(p, lab= as.factor(as.numeric(!(sim$SimCnt$IndSample%in%sim$SimCnt$EEGene$Line)))){
#     roc.out <- roc(1-p, lab) # plot(roc.out)
#     roc.ind <- sum(roc.out$fpr<=.05)
#     roc.min <- roc.out$cutoffs[roc.ind]
#     pauc <- auc(roc.out, min =roc.min)
#     qv <- jabes.q(p)
#     R<- sum(qv <=.05)
#     V <- sum((qv <= .05)*(!as.numeric(as.character(lab))))
#     FDR <- V/max(1, R)
#     S = R - V
#     out <- c(NTP = S, R = R, V = V,FDR = FDR, PAUC = pauc)
#     return(out)
#   }
#   
#   simVal <- pauc_out_general(sim$SimCntOut$res$OWN$VoomPvOutER$pvs[, "Line"])
#   names(simVal) <- paste0("sim.", names(simVal))
#   svaVal0 <- llply(VoomPvOut0, function(x) pauc_out_general(x$pvs[, "Line"]))
#   svaVal0 <- do.call("c", svaVal0)
#   svaVal1 <- llply(VoomPvOut1, function(x) pauc_out_general(x$pvs[, "Line"]))
#   svaVal1 <- do.call("c", svaVal1)
#   
#   corout0 <- llply(svacov0, function(x){
#     corout <- laply(AllCov, function(i){
#       dat <- cbind(y = as.numeric(i),x)
#       lmout <- lm(y~., data = dat)
#       out <- cor(lmout$fitted.values, dat$y)
#       out
#     }
#     )
#     names(corout) <- names(AllCov)
#     if(length(corout)==0) corout <- NA
#     corout
#   })
#   corout0 <- do.call("c", corout0)
#   
#   corout1 <- llply(svacov1, function(x){
#     corout <- laply(AllCov, function(i){
#       dat <- cbind(y = as.numeric(i),x)
#       lmout <- lm(y~., data = dat)
#       out <- cor(lmout$fitted.values, dat$y)
#       out
#     }
#     )
#     names(corout) <- names(AllCov)
#     if(length(corout)==0) corout <- NA
#     corout
#   })
#   corout1 <- do.call("c", corout1)
#   
#   
#   
#   # c(svaVal0, svaVal1, simVal, corout0, numsva0, numsva1)
#   c(simVal, svaVal0, svaVal1,  numsva0, numsva1, corFixCov0,corFixCov1, corout0, corout1)
# }


sva_only_compare <- function(sim, TrueSelCov,TrueCov0 = NULL){
  FixCov <- sim$SimCntOut$FixCov
  TrueCov <- sim$SimCntOut$VarCov[TrueSelCov]
  AllCov <- cbind(sim$SimCntOut$FixCov, sim$SimCntOut$VarCov)
  # DataBSOut <- DataBS_sva_only(counts = sim$SimCnt$SimCnt, FixCov = FixCov, VarCov = FixCov[, 0], alpha = seq(1, 3, length = 10), print.progress = F, svamethod = T)
  # VoomPvOut <- VoomPv(counts = sim$SimCnt$SimCnt, AllCov = cbind(DataBSOut$FixCov, DataBSOut$VarCov))
  # PlotVoomPvOut <- PlotVoomPv(VoomPvOut)
  svacov0 <- cov_sva_only(counts = sim$SimCnt$SimCnt,FixCov = FixCov, VarCov = FixCov[,0])
  names(svacov0) <- paste0(names(svacov0), "0")
  svacov1 <- cov_sva_only(counts = sim$SimCnt$SimCnt,FixCov = FixCov, VarCov = TrueCov)
  names(svacov1) <- paste0(names(svacov1), "1")
  VoomPvOut0 <- llply(svacov0, function(i)VoomPv(counts = sim$SimCnt$SimCnt, AllCov = cbind(FixCov, i)))
  VoomPvOut1 <- llply(svacov1, function(i)VoomPv(counts = sim$SimCnt$SimCnt, AllCov = cbind(FixCov, TrueCov, i)))
  numsva0 <- llply(svacov0, function(i)ncol(i))
  names(numsva0) <- paste0(names(numsva0), ".numsv")
  numsva0 <- do.call("c", numsva0)
  corFixCov0 <- llply(svacov0, function(i)cor_FixCov_resid_sv(FixCov = FixCov, 
                                                              sv = i, 
                                                              TrueCov = sim$SimCntOut$VarCov[TrueCov0]))
  names(corFixCov0) <- paste0(names(corFixCov0), ".corFixCov")
  corFixCov0 <- do.call("c", corFixCov0)
  
  corFixCov1 <- llply(svacov1, function(i)cor_FixCov_resid_sv(FixCov = FixCov, 
                                                              sv = cbind(sim$SimCntOut$VarCov[TrueSelCov],i), 
                                                              TrueCov = sim$SimCntOut$VarCov[TrueCov0]))
  names(corFixCov1) <- paste0(names(corFixCov1), ".corFixCov")
  corFixCov1 <- do.call("c", corFixCov1)
  
  sim.corFixCov <- cor_FixCov_resid_sv(FixCov = FixCov, 
                      sv = sim$SimCntOut$VarCov[TrueSelCov], 
                      TrueCov = sim$SimCntOut$VarCov[TrueCov0])
  names(sim.corFixCov) <- "sim.corFixCov"
  
  numsva1 <- llply(svacov1, function(i)ncol(i))
  names(numsva1) <- paste0(names(numsva1), ".numsv")
  numsva1 <- do.call("c", numsva1)
  
  pauc_out_general <- function(p, lab= as.factor(as.numeric(!(sim$SimCnt$IndSample%in%sim$SimCnt$EEGene$Line)))){
    roc.out <- roc(1-p, lab) # plot(roc.out)
    roc.ind <- sum(roc.out$fpr<=.05)
    roc.min <- roc.out$cutoffs[roc.ind]
    pauc <- auc(roc.out, min =roc.min)
    qv <- jabes.q(p)
    R<- sum(qv <=.05)
    V <- sum((qv <= .05)*(!as.numeric(as.character(lab))))
    FDR <- V/max(1, R)
    S = R - V
    out <- c(NTP = S, R = R, V = V,FDR = FDR, PAUC = pauc)
    return(out)
  }
  
  simVal <- pauc_out_general(sim$SimCntOut$res$OWN$VoomPvOutER$pvs[, "Line"])
  names(simVal) <- paste0("sim.", names(simVal))
  svaVal0 <- llply(VoomPvOut0, function(x) pauc_out_general(x$pvs[, "Line"]))
  svaVal0 <- do.call("c", svaVal0)
  svaVal1 <- llply(VoomPvOut1, function(x) pauc_out_general(x$pvs[, "Line"]))
  svaVal1 <- do.call("c", svaVal1)
  
  corout0 <- llply(svacov0, function(x){
    corout <- laply(AllCov, function(i){
      dat <- cbind(y = as.numeric(i),x)
      lmout <- lm(y~., data = dat)
      out <- cor(lmout$fitted.values, dat$y)
      out
    }
    )
    names(corout) <- names(AllCov)
    if(length(corout)==0) corout <- NA
    corout
  })
  corout0 <- do.call("c", corout0)
  
  corout1 <- llply(svacov1, function(x){
    corout <- laply(AllCov, function(i){
      dat <- cbind(y = as.numeric(i),x)
      lmout <- lm(y~., data = dat)
      out <- cor(lmout$fitted.values, dat$y)
      out
    }
    )
    names(corout) <- names(AllCov)
    if(length(corout)==0) corout <- NA
    corout
  })
  corout1 <- do.call("c", corout1)
  
  
  
  # c(svaVal0, svaVal1, simVal, corout0, numsva0, numsva1)
  c(simVal, svaVal0, svaVal1,  numsva0, numsva1, sim.corFixCov, corFixCov0,corFixCov1, corout0, corout1)
}
# 
# PlotVoomPv(VoomPvOut0[[1]])
# PlotVoomPv(VoomPvOut0[[2]])

# Function wrap all everything---------
FSRnSimAllBS_sva_only <- function(BData, mData, nGene, option, ErrType,  BSim, mSim, amax, gam0Data, gam0, ncores, print.progress, RealDataOutPath){
  RealDataOut <- readRDS(file = paste0(RealDataOutPath, "/", RealDataOutPath,   "_BData", BData, "_mData", mData, "_gam0_", gam0Data, ".rds"))
  TrueSelCov <- FSRSelCov(FSRAnalysisAllOut=RealDataOut, option=option, ErrType = ErrType)
  SimPath <- paste0("SimBS_", "BData", BData, "_mData", mData, "_",option, "_gamma", ErrType, "_gam0_", gam0, "_nGene", nGene, "_Bsim", BSim, "_mSim", mSim)  
  dir.create(path = SimPath, showWarnings = F)
  nf <- list.files(path = SimPath, pattern = "nrep")
  out <- mclapply(nf, function(i){
    print(i)
    sim <- readRDS(file = paste0(SimPath, "/", i))
    simout <- sva_only_compare(sim, TrueSelCov)
    simout
  }, mc.cores = ncores)
  
  oo <- do.call("rbind", out)
  
  # oo <- list(res = laply(out, function(x)x[1:10]), 
  #      svacor = lapply(out, function(x)x[11:length(x)]))
  saveRDS(oo, file = paste0(SimPath, "/sva_only_nocov2.rds"))
  saveRDS(oo, file = paste0("AllResults/",SimPath, "/sva_only_nocov2.rds"))
  dir.create(path = "AllResults_sva_only/", showWarnings = F)
  dir.create(path = paste0("AllResults_sva_only/",SimPath), showWarnings = F)
  saveRDS(oo, file = paste0("AllResults_sva_only/",SimPath, "/sva_only_nocov2.rds"))
  oo
}



##Function to calculate NTP, R, V, PAUC of Line main effect for each simulated data---------
pauc_out_sim_Line <- function(sim){
  option <- c("RX", "ORX", "WN", "OWN")
  ErrType <- c("ER", "RE")
  oE <- expand.grid(option = option, ErrType = ErrType)
  out <- lapply(1:nrow(oE), function(i){
    option <-as.character(oE$option[i]); ErrType <- as.character(oE$ErrType[i])
    p <- sim$SimCntOut$res[[option]][[paste0("VoomPvOut", ErrType)]]$pvs[, "Line"]
    lab= as.factor(as.numeric(!(sim$SimCnt$IndSample%in%sim$SimCnt$EEGene$Line)))
    roc.out <- roc(1-p, lab) # plot(roc.out)
    roc.ind <- sum(roc.out$fpr<=.05)
    roc.min <- roc.out$cutoffs[roc.ind]
    pauc <- auc(roc.out, min =roc.min)
    qv <- jabes.q(p)
    R<- sum(qv <=.05)
    V <- sum((qv <= .05)*(!as.numeric(as.character(lab))))
    FDR <- V/max(1, R)
    S = R - V
    out <- c(NTP = S, R = R, V = V,FDR = FDR, PAUC = pauc)
    names(out) <- paste0(names(out), paste0(".", option, ".", ErrType))
    out
  })
  out <- do.call("c", out)
  return(out)
}


FSRnSimAllBS_LineEffect <- function(BData, mData, nGene, option, ErrType,  BSim, mSim, amax, gam0Data, gam0, ncores, print.progress, RealDataOutPath){
  SimPath <- paste0("SimBS_", "BData", BData, "_mData", mData, "_",option, "_gamma", ErrType, "_gam0_", gam0, "_nGene", nGene, "_Bsim", BSim, "_mSim", mSim)  
  dir.create(path = SimPath, showWarnings = F)
  nf <- list.files(path = SimPath, pattern = "nrep")
  out <- mclapply(nf, function(i){
    print(i)
    sim <- readRDS(file = paste0(SimPath, "/", i))
    simout <- pauc_out_sim_Line(sim)
    simout
  }, mc.cores = ncores)
  oo <- do.call("rbind", out)
  saveRDS(oo, file = paste0("AllResults/",SimPath, "/LineEffect.rds"))
  oo
}

## Function to calculate pacu_out of OnlyLine, AllCov, Oracle, OldBS (published in jabes 2015)
pauc_out <- function(p, lab,  method){
  roc.out <- roc(1-p, lab) # plot(roc.out)
  roc.ind <- sum(roc.out$fpr<=.05)
  roc.min <- roc.out$cutoffs[roc.ind]
  pauc <- auc(roc.out, min =roc.min)
  qv <- jabes.q(p)
  R<- sum(qv <=.05)
  V <- sum((qv <= .05)*(!as.numeric(as.character(lab))))
  FDR <- V/max(1, R)
  S = R - V
  out <- c(NTP = S, R = R, V = V,FDR = FDR, PAUC = pauc)
  names(out) <- paste0(names(out), paste0(".", method))
  out
}

pauc_out_sim_Line_All <- function(sim, BData){
  counts <- sim$SimCnt$SimCnt
  lab <-as.factor(as.numeric(!(sim$SimCnt$IndSample%in%sim$SimCnt$EEGene$Line))) 
   FixCov <- sim$SimCntOut$FixCov
  VarCov <- sim$SimCntOut$VarCov
  # OnlyLine 
  OnlyLine <- pauc_out(p = VoomPv(counts = counts, AllCov = FixCov)$pvs[,"Line"], lab = lab,  method = "OnlyLine")
  
  # AllCov
  All <- pauc_out(p = VoomPv(counts = counts, AllCov = cbind(FixCov, VarCov))$pvs[,"Line"], lab = lab,  method = "All")

  #Oracle
  RealOutData <- readRDS(paste0("RealDataOutBS/RealDataOutBS_BData", BData, "_mData1_gam0_0.001.rds"))
  TrueCov <- names(RealOutData$res$OWN$BestCovOut$BestRE)
  TrueCovSet <- VarCov[TrueCov]
  Oracle <- pauc_out(p = VoomPv(counts = counts, AllCov = cbind(FixCov,TrueCovSet))$pvs[,"Line"], lab = lab,  method = "Oracle")
  # jabes.bs
  BestCov <- jabes.bs(counts = counts, FixCov = FixCov, VarCov = VarCov, print.progress = F)$BestCov
  
  OldBS <- pauc_out(p = VoomPv(counts = counts, AllCov = cbind(FixCov,BestCov))$pvs[,"Line"], lab = lab, method = "OldBS")
  option <- c("RX", "ORX", "WN", "OWN")
  ErrType <- c("ER", "RE")
  oE <- expand.grid(option = option, ErrType = ErrType)
  out <- lapply(1:nrow(oE), function(i){
    option <-as.character(oE$option[i]); ErrType <- as.character(oE$ErrType[i])
    p <- sim$SimCntOut$res[[option]][[paste0("VoomPvOut", ErrType)]]$pvs[, "Line"]
    pauc_out(p = p,lab = lab,  method = paste0( option, ".", ErrType))
    
  })
  out <- do.call("c", out)
  out <- c(out, OnlyLine, All, Oracle, OldBS)
  return(out)
}


FSRnSimAllBS_LineEffect_All <- function(BData, mData, nGene, option, ErrType,  BSim, mSim, amax, gam0Data, gam0, ncores, print.progress, RealDataOutPath){
  SimPath <- paste0("SimBS_", "BData", BData, "_mData", mData, "_",option, "_gamma", ErrType, "_gam0_", gam0, "_nGene", nGene, "_Bsim", BSim, "_mSim", mSim)  
  dir.create(path = SimPath, showWarnings = F)
  nf <- list.files(path = SimPath, pattern = "nrep")
  out <- mclapply(nf, function(i){
    print(i)
    sim <- readRDS(file = paste0(SimPath, "/", i))
    simout <- pauc_out_sim_Line_All(sim, BData)
    simout
  }, mc.cores = ncores)
  oo <- do.call("rbind", out)
  saveRDS(oo, file = paste0("AllResults/",SimPath, "/LineEffectAllLineOnlyOracle.rds"))
  oo
}


## Calculate FSR of the old BS selection in Nguyen et al (JABES 2015)----


FSROldBS <- function(sim, BData){
  counts <- sim$SimCnt$SimCnt
  FixCov <- sim$SimCntOut$FixCov
  VarCov <- sim$SimCntOut$VarCov
  RealOutData <- readRDS(paste0("RealDataOutBS/RealDataOutBS_BData", BData, "_mData1_gam0_0.001.rds"))
  TrueCov <- names(RealOutData$res$OWN$BestCovOut$BestRE)
  BestCov <- names(jabes.bs(counts = counts, FixCov = FixCov, VarCov = VarCov, print.progress = F)$BestCov)
  R <- length(BestCov) # S in the Wu2007 paper
  S <- length(intersect(TrueCov, BestCov))# I in the Wu2007 paper
  U <- R - S
  FSR <- U/(R+1)
  out <- c(S.OldBS=S, R.OldBS=R, U.OldBS=U, FSR.OldBS=FSR)
  out
}


FSROldBSAll <- function(BData, mData, nGene, option, ErrType,  BSim, mSim, amax, gam0Data, gam0, ncores, print.progress, RealDataOutPath){
  SimPath <- paste0("SimBS_", "BData", BData, "_mData", mData, "_",option, "_gamma", ErrType, "_gam0_", gam0, "_nGene", nGene, "_Bsim", BSim, "_mSim", mSim)  
  dir.create(path = SimPath, showWarnings = F)
  nf <- list.files(path = SimPath, pattern = "nrep")
  out <- mclapply(nf, function(i){
    print(i)
    sim <- readRDS(file = paste0(SimPath, "/", i))
    simout <- FSROldBS(sim, BData)
    simout
  }, mc.cores = ncores)
  oo <- do.call("rbind", out)
  saveRDS(oo, file = paste0("AllResults/",SimPath, "/FSRoutOldBS.rds"))
  oo
}




# le <- readRDS(file = paste0("AllResults/",SimPath, "/LineEffect.rds"))
# mle <- apply(le, 2, mean)
# se.mle <- apply(le, 2, sd)/10
# se.mle
# mle
# FDR <- mle


### Simulation with unknown relevant covariates ----------

## First we use simcounts from ModelSize 8 and 7, 
## Let's suppose the unknown covariates is RFI first and see what happen

FSRnSimAllBSUnknownCov <- function(BData, mData, nGene, option, ErrType,  BSim, mSim, amax,  gam0, ncores, print.progress, UnknownCov){
  SimPath <- paste0("SimBS_", "BData", BData, "_mData", mData, "_",option, "_gamma", ErrType, "_gam0_", gam0, "_nGene", nGene, "_Bsim", BSim, "_mSim", mSim)  
  dir.create(path = SimPath, showWarnings = F)
  nf <- list.files(path = SimPath, pattern = "nrep")
  if (BData ==6){
    TrueCov0 <- c( "RFI", "lymp", "baso", "RINa", "Block", "neut", "Concb", "mono")
  }else if(BData ==5){
    TrueCov0 <- c("lymp", "baso", "RINa", "Block", "neut", "Concb", "mono")
  } else if (BData ==4){
    TrueCov0 <- c("baso", "RINa", "Block", "neut", "Concb", "mono")
  }else if(BData == 3){
    TrueCov0 <- c( "Concb", "mono")
  }else if(BData ==2){
    TrueCov0 <- c( "mono")
  }else  TrueCov0 <- NULL
  
  pauc_out0 <- function(p, lab,  method){
    roc.out <- roc(1-p, lab) # plot(roc.out)
    roc.ind <- sum(roc.out$fpr<=.05)
    roc.min <- roc.out$cutoffs[roc.ind]
    pauc <- auc(roc.out, min =roc.min)
    qv <- jabes.q(p)
    R<- sum(qv <=.05)
    V <- sum((qv <= .05)*(!as.numeric(as.character(lab))))
    FDR <- V/max(1, R)
    S = R - V
    out <- c(NTP = S, R = R, V = V,FDR = FDR, PAUC = pauc)
    names(out) <- paste0(paste0(method, "."), names(out))
    out
  }
  
  out <- mclapply(nf, function(i){ # i <- nf[3]
    print(i)
    sim <- readRDS(file = paste0(SimPath, "/", i))
    
    simout <- FSRAnalysisBS(counts=sim$SimCnt$SimCnt, FixCov = sim$SimCntOut$FixCov, 
                            VarCov = sim$SimCntOut$VarCov[-1][setdiff(names(sim$SimCntOut$VarCov[-1]), UnknownCov)], # delete Diet which is the covariate correlated with RFI
                             option = option, B= BSim, m = mSim, 
                             amax= amax, gam0 = gam0, ncores=ncores, 
                             print.progress=print.progress, svamethod = F)
    lab <-as.factor(as.numeric(!(sim$SimCnt$IndSample%in%sim$SimCnt$EEGene$Line))) 
    FSRres <- pauc_out0(p = simout$VoomPvOutER$pvs[, "Line"], lab = lab, method = "FSR")
    R <- length(names(simout$BestCovOut$BestER))
    S <- length(intersect(names(simout$BestCovOut$BestER), TrueCov0))
    U <- R - S
    FSP <-U/(R+1) 
    SVAres <- sva_only_compare(sim = sim, TrueSelCov = names(simout$BestCovOut$BestER), TrueCov0 = TrueCov0)
    corFixFSR <- cor_FixCov_resid_sv(FixCov = sim$SimCntOut$FixCov, 
                                     sv = sim$SimCntOut$VarCov[names(simout$BestCovOut$BestER)], 
                                     TrueCov0 = sim$SimCntOut$VarCov[TrueCov0])
   c(FSR.S = S, FSR.U = U, FSR.R = R,FSR.FSP = FSP,FSR.corFixCov =corFixFSR,  FSRres, SVAres)
    }, mc.cores = 1)
  oo <- do.call("rbind", out)
  # saveRDS(oo, file = paste0("AllResults/",SimPath, "/FSRUnknownCov_", paste0(UnknownCov, collapse = "_"),".rds"))
  saveRDS(oo, file = paste0("AllResults/",SimPath, "/FSRUnknownCov_noDiet_", paste0(UnknownCov, collapse = "_"),".rds"))
  oo
}


# DataBSOut <- DataBS_sva_only(counts = sim$SimCnt$SimCnt, FixCov = FixCov, VarCov = FixCov[, 0], alpha = seq(1, 3, length = 10), print.progress = F, svamethod = T)
# VoomPvOut <- VoomPv(counts = sim$SimCnt$SimCnt, AllCov = cbind(DataBSOut$FixCov, DataBSOut$VarCov))
# PlotVoomPvOut <- PlotVoomPv(VoomPvOut)
# svacov0 <- cov_sva_only(counts = counts,FixCov = FixCov, VarCov = FixCov[,0])
# names(svacov0) <- paste0(names(svacov0), "0")
# svacov1 <- cov_sva_only(counts = counts,FixCov = FixCov, VarCov = TrueCov)
# names(svacov1) <- paste0(names(svacov1), "1")
# VoomPvOut0 <- llply(svacov0, function(i)VoomPv(counts = counts, AllCov = cbind(FixCov, i)))
# VoomPvOut1 <- llply(svacov1, function(i)VoomPv(counts = counts, AllCov = cbind(FixCov, TrueCov, i)))
# PlotVoomPv(VoomPvOut0[[1]])
# vo <- VoomPv(counts = counts, AllCov = cbind(FixCov, TrueCov[-1]))
# PlotVoomPv(vo)


