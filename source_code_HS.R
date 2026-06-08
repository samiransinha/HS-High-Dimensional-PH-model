# ==============================================================================
# Simulation and estimation script
# Order requested:
#   1. All functions
#   2. Data generation
#   3. Estimation using different methods
#
# Note: This keeps your original mathematical code mostly unchanged. The main
# ordering fix is that X is scaled directly after simulation, before mydata exists.
# ===============================================================================

# ==============================================================================
# 1. Libraries
# ===============================================================================

library(survival)
library(glmnet)
library(MASS)          # ginv(); older version may be needed for R 4.3 compatibility
library(LaplacesDemon)
library(cluster)
library(parallel)

# ==============================================================================
# 2. All functions
# ===============================================================================

scale_to_minus1_1 <- function(x) {
  2 * (x - min(x)) / (max(x) - min(x)) - 1
}

split_vector <- function(p, nparts) {
  v <- 1:p
  groups <- cut(seq_along(v), breaks = nparts, labels = FALSE)
  parts <- split(v, groups)
  sizes <- sapply(parts, length)
  list(parts = parts, sizes = sizes)
}

simulate_survival_data <- function(n, p, s, beta.true, gamma, plot_km = TRUE) {
  p2 <- s * p

  Sigma <- matrix(0, ncol = p, nrow = p)
  Sigma[1:p2, 1:p2] <- 0.5
  Sigma[1:p2, -(1:p2)] <- 0.3
  Sigma[-(1:p2), 1:p2] <- 0.3
  Sigma[-(1:p2), -(1:p2)] <- 0.6
  diag(Sigma) <- 1

  mychol <- chol(Sigma)
  X <- t(t(mychol) %*% matrix(rnorm(n * p), nrow = p))

  ## IMPORTANT: scale X directly. 
 # X <- apply(X, 2, scale_to_minus1_1)

  Z <- cbind(runif(n, -1, 1), rbinom(n, 1, 0.5))
  eta <- X %*% beta.true + Z %*% gamma

  ## For independent features, replace the above X and eta with:
  ## X <- matrix(rnorm(n * p), ncol = p)
  ## eta <- X %*% beta.true
  if(p==1000){
  myshape <- s * p / 4} else{# for p=1000
  myshape= s*p/10 # for p=5000
}
  myt <- (-log(1 - runif(n)) / exp(eta))^(1 / myshape)
  myc <- quantile(myt, prob = c(0.8, 0.9))
  CT <- runif(n, myc[1], myc[2])

  myv <- pmin(myt, CT)
  delta <- 1 * (myt < CT)
  myv <- myv / max(myv)

  mydata <- data.frame(
    V = myv,
    delta = delta,
    X = I(as.matrix(X)),
    Z = I(as.matrix(Z))
  )
  mydata <- mydata[order(mydata$V), ]

  if (plot_km) {
    fit <- survfit(Surv(V, delta) ~ 1, data = mydata)
    plot(fit)
  }

  list(
    mydata = mydata,
    censoring_indicator_mean = mean(delta),
    true.idx = which(beta.true != 0)
  )
}

lasso_initialization <- function(mydata, p, max_selected_fraction = 0.05) {
  if ("Z" %in% names(mydata)) {
    Xfull <- cbind(mydata$X, mydata$Z)

    pX <- ncol(mydata$X)
    pZ <- ncol(mydata$Z)

    mypenalty <- c(
      rep(1, pX),  # X variables: penalized
      rep(0, pZ)   # Z variables: unpenalized
    )

    cvfit <- cv.glmnet(
      x = Xfull,
      y = Surv(mydata$V, mydata$delta),
      family = "cox",
      type.measure = "C",
      alpha = 1,
      penalty.factor = mypenalty
    )

    outnnet <- glmnet(
      x = Xfull,
      y = Surv(mydata$V, mydata$delta),
      family = "cox",
      lambda = cvfit$lambda.min,
      alpha = 1,
      penalty.factor = mypenalty
    )
  } else {
    cvfit <- cv.glmnet(
      x = mydata$X,
      y = Surv(mydata$V, mydata$delta),
      family = "cox",
      type.measure = "C",
      alpha = 1
    )

    outnnet <- glmnet(
      x = mydata$X,
      y = Surv(mydata$V, mydata$delta),
      family = "cox",
      lambda = cvfit$lambda.min,
      alpha = 1
    )
  }

  betahatlasso <- as.numeric(coef(outnnet))
  beta.est <- betahatlasso

  if (length(which(betahatlasso[1:p] != 0)) > max_selected_fraction * p) {
    n_keep <- ceiling(max_selected_fraction * length(betahatlasso[1:p]))
    top_indices <- order(abs(betahatlasso[1:p]), decreasing = TRUE)[1:n_keep]

    beta_filtered <- numeric(length(betahatlasso[1:p]))
    beta_filtered[top_indices] <- betahatlasso[top_indices]

    beta.est <- c(beta_filtered, betahatlasso[-(1:p)])
  }

  list(
    beta.est = beta.est,
    betahatlasso = betahatlasso,
    cvfit = cvfit,
    outnnet = outnnet
  )
}

prepare_global_objects <- function(mydata, npart = 10, varz_value = 1.10) {
  ## This function assigns the global objects used by the original estimation
  ## functions below.

  n <<- nrow(mydata$X)
  p <<- ncol(mydata$X)
  V <<- mydata$V
  delta <<- as.numeric(mydata$delta)
  X <<- mydata$X

  Z <<- mydata$Z
  q_z <<- ncol(Z)
  varz <<- varz_value
  covmat <<- cbind(X, Z)

  mpq <<- p + q_z
  npara <<- p + q_z + 1
  npart <<- npart

  outs <- split_vector(p, npart)
  istarts <<- c(1, cumsum(outs$sizes)[-length(outs$sizes)] + 1)
  iends <<- cumsum(outs$sizes)

  chunks <<- vector("list", npart)
  for (i in 1:npart) {
    chunks[[i]] <<- X[, istarts[i]:iends[i]]
  }

  a.alpha <<- 0
  b.alpha <<- 2
  a.varpi <<- 0
  b.varpi <<- 2

  myknots <<- quantile(V[delta == 1], prob = c(0.2, 0.4, 0.6, 0.8))
  m <<- length(myknots)

  knotmat <<- matrix(rep(myknots, n), ncol = m, byrow = TRUE)
  knotadjv <<- matrix(rep(V, m), byrow = FALSE, ncol = m) - knotmat

  knotadjvsq <<- knotadjv^2
  knotadjvpwr3 <<- knotadjvsq * knotadjv
  knotadjvpwr4 <<- knotadjvsq * knotadjvsq
  myknotspwr3 <<- matrix(rep(myknots^3, n), ncol = m, byrow = TRUE)

  invisible(list(
    n = n, p = p, q_z = q_z, mpq = mpq, npara = npara,
    npart = npart, m = m, varz = varz
  ))
}

lambda.estimation <- function(mydata){
  
 
  # prior parameters 
  a.alpha=0 
  b.alpha=100
  a.varpi=0
  b.varpi=100 
  const1= sqrt(2*pi)
  diag_b.alpha2=diag(1/rep(b.alpha^2, m))
  
  gunction=function(para){
    
    term1.1= exp(0.5*para[m+1]) 
    term1.2= exp(para[m+1]) 
    term2= exp(para[1:m])
    ##
    term7=exp(-0.5*knotadjvsq/term1.2)
    lambda.construct= term2*t(term7) 
    ## 
    lambda=colSums(lambda.construct)
    lambda=pmax(lambda, 1e-300)
    
    term9=(pnorm(knotadjv/term1.1)+pnorm(knotmat/term1.1)-1)
    mat1=const1*term1.1*t(term2* t(term9))
    caplambda=rowSums(mat1)
    caplambda=pmax(caplambda, 1.17e-17)
    
    # Computation of log-posterior density without the proportionality constant 
    logpost=(sum(delta*(log(lambda)))- sum(caplambda)  
             -0.5* sum((para[1:m]-a.alpha)^2)/b.alpha^2
             -0.5*(para[m+1]- a.varpi)^2/b.varpi^2
    )
    
    return(-logpost)
  }
  #hessian(gunction, inipara)
  
  inipara=c(rep(0.1,m), 0.1)
  outm= optim(inipara, gunction, method="L-BFGS-B")#, hessian=TRUE)
  return(outm$par)
}






estimatepsi=function(beta.est, initial.para, varpi_update_psi){

term3  <- as.numeric(covmat %*% beta.est)
expxbeta <- exp(term3)
expxbeta <- pmin(pmax(expxbeta, 1e-100), 1e8)

# parameters
para=initial.para # should be commented out 
for(irep in 1:4){
alpha  <- para[1:m]
varpi  <- para[m+1]

exp_half_varpi <- exp(0.5 * varpi)   # e^{0.5 varpi}
exp_varpi      <- exp(varpi)         # e^{varpi}
exp_alpha   <- exp(alpha)

##############################
# Kernel quantities
##############################

# φ_{ij}
kernel_val <- exp(-0.5 * knotadjvsq / exp_varpi)

# λ_i components
lambda.construct <- exp_alpha * t(kernel_val)    # m x n
lambda.construct <- pmax(lambda.construct, 1e-300)
lambda <- colSums(lambda.construct)
lambda <- pmax(lambda, 1e-300)

##############################
# CDF terms
##############################

scaled_diff  <- knotadjv / exp_half_varpi
scaled_knots <- -knotmat / exp_half_varpi

Phi_diff  <- pnorm(scaled_diff)
Phi_knots <- pnorm(scaled_knots)

Hij <- Phi_diff - Phi_knots   # n x m

##############################
# Initialize gradient
##############################

deriv_vec <- rep(0, (m + 1))

##############################
# (1) Gradient w.r.t alpha
##############################

for (k in 1:m) {
  deriv_vec[k] =
    sum(delta * lambda.construct[k, ] /lambda) -
    sum(expxbeta * (sqrt(2*pi) * exp_half_varpi) * exp_alpha[k] * Hij[, k]) -
    (alpha[k] - a.alpha) / b.alpha^2
}

##############################
# (2) Gradient w.r.t varpi
##############################

## ---- Term A: derivative of log S_i ----
dSi <- colSums(lambda.construct * t(knotadjvsq)) / (2 * exp_varpi)
termA <- sum(delta * dSi / lambda)

## ---- Term B: derivative of cumulative hazard ----

# standard normal pdf terms
pdf_diff  <- dnorm(scaled_diff)
pdf_knots <- dnorm(scaled_knots)

termB <- 0.5*sqrt(2*pi)*sum(((knotadjv * pdf_diff +
                knotmat * pdf_knots)%*%exp_alpha)*expxbeta)
termC <- 0.5*sqrt(2*pi)*exp_half_varpi*sum((Hij%*%exp_alpha)*expxbeta)

## ---- Final varpi gradient ----
deriv_vec[m+1] =(
  termA +  termB-termC
  -(varpi - a.varpi) / b.varpi^2
)

################
myhessian=matrix(0, ncol=(m+1), nrow=(m+1))  
for(k1 in 1: m){
  for( k2 in k1:m){
myhessian[k1, k2]=-sum(delta*lambda.construct[k1, ]*
       lambda.construct[k2, ]/lambda^2)
                              
  }
} 
for(k1 in 1: m){
  myhessian[k1, k1]=myhessian[k1, k1]+(sum(delta*lambda.construct[k1, ]/lambda) 
  -sum(expxbeta * (sqrt(2*pi) * exp_half_varpi) * exp_alpha[k1] * Hij[, k1]) -
    1 / b.alpha^2)
}

for(k2 in 1: (m-1)){
  for( k1 in (k2+1):m){
  myhessian[k1, k2]= myhessian[k2, k1]
     }
}

for(k1 in 1:m){
myhessian[(m+1), k1]=(
   (0.5/exp_varpi)*sum(delta*lambda.construct[k1, ]*
                                            knotadjvsq[, k1]/lambda)
-(0.5/exp_varpi)*sum(delta*lambda.construct[k1, ]*
                       colSums(t(knotadjvsq)*lambda.construct)/lambda^2)
+0.5* sqrt(2*pi)*exp_alpha[k1]*sum((pdf_diff[, k1]*knotadjv[, k1]+
                     pdf_knots[, k1]*knotmat[, k1])*expxbeta)
-0.5*sqrt(2*pi)*exp_half_varpi*sum(Hij[, k1]*exp_alpha[k1]*expxbeta)
)        
############
 myhessian[k1, (m+1)]=myhessian[(m+1), k1]
}
myhessian[(m+1), (m+1)]=( sum((delta/exp_varpi)* colSums(t(0.25*knotadjvpwr4/exp_varpi-
                                   0.5*knotadjvsq)*lambda.construct)/lambda)
- (0.25/exp_varpi^2)*sum(delta*(colSums(t(knotadjvsq)*lambda.construct)/lambda)^2)
-0.25*sqrt(2*pi)*exp_half_varpi*sum((Hij%*%exp_alpha)*expxbeta)   
+0.25*sqrt(2*pi)* sum(( (pdf_diff*knotadjv+pdf_knots*knotmat)%*%exp_alpha)*expxbeta)
+0.25*sqrt(2*pi)*(1/exp_varpi)*sum(((pdf_diff*knotadjvpwr3+
                                    pdf_knots*myknotspwr3)%*%exp_alpha)*expxbeta)   

-1/b.varpi^2
)
################

if(varpi_update_psi)  {
  #cat("old", myhessian, "\n")
  newton_dir <- -ginv(myhessian) %*% deriv_vec
  ascent_check <- as.numeric(t(deriv_vec) %*% newton_dir)
  
  if (ascent_check > 0) {
    direction <- as.vector(newton_dir)
    #cat("Using Newton direction\n")
  } else {
    grad_norm <- sqrt(sum(deriv_vec^2))
    direction <- as.vector(deriv_vec / grad_norm)
   #cat("Newton direction is not ascent. Using normalized gradient direction\n")
  }
  step_size <- 0.1
  para <- para + step_size * direction
  #print(para);
  #cat("new", myhessian, "\n");
  } else{
    #cat("regular", para, "\n")
    para[1:m]=para[1:m]- 
      1*ginv(myhessian[1:m, 1:m])%*%deriv_vec[1:m]
}
  
}
#print(para)
  #para[m+1]=outm$par[m+1]
  return(para)
}

############################################################
############################################################






###########################################################
############################################################
reg.para.est <- function(mydata, npart,  niter=10000, tol=0.0005, 
                         theta=mytheta, psi=mypsi, acns=1, tau02, update_vartheta){
  eps1=eps2=2
  ncount=0
  logpost=5
  store=rep(0, 5000)
  const1= sqrt(2*pi)
  
  
  para=theta[1:npara]
  etaj=theta[-(1:npara)] 
  alpha=psi[1:m]
  varpi=psi[m+1]
  
  exp_half_varpi= exp(0.5*varpi) 
  exp_varpi= exp(varpi) 
  exp_alpha= exp(alpha)
  term7=exp(-0.5*knotadjvsq/exp_varpi)
  lambda= term7%*%exp_alpha 
  ## 
  lambda=pmax(lambda, 1e-300)
  
  term9=(pnorm(knotadjv/exp_half_varpi)+pnorm(knotmat/exp_half_varpi)-1)
  mat1=const1*exp_half_varpi*t(term9%*%exp_alpha)
  caplambda=as.numeric(pmax(mat1, 1.175e-17))
  ncount=0
  while(eps1>tol  & ncount<niter){
    #  for (it in 1: niter){
    myoldpara=para
    ncount=ncount+1
   # print(c(ncount, eps1, eps2))
    oldpara=para
    oldlogpost=logpost
    store[ncount]=logpost
    
    term3= as.numeric(covmat%*%para[1:mpq])
    #term3= pmin(term3, 15)
    expxbeta= exp(term3)
    expxbeta=pmin(expxbeta, 1e8)
    expxbeta=pmax(expxbeta, 1e-100)
    term4.1= exp(0.5*para[mpq+1])
    term4.2=exp(-para[mpq+1])
    #
    prod_expxbeta_caplambda= expxbeta*caplambda
    
    
    
    
    
    deriv_4_beta=lapply(1:npart, function(r){
      r1= istarts[r]
      r2=iends[r] 
      ##
      prod1= chunks[[r]]*prod_expxbeta_caplambda
      ##
      parametersq=para[r1:r2]^2
      
      deriv_order1=( colSums(delta*chunks[[r]])- colSums(prod1)
                     -para[r1:r2]*term4.2*exp(-etaj[r1:r2])
      )
      #### deriv_order2 is the most time consuming chunk, specifically the following matrix 
      #### multiplication 
      deriv_order2= -(npart+1)*crossprod(chunks[[r]], prod1)
      diag(deriv_order2)= (diag(deriv_order2)- ((1+5*parametersq)/(1+parametersq))*
                             term4.2*exp(-etaj[r1:r2]))
      
      
      #change
      return(para[r1:r2]-acns*solve(deriv_order2, deriv_order1) )
      #    
    })
    para[1:p]= unlist(deriv_4_beta, use.names = FALSE)#as.vector(deriv_4_beta)
    
    #########
    prod1_z= covmat[, -(1:p)]*prod_expxbeta_caplambda
    deriv_order1_z=( colSums(delta*Z)- colSums(prod1_z)
                     -para[(p+1):mpq]/varz)
    deriv_order2_z= -(npart+1)*crossprod(covmat[, -(1:p)], prod1_z)-1/varz
    para[(p+1):mpq]=para[(p+1):mpq]-acns*solve(deriv_order2_z, deriv_order1_z)
    
    
    #########
    
    
    
    if (update_vartheta){
    drv_vartheta=(
      -0.5*(p-1)-1/(tau02*term4.2+1 )
      +0.5*term4.2*sum(exp(-etaj)*para[1:p]^2)
    )
    drv_order2_vartheta=( - tau02*term4.2/(tau02*term4.2+1)^2
                          - 1.5*term4.2*sum(exp(-etaj)*(para[1:p]^2+1))
                          
    )
    para[npara]=para[npara]- acns*drv_vartheta/drv_order2_vartheta 
    }
    ############################
    #####
    deriv_etaj= -1/(1+exp(-etaj))+0.5*term4.2*exp(-etaj)*para[1:p]^2
    deriv2_etaj= -exp(-etaj)/(1+exp(-etaj))^2 -1.5*term4.2*exp(-etaj)*(para[1:p]^2+1)
    
    
    etaj= etaj- acns*deriv_etaj/deriv2_etaj
    
    #####
    
    
    #print(drv)
    #para=para+rnorm(npara)*(.1/ncount^(1/2))
    eps1= max(abs((oldpara-para)/
                    apply(cbind(abs(oldpara), 0.1), 1, max)))
    
    
    # Computation of log-posterior density without the proportionality constant 
    logpost <- (
      sum(delta * (term3 + log(lambda))) -
        sum(caplambda * expxbeta) -
        0.5 * (p - 1) * para[npara] -
        log(1 + 1 / (term4.2 * tau02)) -
        0.5 * sum(para[1:p]^2 * term4.2 * exp(-etaj)) -
        sum(log1p(exp(etaj))) -
        0.5 * sum(para[(p + 1):mpq]^2) / varz
    )
    
    
    #  
    eps2= abs((oldlogpost-logpost)/oldlogpost)
    
    if (is.nan(eps1) || is.nan(eps2)) {
      stop()# next
    }
    
  }
  
  return(list(allpara=c(para, etaj), 
              m=m, 
              p=p,
              q_z=q_z, 
              varz=varz,
              logpost=logpost,
              lambda=lambda
  )   
  ) 
  
}






run_proposed_estimation <- function(mydata, beta.est, npart = 10,
                                    outer_tol = 0.001,
                                    initial_tau02 = 2,
                                    update_tau02 = 1) {
  initial.para <- lambda.estimation(mydata)
  a.varpi <<- initial.para[m + 1]

  #myetaj <- rep(-6, p)
  #myetaj[beta.est[1:p] != 0] <- 6#10
  lasso_nonzero <- beta.est[1:p] != 0
  n_nonzero <- sum(lasso_nonzero)
  
  myetaj <- rep(-4, p)  # stronger shrinkage on zeros
  # Scale inclusion confidence by absolute coefficient size
  if (n_nonzero > 0) {
    abs_coef <- abs(beta.est[1:p][lasso_nonzero])
    scaled <- 4 + 4 * (abs_coef - min(abs_coef)) / (max(abs_coef) - min(abs_coef) + 1e-10)
    myetaj[lasso_nonzero] <- scaled  # range [4, 8] instead of flat 6
  }
  s_hat=sum(beta.est[1:p]!=0)
  s_hat <- max(s_hat, 1)
  s_hat <- min(s_hat, p - 1)
  vartheta_initial= log((s_hat/(p-s_hat))^2)
 
  #c0 <- 0.05 #1e-03
  #myetaj <- log(2*beta.est[1:p]^2 + c0) - vartheta_initial
  #myetaj <- pmin(pmax(myetaj, -4), 6)
  
  
  mypsi <- estimatepsi(beta.est = beta.est, initial.para = initial.para, 
                       varpi_update_psi=FALSE)

  mytheta <- reg.para.est(
    mydata,
    npart = npart,
    niter = 10000,
    tol = 0.1,
    theta = c(beta.est[1:p], beta.est[-(1:p)], vartheta_initial, myetaj),
    psi = mypsi,
    acns = 1,
    tau02 = initial_tau02, update_vartheta=FALSE
  )

  
  ################
  ###############
  ########
  ######
  eps <- 20.00001
  while (eps > outer_tol) {
  #for(i00 in 1:5){
    oldbeta <- mytheta$allpara[1:mpq]
    
    mypsi <- estimatepsi(
      beta.est = mytheta$allpara[1:mpq],
      initial.para = mypsi, varpi_update_psi=TRUE
    )
    
    mytheta <- reg.para.est(
      mydata,
      npart = npart,
      niter = 1000,
      tol = 0.1,
      theta = mytheta$allpara,
      psi = mypsi,
      acns = 1,
      tau02 = update_tau02, update_vartheta=TRUE
    )
    
    old.pr <- 1 / (1 + exp(-oldbeta))
    new.pr <- 1 / (1 + exp(-mytheta$allpara[1:mpq]))
    eps <- sum(abs(old.pr - new.pr) / old.pr)
    #print(eps)
  }
  
  ######
  ########
  ###############
  ################
  
  
  
  
  
  list(
    mytheta = mytheta,
    mypsi = mypsi,
    initial.para = initial.para,
    eps = eps
  )
}


####### FDR, FNR, F1
feature_selection=function(selected, beta.true){
true.active <- which(beta.true != 0)
true.null <- which(beta.true == 0)

TP <- length(intersect(selected, true.active))
FP <- length(intersect(selected, true.null))
FN <- length(setdiff(true.active, selected))

FDR <- ifelse(length(selected)==0,
              0,
              FP/length(selected))

FNR <- FN/length(true.active)

F1 <- ifelse((2*TP+FP+FN)==0,
             0,
             2*TP/(2*TP+FP+FN))
return(list(FDR=FDR, FNR=FNR, F1=F1))
}



