# ==============================================================================
# Data simulation
# ===============================================================================
source("source_code_HS.R")
nnrep=1
store.our=store.lasso=rep(0, nnrep)
#for( numrep in 1:nnrep){
numrep=100 
set.seed(numrep)
print(paste("Running seed", numrep))

n <- 600
p <- 1000
s <- 0.05
p2 <- s * p
p1 <- floor(p2 / 2)

beta.true <- rep(0, p)
beta.true[1:p1] <- seq(0.2, 0.69, length = p1)
beta.true[(1 + p1):p2] <- -seq(0.2, 0.69, length = (p2 - p1))

gamma <- c(-0.5, 0.5)
true.idx <- which(beta.true != 0)

sim.out <- simulate_survival_data(
  n = n,
  p = p,
  s = s,
  beta.true = beta.true,
  gamma = gamma,
  plot_km = TRUE
)

mydata <- sim.out$mydata
mean_delta <- sim.out$censoring_indicator_mean
#print(mean_delta)

# ==============================================================================
#  Estimation using different methods
# ===============================================================================

# ------------------------------------------------------------------------------
# Method 1: Cox LASSO initialization
# ------------------------------------------------------------------------------

lasso.out <- lasso_initialization(mydata = mydata, p = p)
betahatlasso <- lasso.out$betahatlasso
beta.est <- lasso.out$beta.est

# ------------------------------------------------------------------------------
# Method 2: Proposed iterative estimation
# ------------------------------------------------------------------------------

npart <- 10
prepare_global_objects(mydata = mydata, npart = npart, varz_value = 2)

proposed.out <- run_proposed_estimation(
  mydata = mydata,
  beta.est = beta.est,
  npart = npart,
  outer_tol = 0.005,
  initial_tau02 = 2,
  update_tau02 = 2
)

mytheta <- proposed.out$mytheta
mypsi <- proposed.out$mypsi
initial.para <- proposed.out$initial.para

# ==============================================================================
# Basic diagnostic plots
# ===============================================================================

plot(mytheta$allpara[1:mpq], main = "Proposed estimates", ylab = "Estimate")
plot(beta.est, main = "LASSO initial estimates", ylab = "Estimate")
plot(beta.true, main = "True values of the regression coefficients", ylab = "Estimate")
store.our[numrep]=sqrt(mean((mytheta$allpara[1:p]-beta.true)^2))
store.lasso[numrep]=sqrt(mean((beta.est[1:p]-beta.true)^2))

# ==========================================================
# Estimation accuracy 
# ==========================================================

cat('RMSE for HS= ', sqrt(mean((mytheta$allpara[1:p]-beta.true)^2)), "\n")
cat('RMSE for LASSO= ', sqrt(mean((beta.est[1:p]-beta.true)^2)), "\n")


#==========================================================
# Feature selection metrics
#==========================================================
lasso.selection=feature_selection(which(beta.est!=0), beta.true)
cat("Feature selection metrics using LASSO:",
    "FDR =", lasso.selection$FDR,
    ", FNR =", lasso.selection$FNR,
    ", F1 =", lasso.selection$F1,
    "\n")
###########
###########
anewdata <- data.frame(V=mydata$V,lambda=mytheta$lambda)
anewdata=anewdata[order(anewdata$V), ]
sigmahat=(sum(c(anewdata$V[1:n]-c(0, anewdata$V[1:(n-1)]))*anewdata$lambda))^0.5
threshold=sigmahat*sqrt(log(p)/sum(mydata$delta));
selected=which(abs(mytheta$allpara[1:p])>threshold)
our.selection.HS1=feature_selection(selected, beta.true)
cat("Feature selection metrics using HS1:",
    "FDR =", our.selection.HS1$FDR,
    ", FNR =", our.selection.HS1$FNR,
    ", F1 =", our.selection.HS1$F1, "\n")
############
theta.abs <- abs(mytheta$allpara[1:p])
set.seed(123)
km <- kmeans(theta.abs, centers = 2, nstart = 50)
cluster.means <- tapply(theta.abs, km$cluster, mean)
zero.cluster <- as.numeric(names(which.min(cluster.means)))
nonzero.cluster <- as.numeric(names(which.max(cluster.means)))
selected <- which(km$cluster == nonzero.cluster)

our.selection.HS2 <- feature_selection(selected, beta.true)
cat("Feature selection metrics using HS2:",
    "FDR =", our.selection.HS2$FDR,
    ", FNR =", our.selection.HS2$FNR,
    ", F1 =", our.selection.HS2$F1, "\n")



