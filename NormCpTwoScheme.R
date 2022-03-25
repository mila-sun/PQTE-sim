library(EnvStats)
library(quantreg)
library(betareg)
library(copula)
library(gamCopula)
library(VineCopula)
library(foreach)
library(doParallel)

#ncores = Sys.getenv("SLURM_CPUS_PER_TASK")
ncores = 1
registerDoParallel(cores=ncores)
print(ncores)

# Data generation from two simulation schemes
sim.ds = function(n, phi, rho.param, fam){
  
  # simulate covariates
  x1 = rbinom(n, size=1, prob = 0.6)
  x1 = c(x1,x1)
  x2 = rbinom(n, size=1, prob = 0.4)
  x2 = c(x2,x2)
  u = c(rnorm(n), (1+rnorm(n)))
  z = rep(c(1,0), each=n)
  id = rep(c(1:n), times=2)
  
  # simulate S0 and S1 from copula with family=fam
  Us = condBiCopSim(family = fam, 
                    calib.fnc = function(x1,x2){
                      cbind(1,x1,x2)%*%rho.param},
                    X = as.matrix(cbind(x1,x2))[1:n,], 
                    return.par = TRUE, tau = T)
  Us.dt = Us$data
  colnames(Us.dt) = c("u1", "u2")
  
  # Beta distribution parameter 
  mu.s0 = 1/(1+exp(0.1+0.25*x1-0.35*x2+0.15*u))
  mu.s1 = 1/(1+exp(-0.2+0.3*x1+0.2*x2+0.2*u))
  
  a0 = mu.s0*phi[1]
  b0 = (1-mu.s0)*phi[1]
  a1 = mu.s1*phi[2]
  b1 = (1-mu.s1)*phi[2]
  
  # error term 
  e = rpareto(n=2*n, location = 1, shape = 3)
  
  #---- (1) comonotone scheme ----
  Us.dt.com = rbind(Us.dt, Us.dt)
  
  s0.com = qbeta(pobs(Us.dt.com[,1]), shape1 = a0, shape2 = b0)
  s1.com = qbeta(pobs(Us.dt.com[,2]), shape1 = a1, shape2 = b1)
  
  # simulate outcome Y
  y.com = 1+z+z*s0.com+z*s1.com+s0.com+s1.com+x1+x2+e
  
  # create data frame
  ds.com = data.frame(id=id, y=y.com, x1=x1, x2=x2, u=u, z=z, s0=s0.com, s1=s1.com, a0=a0, b0=b0, a1=a1, b1=b1)
  
  #---- (1) independent scheme ----
  Us.dt.indep = cbind(Us.dt, Us$par)
  
  # simulate S0 and S1 using conditional copula
  u2.prime = apply(Us.dt.indep, 1, function(x){
    cCopula(cbind(x[1], runif(1)), indices = 2,
            copula = normalCopula(param=x[3], dim=2), 
            inverse = T)})
  u1.prime = apply(Us.dt.indep, 1, function(x){
    cCopula(cbind(x[2], runif(1)), indices = 2,
            copula = normalCopula(param=x[3], dim=2), 
            inverse = T)})
  U = data.frame(u1=c(u1.prime,Us.dt[,1]), u2=c(Us.dt[,2],u2.prime))
  
  s0.indep = qbeta(pobs(U[,1]), shape1 = a0, shape2 = b0)
  s1.indep = qbeta(pobs(U[,2]), shape1 = a1, shape2 = b1)
  
  # simulate outcome Y
  y.indep = 1+z+z*s0.indep+z*s1.indep+s0.indep+s1.indep+x1+x2+e
  
  # create data frame
  ds.indep = data.frame(id=id, y=y.indep, x1=x1, x2=x2, u=u, z=z, s0=s0.indep, s1=s1.indep, a0=a0, b0=b0, a1=a1, b1=b1)
  
  return(list(ds.com, ds.indep))
}

# Copula imputation method
imp.func.cp = function(ds, fam){
  
  # marginal model
  s0.brg = betareg(s0~x1+x2+u, data = ds, link = "logit", type="BC", na.action = na.omit)
  s1.brg = betareg(s1~x1+x2+u, data = ds, link = "logit", type="BC", na.action = na.omit)
  
  # fitted and predicted Beta distribution parameters
  hat.a0.fitted = s0.brg$fitted.values*(s0.brg$coefficients$precision)
  hat.b0.fitted = (1-s0.brg$fitted.values)*(s0.brg$coefficients$precision)
  hat.a0.pred = predict(s0.brg, type="response", newdata=ds[ds$z==1,c("x1","x2","u")])*(s0.brg$coefficients$precision)
  hat.b0.pred = (1-predict(s0.brg, type="response", newdata=ds[ds$z==1,c("x1","x2","u")]))*(s0.brg$coefficients$precision)
  
  hat.a1.fitted = s1.brg$fitted.values*(s1.brg$coefficients$precision)
  hat.b1.fitted = (1-s1.brg$fitted.values)*(s1.brg$coefficients$precision)
  hat.a1.pred = predict(s1.brg, type="response", newdata=ds[ds$z==0,c("x1","x2","u")])*(s1.brg$coefficients$precision)
  hat.b1.pred = (1-predict(s1.brg, type="response", newdata=ds[ds$z==0,c("x1","x2","u")]))*(s1.brg$coefficients$precision)
  
  
  # get CDF U0 and U1
  hat.u0 = pbeta(ds$s0[ds$z==0], shape1 = hat.a0.fitted, shape2 = hat.b0.fitted)
  hat.u1 = pbeta(ds$s1[ds$z==1], shape1 = hat.a1.fitted, shape2 = hat.b1.fitted)
  
  # use the rank of hat.mu r
  rank.u1 = pobs(hat.u1)
  rank.u0 = pobs(hat.u0)
  
  # copula model
  n = nrow(ds)/2
  ds.cp = data.frame(u1=rank.u1, u2=rank.u0, x1=ds$x1[1:n], x2=ds$x2[1:n])
  fit.cp = gamBiCopFit(ds.cp, formula=~1+x1+x2, family = fam, tau = T)$res
  hat.tau = gamBiCopPredict(fit.cp, target = "tau")$tau
  hat.par = BiCopTau2Par(family = fam, tau = hat.tau)
  
  grid.cp = cbind(hat.par, rank.u0, rank.u1)
  
  #use the empirical mean
  u1.imp = apply(grid.cp, 1, function(x){
    mean(cCopula(cbind(x[2], runif(1000)), indices = 2,
                 copula = normalCopula(param=x[1], dim=2), 
                 inverse = T))})
  u0.imp = apply(grid.cp, 1, function(x){
    mean(cCopula(cbind(x[3], runif(1000)), indices = 2,
                 copula = normalCopula(param=x[1], dim=2), 
                 inverse = T))})
  s1.imp = qbeta(u1.imp, shape1 = hat.a1.pred, shape2 = hat.b1.pred)
  s0.imp = qbeta(u0.imp, shape1 = hat.a0.pred, shape2 = hat.b0.pred)
  
  
  # imputed data
  ds.imp  = ds
  ds.imp$s0 = ifelse(ds.imp$z==1, s0.imp, ds.imp$s0)
  ds.imp$s1 = ifelse(ds.imp$z==0, s1.imp, ds.imp$s1)
  
  return(ds.imp)
  
}

# Comonotonic imputation method
imp.func.com = function(ds, fam){
  
  # marginal model
  s0.brg = betareg(s0~x1+x2+u, data = ds, link = "logit", type="BC", na.action = na.omit)
  s1.brg = betareg(s1~x1+x2+u, data = ds, link = "logit", type="BC", na.action = na.omit)
  
  # fitted and predicted Beta distribution parameters
  hat.a0.fitted = s0.brg$fitted.values*(s0.brg$coefficients$precision)
  hat.b0.fitted = (1-s0.brg$fitted.values)*(s0.brg$coefficients$precision)
  hat.a0.pred = predict(s0.brg, type="response", newdata=ds[ds$z==1,c("x1","x2","u")])*(s0.brg$coefficients$precision)
  hat.b0.pred = (1-predict(s0.brg, type="response", newdata=ds[ds$z==1,c("x1","x2","u")]))*(s0.brg$coefficients$precision)
  
  hat.a1.fitted = s1.brg$fitted.values*(s1.brg$coefficients$precision)
  hat.b1.fitted = (1-s1.brg$fitted.values)*(s1.brg$coefficients$precision)
  hat.a1.pred = predict(s1.brg, type="response", newdata=ds[ds$z==0,c("x1","x2","u")])*(s1.brg$coefficients$precision)
  hat.b1.pred = (1-predict(s1.brg, type="response", newdata=ds[ds$z==0,c("x1","x2","u")]))*(s1.brg$coefficients$precision)
  
  # assume that u0 when Z=0 equals to u0 when Z=1
  # get CDF U0 and U1
  hat.u0 = pbeta(ds$s0[ds$z==0], shape1 = hat.a0.fitted, shape2 = hat.b0.fitted)
  hat.u1 = pbeta(ds$s1[ds$z==1], shape1 = hat.a1.fitted, shape2 = hat.b1.fitted)
  
  # impute S(0) when Z=1 using marginal beta distribution 
  s0.imp = qbeta(hat.u0, shape1 = hat.a0.pred, shape2 = hat.b0.pred)
  s1.imp = qbeta(hat.u1, shape1 = hat.a1.pred, shape2 = hat.b1.pred)
  
  
  # imputed data
  ds.imp = ds
  ds.imp$s0 = ifelse(ds.imp$z==1, s0.imp, ds.imp$s0)
  ds.imp$s1 = ifelse(ds.imp$z==0, s1.imp, ds.imp$s1)
  
  return(ds.imp)
  
}

boot.func = function(D.I, B, tau, fam){
  
  coef.imp.boot = foreach(icount(B), .combine = rbind, .export=c('imp.func.cp','imp.func.com','D.I','tau'),.packages = c("gamCopula","VineCopula","copula","quantreg")) %dopar% {
    D.star.id = sample(c(1:max(D.I$id)), size = 0.7*max(D.I$id), replace = T)
    D.star = D.I[c(D.star.id,(D.star.id+max(D.I$id))),]
    
    # set missing values
    D.star$s0 = ifelse(D.star$z==1, NA, D.star$s0)
    D.star$s1 = ifelse(D.star$z==0, NA, D.star$s1)
    
    # impute missing values (copula method)
    D.I.star.cp = imp.func.cp(D.star, fam = fam)
    
    # impute missing values (comonotone method)
    D.I.star.com = imp.func.com(D.star, fam = fam)
    
    # outcome model
    fit.imp.boot.cp = summary(rq(y~z*s0+z*s1+x1+x2, data = D.I.star.cp, tau=tau))$coef
    fit.imp.boot.com = summary(rq(y~z*s0+z*s1+x1+x2, data = D.I.star.com, tau=tau))$coef
    c(as.numeric(fit.imp.boot.cp[c(2,7,8),1]),as.numeric(fit.imp.boot.com[c(2,7,8),1]))
  }
  return(coef.imp.boot)
}

model.func = function(ds.com, ds.indep, fam, tau, B){
  
  ds.null.com = ds.com
  ds.com$s0 = ifelse(ds.com$z==1, NA, ds.com$s0)
  ds.com$s1 = ifelse(ds.com$z==0, NA, ds.com$s1)
  
  ds.null.indep = ds.indep
  ds.indep$s0 = ifelse(ds.indep$z==1, NA, ds.indep$s0)
  ds.indep$s1 = ifelse(ds.indep$z==0, NA, ds.indep$s1)
  
  #---- comonotonic simulation scheme ----
  ds.com.cp = imp.func.cp(ds.com, fam = fam)
  ds.com.com = imp.func.com(ds.com, fam=fam)
  
  fit.com.cp = summary(rq(y~z*s0+z*s1+x1+x2, data = ds.com.cp, tau=tau))$coef
  fit.com.com = summary(rq(y~z*s0+z*s1+x1+x2, data = ds.com.com, tau=tau))$coef
  fit.com.true = summary(rq(y~z*s0+z*s1+x1+x2, data = ds.null.com, tau=tau))$coef
  
  # bootstrap variance
  com.coef.boot = boot.func(D.I=ds.com, B, tau = tau, fam = fam)
  
  #---- independent simulation scheme ----
  ds.indep.cp = imp.func.cp(ds.indep, fam = fam)
  ds.indep.com = imp.func.com(ds.indep, fam=fam)
  
  fit.indep.cp = summary(rq(y~z*s0+z*s1+x1+x2, data = ds.indep.cp, tau=tau))$coef
  fit.indep.com = summary(rq(y~z*s0+z*s1+x1+x2, data = ds.indep.com, tau=tau))$coef
  fit.indep.true = summary(rq(y~z*s0+z*s1+x1+x2, data = ds.null.indep, tau=tau))$coef
  
  # bootstrap variance
  indep.coef.boot = boot.func(D.I=ds.indep, B, tau = tau, fam = fam)
  
  # combine results
  com.coef.est = c(as.numeric(fit.com.true[c(2,7,8),1]), as.numeric(fit.com.cp[c(2,7,8),1]), as.numeric(fit.com.com[c(2,7,8),1]), apply(com.coef.boot, 2, mean))
  com.coef.se = c(as.numeric(fit.com.true[c(2,7,8),2]), as.numeric(fit.com.cp[c(2,7,8),2]),  as.numeric(fit.com.com[c(2,7,8),2]),apply(com.coef.boot, 2, sd))
  
  indep.coef.est = c(as.numeric(fit.indep.true[c(2,7,8),1]), as.numeric(fit.indep.cp[c(2,7,8),1]), as.numeric(fit.indep.com[c(2,7,8),1]), apply(indep.coef.boot, 2, mean))
  indep.coef.se = c(as.numeric(fit.indep.true[c(2,7,8),2]), as.numeric(fit.indep.cp[c(2,7,8),2]), as.numeric(fit.indep.com[c(2,7,8),2]),apply(indep.coef.boot, 2, sd))
  
  return(list(com.coef.est, com.coef.se, indep.coef.est, indep.coef.se))
}

#m=100
m = 1000
n = 3000
fam = 1
tau = 0.9
B = 1000
rho.param=c(1.25,1.5,1)
phi=c(10,5)

com.coef.est = com.coef.sd = matrix(NA, nrow = m, ncol = 21)
indep.coef.est = indep.coef.sd = matrix(NA, nrow = m, ncol = 21)

for(i in 1:m){
  ds.temp = sim.ds(n, phi=phi, rho.param=rho.param, fam=fam)
  model.temp = model.func(ds.com=ds.temp[[1]], ds.indep=ds.temp[[2]], fam, tau, B)
  com.coef.est[i,] = model.temp[[1]]
  com.coef.sd[i,] = model.temp[[2]]
  indep.coef.est[i,] = model.temp[[3]]
  indep.coef.sd[i,] = model.temp[[4]]
}

col_name = c("z_true","zs0_true","zs1_true","z_imp_cp","zs0_imp_cp","zs1_imp_cp","z_imp_com","zs0_imp_com","zs1_imp_com","z_boot_cp","zs0_boot_cp","zs1_boot_cp","z_boot_com","zs0_boot_com","zs1_boot_com")
com.coef.est = as.data.frame(com.coef.est)
colnames(com.coef.est) = col_name

com.coef.sd = as.data.frame(com.coef.sd)
colnames(com.coef.sd) = col_name

indep.coef.est = as.data.frame(indep.coef.est)
colnames(indep.coef.est) = col_name

indep.coef.sd = as.data.frame(indep.coef.sd)
colnames(indep.coef.sd) = col_name

#args = commandArgs(trailingOnly = TRUE)
stopImplicitCluster()
