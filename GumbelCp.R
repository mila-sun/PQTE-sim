library(car)
library(EnvStats)
library(quantreg)
library(betareg)
library(copula)
library(gamCopula)
library(VineCopula)
library(foreach)
library(doParallel)
library(doSNOW)

gFunc_deriv = function(t, theta){
  -(1/theta)*(t^(1/theta-1))*exp(-t^(1/theta))
}
f_u = function(u,theta,v_y,z_i){
  gFunc_deriv((-log(u))^(theta)+v_y, theta)-z_i*gFunc_deriv(v_y, theta)
}


sim.ds = function(n, phi, rho.param){
  
  # simulate covariates
  x1 = rbinom(n, size=1, prob = 0.6)
  x1 = c(x1,x1)
  x2 = rbinom(n, size=1, prob = 0.4)
  x2 = c(x2,x2)
  u = c(rnorm(n), (1+rnorm(n)))
  z = rep(c(1,0), each=n)
  id = rep(c(1:n), times=2)
  
  # simulate S0 and S1 from copula with family=fam
  Us = condBiCopSim(family = 401, 
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
  
  #---- (1) independent scheme ----
  Us.dt.indep = cbind(Us.dt, Us$par)
  
  # simulate S0 and S1 using conditional copula
  u.imp = apply(Us.dt.indep, 1, function(x){
    temp.u2.prime = uniroot(f_u,interval = c(0,1),theta=x[3], 
                    v_y=(-log(x[1]))^(x[3]),z_i=runif(1))$root
    temp.u1.prime = uniroot(f_u,interval = c(0,1),theta=x[3], 
              v_y=(-log(x[2]))^(x[3]),z_i=runif(1))$root
    c(temp.u2.prime, temp.u1.prime)
  })
  
  U = data.frame(u1=c(u.imp[2,],Us.dt[,1]), u2=c(Us.dt[,2],u.imp[1,]))
  
  s0.indep = qbeta(pobs(U[,1]), shape1 = a0, shape2 = b0)
  s1.indep = qbeta(pobs(U[,2]), shape1 = a1, shape2 = b1)
  
  
  # simulate outcome Y
  y.indep = 1+z+z*s0.indep+z*s1.indep+s0.indep+s1.indep+x1+x2+e
  
  # create data frame
  ds.indep = data.frame(id=id, y=y.indep, x1=x1, x2=x2, u=u, z=z, s0=s0.indep, s1=s1.indep, a0=a0, b0=b0, a1=a1, b1=b1)
  
  return(ds.indep)
}

imp.func = function(ds){
  
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
  
  # use the rank of hat.mu (IME) paper
  rank.u1 = pobs(hat.u1)
  rank.u0 = pobs(hat.u0)
  
  # copula model
  n = nrow(ds)/2
  ds.cp = data.frame(u1=rank.u1, u2=rank.u0, x1=ds$x1[1:n], x2=ds$x2[1:n])
  fit.cp = gamBiCopFit(ds.cp, formula=~1+x1+x2, family = 401, tau = T)$res
  hat.tau = gamBiCopPredict(fit.cp, target = "tau")$tau
  hat.par = BiCopTau2Par(family = 4, tau = hat.tau)
  
  grid.cp = data.frame(cbind(hat.par, rank.u0, rank.u1))
  
  # copula imputation using the empirical mean
  gFunc_deriv = function(t, theta){-(1/theta)*(t^(1/theta-1))*exp(-t^(1/theta))}
  f_u = function(u,theta,v_y,z_i){
    gFunc_deriv((-log(u))^(theta)+v_y, theta)-z_i*gFunc_deriv(v_y, theta)
  }
  
  # trace the seed
  initial_seed = as.integer(Sys.time())
  the_seed = initial_seed %% 100000
  set.seed(the_seed)
  
  unif.rm = runif(500)
  
  u.imp = apply(grid.cp, 1, function(x){
    temp1 = mean(sapply(unif.rm, function(y){
      uniroot(f_u,interval = c(0,1),
              theta=x[1], 
              v_y=(-log(x[2]))^(x[1]),
              z_i=y)$root}))
    temp2 = mean(sapply(unif.rm, function(y){
      uniroot(f_u,interval = c(0,1),
              theta=x[1], 
              v_y=(-log(x[3]))^(x[1]),
              z_i=y)$root}))
    c(temp1, temp2)
  })
  
  s1.imp = qbeta(u.imp[1,], shape1 = hat.a1.pred, shape2 = hat.b1.pred)
  s0.imp = qbeta(u.imp[2,], shape1 = hat.a0.pred, shape2 = hat.b0.pred)
  
  
  # imputed data
  ds.imp = ds
  ds.imp$s0 = ifelse(ds.imp$z==1, s0.imp, ds.imp$s0)
  ds.imp$s1 = ifelse(ds.imp$z==0, s1.imp, ds.imp$s1)
  
  #cat(" imp_func runif seed =", the_seed, "\n")
  return(ds.imp)
  
}

boot.func = function(D.I, B, tau){
  #cat("Bootstrap starts: \n")
  # parallel computing
  coef.imp.boot = foreach(icount(B), .combine=rbind, .export=c('imp.func','D.I','tau'),.packages = c("gamCopula","VineCopula","copula","quantreg")) %dopar% {
    
    # trace seed
    the_seed_bt = as.integer(Sys.time()) %% 10000
    set.seed(the_seed_bt)
    #cat("boostrap sample seed =", the_seed_bt,";")
    D.star.id = sample(c(1:max(D.I$id)), size = max(D.I$id), replace = T)
    D.star = D.I[c(D.star.id,(D.star.id+max(D.I$id))),]
    
    # impute missing values
    D.I.star = imp.func(D.star)
    
    # outcome model
    fit.imp.boot = summary(rq(y~z*s0+z*s1+x1+x2, data = D.I.star, tau=tau))$coef
    as.numeric(fit.imp.boot[c(2,7,8),1])
  }

  return(coef.imp.boot)
}

model.func = function(ds, tau, B){
  
  ds.null = ds
  ds$s0 = ifelse(ds$z==1, NA, ds$s0)
  ds$s1 = ifelse(ds$z==0, NA, ds$s1)
  
  ds.imp = imp.func(ds)
  
  # outcome model
  fit.imp = summary(rq(y~z*s0+z*s1+x1+x2, data = ds.imp, tau=tau))$coef
  fit.true = summary(rq(y~z*s0+z*s1+x1+x2, data = ds.null, tau=tau))$coef
  
  # bootstrap variance
  coef.boot = boot.func(D.I=ds, B, tau = tau)
  stopImplicitCluster()
  
  # combine results
  coef.est = c(as.numeric(fit.true[c(2,7,8),1]), as.numeric(fit.imp[c(2,7,8),1]), apply(coef.boot, 2, mean))
  coef.se = c(as.numeric(fit.true[c(2,7,8),2]), as.numeric(fit.imp[c(2,7,8),2]), apply(coef.boot, 2, sd))
  
  return(list(coef.est, coef.se))
}

m = 1
n = 3000
tau = 0.8
B = 500
rho.param=c(1.25,1.5,1)
phi=c(10,5)

coef.est = coef.sd = matrix(NA, nrow = m, ncol = 9)

ncores = Sys.getenv("SLURM_CPUS_PER_TASK")
registerDoParallel(cores=ncores)
print(ncores)
getDoParWorkers()

for(i in 1:m){
  ds.temp = sim.ds(n, phi=phi, rho.param=rho.param)
  model.temp = model.func(ds=ds.temp, tau, B)
  coef.est[i,] = model.temp[[1]]
  coef.sd[i,] = model.temp[[2]]
}

coef.est = as.data.frame(coef.est)
colnames(coef.est) = c("z_true","zs0_true","zs1_true","z_imp","zs0_imp","zs1_imp","z_boot","zs0_boot","zs1_boot")

coef.sd = as.data.frame(coef.sd)
colnames(coef.sd) = c("z_true","zs0_true","zs1_true","z_imp","zs0_imp","zs1_imp","z_boot","zs0_boot","zs1_boot")
