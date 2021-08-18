setwd("~/Desktop/Duke/00 Thesis/01 bridging study/R code")
# MCT simulation
library(dplyr)
library(ggplot2)
library(nlme)
library(tidyr)
library(viridis)
library(multcompView)

########
# TOST #
########
test.tost = function(dat, MSE, n, alpha) {
  k = 6 # sequence
  j = 3 # period
  # BE limits
  lr = 0.8
  up = 1.25
  # contrasts
  C.TA = c(-4,-2,6,5,-2,-3,-1,4,-3,4,2,-6,1,-4,3,-5,2,3)/24
  C.TB = c(-5,2,3,1,-4,3,4,2,-6,-1,4,-3,5,-2,-3,-4,-2,6)/24
  C.AB = c(-1,4,-3,-4,-2,6,5,-2,-3,-5,2,3,4,2,-6,1,-4,3)/24
  sum(C.TA)
  sum(C.TB)
  sum(C.AB)
  S = MSE # intra-subject MSE with DF 2N-6 << N=6*n
  # unbiased estimator for thetas (difference between treatment)
  logmu.T = log(mean(dat[dat$trt=='T', ]$response))
  logmu.A = log(mean(dat[dat$trt=='A', ]$response))
  logmu.B = log(mean(dat[dat$trt=='B', ]$response))
  # theta.hat.ta = 1/(k*n)*(log(sum(dat[dat$trt=='T', ]$response))-log(sum(dat[dat$trt=='A', ]$response)))
  # theta.hat.tb = 1/(k*n)*(log(sum(dat[dat$trt=='T', ]$response))-log(sum(dat[dat$trt=='B', ]$response)))
  # theta.hat.ab = 1/(k*n)*(log(sum(dat[dat$trt=='A', ]$response))-log(sum(dat[dat$trt=='B', ]$response)))
  
  theta.hat.ta = logmu.T - logmu.A
  theta.hat.tb = logmu.T - logmu.B
  theta.hat.ab = logmu.A - logmu.B
  # test statistics
  t.crit = qt(alpha, 2*n*k-6)
  TL.ta = (theta.hat.ta - log(lr))/(S*sqrt(sum(C.TA^2)/n))
  TU.ta = (theta.hat.ta - log(up))/(S*sqrt(sum(C.TA^2)/n))
  TL.ta > t.crit
  TU.ta < -t.crit
  TL.tb = (theta.hat.tb - log(lr))/(S*sqrt(sum(C.TB^2)/n))
  TU.tb = (theta.hat.tb - log(up))/(S*sqrt(sum(C.TB^2)/n))
  TL.tb > t.crit
  TU.tb < -t.crit
  TL.ab = (theta.hat.ab - log(lr))/(S*sqrt(sum(C.AB^2)/n))
  TU.ab = (theta.hat.ab - log(up))/(S*sqrt(sum(C.AB^2)/n))
  TL.ab > t.crit
  TU.ab < -t.crit
  
  if (TL.ta > t.crit & TU.ta < -t.crit & TL.tb > t.crit & TU.tb < -t.crit & TL.ab > t.crit & TU.ab < -t.crit) {
    # if all true, reject the bull hypothesis of NO BE and conclude BE
    return(1)
  } else {
    return(0)
  }
}

#######################
# Confidence Interval #
#######################
test.ci = function(dat, MSE, n, alpha) {
  k = 6 # sequence
  j = 3 # period
  # BE limits
  lr = 0.8
  up = 1.25
  # contrasts
  C.TA = c(-4,-2,6,5,-2,-3,-1,4,-3,4,2,-6,1,-4,3,-5,2,3)/24
  C.TB = c(-5,2,3,1,-4,3,4,2,-6,-1,4,-3,5,-2,-3,-4,-2,6)/24
  C.AB = c(-1,4,-3,-4,-2,6,5,-2,-3,-5,2,3,4,2,-6,1,-4,3)/24
  sum(C.TA)
  sum(C.TB)
  sum(C.AB)
  S = MSE # intra-subject MSE with DF 2N-6 << N=6*n
  # unbiased estimator for thetas (difference between treatment)
  logmu.T = log(mean(dat[dat$trt=='T', ]$response))
  logmu.A = log(mean(dat[dat$trt=='A', ]$response))
  logmu.B = log(mean(dat[dat$trt=='B', ]$response))

  theta.hat.ta = logmu.T - logmu.A
  theta.hat.tb = logmu.T - logmu.B
  theta.hat.ab = logmu.A - logmu.B
  
  t.crit = qt(alpha, 2*n*k-6)
  # T-A
  L.ta = theta.hat.ta - t.crit * (S*sqrt(sum(C.TA^2)/n))
  U.ta = theta.hat.ta + t.crit * (S*sqrt(sum(C.TA^2)/n))
  # T-B
  L.tb = theta.hat.tb - t.crit * (S*sqrt(sum(C.TB^2)/n))
  U.tb = theta.hat.tb + t.crit * (S*sqrt(sum(C.TB^2)/n))
  # A-B
  L.ab = theta.hat.ab - t.crit * (S*sqrt(sum(C.AB^2)/n))
  U.ab = theta.hat.ab + t.crit * (S*sqrt(sum(C.AB^2)/n))
  
  if (L.ta > log(lr) & U.ta < log(up) & L.tb > log(lr) & U.tb < log(up) & L.ab > log(lr) & U.ab < log(up)) {
    return(1)
  } else {
    return(0)
  }
}

############
# Pairwise #
############
test.tukey = function(dat, alpha) {
  dat = dat %>% filter(prd==1)
  model=lm(response ~ trt, data=dat)
  ANOVA=aov(model)
  TUKEY <- TukeyHSD(x=ANOVA, 'trt', conf.level=1-2*alpha)
  if (all(TUKEY$trt[,4]>(1-(1-2*alpha)^(1/3)))) {
    return(1)
  } else {
    return(0)
  }
}

test.t = function(dat, alpha) {
  dat = dat %>% filter(prd==1)
  dat.T = dat %>% filter(trt=='T') %>% select(response)
  dat.A = dat %>% filter(trt=='A') %>% select(response)
  dat.B = dat %>% filter(trt=='B') %>% select(response)
  t1 = t.test(dat.T, dat.A, alternative = "two.sided", var.equal = TRUE, conf.level = 1-alpha)
  t1.p = t1$p.value
  t2 = t.test(dat.T, dat.B, alternative = "two.sided", var.equal = TRUE, conf.level = 1-alpha)
  t2.p = t2$p.value
  t3 = t.test(dat.A, dat.A, alternative = "two.sided", var.equal = TRUE, conf.level = 1-alpha)
  t3.p = t3$p.value
  t.p = c(t1.p, t2.p, t3.p)
  if (all(t.p>alpha)) {
    return(1)
  } else {
    return(0)
  }
}

##############
# Simulation #
##############
run.sim = function(n=10, mu.T=0.5, mu.A=0.5, mu.B=0.5, CV=0.1, n.sim=100, alpha=0.1, seed=100) {
  set.seed(seed)
  k = 6
  j = 3
  N = n * k * j
  sd.T = CV * mu.T
  sd.A = CV * mu.A
  sd.B = CV * mu.B
  
  sim.dat = data.frame(sim=rep(seq(n.sim), each=N), 
                       id=rep(rep(seq(n*k),each=j),n.sim),
                       sqc=rep(rep(seq(k), each=n*j),n.sim),
                       prd=rep(seq(j),k*n*n.sim)
  )
  trt = rep(c(rep(c('T','B','A'), n), rep(c('A','T','B'), n), rep(c('B','A','T'), n),
              rep(c('A','B','T'), n), rep(c('B','T','A'), n), rep(c('T','A','B'), n)), n.sim)
  sim.dat$trt = trt
  sim.dat$prd = as.factor(sim.dat$prd)
  sim.dat$trt = factor(sim.dat$trt, levels = c('T','A','B'))
  sim.dat$response = 0
  sim.dat[sim.dat$trt=='T', 'response'] = rnorm(k*n*n.sim, mean=mu.T, sd=sd.T)
  sim.dat[sim.dat$trt=='A', 'response'] = rnorm(k*n*n.sim, mean=mu.A, sd=sd.A)
  sim.dat[sim.dat$trt=='B', 'response'] = rnorm(k*n*n.sim, mean=mu.B, sd=sd.B)
  
  sim.result = data.frame(n.sim=integer(), MSE=double(), p=double())
  
  p = data.frame(tost=double(), ci=double(), pairwise=double())
  sim.record = data.frame(i = integer(), ci.lb=double(), ci.ub=double())
  for (i in 1:n.sim) {
    dat = sim.dat[sim.dat$sim==i,]
    # fit linear mixed effects model
    tryCatch({
      # model = lme(data=dat, fixed=response~sqc+prd+trt, random=~1|id)
      # MSE = mean(model$residuals^2)
      # MSE = sqrt(var.e)
      aov.out = aov(response~sqc+trt+prd + Error(id/sqc+trt+prd), data=dat)
      MSE = mean(aov.out$Within$residuals^2)
      # sim.record =  rbind(sim.record, data.frame(i=i, ci.lb=ci.lb, ci.ub=ci.ub))
      p.result.tost = test.tost(dat, MSE, n, alpha)
      p.result.ci = test.ci(dat, MSE, n, alpha)
      p.result.t = test.t(dat, alpha)
      
      p = rbind(p, data.frame(tost=p.result.tost, ci=p.result.ci, pairwise=p.result.t))
    }, error=function(e) {
      p = p
      print(c(i, e))}
    )
  }
  pwr = apply(p, 2, sum)/n.sim
  return(pwr)
}

q90 = function(col) quantile(col, c(0.05,0.95)) 

if (dir.exists('output')==F) {
  dir.create('output')
}

##############
# scenario 1 #
##############
# BE between EU and US, test for BE of proposed biosimilar
exp = expand.grid(n=c(1,2,3,4,5,10), 
                  mu.T=c(0.395, 0.40, 0.45, 0.5, 0.5625, 0.625, 0.63), 
                  mu.A=0.5,
                  mu.B=0.5, 
                  CV=c(0.1,0.2,0.3,0.4))
sim.results = data.frame(pwr.tost=double(), pwr.ci=double(), pwr.pairwise=double())
start.time = Sys.time()
for (i in 1:nrow(exp)) {
  n = exp$n[i]
  mu.T = exp$mu.T[i]
  mu.A = exp$mu.A[i]
  mu.B = exp$mu.B[i]
  CV = exp$CV[i]
  print(paste('progress: ', i, round(100*i/nrow(exp), 2), '% time lapsed: ', round(Sys.time()-start.time, 2), 'sec'))
  pwr = run.sim(n=n, mu.T=mu.T, mu.A=mu.A, mu.B=mu.B, CV=CV, n.sim=1000, alpha=0.05, seed=100)
  sim.results = rbind(sim.results, data.frame(pwr.tost=pwr[[1]], pwr.ci=pwr[[2]], pwr.pairwise=pwr[[3]]))
}

sim.results1 = cbind(exp, sim.results)
write.csv(sim.results1, paste('output/simulation 1 ', Sys.time(), '.csv', sep=''))
write.csv(sim.results1, 'output/simulation 1.csv')
# sim.results1 = read.csv('simulation 1.csv')

##############
# scenario 2 #
##############
# No BE between EU and US, test for BE of proposed biosimilar
exp = expand.grid(n=c(1,2,3,4,5,10), 
                  mu.T=c(0.395, 0.40, 0.45, 0.5, 0.5625, 0.625, 0.63), 
                  mu.A=0.5,
                  mu.B=c(0.395, 0.630), 
                  CV=0.2)
sim.results = data.frame(pwr.tost=double(), pwr.ci=double(), pwr.pairwise=double())
start.time = Sys.time()
for (i in 1:nrow(exp)) {
  n = exp$n[i]
  mu.T = exp$mu.T[i]
  mu.A = exp$mu.A[i]
  mu.B = exp$mu.B[i]
  CV = exp$CV[i]
  print(paste('progress: ', i, round(100*i/nrow(exp), 2), '% time lapsed: ', round(Sys.time()-start.time, 2), 'sec'))
  pwr = run.sim(n=n, mu.T=mu.T, mu.A=mu.A, mu.B=mu.B, CV=CV, n.sim=1000, alpha=0.05, seed=100)
  sim.results = rbind(sim.results, data.frame(pwr.tost=pwr[[1]], pwr.ci=pwr[[2]], pwr.pairwise=pwr[[3]]))
}

sim.results2 = cbind(exp, sim.results)
write.csv(sim.results2, paste('output/simulation 2 ', Sys.time(), '.csv', sep=''))
write.csv(sim.results2, 'output/simulation 2.csv')


############
# Plotting #
############

pdf("output/plots.pdf",width = 12, height = 2.5)

# scenario 1
# within BE range, compare power across sample size across method
get.plot1 = function(result, T.val, A.val) {
  dat = result %>% filter(mu.T==T.val) %>% 
    pivot_longer(names_to='power.type', cols = c(pwr.tost, pwr.ci, pwr.pairwise), values_to = 'power') %>%
    filter(power.type!='pwr.ci')
  p = ggplot(data=dat, aes(x=n*6,y=power*100, color=power.type)) +
    geom_line() +
    geom_hline(yintercept=80, linetype='dashed') +
    labs(title=paste('Simulation Result Ratio =', T.val/A.val)) + theme_bw() + 
    scale_color_viridis(discrete=T, begin=0.3, end=0.65, direction=-1, labels=c('Pairwise', 'TOST'))+
    facet_wrap(~CV, ncol=4) +
    xlab('Total Sample Size') + ylab('Power (%)') + ylim(c(0,100)) +
    scale_x_continuous(breaks=seq(10, 60, 5)) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(size=12), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.title=element_blank())
  plot(p)
}

get.plot1(sim.results1, 0.5, 0.5)
get.plot1(sim.results1, 0.4, 0.5)
get.plot1(sim.results1, 0.625, 0.5)
get.plot1(sim.results1, 0.45, 0.5)
get.plot1(sim.results1, 0.5625, 0.5)

# scenario 2
get.plot2 = function(result, T.val, A.val) {
  dat = result %>% filter(mu.T==T.val) %>% 
    pivot_longer(names_to='power.type', cols = c(pwr.tost, pwr.ci, pwr.pairwise), values_to = 'power') %>%
    filter(power.type!='pwr.ci')
  ggplot(data=dat, aes(x=n*6,y=power*100, color=power.type)) +
    geom_line() +
    geom_hline(yintercept=80, linetype='dashed') +
    labs(title=paste('Simulation Result Ratio =', T.val/A.val)) + theme_bw() + 
    scale_color_viridis(discrete=T, begin=0.3, end=0.65, direction=-1, labels=c('Pairwise', 'TOST'))+
    facet_wrap(~mu.B, ncol=2) +
    xlab('Total Sample Size') + ylab('Power (%)') + ylim(c(0,100)) +
    scale_x_continuous(breaks=seq(10,60,5)) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(size=12), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.title=element_blank())
}

get.plot2(sim.results2, 0.5, 0.5)
get.plot2(sim.results2, 0.4, 0.5)
get.plot2(sim.results2, 0.625, 0.5)
get.plot2(sim.results2, 0.45, 0.5)
get.plot2(sim.results2, 0.5625, 0.5)


dev.off()

