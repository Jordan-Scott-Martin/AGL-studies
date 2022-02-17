#Bayesian power analysis for univariate comparison

#set working directory
setwd("")

#load data
copy<-read.csv("copy data.csv")

#sample
  N_d = 15
  N_h = 15
  trials = nrow(copy[copy$Subject==1,])
  
#difference in means (effect size)
  es_d = 0.3 #small
  logit_sd = (pi/sqrt(3))
  beta_h = es_d*logit_sd
  
#intercept
  mu = 1.5 #p ~ 0.80
  
#subject intercepts
  id_sd = 0.5
  re_d = rnorm(N_d, 0, id_sd)
  re_h = rnorm(N_h, 0, id_sd)
  
#latent values
  eta_d = mu + re_d
  eta_h = mu + beta_h + re_h
  
#simulate bernoulli trials
  logistic = function (x) {p <- 1/(1 + exp(-x))
                           p <- ifelse(x == Inf, 1, p)
                           p}
  y_d = rbinom(n = length(eta_d)*trials, size = 1, prob = rep(logistic(eta_d), each = trials)) 
  y_h = rbinom(n = length(eta_h)*trials, size = 1, prob = rep(logistic(eta_h), each = trials)) 
  df_d = data.frame(group = "d", y = y_d, id = rep(seq(1:N_d),each = trials))
  df_h = data.frame(group = "h", y = y_h, id = rep(seq(from = c(N_h+1), to = c(N_h+N_d)),each = trials))
  df = rbind(df_d,df_h)
  
#run model
  library(brms)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  mf = bf(y ~ group + (1|id)) + bernoulli()
  
  m_1 <- brm(formula = mf, data = df,
               prior=c(prior("normal(0,2)",class="Intercept"),
                       prior("normal(0,2)",class="b"),
                       prior("cauchy(0,2)",class="sd")),
               warmup=1000,iter=2000, chains=4, seed=9,
               control=list(adapt_delta=0.80))
  

#list for posterior probabilities (200 datasets)
result = list()

for(i in 1:200){
  #sample
  N_d = 15
  N_h = 15
  trials = nrow(copy[copy$Subject==1,])
  
  #difference in means (effect size)
  es_d = 0.3 #small
  logit_sd = (pi/sqrt(3))
  beta_h = es_d*logit_sd
  
  #intercept
  mu = 1.5 #p ~ 0.80
  
  #subject intercepts
  id_sd = 0.5
  re_d = rnorm(N_d, 0, id_sd)
  re_h = rnorm(N_h, 0, id_sd)
  
  #latent values
  eta_d = mu + re_d
  eta_h = mu + beta_h + re_h
  
  #simulate bernoulli trials
  logistic = function (x) {p <- 1/(1 + exp(-x))
                           p <- ifelse(x == Inf, 1, p)
                           p}
  y_d = rbinom(n = length(eta_d)*trials, size = 1, prob = rep(logistic(eta_d), each = trials)) 
  y_h = rbinom(n = length(eta_h)*trials, size = 1, prob = rep(logistic(eta_h), each = trials)) 
  df_d = data.frame(group = "d", y = y_d, id = rep(seq(1:N_d),each = trials))
  df_h = data.frame(group = "h", y = y_h, id = rep(seq(from = c(N_h+1), to = c(N_h+N_d)),each = trials))
  df = rbind(df_d,df_h)
  
  #run model
  m_1 <- update(m_1, newdata = df)
  
  result[[i]] = hypothesis(m_1, "grouph > 0")$hypothesis$Post.Prob
}

bin_result = apply(matrix(result), 1, FUN = function(x) ifelse(x > 0.95, 1, 0))

sum(bin_result)/length(bin_result)
saveRDS(result, "result_1.RDS")
result = readRDS("result_1.RDS")

