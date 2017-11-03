
# numerical simulations of the 'score' statistics
# for 2017 reading group on complexity


lBinom <- function(x,p,n){
  L = (n-x)*log(1-p)+x*log(p)+log(choose(n,x));
}

# score = dtheta (L(theta,x))
scoreBinom <- function(x,p,n){
  (n*p-x)/(p*(p-1))
}


n=25;
Nsims= 1000;
ptrue=.5;
simdata_binom <- rbinom(Nsims, n, ptrue);
pvals = seq(0,1,.01);

scores_binom = c()

for (j in 1:Nsims){
  obs <- simdata_binom[j]
  s = scoreBinom(obs,pvals,n)
  scores_binom = cbind2(scores_binom, s)
}

meanscores_binom = rowMeans(scores_binom)
#plot(pvals, meanscores_binom)
  


#################

scoreGamma <- function(x,alpha){
  log(x)-digamma(alpha)
}

alphavals=seq(0,1,.1);
alphatrue=1.5;
simdata_gamma <- rgamma(Nsims, alphatrue);

scores_gamma = c()

for (j in 1:Nsims){
  obs <- simdata_gamma[j]
  s = scoreGamma(obs,alphavals)
  scores_gamma = cbind2(scores_gamma, s)
}


meanscores_gamma = rowMeans(scores_gamma)
#plot(alphavals, meanscores_gamma)



## ####################

Nsims <- 25;
mutrue = -0.5;
sigmatrue = 1;
muvals = seq(-2,2,.01)

scoreNormal <- function(x,mu){
    x-mu
}

simdata_n  <- rnorm(Nsims,mutrue,sigmatrue)
scores_n = c()

for (j in 1:Nsims){
    obs <- simdata_n[j]
    s = scoreNormal(obs,muvals)
    scores_n = cbind2(scores_n, s)
}


meanscores_n = rowMeans(scores_n)
plot(muvals, meanscores_n)


