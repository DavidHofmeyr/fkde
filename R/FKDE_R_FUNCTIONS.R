makeBetas = function(ord) 1/factorial(0:ord)

n_const = function(betas){
  2*sum(sapply(1:length(betas), function(k) betas[k]*factorial(k-1)))
}

sig2poly = function(betas){
  beta_use = betas/n_const(betas)
  2*sum(sapply(1:length(betas), function(k) beta_use[k]*factorial(k+1)))
}

RKpoly = function(betas){
  beta_use = betas/n_const(betas)
  betakj = beta_use%*%t(beta_use)
  kpj = log((2^(0:(length(betas)-1)))%*%t(2^(0:(length(betas)-1))), base = 2)
  sum(betakj/2^kpj*factorial(kpj))
}

hbs = function(x, ord, deriv = 0){
  coefs = makeBetas(ord)
  coefs = coefs/n_const(coefs)
  if(deriv==0) sd(x)*(RKpoly(coefs)/sig2poly(coefs)^2*8*sqrt(pi)/3/length(x))^(1/5)
  else if(deriv==1) sd(x)*(3*factorial(2*ord)/2^(2*ord+2)/factorial(ord+1)^2/sig2poly(coefs)^2*16/15*sqrt(pi)/length(x))^(1/7)
  else stop('only the density and its first derivative are currently implemented')
}


normGen = function(mus, sigs, ps, n){
  if(length(mus)==1) return(rnorm(n)*sigs[1]+mus[1])
  ns = c(rmultinom(1, n, ps))
  x = numeric(0)
  for(i in 1:length(mus)){
    x = c(x, rnorm(ns[i])*sigs[i]+mus[i])
  }
  sm = sample(1:n, n)
  x[sm]
}

unifGen = function(mus, sigs, ps, n){
  if(length(mus)==1) return(rnorm(n)*sigs[1]+mus[1])
  ns = c(rmultinom(1, n, ps))
  x = numeric(0)
  for(i in 1:length(mus)){
    x = c(x, (runif(ns[i])-.5)*sigs[i]+mus[i])
  }
  sm = sample(1:n, n)
  x[sm]
}


dataGen = function(n, den_number){
  if(den_number==1){
    mus = c(0)
    sigs = c(1)
    props = c(1)
    return(normGen(mus, sigs, props, n))
  }
  if(den_number==2){
    return(runif(n))
  }
  if(den_number==3){
    mus = c(0, 0)
    sigs = c(1, .1)
    props = c(2/3, 1/3)
    return(normGen(mus, sigs, props, n))
  }
  if(den_number==4){
    mus = c(-1, 1)
    sigs = c(2/3, 2/3)
    props = c(.5, .5)
    return(normGen(mus, sigs, props, n))
  }
  if(den_number==5){
    mus = 3*((2/3)^(0:7)-1)
    sigs = (2/3)^(0:7)
    props = numeric(8) + 1/8
    return(normGen(mus, sigs, props, n))
  }
  if(den_number==6){
    mus = c(-1, 1, ((0:6)-3)/2)
    sigs = c(2/3, 2/3, numeric(7)+1/100)
    props = c(49/100, 49/100, numeric(7)+1/350)
    return(normGen(mus, sigs, props, n))
  }
  if(den_number==7){
    mus = c(0, (0:4)/2-1)
    sigs = c(1, numeric(5)+.1)
    props = c(.5, numeric(5)+.1)
    return(normGen(mus, sigs, props, n))
  }
  if(den_number==8){
    mus = c(3*((2/3)^(0:7)-1), -5)
    sigs = c((2/3)^(0:7), 1.5)
    props = c(numeric(8) + 1/8*4/5, .2)
    return(normGen(mus, sigs, props, n))
  }
  if(den_number==9){
    mus = 3*((2/3)^(0:7)-1)
    sigs = (2/3)^(0:7)
    props = numeric(8) + 1/8
    nnorm = sum(runif(n)<.8)
    return(c(normGen(mus, sigs, props, nnorm), runif(n-nnorm)-1))
  }
  if(den_number==10){
    mus = c(1.6,  0.2,  2.0,  0.7, -1.6, -0.2,  0.1,  0.3,  3.4, -0.2)
    sigs = c(0.5, 0.9, 2, 0.2, 0.4, 0.5, 0.8, 1.7, 1.7, 1.9)
    props = c(3, 1, 1, .2, .2, .6, 1, 1, .3, 2)/10.3
    return(unifGen(mus, sigs, props, n))
  }
}

MISE = function(den_number, xs, ys){
  if(den_number==1){
    mus = c(0)
    sigs = c(1)
    props = c(1)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret
    }
  }
  if(den_number==2){
    f = function(x) dunif(x)
  }
  if(den_number==3){
    mus = c(0, 0)
    sigs = c(1, .1)
    props = c(2/3, 1/3)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret
    }
  }
  if(den_number==4){
    mus = c(-1, 1)
    sigs = c(2/3, 2/3)
    props = c(.5, .5)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret
    }
  }
  if(den_number==5){
    mus = 3*((2/3)^(0:7)-1)
    sigs = (2/3)^(0:7)
    props = numeric(8) + 1/8
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret
    }
  }
  if(den_number==6){
    mus = c(-1, 1, ((0:6)-3)/2)
    sigs = c(2/3, 2/3, numeric(7)+1/100)
    props = c(49/100, 49/100, numeric(7)+1/350)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret
    }
  }
  if(den_number==7){
    mus = c(0, (0:4)/2-1)
    sigs = c(1, numeric(5)+.1)
    props = c(.5, numeric(5)+.1)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret
    }
  }
  if(den_number==8){
    mus = c(3*((2/3)^(0:7)-1), -5)
    sigs = c((2/3)^(0:7), 1.5)
    props = c(numeric(8) + 1/8*4/5, .2)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret
    }
  }
  if(den_number==9){
    mus = 3*((2/3)^(0:7)-1)
    sigs = (2/3)^(0:7)
    props = numeric(8) + 1/8
    nnorm = sum(runif(n)<.8)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret = ret + .2*dunif(x, min = -1, max = 0)
      ret
    }
  }
  if(den_number==10){
    mus = c(1.6,  0.2,  2.0,  0.7, -1.6, -0.2,  0.1,  0.3,  3.4, -0.2)
    sigs = c(0.5, 0.9, 2, 0.2, 0.4, 0.5, 0.8, 1.7, 1.7, 1.9)
    props = c(3, 1, 1, .2, .2, .6, 1, 1, .3, 2)/10.3
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]*dunif(x, min = mus[i]-sigs[i]/2, mus[i]+sigs[i]/2)
      ret
    }
  }
  fvals = sapply(xs, f)
  sum((fvals-ys)^2*c(0,diff(xs)))
}

plotDen = function(den_number){
  if(den_number==1){
    mus = c(0)
    sigs = c(1)
    props = c(1)
    xs = seq(min(mus-4*sigs), max(mus+4*sigs), length = 1000)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret
    }
  }
  if(den_number==2){
    xs = seq(-.2, 1.2, length = 1000)
    f = function(x) dunif(x)
  }
  if(den_number==3){
    mus = c(0, 0)
    sigs = c(1, .1)
    props = c(2/3, 1/3)
    xs = seq(min(mus-4*sigs), max(mus+4*sigs), length = 1000)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret
    }
  }
  if(den_number==4){
    mus = c(-1, 1)
    sigs = c(2/3, 2/3)
    props = c(.5, .5)
    xs = seq(min(mus-4*sigs), max(mus+4*sigs), length = 1000)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret
    }
  }
  if(den_number==5){
    mus = 3*((2/3)^(0:7)-1)
    sigs = (2/3)^(0:7)
    props = numeric(8) + 1/8
    xs = seq(min(mus-4*sigs), max(mus+4*sigs), length = 1000)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret
    }
  }
  if(den_number==6){
    mus = c(-1, 1, ((0:6)-3)/2)
    sigs = c(2/3, 2/3, numeric(7)+1/100)
    props = c(49/100, 49/100, numeric(7)+1/350)
    xs = seq(min(mus-4*sigs), max(mus+4*sigs), length = 1000)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret
    }
  }
  if(den_number==7){
    mus = c(0, (0:4)/2-1)
    sigs = c(1, numeric(5)+.1)
    props = c(.5, numeric(5)+.1)
    xs = seq(min(mus-4*sigs), max(mus+4*sigs), length = 1000)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret
    }
  }
  if(den_number==8){
    mus = c(3*((2/3)^(0:7)-1), -5)
    sigs = c((2/3)^(0:7), 1.5)
    props = c(numeric(8) + 1/8*4/5, .2)
    xs = seq(min(mus-4*sigs), max(mus+4*sigs), length = 1000)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret
    }
  }
  if(den_number==9){
    mus = 3*((2/3)^(0:7)-1)
    sigs = (2/3)^(0:7)
    props = numeric(8) + 1/8
    nnorm = sum(runif(n)<.8)
    xs = seq(min(mus-4*sigs), max(mus+4*sigs), length = 1000)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]/sqrt(2*pi)/sigs[i]*exp(-(x-mus[i])^2/2/sigs[i]^2)
      ret = ret + .2*dunif(x, min = -1, max = 0)
      ret
    }
  }
  if(den_number==10){
    mus = c(1.6,  0.2,  2.0,  0.7, -1.6, -0.2,  0.1,  0.3,  3.4, -0.2)
    sigs = c(0.5, 0.9, 2, 0.2, 0.4, 0.5, 0.8, 1.7, 1.7, 1.9)
    props = c(3, 1, 1, .2, .2, .6, 1, 1, .3, 2)/10.3
    xs = seq(min(mus-1.5*sigs), max(mus+1.5*sigs), length = 1000)
    f = function(x){
      ret = 0
      for(i in 1:length(mus)) ret = ret + props[i]*dunif(x, min = mus[i]-sigs[i]/2, mus[i]+sigs[i]/2)
      ret
    }
  }
  fvals = sapply(xs, f)
  plot(xs, fvals, type = 'l', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
}



fkde = function(x, h = NULL, ord = 1, ngrid = 0, nbin = 0, m = NULL, M = NULL, hnorm = NULL){
  if(!is.null(hnorm)) h = hnorm*hbs(x,ord,0)/(sd(x)*(4/3/length(x))^.2)
  else if(is.null(h)) h = hbs(x, ord, 0)
  coefs = makeBetas(ord)
  coefs = coefs/n_const(coefs)
  if(is.null(m)) m = min(x)-6*h
  if(is.null(M)) M = max(x)+6*h
  if(ngrid==0 && nbin==0){
    list(x = sort(x), y = f_kde_mix(x, h, length(x), ord, coefs))
  }
  else if(nbin==0){
    list(x = seq(m, M, length = ngrid), y = f_kde_mix_grid(x, h, length(x), ord, coefs, ngrid, m, M))
  }
  else{
    list(x = seq(m, M, length = nbin), y = f_kde_mix_binned(x, h, length(x), ord, coefs, nbin, m, M))
  }
}

dfkde = function(x, h = NULL, ord = 1, hnorm = NULL){
  if(!is.null(hnorm)) h = hnorm*hbs(x, ord, 1)/(3/4/sqrt(pi)/(15/16/sqrt(pi)))^(1/7)/sd(x)
  else if(is.null(h)) h = hbs(x, ord, 1)
  list(x = sort(x), y = df_kde(x, h, length(x), ord))
}
