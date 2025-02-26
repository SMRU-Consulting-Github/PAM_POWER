### SIMPLER DOSE-RESPONSE BASED ON SAANA CODE

logistic_fun <- function(x, x0 = 0, k = k, L = 1) {
  p <- L / (1 + exp(-k * (x - x0)))
  return(p)
}

dose_resp_fun <- function(ed) {
  RL <- SL - logfac * log10(ed)
  p <- logistic_fun(x = RL, x0 = TH, k = k, L = 1)
  return(p)
}

### DOSE-RESPONSE BASED ON EXPERT ELICITATION

dose_resp_fun_EE <- function(ed,lb=85, ub=200, mean, sd) {
  loss <- (logfac * log10(ed*1000) + (alphaCoeff * ed))
  loss[loss < 0] <- 0
  RL <- SL - loss
  if(is.na(mean)|is.na(sd)) 
    return(rep(NA, length(RL)))
  resp <- (pnorm(RL, mean, sd) - pnorm(lb, mean, sd)) / 
    (pnorm(ub, mean, sd) - pnorm(lb, mean, sd))
  resp[resp < 0] <- 0; resp[resp > 1] <- 1
  return(resp)
}


# for splitting windfarm footprints into subpolygons
split_poly <- function(sf_poly, n_areas) {
  # create random points
  points_rnd <- st_sample(sf_poly, size = 10000)
  # k-means clustering
  points <- do.call(rbind, st_geometry(points_rnd)) %>%
    as_tibble() %>%
    setNames(c("lon", "lat"))
  k_means <- kmeans(points, centers = n_areas)
  # create voronoi polygons
  voronoi_polys <- dismo::voronoi(k_means$centers, ext = sf_poly)
  # clip to sf_poly
  crs(voronoi_polys) <- crs(sf_poly)
  voronoi_sf <- st_as_sf(voronoi_polys)
  equal_areas <- st_intersection(voronoi_sf, sf_poly)
  equal_areas$area <- st_area(equal_areas)
  return(equal_areas)
}

# to redistribute animals based on either baseline or whatever map
shift_density <- function(mov_d, habitat_d, PResp = NULL, A) {
  # mov_d = density to be redistributed/shifted
  # habitat_d = habitat preference
  # PResp = discounted/added habitat preference
  # A is the accessibility (home range, or maximum distance they can be moved)

  if (is.null(PResp)) {
    PResp <- rep(0, length(habitat_d))
  }

  # Target density: habitat preference discounted by P(response)
  targ_d <- habitat_d * (1 - PResp)
  targ_d <- targ_d / sum(targ_d, na.rm = T) # *sum(mov_d)
  n <- length(habitat_d)
  # Destination density: target density, given accessibility
  dest_d <- A * matrix(targ_d, nrow = n, ncol = n, byrow = T)

  # Scale each row (source locations with different destination options) up by source density
  dest_d <- dest_d / matrix(rowSums(dest_d, na.rm=T), nrow = n, ncol = n)
  #dest_d[is.na(dest_d)] <- 0
  #dest_d[is.nan(dest_d)] <- 0
  dest_d <- dest_d * matrix(mov_d, nrow = n, ncol = n)

  # Sum probabilities across destination locations (columns)
  shifted_d <- apply(dest_d, 2, sum, na.rm=T)
  # shifted_d[is.na(shifted_d)] <- 0

  return(shifted_d)
}

# accessibility function (prob space use as a function of distance)
access_fun <- function(x, hr) {
  # x = distance matrix
  # hr = home range
  p <- x <= hr
  p <- p / matrix(rowSums(p), nrow = dim(p)[1], ncol = dim(p)[2])
  return(p)
}

# this and next function was written by Len and I copied it from mcmc_fixer.R
## Simple post fit mcmc for mgcv.
## (c) Simon N. Wood (2020)

## some useful densities (require mgcv::rmvn)...

rmvt <- function(n,mu,V,df) {
  ## simulate multivariate t variates  
  y <- rmvn(n,mu*0,V)
  v <- rchisq(n,df=df)
  t(mu + t(sqrt(df/v)*y))
}

r.mvt <- function(n,mu,V,df) rmvt(n,mu,V,df)

dmvt <- function(x,mu,V,df,R=NULL) {
  ## multivariate t log density...
  p <- length(mu);
  if (is.null(R)) R <- chol(V)
  z <- forwardsolve(t(R),x-mu)
  k <- - sum(log(diag(R))) - p*log(df*pi)/2 + lgamma((df+p)/2) - lgamma(df/2)
  k - if (is.matrix(z)) (df+p)*log1p(colSums(z^2)/df)/2 else (df+p)*log1p(sum(z^2)/df)/2
}

d.mvt <- function(x,mu,V,df,R=NULL) dmvt(x,mu,V,df,R)

dmvn <- function(x,mu,V,R=NULL) {
  ## multivariate normal density mgcv:::rmvn can be used for generation 
  if (is.null(R)) R <- chol(V)
  z <- forwardsolve(t(R),x-mu)
  -colSums(z^2)/2-sum(log(diag(R))) - log(2*pi)*length(mu)/2
}

## some functions to extract important components of joint density from
## fitted gam...

bSb <- function(b,beta=coef(b)) {
  ## evaluate penalty for fitted gam, possibly with new beta
  bSb <- k <-  0
  sp <- if (is.null(b$full.sp)) b$sp else b$full.sp ## handling linked sp's
  for (i in 1:length(b$smooth)) {
    m <- length(b$smooth[[i]]$S)
    if (m) {
      ii <- b$smooth[[i]]$first.para:b$smooth[[i]]$last.para
      for (j in 1:m) {
        k <- k + 1
        bSb <- bSb + sp[k]*(t(beta[ii])%*%b$smooth[[i]]$S[[j]]%*%beta[ii])
      }
    }  
  }
  bSb
} ## bSb

devg <- function(b,beta=coef(b),X=model.matrix(b)) {
  ## evaluate the deviance of a fitted gam given possibly new coefs, beta
  ## for general families this is simply -2*log.lik
  if (inherits(b$family,"general.family")) {
    -2*b$family$ll(b$y,X,beta,b$prior.weights,b$family,offset=b$offset)$l
  } else { ## exp or extended family
    sum(b$family$dev.resids(b$y,b$family$linkinv(X%*%beta+b$offset),b$prior.weights))
  }
} ## devg

lpl <- function(b,beta=coef(b),X=model.matrix(b)) {
  ## log joint density for beta, to within uninteresting constants
  -(devg(b,beta,X)/b$sig2+bSb(b,beta)/b$sig2)/2
}

gam.mh <- function(b,ns=10000,burn=1000,t.df=40,rw.scale=.25,thin=1) {
  ## generate posterior samples for fitted gam using Metroplois Hastings sampler
  ## alternating fixed proposal and random walk proposal, both based on Gaussian
  ## approximation to posterior...
  if (inherits(b,"bam")) stop("not usable with bam fits")
  beta <- coef(b);Vb <- vcov(b)
  X <- model.matrix(b); burn <- max(0,burn)
  prog <- interactive();iprog <- 0
  di <- floor((ns+burn)/100)
  if (prog) prg <- txtProgressBar(min = 0, max = ns+burn, initial = 0,
                                  char = "=",width = NA, title="Progress", style = 3)
  bp <- rmvt(ns+burn,beta,Vb,df=t.df) ## beta proposals
  bp[1,] <- beta ## Don't change this after density step!!
  lfp <- dmvt(t(bp),beta,Vb,df=t.df) ## log proposal density
  
  rw <- is.finite(rw.scale)&&rw.scale>0
  if (rw) {
    R <- chol(Vb) 
    step <- rmvn(ns+burn,beta*0,Vb*rw.scale) ## random walk steps (mgcv::rmvn)
  }
  u <- runif(ns+burn);us <- runif(ns+burn) ## for acceptance check
  bs <- bp;j <- 1;accept <- rw.accept <- 0
  lpl0 <- lpl(b,bs[1,],X)
  for (i in 2:(ns+burn)) { ## MH loop
    ## first a static proposal...
    lpl1 <- lpl(b,bs[i,],X)
    if (u[i] < exp(lfp[j]-lfp[i]+lpl1-lpl0)) {
      lpl0 <- lpl1;accept <- accept + 1
      j <- i ## row of bs containing last accepted beta
    } else bs[i,] <- bs[i-1,]
    ## now a random walk proposal...
    if (rw) {
      lpl1 <- lpl(b,bs[i,]+step[i,],X)
      if (us[i] < exp(lpl1-lpl0)) { ## accept random walk step
        lpl0 <- lpl1;j <- i
        bs[i,] <- bs[i,] + step[i,]
        rw.accept <- rw.accept+1 
        lfp[i] <- dmvt(bs[i,],beta,Vb,df=t.df,R=R) ## have to update static proposal density
      }
    }  
    if (i==burn) accept <- rw.accept <- 0
    if (prog&&i%%di==0) setTxtProgressBar(prg, i)
  } ## MH loop
  if (burn>0) bs <- bs[-(1:burn),]
  if (thin>1) bs <- bs[seq(1,ns,by=thin),]
  if (prog) {
    setTxtProgressBar(prg, i);cat("\n")
    cat("fixed acceptance = ",accept/ns,"  RW acceptance = ",rw.accept/ns,"\n")
  }  
  list(bs=bs,rw.accept = rw.accept/ns,accept=accept/ns)
} ## gam.mh

# some functions to organise/extract data
lalaFunction <- function(data){
  
  tee <- list()
  for (i in 1:length(data)){
    tee[[i]] <- data[[i]]
    tee[[i]] <- as.data.frame(tee[[i]])
    tee[[i]]$pam <- c(1:nrow(tee[[i]]))
    tee[[i]]<- do.call("rbind",tee[[i]])
    #tee[[i]]$month <- monthsPiling[i]
    
  }
  
  tee <- do.call("rbind",tee)
  
  newList <- list() # here each element list is one realization map
  for (i in 1:(ncol(tee)-2)){
    newList[[i]] <- tapply(tee[,i],tee$pam, mean, na.rm=T)
  }
  newList
}

lalaFunction2 <- function(data){
  
  tee <- list()
  for (i in 1:length(data)){
    tee[[i]] <- data[[i]]
    tee[[i]]<- do.call("cbind",tee[[i]])
    tee[[i]] <- as.data.frame(tee[[i]])
    tee[[i]]$month <- monthsPiling[i]
    tee[[i]]$pam <- c(1:nrow(tee[[i]]))
  }
  
  tee <- do.call("rbind",tee)
  tee
}

