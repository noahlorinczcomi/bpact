#' Trapezoidal density
#'
#' This function calculates the density of the trapezoidal distribution
#' @param x argument to the density function
#' @param a smallest parameter of the trapezoidal density function
#' @param x1 second-smallest parameter of the trapezoidal density function
#' @param x2 second-largest parameter of the trapezoidal density function
#' @param b largest parameter of the trapezoidal density function
#' @param log TRUE if log density should be returned, FALSE otherwise
#' @keywords trapezoidal
#' @export
#' @examples dtrap(1,0,1,2,3)
#' dtrap()
dtrap=function(x,a,x1,x2,b,log=FALSE) {
  # trapezoidal pdf, vectorized version
  y0=2/(b-a-x1+x2)
  res=rep(0,length(x))
  res[x>=a & x<x1]=y0/(x1-a)*(x[x>=a & x<x1]-a)
  res[x>=x1 & x<x2]=y0
  res[x>=x2 & x<b]=y0/(b-x2)*(b-x[x>=x2 & x<b])
  if(log) return(log(res))
  return(res)
}

#' This function calculates the cumulative density of the trapezoidal distribution
#'
#' @param x argument to the cumulative trapezoidal density function
#' @param a smallest parameter of the trapezoidal density function
#' @param x1 second-smallest parameter of the trapezoidal density function
#' @param x2 second-largest parameter of the trapezoidal density function
#' @param b largest parameter of the trapezoidal density function
#' @param log TRUE if log density should be returned, FALSE otherwise
#' @keywords trapezoidal
#' @export
#' @examples ptrap(2,0,1,2,3)
#' ptrap()
ptrap=function(x,a,x1,x2,b,log=FALSE) {
  x=sort(x)
  y0=2/(b-a-x1+x2)
  # known integrals: int_c1^c2
  int_a_x1=y0*((x1^2-a^2)/(2*x1-2*a)-a)
  int_x1_x2=y0*(x2-x1)
  int_x2_b=1-int_a_x1-int_x1_x2
  # int_{c}^x
  int_a_x=y0/((2*x1-a))*(x^2-a^2)-a*y0/(x1-a)*(x-a)
  int_x1_x=y0*(x-x1)
  int_x2_x=b*y0/(b-x2)*(x-x2)-y0/(2*(b-x2))*(x^2-x2^2)
  res=rep(0,length(x))
  # indices
  ix_a_x1=which(x>=a & x<=x1)
  ix_x1_x2=which(x>=x1 & x<=x2)
  ix_x2_b=which(x>=x2 & x<=b)
  # CDF
  res[ix_a_x1]=int_a_x[ix_a_x1]
  res[ix_x1_x2]=int_a_x1+int_x1_x[ix_x1_x2]
  res[ix_x2_b]=int_a_x1+int_x1_x2+int_x2_x[ix_x2_b]
  # end
  if(log) return(log(res)) else (return(res))
}

#' This function numerically calculates the quantile function of the trapezoidal distribution
#'
#' @param p the desired quantile defines this level of probability
#' @param a smallest parameter of the trapezoidal density function
#' @param x1 second-smallest parameter of the trapezoidal density function
#' @param x2 second-largest parameter of the trapezoidal density function
#' @param b largest parameter of the trapezoidal density function
#' @param eps error in numerical approximation of the true quantile
#' @param max_passes number of intervals, of decreasing width, in which the quantile should be approximated
#' @param n number of grid points at which to evaluate the quantile in each interval
#' @keywords trapezoidal
#' @export
#' @examples qtrap(0.5,0,1,2,3)
#' qtrap()
qtrap=function(p,a,x1,x2,b,eps=1e-5,max_passes=5,n=100) {
  # numerical approach to finding the quantile
  # multiple passes of increasing precision
  # increasing precision as more passes
  iter=0;error=eps+1
  bounds=c(a,b)
  while(error>eps & iter<max_passes) {
    iter=iter+1
    s=seq(bounds[1],bounds[2],length.out=n)
    cdf=ptrap(s,a,x1,x2,b)
    errors=abs(cdf-p)
    ix=which.min(errors)
    error=errors[ix]
    bounds=c(s[ix-1],s[ix+1])
    if(length(bounds)==1) bounds=c(bounds/1e3,bounds)
    quant=s[ix]
  }
  quant
}

#' This function geometrically samples from the trapezoidal distribution
#'
#' @param n the desired number of random draws
#' @param a smallest parameter of the trapezoidal density function
#' @param x1 second-smallest parameter of the trapezoidal density function
#' @param x2 second-largest parameter of the trapezoidal density function
#' @param b largest parameter of the trapezoidal density function
#' @keywords trapezoidal
#' @export
#' @examples rtrap(10,0,1,2,3)
#' rtrap()
rtrap=function(n,a,x1,x2,b) {
  y0=2/(b-a-x1+x2)
  x=c()
  for(i in 1:n) {
    r2=Inf
    r1=runif(1,a,b)
    while(r2>dtrap(r1,a,x1,x2,b)) {
      r1=runif(1,a,b)
      r2=runif(1,0,y0)
    }
    x[i]=r1
  }
  x
}

#' This function numerically estimates the parameters of a trapezoidal distribution using observations
#'
#' @param taus vector of observed random variables following a trapezoidal distribution
#' @param subn number of sub-sampled observations to use if length of taus is very large
#' @param atrim the lower quantile of of observed values to consider, such as 0.05
#' @param btrim the upper quantile of of observed values to consider, such as 0.95
#' @param doplot TRUE if a plot should be generated, FALSE otherwise
#' @keywords trapezoidal
#' @export
#' @examples numerical_trap_distribution(c(1,3,4,5,7))
#' numerical_trap_distribution()
numerical_trap_distribution=function(taus,subn=length(taus),atrim=0,btrim=0.975,doplot=FALSE) {
  ## numerical estimation of parameters of trapezoidal distribution
  # actually compare empirical density to theoretical density
  ecdf=function(x) {
    xs=cbind(1:length(x),x)[order(x),]
    xc=cbind(1:length(x),xs)
    xc=cbind(xc[,1]/nrow(xc),xc)
    xc=xc[order(xc[,3]),]
    xc[,1]
  }
  # sub-sampled taus
  staus=sample(taus,min(c(length(taus),subn)),replace=FALSE)
  staus=sort(staus)
  # trim outliers?
  # where to place b?
  q=quantile(staus,prob=c(atrim,btrim))
  a=unname(q[1])
  b=unname(q[2])
  ed=ecdf(staus)
  # grid of possible (x1,x2) pairs
  eg=expand.grid(seq(a,b,length.out=50),seq(a,b,length.out=50))
  colnames(eg)=c('x1','x2')
  eg=eg[eg[,1]<eg[,2],]
  eg=eg[eg[,1]>a,]
  eg=eg[eg[,2]<b,]
  pen=c()
  for(i in 1:nrow(eg)) {
    x1=eg$x1[i]
    x2=eg$x2[i]
    cdf=ptrap(staus,a,x1,x2,b)
    # 2-norm between eCDF and proposal CDF as a penalty
    pen[i]=norm(ed-cdf,'2')^2
  }
  ix=which.min(pen)
  # plot?
  if(doplot) {
    par(mfrow=c(1,3),mar=c(5.1,6,4.1,2.1),mgp=c(1,1/2,0))
    # optimization (2-norm ^2 of CDF-eCDF)
    plot(log(pen),col='#0000004B',xaxt='n',yaxt='n',
         xlab='trapezoidal distribution parameters',
         ylab=expression('||'*italic(F)*'(*'*italic(x)*')-'*italic(F)[e]*'('*italic(x)*')||'[2]^2))
    points(ix,log(min(pen)),bg='red',pch=21)
    # histogram of tau vs inferred density
    par(mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0))
    hist(staus,pr=T,col='gray90',border='gray90',xlab=expression(tau),main='')
    abline(v=c(a,b),lty=2)
    curve(dtrap(x,a,eg[ix,1],eg[ix,2],b),min(staus),max(staus),add=T)
    # CDF vs eCDF plot
    plot(staus,ed,type='l',xlim=c(a,b),
         xlab=expression(tau),
         ylab=expression(italic(P)*'('*italic(X)*'<'*italic(x)*')'))
    curve(ptrap(x,a,eg[ix,1],eg[ix,2],b),a,b,add=T,col='blue')
    legend('topleft',
           legend=c(expression(italic(hat(F))*'('*italic(x)*')'),
                    expression(italic(F)[E]*'('*italic(x)*')')),
           lty=c(1,1),col=c('black','blue'),cex=3/4)
    par(mfrow=c(1,1))
  }
  c(hard_a=min(taus),soft_a=a,x1=eg$x1[ix],x2=eg$x2[ix],soft_b=b,hard_b=max(taus))
}

#' This function uses the method of moments approximation to estimate the parameters of a trapezoidal distribution from observed values when the boundaries of the observable space are known
#'
#' @param taus vector of observed random variables following a trapezoidal distribution
#' @param a smallest parameter of the trapezoidal density function
#' @param b largest parameter of the trapezoidal density function
#' @keywords trapezoidal MoM
#' @export
#' @examples mom_trap_distribution(c(1,3,4,5,7),0,8)
#' numerical_trap_distribution()
mom_trap_distribution=function(taus,a,b) {
  ## approx. method of moments for trapezoidal distribution
  # a: minimum allowable tau value
  # b: maximum allowable tau value
  # not actually a method of moments - searches a grid to find most likely (x1,x2) values
  eg=expand.grid(seq(a,b,length.out=100),seq(a,b,length.out=100))
  colnames(eg)=c('x1','x2')
  eg=eg[eg[,1]<=eg[,2],]
  eg$mean=eg$variance=NA
  for(i in 1:nrow(eg)) {
    x1=eg$x1[i]
    x2=eg$x2[i]
    ex=1/(3*(b+x2-a-x1))*((b^3-x2^3)/(b-x2)-(x1^3-a^3)/(x1-a))
    vx=1/(6*(b+x2-a-x1))*((b^4-x2^4)/(b-x2)-(x1^4-a^4)/(x1-a))-ex^2
    eg$mean[i]=ex
    eg$variance[i]=vx
  }
  EV=matrix(c(mean(logtaus),var(logtaus)),nr=nrow(eg),nc=2,byrow=TRUE)
  ix=as.matrix(eg[,c('mean','variance')])
  ix=(EV[,1]-ix[,1]+EV[,2]-ix[,2])^2
  egdf=eg[which.min(ix),c('x1','x2')]
  # hist(taus,pr=T,col='gray90',border='gray90')
  # curve(dtrap(x,a,egdf$x1,egdf$x2,b),min(taus),max(taus),add=T,lty=2)
  c('x1'=egdf$x1,'x2'=egdf$x2)
}

#' This function uses the method of moments to estimate the parameters of a Beta distribution from observed values
#'
#' @param deltas vector of observed random variables following a Beta distribution
#' @param deltas estimate of the variance of the observed values
#' @keywords Beta MoM
#' @export
#' @examples mom_beta_distribution(c(1,3,4,5,7)/10,1/20)
#' mom_beta_distribution()
mom_beta_distribution=function(deltas,var_deltas=NULL) {
  ## method of moments for Beta distribution
  if(is.null(var_deltas)) var_deltas=var(deltas)
  mom_alpha=mean(deltas)*(mean(deltas)*(1-mean(deltas))/var_deltas-1)
  mom_beta=mom_alpha/mean(deltas)*(1-mean(deltas))
  c(alpha=mom_alpha,beta=mom_beta)
}

#' This function uses the Metropolis-Hastings algorithm to sample from the posterior distributions of scaled heritability and proportion of causal genes
#'
#' @param gentres_chr_ldblock gene-based association test results for a single chromosome for a random set of independent genes
#' @param lddf_chr_ldblock weighted LD scores for a single chromosome for a random set of independent genes
#' @param nk GWAS sample size
#' @param chain_length length of the Markov chain to use
#' @param burnin number of originally drawn samples to exclude from the observed posterior distribution
#' @param madj TRUE if you want to divide likelihoods by the number of genes used, FALSE otherwise
#' @param tau_a smallest parameter of the trapezoidal density function that is the prior for scaled heritability
#' @param tau_x1 second-smallest parameter of the trapezoidal density function that is the prior for scaled heritability
#' @param tau_x2 second-largest parameter of the trapezoidal density function that is the prior for scaled heritability
#' @param tau_b largest parameter of the trapezoidal density function that is the prior for scaled heritability
#' @param beta_shape2 second shape parameter of the prior distribution of the proportion of non-causal genes
#' @param verbose. TRUE if results should be printed as links in the Markov chain are added, FALSE otherwise
#' @keywords prior
#' @export
#' @examples
#' mh()
mh=function(gentres_chr_ldblock,lddf_chr_ldblock,nk,
            chain_length=100000,burnin=100,madj_likelihood=FALSE,
            tau_a=2e-10,tau_x1=2e-9,tau_x2=2e-7,tau_b=2e-4,beta_shape2=1.0001,verbose.=TRUE) {
  # Metropolis-Hastings function with built-in priors for delta and tau
  ## chromosome-specific GenT data set subsetted to one randomly selected gene per LD block
  ## chromosome-specific LD score data set subsetted to one randomly selected gene per LD block
  # extract data
  mu0=gentres_chr_ldblock$m
  sigma0=gentres_chr_ldblock$gent_sigma2_h0
  Qk=gentres_chr_ldblock$gent_test_statistic
  nk=nk
  wldscores=lddf_chr_ldblock$sumwldscore
  wld2scores=lddf_chr_ldblock$sumwld2score
  # function to find estimate of alpha in Beta(alpha,~1) given single observation
  myf=function(par,x) -dbeta(x,par,beta_shape2,log=TRUE)
  # likelihood function
  lf=function(delta,tau,mu0,sigma0,mu1,nk,wldscores,wld2scores,adjm=TRUE) {
    beta0=mu0/sigma0
    alpha0=mu0*beta0
    e1=mu0+tau*nk*wldscores
    v1=sigma0+2*(nk*tau)^2*wld2scores+4*nk*tau*wldscores
    beta1=e1/v1
    alpha1=e1*beta1
    f0=dgamma(Qk,shape=alpha0,rate=beta0)
    f1=dgamma(Qk,shape=alpha1,rate=beta1)
    f0[f0==0]=min(f0[f0>0])
    f1[f1==0]=min(f1[f1>0])
    # this is the conditional likelihood (conditional on params)
    # should not be P(I_k|T_k)=1. that is for the conditional marginalized distribution
    ll=sum(log(delta*f0+(1-delta)*f1))
    if(adjm) ll=ll/length(mu0)
    ll
  }
  # hyperparameters
  tauchain=deltachain=accepted=c()
  # start chains with chosen values
  tauchain[1]=4e-7
  deltachain[1]=0.99
  # start algorithm
  for(iter in 2:chain_length) {
    if(verbose. & iter%%round(chain_length*0.1)==0) cat('  ',iter,'th link in chain\n',sep='')
    delta0=deltachain[iter-1]
    tau0=tauchain[iter-1]
    # draw new delta, tau
    alphahat=optim(100,myf,x=delta0,lower=3,upper=300,method='Brent')$par
    delta1=rbeta(1,alphahat,beta_shape2)
    tau1=rtrap(1,tau_a,tau_x1,tau_x2,tau_b)
    # densities of delta and tau
    d_delta0=dbeta(delta0,100,beta_shape2,log=TRUE)
    d_delta1=dbeta(delta1,100,beta_shape2,log=TRUE)
    d_tau0=dtrap(tau0,tau_a,tau_x1,tau_x2,tau_b,log=TRUE)
    d_tau1=dtrap(tau1,tau_a,tau_x1,tau_x2,tau_b,log=TRUE)
    # likelihood (assumes independence between genes)
    like0=lf(delta0,tau0,mu0,sigma0,mu1,nk,wldscores,wld2scores,adjm=madj_likelihood)
    like1=lf(delta1,tau1,mu0,sigma0,mu1,nk,wldscores,wld2scores,adjm=madj_likelihood)
    # posterior likelihood ratio
    gam=d_delta1+d_tau1+like1-d_delta0-like0-d_tau0
    # acceptance probability
    aprob=exp(gam)
    if(aprob>runif(1)) {
      deltachain[iter]=delta1
      tauchain[iter]=tau1
      accepted[iter-1]=1
    } else {
      deltachain[iter]=delta0
      tauchain[iter]=tau0
      accepted[iter-1]=0
    }
  }
  if(burnin>0) {
    deltachain=deltachain[-(1:burnin)]
    tauchain=tauchain[-c(1:burnin)]
  }
  data.frame(delta=deltachain,tau=tauchain)
}

#' This function refines the prior distributions of scaled heritability and the proportion of noncausal independent genes
#'
#' @param gent.data data frame of full gene-based association test results
#' @param ld.df data frame of weighted LD scores for all genes in gent.data
#' @param gent.Rho list of matrices. Each matrix in this list corresponds to a chromosome and represents correlations between gene-based association test statistics
#' @param gwasn GWAS sample size
#' @param verbose TRUE if progress should be printed, FALSE otherwise
#' @param niter number of iterations for which random sampling of independent genes should be done
#' @param thr_r2 genes correlated beyond this threshold will be placed in the same block; the thr_r2 argument in bigsnpr::snp_ldsplit
#' @param min_size the minimum size of blocks of genes; the min_size argument in bigsnpr::snp_ldsplit
#' @param max_size the maximum size of blocks of genes; the min_size argument in bigsnpr::snp_ldsplit
#' @param max_K the maximum number of gene blocks; the max_K argument in bigsnpr::snp_ldsplit
#' @param max_r2 the maximum allowable correlation between genes in separate blocks; the max_r2 argument in bigsnpr::snp_ldsplit
#' @param ... additional arguments passed to the mh function for prior estimation
#' @keywords prior
#' @export
#' @examples
#' compositemh()
compositemh=function(gent.data,ld.df,gent.Rho,gwasn,verbose=TRUE,niter=10,
                     thr_r2=0.5,min_size=1,max_size=100,max_K=1e6,max_r2=0.75,...) {
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # gent.data: GenT summary statistics, these columns assumed:
  #   - gene (gene symbol)
  #   - chr (chromsome)
  #   - mid (gene base pair start/end midpoint)
  #   - pval (GenT P-value)
  #   - m (number SNPs used)
  #   - gent_sigma2_h0 (null variance)
  #   - gent_test_statistic (GenT test statistic)
  # ld.df: LD data frame for genes, these columns assumed:
  #   - chr, gene, m, sumwldscore, sumwld2score (MAF-weighted; see paper)
  # gent.Rho: list of GenT statistic correlation matrices for each chromosome, these names assumed:
  #   - chr1, chr2, ..., chr22
  # gwasn: GWAS sample size
  # verbose: should (chr,composite likelihood iteration) progress be printed?
  # niter: number of composite likelihood iterations
  # return.fullres: TRUE/FALSE; should chr- and composite likelihood-specific values be returned?
  # doplot: should empirical densities be plotted for each chromosome?
  # ...: additional arguments passed to mh(), including chain_length, burnin, tau_a, tau_x1, tau_x2, tau_b
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # make sure the data have the necessary columns
  # loop over chromosomes
  chrs=unique(gent.data$chr)
  chrs=sort(as.numeric(chrs))
  fullres=list()
  posterior_parameters=data.frame()
  for(cc in 1:length(chrs)) {
    if(verbose) cat('chr ',chrs[cc],'\n',sep='')
    # filter GenT results to just this chromosomes and arrange by base pair start/end midpoint
    gentres_chr=gent.data %>% filter(chr==chrs[cc]) %>% arrange(mid)
    # LD scores
    lddf_chr=ld.df %>% filter(chr==chrs[cc]) %>% filter(gene %in% gentres_chr$gene)
    ix=sapply(gentres_chr$gene,function(h) which(lddf_chr$gene==h))
    lddf_chr=lddf_chr[ix,]
    # GenT stat correlation matrix
    chrmat=gent.Rho[[which(names(gent.Rho)==paste0('chr',chrs[cc]))]]
    # match order of genes in gentres_chr, lddf_chr, and chrmat to be the same
    usegenes=intersect(gentres_chr$gene,rownames(chrmat))
    gentres_chr=gentres_chr %>% filter(gene %in% usegenes)
    lddf_chr=lddf_chr %>% filter(gene %in% usegenes)
    ix=which(rownames(chrmat) %in% usegenes)
    chrmat=chrmat[ix,ix]
    ix=sapply(gentres_chr$gene,function(h) which(rownames(chrmat)==h))
    chrmat=chrmat[ix,ix]
    # define correlation blocks (may need some adjusting)
    Rho=as(chrmat,'sparseMatrix')
    blockres=bigsnpr::snp_ldsplit(Rho,
                                  thr_r2=thr_r2,
                                  min_size=min_size,
                                  max_size=min(c(nrow(Rho),max_size)),
                                  max_K=max_K,
                                  max_r2=max_r2)
    blocks=unlist(tail(blockres$all_last,1))
    if(is.null(blocks)) blocks=0:nrow(Rho) else blocks=c(0,blocks)
    # multiple iterations of composite likelihood
    mh_deltares=mh_taures=list()
    for(iter in 1:niter) {
      ixs=c()
      # loop over each correlation block
      for(i in 1:(length(blocks)-1)) {
        # sample 1 observation in each correlation block
        blockinds=(blocks[i]+1):blocks[i+1]
        if(length(blockinds)==1) {ixs[i]=blockinds;next}
        ixs[i]=sample(blockinds,1)
      }
      # show progress?
      if(verbose) cat(' iteration',iter,'\n')
      # Metroplis-Hastings algorithm to estimate prior distribution hyperparameters for delta and tau
      mh_fit=mh(
        gentres_chr_ldblock=gentres_chr[ixs,],
        lddf_chr_ldblock=lddf_chr[ixs,],
        nk=gwasn,verbose.=verbose,...)
      mh_deltares[[iter]]=mh_fit$delta
      mh_taures[[iter]]=mh_fit$tau
    }
    DELTA=do.call(cbind,mh_deltares)
    TAU=do.call(cbind,mh_taures)
    # estimate parameters of posterior distributions of delta and tau
    ## NOTE: they are in fact correlated 'posterior-ily'; still safe to assume independent priors, though
    # see: diag(cor(DELTA,TAU))
    toadd=data.frame(
      # overall means, same as mean of iteration-specific means
      mean_delta=mean(DELTA), mean_tau=mean(TAU),
      # overall medians, NOT same as mean/median of iteration-specific medians
      median_delta=median(DELTA),median_tau=median(TAU),
      # variances of delta and tau estimates (within + between iterations)
      sd_delta=mean(apply(DELTA,2,var))+sum((colMeans(DELTA)-mean(DELTA))^2)/(niter-1),
      sd_tau=mean(apply(TAU,2,var))+sum((colMeans(TAU)-mean(TAU))^2)/(niter-1),
      # mean quantiles
      quant5_delta=mean(apply(DELTA,2,function(h) quantile(h,prob=0.05))),
      quant95_delta=mean(apply(DELTA,2,function(h) quantile(h,prob=0.95))),
      quant5_tau=mean(apply(TAU,2,function(h) quantile(h,prob=0.05))),
      quant95_tau=mean(apply(TAU,2,function(h) quantile(h,prob=0.95))),
      chr=chrs[cc]
    )
    # record results
    fullres[[cc]]=list(mh_deltares=DELTA,mh_taures=TAU)
    names(fullres)[cc]=paste0('chr',chrs[cc])
    posterior_parameters=rbind(posterior_parameters,toadd)
  }
  # return results
  out=list(posterior_parameters=posterior_parameters,fullres=fullres)
  return(out)
}


#' This function returns gene-level posterior probabilities of causality for a single trait that do not consider correlation between genes
#'
#' @param gent.data data frame of full gene-based association test results
#' @param ld.df data frame of weighted LD scores for all genes in gent.data
#' @param gent.Rho list of matrices. Each matrix in this list corresponds to a chromosome and represents correlations between gene-based association test statistics
#' @param gwasn GWAS sample size
#' @param prior.estimation.list direct output of the compositemh function containing information regarding the prior distributions to use in posterior estimation
#' @param verbose TRUE if progress should be printed, FALSE otherwise
#' @keywords posterior
#' @export
#' @examples
#' posterior_gene()
posterior_gene=function(gent.data,ld.df,gwasn,prior.estimation.list,verbose=TRUE) {
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # gent.data: GenT summary statistics, these columns assumed:
  #   - gene (gene symbol)
  #   - chr (chromsome)
  #   - mid (gene base pair start/end midpoint)
  #   - pval (GenT P-value)
  #   - m (number SNPs used)
  #   - gent_sigma2_h0 (null variance)
  #   - gent_test_statistic (GenT test statistic)
  # ld.df: LD data frame for genes, these columns assumed:
  #   - chr, gene, m, sumwldscore, sumwld2score (MAF-weighted; see paper)
  # gwasn: GWAS sample size
  # prior.estimation.list: direct output of compositemh() (both the parameters and the chains; need chains for prior distribution estimation)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # posterior P(causal) for each gene (a conditional marginalized distribution)
  # posterior probability that gene k is causal for trait t

  # estimate parameters of delta~Beta(alpha,beta) and tau~Trap(a,x1,x2,b)
  ## deltas (chromosome-specific)
  DELTAS=list()
  chrs=unique(names(prior.estimation.list$fullres))
  for(o in 1:length(chrs)) DELTAS[[o]]=prior.estimation.list$fullres[[o]]$mh_deltares
  names(DELTAS)=chrs
  mom_delta=lapply(DELTAS,function(h) mom_beta_distribution(c(h)))
  mom_delta=matrix(unlist(mom_delta),nc=2,byrow=T)
  rownames(mom_delta)=chrs
  ## taus (chromosome-specific)
  TAUS=list()
  for(o in 1:22) TAUS[[o]]=prior.estimation.list$fullres[[o]]$mh_taures
  names(TAUS)=chrs
  mom_tau=lapply(TAUS,function(h) numerical_trap_distribution(c(h),subn=1e5,atrim=0,btrim=0.975,doplot=FALSE))
  mom_tau=matrix(unlist(mom_tau),nc=6,byrow=T)
  rownames(mom_tau)=chrs
  ## store
  toadd=data.frame(
    chr=chrs,
    delta_param1=mom_delta[,1],
    delta_param2=mom_delta[,2],
    tau_param_hard_a=mom_tau[,1],
    tau_param_soft_a=mom_tau[,2],
    tau_param_x1=mom_tau[,3],
    tau_param_x2=mom_tau[,4],
    tau_param_soft_b=mom_tau[,5],
    tau_param_hard_b=mom_tau[,6])
  # prior density functions needed for posterior calculation (ie the distributions integrated over)
  pz=function(m,gent_sigma2_h0,gent_test_statistic,sumwldscore,sumwld2score,nk,delta,tau) {
    mu0=m
    sigma0=gent_sigma2_h0
    Tk=gent_test_statistic
    beta0=mu0/sigma0
    alpha0=mu0*beta0
    e1=mu0+tau*nk*sumwldscore
    v1=sigma0+2*(nk*tau)^2*sumwld2score+4*nk*tau*sumwldscore
    beta1=e1/v1
    alpha1=e1*beta1
    # f0=dgamma(Tk,shape=alpha0,rate=beta0)
    # f1=dgamma(Tk,shape=alpha1,rate=beta1)
    # prob=f1*(1-delta)/(f1*(1-delta)+f0*delta)
    ## log version to avoid numerical instability
    log_f0=dgamma(Tk,shape=alpha0,rate=beta0,log=TRUE)
    log_f1=dgamma(Tk,shape=alpha1,rate=beta1,log=TRUE)
    log_ratio=log_f0-log_f1  # Compute the logarithm of the ratio f0/f1
    delta_ratio=delta/(1-delta)  # Ratio of delta to 1-delta
    prob=1/(1+exp(log_ratio)*delta_ratio) # compute
    return(prob)
  }
  ptau=function(tau,a=tau_a,x1=tau_x1,x2=tau_x2,b=tau_b) {
    # trapezoidal pdf, vectorized version
    y0=2/(b-a-x1+x2)
    res=rep(0,length(tau))
    res[tau>=a & tau<x1]=y0/(x1-a)*(tau[tau>=a & tau<x1]-a)
    res[tau>=x1 & tau<x2]=y0
    res[tau>=x2 & tau<b]=y0/(b-x2)*(b-tau[tau>=x2 & tau<b])
    return(res)
  }
  pdelta=function(delta,delta_alpha,delta_beta) dbeta(delta,delta_alpha,delta_beta)
  f=function(delta,tau,
             m,gent_sigma2_h0,gent_test_statistic, # gene test statistics
             sumwldscore,sumwld2score, # LD statistics
             nk, # gwas N
             a=tau_a,x1=tau_x1,x2=tau_x2,b=tau_b, # priors of tau
             delta_alpha,delta_beta) {
    pz(m,gent_sigma2_h0,gent_test_statistic,
       sumwldscore,sumwld2score,
       nk,
       delta,tau)*
      ptau(tau,a,x1,x2,b)*
      pdelta(delta,delta_alpha,delta_beta)
  } # priors of delta
  # loop over chromosomes
  chrs=sort(as.numeric(unique(gent.data$chr)))
  probdf=data.frame()
  for(cc in 1:length(chrs)) {
    if(verbose) cat('starting chromosome',chrs[cc],'\n')
    # filter GenT results to just this chromosome and sort by gene base pair start/end midpoint
    gentres_chr=gent.data %>% filter(chr==chrs[cc]) %>% arrange(mid)
    # LD scores
    lddf_chr=ld.df %>% filter(chr==chrs[cc]) %>% filter(gene %in% gentres_chr$gene)
    ix=sapply(gentres_chr$gene,function(h) which(lddf_chr$gene==h))
    lddf_chr=lddf_chr[ix,]
    PROBS=data.frame()
    # select row of prior distribution parameters
    phi=toadd %>% filter(chr==paste0('chr',chrs[cc]))
    for(i in 1:nrow(gentres_chr)) {
      if(verbose & i%%round(nrow(gentres_chr)*0.1)==0) cat(' ',round(i/nrow(gentres_chr)*100),'% complete\n',sep='')
      # conditional marginal density (integrate over delta and tau)
      prob=pracma::quad2d(f,
                          xa=0, # lower bound of first integration
                          xb=1, # upper bound of first integration
                          ya=phi$tau_param_hard_a, # lower limit of second integration
                          yb=phi$tau_param_hard_b, # upper limit of second integration
                          m=gentres_chr$m[i],
                          gent_sigma2_h0=gentres_chr$gent_sigma2_h0[i],
                          gent_test_statistic=gentres_chr$gent_test_statistic[i],
                          sumwldscore=lddf_chr$sumwldscore[i],
                          sumwld2score=lddf_chr$sumwld2score[i],
                          nk=gwasn,
                          a=phi$tau_param_hard_a,
                          x1=phi$tau_param_x1,
                          x2=phi$tau_param_x1,
                          b=phi$tau_param_hard_b,
                          delta_alpha=phi$delta_param1,
                          delta_beta=phi$delta_param2,
                          n=100)
      PROBS=rbind(PROBS,data.frame(gene=gentres_chr$gene[i],chr=chrs[cc],posterior_probability=prob))
    }
    probdf=rbind(probdf,PROBS %>% mutate(chr=chrs[cc],mid=gentres_chr$mid))
  }
  # done
  return(probdf)
}

#' This function estimates the counts and proportions of independent causal genes which are shared between two traits
#'
#' @param posteriors1 direct output of posterior_gene applied to trait 1
#' @param posteriors2 direct output of posterior_gene applied to trait 2
#' @param gwasn1 GWAS sample size of trait 1
#' @param gwasn2 GWAS sample size of trait 2
#' @param gent.Rho list of matrices. Each matrix in this list corresponds to a chromosome and represents correlations between gene-based association test statistics
#' @param niter  number of iterations for which random sampling of independent genes should be done
#' @param verbose TRUE if progress should be printed, FALSE otherwise
#' @param thr_r2 genes correlated beyond this threshold will be placed in the same block; the thr_r2 argument in bigsnpr::snp_ldsplit
#' @param min_size the minimum size of blocks of genes; the min_size argument in bigsnpr::snp_ldsplit
#' @param max_size the maximum size of blocks of genes; the min_size argument in bigsnpr::snp_ldsplit
#' @param max_K the maximum number of gene blocks; the max_K argument in bigsnpr::snp_ldsplit
#' @param max_r2 the maximum allowable correlation between genes in separate blocks; the max_r2 argument in bigsnpr::snp_ldsplit
#' @keywords shared
#' @export
#' @examples
#' propshared()
propshared=function(posteriors1,posteriors2,gwasn1,gwasn2,gent.Rho,niter=10,verbose=TRUE,
                    thr_r2=0.5,min_size=1,max_size=30,max_K=1e6,max_r2=0.75) {
  # function to estimate proportion of shared genes for each chromosome and
  probdf=data.frame()
  chrs=intersect(unique(posteriors1$chr),unique(posteriors1$chr))
  chrs=sort(as.numeric(chrs))
  blocks_list=list()
  for(cc in 1:length(chrs)) {
    if(verbose) cat('chr',cc,'\n',sep='')
    # merge gent.data1 and gent.data2
    posteriors1_chr=posteriors1 %>% filter(chr==chrs[cc])
    posteriors2_chr=posteriors2 %>% filter(chr==chrs[cc])
    merged_12=inner_join(
      posteriors1_chr %>% select(gene,chr,mid,pp1=posterior_probability),
      posteriors2_chr %>% select(gene,pp2=posterior_probability),by='gene')
    # filter the correlation matrix to just these genes, in this order
    chrmat=gent.Rho[[which(names(gent.Rho)==paste0('chr',chrs[cc]))]]
    merged_12=merged_12 %>% filter(gene %in% rownames(chrmat)) %>% arrange(mid)
    chrmat=chrmat[merged_12$gene,merged_12$gene]
    # identify independent LD blocks (will only be using one matrix)
    Rho=as(chrmat,'sparseMatrix')
    blockres=bigsnpr::snp_ldsplit(Rho,
                                  thr_r2=thr_r2,
                                  min_size=min_size,
                                  max_size=min(c(nrow(Rho),max_size)),
                                  max_K=max_K,
                                  max_r2=max_r2)
    blocks=unlist(tail(blockres$all_last,1))
    if(is.null(blocks)) blocks=1:nrow(chrmat)
    blocks=c(0,blocks)
    # corrplot::corrplot(chrmat[1:50,1:50],method='color') %>% corrRect(blocks+1,col='red')
    # store blocks for the user
    bll=list();for(o in 2:length(blocks)) bll[[o-1]]=merged_12$gene[(1+blocks[o-1]):blocks[o]]
    names(bll)=paste0('block',1:(length(blocks)-1)) # minus 1 because I artifically added the 0
    blocks_list[[cc]]=bll
    names(blocks_list)[cc]=paste0('chr',chrs[cc])
    # iterate over blocks to calculate joint probabilities
    # blockPRODS=blockPROBS1=blockPROBS2=matrix(nr=length(blocks)-1,nc=niter)
    iterdf=data.frame()
    for(iter in 1:niter) {
      ixs=c()
      for(i in 1:(length(blocks)-1)) {
        # sample 1 observation in each correlation block
        blockinds=(blocks[i]+1):blocks[i+1]
        if(length(blockinds)==1) {ixs[i]=blockinds;next}
        ixs[i]=sample(blockinds,1)
      }
      # ixs are indices of independent genes when tested with GenT
      # sum of products of gene-specific posterior probabilities
      prob1=merged_12$pp1[ixs]
      prob2=merged_12$pp2[ixs]
      joint_prob=merged_12$pp1[ixs]*merged_12$pp2[ixs]
      toadd=data.frame(chr=chrs[cc],block=names(bll),iteration=iter,pp1=prob1,pp2=prob2,joint_prob)
      iterdf=rbind(iterdf,toadd)
    }
    # store in larger data frame (running across chromosomes)
    probdf=rbind(probdf,iterdf)
  }
  # calculate estimated counts, shared counts, exclusive counts, and their SEs
  ## trait-specific estimated counts
  estimated_counts1=probdf %>% group_by(iteration) %>% summarise(xbar=sum(pp1))
  estimated_counts2=probdf %>% group_by(iteration) %>% summarise(xbar=sum(pp2))
  v_counts1=probdf %>% group_by(iteration) %>% summarise(vbar=var(pp1))
  v_counts2=probdf %>% group_by(iteration) %>% summarise(vbar=var(pp2))
  impv1=var(estimated_counts1$xbar)*(1+1/niter)+mean(v_counts1$vbar) # imputation variance
  impv2=var(estimated_counts2$xbar)*(1+1/niter)+mean(v_counts2$vbar) # imputation variance
  counts1=mean(estimated_counts1$xbar)
  counts2=mean(estimated_counts2$xbar)
  ## estimated shared counts
  estimated_counts=probdf %>% group_by(iteration) %>% summarise(xbar=sum(joint_prob))
  v_counts=probdf %>% group_by(iteration) %>% summarise(vbar=var(joint_prob))
  joint_var=var(estimated_counts$xbar)*(1+1/niter)+mean(v_counts$vbar) # imputation variance
  joint_counts=mean(estimated_counts$xbar)
  ## estimate unique counts are implied, as are estimated null counts (think of a Venn diagram)
  countdf=data.frame(count_type=c('trait1','trait2','trait1_and_trait2','trait1_or_trait2'),
                     estimated_count=c(counts1,counts2,joint_counts,counts1+counts2-joint_counts),
                     SD=sqrt(c(impv1,impv2,joint_var,impv1+impv2+joint_var)))
  # end
  out=list(
    iteration_data=probdf,
    blocks=blocks_list,
    count_results=countdf
  )
  return(out)
}

#' This function adjusts the raw estimates of shared causal genes between a pair of traits for GWAS sample overlap
#'
#' @param M number of genes tested for association with each trait
#' @param mcausalgenes1 estimated number of genes causing trait 1
#' @param mcausalgenes2 estimated number of genes causing trait 2
#' @param ngwas1 sample size of trait 1 GWAS
#' @param ngwas2 sample size of trait 2 GWAS
#' @param estimated_overlap_count raw data-estimated number of overlapping causal genes
#' @param upsilon_overlap sample overlap correlation coefficient, ie correlation between SNP effect sizes due to sample overlap. See LeBlanc et al. 2015
#' @param h21 estimated/assumed SNP heritability of trait 1
#' @param h22 estimated/assumed SNP heritability of trait 2
#' @param niter number of iterations in each simulation
#' @param m number of tested SNPs per gene
#' @param mcausal_snps_per_gene assumed number of causal SNPs per causal gene
#' @param LD_type structure of correlations between gene-based test statistics due to shared LD of their SNP sets. AR for first order autoregressive, CS for compound symmetry, or I for independence.
#' @param LD_rho correlation parameter specifying the matrix of correlations between gene-based test statistics
#' @param nlambdas number of simulated sample overlap correlation values to evaluate in SIMEX
#' @param doplot TRUE if a plot should be generated, FALSE otherwise
#' @param verbose TRUE if running results should be printed, FALSE otherwise
#' @keywords SIMEX
#' @export
#' @examples
#' shared_count_simex()
shared_count_simex=function(
    M, # number of jointly tested genes
    mcausalgenes1, # estimated number of genes causing trait 1
    mcausalgenes2, # estimated number of genes causing trait 2
    ngwas1, # sample size of trait 1 GWAS
    ngwas2, # sample size of trait 2 GWAS
    estimated_overlap_count, # data-estimated number of overlapping causal genes
    upsilon_overlap, # n01/sqrt(n1*n2)*Corr(t1,t2)
    # assumed simulation parameters
    h21=0.2, # SNP heritability trait 1
    h22=0.2, # SNP heritability trait 1
    niter=1000, # number of iterations
    m=50, # number of tested SNPs per gene
    mcausal_snps_per_gene=3, # number of causal SNPs per causal gene
    LD_type='AR', # type of LD correlation matrix (AR for autoregressive, CS for compound symmetry, I for independence)
    LD_rho=0.5, # correlation parameter
    nlambdas=10, # number of lambdas to evaluate in SIMEX
    doplot=TRUE, # should a plot be created?
    verbose=TRUE
) {
  # SIMEX adjustment to shared counts.
  # this function is used to estimate the number of genes which are inferred to be shared due only to sample overlap
  tr=function(x) sum(diag(x))
  cs=function(n,rho) matrix(rho,n,n)+diag(1-rho,n)
  # genome-wide parameters
  ixs=1:M # gene indices
  delta1=1-mcausalgenes1/M
  delta2=1-mcausalgenes2/M
  ix1=sample(ixs,mcausalgenes1,replace=FALSE) # indices of causal trait 1 genes
  if(estimated_overlap_count==0) {
    ix2=sample(ixs[-ix1],mcausalgenes2,replace=FALSE) # indices of causal trait 2 genes
  } else {
    ix2_a=sample(ix1,estimated_overlap_count,replace=FALSE)
    ix2_b=sample(ixs[-c(ix1,ix2_a)],mcausalgenes2-estimated_overlap_count,replace=FALSE)
    ix2=c(ix2_a,ix2_b)
  }
  # SNP-level parameters
  causalsnpix=round(seq(1,m,length.out=mcausal_snps_per_gene)) # indices of causal SNPs for each gene
  mcausalsnps1=mcausalgenes1*mcausal_snps_per_gene # total number of causal trait 1 SNPs
  mcausalsnps2=mcausalgenes2*mcausal_snps_per_gene # total number of causal trait 2 SNPs
  R=diag(m)
  if(tolower(LD_type)=='ar') R=LD_rho^toeplitz(0:(m-1))
  if(tolower(LD_type)=='cs') R=cs(LD_rho)
  Z0_1=rep(0,m); Z0_1[causalsnpix]=sqrt(ngwas1*h21/mcausalsnps1) # true Z-stats for trait 1 for each causal gene
  Z0_2=rep(0,m); Z0_2[causalsnpix]=sqrt(ngwas2*h22/mcausalsnps2) # true Z-stats for trait 2 for each causal gene
  # moments of null and non-null test statistic distributions
  e0_1=m # null mean of trait 1
  v0_1=2*tr(R%*%R) # null variance of trait 1
  e1_1=m+c(t(Z0_1)%*%Z0_1) # non-null mean of trait 1
  v1_1=v0_1+4*c(t(Z0_1)%*%R%*%Z0_1) # non-null variance of trait 1
  e0_2=m # null mean of trait 2
  v0_2=2*tr(R%*%R) # null variance of trait 2
  e1_2=m+c(t(Z0_2)%*%Z0_2) # non-null mean of trait 2
  v1_2=v0_2+4*c(t(Z0_2)%*%R%*%Z0_2) # non-null variance of trait 2
  # distribution parameters
  beta0_1=e0_1/v0_1 # rate parameter of null distribution for trait 1
  alpha0_1=e0_1*beta0_1 # shape parameter of null distribution for trait 1
  beta1_1=e1_1/v1_1 # rate parameter of non-null distribution for trait 1
  alpha1_1=e1_1*beta1_1 # shape parameter of non-null distribution for trait 1
  beta0_2=e0_2/v0_2 # rate parameter of null distribution for trait 2
  alpha0_2=e0_2*beta0_2 # shape parameter of null distribution for trait 2
  beta1_2=e1_2/v1_2 # rate parameter of non-null distribution for trait 2
  alpha1_2=e1_2*beta1_2 # shape parameter of non-null distribution for trait 2
  # function to perform simulation - defined down here so I don't need to give it any arguments - they come from above
  simres=function(overlap_corr,verbose.=verbose) {
    SigmaOverlap=matrix(c(1,overlap_corr,overlap_corr,1),2,2) # sample overlap correlation matrix
    K=kronecker(SigmaOverlap,R) # variance-covariance matrix of Z-stats for each trait and gene
    RES1=RES2=matrix(nr=niter,nc=M) # store results in these
    for(i in 1:M) {
      Z1=Z2=rep(0,m)
      if(i %in% ix1) Z1[causalsnpix]=sqrt(ngwas1*h21/mcausalsnps1)
      if(i %in% ix2) Z2[causalsnpix]=sqrt(ngwas2*h22/mcausalsnps2)
      z=mvnfast::rmvn(niter,c(Z1,Z2),K)
      T1=rowSums(z[,1:m]^2)
      T2=rowSums(z[,-c(1:m)]^2)
      # log versions to avoid numerical instability
      ## trait 1
      log_f0_1=dgamma(T1,shape=alpha0_1,rate=beta0_1,log=TRUE)
      log_f1_1=dgamma(T1,shape=alpha1_1,rate=beta1_1,log=TRUE)
      log_ratio_1=log_f0_1-log_f1_1
      delta_ratio_1=delta1/(1-delta1)
      p1=1/(1+exp(log_ratio_1)*delta_ratio_1)
      ## trait 2
      log_f0_2=dgamma(T2,shape=alpha0_2,rate=beta0_2,log=TRUE)
      log_f1_2=dgamma(T2,shape=alpha1_2,rate=beta1_2,log=TRUE)
      log_ratio_2=log_f0_2-log_f1_2
      delta_ratio_2=delta2/(1-delta2)
      p2=1/(1+exp(log_ratio_2)*delta_ratio_2)
      RES1[,i]=p1
      RES2[,i]=p2
      if(i %% round(0.1*M)==0 & verbose.) cat(' ',round(i/M*100),'% complete\n',sep='')
    }
    RES=RES1*RES2 # rows are simulations, columns are genes
    rowSums(RES) # niter-length vector of estimated shared counts
  }
  # perform SIMEX
  corrs=seq(upsilon_overlap,sign(upsilon_overlap)*0.99,length.out=nlambdas)
  RESI=matrix(nr=niter,nc=length(corrs))
  for(i in 1:length(corrs)) {
    if(verbose) cat(i,'\n')
    RESI[,i]=simres(corrs[i],verbose.=verbose)
  }
  # end
  x=cbind(corrs,colMeans(RESI))
  colnames(x)=c('sample_overlap_corr','estimated_shared_count')
  # extrapolate
  s=seq(0,upsilon_overlap,length.out=nlambdas)
  fit=lm(x[,2]~x[,1]+I(x[,1]^2))
  ypred=cbind(1,s,s^2)%*%coef(fit)
  if(doplot) {
    plot(c(0,1*sign(upsilon_overlap)),sort(range(c(ypred,x[,2]))),type='n',
         xlab=expression(upsilon*'*'[tt]),ylab='est. shared count')
    points(x[,1],x[,2],type='b',pch=19,col='gray80',cex=2/3)
    lines(s,ypred,col='blue')
    points(0,ypred[which.min(abs(s))],col='blue',pch=19,cex=2/3)
  }
  adj=ypred[1]/x[1,2]
  list(
    adj_shared_count=estimated_overlap_count*adj,
    init_shared_count_est=estimated_overlap_count,
    adj=adj,
    lambda_res=x,
    RES=RESI)
}












