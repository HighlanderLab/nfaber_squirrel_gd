##################
## ---- Setup ----
##################

rm(list = ls())
# install.packages(pkg = "ggplot2")
library(ggplot2)
# install.packages(pkg = "reshape2")
library(reshape2)
# install.packages(pkg = "compiler")
library(compiler)
# install.packages(pkg = "iterators")
library(iterators)
# install.packages(pkg = "snow")
library(snow)
# install.packages(pkg = "doSNOW")
library(doSNOW)
# install.packages(pkg = "foreach")
library(foreach)
# install.packages(pkg = "tidyverse")
library(tidyverse)
# install.packages(pkg = "viridis")
library(viridis)

###########################################
## ---- Gene-drive simulation function ----
###########################################

gd <- function(input) {

  ##########################################
  ## ---- Assign all variables in input ----
  ##########################################

  for (i in 1:ncol(input)) {assign(names(input)[i], input[,i])}

  ######################
  ## ---- Functions ----
  ######################

  ## calculate population growth rate from population projection matrix
  r.func <- function(x) log(Re((eigen(x)$values)[1]))

  ## calculate stable age distibution from population projection matrix
  stable.age.dist <- function(x) ((x %*% (Re((eigen(x)$vectors)[,1])))/(sum((x %*% (Re((eigen(x)$vectors)[,1]))))))[,1]

  ## Classify phenotypic sex and reproductive status
  reproductive <- function(pop, strategy) {

    if (nrow(pop)>0) {

      if (strategy==1) {
        ## sex reversal for females, only 1 gene-drive copy needed
        pop$p.sex <- ifelse(pop$g.sex=='f' & (pop$aut.m==-999 | pop$aut.p==-999),'m',pop$g.sex)
        ## only females are made sterile by gene drive
        pop$reprod <- ifelse(pop$g.sex=='f' & (pop$aut.m==-999 | pop$aut.p==-999),0,1)
      } else if (strategy==2) {
        ## sex reversal for females, only 1 gene-drive copy needed
        pop$p.sex <- ifelse(pop$g.sex=='f' & (pop$aut.m==-999 | pop$aut.p==-999),'m',pop$g.sex)
        ## no sterility
        pop$reprod <- 1
      } else if (strategy==3) {
        ## no sex-reversal
        pop$p.sex <- pop$g.sex
        ## no sterility
        pop$reprod <- 1
      } else if (strategy==4) {
        ## no sex-reversal
        pop$p.sex <- pop$g.sex
        ## females homozygous for gene-drive allele are infertile
        if (position=='intron') {
          pop$reprod <- ifelse(pop$g.sex=='f' & pop$aut.m==-999 & pop$aut.p==-999,0,1)
        } else {
          pop$reprod <- ifelse(pop$g.sex=='f' & pop$aut.m.func==0 & pop$aut.p.func==0,0,1)
        }
      } else if (strategy==5) {
        ## no sex-reversal
        pop$p.sex <- pop$g.sex
        ## females heterozygous for gene-drive allele are infertile
        if (position=='intron') {
          pop$reprod <- ifelse(pop$g.sex=='f' & (pop$aut.m==-999 | pop$aut.p==-999),0,1)
        } else {
          pop$reprod <- ifelse(pop$g.sex=='f' & (pop$aut.m.func==0 | pop$aut.p.func==0),0,1)
        }
      }
    }
    pop
  }

  ## Create population
  create <- function(nInit, nInit.vec, nG, Pr){

    pop <- data.frame(id=1:nInit, fid=0, mid=0)
    pop$age <- rep(rep(c(genTime,genTime*2),2),nInit.vec)

    ## genotypic and phenotypic sex
    pop$g.sex <- c('f','m')[rep(1:2,c(sum(nInit.vec[1:2]),sum(nInit.vec[3:4])))] ## genotypic sex
    pop$p.sex <- pop$g.sex ## phenotypic sex

    ## Autosomes
    pop$aut.p <- pop$aut.m <- nGuides ## start with the full complement of cutting sites

    ## Allosomes
    pop$all.m <- 'X' ## maternal allosome (sex chromosome)
    pop$all.p <- ifelse(pop$g.sex=='f','X','Y') ## paternal allosome (sex chromosome)
    
    ## Autosome functionality for strategies 3:5 when positioned in exon (0 - nonfunctional, 1 - functional)
    if (strategy %in% 3:5 & position=='exon') {
      pop$aut.m.func <- ifelse(pop$aut.m<0,0,ifelse(runif(nrow(pop))>(1-pRnonFunc)^(nGuides-pop$aut.m),0,1))
      pop$aut.p.func <- ifelse(pop$aut.p<0,0,ifelse(runif(nrow(pop))>(1-pRnonFunc)^(nGuides-pop$aut.p),0,1))
      pop$aut.m <- ifelse(pop$aut.m%in%-(1:(nGuides+1)),-pop$aut.m-1,pop$aut.m) ## correct back to number of cutting sites
      pop$aut.p <- ifelse(pop$aut.p%in%-(1:(nGuides+1)),-pop$aut.p-1,pop$aut.p) ## correct back to number of cutting sites
    }

    ## Gene drive inoculation (additional to starting population size)
    if (nG>0) {
      pop <- supplementation(pop=pop, amount=nG)
    }

    ## correction for resistant sites through polymorphism
    if (Pr>0) {
      # at this stage alleles are here either wild-type (>0) or gene-drive (-999)
      sel <- pop$aut.m!=-999
      pop$aut.m[sel] <- rbinom(nrow(pop),pop$aut.m[sel],1-Pr)
      sel <- pop$aut.p!=-999
      pop$aut.p[sel] <- rbinom(nrow(pop),pop$aut.p[sel],1-Pr)
    }

    ## Classification of phenotypic sex and reproductive individuals
    pop <- reproductive(pop, strategy=strategy)

    ## return
    pop
  }

  ## Reproduction (breeding, inheritance)
  reproduce <- function(pop, m){

    ## select breeders
    rep.f.ind <- which(pop$p.sex=='f' & pop$reprod==1) ## females
    rep.m.ind <- which(pop$p.sex=='m' & pop$reprod==1) ## males
    
    ## pairing
    if (length(rep.f.ind)>0 & length(rep.m.ind)>0) {
      
      rep.m.ind <- sample(rep.m.ind,length(rep.f.ind),replace=T) ## update males that breed by sampling with replacement
      
      ## final numbers of females and males mating
      (n.rep.f <- length(rep.f.ind))
      (n.rep.m <- length(rep.m.ind))
      n.offspr.vec <- rpois(n.rep.f, m) # how many offspring per reproducing female
      (n.offspr <- sum(n.offspr.vec))

    } else  {
      n.offspr <- 0
    }

    ## births
    if (n.offspr > 0){

      rep.f.df <- offspring <- pop[rep(rep.f.ind,n.offspr.vec),] ## replicating females
      rep.m.df              <- pop[rep(rep.m.ind,n.offspr.vec),] ## replicating males
      f.aut.m.start <- rep.f.df$aut.m ## maternal
      f.aut.p.start <- rep.f.df$aut.p ## paternal
      m.aut.m.start <- rep.m.df$aut.m ## maternal
      m.aut.p.start <- rep.m.df$aut.p ## paternal
      offspring$id <- max(pop$id)+seq_len(n.offspr)
      offspring$age <- 0
      offspring$fid <- rep.f.df$id
      offspring$mid <- rep.m.df$id

      ## Germline homing mothers
      homing.ind.m <- which(rep.f.df$aut.m>0 & rep.f.df$aut.p==-999)
      homing.ind.p <- which(rep.f.df$aut.m==-999 & rep.f.df$aut.p>0)
      
      if(length(homing.ind.m)>0) {
        rmulti <- t(apply(transit[match(rep.f.df$aut.m[homing.ind.m], rownames(transit)),,drop=F],
                          1,function(x) rmultinom(n=1,size=1,prob=x)))
        rep.f.df$aut.m[homing.ind.m] <- as.numeric(colnames(transit)[apply(rmulti,1,function(x) which(x==1))])
      }
      if(length(homing.ind.p)>0) {
        rmulti <- t(apply(transit[match(rep.f.df$aut.p[homing.ind.p], rownames(transit)),,drop=F],
                          1,function(x) rmultinom(n=1,size=1,prob=x)))
        rep.f.df$aut.p[homing.ind.p] <- as.numeric(colnames(transit)[apply(rmulti,1,function(x) which(x==1))])
      }
      
      ## Germline homing fathers
      homing.ind.m <- which(rep.m.df$aut.m>0 & rep.m.df$aut.p==-999)
      homing.ind.p <- which(rep.m.df$aut.m==-999 & rep.m.df$aut.p>0)
      
      if(length(homing.ind.m)>0) {
        rmulti <- t(apply(transit[match(rep.m.df$aut.m[homing.ind.m], rownames(transit)),,drop=F],
                          1,function(x) rmultinom(n=1,size=1,prob=x)))
        rep.m.df$aut.m[homing.ind.m] <- as.numeric(colnames(transit)[apply(rmulti,1,function(x) which(x==1))])
      }
      if(length(homing.ind.p)>0) {
        rmulti <- t(apply(transit[match(rep.m.df$aut.p[homing.ind.p], rownames(transit)),,drop=F],
                          1,function(x) rmultinom(n=1,size=1,prob=x)))
        rep.m.df$aut.p[homing.ind.p] <- as.numeric(colnames(transit)[apply(rmulti,1,function(x) which(x==1))])
      }
      
      ## Autosome and allosome allocation
      aut.m.rand <- runif(n.offspr)
      aut.p.rand <- runif(n.offspr)
      all.m.rand <- runif(n.offspr)
      all.p.rand <- runif(n.offspr)
      offspring$aut.m <- ifelse(aut.m.rand<0.5, rep.f.df$aut.m, rep.f.df$aut.p)
      offspring$aut.p <- ifelse(aut.p.rand<0.5, rep.m.df$aut.m, rep.m.df$aut.p)
      offspring$all.m <- ifelse(all.m.rand<0.5, rep.f.df$all.m, rep.f.df$all.p)
      offspring$all.p <- ifelse(all.p.rand<0.5, rep.m.df$all.m, rep.m.df$all.p) ## Normal paternal allosome allocation

      ## Transfer functionality of autosomes
      ## and recalculate functionality of autosomes for alleles that were functional in the parents
      if (strategy %in% 3:5 & position=='exon') {
        ## transfer
        offspring$aut.m.func <- ifelse(aut.m.rand<0.5, rep.f.df$aut.m.func, rep.f.df$aut.p.func)
        offspring$aut.p.func <- ifelse(aut.p.rand<0.5, rep.m.df$aut.m.func, rep.m.df$aut.p.func)
        aut.m.start <- ifelse(aut.m.rand<0.5, f.aut.m.start, f.aut.p.start)
        aut.p.start <- ifelse(aut.p.rand<0.5, m.aut.m.start, m.aut.p.start)
        aut.m.diff <- aut.m.start-offspring$aut.m ## loss of cutting sites
        aut.p.diff <- aut.p.start-offspring$aut.p ## loss of cutting sites

        ## recalculate
        offspring$aut.m.func <- ifelse(offspring$aut.m<0 | offspring$aut.m.func==0,0,
                                       ifelse(runif(nrow(offspring))>(1-pRnonFunc)^aut.m.diff,0,1))
        offspring$aut.p.func <- ifelse(offspring$aut.p<0 | offspring$aut.p.func==0,0,
                                       ifelse(runif(nrow(offspring))>(1-pRnonFunc)^aut.p.diff,0,1))
      }

      ## convert back to true number of cutting sites
      offspring$aut.m <- ifelse(offspring$aut.m%in%-(1:(nGuides+1)),-offspring$aut.m-1,offspring$aut.m) ## correct back to number of cutting sites
      offspring$aut.p <- ifelse(offspring$aut.p%in%-(1:(nGuides+1)),-offspring$aut.p-1,offspring$aut.p) ## correct back to number of cutting sites

      ## genetic sex allocation
      offspring$g.sex <- ifelse(offspring$all.p%in%'Y','m','f')

      ## if strategy 3, kill homozygous gene drive
      if (strategy==3) {
        if (position=='intron') {
          mort <- which(offspring$aut.m==-999 & offspring$aut.p==-999)
        } else {
          ## if drive is in an exon, kill if no functional allele present
          mort <- which(offspring$aut.m.func==0 & offspring$aut.p.func==0)
        }
        if (length(mort)>0) {
          offspring <- offspring[-mort,]  # kill these ones
        }
      }

      ## if strategy 3, and modelling compensation in utero, correct survival probability based on # implants
      if (strategy==3 & compensation) {
        offspring <- offspring[order(offspring$id),]
        implants <- as.numeric(table(offspring$id))
        implant.s <- implantSmax - (implantSmax-implantSmin)*((implants/m)^implantTheta)
        implant.s.vec <- rep(implant.s, implants)
        survivors <- rbinom(nrow(offspring),1,implant.s.vec)
        offspring <- offspring[survivors==1,]
      }

      ## fix ids and rownames
      offspring$id <- max(pop$id)+seq_len(nrow(offspring))
      rownames(offspring) <- offspring$id

      ## update reproduction status
      offspring <- reproductive(pop=offspring,strategy=strategy)
      
      # add offspring to population
      if (nrow(offspring)>0) {
        pop <- rbind(pop, offspring)
      }
    }
    pop
  }

  ## Mortality
  mortality <- function(pop, s){
    mort <- which(runif(nrow(pop)) > s) # kill these ones
    if (length(mort)>0) {
      pop <- pop[-mort,]
    }
    pop
  }

  # Ageing
  age <- function(pop){
    pop$age <- pop$age + genTime # advance age
    pop
  }

  ## supplementation
  supplementation <- function(pop, amount) {
    supp.pop <- pop[rep(nrow(pop),amount),]
    supp.pop$age <- 2 * genTime
    
    ## Allosomes
    supp.pop$g.sex <- supp.pop$p.sex <- 'm'
    supp.pop$all.m <- 'X'
    supp.pop$all.p <- 'Y'
    
    ## Autosomes
    supp.pop$aut.m <- rbinom(amount,nGuides,1-Pr)
    supp.pop$aut.p <- -999 ## -999 for gene drive
    
    ## Autosome functionality for strategies 3:5,7 when positioned in exon (0 - nonfunctional, 1 - functional)
    if (strategy %in% 3:5 & position=='exon') {
      supp.pop$aut.m.func <- ifelse(supp.pop$aut.m<0,0,ifelse(runif(nrow(supp.pop))>(1-pRnonFunc)^(nGuides-supp.pop$aut.m),0,1))
      supp.pop$aut.p.func <- ifelse(supp.pop$aut.p<0,0,ifelse(runif(nrow(supp.pop))>(1-pRnonFunc)^(nGuides-supp.pop$aut.p),0,1))
      supp.pop$aut.m <- ifelse(supp.pop$aut.m%in%-(1:(nGuides+1)),-supp.pop$aut.m-1,supp.pop$aut.m) ## correct back to number of cutting sites
      supp.pop$aut.p <- ifelse(supp.pop$aut.p%in%-(1:(nGuides+1)),-supp.pop$aut.p-1,supp.pop$aut.p) ## correct back to number of cutting sites
    }
    
    supp.pop$id <- tail(pop$id,1)+seq_len(amount)
    rownames(supp.pop) <- supp.pop$id
    pop <- rbind(pop, supp.pop)
    
    pop
  }
  
  #########################################################
  ## ---- Transition matrix for germline homing events ----
  #########################################################

  GeneDriveSeqRec <- function(S, Pc, Pn) {

    mat <- matrix(0, S+1, S+1)
    rownames(mat) <- colnames(mat) <- 0:S
    diag(mat) <- sapply(0:S, function(v) (1-Pc)^v)
    n.vec <- unlist(apply(lower.tri(mat)[-1,,drop=F],2,function(x) which(x))) ## cutting sites remaining
    i.vec <- rep(0:(S-1), apply(lower.tri(mat)[-1,,drop=F],2,sum)[-(S+1)]) ## number of resulting sites
    mat[lower.tri(mat)] <- sapply(1:length(n.vec), function(x) choose(n.vec[x],i.vec[x])*((1-Pc)^i.vec[x])*(Pc^(n.vec[x]-i.vec[x]))*Pn^(n.vec[x]-i.vec[x]))

    ## add gene drive column
    mat <- cbind(mat, rep(0,S+1))
    colnames(mat)[ncol(mat)] <- -999
    n.vec2 <- 0:S ## cutting sites remaining
    mat[,'-999'] <- sapply(0:S,function(n) sum(sapply(1:n, function(i) choose(n,i)*(Pc^i)*((1-Pc)^(n-i))*(1-(Pn^i)))))
  
    ## round matrix to 15 digits to avoid finite floating-point precision problems
    mat <- round(mat,digits=15)
    
    return(mat)
  }

  GeneDriveSimRec <- function(S,Pc,Pn) {

    mat <- matrix(0,S,S+2)
    mat[1,1] <- Pc*Pn
    mat[1,2] <- (1-Pc)

    if (S>1) { # if there is more than 1 gRNA
      for (nS in 2:S) { # for each gRNA starting number
        mat[nS,2:(nS+1)] <- (1-Pc)*mat[nS-1,1:nS]
        mat[nS,1:nS] <- mat[nS,1:nS] + Pc*Pn*mat[nS-1,1:nS]
        for (ni in 2:nS) { # for each gRNA ending number
          prob <- 0 # the probability that more than one gRNA was cut at the same time
          for (j in 2:ni) { # for each amount of gRNA's that could be cut at the same time
            prob <- prob + choose(ni-2,j-2)*Pc^j*(1-Pc)^(ni-j)*Pn*(1-Pn)^j
          }
          for (nl in seq(nS-ni,0,by=-1)) {
            if ((nS-ni)>0) {
              mat[nS,nS-ni-nl+1] <- mat[nS,nS-ni-nl+1] + prob*mat[nS-ni,nS-ni-nl+1]
            } else {
              mat[nS,1] = mat[nS,1] + prob
            }
          }
        }
      }
    }
    mat[,S+2] = 1 - rowSums(mat[,1:(S+1),drop=F])
    mat <- rbind(c(1,rep(0,S+1)), mat)
    rownames(mat) <- 0:S
    colnames(mat) <- c(0:S,'-999')
    
    ## round matrix to 15 digits to avoid finite floating-point precision problems
    mat <- round(mat,digits=15)
    
    return(mat)
  }

  ## create transition matrix depending on whether simultaneous/sequential cutting and position of gene drive
  ## when transType=='seq', intervening sequence deletion is indicated by -(s + 1)
  if (transType=='seq') {
    transit <- GeneDriveSeqRec(S=nGuides,Pc=Pc,Pn=Pn)
  } else if (transType=='sim' & position=='intron') {
    transit <- GeneDriveSimRec(S=nGuides,Pc=Pc,Pn=Pn)
  } else {
    transit.seq <- GeneDriveSeqRec(S=nGuides,Pc=Pc,Pn=Pn)
    transit.sim <- GeneDriveSimRec(S=nGuides,Pc=Pc,Pn=Pn)
    diff.mat <- (transit.sim-transit.seq)[,-(nGuides+2)]
    colnames(diff.mat) <- -(as.numeric(colnames(diff.mat))+1)
    transit <- cbind(diff.mat,transit.seq)
    transit[,'-999'] <- 1-rowSums(transit[,-ncol(transit)])
  }

  ###########################################################################
  ### ---- Survival for stable population and starting population vector ----
  ###########################################################################

  ## calculate fertility rate for r=0
  r0.func <- function(s) {
    a <- matrix(0,4,4)
    a[1,1:2] <- s*m*0.5
    a[3,1:2] <- s*m*0.5
    a[2,1:2] <- a[4,3:4] <- s
    minim <- abs(r.func(a))
    return(minim)
  }
  (s <- optimise(r0.func,interval=c(1e-10,100))$minimum)

  ## get stable age distribution
  a <- matrix(0,4,4)
  a[1,1:2] <- s*m*0.5
  a[3,1:2] <- s*m*0.5
  a[2,1:2] <- a[4,3:4] <- s
  (sad <- stable.age.dist(a))

  ### adjust rmax for generation time
  (rmax.adj <- rmax*(genTime/52))

  ## calculate survival rate for r=rmax
  rmax.func <- function(s) {
    a <- matrix(0,4,4)
    a[1,1:2] <- s*m*0.5
    a[3,1:2] <- s*m*0.5
    a[2,1:2] <- a[4,3:4] <- s
    minim <- abs(r.func(a)-rmax.adj)
    return(minim)
  }
  (s.max <- optimise(rmax.func,interval=c(1e-10,100))$minimum)

  ## intercept and slope of logistic survival equation
  (alpha.s <- -log((1/s.max)-1))
  (beta.s <- -alpha.s - log((1/s)-1))

  ## starting n vector
  nInit <- K
  (nInit.vec <- round(sad*nInit))

  #######################
  ## ---- Simulation ----
  #######################
  
  # Object for storing results (population size in time)
  year <- seq(0,nYear,by=genTime/52)
  gen <- 0:(nYear*(52/genTime))
  nT <- length(gen)-1
  var <- c('nInd', 'nIndF', 'nIndM',
           'nIndWW', 'nIndWG', 'nIndWN', 'nIndGG', 'nIndGN', 'nIndNN',
           'pIndWW', 'pIndWG', 'pIndWN', 'pIndGG', 'pIndGN', 'pIndNN',
           'nIndWWF', 'nIndWGF', 'nIndWNF', 'nIndGGF', 'nIndGNF', 'nIndNNF',
           'pIndWWF', 'pIndWGF', 'pIndWNF', 'pIndGGF', 'pIndGNF', 'pIndNNF',
           'nIndWWM', 'nIndWGM', 'nIndWNM', 'nIndGGM', 'nIndGNM', 'nIndNNM',
           'pIndWWM', 'pIndWGM', 'pIndWNM', 'pIndGGM', 'pIndGNM', 'pIndNNM',
           'nAllW', 'nAllG', 'nAllN',
           'pAllW', 'pAllG', 'pAllN',
           'nAllWF', 'nAllGF', 'nAllNF',
           'pAllWF', 'pAllGF', 'pAllNF',
           'nAllWM', 'nAllGM', 'nAllNM',
           'pAllWM', 'pAllGM', 'pAllNM')
  results <- array(NA,
                   dim=c(length(year),length(var),nIter),
                   dimnames=list(year=year,var=var,iter=1:nIter))
  results[1,,] <- 0
  results[1,c('nInd','nIndF','nIndM'),] <- c(nInit+nG, (nInit/2), (nInit/2)+nG)

  # Run
  for (iter in 1:nIter) {
    # initial population
    nInd <- nInit+nG
    pop <- create(nInit=nInit, nInit.vec=nInit.vec, nG=nG, Pr=Pr)
    
    for(t in 2:(nT+1)) {

      # supplementation each year
      if (supplement > 0 & floor(nrow(pop)*supplement) > 0){
        decTime <- year[t]
        if (round(decTime/suppInterval,10)==ceiling(round(decTime/suppInterval,10))) { 
          amount <- floor(nrow(pop)*supplement)
          pop <- supplementation(pop=pop, amount=amount)
        }
      }
      
      # demography
      pop <- reproduce(pop, m=m)
      (N.now <- results[match(floor(gen[t-1]), gen),'nInd',iter])
      (s.now <- ifelse(is.na(rmax),s,1/(1+exp(-(alpha.s+beta.s*N.now/K)))))
      pop <- mortality(pop, s=s.now)
      pop <- age(pop)

      ## calculate summaries
      (nInd <- nrow(pop))

      femSel <- pop$p.sex == 'f'
      malSel <- pop$p.sex == 'm'
      nIndF <- sum(femSel)
      nIndM <- sum(malSel)
      
      # gene drive resistance
      aut.mG <- pop$aut.m == -999 # gene drive
      aut.mN <- pop$aut.m == 0    # resitant
      aut.mW <- !aut.mG | !aut.mN # wild type

      aut.pG <- pop$aut.p == -999 # gene drive
      aut.pN <- pop$aut.p == 0    # resitant
      aut.pW <- !aut.pG | !aut.pN # wild type

      nIndWW <- sum(aut.mW & aut.pW)
      nIndWG <- sum((aut.mW & aut.pG) | (aut.mG & aut.pW))
      nIndWN <- sum((aut.mW & aut.pN) | (aut.mN & aut.pW))
      nIndGG <- sum(aut.mG & aut.pG)
      nIndGN <- sum((aut.mG & aut.pN) | (aut.mN & aut.pG))
      nIndNN <- sum(aut.mN & aut.pN)

      pIndWW <- nIndWW / nInd
      pIndWG <- nIndWG / nInd
      pIndWN <- nIndWN / nInd
      pIndGG <- nIndGG / nInd
      pIndGN <- nIndGN / nInd
      pIndNN <- nIndNN / nInd

      nIndWWF <- sum(aut.mW[femSel] & aut.pW[femSel])
      nIndWGF <- sum((aut.mW[femSel] & aut.pG[femSel]) | (aut.mG[femSel] & aut.pW[femSel]))
      nIndWNF <- sum((aut.mW[femSel] & aut.pN[femSel]) | (aut.mN[femSel] & aut.pW[femSel]))
      nIndGGF <- sum(aut.mG[femSel] & aut.pG[femSel])
      nIndGNF <- sum((aut.mG[femSel] & aut.pN[femSel]) | (aut.mN[femSel] & aut.pG[femSel]))
      nIndNNF <- sum(aut.mN[femSel] & aut.pN[femSel])

      pIndWWF <- nIndWWF / nIndF
      pIndWGF <- nIndWGF / nIndF
      pIndWNF <- nIndWNF / nIndF
      pIndGGF <- nIndGGF / nIndF
      pIndGNF <- nIndGNF / nIndF
      pIndNNF <- nIndNNF / nIndF

      nIndWWM <- sum(aut.mW[malSel] & aut.pW[malSel])
      nIndWGM <- sum((aut.mW[malSel] & aut.pG[malSel]) | (aut.mG[malSel] & aut.pW[malSel]))
      nIndWNM <- sum((aut.mW[malSel] & aut.pN[malSel]) | (aut.mN[malSel] & aut.pW[malSel]))
      nIndGGM <- sum(aut.mG[malSel] & aut.pG[malSel])
      nIndGNM <- sum((aut.mG[malSel] & aut.pN[malSel]) | (aut.mN[malSel] & aut.pG[malSel]))
      nIndNNM <- sum(aut.mN[malSel] & aut.pN[malSel])

      pIndWWM <- nIndWWM / nIndM
      pIndWGM <- nIndWGM / nIndM
      pIndWNM <- nIndWNM / nIndM
      pIndGGM <- nIndGGM / nIndM
      pIndGNM <- nIndGNM / nIndM
      pIndNNM <- nIndNNM / nIndM

      nAllW <- nIndWW*2 + nIndWG + nIndWN
      nAllG <- nIndGG*2 + nIndWG + nIndGN
      nAllN <- nIndNN*2 + nIndWN + nIndGN
      nAll <- nAllW + nAllG + nAllN

      pAllW <- nAllW / nAll
      pAllG <- nAllG / nAll
      pAllN <- nAllN / nAll
      
      nAllWF <- nIndWWF*2 + nIndWGF + nIndWNF
      nAllGF <- nIndGGF*2 + nIndWGF + nIndGNF
      nAllNF <- nIndNNF*2 + nIndWNF + nIndGNF
      nAllF <- nAllWF + nAllGF + nAllNF

      pAllWF <- nAllWF / nAllF
      pAllGF <- nAllGF / nAllF
      pAllNF <- nAllNF / nAllF

      nAllWM <- nIndWWM*2 + nIndWGM + nIndWNM
      nAllGM <- nIndGGM*2 + nIndWGM + nIndGNM
      nAllNM <- nIndNNM*2 + nIndWNM + nIndGNM
      nAllM <- nAllWM + nAllGM + nAllNM

      pAllWM <- nAllWM / nAllM
      pAllGM <- nAllGM / nAllM
      pAllNM <- nAllNM / nAllM

      # store results (just for each whole generation time)
      for (v in var) {
        results[t,v,iter] <- get(x=v)
      }
      
      # break when one sex remains or N = 0
      if (nInd==0 | nIndF==0 | nIndM==0 ) {
        break
      }

    }

    print(iter)
  }

  ## melt into required format
  res <- melt(results[,'nInd',],varnames=c('year','iter'),value.name='nInd')
  for (v in var[-1]) {
    tmp <- paste(v, ".res", sep="")
    assign(x=tmp, value=melt(results[,v,],varnames=c('year','iter'),value.name=v))
    res <- cbind(res,get(x=tmp)[,v,drop=F])
  }

  ## add all the parameters to the results data frame too
  for (i in 1:ncol(input)) {
    res[[names(input)[i]]] <- input[,i]
  }
  return(res)
}

############################################
## ---- Compile the simulation function ----
############################################

gd.comp <- cmpfun(gd)
