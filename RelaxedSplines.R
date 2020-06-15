# This file contains all the functions necessary to compute splines with level shifts 
# the main function for the fitting is reljumps
# and for computing the fit using the information criteria it is ICjumps


library(HDPenReg)     # c implementation of lars
library(splines)
library(Matrix)

# for plotting 
library(ggplot2)
library(gridExtra)
library(grid)

# for parallelisation
library(foreach)
library(doParallel)
library(doSNOW)
library(iterators)


# function debiasing the first difference term 
pathd1 <- function(X, y, maxsteps){
  
  # initialize
  p1 = dim(X)[2]
  
  # inverse matrix for first differences
  diags = lapply(p1:1, function(i){rep(1,i)})
  Bc    = bandSparse(p1, p1, k = -(0:(p1-1)), diagonals = diags)
 
  # fit model 
  path     = HDlars(as.matrix(X %*% Bc[,-1]), as.vector(y), intercept = FALSE, maxSteps = maxsteps)
  betaCurr = listToMatrix(path)
  beta     = Matrix(0, nrow = p1, ncol = min(maxsteps, length(path@lambda)) + 1)             # matrix of coeffs
  
  #put coeffs into beta dataframe 
  beta[2:p1,] = betaCurr[1:(p1-1),]    # coeffs for jumps
  
  return(list(lambda = path@lambda[1: min(maxsteps, length(path@lambda))], beta = beta, Bc = Bc)) 
}

# function for calculating smoother matrix for given mu
SMats <- function(Xs,P,mu){

  # invert matrix usig qr and if not pd then use svd
  Xinv = tryCatch({
    # invert using qr
    Xinv  =  qr.solve(t(Xs) %*% Xs + mu * t(P) %*% P)     # inverse of matrix using qr decomposition

  }, error = function(e){
  
    eps = 1e-8
    
    svdXs              = svd(t(Xs) %*% Xs + mu * t(P) %*% P)
    d                  = svdXs$d
    dInv               = rep(eps,length(d))
    dInv[abs(d) > eps] = 1/(d[abs(d) > eps])
    Xinv               = svdXs$v %*% Diagonal(length(d),dInv) %*% t(svdXs$u)
  
    })
  
  Smu   = Xinv %*% t(Xs)
  
  
  return(list(Amu = Xs %*% Smu, Smu = Smu))
}


# function calculating the spline penalty matrix and the model matrix XS
XDs <- function(p, xNew, ord){
  
  eps = 1e-10
  
  # initialize 
  p2 = p
  if(ord == 4){ knots = c(min(xNew), min(xNew), min(xNew), seq(min(xNew),max(xNew),length.out = p2-2) , max(xNew), max(xNew), max(xNew)) }
    #knots = c(min(xNew), min(xNew), min(xNew), as.vector(quantile(xNew,probs = seq(0,1,length.out = p2-2))), max(xNew), max(xNew), max(xNew))  }
  if(ord == 2){ knots = c(min(xNew),seq(min(xNew),max(xNew),length.out = p2), max(xNew)) } 
    #knots = c(min(xNew), as.vector(quantile(xNew,probs = seq(0,1,length.out = p2))), max(xNew))}
  
  # create data matrix for smoothing and its penalty matrix
  if(ord == 4){
    diags = list(rep(1,p2-2),rep(-2,p2-2),rep(1,p2-2))
    D2    = bandSparse(p2-2,p2,k = c(0,1,2),diagonals = diags) 
    #D2   = rbind(cbind(Matrix(c(eps, 0, 0, eps), 2, 2),Matrix(0, 2, p2-2)), D2)          # adding small ridge penalty to smoother 
  }else if(ord == 2){
    diags = list(rep(-1,p2-1),rep(1,p2-1))
    D2    = bandSparse(p2-1,p2,k = c(0,1),diagonals = diags) 
    #D2    = rbind(cbind(eps,Matrix(0, 1, p2-1)), D2)          # adding small ridge penalty to smoother 
  }else if(ord == -2){
    diags = list(rep(-1,p2-1),rep(1,p2-1))
    D2    = bandSparse(p2-1,p2,k = c(0,1),diagonals = diags) 
  }else if(ord == -4){
    diags = list(rep(1,p2-2),rep(-2,p2-2),rep(1,p2-2))
    D2    = bandSparse(p2-2,p2,k = c(0,1,2),diagonals = diags) 
  }
  
  # now create spline model matrix 
  if(ord == 2 | ord == 4){
    Xs = splineDesign(knots, xNew, ord = ord, derivs = 0, outer.ok = TRUE, sparse = TRUE)
  }else if(ord == -2 | ord == -4){
    Xs = Diagonal(p2,rep(1,p2))
  }

  
  return(list(Xs = Xs, Ds = D2, knots = knots))
}


# function creating the model matrix for the jumps
Xj <- function(p1, jumps, xNew){
  
  N    = length(xNew)
  
  # data matrix for jumps 
  Xj       = Matrix(0, N, p1)
  zeroCols = vector(mode = "numeric",length = p1)                                        # indices of columns only containing zero's
  for(s in 1:p1){
    if(s != p1){
      currInt = which( xNew >= jumps[s] & xNew < jumps[s+1] )
      if(length(currInt) != 0){ Xj[currInt,s] = 1 }else{ zeroCols[s] = s }
    }else{
      currInt = which( xNew >= jumps[s] & xNew <= jumps[s+1] )
      if(length(currInt) != 0){ Xj[currInt,s] = 1 }
    }
  }
  
  # remove the columns being zero and redifine data matrix and dimension
  zeroCols = zeroCols[zeroCols != 0]
  if(length(zeroCols) != 0){ 
    Xj = Xj[,-zeroCols]
    p1 = p1 -length(zeroCols)
    jumps = jumps[-zeroCols]
  }
  
  
  return(list(Xj = Xj, jumps = jumps))
}

#function creating the model matrix for the jumps for the prediction
Xjpred <- function(p1, jumps, xNew){
  
  N    = length(xNew)
  
  # data matrix for jumps 
  Xj       = Matrix(0, N, p1)
  zeroCols = vector(mode = "numeric",length = p1)                                        # indices of columns only containing zero's
  for(s in 1:p1){
    if(s != p1){
      currInt = which( xNew >= jumps[s] & xNew < jumps[s+1] )
      if(length(currInt) != 0){ Xj[currInt,s] = 1 }else{ zeroCols[s] = s }
    }else{
      currInt = which( xNew >= jumps[s] & xNew <= jumps[s+1] )
      if(length(currInt) != 0){ Xj[currInt,s] = 1 }
    }
  }
  
  return(list(Xj = Xj, jumps = jumps))
}

# function which gets all the model information for all the models at question
getModelsD1 <- function(solPath, p1, H, Amu, X1, yNew){
  
  EPS         = 1e-8
  basisChange =  X1 %*% solPath$Bc 
  
  M = length(solPath$lambda)

  models = lapply(1:M, function(m){
    
    value = solPath$beta[1:p1,m]
    value[abs(value) <= EPS] = 0
    
    # check if there is a jump or not 
    ind = which(value != 0)
    if(m == 1 | (length(ind) != 0 & length(ind) <= H)){  # we always include the first model - with no jumps 
      
      fitC = as.vector(basisChange %*% (solPath$beta[1:p1,m]))
      fitS = as.vector(Amu %*% (yNew - fitC))
      
      return(list(lambda = solPath$lambda[m], fitC = fitC, fitS = fitS, fit = fitS + fitC, ind = ind, indexOfModel = m))
    }
    
  })
  
  models = models[lengths(models) != 0]
  
  indices = lapply(models, function(x) x$ind )
  match( unique(indices),indices ) 
  
  return(models)
}


# function which computes the debiased solution - returns coefficient vector 
debSolD1 <- function(Xmu, ymu, Bc, p1, models, solPath){
  
  # get indices of models 
   M = length(models)
   indices = as.vector(unlist(lapply(1:M,function(m){models[[m]]$indexOfModel})))
   indices  = unique(c(1,indices))
   
   debCoeffs = lapply(1:M,function(m){
   
     if(m == 1){
       
       sol = Bc %*% solPath$beta[,m]
       
     }else{
       
       # compute vectors for checking relaxed lasso properties 
       aMul  = solPath$beta[,indices[m]]
       aMul1 = aMul + solPath$lambda[indices[m]] * (aMul - solPath$beta[,indices[m]-1]) / (solPath$lambda[indices[m]-1] - solPath$lambda[indices[m]] )
       
       if(all(sign(aMul) == sign(aMul1))){  #if relaxed lasso assumptions hold, everything is easy
         
         sol = Bc %*% aMul1

       }else{  #if relaxed lasso assumption does not hold recompute solution 
         
         print("grrr")
         indRemove       = setdiff(1:p1, models[[m]]$ind)              # indices to remove 
         sol             = rep(0,dim(Xmu)[2])
         model           = lm.fit(as.matrix((Xmu %*% Bc)[,-indRemove]), ymu)
         sol[-indRemove] = model$coefficients  
         sol             = Bc %*% sol
         
       }
     }
       
     
     return(sol)
   })
  
  debCoeffs = do.call("cbind", debCoeffs)
  
  return(debCoeffs)
}



# predict function for relJumps object - ind can be a vector indicicating the indices of the models to predict at
predict.jumps <- function(obj, xNew, ind){
  
  EPS    = 1e-8
  p1     = obj$p1
  p      = obj$p
  ord    = obj$ord
  debSol = obj$debSol
  jumps  = obj$jumps
  kinks  = obj$kinks
  knots  = obj$knots

 # 
 # if( (min(xNew) < min(jumps)) |  (max(xNew) > max(jumps))  ){
 #   stop("The new x where to predict must be inside the interval of fit!")
 # }
  
  
  # build model matrices for new data
  X1 = Xjpred(p1, jumps, xNew)$Xj
  if( ord == 4 | ord == 2){
    
    Xs = splineDesign(knots, xNew, ord = ord, derivs = 0, outer.ok = TRUE, sparse = TRUE)
    
    fitJ = X1 %*% debSol[1:p1,ind]
    fitS = Xs %*% debSol[(p1+1):(p1+p),ind]
    
  }else if(ord == -2 | ord == -4){
    
    Xs = Diagonal(p, rep(1,p))
    fitJ = X1 %*% debSol[1:p1,ind]
    
    # per linear interpolation
    fitS = Xs %*% debSol[(p1+1):(p1+p),ind]
    xInd = findInterval(xNew, jumps, left.open = FALSE, rightmost.closed = TRUE)
    fitS = unlist(lapply(1:length(xInd), function(i){
      
      fi  = fitS[xInd[i]]
      fi1 = fitS[xInd[i]+1]
      xi  = jumps[xInd[i]]
      xi1 = jumps[xInd[i]+1]
      
      # calculate value 
      if(abs(xi1 - xi) > EPS){  value = fi * ( xi1 - xNew[i]) / ((xi1 - xi)) + fi1 * (xi - xNew[i]) / (xi - xi1)  }else{ value = fi } 
      
      return(value)
    }))
    
    
    
  }
  

  
  return(list(xNew = xNew, fitJ = fitJ, fitS = fitS, fit = fitJ + fitS))
}


# this function is plotting the fitted function with the original data points - scaled 
# code = jump, smooth, fit ,all 
plot.relJumps <- function(obj, ind, code){
  
  N     = 2000   # for plotting 
  x     = obj$x
  y     = obj$y
  xPlot = seq(min(x), max(x), length.out = N) 
  fPlot = predict.jumps(obj, xPlot, ind)
  fJ    = as.vector(fPlot$fitJ)
  fS    = as.vector(fPlot$fitS)
  f     = as.vector(fPlot$fit)
  
  # construct ggplot objects
  p   = ggplot(data.frame(x = x, y = y), aes(x,y)) + geom_point(alpha = 0.2) + scale_x_continuous(limits = c(min(x), max(x)))
  p   = p + labs(title = "Fitted function", x = "x", y = "y") + theme(axis.title = element_text( size = (15)))
  pj  = ggplot(data.frame(x = xPlot, u = fJ), aes(x, u)) + geom_line() + scale_x_continuous(limits = c(min(x), max(x)))
  pj  = pj + labs(title = "Level shift", x = "x", y = "y") + theme(axis.title = element_text( size = (15)))
  ps  = ggplot(data.frame(x = xPlot, u = fS), aes(x, u)) + geom_line()  + scale_x_continuous(limits = c(min(x), max(x)))
  ps  = ps + labs(title = "Smooth part", x = "x", y = "y") + theme(axis.title = element_text( size = (15)))
  pf  = p + geom_line(data.frame(x = xPlot, u = f), mapping = aes(x, u),col="red")
  #blank  = grid.rect(gp=gpar(fill="white",col="white"),draw = FALSE)
  
  if(code == "jump"){ return(pj) }
  if(code == "smooth"){ return(ps) }
  if(code == "fit"){ return(pf) }
  if(code == "all"){ 
    
    yScalej = c(min(min(y), min(fJ), min(fS), min(f)), max(max(y), max(fJ), max(fS), max(f)))
    pj = pj  #+ scale_y_continuous(limits = yScale)
    pf = pf  #+ scale_y_continuous(limits = yScale)
    ps = ps #+ scale_y_continuous(limits = yScale)
    return(grid.arrange(pf, pj, ps, layout_matrix = matrix(c(1,2,3), 3, 1))) }
  # return(grid.arrange(pf, pj, ps, blank, layout_matrix = matrix(c(1,2,1,2,4,3,4,3), 2, 4)))
}

# function which calculates the coefficients for the debiase model with jumps - ord is the order of the spline (only linear and cubic so far)

# mu is the spline penalty parameter
# p the number of basis functions for the spline
ord = 2, or = 4 if linear or cubic splines desired
# jumps is a vector of possible points where jumps can occur ( usually jumps = x)
# H is a natural number indicating till which step lars should be computed
relJumps <- function(x, y, H, jumps, p = 10, mu, ord){
  
  if(ord == -2 | ord == -4){
    warning("When using Trendfilter then p must be equal to observations")
    p = length(y)
  }
  
  # order x and y accordingly
  yNew = y[order(x)]
  xNew = sort(x)
  N    = length(yNew)
  
  # initialze
  EPS      = 1e-8    # for deciding wheter zero or not
  eps      = 1e-2    # for ridge penalty
  MAXITER  = 200
  MAXSTEPS = H + 1     # for the moment stop after 2 * H steps
  
  # if jumps not given take all the points 
  if(length(jumps) == 0){ jumps = xNew }
  p1  = length(jumps) - 1
  
  # create model matrix and penalty matrix for smooth
  XDsmatrix = XDs(p, xNew, ord)
  Ds        = XDsmatrix$Ds
  Xs        = XDsmatrix$Xs
  
  # data matrix for jumps 
  XjMat = Xj(p1, jumps, xNew)
  X1    = XjMat$Xj
  p1    = dim(X1)[2]
  jumps = XjMat$jumps
  
  # now solve path of the problem with current mu
  ASmats  = SMats(Xs,Ds,mu)
  Amu     = ASmats$Amu
  Smu     = ASmats$Smu
  Xmu     = rbind(X1 - Amu %*% X1,  sqrt(mu) * Ds  %*% (Smu %*% X1) )
  ymu     = c(as.vector(yNew - Amu %*% yNew), as.vector(sqrt(mu) * Ds  %*% (Smu %*%  yNew)))
  solPath = pathd1(Xmu , ymu, maxsteps = MAXSTEPS)
  
  # get the selected models up to H jumps
  models = getModelsD1(solPath, p1, H, Amu, X1, yNew)  

  # compute the debiased solutions 
  Bc     = solPath$Bc
  debSol = debSolD1(Xmu, ymu, Bc, p1, models, solPath)
  debSol = rbind(debSol, Smu %*% (yNew - X1 %*% debSol))
  
  # return relaxed GAM class 
  relaxedJumps = list(x = xNew, y = yNew, jumps = jumps, knots = XDsmatrix$knots,mu = mu, Amu = Amu, Smu = Smu, Xmu = Xmu, Bc = Bc, X1 = X1, Xs = Xs, Ds = Ds, lambda = solPath$lambda, models = models,
                      debSol = debSol, p1 = p1, p = p, ord = ord)
  class(relaxedJumps) = "relaxedJumps"
  return(relaxedJumps)
}


# function for estimation of h - nbr of jumps : information criterion; only for jumps so far
# can also deal with - mus will be for mu values to test and hs for the grid of hs - can also be just one h if known.
# ic = 1 or 2 stands for which mic to use
# p is the order for the spline, ord = 2 or 4 depending if linear or cubic splines chosen
# hs is a vector of number of jumps we want to test
# parallel = TRUE if parallel computation what is wanted
# IC = 1,2 or 3 for desired information criterion
# mus is a vector of possible smoothing parameters
# CLNBR is number of clusters one wants to use

# returns ic scores, the best h (attention h-1 is number of jumps as it starts at h = 1 - no jumps)
# and the best mu found by ic chosen
ICJumps <- function(x, y, jumps, p, ord, mus, hs, parallel = FALSE, ic = 1, gamma = 1.5, CLNBR){
  
  # set cluster number  
  CLNBR = CLNBR
  
  # initialize iterator 
  grid  = expand.grid(mus, hs)            
  param = iter(grid, by = "row")    
  
  if(parallel == TRUE){
    
    # load cluster
    cl = makeCluster(CLNBR) #not to overload your computer
    registerDoSNOW(cl)
    clusterExport(cl = cl,c("relJumps", "XDs", "Xj", "SMats", "pathd1",
                            "getModelsD1", "debSolD1", "predict.jumps",
                            "Xjpred", "mic1", "mic2","mic3"))
    
    ics = foreach(mu = mus, .packages = c('lars', 'HDPenReg', 'splines', 'Matrix', 'glmnet')) %dopar% {
      
      # fit the model for fixed mu and calculate dfs
      model = relJumps(x, y, max(hs)+1, jumps, p, mu, ord)
      dfs   = sum(diag(model$Amu))
      
      icHValues = unlist(lapply(hs, function(h){
        
        fit = predict.jumps(model, x, h)$fit
        dfj = length(model$models[[h]]$ind)
        
        if(ic == 1){ return(mic1(y, fit, dfs, dfj, gamma)) }
        if(ic == 2){ return(mic2(y, fit, dfs, dfj, gamma)) }
        if(ic == 3){ return(mic3(y, fit, dfs, dfj, gamma)) }
      }))
      
      return(icHValues)
    }
    
    stopCluster(cl)  # stopping cluster
    
    # bind all results into matrix together 
    ics = do.call("rbind", ics)
    
    
  }else if(parallel == FALSE){
    
   ics = foreach(mu = mus) %do% {
     
     # fit the model for fixed mu and calculate dfs
     model = relJumps(x, y, max(hs)+1, jumps, p, mu, ord)
     dfs   = sum(diag(model$Amu))
     
     icHValues = unlist(lapply(hs, function(h){
       
       fit = predict.jumps(model, x, h)$fit
       dfj = length(model$models[[h]]$ind)
       
       if(ic == 1){ return(mic1(y, fit, dfs, dfj, gamma)) }
       if(ic == 2){ return(mic2(y, fit, dfs, dfj, gamma)) }
       if(ic == 3){ return(mic3(y, fit, dfs, dfj, gamma)) }
       }))
   
     return(icHValues)
   }
   
   # bind all results into matrix together 
   ics = do.call("rbind", ics)
    
  }

  # get indices where minimal   
  minInd = which(ics == min(ics), arr.ind = TRUE)
  mu     = mus[minInd[1,1]]
  h      = hs[minInd[1,2]]
  
  return(list(ics = ics, h = h, mu =  mu))
}



# information criteria 
mic1 <- function(y, fit, dfs, dfj,gamma){
  
  n     = length(y)
  value = log(sum((y-fit)^2)) - log((1 - gamma * dfs / n)^2) +  dfj * log(n) / n  #3 * dfj *log(n) / n
  
  return(value)
  
  
  
}

mic2 <- function(y, fit, dfs, dfj,gamma){
  
  n     = length(y)
  value = log(sum((y - fit)^2)) -   log((1 -   gamma*dfs/n )^2) +  0.3 * dfj * log(n)^(2.1) / n 

  return(value)
}
# 
 mic3 <- function(y, fit, dfs, dfj,gamma){
   
   n     = length(y)
   value = log(sum((y-fit)^2)) - log((1 - gamma * dfs / n)^2) +  3 * dfj * log(n) / n  #3 * dfj *log(n) / n
   
   return(value)
   
    
   
 }
