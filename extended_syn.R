# library entire packages
library(rpart)
library(MASS)
library(methods)
library(nnet)
library(lattice)
library(foreign)
library(ggplot2)
library(plyr)
library(proto)
library(survival)
library(stringr)
# Functions for synthesising data, some of which are adapted
# from mice package by S. van Buuren and K. Groothuis-Oudshoorn,
# TNO Quality of Life


###-----.norm.fix.syn------------------------------------------------------

.norm.fix.syn <- function(y, x, ridge = 0.00001, ...)
{
  # Calculates regression coefficients + error estimate
  
  xtx <- t(x) %*% x
  pen <- ridge * diag(xtx)
  if (length(pen)==1) pen <- matrix(pen)
  v           <- solve(xtx + diag(pen))
  coef        <- t(y %*% x %*% v)
  residuals   <- y - x %*% coef
  sigma       <- sqrt((sum(residuals^2))/(length(y)-ncol(x)-1))
  parm        <- list(coef, sigma)
  names(parm) <- c("beta","sigma")
  return(parm)
}


###-----.norm.draw.syn-----------------------------------------------------

.norm.draw.syn <- function(y, x, ridge = 0.00001, ...)
{
  # Draws values of beta and sigma for Bayesian linear regression synthesis 
  # of y given x according to Rubin p.167
  
  xtx <- t(x) %*% x
  pen <- ridge * diag(xtx)
  if (length(pen)==1) pen <- matrix(pen)
  v           <- solve(xtx + diag(pen))
  coef        <- t(y %*% x %*% v)
  residuals   <- y - x %*% coef
  sigma.star  <- sqrt(sum((residuals)^2)/rchisq(1, length(y) - ncol(x)))
  beta.star   <- coef + (t(chol((v + t(v))/2)) %*% rnorm(ncol(x))) * sigma.star
  parm        <- list(coef, beta.star, sigma.star)      
  names(parm) <- c("coef","beta","sigma")      
  return(parm)
}


###-----syn.norm-----------------------------------------------------------

syn.norm <- function(y, x, xp, proper = FALSE, ...)
{
  x   <- cbind(1, as.matrix(x))
  xp  <- cbind(1, as.matrix(xp))
  if (proper == FALSE){
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }  
  res <- xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma
  res <- round(res, max(sapply(y, decimalplaces)))
  
  return(list(res = res, fit = parm))
}


###-----syn.lognorm--------------------------------------------------------

syn.lognorm <- function(y, x, xp, proper = FALSE, ...) 
{
  addbit <- FALSE
  if (any(y < 0)) stop("Log transformation not appropriate for negative values.\n", call. = FALSE)
  if (any(y == 0)) {y <- y + .5*min(y[y != 0]); y <- log(y); addbit <- TRUE}  ##  warning about this and above should be in check model
  else y <- log(y)
  x   <- cbind(1, as.matrix(x))
  xp  <- cbind(1, as.matrix(xp))
  if (proper == FALSE) {
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }
  res <- xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma
  if (addbit) {res <- res - .5 * min(y[y != 0]); res[res <= 0] <- 0}
  res <- exp(res)
  res <- round(res, max(sapply(y, decimalplaces)))
  
  return(list(res = res, fit = parm))
}


###-----syn.sqrtnorm-------------------------------------------------------

syn.sqrtnorm <- function(y, x, xp, proper = FALSE, ...) 
{
  addbit <- FALSE
  if (any(y < 0)) stop("Square root transformation not appropriate for negative values.\n", call. = FALSE)   ##  needs check in checkmodel
  else y <- sqrt(y)
  x   <- cbind(1, as.matrix(x))
  xp  <- cbind(1, as.matrix(xp))
  if (proper == FALSE) {
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }
  res <- xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma
  res <- res^2
  res <- round(res, max(sapply(y, decimalplaces)))
  
  return(list(res = res, fit = parm))
}


###-----syn.cubertnorm-----------------------------------------------------

syn.cubertnorm <- function(y, x, xp, proper = FALSE, ...) 
{
  addbit <- FALSE
  y <- sign(y)*abs(y)^(1/3)
  x   <- cbind(1, as.matrix(x))
  xp  <- cbind(1, as.matrix(xp))
  if (proper == FALSE) {
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }
  res <- xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma
  res <- res^3
  res <- round(res, max(sapply(y, decimalplaces)))
  
  return(list(res = res, fit = parm))
}


###-----syn.normrank-------------------------------------------------------

syn.normrank <- function(y, x, xp, smoothing = "", proper = FALSE, ...)
{
  # Regression synthesis of y given x, with a fixed regression
  # line, and with random draws of the residuals around the line.
  # Adapted from norm by carrying out regression on Z scores from ranks
  # predicting new Z scores and then transforming back
  # similar to method by ? and ?
  #
  # First get approx rank position of vector in one of another length
  # so that result returned has correct length for xp
  # matters for sub-samples and missing data
  
  z  <- qnorm(rank(y)/(length(y) + 1))
  x  <- cbind(1, as.matrix(x))
  xp <- cbind(1, as.matrix(xp))
  
  if (proper == FALSE) {
    parm <- .norm.fix.syn(z, x, ...)
  } else {
    parm <- .norm.draw.syn(z, x, ...)
  }
  
  pred <- (xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma)
  res  <- round(pnorm(pred)*(length(y) + 1))
  res[res < 1] <- 1
  res[res > length(y)] <- length(y)
  res  <- sort(y)[res]
  
  if (smoothing != "") {
    res <- syn.smooth(res, y, smoothing = smoothing)
  }
  
  # if (smoothing == "") res  <- sort(y)[res]
  # 
  # if (smoothing == "density") {
  #   ydsamp <- y
  #   ys     <- 1:length(y)
  #   maxfreq <- which.max(table(y))
  #   maxcat  <- as.numeric(names(table(y))[maxfreq])
  #   if (table(y)[maxfreq]/sum(table(y)) > .7) ys <- which(y != maxcat)
  #   if (10 * table(y)[length(table(y)) - 1] < 
  #     tail(table(y), n = 1) - table(y)[length(table(y)) - 1]) {
  #     ys <- ys[-which(y == max(y))]  
  #     maxy <- max(y)
  #   }   
  #   densbw <- density(y[ys], width = "SJ")$bw
  #   ydsamp[ys] <- rnorm(length(ydsamp[ys]), 
  #     mean = sample(ydsamp[ys], length(ydsamp[ys]), replace = TRUE), sd = densbw)
  #   if (!exists("maxy")) maxy <- max(y) + densbw
  #   ydsamp[ys] <- pmax(pmin(ydsamp[ys],maxy),min(y))
  #   res <- sort(ydsamp)[res]
  # }
  
  return(list(res = res, fit = parm))
}


###-----.pmm.match---------------------------------------------------------

.pmm.match <- function(z, yhat = yhat, y = y, donors = 3, ...)
{
  # Auxilary function for syn.pmm.
  # z    = target predicted value (scalar)
  # yhat = array of fitted values, to be matched against z
  # y    = array of donor data values
  
  # Finds the three cases for which abs(yhat-z) is minimal,
  # and makes a random draw from these.
  
  d <- abs(yhat - z)
  m <- sample(y[rank(d, ties.method = "random") <= donors], 1)
  
  return(m)
}


###-----syn.pmm------------------------------------------------------------

syn.pmm <- function(y, x, xp, smoothing = "", proper = FALSE, ...)
{
  # Synthesis of y by predictive mean matching
  # Warning: can be slow for large data sets 
  # for which syn.normrank may be a better choice
  x       <- cbind(1, as.matrix(x))
  xp      <- cbind(1, as.matrix(xp))
  if (proper == FALSE) {
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }
  yhatobs <- x  %*% parm$coef
  yhatmis <- xp %*% parm$beta
  res <- apply(as.array(yhatmis), 1, .pmm.match, yhat = yhatobs, y = y, ...)
  
  if (smoothing != "") {
    res <- syn.smooth(res, y, smoothing = smoothing)
  }
  
  return(list(res = res, fit = parm))
}


###-----augment.syn--------------------------------------------------------

augment.syn <- function(y, x, ...)
{
  # define augmented data for stabilizing logreg and polyreg
  # by the ad hoc procedure of White, Daniel & Royston, CSDA, 2010
  # This function will prevent augmented data beyond the min and
  # the max of the data
  # Input:
  # x: numeric data.frame (n rows)
  # y: factor or numeric vector (length n)
  # Output:
  # return a list with elements y, x, and w with length n+2*(ncol(x))*length(levels(y))
  
  x    <- as.data.frame(x)
  icod <- sort(unique(unclass(y)))
  ki   <- length(icod)
  # if (ki>maxcat) stop(paste("Maximum number of categories (",maxcat,") exceeded", sep=""))
  p    <- ncol(x)
  
  # skip augmentation if there are no predictors
  if (p == 0) return(list(y = y, x = x, w = rep(1, length(y))))
  
  # skip augmentation if there is only 1 missing value  
  if (length(y) == 1) return(list(y = y, x = x, w = rep(1, length(y))))
  
  # calculate values to augment
  mean <- apply(x,2,mean)
  sd   <- sqrt(apply(x,2,var))
  minx <- apply(x,2,min)
  maxx <- apply(x,2,max)
  nr   <- 2 * p * ki
  a    <- matrix(mean, nrow = nr, ncol = p, byrow = TRUE)
  b    <- matrix(rep(c(rep(c(0.5, -0.5), ki), rep(0,nr)), length = nr*p), 
                 nrow = nr, ncol = p, byrow = FALSE)
  c    <- matrix(sd, nrow = nr, ncol = p, byrow = TRUE)
  d    <- a + b * c
  d    <- pmax(matrix(minx, nrow = nr, ncol = p, byrow = TRUE), d)
  d    <- pmin(matrix(maxx, nrow = nr, ncol = p, byrow = TRUE), d)
  e    <- rep(rep(icod, each = 2), p)
  dimnames(d) <- list(paste("AUG", 1:nrow(d), sep = ""), dimnames(x)[[2]])
  xa   <- rbind.data.frame(x, d)
  # beware, concatenation of factors
  # this change needed to avoid reordering of factors                           
  # if (is.factor(y)) ya <- as.factor(levels(y)[c(y,e)]) else ya  <- c(y, e)
  if (is.factor(y)) ya <- addNA(factor(levels(y)[c(y, e)], 
                                       levels = levels(y)), ifany = TRUE) else ya <- c(y, e)   
  wa <- c(rep(1, length(y)),rep((p + 1)/nr, nr))
  
  return(list(y = ya, x = xa, w = wa))
}


###-----syn.logreg---------------------------------------------------------

syn.logreg <- function(y, x, xp, denom = NULL, denomp = NULL, 
                       proper = FALSE, ...)            
{
  # Synthesis for binary or binomial response variables by
  # logistic regression model. See Rubin (1987, p. 169-170) for
  # a description of the method.
  
  # The method consists of the following steps:
  # 1. Fit a logit, and find (bhat, V(bhat))
  # 2. Draw BETA from N(bhat, V(bhat))
  # 3. Compute predicted scores for m.d., i.e. logit-1(X BETA)
  # 4. Compare the score to a random (0,1) deviate, and synthesise.
  
  xmeans <- lapply(x, mean)                      ## x matrix centred
  x  <- mapply(function(x, y) x - y, x, xmeans)
  xp <- mapply(function(x, y) x - y, xp, xmeans) ## also xp to match
  
  if (is.null(denom)) {
    aug <- augment.syn(y, x, ...)
    # when no missing data must set xf to augmented version
    xf   <- aug$x
    y    <- aug$y
    w    <- aug$w
    xf   <- cbind(1, as.matrix(xf))
    xp   <- cbind(1, as.matrix(xp))
    expr <- expression(glm.fit(xf, y, family = binomial(link = logit), weights = w))
    fit  <- suppressWarnings(eval(expr))
    fit.sum <- summary.glm(fit)
    beta <- coef(fit)
    if (proper == TRUE) {
      rv   <- t(chol(fit.sum$cov.unscaled))
      beta <- beta + rv %*% rnorm(ncol(rv))  
    }
    p   <- 1/(1 + exp(-(xp %*% beta)))  
    vec <- (runif(nrow(p)) <= p)
    if (!is.logical(y)) vec <- as.numeric(vec)          
    if (is.factor(y)) vec <- factor(vec,c(0,1), labels = levels(y))
  } else {
    aug <- augment.syn(y, x, ...)
    # when no missing data must set xf to augmented version
    xf   <- aug$x
    y    <- aug$y
    w    <- aug$w
    xf   <- cbind(1, as.matrix(xf))
    xp   <- cbind(1, as.matrix(xp))
    den  <- w
    denind <- which(den == 1)
    den[denind] <- denom
    yy   <- y/den        #denom give then average response
    yy[den < 1]   <- mean(yy[denind]) 
    expr <- expression(glm.fit(xf, yy, family = binomial(link = logit), weights = den))
    fit  <- suppressWarnings(eval(expr))
    fit.sum <- summary.glm(fit)
    beta <- coef(fit.sum)[, 1]
    if (proper == TRUE) {
      rv   <- t(chol(fit.sum$cov.unscaled))
      beta <- beta + rv %*% rnorm(ncol(rv))  
    }
    p <- 1/(1 + exp(-(xp %*% beta)))  
    vec <- rbinom(nrow(p),denomp, p) 
  }
  return(list(res = vec, fit = fit.sum))
}


###-----syn.polyreg--------------------------------------------------------   

syn.polyreg <- function(y, x, xp, proper = FALSE, maxit = 1000, 
                        trace = FALSE, MaxNWts = 10000, ...)
{
  # synthesis for categorical response variables by the Bayesian
  # polytomous regression model. See J.P.L. Brand (1999), Chapter 4,
  # Appendix B.
  #
  # The method consists of the following steps:
  # 1. Fit categorical response as a multinomial model
  # 2. Compute predicted categories
  # 3. Add appropriate noise to predictions.
  #
  # This algorithm uses the function multinom from the libraries nnet and MASS
  # (Venables and Ripley).
  
  x   <- as.matrix(x)
  xp  <- as.matrix(xp)
  
  if (proper == TRUE) { # bootstrap to make proper
    s   <- sample(length(y), replace = TRUE)
    x   <- x[s, , drop = FALSE]
    y   <- y[s]  
    y   <- factor(y)
  }
  aug <- augment.syn(y, x, ...)
  # yf and xf needed for augmented data to save x as non augmented  not now needed can tidy
  xf  <- aug$x
  yf  <- aug$y
  w   <- aug$w
  
  ### rescaling numeric to [0,1]
  toscale <- sapply(xf, function(z) (is.numeric(z) & (any(z < 0) | any(z > 1))))
  rsc <- sapply(xf[, toscale, drop = FALSE], range)
  xf_sc <- xf
  for (i in names(toscale[toscale == TRUE])) xf_sc[, i] <- (xf_sc[, i] - rsc[1,i])/(rsc[2,i] - rsc[1,i])
  for (i in names(toscale[toscale == TRUE])) xp[, i] <- (xp[, i] - rsc[1,i])/(rsc[2,i] - rsc[1,i])
  ###
  
  xfy <- cbind.data.frame(yf, xf_sc)  
  fit <- multinom(formula(xfy), data = xfy, weights = w,
                  maxit = maxit, trace = trace, MaxNWts = MaxNWts, ...)
  if (fit$convergence == 1) cat("\nReached max number of iterations for a multinomial model\nsuggest rerunning with polyreg.maxit increased (default 1000)\n")             
  post <- predict(fit, xp, type = "probs") 
  if (length(y) == 1) post <- matrix(post, nrow = 1, ncol = length(post)) 
  if (!is.factor(y)) y <- as.factor(y)
  nc <- length(levels(yf))                    
  un <- rep(runif(nrow(xp)), each = nc)
  if (is.vector(post)) post <- matrix(c(1 - post, post), ncol = 2)
  draws <- un > apply(post, 1, cumsum)
  idx   <- 1 + apply(draws, 2, sum)
  res <- levels(yf)[idx]
  if (length(table(res)) == 1) {
    cat("\n***************************************************************************************")
    cat("\nWarning the polyreg fit produces only one category for the variable being synthesised." )
    cat("\nThis may indicate that the function multinom used in polyreg failed to iterate, possibly")
    cat("\nbecause the variable is sparse. Check results for this variable carefully.")
    cat("\n****************************************************************************************\n")
  } 
  fitted <- summary(fit)
  return(list(res = res, fit = fitted)) 
}


###-----syn.polr-----------------------------------------------------------

syn.polr <- function(y, x, xp, proper = FALSE, maxit = 1000,
                     trace = FALSE, MaxNWts = 10000, ...)
{
  x   <- as.matrix(x)
  xp  <- as.matrix(xp)
  
  if (proper == TRUE) {  # bootstrap to make proper
    s   <- sample(length(y), replace = TRUE)
    x   <- x[s,]
    y   <- y[s]
    y   <- factor(y)
  }
  
  aug <- augment.syn(y, x, ...)
  # yf, wf and xf needed for augmented data to save x as non augmented  GR
  xf  <- aug$x
  yf  <- aug$y
  wf  <- aug$w
  #xy  <- cbind.data.frame(y = y,  x = xp)
  xfy <- cbind.data.frame(yf, xf)
  
  ## polr may fail on sparse data. We revert to multinom in such cases. 
  fit <- try(suppressWarnings(polr(formula(xfy), data = xfy, Hess = TRUE, weights = wf, ...)), silent = TRUE)
  
  if (inherits(fit, "try-error")) {
    fit <- multinom(formula(xfy), data = xfy, weights = wf,
                    maxit = maxit, trace = trace, Hess = TRUE, MaxNWts = MaxNWts, ...)
    cat("\tMethod changed to multinomial")
    if (fit$convergence == 1) cat("\nReached max number of iterations for a multinomial model\nRerun with polyreg.maxit increased (default 100)\n")
  }
  post  <- predict(fit, xp, type = "probs")
  
  if (length(y) == 1) post <- matrix(post, nrow = 1, ncol = length(post))
  y     <- as.factor(y)
  nc    <- length(levels(yf))                       
  un    <- rep(runif(nrow(xp)), each = nc)
  if (is.vector(post)) post <- matrix(c(1 - post, post), ncol = 2)
  draws <- un > apply(post, 1, cumsum)
  idx   <- 1 + apply(draws, 2, sum)
  # this slightly clumsy code needed to ensure y retains its labels and levels
  #  y[1:length(y)]<-(levels(y)[idx])
  res <- levels(yf)[idx]
  fitted <- summary(fit)
  return(list(res = res, fit = fitted)) 
}


###-----syn.sample---------------------------------------------------

syn.sample <- function(y, xp, smoothing = "", cont.na = NA, proper = FALSE, ...) 
{
  # Generates random sample from the observed y's
  # with bootstrap if proper == TRUE
  if (proper == TRUE) y <- sample(y, replace = TRUE)
  yp <- sample(y, size = xp, replace = TRUE)
  
  if (smoothing != "") yp[!(yp %in% cont.na)] <- 
    syn.smooth(yp[!(yp %in% cont.na)], y[!(y %in% cont.na)], 
               smoothing = smoothing)
  
  return(list(res = yp, fit = "sample"))
}


###-----syn.passive--------------------------------------------------------
syn.passive <- function(data, func)
{
  # Special elementary synthesis method for transformed data.
  # SuppressWarnings to avoid message 'NAs by coercion for NAtemp  
  res <- suppressWarnings(model.frame(as.formula(func), data, 
                                      na.action = na.pass))	
  
  return(list(res = res, fit = "passive"))
}


###-----syn.cart-----------------------------------------------------------

syn.cart <- function(y, x, xp, smoothing = "", proper = FALSE, 
                     minbucket = 5, cp = 1e-08, ...)
{
  ylogical <- is.logical(y)
  
  if (proper == TRUE) {
    s <- sample(length(y), replace = TRUE)
    x <- x[s,,drop = FALSE]
    y <- y[s]
  }
  
  #for (j in 1:ncol(x)){
  #  if(is.factor(x[,j])) { 
  #    attributes(x[,j])$contrasts <- NULL
  #    attributes(xp[,j])$contrasts <- NULL
  #  }
  #}
  minbucket <- max(1, minbucket)  # safety
  if (!is.factor(y) & !is.logical(y)) {
    fit <- rpart(y ~ ., data = as.data.frame(cbind(y, x)), method = "anova",
                 minbucket = minbucket, cp = cp, ...)
    # get leaf number for observed data
    leafnr  <- floor(as.numeric(row.names(fit$frame[fit$where,])))
    # replace yval with leaf number in order to predict later node number 
    # rather than yval (mean y for observations classified to a leaf) 
    fit$frame$yval <- as.numeric(row.names(fit$frame))
    # predict leaf number
    nodes       <- predict(object = fit, newdata = xp)
    # BN:16/06/20
    # node numbering: node * 2 + 0:1    
    notleaf <- setdiff(nodes, leafnr)
    # if (length(notleaf) > 0) {
    #   for (i in notleaf){
    #     nodes[which(nodes == i)] <- 2 * i + sample(0:1, 1)
    #   }
    # }
    if (length(notleaf) > 0) {
      for (i in notleaf){
        j <- i
        while(!(j %in% leafnr)){
          j <- 2 * j + sample(0:1, 1)
        }
        nodes[which(nodes == i)] <- j
      }
    }
    
    uniquenodes <- unique(nodes)
    new  <- vector("numeric",nrow(xp))
    for (j in uniquenodes) {
      donors <- y[leafnr == j] # values of y in a leaf
      new[nodes == j] <- resample(donors, size = sum(nodes == j), 
                                  replace = TRUE)
    }
    
    if (smoothing != "") new <- syn.smooth(new, y, smoothing = smoothing)
    
    #donor <- lapply(nodes, function(s) y[leafnr == s])
    #new   <- sapply(1:length(donor),function(s) resample(donor[[s]], 1))
  } else {
    y     <- factor(y)
    fit   <- rpart(y ~ ., data = as.data.frame(cbind(y, x)), method = "class",
                   minbucket = minbucket, cp = cp, ...)
    nodes <- predict(object = fit, newdata = xp)
    new   <- apply(nodes, MARGIN = 1, FUN = function(s) resample(colnames(nodes), 
                                                                 size = 1, prob = s))
    if (ylogical) {
      new   <- as.logical(new)
    } else {
      new   <- factor(new, levels = levels(y)) 
    }
  }
  
  return(list(res = new, fit = fit))
}


###-----syn.ctree----------------------------------------------------------

syn.ctree <- function(y, x, xp, smoothing = "", proper = FALSE, minbucket = 5, 
                      mincriterion = 0.9, ...)
{ 
  if (proper == TRUE) {
    s <- sample(length(y), replace = truehist())
    y <- y[s]
    x <- x[s, , drop = FALSE]
  }
  
  for (i in which(sapply(x, class) != sapply(xp,class))) xp[,i] <-
      eval(parse(text = paste0("as.", class(x[,i]), "(xp[,i])", sep = "")))
  
  # Fit a tree
  #  datact     <- partykit::ctree(y ~ ., data = as.data.frame(cbind(y, x)), 
  #                  control = partykit::ctree_control(minbucket = minbucket, 
  #                                                    mincriterion = mincriterion, ...))
  datact <- ctree(y ~ ., data = as.data.frame(cbind(y,x)), 
                  controls = ctree_control(minbucket = minbucket, 
                                           mincriterion = mincriterion, ...))
  
  # fit.nodes  <- predict(datact, type = "node")
  fit.nodes  <- where(datact)
  nodes      <- unique(fit.nodes)
  no.nodes   <- length(nodes)
  # pred.nodes <- predict(datact, type = "node", newdata = xp)
  pred.nodes <- where(datact, newdata = xp)
  # Get row numbers for predicted by sampling with replacement from existing data
  rowno      <- 1:length(y)
  newrowno   <- vector("integer", nrow(xp))
  
  for (i in nodes) {
    newrowno[pred.nodes == i] <- sample(rowno[fit.nodes == i],
                                        length(newrowno[pred.nodes == i]),
                                        replace = TRUE)
  }
  new <- y[newrowno]
  
  if (!is.factor(y) & smoothing != "") new <- 
    syn.smooth(new, y, smoothing = smoothing )
  
  return(list(res = new, fit = datact))
}


###-----syn.survctree------------------------------------------------------

syn.survctree <- function(y, yevent, x, xp, proper = FALSE, minbucket = 5, ...)
  # time, event - data column numbers
{      
  if (proper == TRUE) {
    s <- sample(length(y), replace = TRUE)
    y <- y[s]
    x <- x[s, , drop = FALSE]                        
    yevent <- yevent[s]
  }
  
  for (i in which(sapply(x, class) != sapply(xp,class))) xp[,i] <-
      eval(parse(text = paste0("as.", class(x[,i]), "(xp[,i])", sep = "")))
  
  if (is.factor(yevent)) {
    yevent0 <- as.numeric(yevent) - 1
  } else {
    yevent0 <- yevent 
  }
  
  # Fit a tree  
  datact     <- ctree(Surv(y, yevent0) ~ ., 
                      data = as.data.frame(cbind(y, yevent0, x)),
                      controls = ctree_control(minbucket = minbucket, ...))
  # fit.nodes  <- predict(datact, type = "node")
  fit.nodes  <- where(datact)
  nodes      <- unique(fit.nodes)
  no.nodes   <- length(nodes)
  # pred.nodes <- predict(datact, type = "node", newdata = xp)
  pred.nodes <- where(datact, newdata = xp)
  # Get row numbers for predicted by sampling
  # with replacement from existing data
  rowno      <- 1:length(y)
  newrowno   <- rep(0,nrow(xp))
  for (i in nodes) {
    newrowno[pred.nodes == i] <- sample(rowno[fit.nodes == i],
                                        length(newrowno[pred.nodes == i]), replace = TRUE)
  }
  #Predicte node & sample time+event
  faketime  <- y[newrowno]
  fakeevent <- yevent[newrowno]
  
  return(list(syn.time = faketime, syn.event = fakeevent, fit = datact))
}


###-----syn.rf-------------------------------------------------------------
# bagging when mtry = ncol(x) - using all predictors
syn.rf <- function(y, x, xp, smoothing = "", proper = FALSE, ntree = 10, ...) 
{ 
  
  #nodesize <- max(1, nodesize)  # safety
  #if (proper == TRUE) {
  #  s <- sample(length(y), replace = T); y <- y[s]
  #  x <- x[s, , drop = FALSE]
  #}  
  
  for (i in which(sapply(x,class) != sapply(xp, class))) xp[,i] <-
      do.call(paste0("as.", class(x[,i])[1]), unname(xp[, i]))
  
  if (is.factor(y)) {
    obslevels <- levels(y)
    y <- droplevels(y)
  }
  
  # fit a random forest
  # regression (mtry = p/3), classification (mtry = sqrt(p))
  rf.fit <- randomForest(y ~ ., data = cbind.data.frame(y,x), ntree = ntree, ...)
  nodessyn <- attr(predict(rf.fit, newdata = xp, nodes = T), "nodes")
  nodesobs <- attr(predict(rf.fit, newdata = x, nodes = T), "nodes")
  
  ndonors <- vector("list", nrow(xp))
  n0      <- vector("list", ntree)
  for (j in 1:nrow(xp)) {
    for (i in 1:ntree) {
      n0[[i]] <- y[nodesobs[,i] == nodessyn[j,i]]
    }
    empty <- sapply(n0, length)
    ndonors[[j]] <- unlist(n0[empty != 0])
  }
  
  yhat <- sapply(ndonors, sample, size = 1)          
  
  if (is.factor(y)) yhat <- factor(yhat, levels = obslevels) 
  if (!is.factor(y) & smoothing != "") yhat <- 
    syn.smooth(yhat, y, smoothing = smoothing)
  
  return(list(res = yhat, fit = rf.fit)) 
}


###-----syn.ranger---------------------------------------------------------
# bagging when mtry = ncol(x) - using all predictors
# contributed by Caspar J. van Lissa

syn.ranger <- function(y, x, xp, smoothing = "", proper = FALSE, ...) 
{ 
  dots <- list(...)
  dots[c("formula", "data")] <- NULL
  if("min.node.size" %in% names(dots)){
    dots[["min.node.size"]] <- max(1, dots[["min.node.size"]])  # safety
  }
  
  if (proper == TRUE) {
    s <- sample(length(y), replace = TRUE)
    y <- y[s]
    x <- x[s, , drop = FALSE]
  }
  
  for (i in which(sapply(x,class) != sapply(xp, class))) xp[,i] <-
    do.call(paste0("as.", class(x[,i])[1]), unname(xp[, i]))
  
  if (is.factor(y)) {
    obslevels <- levels(y)
    y <- droplevels(y)
  }
  
  # fit a random forest
  Args     <- c(list(formula = y ~ ., data = cbind.data.frame(y,x)), dots)
  rf.fit   <- do.call(ranger, Args)
  nodessyn <- predict(rf.fit, data = xp, type = "terminalNodes")$predictions
  nodesobs <- predict(rf.fit, data = x, type = "terminalNodes")$predictions
  ntree    <- rf.fit$num.trees
  ndonors  <- vector("list", nrow(xp))
  n0       <- vector("list", ntree)
  for (j in 1:nrow(xp)) {
    for (i in 1:ntree) {
      n0[[i]] <- y[nodesobs[,i] == nodessyn[j,i]]
    }
    empty <- sapply(n0, length)
    ndonors[[j]] <- unlist(n0[empty != 0])
  }
  
  yhat <- sapply(ndonors, sample, size = 1)          
  
  if (is.factor(y)) yhat <- factor(yhat, levels = obslevels) 
  if (!is.factor(y) & smoothing != "") yhat <- 
    syn.smooth(yhat, y, smoothing = "smoothing")
  
  return(list(res = yhat, fit = rf.fit))
}


###-----syn.bag-------------------------------------------------------------
# bagging when mtry = ncol(x) - using all predictors
syn.bag <- function(y, x, xp, smoothing = "", proper = FALSE, ntree = 10, ...) 
{ 
  #nodesize <- max(1, nodesize)  # safety
  #if (proper == TRUE) {
  #  s <- sample(length(y), replace = T); y <- y[s]
  #  x <- x[s, , drop = FALSE]
  #}  
  
  for (i in which(sapply(x,class) != sapply(xp, class))) xp[,i] <-
      do.call(paste0("as.", class(x[,i])[1]), unname(xp[, i]))
  
  if (is.factor(y)) {
    obslevels <- levels(y)
    y <- droplevels(y)
  }
  
  # fit a random forest
  # regression (mtry = p/3), classification (mtry = sqrt(p))
  rf.fit <- randomForest(y ~ ., data = cbind.data.frame(y,x), 
                         ntree = ntree, mtry = ncol(x), ...)
  nodessyn <- attr(predict(rf.fit, newdata = xp, nodes = T), "nodes")
  nodesobs <- attr(predict(rf.fit, newdata = x, nodes = T), "nodes")
  
  ndonors <- vector("list", nrow(xp))
  n0      <- vector("list", ntree)
  for (j in 1:nrow(xp)) {
    for (i in 1:ntree) {
      n0[[i]] <- y[nodesobs[,i] == nodessyn[j,i]]
    }
    empty <- sapply(n0, length)
    ndonors[[j]] <- unlist(n0[empty != 0])
  }
  
  yhat <- sapply(ndonors, sample, size = 1)
  
  if (is.factor(y)) yhat <- factor(yhat, levels = obslevels) 
  if (!is.factor(y) & smoothing != "") yhat <- 
    syn.smooth(yhat,y, smoothing = smoothing)
  
  return(list(res = yhat, fit = rf.fit))
}


###-----syn.nested---------------------------------------------------------
# function for allocating to subcategories (random sampling within groups)

syn.nested <- function(y, x, xp, smoothing = "", cont.na = NA, ...)
{
  xr   <- x[,1]
  xpr  <- xp[,1]
  uxpr <- sort(unique(xpr))
  
  index  <- 1:length(y)
  indexp <- rep(0, nrow(xp))
  for (i in uxpr) {
    indexp[xpr == i] <- sample(index[xr == i], sum(xpr == i), TRUE)
  }
  yp <- y[indexp]
  
  if (smoothing != "") yp[!(yp %in% cont.na)] <-
    syn.smooth(yp[!(yp %in% cont.na)], y[!(y %in% cont.na)], 
               smoothing = smoothing)
  
  return(list(res = yp, fit = "nested"))
}


###-----syn.satcat---------------------------------------------------------

syn.satcat <- function(y, x, xp, proper = FALSE, ...)
{
  # Fits a saturated model to combinations of variables.
  # Method fails if the predictor variables generate
  # a combination of variables not found in the original data.
  if (proper == TRUE) {
    s <- sample(length(y), replace = TRUE)
    y <- y[s]
    x <- x[s, , drop = FALSE]
  }
  
  xr  <- apply(x, 1, function(x) paste(x, collapse = "-"))
  syn.categories <- apply(xp, 1, function(x) paste(x, collapse = "-"))
  
  if (!all(names(table(syn.categories)) %in% names(table(xr)))) {
    cat("\n\n")
    print(table(syn.categories)[!names(table(syn.categories)) %in% names(table(xr))])
    stop('The combined groups above for "satcat" have no records in original data.\n
         Consider using grouped synthesis with "syn.catall" to overcome this', call. = FALSE)
  }
  
  uxpr <- sort(unique(syn.categories))
  
  index  <- 1:length(y)
  indexp <- rep(0, nrow(xp))
  for (i in uxpr) {
    indexp[syn.categories == i] <- sample(index[xr == i], sum(syn.categories == i), TRUE)
  }
  yp <- y[indexp]
  fit <- table(xr)
  return(list(res = yp, fit = fit))
}


###-----syn.constant-------------------------------------------------------

syn.constant <- function(y, xp, ...) 
{
  yp <- y
  length(yp) <- xp
  if (xp > length(y)) {
    yp[(length(y) + 1):xp] <- names(which.max(table(y, exclude = NULL)))  # in case y is 'almost' constant 
    if (is.numeric(y)) yp <- as.numeric(yp)
    else if (is.logical(y)) yp <- as.logical(yp)
  }
  return(list(res = yp, fit = "constant"))
}


###-----syn.collinear------------------------------------------------------

syn.collinear <- function(y, x, xp, ...)
{
  x <- x[,1]                               #!BN to check
  xp <- xp[,1]                             #!BN to check
  indexp  <- match(xp, x)
  yp <- y[indexp]
  return(list(res = yp, fit = "collinear"))
}


###-----syn.catall---------------------------------------------------------

syn.catall <- function(x, k, proper = FALSE, priorn = 1, structzero = NULL, 
                       maxtable = 1e8, epsilon = 0, rand = TRUE, ...)
{
  # Fits a saturated model to combinations of variables
  # xp just holds number of synthetic records required
  
  levs <- sapply(x, function(x) {length(levels(x)) + any(is.na(x))})  # all NAtemp here already
  table.size <- prod(levs)   # exp(sum(log(levs)))
  if (table.size > maxtable) stop("Table has more than ", maxtable/1e6,
                                  " million cells (", round(table.size/1e6, 2),
                                  " millions),\nwhich may lead to memory problems.\nYou can rerun syn() with catall.maxtable increased (default: ", maxtable/1e6,
                                  " millions).\nAlternatively use a smaller group of variables for catall.", 
                                  sep = "", call. = FALSE)
  
  N <- dim(x)[1]
  if (proper == TRUE) x <- x[sample(1:N, replace = TRUE), ]
  tab <- table(x)
  n <- length(tab) 
  addon <- priorn/n
  # Set structural zero cells to zeros 
  if (!is.null(structzero)) {
    sz <- checksz(structzero, x) ## checks and converts var names to numbers
    if (sum(tab[sz]) > 0) cat("
\n************************************************************************
WARNING: Total of ", sum(tab[sz])," counts of original data in structural zero cells.
************************************************************************\n", sep = "")
  }
  # Add extra to prior
  tab <- (tab + addon)
  if (!is.null(structzero)) tab[sz] <- 0
  dt  <- dim(tab)
  dn  <- dimnames(tab)
  if (epsilon > 0) {
    if (rand == TRUE) {
      if (!is.null(structzero)) tab[!sz] <- addlapn(tab[!sz], epsilon) 
      else tab <- addlapn(tab, epsilon) 
      fit <- tab
      tab <- tab/sum(tab)   # get it as proportions
      tab <- rmultinom(1, k, tab)
    } else {
      if (!is.null(structzero)) tab[!sz] <- addlapn(tab[!sz], epsilon) 
      else tab <- addlapn(tab, epsilon ) #GR2022
      fit <- tab
      tab <- roundspec(tab*k/sum(tab))
    }
  } else {
    if (rand == TRUE) {
      fit <- tab
      tab <- tab/sum(tab)   # get it as proportions
      tab <- rmultinom(1, k, tab)
    }
    else stop("If you set rand = FALSE when epsilon = 0 synthpop will return the data unchanged", call. = FALSE) 
  }
  tab <- array(tab, dt)
  dimnames(tab) <- dn
  res <- array.to.frame(tab)
  for (i in 1:dim(res)[2]) res[, i] <- addNA(res[, i], ifany = TRUE)
  res <- res[sample(1:nrow(res)), ]
  return(list(res = res, fit = fit))
}


###-----syn.ipf------------------------------------------------------------

syn.ipf <- function(x, k, proper = FALSE, priorn = 1, structzero = NULL, 
                    gmargins = "twoway", othmargins = NULL, tol = 1e-3, max.its = 5000,
                    maxtable = 1e8, print.its = FALSE, epsilon= 0, rand = TRUE,...)
{
  # Fits log-linear model to combinations of variables
  # k just holds number of synthetic records required
  
  levs <- sapply(x, function(x) {length(levels(x)) + any(is.na(x))}) # all NAtemp here already
  table.size <- prod(levs)   # exp(sum(log(levs))) 
  if (table.size > maxtable) stop("Table has more than ", maxtable/1e6,
                                  " million cells (", round(table.size/1e6, 2),
                                  " millions),\nwhich may lead to memory problems.\nYou can rerun syn() with ipf.maxtable increased (default: ", maxtable/1e6,
                                  " millions).\nAlternatively use a smaller group of variables for ipf.", 
                                  sep = "", call. = FALSE)
  
  N  <- dim(x)[1]
  nv <- dim(x)[2] # number of variables
  if (proper == TRUE) x <- x[sample(1:N, replace = TRUE),]
  tab <- table(x, useNA = "ifany")
  n <- length(tab)
  # Add extra to prior
  addon <- priorn/n
  # Set structural zero cells to zeros 
  if (!is.null(structzero)) {
    sz <- checksz(structzero, x) # checks and converts var names to numbers
    if (sum(tab[sz]) > 0) cat("
\n************************************************************************
WARNING: Total of ", sum(tab[sz])," counts of original data in structural zero cells.
************************************************************************\n", sep = "")
  }
  if (!is.null(gmargins)) {
    if (gmargins == "twoway") {
      n_margins  <- nv*(nv - 1)/2
      mx_margins <- combn(1:nv, 2)
      margins <- split(mx_margins, col(mx_margins))
    } else if (gmargins == "oneway") {
      margins <- as.list(1:nv)
    } else stop("Only 'oneway' or 'twoway' are implemented for gmargins.\n", 
                call. = FALSE)
    if (!is.null(othmargins)) for (i in 1:length(othmargins)) {
      margins[[length(margins) + 1]] <- othmargins[[i]]
    }
  } else {
    if (!is.null(othmargins)) margins <- othmargins
    else stop("Need to specify some margins.\n", call. = FALSE)
  }
  
  umar <- unique(unlist(margins))
  missed <- (1:nv)[!(1:nv %in% umar)]
  
  if (length(missed) > 0) cat("\n
**************************************************
SEVERE WARNING: Margins ", missed, " not fitted.
This means they will be fitted as having
the same proportion in each level.
**************************************************\n", sep = "")
  if (epsilon > 0) {eps <- epsilon / length(margins)
  cat("Overall epsilon for DP is ",epsilon," divided equally between ", 
      length(margins)," margins to give ", eps, " each,\n")
  } 
  
  # Get data for margins
  margins.data <- vector("list", length(margins))
  for (i in 1:length(margins)) {
    margins.data[[i]] <- table(x[, margins[[i]]], useNA = "ifany")
    margins.data[[i]] <- margins.data[[i]] + priorn/length(margins.data[[i]])
    if (epsilon > 0) {
      margins.data[[i]] <- addlapn(margins.data[[i]], eps)
    }
  }
  start <- array(1, dim(tab))
  
  results.1 <- suppressWarnings(Ipfp(start, margins, margins.data, 
                                     iter = max.its, print = print.its, tol = tol, ...))  # note larger than default to speed things up needs checking
  
  if (results.1$conv == TRUE) cat("\n['ipf' converged in ", 
                                  length(results.1$evol.stp.crit), " iterations]\n", sep = "")
  else cat("\n['ipf' failed to converge in ", length(results.1$evol.stp.crit),
           " iterations]\n\nYou can try to change parameters of Ipfp() function,\ne.g. syn(..., ipf.iter = 2000).", sep = "")
  
  exp1 <- array(N * results.1$p.hat , dim(tab))
  if (!is.null(structzero)) exp1[sz] <- 0
  
  exp1 <- exp1 / sum(exp1)   # get it as proportions
  if (epsilon == 0 & rand == FALSE) cat("WARNING: No DP noise or random noise added,\nData returned will be close to NULL expectation for model defined by margins.\n")
  if (rand == TRUE)  z    <- rmultinom(1, k, exp1)
  else z <- roundspec(exp1*k/sum(exp1))
  res  <- array(z, dim(exp1))
  dimnames(res) <- dimnames(tab)
  res  <- array.to.frame(res)
  for (i in 1:dim(res)[2]) res[,i] <- addNA(res[,i], ifany = TRUE)
  res <- res[sample(1:nrow(res)),]
  fit <- list(margins = margins, margins.data = margins.data)
  return(list(res = res, fit = fit))
}


# O T H E R   A U X I L I A R Y   F U N C T I O N S  
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 

###-----is.passive---------------------------------------------------------
is.passive <- function(string) return("~" == substring(string,1,1))


###-----resample-----------------------------------------------------------
# used in syn.cart() and syn.cartbboot() instead of sample() 
# for safety reasons, i.e. to avoid sampling from 1:x when length(x)==1
resample <- function(x, ...) x[sample.int(length(x), ...)]


###-----decimalplaces------------------------------------------------------
# counts number of decimal places - used for rounding smoothed values
# approximate in some cases (as.character -> 15 significant digits; 
# scientific notation)
decimalplaces <- function(x) 
{
  x <- x - floor(x) # -> more digit numbers 
  if (!is.na(x) & (x %% 1) != 0 & (round(x, 15) %% 1 != 0)) {
    nchar(strsplit(sub("0+$", "", as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}


###-----get.names----------------------------------------------------------
get.names <- function(formula, names)
{
  res <- match(all.vars(formula)[-1], names)
  return(res)
}


###-----addXfac------------------------------------------------------------
# function to add factors
addXfac <- function(x,...)
{
  df  <- cbind.data.frame(x,...)
  if (any(sapply(df,is.factor) + sapply(df,is.numeric) != 1)) 
    stop("All arguments must be factors or numeric", call. = FALSE)
  fn  <- function(f){
    if (is.factor(f)) f <- as.numeric(levels(f))[f]
    else f <- f
  }
  df  <- as.data.frame(sapply(df, fn))
  add <- factor(rowSums(df))
  return(add)
}


###-----syn.smooth--------------------------------------------------------- 

syn.smooth <- function(ysyn, yobs = NULL, smoothing = "spline", window = 5, ...)
{
  if (!(smoothing %in% c("", "spline", "density", "rmean"))) 
    cat('Smoothing must be one of "spline", "density" or "rmean". No smoothing done.\n')
  if (any(is.na(ysyn))) stop("ysyn cannot contain missing values", call. = FALSE)
  
  else if (smoothing == "density") {
    ys <- 1:length(ysyn)
    # exclude from smoothing if freq for a single value higher than 70% 
    maxfreq <- which.max(table(ysyn))
    maxcat  <- as.numeric(names(table(ysyn))[maxfreq])
    if (table(ysyn)[maxfreq]/sum(table(ysyn)) > .7) ys <- which(ysyn != maxcat)
    # exclude from smoothing if data are top-coded - approximate check
    if (10 * table(ysyn)[length(table(ysyn)) - 1] <
        tail(table(ysyn), n = 1) - table(ysyn)[length(table(ysyn)) - 1]) {
      ys   <- ys[-which(ysyn == max(yobs))]
      maxy <- max(yobs)
    }
    
    densbw  <- density(ysyn[ys], width = "SJ")$bw
    ysyn[ys] <- rnorm(n = length(ysyn[ys]), mean = ysyn[ys], sd = densbw)
    if (!exists("maxy")) maxy <- max(yobs) + densbw
    ysyn[ys] <- pmax(pmin(ysyn[ys], maxy), min(yobs))
    ysyn[ys] <- round(ysyn[ys], max(sapply(yobs, decimalplaces))) 
  }    
  else if (smoothing == "rmean") {
    ord <- order(ysyn)
    ysyn <- runningmean(1:length(ysyn), ysyn[ord], window = window)
    ysyn[ord] <- round(ysyn, max(sapply(yobs, decimalplaces))) 
  }
  else if (smoothing == "spline") {
    ord <- order(ysyn)
    ysyn <- smooth.spline(sort(ysyn), all.knots = FALSE)$y
    ysyn[ord] <- round(ysyn, max(sapply(yobs, decimalplaces))) 
  }
  
  return(ysyn)
}


###-----checksz------------------------------------------------------------

checksz <- function(sz, x) 
{
  if (!is.list(sz)) stop("structzero needs to be a list.\n", call. = FALSE)
  if (!is.character(names(sz)) || !all(grepl("_", names(sz)))) stop("\nstructzero list elements must be named using variable names\nseperated by an underscore, e.g. sex_edu", call. = FALSE)
  
  list_sub_names <- names(do.call("c", sz))
  allvars <- unique(sub(".*\\.", "", list_sub_names))  
  
  if (!all(allvars %in% names(x))) stop("\nStructural zero variables must match names of variables in the data.", call. = FALSE)
  
  dd  <- as.data.frame(table(x, useNA = "ifany"))
  res <- rep(FALSE, nrow(dd))
  
  for (i in 1:length(sz)) {
    
    tempz <- rep(TRUE, nrow(dd))
    if (names(sz)[i] != paste0(names(sz[[i]]), collapse = "_")) stop("\nNames of structzero list elements must correspond to names of their sublists.", call. = FALSE)
    vars  <- names(sz[[i]])    
    vars  <- match(vars, names(x))
    nvars <- length(vars)
    if (!(is.list(sz[[i]]) & length(sz[[i]]) == nvars)) stop("Each element of structzero list must be a list of length\nequal to the number of variables used to define structural zeros.\n", call. = FALSE)
    
    for (j in 1:nvars) {
      if (!is.numeric(sz[[i]][[j]])) {
        if (!all(sz[[i]][[j]] %in% levels(x[,vars[j]]))) 
          stop("Structural zeros (element ", i, ", variable ", j, 
               "): level(s) of variable not in data.\n", call. = FALSE)
        sz[[i]][[j]] <- match(sz[[i]][[j]], levels(x[, vars[j]]))
      }
      if (!all(sz[[i]][[j]] %in% 1:nlevels(x[,vars[j]]))) 
        stop("Structural zeros (element ", i, ", variable ", j, 
             "): numeric level(s) of variable not in data.\n", sep = "", 
             call. = FALSE)
      tempz <- tempz & (as.numeric(dd[, vars[j]]) %in% sz[[i]][[j]])
    }
    res[tempz] <- TRUE
  }
  return(res)
}


###-----array.to.frame-----------------------------------------------------

array.to.frame <- function(x)
{
  df1 <- as.data.frame.table(x)
  # df1 <- df1[df1$Freq != 0,]
  res <- df1[rep(1:nrow(df1), df1$Freq), -ncol(df1)]
  dimnames(res)[[1]] <- 1:sum(x)
  return(res)
}


###-----roundspec----------------------------------------------------------

roundspec <- function(tab) {  ## special rounding to preserve total
  if (abs(round(sum(tab)) - sum(tab)) > 1e-6) 
    stop("This function assumes sum of tab is a whole number\n", .call= FALSE)
  diff <- round(sum(tab) - sum(round(tab)))
  
  if (diff == 0 & all(tab >= 0)) { 
    result <- round(tab)
  } else if (any(round(tab) <= 0)){
    if (diff < 0 ) {
      inds <- (1:length(tab))[round(tab) >0]
      vals <- tab[round(tab) > 0]
      ordinds <- inds[order(vals)]
      newtab <- round(tab)
      newtab[ordinds[1:(-diff)]]<- newtab[ordinds[1:(-diff)]]  - 1
    } else {
      inds <- (1:length(tab))
      newtab <- round(tab)
      newtab[1:(diff)]<- newtab[1:(diff)] + 1
      
    }
    result <- newtab
  } else {
    if (abs(diff) > length(tab)) cat("Unexpected problem with rounding.\n
Please report to maintainer of synthpop with example\n
including this output\n", tab, "\n")
    result <- round(tab)
    inds <- (1:length(result))[order(result)][1:abs(diff)]
    if (diff >0 ) result[inds] <- result[inds] + 1
    if (diff <0 ) result[inds] <- result[inds] - 1
  }
  return(result)
}


###-----addlapn------------------------------------------------------------

addlapn <- function(x, eps){
  # add Laplace noise with re-scaling to a total or sample size
  
  if (eps <= 0) stop("eps must be > 0\n")
  res <- x + rlaplace(length(x), 0, 1/eps)
  res <- makepos(res, sum(x))
  return(res)
}


###---------------makepos--------------------------------------------------

makepos <- function(lap, tot) {  ## makes positive summing to total
  olap <- order(-lap)
  lapo <- lap[olap]
  lap[olap[cumsum(abs(lapo)) >= tot]] <- 0
  abs(lap)
}






is.passive <- function(string) return("~" == substring(string,1,1))

###-----replicated.uniques-------------------------------------------------
# unique units in the synthesised data that replicates unique real units
# (+number +percent of all observations)

replicated.uniques <- function(object, data, exclude = NULL){
  # check names of vars to be excluded
  if (!is.null(exclude)) {
    exclude.cols <- match(exclude, colnames(data))
    if (any(is.na(exclude.cols))) stop("Unrecognized variable(s) in exclude parameter for uniques: ",
                                       paste(exclude[is.na(exclude.cols)],collapse=", "), call. = FALSE)
  }
  
  # exclude vars from the original data
  if (!is.null(exclude)) data <- data[,-exclude.cols, drop = FALSE] else data <- data
  # extract uniques from the original data
  uReal <- data[!(duplicated(data) | duplicated(data,fromLast=TRUE)), , drop = FALSE]
  no.uniques <- nrow(uReal)
  
  if (is.null(no.uniques) || no.uniques == 0) {
    no.uniques <- 0
    no.duplicates <- per.duplicates <- rep(0,object$m)
    if (object$m == 1) rm.Syn <- rep(FALSE,nrow(object$syn))
    if (object$m > 1) rm.Syn <- matrix(FALSE,nrow=nrow(object$syn[[1]]),ncol=object$m)
    
  } else {
    
    if (object$m == 1){
      if (!is.null(exclude)) Syn <- object$syn[,-exclude.cols, drop = FALSE] else Syn <- object$syn
      rm.Syn <- rep(FALSE,nrow(Syn))
      i.unique.Syn <- which(!(duplicated(Syn) | duplicated(Syn,fromLast=TRUE)))
      if (length(i.unique.Syn)!=0) {
        uSyn <- Syn[i.unique.Syn, , drop = FALSE]
        uAll <- rbind.data.frame(uReal,uSyn)
        dup.of.unique <- duplicated(uAll)[(nrow(uReal)+1):nrow(uAll)]
        rm.Syn[i.unique.Syn] <- dup.of.unique
      }
      no.duplicates <- sum(rm.Syn)
    }
    
    if (object$m > 1){
      rm.Syn <- matrix(FALSE,nrow=nrow(object$syn[[1]]),ncol=object$m)
      for (i in 1:object$m){
        if (!is.null(exclude)) Syn <- object$syn[[i]][,-exclude.cols,drop = FALSE] else Syn <- object$syn[[i]]
        i.unique.Syn <- which(!(duplicated(Syn) | duplicated(Syn,fromLast=TRUE)))
        if (length(i.unique.Syn)!=0) {
          uSyn <- Syn[i.unique.Syn, , drop = FALSE]
          uAll <- rbind.data.frame(uReal,uSyn)
          dup.of.unique <- duplicated(uAll)[(nrow(uReal)+1):nrow(uAll)]
          rm.Syn[i.unique.Syn,i] <- dup.of.unique
        }
      }
      no.duplicates <- colSums(rm.Syn)
    }
    per.duplicates <- no.duplicates/nrow(data)*100
  }
  
  return(list(replications = rm.Syn, no.uniques = no.uniques,
              no.replications = no.duplicates, per.replications = per.duplicates))
}


###-----sdc----------------------------------------------------------------
# sdc - statistical disclosure control:
# labeling, removing unique replicates of unique real individuals

sdc <- function(object, data, label = NULL, rm.replicated.uniques = FALSE, 
                uniques.exclude = NULL, recode.vars = NULL, bottom.top.coding = NULL,
                recode.exclude = NULL, smooth.vars = NULL){
  
  if (!is.null(smooth.vars)) {
    if (object$m == 1) { 
      if (any(!smooth.vars %in% names(object$syn))) stop("Some of smooth.vars not in the data", call. = FALSE)  
      if (any(!(sapply(object$syn[, smooth.vars], function(x) is.numeric(x) | is.integer(x))))) stop("Some of smooth.vars not numeric", call. = FALSE)  
    } else {
      if (any(!smooth.vars %in% names(object$syn[[1]]))) stop("Some of smooth.vars not in the data", call. = FALSE)  
      if (any(!(sapply(object$syn[[1]][, smooth.vars], function(x) is.numeric(x) | is.integer(x))))) stop("Some of smooth.vars not numeric", call. = FALSE)  
    }  
  }
  
  if (!is.null(recode.vars)) {
    if (!is.null(bottom.top.coding) && !is.list(bottom.top.coding)) 
      bottom.top.coding <- list(bottom.top.coding)
    if (!is.null(recode.exclude) && !is.list(recode.exclude)) 
      recode.exclude <- list(recode.exclude)
    if (length(bottom.top.coding) != length(recode.vars) | 
        any(sapply(bottom.top.coding,length) != 2)) 
      stop("Bottom and top codes have to be provided for each variable in recode.vars.\nUse NA if there is no need for bottom or top recoding.\nFor more than one variable to be recoded provide a list of two-element vectors, e.g. list(c(0,60),c(NA,5000))",
           call. = FALSE)
    if (!is.null(recode.exclude) && length(bottom.top.coding) != length(recode.exclude))
      stop("recode.exclude have to include codes for each variable in recode.vars.\nUse NA if all values should be considered for recoding.\nFor more than one variable to be recoded provide a list, e.g. list(NA,c(NA,-8)).",
           call. = FALSE)
  }
  
  if (object$m == 1) {
    if (!is.null(recode.vars)) {
      cols <- match(recode.vars,colnames(object$syn)) 
      for (i in cols) {
        j <- match(i,cols) 
        recoded <- bottom.top.recoding(object$syn[,i],bottom.top.coding[[j]][1],
                                       bottom.top.coding[[j]][2],recode.exclude[[j]])
        object$syn[,i] <- recoded$x
        cat("\n",recode.vars[j],": no. of bottom-coded values - ",
            recoded$no.recoded.bottom,", no. of top-coded values - ",
            recoded$no.recoded.top, sep = "")
      }
      cat("\n")
    }
    if (rm.replicated.uniques) {
      du <- replicated.uniques(object, data, exclude = uniques.exclude) 
      object$syn <- object$syn[!du$replications,]
      cat("no. of replicated uniques: ", du$no.replications, "\n", sep = "")
    }
    if (!is.null(label)) object$syn <- cbind.data.frame(flag = label, object$syn)
    
    if (!is.null(smooth.vars)) {
      numindx  <- which(names(object$syn) %in% smooth.vars)
      for (i in numindx) {
        yy <- object$syn[,i][!(object$syn[,i] %in% object$cont.na[[i]])]  
        yyrank <- rank(yy)
        yyforsmooth <- sort(yy)
        yysmoothed  <- smooth.spline(yyforsmooth, all.knots = FALSE)
        object$syn[,i][!(object$syn[,i] %in% object$cont.na[[i]])] <- yysmoothed$y[yyrank]  
      }     
    } 
  }
  
  if (object$m > 1) {
    if (!is.null(recode.vars)) {
      cols <- match(recode.vars,colnames(object$syn[[1]])) 
      for (k in 1:object$m) {
        cat("\nm =",k)
        for (i in cols) {
          j <- match(i,cols) 
          recoded <- bottom.top.recoding(object$syn[[k]][,i],bottom.top.coding[[j]][1],
                                         bottom.top.coding[[j]][2],recode.exclude[[j]])
          object$syn[[k]][,i] <- recoded$x
          cat("\n",recode.vars[j], ": no. of bottom-coded values - ",
              recoded$no.recoded.bottom, ", no. of top-coded values - ",
              recoded$no.recoded.top, sep = "")
        }
      }
      cat("\n")
    }
    if (rm.replicated.uniques) {
      du <- replicated.uniques(object, data, exclude = uniques.exclude) 
      for (i in 1:object$m) { 
        object$syn[[i]] <- object$syn[[i]][!du$replications[,i],]
      }
      cat("no. of replicated uniques: ", 
          paste0(du$no.replications, collapse = ", "),"\n", sep="")
    }
    if (!is.null(label)) object$syn <- mapply(cbind.data.frame, flag=label,
                                              object$syn, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    
    if (!is.null(smooth.vars)){
      numindx  <- which(names(object$syn[[1]]) %in% smooth.vars)
      for (k in 1:object$m){
        for (i in numindx){
          yy <- object$syn[[k]][,i][!(object$syn[[k]][,i] %in% object$cont.na[[i]])]  
          yyrank <- rank(yy)
          yyforsmooth <- sort(yy)
          yysmoothed <- smooth.spline(yyforsmooth, all.knots = FALSE)
          object$syn[[k]][,i][!(object$syn[[k]][,i] %in% object$cont.na[[i]])] <- yysmoothed$y[yyrank]  
        }
      }
    }
  }
  return(object) 
}


###---- bottom.top.recoding -----------------------------------------------

bottom.top.recoding <- function(x,bottom,top,exclude=NULL){
  below <- which(x < bottom & !x%in%exclude); no.below <- length(below)
  above <- which(x > top & !x%in%exclude); no.above <- length(above)
  x[below] <- bottom
  x[above] <- top
  return(list(x=x,no.recoded.bottom=no.below,no.recoded.top=no.above))
}


###---- read.obs ----------------------------------------------------------

read.obs <- function(file, convert.factors = TRUE, lab.factors = FALSE, 
                     export.lab = FALSE, ...){
  
  pos <- regexpr("\\.([[:alnum:]]+)$", file)
  ext <- ifelse(pos > -1L, substring(file, pos + 1L), "")
  
  if (ext=="sav") {
    real.data <- read.spss(file, to.data.frame = FALSE, 
                           use.value.labels = convert.factors, 
                           trim.factor.names = TRUE, ...)
    # trim.factor.names=T - trim trailing spaces from factor levels
    # use.value.labels=F -> to prevent combining factor levels with missing labels
    # for read.spss -> cbind(value=attributes(data)$label.table$...)
    # {Hmisc} real.data <- spss.get(file,use.value.labels=TRUE,max.value.labels=20)
    varlab <- attr(real.data, "variable.labels")
    vallab <- attr(real.data, "label.table")
    vallab <- lapply(vallab,rev)
    codepage <- attr(real.data,"codepage")
    real.data <- as.data.frame(real.data)
    attr(real.data, "variable.labels") <- varlab
    attr(real.data, "label.table") <- vallab
    if (codepage > 500) attr(real.data, "codepage") <- codepage
    
    labs <- list(var.lab=varlab,val.lab=vallab)
    if (export.lab) dput(labs,"SPSSlabels")
    
    if (convert.factors == FALSE & lab.factors == TRUE){  
      # convert completly labeled variables into factors
      ff <- !sapply(vallab,is.null)
      suppressWarnings(llall <-
                         sapply(vallab, function(x) !(sum(!is.na(as.numeric(names(x))))!=0)))
      factors <- which(ff & llall)
      for (i in factors){
        real.data[,i] <- factor(real.data[,i],levels=vallab[[i]])
      }
    }
    
  } else if (ext=="dta") {
    real.data <- read.dta(file, convert.factors = convert.factors, ...)
    
    varlab <- attr(real.data, "val.labels")
    vallab <- attr(real.data, "label.table")
    labs <- list(var=varlab,val=vallab)
    if (export.lab) dput(labs,"Statalabels")
    
    if (convert.factors == FALSE & lab.factors == TRUE){  
      # convert completly labeled variables into factors
      ff <- which(varlab != "")
      suppressWarnings(varlaball <-
                         sapply(vallab, function(x) !(sum(!is.na(as.numeric(names(x))))!=0)))
      factors <- ff[varlaball[varlab[ff]]]
      for (i in factors){
        real.data[,i] <- factor(real.data[,i],levels=vallab[[varlab[i]]])
      }
    }
    
    # for read.dta -> cbind(value=attributes(data)$label.table$...)
    # convert.factors=T, convert.underscore=F, warn.missing.labels=T, missing.type=T
    # {Hmisc} real.data <- stata.get(file,...)
    
  } else if (ext=="xpt") {
    real.data <- read.xport(file)
    # {Hmisc} real.data <- sasxport.get(file, ...)
    
  } else if (ext=="csv") {
    real.data <- read.csv(file, header=TRUE, ...)
    
  } else if (ext=="txt") {
    real.data <- read.table(file, header=TRUE, ...)
    
  } else if (ext=="tab") {
    real.data <- read.table(file, header=TRUE, sep="\t")
    
  } else {
    stop(".",ext," is an unrecognized data format",call.=FALSE)
  }
  
  attr(real.data,"filetype") <- ext
  return(real.data)
}
# R files (*.RData, *.rda)
# load(".rda") - don't assign to an object!


###---- write.syn ---------------------------------------------------------

write.syn <- function(object, filename,
                      filetype = c("SPSS","Stata","SAS","csv","tab","rda","RData","txt"),
                      convert.factors = "numeric", data.labels = NULL,
                      save.complete = TRUE, extended.info = TRUE, ...){
  
  # pos <- regexpr("\\.([[:alnum:]]+)$", file)
  # ext <- ifelse(pos > -1L, substring(file, pos + 1L), "")
  # without.ext <- strsplit(file,"\\.")[[1]][1]      #! will include path
  
  m <- object$m
  if (m == 1) object$syn <- list(object$syn)
  filetype <- match.arg(filetype)
  if (is.null(data.labels)) data.labels <- list(var.lab = object$var.lab,
                                                val.lab = object$val.lab)
  
  if (filetype=="SPSS"){
    if (m==1) {
      f1 <- paste0(filename,".sps"); f2 <- paste0(filename,".txt")
    } else {
      f1 <- paste0(filename,"_",1:m,".sps"); f2 <- paste0(filename,"_",1:m,".txt")
    }
    for (i in 1:m) {
      #write.foreign(object$syn[[i]], codefile=f1, datafile=f2, package=filetype, ...)
      write.syn.SPSS(object$syn[[i]], codefile = f1[i], datafile = f2[i],
                     varnames = names(object$syn[[i]]), data.labels = data.labels, ...)
    }
  } else if (filetype == "Stata"){
    if (m==1) f1 <- paste0(filename,".dta") else f1 <- paste0(filename,"_",1:m,".dta")
    for (i in 1:m){
      write.dta(object$syn[[i]], file = f1[i], convert.factors = convert.factors, ...)
      #!### check why default convert.factors="labels" cuts the names
    }
  } else if (filetype == "SAS"){
    if (m==1) {
      f1 <- paste0(filename,".sas"); f2 <- paste0(filename,".txt")
    } else {
      f1 <- paste0(filename,"_",1:m,".sas"); f2 <- paste0(filename,"_",1:m,".txt")
    }
    for (i in 1:m){
      write.foreign(object$syn[[i]], codefile = f1[i], datafile = f2[i], package = filetype, ...)
    }
  } else if (filetype == "rda" | filetype == "RData"){
    if (m==1) f1 <- paste(filename,filetype,sep=".") else f1 <- paste0(filename,"_",1:m,".",filetype)
    for (i in 1:m){
      syn <- object$syn[[i]]
      save(syn, file = f1[i],...)
    }
  } else if (filetype == "csv"){
    if (m==1) f1 <- paste0(filename,".csv") else f1 <- paste0(filename,"_",1:m,".csv")
    for (i in 1:m){
      write.csv(object$syn[[i]], file = f1[i], row.names = FALSE, ...)
    }
  } else if (filetype == "txt"){
    if (m==1) f1 <- paste0(filename,".txt") else f1 <- paste0(filename,"_",1:m,".txt")
    for (i in 1:m){
      write.table(object$syn[[i]], file = f1[i], row.names = FALSE, ...)
    }
  } else if (filetype == "tab"){
    if (m==1) f1 <- paste0(filename,".tab") else f1 <- paste0(filename,"_",1:m,".tab")
    for (i in 1:m){
      write.table(object$syn[[i]],file=f1[i], row.names=FALSE, sep="\t", ...)
    }
  }
  
  call <- match.call()
  infofile <- paste("synthesis_info_",filename,".txt",sep="")
  #---
  sink(infofile)
  cat("Date saved:",format(Sys.time(), "%d %b %Y, %H:%M"), "\n")
  cat("Data frame with original data:", deparse(object$call$data), "\n")
  cat("Number of synthetic data sets:", m, "\n")
  cat("Output file(s):")
  if (filetype=="SPSS" | filetype=="SAS"){
    cat(paste0("(",filetype,") "))
    cat(paste(f1, f2, collapse=" "))
  } else {
    cat(paste0("(",filetype,") ",f1))
  }
  if (save.complete) {                            
    addname <- paste0("synobject_",filename,".RData")
    save(object,file=addname)
    cat("\nAdditional file: ", addname, sep="")
  }
  if (extended.info) {
    cat("\nMethods used:\n")
    print(object$method)
    cat("Seed used:",object$seed,"\n")
  }
  sink()
  #---
  
  cat("Synthetic data exported as ",filetype," file(s).",sep="")
  cat("\nInformation on synthetic data written to\n ",
      paste(getwd(),"/",infofile,sep=""),"\n")
  
}


###---- write.syn.SPSS ----------------------------------------------------

write.syn.SPSS <- function (df, datafile, codefile, varnames = names(df),
                            data.labels = NULL, ...)
{
  varlabels <- data.labels$var.lab
  vallabels <- data.labels$val.lab
  if (is.null(varnames)) varnames <- names(df)
  
  dfn <- df
  for (i in 1:ncol(dfn)){
    if (!is.null(vallabels[[varnames[i]]])){
      dfn[,i] <- mapvalues(dfn[,i], from = names(vallabels[[varnames[i]]]), 
                           to = vallabels[[varnames[i]]])
    }
  }
  write.table(dfn, file = datafile, row.names = FALSE, col.names = FALSE,
              sep = ",", quote = FALSE, na = "", eol = ",\n")
  varnames  <- gsub("[^[:alnum:]_\\$@#]", "\\.", varnames)
  if (is.null(varlabels)) varlabels <- varnames 
  dl.varnames <- varnames
  if (any(chv <- sapply(df, is.character))) {
    lengths <- sapply(df[chv], function(v) max(nchar(v)))
    if (any(lengths > 255L))
      stop("Cannot handle character variables longer than 255")
    lengths <- paste0("(A", lengths, ")")
    star <- ifelse(c(TRUE, diff(which(chv) > 1L)), " *"," ")
    dl.varnames[chv] <- paste(star, dl.varnames[chv], lengths)
  }
  cat("DATA LIST FILE=", adQuote(paste(getwd(), datafile, sep = "/")),
      " free (\",\")\n", file = codefile)
  cat("/", dl.varnames, " .\n\n", file = codefile, append = TRUE)
  cat("VARIABLE LABELS\n", file = codefile, append = TRUE)
  cat(paste(varnames, adQuote(varlabels[varnames]), "\n"), ".\n", file = codefile,
      append = TRUE)
  factors <- sapply(df, is.factor) & !sapply(vallabels[varnames], is.null)
  if (any(factors)) {
    cat("\nVALUE LABELS\n", file = codefile, append = TRUE)
    for (v in which(factors)) {
      cat("/\n", file = codefile, append = TRUE)
      cat(varnames[v], " \n", file = codefile, append = TRUE, sep = "")
      levs     <- vallabels[[varnames[v]]]
      levslabs <- names(vallabels[[varnames[v]]])
      cat(paste((levs), adQuote(levslabs), "\n", sep = " "),
          file = codefile, append = TRUE)
    }
    cat(".\n", file = codefile, append = TRUE)
  }
  cat("\nEXECUTE.\n", file = codefile, append = TRUE)
}


###---- adQuote -----------------------------------------------------------

adQuote <- function (x) paste("\"", x, "\"", sep = "")


###---- maplabs -----------------------------------------------------------
maplabs <- function (x, from, to) 
{
  if (is.factor(x)) {
    levels(x) <- maplabs(levels(x), from, to)
    return(x)
  }
  mapidx   <- match(x, from)
  mapidxNA <- is.na(mapidx)
  x[!mapidxNA] <- to[mapidx[!mapidxNA]]
  return(x)
}


###---- mapvalues -----------------------------------------------------------
mapvalues <- function (x, from, to, warn_missing = TRUE) 
{
  if (length(from) != length(to)) {
    stop("`from` and `to` vectors are not the same length.")
  }
  if (!is.atomic(x)) {
    stop("`x` must be an atomic vector.")
  }
  if (is.factor(x)) {
    levels(x) <- mapvalues(levels(x), from, to, warn_missing)
    return(x)
  }
  mapidx <- match(x, from)
  mapidxNA <- is.na(mapidx)
  from_found <- sort(unique(mapidx))
  if (warn_missing && length(from_found) != length(from)) {
    message("The following `from` values were not present in `x`: ", 
            paste(from[!(1:length(from) %in% from_found)], collapse = ", "))
  }
  x[!mapidxNA] <- to[mapidx[!mapidxNA]]
  x
}


collinear.out <- function(x, threshold = 0.99999) {
  nvar     <- ncol(x)
  x        <- data.matrix(x)
  varnames <- dimnames(x)[[2]]
  z        <- suppressWarnings(cor(x, use = "pairwise.complete.obs"))
  z[is.na(z)] <- 0    #GR-09/2016 to handle all NAs
  hit      <- abs(z) >= threshold
  mysets   <- !duplicated(hit) & rowSums(hit) > 1
  if (sum(mysets) == 0) setlist <- NULL
  else {
    setlist  <- vector("list",sum(mysets))
    for (i in 1:sum(mysets)) setlist[[i]] <- colnames(hit)[hit[mysets, , drop = FALSE][i,]]  
  } 
  return(setlist)
}

###-----compare------------------------------------------------------------
compare <- function(object, data, ...) UseMethod("compare")

###-----compare.default----------------------------------------------------
compare.default <- function(object, ...)
  stop("No compare method associated with class ", class(object), call. = FALSE)


###-----compare.synds------------------------------------------------------
compare.synds <- function(object, data, vars = NULL, msel = NULL,  
                          stat = "percents", breaks = 20, 
                          nrow = 2, ncol = 2, rel.size.x = 1,
                          utility.stats = c("pMSE", "S_pMSE", "df"),
                          utility.for.plot = "S_pMSE",  
                          cols = c("#1A3C5A","#4187BF"),
                          plot = TRUE, table = FALSE, ...){    
  
  if (is.null(data)) stop("Requires parameter data to give name of the real data.\n", call. = FALSE)
  if (!is.data.frame(data)) stop("Argument data must be a data frame.\n", call. = FALSE) 
  
  if (any(class(data) %in% c("tbl", "tbl_df"))) data <- as.data.frame(data)  
  
  if (!inherits(object, "synds")) stop("Object must have class synds.\n", call. = FALSE )                                                                    
  if (!is.null(msel) & !all(msel %in% (1:object$m))) stop("Invalid synthesis number(s).", call. = FALSE)
  if (!all(utility.stats %in% c("VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", "MabsDD", "dBhatt",
                                "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE", "df", "all"))) 
    stop('utility.stats must be set to "all" or selected from "VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", "MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE" or "df".\n', call. = FALSE)  
  
  if (!is.null(utility.for.plot) &&
      !(length(utility.for.plot) == 1 & utility.for.plot %in% c("VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", "MabsDD", "dBhatt",
                                                                "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE"))) 
    stop('utility.for.plot must be one of "VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", "MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G" or "S_pMSE" .\n', call. = FALSE)  
  if (!is.null(utility.for.plot) && !utility.for.plot %in% utility.stats ) cat(utility.for.plot, "used on plots added to table of results.\n")   
  utility.stats <- unique(c(utility.stats, utility.for.plot))
  
  if (!(length(stat) == 1 & stat %in% c("percents", "counts"))) { 
    cat('Parameter stat must be "percents" or "counts".\n' )
    stat <- "percents"
    cat('Changed to default value of "percents".\n')
  }
  
  # single / pooled synthetic data sets                                  
  if (object$m == 1) {
    syndsall <- object$syn 
  } else if (length(msel) == 1) {
    syndsall <- object$syn[[msel]]
  } else if (length(msel) > 1 | is.null(msel)) {                               
    syndsall <- do.call(rbind,object$syn)
  }
  # list of synthetic data sets for non-pooled results 
  if (length(msel) > 1) {
    synds <- vector("list",length(msel))
    for (i in 1:length(msel)) synds[[i]] <- object$syn[[msel[i]]]
  }
  synnames    <- names(syndsall)
  realnames   <- names(data)
  commonnames <- synnames[match(realnames,synnames)]
  commonnames <- commonnames[!is.na(commonnames)]
  
  if (!is.null(vars)) {
    if (!(all(vars %in% synnames))) stop("Variable(s) ", 
                                         paste0(vars[is.na(match(vars,synnames))], collapse = ", "),
                                         " not in synthetic data \n", call. = FALSE)
    if (!(all(vars %in% realnames))) stop("Variable(s) ", 
                                          paste0(vars[is.na(match(vars,realnames))], collapse = ", "),
                                          " not in observed data \n", call. = FALSE)
    commonnames <- commonnames[match(vars,commonnames)]
  }
  
  if (!(all(synnames %in% realnames))) cat("Warning: Variable(s)", 
                                           paste0(synnames[is.na(match(synnames,realnames))], collapse = ", "),
                                           "in synthetic object but not in observed data\n",
                                           " Looks like you might not have the correct data for comparison\n")
  
  if ((length(commonnames) == 0) && (typeof(commonnames) == "character"))        #! when would it apply?
    stop("None of variables selected for comparison in data", call. = FALSE)
  
  df.obs    <- data[, commonnames, drop = FALSE]
  df.synall <- syndsall[, commonnames, drop = FALSE]
  if (length(msel) > 1) {
    df.syn <- vector("list",length(msel))
    for (i in 1:length(msel)) df.syn[[i]] <- synds[[i]][, commonnames, drop = FALSE]
  }
  
  # change any numeric variables with < 6 distinct values to factors
  for (i in 1:length(commonnames)) {
    if (is.numeric(df.obs[,i]) && length(table(df.obs[,i])) < 6){
      df.obs[,i] <- as.factor(df.obs[,i])
      df.synall[,i] <- as.factor(df.synall[,i])
    }
  }
  
  num <- sapply(df.synall, is.numeric) | sapply(df.synall, is.integer)  
  fac <- sapply(df.synall, function(x) is.factor(x) | is.logical(x))  
  
  if (is.null(vars)) {
    if (object$m == 1) vars <- names(object$syn)
    else vars <- names(object$syn[[1]])
  }
  
  if (!is.null(utility.stats) | !is.null(utility.for.plot)) {
    utility.list <- as.list(1:length(vars))
    names(utility.list) <- vars
    if (utility.stats[1] == "all") utility.stats <- 
      c("VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", 
        "MabsDD", "dBhatt","S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", 
        "S_pMSE", "df")
    tab.utility <- matrix(NA, length(vars), length(unique(c(utility.stats, utility.for.plot))))
    dimnames(tab.utility) <- list(vars, unique(c(utility.stats, utility.for.plot)))
    
    for (i in 1:length(vars)) {
      utility.list[[i]] <- utility.tab(object, data, vars = vars[i])
      if (i == 1) tab.ind <- match(utility.stats, names(utility.list[[i]]))
      tab.utility[i, ] <- sapply(utility.list[[i]][tab.ind], mean)
    }
    if (!is.null(utility.for.plot)) utilvals.for.plot <- tab.utility[, match(utility.for.plot, dimnames(tab.utility)[[2]])]
  } else {
    if (is.null(utility.stats)) tab.utility <- NULL
    if (is.null(utility.for.plot)) utilvals.for.plot <- NULL
  }
  # frequency tables for factors
  if (sum(fac) > 0) {
    any.fac.na <- unlist(apply(df.obs[,fac,drop = FALSE],2,function(x) any(is.na(x))))    
    per.obs.fac <- ggfac(df.obs[, fac, drop = FALSE], anyna = any.fac.na, stat = stat)  
    if (length(msel) <= 1) {
      per.syn.facall <- ggfac(df.synall[, fac, drop = FALSE], 
                              name = "synthetic", anyna = any.fac.na, stat = stat)
      if (stat == "counts") {
        per.syn.facall$perdf$Count <- per.syn.facall$perdf$Count/object$m
        per.syn.facall$perlist <- lapply(per.syn.facall$perlist,"/",object$m)
      }
    }
    if (length(msel) > 1) {
      per.syn.fac <- vector("list",length(msel))
      for (i in 1:length(msel)) per.syn.fac[[i]] <- ggfac(df.syn[[i]][, fac, drop = FALSE], 
                                                          name = paste0("syn=", msel[i]), anyna = any.fac.na, stat = stat)  
    }
  } else {                                       
    per.obs.fac    <- NULL
    per.syn.facall <- NULL
  }
  
  # frequency tables for numeric variables
  if (sum(num) > 0) {
    cont.index <- match(colnames(df.obs[, num, drop = FALSE]), colnames(syndsall))    
    na <- object$cont.na[cont.index] 
    # to exclude from summaries if no missing in data
    any.na <- unlist(apply(df.obs[, num, drop = FALSE], 2, function(x) any(is.na(x))))  
    
    lbreaks    <- as.list(rep(breaks, length(na)))  
    
    ## get limits(=breaks) from both observed and synthetic
    #--
    df.both  <- rbind.data.frame(df.obs, df.synall)                
    per.both <- ggnum(df.both[, num, drop = FALSE], na = na, 
                      breaks = lbreaks, anyna = any.na, stat = stat) ##GR stat added 
    #--  
    per.obs.num <- ggnum(df.obs[, num, drop = FALSE], na = na, 
                         breaks = per.both$hbreaks, anyna = any.na, stat = stat) ##GR stat added
    
    if (length(msel) <= 1) {
      per.syn.numall <- ggnum(df.synall[, num, drop = FALSE], 
                              name = "synthetic", na = na, breaks = per.both$hbreaks, 
                              anyna = any.na, stat = stat) 
      if (stat == "counts") {
        per.syn.numall$perdf$Count <- per.syn.numall$perdf$Count/object$m
        per.syn.numall$perlist <- lapply(per.syn.numall$perlist,"/",object$m)
      }
    } 
    if (length(msel) > 1) {
      per.syn.num <- vector("list",length(msel))
      for (i in 1:length(msel)) per.syn.num[[i]] <- ggnum(df.syn[[i]][, num, drop = FALSE], 
                                                          name =  paste0("syn=",msel[i]), na = na, breaks = per.both$hbreaks, 
                                                          anyna = any.na, stat = stat) 
    } 
  } else {
    per.obs.num <- NULL
    per.syn.numall <- NULL
  }
  
  # data frame for plotting 
  if (length(msel) <= 1) per.fac <- rbind.data.frame(per.obs.fac$perdf, 
                                                     per.obs.num$perdf, per.syn.facall$perdf, per.syn.numall$perdf )
  
  if (length(msel) > 1) {
    per.fac <- rbind.data.frame(per.obs.fac$perdf, per.obs.num$perdf)
    for (i in 1:length(msel)) {
      if (sum(fac) > 0) temp.fac <- per.syn.fac[[i]]$perdf else temp.fac <- NULL
      if (sum(num) > 0) temp.num <- per.syn.num[[i]]$perdf else temp.num <- NULL  
      per.fac <- rbind.data.frame(per.fac, temp.fac, temp.num)
    }
  }
  per.fac$Variable <- factor(per.fac$Variable, levels = commonnames, 
                             ordered = T, exclude = NULL)
  per.fac$Value <- factor(per.fac$Value, levels = unique(per.fac$Value), 
                          ordered = T, exclude = NULL)
  
  # list of result tables
  if (length(msel) <= 1) {
    os.table.fac <- mapply(rbind, obs = per.obs.fac$perlist, 
                           syn = per.syn.facall$perlist, SIMPLIFY = FALSE)  
    os.table.num <- mapply(rbind, obs = per.obs.num$perlist, 
                           syn = per.syn.numall$perlist, SIMPLIFY = FALSE)  
  }
  if (length(msel) > 1) {
    os.table.fac <- per.obs.fac$perlist 
    os.table.num <- per.obs.num$perlist 
    for (i in 1:length(msel)) {
      if (sum(fac) > 0) temp.fac <- per.syn.fac[[i]]$perlist else temp.fac <- NULL
      if (sum(num) > 0) temp.num <- per.syn.num[[i]]$perlist else temp.num <- NULL
      os.table.fac <- mapply(rbind, os.table.fac, temp.fac, SIMPLIFY = FALSE)    
      os.table.num <- mapply(rbind, os.table.num, temp.num, SIMPLIFY = FALSE)     
    }
  }
  os.table <- c(os.table.fac, os.table.num)
  if (is.null(msel)) {
    for (i in 1:length(os.table)) dimnames(os.table[[i]])[[1]] <- c("observed","synthetic")
  } else {
    for (i in 1:length(os.table)) dimnames(os.table[[i]])[[1]] <- c("observed", paste0("syn=", msel))
  }
  Value <- Percent <- Count <- Data <- NULL    ## otherwise 'no visible binding for global variables'
  # sorts the factor labels in the right order for numeric vars
  per.fac$Value <- as.character(per.fac$Value)
  vals <- unique(per.fac$Value)
  valsnum <- unique(per.fac$Value[per.fac$Variable %in% names(num[num == TRUE])])
  valsnum.nonmiss <- sort(as.numeric(vals[vals %in% valsnum & substr(vals, 1, 4) != "miss"]))
  valsnum.nonmiss <- format(valsnum.nonmiss, scientific = FALSE,
                            trim = TRUE, drop0trailing = TRUE)
  valsnum.miss <- sort(vals[vals %in% valsnum & substr(vals, 1, 4) == "miss"])
  vals[vals %in% valsnum] <- c(valsnum.nonmiss, valsnum.miss)
  per.fac$Value <- factor(as.character(per.fac$Value), levels = vals)
  
  if (!is.null(utility.for.plot)){
    levels(per.fac$Variable) <- paste0(levels(per.fac$Variable), ": ", 
                                       utility.for.plot, " = ", round(utilvals.for.plot, 2))
    commonnames_lab <- paste0(commonnames, ": ",
                              utility.for.plot, " = ", round(utilvals.for.plot, 2))
  } else {
    commonnames_lab <- commonnames
  }
  
  # get different plots in order of data
  nperplot <- nrow*ncol
  nplots   <- ceiling(length(commonnames)/nperplot)
  plots    <- vector("list", nplots)
  tables   <- vector("list", nplots)
  
  for (i in 1:nplots) {
    min <- (i - 1)*nperplot + 1
    max <- min(length(commonnames), (i - 1)*nperplot + nperplot)
    # tables
    ttables <- vector("list", (max - min + 1))
    names(ttables) <- commonnames[min:max]
    for (j in commonnames[min:max]) {
      ttables[[j]] <- os.table[[j]]
    } 
    tables[[i]] <- ttables
    
    # plots
    per.fact <- per.fac[per.fac$Variable %in% commonnames_lab[min:max],]
    if (stat == "percents") p <- ggplot(data = per.fact, aes(x = Value, y = Percent, fill = Data))
    else p <- ggplot(data = per.fact, aes(x = Value, y = Count, fill = Data))
    p <- p + geom_bar(position = "dodge", colour = cols[1], stat = "identity") + 
      facet_wrap(~ Variable, scales = "free", ncol = ncol)  
    p <- p + guides(fill = guide_legend(override.aes = list(colour = NULL))) + 
      theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1, size = rel(rel.size.x)), 
            legend.position = "top", 
            legend.key = element_rect(colour = cols[1])) 
    p <- p + theme(legend.title = element_blank())
    if (length(msel) > 1) p <- p + scale_fill_manual(values = c(cols[1], rep(cols[2], length(msel))))
    if (length(msel) <= 1) p <- p + scale_fill_manual(values = cols)
    plots[[i]] <- p
  }
  
  if (length(tables) == 1) {
    tables <- tables[[1]]
    plots  <- plots[[1]]
  }
  
  res <- list(tables = tables, plots = plots, stat = stat, vars = vars,
              tab.utility = tab.utility, table = table, plot = plot)
  
  class(res) <- "compare.synds"
  return(res)
}


###-----compare.data.frame---compare.list----------------------------------    
compare.data.frame <- compare.list <- function(object, data, vars = NULL, cont.na = NULL,         
                                               msel = NULL, stat = "percents", breaks = 20, 
                                               nrow = 2, ncol = 2, rel.size.x = 1,
                                               utility.stats = c("pMSE", "S_pMSE", "df"),
                                               utility.for.plot = "S_pMSE",
                                               cols = c("#1A3C5A","#4187BF"),   
                                               plot = TRUE, table = FALSE, ...){
  
  if (is.null(data)) stop("Requires parameter 'data' to give name of the real data.\n\n",  call. = FALSE)
  if (is.null(object)) stop("Requires parameter 'object' to give name of the synthetic data.\n\n",  call. = FALSE)   
  
  if (is.list(object) & !is.data.frame(object)) m <- length(object)
  else if (is.data.frame(object)) m <- 1
  else stop("object must be a data frame or a list of data frames.\n", call. = FALSE)
  
  # sort out cont.na to make it into a complete named list
  cna <- cont.na
  cont.na <- as.list(rep(NA, length(data)))
  names(cont.na) <- names(data)
  if (!is.null(cna)) {
    if (!is.list(cna) | any(names(cna) == "") | is.null(names(cna))) 
      stop("Argument 'cont.na' must be a named list with names of selected variables.", call. = FALSE)  
    if (any(!names(cna) %in% names(data))) stop("Names of the list cont.na must be variables in data.\n", call. = FALSE)
    for (i in 1:length(cna)) {
      j <- (1:length(data))[names(cna)[i] == names(data)]
      cont.na[[j]] <- unique(c(NA,cna[[i]]))
    }
  }
  
  object <- list(syn = object, m = m, cont.na = cont.na)
  class(object ) <- "synds"
  
  res <- compare.synds(object = object, data = data, vars = vars, 
                       msel = msel, stat = stat, breaks = breaks, 
                       nrow = nrow, ncol = ncol, rel.size.x = rel.size.x,
                       utility.stats = utility.stats,
                       utility.for.plot = utility.for.plot,
                       cols = cols, plot = plot, table = table, ...) 
  res$call <- match.call()
  return(res)
}


###-----pertable-----------------------------------------------------------
# calculate counts/percentages
pertable <- function(x, stat, ...) { ##GR stat added
  res <- table(x, useNA = "always")   
  if (stat == "percents") res <- res*100/sum(res) 
  return(res)
}


###-----ggfac--------------------------------------------------------------
# calculate percentages for factors and store in a data frame (a long format)
ggfac <- function(data, name = "observed", anyna, stat){ 
  data <- as.data.frame(data)
  perlist  <- lapply(data, pertable, stat = stat)
  for (i in 1:length(perlist)) {
    if (anyna[i] == FALSE) perlist[[i]] <- perlist[[i]][-length(perlist[[i]])]  
  }
  if (stat == "percents") Percent  <- unlist(perlist, use.names = FALSE) 
  else Count <- unlist(perlist, use.names = FALSE)
  Value <- unlist(lapply(perlist,names), use.names = FALSE)
  Variable <- rep(names(perlist), sapply(perlist,length))
  if (stat == "percents") perdf <- cbind.data.frame(Percent, Value, Variable) 
  else  perdf <- cbind.data.frame(Count, Value, Variable)
  
  perdf$Data <- name
  return(list(perdf = perdf, perlist = perlist))
} 

###-----ggnum--------------------------------------------------------------
# calculate percentages for numeric variables (non-missing values and 
# missing data categories seperately) and store in a data frame (a long format)
ggnum <- function(data, name = "observed", na = as.list(rep(NA,ncol(data))), 
                  breaks = as.list(rep(30,ncol(data))), anyna, 
                  stat = stat){ 
  data <- as.data.frame(data)
  
  # non-missing values  
  nvar <- ncol(data)
  perlist <- vector("list", nvar)
  hbreaks <- vector("list", nvar)
  
  for (i in 1:nvar) {
    ## counts for non-missing values  
    vardata <- data[!(data[,i] %in% na[[i]]),i]                  
    hh      <- hist(vardata, breaks = breaks[[i]], plot = FALSE)
    counts  <- hh$counts
    names(counts) <- format(hh$breaks[-length(hh$breaks)],
                            scientific = FALSE, trim = TRUE,
                            drop0trailing = TRUE)
    hbreaks[[i]]  <- hh$breaks
    ## counts for missing values  
    dataNA <- data[data[,i] %in% c(NA,na[[i]]),i]
    
    if (length(dataNA) == 0) {
      if (anyna[i] == TRUE) {
        NAcounts <- 0   
        names(NAcounts) <- NA 
      } else {
        NAcounts <- NULL
      }
    } else {
      if (anyna[i] == TRUE) {
        NAcounts <- table(dataNA, useNA = "always")
      } else {
        NAcounts <- table(dataNA)  
      }
    } 
    if (!is.null(NAcounts)) names(NAcounts) <- paste("miss", names(NAcounts), sep = ".") 
    counts <- c(counts, NAcounts)
    if (stat == "percents") perlist[[i]] <- counts*100/length(data[,i])
    else perlist[[i]] <- counts
  }
  names(perlist) <- colnames(data)
  
  # create data frame in a long format  
  if (stat == "percents") Percent  <- unlist(perlist, use.names = F) 
  else Count  <- unlist(perlist, use.names = F)
  Value    <- unlist(lapply(perlist,names), use.names = F)
  Variable <- rep(names(perlist),sapply(perlist,length))
  if (stat == "percents") perdf    <- cbind.data.frame(Percent, Value, Variable)
  else perdf    <- cbind.data.frame(Count, Value, Variable)
  perdf$Data <- name
  
  return(list(perdf = perdf, perlist = perlist, hbreaks = hbreaks))
}


###-----dfNA---------------------------------------------------------------
dfNA <- function(data, na){
  all  <- length(data)
  data <- data[data %in% unlist(na)]
  NA.counts <- table(data, exclude = NULL)*100/all
  return(NA.counts)
}


###-----compare.fit.synds--------------------------------------------------
compare.fit.synds <- function(object, data, plot = "Z",
                              print.coef = FALSE, return.plot = TRUE, plot.intercept = FALSE, 
                              lwd = 1, lty = 1, lcol = c("#1A3C5A","#4187BF"),
                              dodge.height = .5, point.size = 2.5, 
                              population.inference = FALSE, ci.level = 0.95, ...) {   # c("#132B43", "#56B1F7")
  
  # Compares and plots fits to synthetic and original data
  # First parameter must be a fit to synthetic data from glm.synds(),
  # lm.synds(), polr.synds() or multinom.synds()
  
  value <- "Value"
  coefficient <- c("Coefficient", "Model")  
  
  if (!inherits(object, "fit.synds")) stop("Object must have class fit.synds.\n")
  if (!is.data.frame(data)) stop("Data must be a data frame.\n")  # theoretically can be a matrix (?)
  if (ci.level <= 0 | ci.level > 1) stop("ci.level must be beteen 0 and 1.\n")
  
  m <- object$m
  n <- object$n
  k <- object$k
  fitting.function <- object$fitting.function
  
  ##?? syn.coef    <- object$mcoefavg
  
  call <- match.call()
  
  # get fit to real data
  if (fitting.function %in% c("multinom", "polr")) {
    real.fit.0 <- do.call(object$fitting.function,
                          args = list(formula = formula(object),
                                      Hess = TRUE, data = call$data))
    real.fit <- summary(real.fit.0)
  } else if (fitting.function %in% c("lm")) {
    real.fit <- summary(do.call(object$fitting.function,
                                args = list(formula = formula(object), 
                                            data = call$data)))
  } else { ##  for glm
    real.fit <- summary(do.call(object$fitting.function,
                                args = list(formula = formula(object),
                                            family = object$call$family, 
                                            data = call$data)))
  }
  
  if (object$fitting.function == "multinom") { 
    real.varcov <- vcov(real.fit.0) 
    dd <- dimnames(t(real.fit$coefficients))
    real.fit$coefficients <- cbind(as.vector(t(real.fit$coefficients)), 
                                   as.vector(t(real.fit$standard.errors)),
                                   as.vector(t(real.fit$coefficients))/as.vector(t(real.fit$standard.errors))) 
    dimnames(real.fit$coefficients) <- list(
      paste(rep(dd[[2]], each = length(dd[[1]])), 
            rep(dd[[1]], length(dd[[2]])), sep = " : "), c("Estimate","Std error","t stat"))
  } else if (object$fitting.function == "polr") { 
    real.varcov <- vcov(real.fit.0) 
  } else if (object$fitting.function == "lm") {
    real.varcov <- real.fit$cov.unscaled * real.fit$sigma^2  ##??
  } else { # for "glm"
    real.varcov <- real.fit$cov.scaled   
  }
  
  syn.fit  <- summary.fit.synds(object, real.varcov = real.varcov,
                                population.inference = population.inference)
  incomplete <- syn.fit$incomplete
  
  # detailed results
  res.obs <- real.fit$coefficients[,1:3]
  colnames(res.obs) <- c("Beta","se(Beta)","Z")
  res.syn  <- syn.fit$coefficients[,1:3] 
  res.syn  <- res.syn[order(match(rownames(res.syn), rownames(res.obs))), ]
  res.overlap <- compare.CI(res.syn, res.obs, ci.level = ci.level, intercept = TRUE)
  ncoef <- nrow(res.obs) 
  
  res.diff <-  cbind(res.syn[,1], res.obs[,1], 
                     res.syn[,1] - res.obs[,1], 
                     (res.syn[,1] - res.obs[,1]) / res.obs[,2])
  dimnames(res.diff)[[2]] <- c("Synthetic", "Observed", "Diff", "Std. coef diff")
  #? pval <- round(2 * (1 - pnorm(abs(res.syn[,1] - res.obs[,1])/sqrt(diag(lof.varcov)))), 3)  # "p value"
  
  if (incomplete == TRUE) {   
    if (object$m < ncoef) {
      cat("\n\nWnen some variables are not synthesised m  (= ",m,") must exceed number of",
          "\ncoefficients (= ",ncoef,") for lack of fit test. No test can be reported.\n\n")  
      lack.of.fit <- NULL
      lof.pvalue  <- NULL
    } else {
      QB <- matrix(NA, m, length(object$mcoefavg))
      for (i in 1:m) {
        QB[i,] <- object$mcoef[i,] - object$mcoefavg
      }
      lof.varcov <- t(QB) %*% QB/(m - 1)/m
      lack.of.fit <- t(res.diff[,3]) %*% ginv(lof.varcov) %*% res.diff[,3]
      lack.of.fit <- lack.of.fit * (object$m - ncoef)/ncoef/(object$m - 1) # Hotellings T square
      lof.pvalue  <- 1 - pf(lack.of.fit, ncoef, object$m - ncoef)          
    }
  } else {
    lof.varcov <- real.varcov*n/k/m
    lack.of.fit <- t(res.diff[,3]) %*% ginv(lof.varcov) %*% res.diff[,3]  #! multiply by m for combined test
    lof.pvalue  <- 1 - pchisq(lack.of.fit, ncoef)
  }
  
  ##?? if (object$proper == TRUE) lof.varcov <- real.varcov * (1 + n/k)/m
  
  # Calculate summary measures
  mean.ci.overlap   <-  mean(res.overlap[,1])
  mean.abs.std.diff <-  mean(abs(res.diff[,4]))
  
  if (return.plot == TRUE) {
    yvar <- as.character(formula(object)[[2]])
    
    if (plot == "Z") {
      BetaCI <- dfCI(real.fit, Z = TRUE, ci.level = ci.level)
      
      if (population.inference == FALSE) {  ## get interval from real var
        BsynCI <- BetaCI
        for (i in c(1,3,4)) BsynCI[,i] <- BsynCI[,i] + (res.syn[,1] - res.obs[,1])/res.obs[,2]
        BsynCI[,5] <- "synthetic"
      } else {
        BsynCI <- dfCI(syn.fit, Z = TRUE, name.Z = "Z.syn", 
                       model.name = "synthetic", ci.level = ci.level)
      }
      xlab = "Z value"
      title = paste0("Z values for fit to ",yvar)
      
    } else {
      BetaCI <- dfCI(real.fit, Z = FALSE, ci.level = ci.level)  
      
      if (population.inference == FALSE) {  ## get interval from real var
        BsynCI <- BetaCI
        for (i in c(1,3,4)) BsynCI[,i] <- BsynCI[,i] + (res.syn[,1] - res.obs[,1])
        BsynCI[,5] <- "synthetic"
      } else {
        BsynCI <- dfCI(syn.fit, Z = FALSE, name.Z = "syn.coef", 
                       model.name = "synthetic", ci.level = ci.level)
      }
      
      xlab = "Value"
      title = paste0("Coefficients for fit to ",yvar)
    }
    
    modelCI <- rbind.data.frame(BetaCI, BsynCI)
    rownames(modelCI) <- 1:nrow(modelCI)
    
    if (!plot.intercept) modelCI <- modelCI[modelCI$Coefficient != "(Intercept)",]
    
    CI.geom <- geom_errorbar(aes_string(ymin = "LowCI", ymax = "HighCI",
                                        color = "Model", linetype = "Model"), data = modelCI, width = 0, 
                             lwd = lwd, lty = lty, position = position_dodge(width = dodge.height))
    
    point.geom <- geom_point(aes_string(#ymin = value, ymax = value,  #BN-03/02/2017 commented
      color = "Model", shape = "Model"), data = modelCI, 
      size = point.size, position = position_dodge(width = dodge.height))
    
    p <- ggplot(data = modelCI, aes_string(x = "Coefficient", y = "Value"))
    p <- p + geom_hline(yintercept = 0, colour = "grey", linetype = 2, lwd = 1)
    p <- p + CI.geom + point.geom + labs(title = title, y = xlab)
    p <- p + scale_shape_manual(values = c(17:16), breaks = c("synthetic","observed")) +
      scale_colour_manual(values = lcol[2:1], breaks = c("synthetic", "observed"))
    p <- p + coord_flip()
    # p <- p + theme_bw()
    # scale_colour_manual(values = rev(brewer.pal(3,"Blues")))
    # scale_colour_grey(start = 0, end = .6)
    p
  } else p <- NULL
  
  res <- list(call = object$call, coef.obs = res.obs, coef.syn = res.syn, 
              coef.diff = res.diff, mean.abs.std.diff = mean.abs.std.diff,
              ci.overlap = res.overlap, mean.ci.overlap =  mean.ci.overlap, 
              lack.of.fit = lack.of.fit, lof.pvalue = lof.pvalue, 
              ci.plot = p, print.coef = print.coef,       
              m = object$m, ncoef = ncoef,
              incomplete = incomplete,
              population.inference = population.inference)  
  
  class(res) <- "compare.fit.synds"
  return(res)
}


###-----dfCI---------------------------------------------------------------
# extract info for plotting confidence intervals
dfCI <- function(modelsummary, names.est.se = c("Estimate","Std. Error"),
                 model.name = "observed",  ci.level = 0.95, Z = FALSE, 
                 name.Z = colnames(modelsummary$coefficients)[3]){
  
  CI <- qnorm(1 - (1 - ci.level)/2)
  if (!Z) {
    #msCI <- as.data.frame(modelsummary$coefficients[,names.est.se])
    msCI <- as.data.frame(modelsummary$coefficients[,1:2]) 
    names(msCI) <- c("Value", "SE")
  } else {
    # msCI <- as.data.frame(modelsummary$coefficients[,name.Z])
    msCI <- as.data.frame(modelsummary$coefficients[,3]) 
    names(msCI) <- c("Value")
    msCI$SE <- 1
  }  
  msCI$Coefficient <- rownames(msCI)
  
  msCI$HighCI <- msCI$Value + CI*msCI$SE
  msCI$LowCI  <- msCI$Value - CI*msCI$SE
  msCI$SE     <- NULL
  msCI$Model  <- model.name
  msCI$Coefficient <- factor(msCI$Coefficient, levels = rev(msCI$Coefficient))  #!BN290416, rev added 
  
  return(msCI)
}



###-----print.synds--------------------------------------------------------

print.synds <- function(x, ...){
  cat("Call:\n($call) ")
  print(x$call)
  cat("\nNumber of synthesised data sets: \n($m) ",x$m,"\n")  
  if (x$m == 1) {
    cat("\nFirst rows of synthesised data set: \n($syn)\n")
    print(head(x$syn))
  } else {
    cat("\nFirst rows of first synthesised data set: \n($syn)\n")
    print(head(x$syn[[1]]))
  }    
  cat("...\n")
  cat("\nSynthesising methods: \n($method)\n")
  print(x$method)
  cat("\nOrder of synthesis: \n($visit.sequence)\n")
  print(x$visit.sequence)
  cat("\nMatrix of predictors: \n($predictor.matrix)\n")
  print(x$predictor.matrix)     
  invisible(x)
}


###-----summary.synds------------------------------------------------------

summary.synds <- function(object, msel = NULL, 
                          maxsum = 7, digits = max(3, getOption("digits") - 3), ...){
  if (!is.null(msel) & !all(msel %in% (1:object$m))) 
    stop("Invalid synthesis number(s)", call. = FALSE)
  
  sy <- list(m = object$m, msel = msel, method = object$method)
  
  if (object$m == 1) {
    sy$result <- summary(object$syn,...)
  } else if (is.null(msel)) {
    zall <- vector("list",object$m) 
    for (i in 1:object$m) zall[[i]] <- lapply(object$syn[[i]], summary,
                                              maxsum = maxsum, digits = digits, ...)
    zall.df <- Reduce(function(x,y) mapply("rbind",x,y),zall)
    meanres <- lapply(zall.df, function(x) apply(x,2,mean))
    sy$result <- summary.out(meanres)
  } else if (length(msel) == 1) {
    sy$result <- summary(object$syn[[msel]],...)
  } else {
    for (i in (1:length(msel))) {
      sy$result[[i]] <- summary(object$syn[[msel[i]]],...)
    }
  }
  class(sy) <- "summary.synds"
  return(sy)
}


###-----print.summary.synds------------------------------------------------

print.summary.synds <- function(x, ...){
  
  if (x$m == 1) {
    cat("Synthetic object with one synthesis using methods:\n")
    print(x$method)
    cat("\n")
    print(x$result)
  } else if (is.null(x$msel)) {
    cat("Synthetic object with ",x$m," syntheses using methods:\n",sep = "")
    print(x$method)
    cat("\nSummary (average) for all synthetic data sets:\n",sep = "")
    print(x$result)  
  } else if (length(x$msel) == 1) {
    cat("Synthetic object with ",x$m," syntheses using methods:\n",sep = "")
    print(x$method)
    cat("\nSummary for synthetic data set ",x$msel,":\n",sep = "")
    print(x$result)
  } else {
    cat("Synthetic object with ",x$m," syntheses using methods:\n",sep = "")
    print(x$method)
    for (i in (1:length(x$msel))) {
      cat("\nSummary for synthetic data set ",x$msel[i],":\n",sep = "")
      print(x$result[[i]])
    }
  }
  invisible(x)
}


###-----mcoefvar--------------------------------------------------
# Arrange coefficients from all m syntheses in a matrix
# (same with their variances). 
# [used in lm.synds and glm.synds function]

mcoefvar <- function(analyses, ...) {
  m <- length(analyses)
  if (m == 1) {
    matcoef <- mcoefavg <- analyses[[1]]$coefficients[,1]
    matvar  <- mvaravg  <- analyses[[1]]$coefficients[,2]^2
  } else {
    namesbyfit <- lapply(lapply(analyses,coefficients),rownames)
    allnames <- Reduce(union,namesbyfit)
    matcoef <- matvar <- matrix(NA, m, length(allnames))
    dimnames(matcoef)[[2]] <- dimnames(matvar)[[2]] <- allnames
    for (i in 1:m) {
      pos <- match(namesbyfit[[i]],allnames)
      matcoef[i,pos] <- analyses[[i]]$coefficients[,1]
      matvar[i,pos] <- analyses[[i]]$coefficients[,2]^2
    }
    mcoefavg <- apply(matcoef, 2, mean, na.rm = TRUE)
    mvaravg  <- apply(matvar,  2, mean, na.rm = TRUE)
    #bm <- apply(matcoef,2,var) not needed xpt for partial synthesis
  }
  if (m > 1) rownames(matcoef) <- rownames(matvar) <- paste0("syn=", 1:m)
  return(list(mcoef    = matcoef,  mvar    = matvar, 
              mcoefavg = mcoefavg, mvaravg = mvaravg))
}


###-----lm.synds-----------------------------------------------------------

lm.synds <- function(formula, data, ...)
{
  if (!inherits(data, "synds")) stop("Data must have class synds\n", call. = FALSE)
  if (is.matrix(data$method)) data$method <- data$method[1,]
  if (is.matrix(data$visit.sequence)) data$visit.sequence <- data$visit.sequence[1,]
  if (data$m > 1) vars <- names(data$syn[[1]])  else  vars <- names(data$syn)  
  n <- sum(data$n)
  if (is.list(data$k)) k <- sum(data$k[[1]]) else k <- sum(data$k)  
  
  call <- match.call()
  fitting.function <- "lm"
  analyses <- as.list(1:data$m)
  
  # Do the repeated analysis, store the result without data
  if (data$m == 1) {
    analyses[[1]] <- summary(lm(formula, data = data$syn, ...))
  } else {
    for (i in 1:data$m) {
      analyses[[i]] <- summary(lm(formula, data = data$syn[[i]], ...))
    }
  }
  
  # Check validity of inference from vars not in visit sequence or with method ""
  incomplete <- checkcomplete(vars, formula, data$visit.sequence, data$method) 
  
  # Get matrices from coefficients
  allcoefvar <- mcoefvar(analyses = analyses)
  
  # Return the complete data analyses as a list of length m
  object <- list(call = call, mcoefavg = allcoefvar$mcoefavg, 
                 mvaravg = allcoefvar$mvaravg, analyses = analyses,  
                 fitting.function = fitting.function,
                 n = n, k = k, proper = data$proper, 
                 m = data$m, method = data$method, incomplete = incomplete,
                 mcoef = allcoefvar$mcoef, mvar = allcoefvar$mvar)
  class(object) <- "fit.synds"
  return(object)
}


###-----glm.synds----------------------------------------------------------

glm.synds <- function(formula, family = "binomial", data, ...)
{
  if (!inherits(data, "synds")) stop("Data must have class synds\n", call. = FALSE)
  if (is.matrix(data$method)) data$method <- data$method[1,]
  if (is.matrix(data$visit.sequence)) data$visit.sequence <- data$visit.sequence[1,]
  if (data$m > 1) vars <- names(data$syn[[1]])  else  vars <- names(data$syn)  
  n <- sum(data$n)
  if (is.list(data$k)) k <- sum(data$k[[1]]) else k <- sum(data$k)  
  
  call <- match.call()
  fitting.function <- "glm"
  analyses <- as.list(1:data$m)
  
  # Do the repeated analysis, store the result without data
  if (data$m == 1) {
    analyses[[1]] <- summary(glm(formula,data = data$syn, family = family, ...))
  } else {
    for (i in 1:data$m) {
      analyses[[i]] <- summary(glm(formula,data = data$syn[[i]], family = family, ...))
    }
  }
  
  # Check completeness for inference from vars not in visit sequence or with method ""
  incomplete <- checkcomplete(vars, formula, data$visit.sequence, data$method) 
  
  # Get matrices from coefficients
  allcoefvar <- mcoefvar(analyses = analyses)
  
  # Return the complete data analyses as a list of length m
  object <- list(call = call, mcoefavg = allcoefvar$mcoefavg, 
                 mvaravg = allcoefvar$mvaravg, analyses = analyses,  
                 fitting.function = fitting.function,
                 n = n, k = k, proper = data$proper, 
                 m = data$m, method = data$method, incomplete = incomplete,
                 mcoef = allcoefvar$mcoef, mvar = allcoefvar$mvar)
  class(object) <- "fit.synds"
  return(object)
}


###-----polr.synds-----------------------------------------------------

polr.synds <- function(formula, data, ...)
{
  if (!inherits(data, "synds")) stop("Data must have class 'synds'.\n", call. = FALSE)
  if (is.matrix(data$method)) data$method <- data$method[1,]
  if (is.matrix(data$visit.sequence)) data$visit.sequence <- data$visit.sequence[1,]
  if (data$m > 1) vars <- names(data$syn[[1]]) else  vars <- names(data$syn)  
  n <- sum(data$n)
  if (is.list(data$k)) k <- sum(data$k[[1]]) else k <- sum(data$k)  
  
  call <- match.call()
  fitting.function <- "polr"
  analyses <- as.list(1:data$m)
  
  # Do the repeated analysis, store the result without data
  for (i in 1:data$m) {
    if (data$m == 1) fit <- polr(formula, data = data$syn, Hess = TRUE, ...)
    else fit <- polr(formula, data = data$syn[[i]], Hess = TRUE, ...)
    ss <- summary(fit)
    analyses[[i]] <- ss
  }
  
  # Check validity of inference from vars not in visit sequence or with method ""
  incomplete <- checkcomplete(vars, formula, data$visit.sequence, data$method) 
  
  # Get matrices from coefficients
  allcoefvar <- mcoefvar(analyses = analyses)
  
  # Return the complete data analyses as a list of length m
  object <- list(call = call, mcoefavg = allcoefvar$mcoefavg, 
                 mvaravg = allcoefvar$mvaravg, analyses = analyses,  
                 fitting.function = fitting.function,
                 n = n, k = k, proper = data$proper, 
                 m = data$m, method = data$method, incomplete = incomplete,
                 mcoef = allcoefvar$mcoef, mvar = allcoefvar$mvar)
  class(object) <- "fit.synds"
  return(object)
}


###-----multinom.synds-----------------------------------------------------

multinom.synds <- function(formula, data, ...)
{
  if (!inherits(data, "synds")) stop("Data must have class 'synds'.\n", call. = FALSE)
  if (is.matrix(data$method)) data$method <- data$method[1,]
  if (is.matrix(data$visit.sequence)) data$visit.sequence <- data$visit.sequence[1,]
  if (data$m > 1) vars <- names(data$syn[[1]]) else  vars <- names(data$syn)  
  n <- sum(data$n)
  if (is.list(data$k)) k <- sum(data$k[[1]]) else k <- sum(data$k)  
  
  call <- match.call()
  fitting.function <- "multinom"
  analyses <- as.list(1:data$m)
  
  # Do the repated analysis, store the result without data
  for (i in 1:data$m) {
    if (data$m == 1) fit <- multinom(formula, data = data$syn, Hess = TRUE, ...)
    else fit <- multinom(formula, data = data$syn[[i]], Hess = TRUE, ...)
    ss <- summary(fit)
    analyses[[i]] <- list(coefficients = cbind(as.vector(t(ss$coefficients)),
                                               as.vector(t(ss$standard.errors)),
                                               as.vector(t(ss$coefficients)/t(ss$standard.errors))))
    dd <- dimnames(t(ss$coefficients))
    dimnames(analyses[[i]]$coefficients) <- list(paste(rep(dd[[2]], each = length(dd[[1]])),
                                                       rep(dd[[1]], length(dd[[2]])), sep = ":"),
                                                 c("Estimate", "se", "z value"))
  }
  
  # Check validity of inference from vars not in visit sequence or with method ""
  incomplete <- checkcomplete(vars, formula, data$visit.sequence, data$method) 
  
  # Get matrices from coefficients
  allcoefvar <- mcoefvar(analyses = analyses)
  
  # Return the complete data analyses as a list of length m
  object <- list(call = call, mcoefavg = allcoefvar$mcoefavg, 
                 mvaravg = allcoefvar$mvaravg, analyses = analyses,  
                 fitting.function = fitting.function,
                 n = n, k = k, proper = data$proper, 
                 m = data$m, method = data$method, incomplete = incomplete,
                 mcoef = allcoefvar$mcoef, mvar = allcoefvar$mvar)
  class(object) <- "fit.synds"
  return(object)
}


###-----checkcomplete------------------------------------------------------
# Used in lm.synds, glm.synds, multinom.synds, and polr.synds

checkcomplete <- function(vars, formula, vs, method)
{  
  inform <- all.vars(formula) # get all variables in formula
  if ("." %in% inform) inform <- vars
  if (any(!inform %in% vars)) stop("Variable(s) in formula (model to be fitted) are not in synthetic data: ",
                                   paste(inform[!inform %in% vars], collapse = ", "), call. = FALSE)
  methin <- method[names(method) %in% inform]
  if (all(methin == "")) cat("No variables in your formula (model to be fitted) have been synthesised.\nIf the data contain exactly the same observations as the original ones,\nresults will be identical.\n")
  
  order_vs  <- match(names(vs), names(methin))
  order_vs  <- order_vs[!is.na(order_vs)]
  order_oth <- setdiff(1:length(methin), order_vs)
  methin_order <- methin[c(order_vs, order_oth)] 
  blankmeths <- (1:length(methin_order))[methin_order == ""]
  
  if (!all(blankmeths == (1:length(blankmeths)))){ 
    cat(
      "**********************************************************",
      "\nWARNING: Some variable(s) in formula (model to be fitted)  
are not synthesised, so not used in synthesising models
for other variables:", 
      paste(names(methin_order)[blankmeths][!(blankmeths == (1:length(blankmeths)))], collapse = ", "), 
      "\nMethods in synthesis order are:\n")
    print(methin_order)
    cat("Results may not be correct.
**********************************************************\n")
  }
  
  incomplete <- methin[1] == ""
  return(incomplete)
}


###-----print.fit.synds----------------------------------------------------

print.fit.synds <- function(x, msel = NULL, ...)
{
  if (!is.null(msel) & !all(msel %in% (1:x$m))) stop("Invalid synthesis number(s): `msel` must be selected from 1:", x$m, call. = FALSE, sep = "")
  
  #if (x$n != x$k | x$m > 1) cat("Note: To get a summary of results you would expect from the original data\nor for population inference use the summary() function on your fit.\nSee vignette on inference to get more details.\n") 
  if (x$n != x$k | x$m > 1) cat("Note: To get more details of the fit see vignette on inference.\n") 
  
  cat("\nCall:\n")
  print(x$call)
  if (is.null(msel) & x$m > 1) {
    cat("\nAverage coefficient estimates from", x$m, "syntheses:\n")
    print(x$mcoefavg)
  } else if (x$m == 1) {
    cat("\nCoefficient estimates from a single synthesis:\n")
    print(x$mcoefavg)
  } else {
    cat("\nCoefficient estimates for selected synthetic data set(s):\n")
    print(x$mcoef[msel, , drop = FALSE])
  }
  invisible(x)
}


###-----summary.fit.synds--------------------------------------------------

summary.fit.synds <- function(object, population.inference = FALSE, msel = NULL, 
                              real.varcov = NULL, incomplete = NULL, ...)
{ # df.residual changed to df[2] because didn't work for lm 
  if (!inherits(object, "fit.synds")) stop("Object must have class fit.synds\n", call. = FALSE)
  m <- object$m
  n <- object$n
  k <- object$k
  if (is.null(incomplete)) incomplete <- object$incomplete
  
  coefficients <- object$mcoefavg  # mean of coefficients (over m syntheses)
  if (!is.null(real.varcov))  vars <- diag(real.varcov)
  else  vars <- object$mvaravg * k/n  # mean of variances (over m syntheses) * adjustment
  
  ## Checks, messages and warnings for population inference
  #---
  if (population.inference == TRUE) {
    if (incomplete == TRUE & m == 1) {
      cat("Warning: You have selected population inference when your dependent variable is not synthesised",
          "\nor incomplete is set to TRUE and when only a single synthetic data set has been created (m = 1).",
          "\nThe correct method for this case requires m > 1, ideally m > 5.",
          "\nTo provide some results calculations proceed as if all variables had been synthesised.\n\n")
      incomplete <- FALSE}
    else if (incomplete == TRUE & m < 5) {
      cat("Note: You have selected population inference for incompletely synthesised data with m = ", m,
          ",\nwhich is smaller than the minimum of 5 recommended. The estimated standard errors of your\ncoefficients may be inaccurate.\n\n", sep = "")
    }
  }
  #--- 
  
  ## Inference to Q hat
  #---
  if (population.inference == FALSE) { 
    result <- cbind(coefficients,
                    sqrt(vars),
                    coefficients/sqrt(vars),
                    2*pnorm(-abs(coefficients/sqrt(vars))))
    colnames(result) <- c("xpct(Beta)", "xpct(se.Beta)", "xpct(z)", "Pr(>|xpct(z)|)")
    #--- 
    
    ## Population inference to Q
    #---   
  } else { 
    
    ## incomplete method  
    if (incomplete == TRUE) {
      bm <- apply(object$mcoef, 2, var)
      result <- cbind(coefficients,
                      sqrt(bm/m + vars),
                      coefficients/sqrt(bm/m + vars),
                      2*pnorm(-abs(coefficients/sqrt(bm/m + vars))))
      
      ## simple synthesis   
    } else {
      if (object$proper == FALSE) Tf <- vars*(1 + n/k/m) else Tf <- vars*(1 + (n/k + 1)/m)
      result <- cbind(coefficients, 
                      sqrt(Tf), 
                      coefficients/sqrt(Tf),
                      2*pnorm(-abs(coefficients/sqrt(Tf)))) 
    }
    colnames(result) <- c("Beta.syn","se.Beta.syn","z.syn","Pr(>|z.syn|)")
  }
  #---
  
  res <- list(call = object$call, proper = object$proper,
              population.inference = population.inference,
              incomplete = incomplete, 
              fitting.function = object$fitting.function,
              m = m, coefficients = result, n = n, k = k, 
              analyses = object$analyses, msel = msel)
  class(res) <- "summary.fit.synds"
  return(res)
}


###-----print.summary.fit.synds--------------------------------------------

print.summary.fit.synds <- function(x, ...) {
  
  if (!is.null(x$msel) & !all(x$msel %in% (1:x$m))) stop("Invalid synthesis number(s)", call. = FALSE)
  #cat("\nNote that all these results depend on the synthesis model being correct.\n")  
  #cat("\nFor details, see package vignette on inference.\n")
  
  if (x$m == 1) {
    cat("Fit to synthetic data set with a single synthesis. ")
  } else {
    cat("Fit to synthetic data set with ", x$m, " syntheses. ", sep = "")
  }
  
  if (x$population.inference) {
    if (x$incomplete == TRUE) cat("Inference to population coefficients when\nsome variables in the model are not synthesised. Methods for incomplete/partial\nsynthesis are used.\n")
    else cat("Inference to population\ncoefficients when all variables in the model are synthesised.\n")
  } else {
    cat("Inference to coefficients\nand standard errors that would be obtained from the original data.\n")
    if (x$k != x$n) 
      cat("\nThe synthetic data have a different size (", x$k, ") from the original data (", x$n, "),",
          "\nso the standard errors of the coefficients have been adjusted to estimate",
          "\nthe standard errors from the original data.\n", sep = "")
  }
  
  cat("\nCall:\n")
  print(x$call)
  cat("\nCombined estimates:\n")
  printCoefmat(x$coefficients)
  
  if (!is.null(x$msel)) {
    allcoef <- lapply(lapply(x$analyses[x$msel], "[[", "coefficients"), as.data.frame)
    
    estimates <- lapply(allcoef, "[", "Estimate")
    allestimates <- do.call(cbind, estimates)
    
    zvalues <- lapply(allcoef, "[", "z value")
    allzvalues <- do.call(cbind, zvalues)
    
    colnames(allestimates) <- colnames(allzvalues) <- paste0("syn=",x$msel)
    
    cat("\nEstimates for selected syntheses contributing to the combined estimates:\n")
    
    cat("\nCoefficients:\n")
    print(allestimates)
    cat("\nz values:\n")
    print(allzvalues)
    
    # for(i in x$msel) {          
    #   cat("\nsyn=",i,"\n",sep = "")
    #   print(x$analyses[[i]]$coefficients)
    # }
  }      
  invisible(x)
}


###-----print.compare.fit.synds--------------------------------------------

print.compare.fit.synds <- function(x, print.coef = x$print.coef, ...){
  
  cat("\nCall used to fit models to the data:\n")
  print(x$call)
  if (print.coef == TRUE) {
    cat("\nEstimates for the observed data set:\n")
    print(x$coef.obs)
    cat("\nCombined estimates for the synthesised data set(s):\n")
    print(x$coef.syn)
  }  
  
  cat("\nDifferences between results based on synthetic and observed data:\n")
  print(cbind.data.frame(x$coef.diff,x$ci.overlap))                                 
  if (x$m == 1) {
    cat("\nMeasures for one synthesis and ", x$ncoef, " coefficients", sep = "") 
  } else {
    cat("\nMeasures for ", x$m, " syntheses and ", x$ncoef, " coefficients", sep = "") 
  }   
  cat("\nMean confidence interval overlap: ", x$mean.ci.overlap)
  cat("\nMean absolute std. coef diff: ", x$mean.abs.std.diff)
  if (!is.null(x$lack.of.fit)){
    cat("\n\nMahalanobis distance ratio for lack-of-fit (target 1.0): ", 
        round(x$lack.of.fit/x$ncoef, 2), sep = "")
    cat("\nLack-of-fit test: ", x$lack.of.fit,"; p-value ", round(x$lof.pval, 4), 
        " for test that synthesis model is\ncompatible ", sep = "")
    if (x$incomplete == FALSE) cat("with a chi-squared test with ", 
                                   x$ncoef, " degrees of freedom.\n", sep = "")
    else cat("with an F distribution with ", x$ncoef, " and ", 
             x$m - x$ncoef, " degrees of freedom.\n", sep = "") 
  }
  if (!is.null(x$ci.plot)) {
    cat("\nConfidence interval plot:\n")
    print(x$ci.plot)
  }
  invisible(x)
}


###-----print.compare.synds------------------------------------------------

print.compare.synds <- function(x, ...) {
  if (x$table | x$plot) {
    if (x$stat == "counts") cat("\nComparing counts observed with synthetic\n\n") 
    else cat("\nComparing percentages observed with synthetic\n\n")
  }
  
  if (class(x$plots)[1] == "gg") {
    if (x$table) print(x$tables) 
    if (x$plot) print(x$plots)
  } else {
    for (i in 1:length(x$tables)) {
      if (x$table) print(x$tables[[i]]) 
      if (x$plot) {
        print(x$plots[[i]])
        if (i < length(x$tables)) {
          cat("Press return for next variable(s): ")
          ans <- readline()
        }
      }
    }
  }
  if (!is.null(x$tab.utility)) {
    cat("\nSelected utility measures:\n")  
    print(round(x$tab.utility, 6))
  }
  #if (x$print.tab.utility) print(x$tab.utility)
  invisible(x)
}


###-----print.utility.gen--------------------------------------------------

print.utility.gen <- function(x, digits = NULL, zthresh = NULL,    
                              print.zscores = NULL, print.stats = NULL,
                              print.ind.results = NULL,                               
                              print.variable.importance = NULL, ...){
  
  if (x$method == "logit"){  
    converged <- ifelse(x$m == 1, x$fit$converged, x$fit[[1]]$converged)
    length.y <- ifelse(x$m == 1, length(x$fit$y)/2, length(x$fit[[1]]$y)/2)
    length.coef <- ifelse(x$m == 1, length(x$fit$coefficients), length(x$fit[[1]]$coefficients)) 
    if (converged == FALSE)  cat("Warning: Logistic model did not converge in ", x$maxit,
                                 " iterations. Check results.\nThis can be due to cells with small numbers or too complicated a model.
In first case utility measures are probably fine but coefficients may be unstable. 
In the latter case you may be able to overcome this by increasing the parameter 'maxit'.
Your combined data has ", length.y, " observations and your model has ",
                                 length.coef, " parameters.\n", sep = "")
  }
  
  if (is.null(digits)) digits <- x$digits
  if (is.null(print.stats)) print.stats <- x$print.stats
  if (is.null(zthresh)) zthresh <- x$zthresh
  if (is.null(print.zscores)) print.zscores <- x$print.zscores
  if (is.null(print.ind.results)) print.ind.results <- x$print.ind.results   
  if (is.null(print.variable.importance)) print.variable.importance <- x$print.variable.importance
  
  cat("\nUtility score calculated by method: ", x$method, "\n", sep = "")
  cat("\nCall:\n")
  print(x$call)
  #?  
  if (!is.null(x$resamp.method) && !x$resamp.method == "none") {
    if (x$resamp.method == "perm") cat("\nNull utilities simulated from a permutation test with ", x$nperm," replications.\n", sep = "")
    else if (x$resamp.method == "pairs") cat("\nNull utilities simulated from ", x$m*(x$m - 1)/2," pair(s) of syntheses.\n", sep = "")
    
    if (!is.list(x$nnosplits)) { 
      if (x$nnosplits[1] > 0) cat(
        "\n***************************************************************
Warning: null pMSE resamples failed to split ", x$nnosplits[1], " time(s) from ", x$nnosplits[2],
        "\n***************************************************************\n", sep = "")
    } else {
      for (ss in 1:x$m) {
        if (x$nnosplits[[ss]][1] > 0) cat("\nSynthesis ", ss, 
                                          "   null pMSE resamples failed to split ", x$nnosplits[[ss]][1],
                                          " time(s) from ", x$nnosplits[[ss]][2], sep = "")
      }
      cat("\n")
    }
  }
  
  if (print.stats[1] == "all") print.stats <- 
    c("pMSE", "S_pMSE", "SPECKS", "S_SPECKS", "PO50", "S_PO50", "U", "S_U")
  
  if (x$m == 1) { 
    cat("\nSelected utility measures\n")  
    tabres <- c(round(x$pMSE, digits), round(x$S_pMSE, digits), 
                round(x$SPECKS, digits), round(x$S_SPECKS, digits),
                round(x$PO50, digits), round(x$S_PO50, digits),
                round(x$U, digits), round(x$S_U, digits))
    names(tabres) <- c("pMSE", "S_pMSE", "SPECKS", "S_SPECKS", "PO50", "S_PO50", "U", "S_U")
    tabres <- tabres[match(print.stats, names(tabres))]
    print(tabres)
  }
  else {
    cat("\nMean utility results from ", x$m, " syntheses:\n", sep = "")
    tabres <- data.frame(round(x$pMSE, digits), round(x$S_pMSE, digits), 
                         round(x$SPECKS, digits), round(x$S_SPECKS, digits),
                         round(x$PO50, digits), round(x$S_PO50, digits),
                         round(x$U, digits), round(x$S_U, digits))
    names(tabres) <- c("pMSE", "S_pMSE", "SPECKS", "S_SPECKS", "PO50", "S_PO50", "U", "S_U")
    tabres <- tabres[, match(print.stats, names(tabres)), drop = FALSE]
    print(round(apply(tabres, 2, mean), digits))
    if (print.ind.results == TRUE) {
      cat("\nIndividual utility results from ", x$m, " syntheses:\n", sep = "")
      print(tabres)
    }
  }
  
  if (print.zscores == TRUE) {
    if (x$method == "cart") {
      cat("\nz-scores are not available for CART models.\n")
    } else {
      if (x$m > 1) {
        allzscores <- vector("list", x$m) 
        for (i in 1:x$m) allzscores[[i]] <- summary(x$fit[[i]])$coefficients[ ,3] 
        allnames <- unique(unlist(lapply(allzscores, names)))
        allzscores.NA <- lapply(allzscores, "[", allnames) 
        allzscores.NA.df <- do.call(cbind, allzscores.NA)
        zscores <- apply(allzscores.NA.df, 1, mean, na.rm = TRUE)  
        names(zscores) <- allnames
      } else {
        zscores <- summary(x$fit)$coefficients[, 3]
      }
      
      if (!is.na(zthresh)) { 
        zscores <- zscores[abs(zscores) > zthresh]
        if (length(zscores) == 0) {
          cat("\nNo z-scores (or mean z-scores if m > 1) above threshold of +/-", zthresh,"\n", sep = "")
        } else {
          cat("\nz-scores (or mean z-scores if m > 1) greater than the threshold of +/- ", zthresh, "\n", sep = "")
          print(round(zscores, digits))
        }  
      } else {
        cat("\nAll z-scores (or mean z-scores if m > 1)\n")
        print(round(zscores, digits))
      }
    }
  }
  
  if (print.variable.importance == TRUE) {
    if (x$method != "cart" | x$tree.method != "rpart") {
      cat("\nVariable importance only available for CART models using function 'rpart' (tree.method).\n")
    } else {
      cat("\nRelative importance of each variable scaled to add to 100\n" )
      if (x$m == 1) {
        variable.importance <- x$fit$variable.importance
        variable.importance <- round(variable.importance/sum(variable.importance)*100, digits)
        print(round(variable.importance), digits)
      } else {
        cat("(results for ", x$m, " syntheses)\n", sep = "")
        variable.importance <- vector("list", x$m)
        for (i in 1:x$m) {
          if (is.null(x$fit[[i]]$variable.importance)) x$fit[[i]]$variable.importance <- NA
          variable.importance[[i]] <- x$fit[[i]]$variable.importance
          variable.importance[[i]] <- round(variable.importance[[i]]/sum(variable.importance[[i]])*100, digits)
        }
        allnames <- unique(unlist(lapply(variable.importance, names)))
        all.vars <- lapply(variable.importance, "[", allnames) 
        all.vars.importance <- do.call(rbind, all.vars)
        colnames(all.vars.importance) <- allnames
        rownames(all.vars.importance) <- 1:x$m
        print(round(all.vars.importance, digits))
      }
    }
  }
  invisible(x)
  
}


###-----print.utility.tab--------------------------------------------------

print.utility.tab <- function(x, print.tables = NULL,
                              print.zdiff = NULL, 
                              print.stats = NULL, 
                              digits = NULL, ...){
  
  if (is.null(print.stats)) print.stats <- x$print.stats
  if (is.null(print.zdiff)) print.zdiff <- x$print.zdiff
  if (is.null(print.tables)) print.tables <- x$print.tables
  if (is.null(digits)) digits <- x$digits
  
  if (print.tables == TRUE) {
    if (x$k.syn) cat("\nObserved not adjusted to match the size of the synthetic data since 'k.syn' set to TRUE:\n($tab.obs)\n")
    else if (sum(x$tab.obs) != x$n & !x$k.syn) { 
      cat("\nObserved adjusted to match the size of the synthetic data:\n($tab.obs)\n")
      print(round(x$tab.obs, digits))
    } else {
      cat("\nObserved:\n($tab.obs)\n")
      print(x$tab.obs)
    }  
    
    if (x$m == 1) {
      cat("\nSynthesised: \n($tab.syn)\n")
      print(x$tab.syn) 
    } else {
      if (length(dim(simplify2array(x$tab.syn))) == 2) {
        meantab <- apply(simplify2array(x$tab.syn), 1, mean)
      } else {
        meantab <- apply(simplify2array(x$tab.syn), c(1, 2), mean)
      }
      cat("\nMean of ", x$m, " synthetic tables ($tab.syn):\n", sep = "")
      print(round(meantab, digits))
    }
  }
  
  if (print.zdiff == TRUE) {
    if (x$m == 1) {
      cat("\nTable of z-scores for differences: \n($tab.zdiff)\n")      
      print(round(x$tab.zdiff, digits)) 
    } else {
      if (length(dim(simplify2array(x$tab.zdiff))) == 2) {
        meanzdiff <- apply(simplify2array(x$tab.zdiff), 1, mean)
      } else {
        meanzdiff <- apply(simplify2array(x$tab.zdiff), c(1, 2), mean)
      }
      cat("\nMean of ", x$m, " z-score tables for differences ($tab.zdiff):\n", sep = "")
      print(round(as.table(meanzdiff), digits))
    }
  }
  if (print.stats[1] == "all") print.stats <- c("VW", "FT", "JSD", 
                                                "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", "MabsDD", "dBhatt",
                                                "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE", "df", "dfG")
  
  if (x$k.syn) celldiff <- 0 else celldiff <- 1
  tab.res <- matrix(NA, x$m, length(print.stats))
  dimnames(tab.res) <- list(1:x$m, print.stats)
  if (x$m > 1) cat("\nSelected utility measures from ", x$m, " syntheses:\n", sep = "")
  else cat("\nSelected utility measures:\n")
  for (i in 1:length(print.stats)) {
    tab.res[, i] <- unlist(x[match(print.stats[i], names(x))])
  } 
  print(round(tab.res, digits))
  invisible(x)
}


###-----summary.out--------------------------------------------------------
summary.out <- function(z, digits = max(3L, getOption("digits") - 3L), ...)
{
  ncw <- function(x) {
    zz <- nchar(x, type = "w")
    if (any(na <- is.na(zz))) {
      zz[na] <- nchar(encodeString(zz[na]), "b")
    }
    zz
  }
  nv <- length(z)
  nm <- names(z)
  lw <- numeric(nv)
  nr <- if (nv)
    max(unlist(lapply(z, NROW)))
  else 0
  for (i in seq_len(nv)) {
    sms <- z[[i]]
    if (is.matrix(sms)) {
      cn <- paste(nm[i], gsub("^ +", "", colnames(sms),
                              useBytes = TRUE), sep = ".")
      tmp <- format(sms)
      if (nrow(sms) < nr)
        tmp <- rbind(tmp, matrix("", nr - nrow(sms),
                                 ncol(sms)))
      sms <- apply(tmp, 1L, function(x) paste(x, collapse = "  "))
      wid <- sapply(tmp[1L, ], nchar, type = "w")
      blanks <- paste(character(max(wid)), collapse = " ")
      wcn <- ncw(cn)
      pad0 <- floor((wid - wcn)/2)
      pad1 <- wid - wcn - pad0
      cn <- paste0(substring(blanks, 1L, pad0), cn, substring(blanks,
                                                              1L, pad1))
      nm[i] <- paste(cn, collapse = "  ")
      z[[i]] <- sms
    }
    else {
      sms <- format(sms, digits = digits)
      lbs <- format(names(sms))
      sms <- paste0(lbs, ":", sms, "  ")
      lw[i] <- ncw(lbs[1L])
      length(sms) <- nr
      z[[i]] <- sms
    }
  }
  if (nv) {
    z <- unlist(z, use.names = TRUE)
    dim(z) <- c(nr, nv)
    if (anyNA(lw))
      warning("Probably wrong encoding in names(.) of column ",
              paste(which(is.na(lw)), collapse = ", "))
    blanks <- paste(character(max(lw, na.rm = TRUE) + 2L),
                    collapse = " ")
    pad <- floor(lw - ncw(nm)/2)
    nm <- paste0(substring(blanks, 1, pad), nm)
    dimnames(z) <- list(rep.int("", nr), nm)
  }
  else {
    z <- character()
    dim(z) <- c(nr, nv)
  }
  attr(z, "class") <- c("table")
  z
}

###-----multi.compare------------------------------------------------------
multi.compare <- function(object, data, var = NULL, by = NULL, msel = NULL, 
                          barplot.position = "fill", cont.type = "hist", y.hist = "count", 
                          boxplot.point = TRUE, binwidth = NULL, ...) {
  
  # CHECKS
  #---  
  if (is.null(data)) stop("Requires parameter data to give name of the real data.\n", call. = FALSE)
  if (!is.data.frame(data)) stop("Argument data must be a data frame.\n", call. = FALSE)
  if (!is.null(msel) & !all(msel %in% (1:object$m))) stop("Invalid synthesis number(s).", call. = FALSE)
  if (is.null(var)) stop("Argument var is missing.\n", call. = FALSE)
  
  
  if (length(var) > 1) {
    cat("\nParameter var set to ", var, ", should be a single variable.\nOnly first variable used.\n") 
    var <- var[1]
  }
  if (!(var %in% names(data))) stop("\nArgument var set to ",var,
                                    ", should be the name of a variable in data\n", call. = FALSE)
  if (!all(by %in% names(data))) stop("\nArgument by set to ", by,
                                      ", should all be names of variables in data.\n", call. = FALSE)
  if (any(var %in% by)) {  
    cat("\nArgument var: ",var,"included in argument by: ", by,
        " now removed from by.\n")
    by <- by[-match(var, by)]
  }
  #--- 
  
  for (i in match(by,names(data))) data[,i] <- addNA(data[,i], ifany = TRUE)     # numeric will be changed to factor (create seperate data for this?)
  bytable <- xtabs(as.formula(paste("~", paste(by, collapse = " + "))), 
                   data = data, exclude = NULL, na.action = na.pass) 
  cat("\nPlots of ",var," by ",by,"\nNumbers in each plot (observed data):\n\n")
  print(bytable)
  nplots <- length(bytable) 
  if (nplots > 100) cat("\nCAUTION: You have ", nplots, " sections in your plot.\n")
  
  add <- NULL 
  
  # data prep
  #----
  # if (is.null(msel) & object$m > 1) msel <- 1:object$m
  if (object$m == 1) {
    synall <- cbind(object$syn, source = "syn")
  } else if (length(msel) == 1) {
    synall <- cbind(object$syn[[msel]], source = "syn")
  } else if (is.null(msel) & object$m > 1) {
    synall <- cbind(do.call(rbind,object$syn), source = "syn")
  } else if (object$m > 1 & length(msel) > 1) {
    synall <- Map(cbind, object$syn[msel], source = paste0("syn=", msel))
    synall <- do.call(rbind, synall)
  }
  obssyn <- rbind(cbind(data, source = "obs"), synall)
  if ("syn" %in% unique(obssyn$source)){
    obssyn$source <- factor(obssyn$source, levels = c("obs", "syn")) 
  } else {
    obssyn$source <- factor(obssyn$source, levels = c("obs", paste0("syn=", msel))) 
  }
  #----
  
  # change any numeric variables with < 6 distinct values to factors
  if (is.numeric(data[,var]) && length(table(data[,var])) < 6) {
    obssyn[, var] <- as.factor(obssyn[, var])
    data[, var] <- as.factor(data[, var])
  }
  
  ..count.. <- ..density.. <- NULL
  
  if (is.numeric(data[, var])) {
    if (cont.type == "hist") {
      p <- ggplot(data = obssyn, aes(x = eval(parse(text = var))))
      plabs <- labs(x = var, fill = "")
      if (y.hist == "count") {
        ptype <- geom_histogram(aes(y = ..count.., fill = source), 
                                position = "dodge", binwidth = binwidth)
      } else if (y.hist == "density"){
        ptype <- geom_histogram(aes(y = ..density.., fill = source), 
                                position = "dodge", binwidth = binwidth)
      }
    } else if (cont.type == "boxplot") {
      p <- ggplot(data = obssyn, aes(x = source, y = eval(parse(text = var))))
      ptype <- geom_boxplot(aes(colour = source), alpha = 0.7)
      plabs <- labs(y = var, x = "", colour = "")
      if (boxplot.point) add <- geom_jitter(size = 0.2, alpha = 0.2)
    }
  } else {
    if (barplot.position == "dodge") {
      p <- ggplot(data = obssyn, aes(x = eval(parse(text = var))))
      ptype <- geom_bar(aes(fill = source), position = barplot.position)
      plabs <- labs(x = var, fill = "")      
    } else {
      p <- ggplot(data = obssyn, aes(x = source))
      ptype <- geom_bar(aes(fill = eval(parse(text = var))), position = barplot.position)
      plabs <- labs(x = "", fill = var)
    }
  }
  
  if (length(by) == 1){
    playaout <- facet_wrap(by)   # nrow = 1
  } else if (length(by) > 1) {
    form <- paste(by[1], "~", paste0(by[-1], collapse = "+"))
    playaout <- facet_grid(eval(parse(text = form)))
  }
  
  p <- p + add + ptype + playaout + plabs 
  
  if (length(msel > 5)) {
    p <- p + theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1))
  } 
  # + scale_fill_brewer(palette = "Set1") + scale_colour_brewer(palette = "Set1")
  
  return(p)
}

###-----numtocat.syn-------------------------------------------------------
# group numeric variables in a data frame into categories

numtocat.syn <- function(data, numtocat = NULL, 
                         print.flag = TRUE, cont.na = NULL, 
                         catgroups = 5, style.groups = "quantile")
{
  if (!is.data.frame(data)) stop("Data must be a data.frame.\n", call. = FALSE)
  varnames <- names(data)
  
  # checks on numtocat
  if (!is.null(numtocat)) {
    if (is.numeric(numtocat)) {
      if (!all(numtocat %in% 1:length(data))) stop("Column indices must be between 1 and ", 
                                                   length(data), ".", sep = "", call. = FALSE)  
      varnos <- numtocat
      numtocat <- names(data)[varnos]
    } else {
      if (!all(numtocat %in% varnames)) stop("Variable(s) ", 
                                             paste(numtocat[!numtocat %in% varnames], collapse = ", "),
                                             " not in data.\n", sep = "", call. = FALSE)
      varnos <- match(numtocat,varnames)
    }
    vclass <- sapply(data[, varnos, drop = FALSE], is.numeric)
    if (!all(vclass)) stop("Variable(s) in numtocat (", 
                           paste(numtocat[!vclass], collapse = ", "), 
                           ") not numeric.\n", sep = "", call. = FALSE)
  } else { 
    ## if NULL use all numeric variables
    varnos   <- (1:length(data))[sapply(data, is.numeric)]
    numtocat <- names(data)[varnos]
  }
  
  # checks on catgroups
  if (length(catgroups) == 1) catgroups <- rep(catgroups,length(numtocat))
  else if (length(catgroups) != length(numtocat)) stop("catgroups must be a single number or a vector of the same length as numtocat.\n", call. = FALSE)
  
  if (any(catgroups < 2)) stop("catgroups must all be > 1.", call. = FALSE)
  # checks on cont.na
  if (!is.null(cont.na)) {
    if (!is.list(cont.na)) stop("cont.na must be a  list.\n", call. = FALSE)
    if (!all(names(cont.na) %in% numtocat)) stop("Variable(s): ",
                                                 paste(names(cont.na)[!names(cont.na) %in% numtocat],collapse = ", "),
                                                 " in cont.na not in the variables being grouped.\n", call. = FALSE)
    cna <- as.list(rep(NA,length(numtocat)))
    for (i in 1:length(cont.na)) cna[[match(names(cont.na)[i],numtocat)]] <- cont.na[[i]]
  } else {
    cna <- as.list(rep(NA,length(numtocat)))
  }
  names(cna) <- numtocat
  
  if (print.flag == TRUE) cat("Variable(s) ", paste(numtocat, collapse = ", "),
                              " grouped into categories.\n", sep = "")
  breaks <- vector("list", length(varnos)); names(breaks) <- numtocat
  levels <- vector("list", length(varnos)); names(levels) <- numtocat
  orig <- data[, varnos, drop = FALSE]
  names(orig) <- paste("orig", names(orig), sep = ".")
  for (i in 1:length(varnos)) {
    grpd <- group_var(data[, varnos[i]], cont.na = cna[[i]], n = catgroups[i], style = style.groups)
    data[, varnos[i]] <- grpd$x
    breaks[[i]] <- grpd$breaks
    levels[[i]] <- grpd$levels
  } 
  
  return(list(data = data, breaks = breaks, 
              levels = levels, orig = orig, cont.na = cna, 
              numtocat = numtocat, ind = varnos))
}


###-----group_var----------------------------------------------------------

group_var <- function(x,  n = 5, style = "quantile", cont.na = NA) {
  # categorise one continous variable into groups
  if (!is.numeric(x)) stop("x must be numeric.\n", call. = FALSE)
  # select non-missing(nm) values 
  xnm <- x[!(x %in% cont.na) & !is.na(x)]
  my_breaks   <- unique(classIntervals(xnm, n = n, style = style)$brks)
  xnm_grouped <- cut(xnm, breaks = my_breaks, dig.lab = 8, 
                     right = FALSE, include.lowest = TRUE)
  my_levels   <- c(levels(xnm_grouped), cont.na[!is.na(cont.na)])
  # assign groupings to non-missing data
  x[!(x %in% cont.na) & !is.na(x)] <- as.character(xnm_grouped)
  x <- as.factor(x)
  list(x = x, breaks = my_breaks, levels = my_levels)  
}  


padMis.syn <- function(data, method, predictor.matrix, visit.sequence,
                       nvar, rules, rvalues, default.method, cont.na, 
                       smoothing, event, denom) {
  
  # Function called by syn to make dummy/factor variable for missing values
  # in continuous variables. Data is augmented by columns for dummy/factor 
  # variables when they are used in synthesis. 
  
  # Check presence of missing values not covered by missing rules
  # (missing values for non-numeric variables are not counted)   
  No.NA <- vector("list", nvar)
  yes.rules <- sapply(rules, function(x) any(x != ""))
  com.rules <- lapply(rules, paste, collapse = " | ")
  for (j in 1:nvar) {
    if (yes.rules[j]) {
      No.NA[j] <- with(data,sum(data[!eval(parse(text = com.rules[[j]])), j] %in% cont.na[[j]]))      
    } else {
      No.NA[j] <- sum(data[,j] %in% cont.na[[j]]) 
    }
  }
  No.NA    <- sapply(No.NA,function(x) x > 0)
  inpred   <- apply(predictor.matrix != 0, 1, any) | apply(predictor.matrix != 0, 2, any)
  factorNA <- rep(FALSE, nvar)
  
  for (j in 1:nvar) {
    # if (No.NA[j] & is.numeric(data[,j]) & inpred[j]==TRUE){  #GR1.7-1
    if (No.NA[j] & is.numeric(data[,j]) & inpred[j] == TRUE & 
        !is.passive(method[j]) & !method[j] %in% c("nested", "ipf", "catall")) {
      
      # augment the data with a column for the original continuous variable with 
      # missing values replaced by zeros and a column for a new factor for 
      # missing values 
      nonmiscode <- 10^(nchar(round(max(data[,j], na.rm = TRUE))) + 1) - 1               
      y.0  <- ifelse(data[,j] %in% c(cont.na[[j]], rvalues[[j]]), 0, data[,j])
      y.NA <- ifelse(data[,j] %in% c(cont.na[[j]], rvalues[[j]]), data[,j], nonmiscode) #BN13/11 0 changed with nonmiscode
      y.NA <- addNA(y.NA, ifany = TRUE) 
      levels(y.NA)[is.na(levels(y.NA))] <- "NAtemp"   # to allow random forest
      data <- cbind(data,y.0,y.NA)           
      name.0  <- paste(attr(data,"names")[j], 0, sep = ".")
      name.NA <- paste(attr(data,"names")[j], NA, sep = ".")
      
      names(data)[(ncol(data) - 1):ncol(data)] <- c(name.0, name.NA)
      factorNA[(ncol(data) - 1):ncol(data)] <- c(FALSE,TRUE) 
      
      # predictor.matrix is given two extra rows and columns for the new variables
      # rows and columns are copied from an original variable j in predictor.matrix
      predictor.matrix <- rbind(predictor.matrix,  matrix(rep(predictor.matrix[j,], 
                                                              times = 2), byrow = TRUE, nrow = 2))
      predictor.matrix <- cbind(predictor.matrix, matrix(rep(predictor.matrix[,j], 
                                                             times = 2), ncol = 2))
      # the original variable is removed from predictors (=insert zeros) 
      predictor.matrix[,j] <- 0
      # the original variable is synthesized passively so its predictors can be removed as well
      predictor.matrix[j,] <- 0
      
      # add methods for new variables
      method[ncol(data) - 1] <- method[j]
      if (method[j] %in% c("ctree", "ctree.proper", "cart", "cart.proper",
                           "rf", "ranger", "bag",
                           "sample", "")) {
        method[ncol(data)] <- method[j] 
      } else {
        method[ncol(data)] <- ifelse(nlevels(data[,ncol(data)]) == 2,
                                     default.method[2], default.method[3])
      }   
      
      # pass smoothing to new variable and remove from original one
      smoothing[ncol(data) - 1] <- smoothing[j]
      smoothing[ncol(data)]   <- ""
      smoothing[j]            <- ""  
      
      cont.na[ncol(data) - 1] <- cont.na[j]
      cont.na[ncol(data)]     <- NA
      cont.na[j]              <- NA  
      
      # pass denom and event for new variable and remove from original one
      denom[ncol(data) - 1] <- denom[j]    #!!!!!!!!!!!!!! check if correct
      denom[ncol(data)]     <- 0 
      denom[j]              <- 0
      event[ncol(data) - 1] <- event[j]    
      event[ncol(data)]     <- 0 
      event[j]              <- 0
      
      # insert the column numbers for the new variables into the visit sequence 
      # before the jth column
      if (any(visit.sequence == j)) {                    
        newcols <- c(ncol(data), ncol(data) - 1)
        idx <- (1:length(visit.sequence))[visit.sequence == j] - 1
        visit.sequence <- append(visit.sequence,newcols,idx)
        # modify method for the original variable
        method[j] <- paste0("~(ifelse(",name.NA,"!=", nonmiscode," | is.na(",name.0,
                            "),as.numeric(levels(",name.NA,"))[",name.NA,"],",name.0,"))")    
      }
      
      # update missing rules and values for the new variables
      if (any(rules[[j]] != "")) rules[[ncol(data) - 1]] <-
        c(rules[[j]][rules[[j]] != ""],paste(name.NA, "!=", nonmiscode, sep = ""))  #BN13/11
      else rules[[ncol(data) - 1]] <- paste(name.NA, "!=", nonmiscode, sep = "")    #BN13/11 
      rules[[ncol(data)]]        <- rules[[j]]                      
      rules[[j]]                 <- ""
      #!BN1513 rule "year_death.NA!=0" should have only one correspnding 
      #! rvalue equal to 0; before a vector c(NA,0) was assigned instead of 0 
      if (length(rules[[j]]) == 1) rvalues[[ncol(data) - 1]] <- 0           
      else rvalues[[ncol(data) - 1]] <- c(rvalues[[j]], 0)                 
      rvalues[[ncol(data)]]        <- rvalues[[j]]
    }
  }
  
  varnames <- dimnames(data)[[2]]  
  dimnames(predictor.matrix) <- list(varnames,varnames)
  names(method) <- varnames
  names(visit.sequence) <- varnames[visit.sequence]
  return(list(data = as.data.frame(data), 
              nvar = ncol(data),
              predictor.matrix = predictor.matrix, 
              method = method, 
              visit.sequence = visit.sequence, 
              rules = rules,
              rvalues = rvalues,
              factorNA = factorNA,
              smoothing = smoothing,
              event = event,
              denom = denom,
              cont.na = cont.na))
}

padModel.syn <- function(data, method, predictor.matrix, visit.sequence,
                         nvar, rules, rvalues, factorNA, smoothing, event, denom) {
  
  # Function called by syn to make dummy variables data frame data is 
  # augmented by columns for dummy variables when they are used as predictors. 
  # This is returned as part of a list which also contains structures 
  # to tell sampler.syn how to select columns
  
  categories <- data.frame(yes.no.categorical = factor(rep(FALSE, nvar), 
                                                       levels = c("TRUE", "FALSE")), 
                           number.of.dummies  = rep(0, nvar), 
                           yes.no.dummy       = factor(rep(FALSE, nvar), 
                                                       levels = c("TRUE", "FALSE")), 
                           corresponding.column.dummy = rep(0, nvar))
  
  # This is a data frame with a row for each variable and extra rows for 
  # each of the dummy variables added at the end of the j loop with
  # col1 = T/F if category and predictor in parametric model                   #!BN1605
  # col2 = number of dummy variables for factors with col1==TRUE, else 0       #!BN1605
  # col3 = yes.no.dummy TRUE for dummy variables
  # col4 = corresponding.column.dummy original column number for dummies, else 0
  
  # pred.with.cart <- method %in% c("ctree", "ctree.proper", "cart", "cart.proper", "collinear", "satcat") #!BN to check if collinear needed
  
  pred.with.cart <- !method %in% c("norm", "normrank", "logreg", "lognorm", 
                                   "polr", "polyreg", "cubertnorm", "sqrtnorm")  #!GR050318
  for (j in 1:nvar) {
    if ((is.factor(data[,j]) & any(predictor.matrix[1:nvar,j] != 0 & !pred.with.cart)) |  #!BN-16/05/2016
        (factorNA[j] == TRUE & !pred.with.cart[j])) {                                      #!BN-16/05/2016
      categories[j, 1] <- TRUE
      
      # all factors defined to have treatment contrasts
      #!data[, j] <- C(data[, j], contr.treatment)                              UNCOMMENT???? BN-28/04/2016
      n.dummy   <- length(levels(data[, j])) - 1
      categories[j, 2] <- n.dummy
      
      # predictor.matrix is given extra rows and columns for the dummy variables
      # rows are set to zero initially
      predictor.matrix <- rbind(predictor.matrix, matrix(0,
                                                         ncol = ncol(predictor.matrix), nrow = n.dummy))
      
      # columns are set to zero and then for vars with non-CART method        #!BN1605
      # copied from an original variable j in predictor.matrix for            #!BN1605
      # -> 1 for all the rows for which this variable is being used as        #!BN1605
      # a predictor in a non-CART model, 0 otherwise                          #!BN1605
      predictor.matrix <- cbind(predictor.matrix, matrix(0, ncol = n.dummy,     #!BN1605
                                                         nrow = nrow(predictor.matrix)))                  #!BN1605
      predictor.matrix[!pred.with.cart,(ncol(predictor.matrix) - n.dummy + 1):ncol(predictor.matrix)] <- matrix(rep(predictor.matrix[!pred.with.cart,j], times = n.dummy)) #!BN1605
      
      # the original categorical variable is removed from predictors (=insert zeros)
      # for variables with non-CART method 
      predictor.matrix[!pred.with.cart,j] <- 0                                   #!BN1605
      
      
      # insert the column number for first of this set of dummies into
      # the visit sequence immediately after the jth column is predicted
      if (any(visit.sequence == j)) {
        # set an original categorical variable as predictor for its dummies  
        predictor.matrix[(ncol(predictor.matrix) - n.dummy + 1):ncol(predictor.matrix), j] <- rep(1, times = n.dummy)
        # insert dummies into visit sequence 
        newcol <- ncol(predictor.matrix) - n.dummy + 1
        nloops <- sum(visit.sequence == j)
        for (ii in 1:nloops) {
          idx <- (1:length(visit.sequence))[visit.sequence == j][ii]
          visit.sequence <- append(visit.sequence, newcol, idx)
        }
      }
      
      # augment the data with columns for the new dummies
      data <- (cbind(data, matrix(0, ncol = n.dummy, nrow = nrow(data))))
      
      # set dummies to missing when variable is missing
      data[is.na(data[, j]), (ncol(predictor.matrix) - n.dummy + 1):ncol(predictor.matrix)] <- NA
      cat.column <- data[!is.na(data[, j]), j]  # these are the non missing values of this factor
      
      # next bit sets the colums for the dummies to the dummy variables 
      # when data are not missing and labels columns
      data[!is.na(data[, j]), (ncol(predictor.matrix) - n.dummy + 1):ncol(predictor.matrix)] <- model.matrix(~cat.column - 1)[,-1]
      names(data)[(ncol(predictor.matrix) - n.dummy + 1):ncol(predictor.matrix)] <- paste(attr(data,"names")[j],(1:n.dummy), sep = ".")
      method     <- c(method, rep("dummy", n.dummy))
      rules      <- c(rules,rep(rules[j],n.dummy))
      rvalues    <- c(rvalues,rep(rvalues[j],n.dummy))
      categories <- rbind(categories, 
                          data.frame(yes.no.categorical = rep(FALSE,n.dummy),
                                     number.of.dummies = rep(0, n.dummy), 
                                     yes.no.dummy      = rep(TRUE, n.dummy), 
                                     corresponding.column.dummy = rep(j,n.dummy)))
    }
  }                                                       
  
  varnames <- dimnames(data)[[2]]  # now includes dummy names
  dimnames(predictor.matrix) <- list(varnames, varnames)
  names(method) <- varnames
  names(visit.sequence) <- varnames[visit.sequence]
  dimnames(categories)[[1]] <- dimnames(data)[[2]]
  
  return(list(data = as.data.frame(data), 
              syn = as.data.frame(data),
              predictor.matrix = predictor.matrix, 
              method = method, 
              visit.sequence = visit.sequence, 
              rules = rules,
              rvalues = rvalues, 
              categories = categories,
              smoothing = smoothing,
              event = event,
              denom = denom))
}

#-----------------------------sampler.syn-------------------------------
# The sampler controls the generation of conditional distributions
# This function is called by syn()

sampler.syn <- function(p, data, m, syn, visit.sequence,
                        rules, rvalues, event, proper,
                        print.flag, k, pred.not.syn, 
                        models, numtocat,  ...)
{
  #--- Assign optional parameters (...) to appropriate synthesising function   
  dots  <- as.list(substitute(list(...)))[-1L]         
  meth.with.opt <- paste(c("cart", "cartbboot", "ctree", "survctree", "polyreg", 
                           "norm", "lognorm", "sqrtnorm", "cubertnorm", "normrank", "pmm",
                           "polr", "rf", "ranger", "bag", "ipf", "catall"), collapse = "\\.|")
  meth.check <- grep(meth.with.opt, names(dots), value = TRUE)
  args.err <- !(names(dots) %in% meth.check)
  if (any(args.err)) stop("Unknown optional parameter(s): ", 
                          paste(names(dots)[args.err], collapse = ", "),
                          "\nNote that they have to be method specific, e.g. 'ctree.minbucket' and NOT 'minbucket'\n", 
                          call. = FALSE)
  if (length(dots) == 0) {
    mth.args <- NULL
  } else {  
    #mth.args.dots <- strsplit(names(dots), "\\.")
    mth.args.dots <- regmatches(names(dots), regexpr("\\.", names(dots)), invert = TRUE)
    mth.dots  <- unique(lapply(mth.args.dots, "[[", 1))
    args.dots <- lapply(mth.args.dots, "[[", -1)
    mth.args  <- setNames(vector("list", length(mth.dots)), unlist(mth.dots))
    
    for (i in 1:length(mth.dots)) { 
      ind <- grep(mth.dots[[i]], names(dots))
      mth.args[[i]] <- setNames(dots[ind], args.dots[ind])
    } 
  } 
  #---
  
  fits <- NULL 
  
  if (m > 0) {
    if (models) fits <- rep(list(setNames(vector("list", length(p$method)),
                                          names(p$method))), m) 
    for (i in 1:m) {  # Synthesising loop
      if (print.flag & m > 1) cat("\nSynthesis number ", i, 
                                  "\n--------------------\n", sep = "")  
      if (print.flag & m == 1) cat("\nSynthesis\n-----------\n", sep = "")  
      
      # Code for methods that take more than one variable together: ipf & catall      
      #--------------------------------------------------------------------------
      rest.visit.sequence <- p$visit.sequence  # when no grouped methods used
      
      if (any(p$method %in% c("catall", "ipf"))) {
        ordmethod <- p$method[p$visit.sequence]
        grind <- (1:length(p$visit.sequence))[ordmethod %in% ordmethod[1]]
        
        ## to reorder any dummies for grouped variables
        if (any(names(p$visit.sequence) %in% 
                paste(names(p$visit.sequence[grind]), "1", sep = "."))) {  
          dumind <- (1:length(p$visit.sequence))[names(p$visit.sequence) %in% 
                                                   paste(names(p$visit.sequence[grind]), "1", sep = ".")]
          othind <- (1:length(p$visit.sequence))[-c(grind, dumind)]
          p$visit.sequence <- p$visit.sequence[c(grind, dumind, othind)]
          ordmethod <- p$method[p$visit.sequence]
        }
        
        grouped <- p$visit.sequence[ordmethod %in% ordmethod[1]]
        
        if (print.flag == TRUE) {
          if (length(rest.visit.sequence) > 0  && 
              ncol(data) - length(numtocat) > length(grouped)) {
            
            cat("First ", length(grouped), " variables (", 
                paste(names(grouped), collapse = ", "),
                ") synthesised together by method '", ordmethod[1], "'\n", sep = "")
            if (ordmethod[1] == "catall" && !is.null(mth.args) && 
                "epsilon" %in% names(mth.args$catall) && mth.args$catall$epsilon > 0)
              cat("Synthesis made differentially private with parameter epsilon of ",
                  mth.args$catall$epsilon,"\n",
                  "Note that only these first variables will be made differentially private.\n")
            if (ordmethod[1] == "ipf" && !is.null(mth.args) && 
                "epsilon" %in% names(mth.args$ipf) && mth.args$ipf$epsilon > 0) 
              cat("Synthesis made differentially private with parameter epsilon of ",
                  mth.args$ipf$epsilon,"\n",
                  "Note that only these first variables will be made differentially private.\n")
          } else {
            cat("All ", length(grouped), 
                " variables in the data synthesised together by method '", 
                ordmethod[1], "'\n", sep = "")
            
            if (ordmethod[1] == "catall" && !is.null(mth.args) && 
                "epsilon" %in% names(mth.args$catall) && 
                mth.args$catall$epsilon > 0) 
              cat("Synthesis made differentially private with parameter epsilon of ",
                  mth.args$catall$epsilon,"\n")
            if (ordmethod[1] == "ipf" && !is.null(mth.args) && 
                "epsilon" %in% names(mth.args$ipf) && 
                mth.args$ipf$epsilon > 0) 
              cat("Synthesis made differentially private with parameter epsilon of ",
                  mth.args$ipf$epsilon,"\n")
          }   
        }   
        x <- p$data[, grouped]
        if (!(ordmethod[1] %in% names(mth.args))) fun.args <- NULL else
          fun.args  <- mth.args[[ordmethod[1]]]
        f <- paste("syn", ordmethod[1], sep = ".")
        synfun <- do.call(f, args = c(list(x = x, k = k, 
                                           proper = proper), fun.args))
        p$syn[, grouped] <- synfun$res
        if (models) {
          fits[[i]][[grouped[1]]] <- synfun$fit
          for (j in 2:length(grouped)) fits[[i]][[grouped[j]]] <- 
              paste("See first in group:", names(grouped)[1])
        }
        
        rest.visit.sequence <- p$visit.sequence[-(1:length(grouped))]
        if (length(rest.visit.sequence) > 0 & print.flag & 
            ncol(data) - length(numtocat) > length(grouped)) cat("\nRemaining variables:\n")
      }
      
      # Other variables 
      #--------------------------------------------------------------------------
      if (length(rest.visit.sequence) > 0) {           
        prcount <- 0 # to get new lines come out right
        for (j in rest.visit.sequence) {
          
          theMethod <- p$method[j]
          # get optional parameters for theMethod if they are provided
          if (!(theMethod %in% names(mth.args))) fun.args <- NULL else            
            fun.args  <- mth.args[[theMethod]]                                     
          
          vname <- dimnames(p$data)[[2]][j]
          if (print.flag & theMethod != "dummy"  
              & j <= (ncol(data) - length(numtocat))) { 
            cat(" ", vname, sep = "")
            prcount <- prcount + 1
          }  
          if (print.flag & prcount %% 10 == 0 & 
              j <= (ncol(data) - length(numtocat))) cat("\n")                                                  
          
          ya  <-  1:nrow(p$data) 
          ypa <- 1:k    
          
          # ya = yavailable, ym = ymissing                                            
          if (any(p$rules[[j]] != "")) {
            com.rules  <- paste(p$rules[[j]], collapse = " | ")
            evalrul.y  <- with(p$data,eval(parse(text = com.rules)))
            ym         <- which(evalrul.y == TRUE & !is.na(evalrul.y))
            ya         <- setdiff(1:nrow(p$data), ym)                                  
            evalrul.yp <- with(p$syn,eval(parse(text = com.rules)))         
            ypm        <- which(evalrul.yp == TRUE & !is.na(evalrul.yp))        
            ypa        <- setdiff(1:nrow(p$syn), ypm)       
          }                                                                       
          
          # != "", != "dummy", != "passive"
          if (theMethod != "" & (!is.passive(theMethod)) & theMethod != "dummy" ) {
            
            if (theMethod %in% c("sample", "sample.proper", "constant")) {
              
              y   <- p$data[ya, j]
              if (is.factor(y)) y <- y[, drop = TRUE]
              xp  <- length(ypa)
              x   <- length(ya)
              nam <- vname
              f   <- paste("syn", theMethod, sep = ".")
              if (theMethod == "constant") {
                synfun <- do.call(f, args = list(y = y, xp = xp, ...))    
              } else if (is.numeric(y)) {
                synfun <- do.call(f, args = list(y = y, xp = xp,       
                                                 smoothing = p$smoothing[j], cont.na = p$cont.na[[j]], 
                                                 proper = proper, ...)) 
              } else {
                synfun <- do.call(f, args = list(y = y, xp = xp, 
                                                 proper = proper, ...)) 
              }
              p$syn[ypa, j]  <- synfun$res
              if (models) fits[[i]][[j]] <- synfun$fit
              
            } else {
              
              x    <- p$data[ya, p$predictor.matrix[j, ] == 1, drop = FALSE]
              xp   <- p$syn[ypa, p$predictor.matrix[j, ] == 1, drop = FALSE]
              y    <- p$data[ya, j]
              if (is.factor(y)) y <- y[, drop = TRUE]
              nam  <- vname
              f    <- paste("syn", theMethod, sep = ".") 
              if (!theMethod %in% c("collinear", "nested")) {   # nested needs added to allow missing values 
                #if(theMethod!="collinear"){                  
                keep <- remove.lindep.syn(x, y, ...)
                x    <- x[, keep, drop = FALSE]
                xp   <- xp[, keep, drop = FALSE]
              }                                            
              if (theMethod == "survctree") {
                if (p$event[j] == -1) yevent <- rep(1,length(y))                   
                else yevent  <- p$data[ya,p$event[j]]
                survres      <- do.call(f, args = c(list(y = y, yevent = yevent,
                                                         x = x, xp = xp, proper = proper), 
                                                    fun.args))
                p$syn[ypa, j] <- survres[[1]]                                # synthetic data survival goes to p$syn
                if (p$event[j] != -1) p$syn[ypa,p$event[j]] <- survres[[2]]  # synthetic data event goes to p$syn
                if (models) fits[[i]][[j]] <- survres$fit 
              } else if (theMethod == "logreg" & p$denom[j] != 0) {                 
                synfun <- do.call(f, args = list(y = y, x = x, xp = xp,
                                                 denom = p$data[ya,p$denom[j]], 
                                                 denomp = p$syn[ypa, p$denom[j]],       
                                                 proper = proper, ...))
                p$syn[ypa, j] <- synfun$res
                if (models) fits[[i]][[j]] <- synfun$fit           
                
              } else if (theMethod == "nested") {
                if (is.numeric(y)) {
                  synfun <- do.call(f, args = c(list(y = y, x = x, xp = xp,
                                                     smoothing = p$smoothing[j], cont.na = p$cont.na[[j]], 
                                                     proper = proper), fun.args))
                } else {
                  synfun <- do.call(f, args = c(list(y = y, x = x, xp = xp,
                                                     proper = proper), fun.args))
                }
                p$syn[ypa, j] <- synfun$res
                if (models) fits[[i]][[j]] <- synfun$fit
                
              } else {
                if (is.numeric(y)) {
                  synfun <- do.call(f, args = c(list(y = y, x = x, xp = xp,
                                                     smoothing = p$smoothing[j],
                                                     proper = proper), fun.args))
                } else {
                  synfun <- do.call(f, args = c(list(y = y, x = x, xp = xp,
                                                     proper = proper), fun.args))
                }
                p$syn[ypa, j] <- synfun$res
                if (models) fits[[i]][[j]] <- synfun$fit
              }
            }
            
            if (any(p$rules[[j]] != "")) {
              if (length(p$rules[[j]]) == 1 & length(ypm) > 0) {
                p$syn[ypm,j] <- p$rvalues[[j]] 
              } else {
                for (r in 1:length(p$rules[[j]])) {
                  revalrul.yp  <- with(p$syn,eval(parse(text = p$rules[[j]][r])))  
                  rypm <- which(revalrul.yp == TRUE & !is.na(revalrul.yp))
                  if (length(rypm) > 0) p$syn[rypm,j] <- p$rvalues[[j]][r]
                }
              }                 
            }  
          } # end of !="", !="dummy", !="passive"
          
          else if (is.passive(theMethod)) {
            class0 <- class(p$syn[,j])
            synfun <- syn.passive(data = p$syn, func = theMethod) 
            
            if (is.factor(synfun$res[[1]]) & any(is.na(synfun$res[[1]]))) {
              synfun$res[[1]] <- addNA(synfun$res[[1]], ifany = TRUE)
              levels(synfun$res[[1]])[is.na(levels(synfun$res[[1]]))] <- "NAtemp"
            }
            
            p$syn[, j] <- synfun$res
            class(p$syn[,j]) <- class0
            if (models) fits[[i]][[j]] <- synfun$fit 
          }
          
          else if (theMethod == "dummy") {    # replace dummy variables in p$syn
            # getting dummy values from a synthesised categorical variable
            cat.columns <- p$syn[, p$categories[j, 4]]  # this is the single column with the data for which this is the dummy
            model.frame(~cat.columns - 1, data = p$syn) 
            p$syn[, (j:(j + p$categories[p$categories[j, 4], 2] - 1))] <- # replaces all the dummies for this variable with
              matrix((model.matrix(~cat.columns - 1)[, -1]),                # dummies calculated from the synthesised data
                     ncol = p$categories[p$categories[j, 4], 2],
                     nrow = nrow(p$syn))
            p$syn[,j] <- as.numeric(p$syn[, j])
            remove("cat.columns")
            if (models) fits[[i]][[j]] <- "dummy"    
          }
        } # end j loop 
      } # end other variables
      if (print.flag) cat("\n")  
      
      #if (k==dim(data)[1]) syn[[i]] <- p$syn[,1:dim(data)[2]]
      #else syn[[i]] <- p$syn[sample(1:dim(data)[1],k),1:dim(data)[2]]
      
      syn[[i]] <- p$syn[, 1:dim(data)[2], drop = FALSE]
      nms <- names(data)
      # exclude unsynthesised if drop.pred.only set to true
      if (sum(pred.not.syn ) > 0) {
        syn[[i]] <- syn[[i]][, !pred.not.syn]
        nms <- nms[!pred.not.syn]  # GR save names to use below if data just one column
      }
      # Prevent a single character column being changed to a factor
      chgetochar <- (sum(!pred.not.syn) == 1 & any(class(syn[[i]][, 1]) == "character"))       
      syn[[i]] <- as.data.frame(syn[[i]])
      if (chgetochar) {
        syn[[i]][, 1] <- as.character(syn[[i]][, 1])
        names(syn[[i]]) <- nms
      }
      
      #turn NA level in factors / logical to missing NA's
      # and remove contrasts 
      for (j in (1:ncol(syn[[i]]))) {
        if (is.factor(syn[[i]][,j])) {                                    #!BN-20/04/16
          if ("NAlogical" %in% levels(syn[[i]][,j])) {
            levels(syn[[i]][,j])[levels(syn[[i]][,j]) == "NAlogical"] <- NA
            syn[[i]][,j] <- as.logical(syn[[i]][,j])
          } else {                                                            
            # syn[[i]][,j] <- factor(syn[[i]][,j],exclude=NA,levels=levels(syn[[i]][,j]))
            levels(syn[[i]][,j])[levels(syn[[i]][,j]) == "NAtemp"] <- NA  #!BN 10/08/15 
          }
          #! attributes(syn[[i]][,j])$contrasts <- NULL                    #!BN-28/04/16   UNCOMMENT???? 
        }                                                                       
      }
    } # end i loop (m)
  } # end synthesising (m > 0)
  
  return(list(syn = syn, fits = fits))
}


###-----remove.lindep.syn--------------------------------------------------

remove.lindep.syn <- function(x, y, eps = 0.00001, maxcor = 0.99999, 
                              allow.na = FALSE, ...) 
{
  if (ncol(x) == 0) return(NULL) 
  if (eps <= 0) stop("\n Argument 'eps' must be positive.", call. = FALSE)
  xobs <- sapply(x, as.numeric)                                       
  yobs <- as.numeric(y)
  keep <- unlist(apply(xobs, 2, var) > eps)
  keep[is.na(keep)] <- FALSE
  keep <- keep & suppressWarnings((unlist(apply(xobs, 2, cor, yobs)) < maxcor)) # if y includes NA -> NAs error
  if (all(!keep)) warning("\nAll predictors are constant or have too high correlation.\n")
  ksum <- sum(keep)
  cx   <- cor(xobs[, keep, drop = FALSE], use = "all.obs")
  eig  <- eigen(cx, symmetric = TRUE)
  ncx  <- cx
  while (eig$values[ksum]/eig$values[1] < eps) {
    j   <- (1:ksum)[order(abs(eig$vectors[, ksum]), decreasing = TRUE)[1]]
    keep[keep][j] <- FALSE
    ncx  <- cx[keep[keep], keep[keep], drop = FALSE]
    ksum <- ksum - 1
    eig  <- eigen(ncx)
  }
  # if (!all(keep)) cat("\tVariable(s): ", paste(dimnames(x)[[2]][!keep], collapse = ", "),
  #   " removed due to linear dependency",sep="")
  return(keep)
}

# Strata can be provided as a vector (=flag) or as variable names 
# If all stratifying values have the same values within strata they will be
# automatically removed; if not they will be taken into account

# ????? 
#! probably NOT possible and syn.strata needs to be kept as a separate FUN
# creating generic syn function
# changing syn to syn.default
# and if strata != NULL calling syn.strata
# ?????

# Other issues:
# - multiple printout (per strata, per m)
# - compare() by strata?
# - k as a sum?

syn.strata <- function(data, strata = NULL, 
                       minstratumsize = 10 + 10 * length(visit.sequence), 
                       tab.strataobs = TRUE, tab.stratasyn = FALSE,
                       method = "cart", visit.sequence = (1:ncol(data)),
                       predictor.matrix = NULL,
                       m = 1, k = nrow(data), proper = FALSE,
                       minnumlevels = 1, maxfaclevels = 60,
                       rules = NULL, rvalues = NULL,
                       cont.na = NULL, semicont = NULL,
                       smoothing = NULL, event = NULL, denom = NULL,
                       drop.not.used = FALSE, drop.pred.only = FALSE,
                       default.method = c("normrank","logreg","polyreg","polr"),
                       numtocat = NULL, catgroups = rep(5,length(numtocat)), 
                       models = FALSE,
                       print.flag = TRUE,
                       seed = "sample",
                       ...){
  
  m0 <- max(1, m)
  
  if (!is.na(seed) & seed == "sample") {
    seed <- sample.int(1e9,1)
  }
  if (!is.na(seed)) set.seed(seed)
  
  # CHECKS
  #--------
  if (is.null(strata)) stop("Argument strata is missing.", call. = FALSE) 
  # If strata given as variable names (check if they exist) 
  # -> change into one factor with strata names 
  if (is.character(strata) & any(!duplicated(strata))) {
    varindex <- match(strata, colnames(data))
    if (any(is.na(varindex))) stop("Unrecognized variable(s) in strata: ", 
                                   paste(strata[is.na(varindex)],collapse=", "), call. = FALSE)
    else {
      dstrataNA  <- lapply(data[, strata, drop = FALSE], addNA, ifany = TRUE) # change NA to a factor level
      strata.lab <- interaction(dstrataNA, drop = TRUE, sep = "_")
      #strata.varnames <- paste0(strata, collapse="_") 
    }
  } else {
    # check length of strata vector if given as vector; check for missing
    if (length(strata) != nrow(data)) stop(paste("The length of strata index (",
                                                 length(strata), ") does not match the number of rows in the data (",
                                                 nrow(data),").",sep=""), call. = FALSE)
    if (any(is.na(strata))) stop("Strata indicator cannot have missing values.", 
                                 call. = FALSE)
    strata.lab <- factor(strata)
  }
  #--------
  
  # make sure stratification variables are included in visit.sequence
  # important when drop.not.used==T
  if (is.character(strata) & any(!duplicated(strata))){   #GR-20/10/2016 drop.not.used == TRUE removed from the condition
    strata <- match(strata, colnames(data))
    if (is.character(visit.sequence)) visit.sequence <- match(visit.sequence, colnames(data))
    if (any(is.na(visit.sequence))) stop("Unrecognized variable(s) in visit.sequence.", call. = FALSE)
    visit.sequence <- c(visit.sequence, strata[!(strata %in% visit.sequence)])
  }
  
  stratalev.lab <- levels(strata.lab) 
  strata.n.obs  <- table(strata.lab)   
  nstrata       <- dim(strata.n.obs)
  
  if (tab.strataobs == TRUE) {
    cat("Number of observations in strata (original data):")
    print(table(strata.lab, deparse.level = 0))
  }
  
  # check min number of observations in strata
  stratasize.stop <- minstratumsize         
  stratasize.warn <- 100 + 10*length(visit.sequence)
  
  smallstrata <- sum(strata.n.obs < stratasize.stop)
  if (smallstrata > 5) stop("In the original data multiple strata have fewer than the recommended\nnumber of observations. We advise that each should have at least ",
                            stratasize.stop, " observations\n('minstratumsize' which by default is 10 + 10 * no. of variables used in prediction).\n",
                            "You can override this by setting the parameter 'minstratumsize' to a lower value.\n", 
                            sep = "", call. = FALSE)
  if (smallstrata > 0) stop("In the original data some strata (", 
                            paste(stratalev.lab[strata.n.obs < stratasize.stop], collapse = ", "), 
                            ") have fewer than the recommended\nnumber of observations. We advise that each should have at least ",
                            stratasize.stop, " observations\n('minstratumsize' which by default is 10 + 10 * no. of variables used in prediction).\n",
                            "You can override this by setting the parameter 'minstratumsize' to a lower value.\n", 
                            sep = "", call. = FALSE)
  if (any(strata.n.obs < stratasize.warn) & print.flag == TRUE) {
    cat("CAUTION: In the original data some strata (", 
        paste(stratalev.lab[strata.n.obs < stratasize.warn], collapse=", "), 
        ") have limited numbers of observations.\nWe advise that there should be at least ", 
        stratasize.warn, 
        " observations (100 + 10 * no. of variables\nused in prediction).\n", sep = "")
  }
  
  synds.names <- c("call", "m", "syn", "method", "visit.sequence", 
                   "predictor.matrix", "smoothing", "event", "denom", "proper", "n", "k", 
                   "rules", "rvalues", "cont.na", "semicont", "drop.not.used", "drop.pred.only", 
                   "models", "seed", "var.lab", "val.lab", "obs.vars", "strata.syn", "strata.lab","numtocat","catgroups")
  synds <- list(setNames(vector("list",length(synds.names)),synds.names)) 
  synds <- rep(synds, m0)
  sel.names <- match(c("call", "m", "predictor.matrix", "proper", "strata.syn",
                       "strata.lab", "seed", "numtocat", "catgroups", "models"), synds.names)
  same.by.m <- c("call", "m", "method", "visit.sequence", "predictor.matrix", 
                 "smoothing", "event", "denom", "proper", "n", "rules", "rvalues", "cont.na", 
                 "semicont", "drop.not.used", "drop.pred.only",  "seed", "var.lab", 
                 "val.lab", "obs.vars", "strata.lab","numtocat","catgroups")
  same.by.m.idx <- match(same.by.m, synds.names) 
  
  syn.args <- as.list(match.call()[-1])
  
  strata.n.syn <- vector("list", m0) 
  for (j in 1:m0){ 
    synds[[j]]$strata.syn <- sort(factor(sample(stratalev.lab, k, replace = TRUE, 
                                                prob = strata.n.obs), levels = stratalev.lab, exclude = NULL))   
    synds[[j]]$strata.lab <- stratalev.lab
    strata.n.syn[[j]] <- table(synds[[j]]$strata.syn, deparse.level = 0)
    if (tab.stratasyn == TRUE) {  
      cat("\nNumber of observations in strata (synthetic data, m = ", j, "):", sep="")
      print(strata.n.syn[[j]])
    }
  }
  #Different way of printing (all syn in one table)  
  #cat("\nNumber of observations in strata (synthetic data):\n")
  #starta.n.syn.df <- do.call("rbind", strata.n.syn) 
  #rownames(starta.n.syn.df) <- paste0("m = ", 1:m)
  #print(starta.n.syn.df)
  
  for (j in 1:m0){  
    synds.ind <- vector("list", nstrata) # results by stratum  
    # exclude args that are not in syn()
    syn.args$strata <- syn.args$tab.stratasyn <- syn.args$tab.strataobs <- 
      syn.args$minstratumsize <- NULL
    syn.args$m <- 1
    syn.args$visit.sequence <- visit.sequence
    # vs <- visit.sequence
    # syn.args$visit.sequence <- c(vs[(vs %in% strata)],vs[!(vs %in% strata)])  #!GR 9/17 move strata to start of visit sequence
    
    
    for (i in 1:nstrata) {
      if (print.flag) cat("\nm = ",j,", strata = ", stratalev.lab[i],
                          "\n-----------------------------------------------------\n", sep="")
      if (!is.na(stratalev.lab[i])) {
        syn.args$data <- data[strata.lab == stratalev.lab[i],]  
      } else {
        syn.args$data <- data[is.na(as.character(strata.lab)),]
      }
      syn.args$k <- strata.n.syn[[j]][i]; names(syn.args$k) <- "strata"
      syn.args$seed <- NA 
      if (syn.args$k == 0 | m==0) syn.args$m <- 0 else syn.args$m <- 1
      synds.ind[[i]] <- do.call("syn", syn.args)
    }
    
    synds[[j]]$call <- match.call()
    synds[[j]]$m <- m
    synds[[j]]$proper <- proper
    synds[[j]]$predictor.matrix <- eval(parse(text = paste0("list(", 
                                                            paste0("synds.ind[[", 1:nstrata, "]]$predictor.matrix", collapse = ", "),")")))
    synds[[j]]$models <- eval(parse(text = paste0("list(", 
                                                  paste0("synds.ind[[", 1:nstrata, "]]$models", collapse = ", "),")")))
    synds[[j]]$seed <- seed
    synds[[j]][-sel.names] <- eval(parse(text = paste0("Map(rbind, ", 
                                                       paste0("synds.ind[[", 1:nstrata, "]][-sel.names]", collapse = ", "),")")))
    if (k > 0 & m > 0) rownames(synds[[j]]$syn) <- 1:k
  }
  
  if (m==1 | m==0) {
    synds <- synds[[1]] 
  } else {
    synds <- eval(parse(text = paste0("Map(list, ", paste0("synds[[", 1:m, "]]", 
                                                           collapse = ", "),")")))
    synds[same.by.m.idx] <- lapply(same.by.m.idx, function(x) synds[[x]][[1]])
  }
  
  if (m==0) cat("\nCAUTION: method, visit.sequence and predictor.matrix are lists that may vary by stratum and\nshould not be used to initialise syn. For initialising rerun with m = 0 without stratification.\n\n") 
  
  class(synds) <- "synds"
  return(synds) 
}

###-----utility.gen--------------------------------------------------------
utility.gen <- function(object, data, ...) UseMethod("utility.gen")


###-----utility.gen.default------------------------------------------------
utility.gen.default <- function(object, ...)
  stop("No compare method associated with class ", class(object), call. = FALSE)


###-----utility.gen.data.frame---utility.gen.list--------------------------
utility.gen.data.frame <- utility.gen.list <-
  function(object, data, 
           not.synthesised = NULL, cont.na = NULL,
           method = "cart", maxorder = 1,
           k.syn = FALSE, tree.method = "rpart",
           max.params = 400, print.stats = c("pMSE", "S_pMSE"),
           resamp.method = NULL, nperms = 50, cp = 1e-3, 
           minbucket = 5, mincriterion = 0, vars = NULL,
           aggregate = FALSE, maxit = 200, ngroups = NULL, 
           print.flag = TRUE, print.every = 10, 
           digits = 6, print.zscores = FALSE, zthresh = 1.6,
           print.ind.results = FALSE,
           print.variable.importance = FALSE, ...)
  {
    
    if (is.null(data)) stop("Requires parameter 'data' to give name of the real data.\n\n",  call. = FALSE)
    if (is.null(object)) stop("Requires parameter 'object' to give name of the synthetic data.\n\n",  call. = FALSE)
    
    if (is.list(object) & !is.data.frame(object)) m <- length(object)
    else if (is.data.frame(object)) m <- 1
    else stop("object must be a data frame or a list of data frames.\n", call. = FALSE)
    
    # sort out cont.na to make it into a complete named list
    cna <- cont.na
    cont.na <- as.list(rep(NA, length(data)))
    names(cont.na) <- names(data)
    if (!is.null(cna)) {
      if (!is.list(cna) | any(names(cna) == "") | is.null(names(cna)))
        stop("Argument 'cont.na' must be a named list with names of selected variables.", call. = FALSE)
      if (any(!names(cna) %in% names(data))) stop("Names of the list cont.na must be variables in data.\n", call. = FALSE)
      for (i in 1:length(cna)) {
        j <- (1:length(data))[names(cna)[i] == names(data)]
        cont.na[[j]] <- unique(c(NA,cna[[i]]))
      }
    }
    
    syn.method = rep("ok", length(data))
    if (!is.null(not.synthesised)) {
      if (!is.null(not.synthesised) && !all(not.synthesised %in% names(data))) stop("not.synthesised must be names of variables in data.\n", call. = FALSE)
      syn.method[names(data) %in% not.synthesised] <- ""
    }
    
    object <- list(syn = object, m = m, strata.syn = NULL, method = syn.method, cont.na = cont.na)
    class(object ) <- "synds"
    
    res <- utility.gen.synds(object = object, data = data, 
                             method = method, maxorder = maxorder, 
                             k.syn = k.syn, tree.method = tree.method,
                             max.params = max.params, print.stats = print.stats,
                             resamp.method = resamp.method, nperms = nperms, cp = cp, 
                             minbucket = minbucket, mincriterion = mincriterion, 
                             vars = vars, aggregate = aggregate, maxit = maxit, 
                             ngroups = ngroups, print.flag = print.flag, 
                             print.every = print.every, digits = digits, 
                             print.zscores = print.zscores, zthresh = zthresh, 
                             print.ind.results = print.ind.results, 
                             print.variable.importance = print.variable.importance)
    res$call <- match.call()
    return(res)
  }


###-----utility.gen--------------------------------------------------------
utility.gen.synds <- function(object, data, 
                              method = "cart", maxorder = 1,
                              k.syn = FALSE, tree.method = "rpart", 
                              max.params = 400, print.stats = c("pMSE", "S_pMSE"),
                              resamp.method = NULL, nperms = 50, cp = 1e-3,
                              minbucket = 5, mincriterion = 0, vars = NULL,
                              aggregate = FALSE, maxit = 200, ngroups = NULL,
                              print.flag = TRUE, print.every = 10,
                              digits = 6, print.zscores = FALSE, 
                              zthresh = 1.6, print.ind.results = FALSE,
                              print.variable.importance = FALSE, ...)
{
  m  <- object$m
  
  # Check input parameters
  if (is.null(method) || length(method) != 1 || is.na(match(method, c("cart", "logit"))))
    stop("Invalid 'method' type - must be either 'logit' or 'cart'.\n", call. = FALSE)
  if (is.null(print.stats) || any(is.na(match(print.stats, c("pMSE", "SPECKS", "PO50", "U", "S_pMSE", "S_SPECKS", "S_PO50", "S_U", "all")))))
    stop("Invalid 'print.stats'. Can only include 'pMSE', 'SPECKS', 'PO50', 'U', 'S_pMSE', 'S_SPECKS', 'S_PO50', 'S_U'.\nAternatively it can be set to 'all'.\n", call. = FALSE)
  if (!is.null(resamp.method) && is.na(match(resamp.method, c("perm", "pairs", "none"))))
    stop("Invalid 'resamp.method' type - must be NULL, 'perm', 'pairs' or 'none'.\n", call. = FALSE)
  if (aggregate == TRUE & method != "logit") stop("Aggregation only works for 'logit' method.\n", call. = FALSE)
  if (is.null(data)) stop("Requires parameter 'data' to give name of the real data.\n",  call. = FALSE)
  if (!inherits(object, "synds")) stop("Object must have class 'synds'.\n", call. = FALSE)
  if (k.syn & !is.null(resamp.method) && resamp.method == "pairs") stop('\nresamp.method = "pairs" will give the wrong answer when k.syn is TRUE.\n', call. = FALSE)
  if (is.null(tree.method) || length(tree.method) != 1 || is.na(match(tree.method, c("rpart", "ctree"))))
    stop("Invalid 'tree.method' - must be either 'rpart' or 'ctree'.\n", call. = FALSE)
  
  
  # Check selected variables and make observed and synthetic comparable
  if (!(is.null(vars))) {
    if (is.numeric(vars)){
      if (!(all(vars %in% 1:length(data)))) stop("Column indices of 'vars' must be in 1 to length(data).\n", call. = FALSE)
    } else if (!(all(vars %in% names(data)))) stop("Some 'vars' specified not in data.\n", call. = FALSE)
    data <- data[, vars, drop = FALSE]
    if (m == 1) {
      if (!all(vars %in% names(object$syn))) stop("Some 'vars' specified not in synthetic data.\n", call. = FALSE)
      else object$syn <- object$syn[, vars, drop = FALSE ]
    } else {
      if (!all(vars %in% names(object$syn[[1]]))) stop("Some 'vars' specified not in synthetic data.\n", call. = FALSE)
      else object$syn <- lapply(object$syn, "[", vars)
    }
  } else {
    if (m == 1) vars <- names(object$syn) else vars <- names(object$syn[[1]])
    if (!all(vars %in% names(data))) stop("Some variables in synthetic data not in original data.\n", call. = FALSE)
    else data <- data[, vars]  # make data match synthetic
  }
  
  # get cont.na and method parameters for stratified synthesis
  if (!is.null(object$strata.syn)) {
    cna <- object$cont.na[1,]
    syn.method <- object$method[1,]
  } else {
    cna <- object$cont.na
    syn.method <- object$method
  }
  
  cna <- cna[names(cna) %in% vars]
  
  for ( i in 1:length(cna)) {
    nm <- names(cna)[i]
    vals <- unique(cna[[i]][!is.na(cna[[i]])])  # get variables with cont.na other than missing
    if (length(vals) > 0){
      for (j in 1:length(vals))
        n_cna <- sum(vals[j] == data[,nm] & !is.na(data[,nm]))
      if (n_cna == 0) stop("\nValue ", vals[j], " identified as denoting a special or missing in cont.na for ",nm, " is not in data.\n",sep = "", call. = FALSE)
      else if (n_cna < 10 & print.flag) cat ("\nWarning: Only ",n_cna ," record(s) in data with value ",vals[j]," identified as denoting a missing value in cont.na for ",nm, "\n\n", sep = "")
    }
  }
  # Check whether some variables are unsynthesised
  incomplete <- FALSE
  nsynthd <- length(vars)
  unsyn.vars <- names(syn.method)[syn.method == ""]  # identify unsynthesised
  if (any(vars %in% unsyn.vars) & !is.null(unsyn.vars)) {
    notunsyn <- vars[!vars %in% unsyn.vars]  # synthesised vars
    if (!all(unsyn.vars %in% vars)) stop("Unsynthesised variables must be a subset of variables contributing to the utility measure.\n", call. = FALSE)
    if ( all(vars %in% unsyn.vars)) stop("Utility measure impossible if all in vars are unsynthesised.\n", call. = FALSE)
    incomplete <- TRUE
  }
  
  # Set default resampling according to completeness and print.stats (incl. S_SPECKS or S_PO50 or S_U)
  if (is.null(resamp.method)) {
    if ("S_SPECKS" %in% print.stats || "S_PO50" %in% print.stats || "S_U" %in% print.stats || incomplete) {
      resamp.method <- "pairs"
      cat('Resampling method set to "pairs" because S_SPECKS or S_PO50 or S_U in print.stats or incomplete = TRUE.\n') 
    } else if (method == "cart") resamp.method <- "perm"
  } else {
    if (incomplete & resamp.method == "perm")
      stop('Incomplete synthesis requires resamp.method = "pairs".\n', call. = FALSE)
    if (any(c("S_SPECKS", "S_PO50", "S_U") %in% print.stats) & resamp.method == "perm")
      stop('Stat SPECKS, PO50, and U requires resamp.method = "pairs" to get S_SPECKS, S_PO50, and S_U respectively.\n', call. = FALSE)
    if (resamp.method == "pairs" & m == 1) 
      stop('resamp.method = "pairs" needs a synthesis with m > 1, m = 10 suggested.\n', call. = FALSE)
  }
  
  # Drop any single value columns
  leneq1 <- function(x) length(table(as.numeric(x[!is.na(x)]), useNA = "ifany")) %in% (0:1)
  
  dchar <- sapply(data,is.character)
  if (any(dchar == TRUE)) for ( i in 1:dim(data)[2]) if (dchar[i] == TRUE) data[,i] <- factor(data[,i])
  dout <- sapply(data,leneq1)
  if (m == 1) sout <- sapply(object$syn,leneq1)
  else  sout <- sapply(object$syn[[1]],leneq1)
  dout <- dout & sout
  if (any(dout == TRUE) & print.flag) {
    cat("Some columns with single values or all missing values in original and synthetic\nexcluded from utility comparisons (excluded variables: ",
        paste(names(data)[dout], collapse = ", "), ").\n", sep = "")
    data <- data[,!dout]
    if (m == 1) object$syn <- object$syn[, !dout, drop = FALSE]
    else object$syn <- lapply(object$syn, "[", !dout)
  }
  
  # Numeric variables
  numvars <- (1:dim(data)[2])[sapply(data, is.numeric)]
  names(numvars) <- names(data)[numvars]
  # If ngroups != NULL divide numeric variables into ngroups
  data0 <- data  # to save if m > 1
  
  if (!is.null(ngroups)) {
    for (i in numvars) {
      if (m == 1) {
        groups <- group_num(data[,i], object$syn[,i], object$syn[,i], 
                            ngroups, cont.na = cna, ...)
        data[,i] <- groups[[1]]
        object$syn[,i] <- groups[[2]]
      } else {
        syn0 <- c(sapply(object$syn, '[[', i)) 
        for (j in 1:m) {
          groups <- group_num(data0[,i], object$syn[[j]][,i], syn0,
                              ngroups, cont.na = cna[[i]], ...)
          data[,i] <- groups[[1]]
          object$syn[[j]][,i] <- groups[[2]]
        }
      }
    }
  }
  
  # Categorical vars: make missing data part of factor
  catvars <- (1:dim(data)[2])[sapply(data, is.factor)]
  for (i in catvars) {
    data[,i] <- factor(data[,i])
    if (m == 1) object$syn[,i] <- factor(object$syn[,i])
    else for (j in 1:m) object$syn[[j]][,i] <- factor(object$syn[[j]][,i])
    if (any(is.na(data[,i]))) {
      data[,i] <- addNA(data[,i])
      if (m == 1) object$syn[,i] <- addNA(object$syn[,i])
      else for (j in 1:m) object$syn[[j]][,i] <- addNA(object$syn[[j]][,i])
    }
  }
  
  for (i in numvars) {
    if (anyNA(data[,i]) & is.null(ngroups)) {
      newname <- paste(names(data)[i], "NA", sep = "_")
      data <- data.frame(data, 1*(is.na(data[,i])))
      names(data)[length(data)] <- newname
      data[is.na(data[,i]), i] <- 0
      if (m == 1) {
        object$syn <- data.frame(object$syn, 1*(is.na(object$syn[,i])))
        names(object$syn)[length(object$syn)] <- newname
        object$syn[is.na(object$syn[,i]), i] <- 0
      } else {
        for (j in 1:m) {
          object$syn[[j]] <- data.frame(object$syn[[j]], 1*(is.na(object$syn[[j]][,i])))
          names(object$syn[[j]])[length(object$syn[[j]])] <- newname
          object$syn[[j]][is.na(object$syn[[j]][,i]),i] <- 0
        }
      }
    }
    if (any(!is.na(cna[[i]]))  & is.null(ngroups)) {
      cna[[i]] <- cna[[i]][!is.na(cna[[i]])]
      for (j in 1:length(cna[[i]])) {
        newname <- paste(names(data)[i], "cna",j, sep = "_")
        data <- data.frame(data, 1*(data[,i] == cna[[i]][j]))
        data[data[,i] == cna[[i]][j], i] <- 0
        names(data)[length(data)] <- newname
      }
      if (m == 1) {
        for (j in 1:length(cna[[i]])) {
          newname <- paste(names(object$syn)[i], "cna",j, sep = "_")
          object$syn <- data.frame(object$syn, 1*(object$syn[,i] == cna[[i]][j]))
          object$syn[object$syn[,i] == cna[[i]][j], i] <- 0
          names(object$syn)[length(object$syn)] <- newname
        }
      } else {
        for (k in 1:m) {
          for (j in 1:length(cna[[i]])) {
            newname <- paste(names(object$syn[[k]])[i], "cna",j, sep = "_")
            object$syn[[k]] <- data.frame(object$syn[[k]], 1*(object$syn[[k]][,i] == cna[[i]][j]))
            object$syn[[k]][object$syn[[k]][,i] == cna[[i]][j], i] <- 0
            names(object$syn[[k]])[length(object$syn[[k]])] <- newname
          }
        }
      }
    }
  }
  
  # Function for getting propensity scores
  # --------------------------------------
  propcalcs <- function(syndata, data) {
    
    n1 <- dim(data)[1]
    n2 <- dim(syndata)[1]
    N <- n1 + n2
    cc <- n2 / N
    if (k.syn) cc <- 0.5
    
    df.prop <- rbind(syndata, data)  # make data frame for calculating propensity score
    df.prop <- data.frame(df.prop, t = c(rep(1,n2), rep(0,n1)))
    
    # remove any levels of factors that don't exist in data or syndata
    catvars <- (1:(dim(df.prop)[2]))[sapply(df.prop,is.factor)]
    for (i in catvars) {
      if (any(table(df.prop[,i]) == 0)) {
        df.prop[,i] <- as.factor(as.character(df.prop[,i]))
        if (print.flag) cat("Empty levels of factor(s) for variable ", names(df.prop)[i]," removed.\n" )
      }
    }
    
    if (aggregate == TRUE) {
      aggdat <- aggregate(df.prop[,1], by = df.prop, FUN = length)
      wt <- aggdat$x
      aggdat <- aggdat[, -dim(aggdat)[2]]
    }
    
    if (method == "logit" ) {
      
      if (maxorder >= dim(data)[2])
        stop("maxorder cannot be greater or equal to the number of variables.\n", call. = FALSE)
      
      # cheking for large models
      levs <- sapply(data, function(x) length(levels(x)))
      levs[levs == 0] <- 2
      tt1 <- apply(combn(length(levs), 1), 2, function(x) {prod(levs[x] - 1)})
      if (maxorder == 0) nparams <- 1 + sum(tt1)
      else {
        tt2 <- apply(combn(length(levs), 2), 2, function(x) {prod(levs[x] - 1)})
        if (maxorder == 1) nparams <- 1 + sum(tt1) + sum(tt2)
        else {
          tt3 <- apply(combn(length(levs), 3), 2, function(x) {prod(levs[x] - 1)})
          if (maxorder == 2) nparams <- 1 + sum(tt1) + sum(tt2) + sum(tt3)
          else {
            tt4 <- apply(combn(length(levs), 4), 2, function(x) {prod(levs[x] - 1)})
            if (maxorder == 3) nparams <- 1 +  sum(tt1) + sum(tt2) + sum(tt3) + sum(tt4)
            else {
              tt5 <- apply(combn(length(levs), 5), 2, function(x) {prod(levs[x] - 1)})
              if (maxorder == 4) nparams <- 1 +  sum(tt1) + sum(tt2) + sum(tt3) + sum(tt4) + sum(tt5)
            }  
          }  
        }  
      }
      if (nparams > max.params) stop("You will be fitting a large model with ", nparams, 
                                     " parameters that may take a long time and fail to converge.
Have you selected variables with vars?
You can try again, if you really want to, by increasing max.params.\n", sep = "", call. = FALSE)
      else if (nparams > dim(data)[[1]]/5) cat("You will be fitting a large model with ", nparams,
                                               " parameters and only ", dim(data)[[1]], " records
that may take a long time and fail to converge.
Have you selected variables with vars?\n")
      
      if (maxorder >= 1) logit.int <- as.formula(paste("t ~ .^", maxorder + 1))
      else logit.int <- as.formula(paste("t ~ ."))
      
      if (aggregate == TRUE) fit <- glm(logit.int, data = aggdat, family = "binomial",
                                        control = list(maxit = maxit), weights = wt)
      else fit <- suppressWarnings(glm(logit.int, data = df.prop, family = "binomial",
                                       control = list(maxit = maxit)))
      #if (fit$converged == FALSE) cat("\nConvergence failed.\n")
      
      # Get number of parameters that involve synthesised variables
      score <- predict(fit, type = "response")
      if (incomplete == FALSE) km1 <- length(fit$coefficients[!is.na(fit$coefficients)]) - 1  # To allow for non-identified coefficients
      else {
        namescoef <- names(fit$coefficients)
        coefOK <- rep(FALSE, length(namescoef))
        for (nn in notunsyn) coefOK[grepl(nn, namescoef)] <- TRUE
        km1 <- sum(coefOK & print.flag)
        if (m == 1 || (m > 1 & j == 1)) cat("Expectation of utility uses only coefficients involving synthesised variables: ",
                                            km1, " from ", length(fit$coefficients) - 1, "\n", sep = "")
      }
      # one more coefficient (intercept needed if k.syn TRUE)
      if (k.syn) km1 <- km1 + 1
      if (aggregate == TRUE) {
        pMSE <- (sum(wt*(score - cc)^2, na.rm = T)) / N
        KSt <- suppressWarnings(ks.test(rep(score[aggdat$t == 1], wt[aggdat$t == 1]),
                                        rep(score[aggdat$t == 0], wt[aggdat$t == 0])))
        SPECKS <- KSt$statistic
        PO50 <- sum(wt[(score > 0.5 & df.prop$t == 1) | ( score <= 0.5 & df.prop$t == 0)])/N*100 - 50
        U      <- suppressWarnings(wilcox.test(rep(score[aggdat$t == 1], wt[aggdat$t == 1]),
                                               rep(score[aggdat$t == 0], wt[aggdat$t == 0]))$statistic) 
      } else {
        pMSE <- (sum((score - cc)^2, na.rm = T)) / N
        KSt <- suppressWarnings(ks.test(score[df.prop$t == 1], score[df.prop$t == 0]))
        SPECKS <- KSt$statistic
        PO50 <- sum((score > 0.5 & df.prop$t == 1) | ( score <= 0.5 & df.prop$t == 0))/N*100 - 50
        U      <- suppressWarnings(wilcox.test(score[df.prop$t == 1], score[df.prop$t == 0])$statistic) 
      }
      pMSEExp <- km1 * (1 - cc)^2 * cc / N
      S_pMSE  <- pMSE / pMSEExp
      
      # to save space
      fit$data <- NULL
      # fit$model <- fit$residuals <- fit$y <- NULL ?
      
    } else if (method == "cart") {
      km1 <- NA
      if (tree.method == "rpart") {
        fit <- rpart(t ~ ., data = df.prop, method = 'class',
                     control = rpart.control(cp = cp, minbucket = minbucket))
        score <- predict(fit)[, 2]
      } else if (tree.method == "ctree") {
        fit <- ctree(t ~ ., data = df.prop,
                     controls = ctree_control(mincriterion = mincriterion, minbucket = minbucket))
        score <- predict(fit)
      }
      pMSE <- sum((score - cc)^2, na.rm = T) / N
      KSt <- suppressWarnings(ks.test(score[df.prop$t == 1], score[df.prop$t == 0]))
      SPECKS <- KSt$statistic
      PO50 <- sum((score > 0.5 & df.prop$t == 1) | ( score <= 0.5 & df.prop$t == 0))/N*100 - 50
      U <- suppressWarnings(wilcox.test(score[df.prop$t == 1], score[df.prop$t == 0])$statistic)
    }
    
    # Permutations
    if (!is.null(resamp.method) && resamp.method == "none") S_pMSE <- NA
    else if (!is.null(resamp.method) && resamp.method == "perm") { # to allow resamp for logit models
      S_pMSE <- rep(NA, m)
      simpMSE <- rep(0, nperms)
      if (m == 1) j <- 1
      if (j == 1 & print.flag) {
        if (print.every == 0 | print.every >= nperms) cat("Running ", nperms, " permutations to get NULL utilities.", sep = "")
        else cat("Running ", nperms, " permutations to get NULL utilities and printing every ", print.every, "th.", sep = "")
      }
      #if (print.flag) cat("\nsynthesis ", j, "   ", sep = "")
      if (print.flag) cat("\nsynthesis ")
      
      for (i in 1:nperms) {
        if (print.every > 0 & nperms > print.every & floor(i/print.every) == i/print.every & print.flag)  cat(i, " ", sep = "")
        pdata <- df.prop
        if (!k.syn) pdata$t <- sample(pdata$t)
        else pdata$t <- rbinom(N, 1, 0.5)
        
        if (method == "cart") {
          if (tree.method == "rpart") {
            sfit <- rpart(t ~ ., data = pdata, method = 'class', control = rpart.control(cp = cp, minbucket = minbucket))
            score <- predict(sfit)[,2]
          } else if (tree.method == "ctree") {
            sfit <- ctree(t ~ ., data = pdata,
                          controls = ctree_control(mincriterion = mincriterion, minbucket = minbucket))
            score <- predict(sfit)
          }
          simpMSE[i] <- (sum((score - cc)^2, na.rm = T)) / N / 2
          
        } else if (method == "logit") {
          if (maxorder >= 1) logit.int <- as.formula(paste("t ~ .^", maxorder + 1))
          else logit.int <- as.formula(paste("t ~ ."))
          
          if (aggregate == TRUE) {
            aggdat1 <- aggregate(pdata[,1], by = pdata, FUN = length)
            wt <- aggdat1$x
            aggdat1 <- aggdat1[, -dim(aggdat1)[2]]
            sfit <- glm(logit.int, data = aggdat1, family = "binomial",
                        control = list(maxit = maxit), weights = wt)
          } else sfit <- glm(logit.int, data = pdata, family = "binomial", 
                             control = list(maxit = maxit))
          
          if (sfit$converged == FALSE & print.flag) cat("Warning: Logistic model did not converge in ",
                                                        maxit, " iterations.\nYou could try increasing parameter 'maxit'.\n", sep = "")
          score <- predict(sfit, type = "response")
          if (aggregate == TRUE) {
            simpMSE[i] <- sum(wt*(score - cc)^2, na.rm = T) / N / 2  # reduced by factor of 2
          } else {
            simpMSE[i] <- sum((score - cc)^2, na.rm = T) / N / 2  # reduced by factor of 2
          }
        }
      }
      nnosplits <- c(sum(simpMSE < 1e-8), length(simpMSE))
      S_pMSE <- pMSE/mean(simpMSE)
    }
    if (!is.null(resamp.method) && resamp.method == "pairs") 
      res.ind <- list(pMSE = pMSE, SPECKS = SPECKS, PO50 = PO50, U = U, 
                      S_pMSE= NA, S_SPECKS = NA,  S_PO50 = NA, S_U = NA, 
                      fit = fit, nnosplits = NA, df = NA)
    else if (!is.null(resamp.method) && resamp.method == "perm") 
      res.ind <- list(pMSE = pMSE, SPECKS = SPECKS, PO50 = PO50,U = U,
                      S_pMSE= S_pMSE, S_SPECKS = NA, S_PO50 = NA, S_U = NA,
                      fit = fit, nnosplits = nnosplits, df = NA)
    else res.ind <- list(pMSE = pMSE, SPECKS = SPECKS, PO50 = PO50, U =U,
                         S_pMSE = S_pMSE, S_SPECKS = NA, S_PO50 = NA, S_U = NA,
                         fit = fit, nnosplits = NA, df = km1) ## changed to NA
    return(res.ind)
  }
  # --------------------------------------
  # end propcalcs
  
  n1 <- nrow(data)
  
  if (m == 1) {
    n2 <- nrow(object$syn)
    res.ind <- propcalcs(object$syn, data)
    res <- list(call = match.call(), m = m, method = method, tree.method = tree.method,
                resamp.method = resamp.method, maxorder = maxorder, vars = vars,
                k.syn = k.syn, aggregate = aggregate, maxit = maxit,
                ngroups = ngroups, mincriterion = mincriterion,
                nperms = nperms, df = res.ind$df, incomplete = incomplete,
                pMSE = res.ind$pMSE, S_pMSE = res.ind$S_pMSE,
                S_SPECKS = res.ind$S_SPECKS, S_PO50 = res.ind$S_PO50,S_U = res.ind$S_U,
                SPECKS = res.ind$SPECKS, PO50 = res.ind$PO50, U = res.ind$U,
                print.stats = print.stats,
                fit = res.ind$fit, nnosplits = res.ind$nnosplits,
                digits = digits, print.ind.results = print.ind.results,
                print.zscores = print.zscores, zthresh = zthresh,
                print.variable.importance = print.variable.importance)
  } else {
    n2 <- nrow(object$syn[[1]])
    pMSE <- SPECKS <- PO50 <- U <- S_pMSE <- S_SPECKS <- S_PO50 <- S_U <- rep(NA, m)
    fit <- nnosplits <- as.list(1:m)
    if (!is.null(resamp.method) && !(resamp.method == "none") && resamp.method == "pairs") {
      kk <- 0
      simpMSE <- simKS <- simPO50 <- simU <- rep(NA, m*(m - 1)/2)
    }
    for (j in 1:m) {
      res.ind <- propcalcs(object$syn[[j]], data)
      pMSE[j] <- res.ind$pMSE
      SPECKS[j] <- res.ind$SPECKS
      PO50[j] <- res.ind$PO50
      U[j] <- res.ind$U
      fit[[j]] <- res.ind$fit
      
      if (resamp.method == "none" || (method == "logit" & (is.null(resamp.method)))) {
        if (j == 1 & print.flag) cat("Fitting syntheses: ")
        if (print.flag) {
          cat(j, " ", sep = "")
          if (res.ind$fit$converged == FALSE) cat("Convergence failed.\n")
        }
        if (j == m ) cat("\n")
        S_pMSE[j] <- res.ind$S_pMSE
      }
      
      if (!is.null(resamp.method) && resamp.method == "pairs") {
        if (j == 1 & print.flag) {
          if (print.every == 0 | m*(m - 1)/2 <= print.every) cat("Simulating NULL pMSE from ", m*(m - 1)/2, " pair(s).", sep = "")
          else cat("Simulating NULL pMSE from ", m*(m - 1)/2, " pairs, printing every ", print.every, "th:\n", sep = "")
          if (m*(m - 1)/2 < 6 ) cat("\nNumber of pairs too low, we suggest increasing number of syntheses (m).\n")
        }
        if (j < m) {
          for (jj in (j + 1):(m)) {
            kk <- kk + 1
            if (print.every > 0 & print.every < m*(m - 1)/2 & floor(kk/print.every) == kk/print.every & print.flag) cat(kk," ",sep = "")
            simvals <- propcalcs(object$syn[[j]], object$syn[[jj]])
            simpMSE[kk] <- simvals$pMSE
            simKS[kk] <- simvals$SPECKS
            simPO50[kk] <- simvals$SPECKS
            simU[kk] <- simvals$U
          }
        }
        nnosplits<- c(sum(simpMSE < 1e-8), length(simpMSE))
        for (j in 1:m) {
          S_pMSE[j] <- pMSE[j] *2 /mean(simpMSE)
          S_SPECKS[j] <- SPECKS[j] *2 /mean(simKS)
          S_PO50[j] <- PO50[j] *2 /mean(simPO50)
          S_U[j] <- U[j] *2 /mean(simU)
        }
        
      } else {
        nnosplits[[j]] <- res.ind$nnosplits  
        S_pMSE[j] <- res.ind$S_pMSE
        S_SPECKS[j] <- res.ind$S_SPECKS
        S_PO50[j] <- res.ind$S_PO50
        S_U[j] <- res.ind$S_U
      }
    }   
    res <- list(call = match.call(), m = m, method = method, tree.method = tree.method,
                resamp.method = resamp.method, maxorder = maxorder, vars = vars,
                k.syn = k.syn, aggregate = aggregate, maxit = maxit, 
                ngroups = ngroups, mincriterion = mincriterion,
                nperms = nperms, df = res.ind$df, incomplete = incomplete,
                pMSE = pMSE,  S_pMSE = S_pMSE,
                S_SPECKS = S_SPECKS, S_PO50 = S_PO50, S_U = S_U,
                SPECKS = SPECKS, PO50 = PO50, U = U,
                print.stats = print.stats,
                fit = fit, nnosplits = nnosplits,
                digits = digits, print.ind.results = print.ind.results,
                print.zscores = print.zscores, zthresh = zthresh,
                print.variable.importance = print.variable.importance)
    
  }
  class(res) <- "utility.gen"
  res$call <- match.call()
  return(res)
}


###-----utility.tab--------------------------------------------------------
utility.tab <- function(object, data, ...) UseMethod("utility.tab")


###-----utility.tab.default------------------------------------------------
utility.tab.default <- function(object, ...)
  stop("No compare method associated with class ", class(object), call. = FALSE)


###-----utility.tab.data.frame---utility.tab.list--------------------------
utility.tab.data.frame <- utility.tab.list <-
  function(object, data, vars = NULL, cont.na = NULL,
           ngroups = 5, useNA = TRUE, max.table = 1e6,
           print.tables = length(vars) < 4,
           print.stats = c("pMSE", "S_pMSE", "df"), 
           print.zdiff = FALSE, print.flag = TRUE,
           digits = 4, k.syn = FALSE, ...)
  {
    
    if (is.null(data)) stop("Requires parameter 'data' to give name of the real data.\n",  call. = FALSE)
    if (is.null(object)) stop("Requires parameter 'object' to give name of the synthetic data.\n",  call. = FALSE)
    
    if (is.list(object) & !is.data.frame(object)) m <- length(object)
    else if (is.data.frame(object)) m <- 1
    else stop("object must be a data frame or a list of data frames.\n", call. = FALSE)
    
    # sort out cont.na to make it into a complete named list
    cna <- cont.na
    cont.na <- as.list(rep(NA, length(data)))
    names(cont.na) <- names(data)
    if (!is.null(cna)) {
      if (!is.list(cna) | any(names(cna) == "") | is.null(names(cna)))
        stop("Argument 'cont.na' must be a named list with names of selected variables.", call. = FALSE)
      if (any(!names(cna) %in% names(data))) stop("Names of the list cont.na must be variables in data.\n", call. = FALSE)
      for (i in 1:length(cna)) {
        j <- (1:length(data))[names(cna)[i] == names(data)]
        cont.na[[j]] <- unique(c(NA,cna[[i]]))
      }
    }
    
    object <- list(syn = object, m = m, cont.na = cont.na)
    class(object ) <- "synds"
    
    res <- utility.tab.synds(object = object, data = data, vars = vars, 
                             ngroups = ngroups, useNA = useNA, 
                             print.tables = print.tables,
                             print.stats = print.stats, 
                             print.zdiff = print.zdiff,
                             print.flag = print.flag,
                             digits = digits, k.syn = k.syn, ...)
    return(res)
  }


###-----utility.tab--------------------------------------------------------
utility.tab.synds <- function(object, data, vars = NULL, ngroups = 5,
                              useNA = TRUE, max.table = 1e6, 
                              print.tables = length(vars) < 4,
                              print.stats = c("pMSE", "S_pMSE", "df"), 
                              print.zdiff = FALSE, print.flag = TRUE,
                              digits = 4, k.syn = FALSE, ...)
{
  vars <- unique(vars)
  
  # CHECKS
  #---------
  if (is.null(data))
    stop("Requires parameter 'data' to give name of the real data.\n", call. = FALSE)
  if (!is.data.frame(data))
    stop("Data must have class 'data.frame'.\n", call. = FALSE)
  if (!inherits(object, "synds"))
    stop("Object must have class 'synds'.\n", call. = FALSE)
  if (is.null(vars)) stop("You need to set variables with vars parameter.\n", call. = FALSE) else if
  (!(all(vars %in% names(data)))) stop("Unrecognized variable(s) in vars parameter: ",
                                       paste(vars[!(vars %in% names(data))], collapse = ", "), call. = FALSE)
  if (!all(print.stats %in% c("VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", "MabsDD", "dBhatt","S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE", "df", "dfG", "all")))
    stop('print.stats must be set to "all" or selected from "VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", "MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE", "df" or "dfG".\n', call. = FALSE)
  #---------
  
  data <- data[, vars, drop = FALSE]
  nvars <- ncol(data)
  data.orig <- data
  
  # get cont.na parameters for stratified synthesis
  # --------
  if (!is.null(object$strata.syn)) {
    # cna <- apply(object$cont.na, 2, function(y) {unlist(unique(y))})
    cna <- object$cont.na[1, ]
  } else {
    cna <- object$cont.na
  }
  cna <- cna[vars]
  
  m <- object$m
  if (m == 1) syndata <- list(object$syn) else syndata <- object$syn
  syndata <- lapply(syndata, '[', vars)
  
  pMSE <- S_pMSE  <- df <- dfG <- VW <- S_VW  <- FT <- S_FT  <- G  <-  S_G  <- 
    JSD  <- U <- S_JSD <- MabsDD <- WMabsDD <- S_WMabsDD <- SPECKS <-  
    dBhatt <- PO50 <- nempty <- vector("numeric", m)
  tab.syn <- tab.obs <- tab.zdiff <- vector("list", m)
  
  syn.mvar <- vector("list", nvars)
  for (j in 1:nvars) {
    if (is.numeric(syndata[[1]][, j])) syn.mvar[[j]] <- c(sapply(syndata, '[[', j))
  }  
  
  for (i in 1:m) {
    data <- data.orig
    # make all variables into factors
    for (j in 1:nvars) {
      if (is.numeric(data[, j])) {
        grpd <- group_num(data[, j], syndata[[i]][, j], syn.mvar[[j]],
                          n = ngroups, cont.na = cna[[j]], ...)
        data[, j] <- grpd[[1]]; syndata[[i]][, j] <- grpd[[2]]
      } else if (is.character(data[, j])) {
        data[, j] <- factor(data[, j])
        syndata[[i]][, j] <- factor(syndata[[i]][, j],
                                    levels = levels(data[, j]))
      }
      if (any(is.na(data[, j])) & useNA) {
        # makes missings into part of factors if present
        data[, j] <- addNA(data[, j])
        syndata[[i]][, j] <- addNA(syndata[[i]][, j])
      }
    }
    
    ## check table size
    table.size <- prod(sapply(data, function(x) length(levels(x))))
    if (table.size > max.table)
      stop("Table size ", round(table.size), " exceeds max.table limit of ", round(max.table),".",
           "\nYou could try increasing max.table but memory problems are likely.\n", call. = FALSE)
    else if (i == 1 & table.size > dim(data)[1]/2 & print.flag) cat("Warning: You are creating tables with ", table.size, 
                                                                    " cells from ", dim(data)[1], " observations.\nResults from sparse tables may be unreliable.\n", sep = "")
    ## make tables
    if (useNA){
      tab.obs[[i]] <- table(data, useNA = "ifany", deparse.level = 0)
      tab.syn[[i]] <- table(syndata[[i]], useNA = "ifany", deparse.level = 0)
    } else {
      tab.obs[[i]] <- table(data, useNA = "no", deparse.level = 0)
      tab.syn[[i]] <- table(syndata[[i]], useNA = "no", deparse.level = 0)
    }
    
    ## remove cells all zeros
    nempty[i] <-   sum(tab.obs[[i]] + tab.syn[[i]] == 0)
    td <- tab.obs[[i]][tab.obs[[i]] + tab.syn[[i]]  > 0]
    ts <- tab.syn[[i]][tab.obs[[i]] + tab.syn[[i]]  > 0]
    totcells <- length(td)
    
    ## calculate utility measures
    if (!k.syn) df[i] <- totcells - 1 else df[i] <- totcells
    cc      <- sum(ts) / sum(ts + td)
    N       <- sum(ts + td)
    sumos   <- ts + td
    expect  <- sumos * cc
    diff    <- ts - td * cc / (1 - cc)
    VW[i]   <- sum(diff^2 / expect)
    FT[i]   <- 4*sum((ts^(0.5) - (cc / (1 - cc) * td)^(0.5))^2)
    S_FT[i] <- FT[i] / df[i]
    S_VW[i] <- S_pMSE[i] <- VW[i] / df[i]
    pMSE[i] <- VW[i] * cc * (1 - cc)^2 / N
    ## standardized difference (diff/sqrt(expect))
    tab.zdiff[[i]] <- suppressWarnings((tab.syn[[i]] - tab.obs[[i]] * cc/(1-cc)) / 
                                         sqrt((tab.syn[[i]] + tab.obs[[i]]) * cc))
    ## Jensen-Shannon divergence 
    ptabd    <- td / sum(td)
    ptabs    <- ts / sum(ts)
    phalf    <- (ptabd + ptabs) *0.5
    JSD[i]   <- sum((ptabd * log2(ptabd/phalf))[ptabd > 0])/2 +        
      sum((ptabs * log2(ptabs/phalf))[ptabs > 0])/2
    S_JSD[i] <- JSD[i]*2*N/df[i]/log(2)
    ## Symmetric likelihood ratio chisq
    sok     <- ts[ts > 1e-8 & td > 1e-8]                              
    dok     <- td[ts > 1e-8 & td > 1e-8]
    if (!k.syn) dfG[i] <- length(dok) - 1 else dfG[i] <- length(dok)
    G[i]    <-  2 *sum(sok*log(sok/sum(sok)/dok*sum(dok)))
    S_G[i]  <- G[i] / dfG[i]
    ## Kolmogorov-Smirnov
    score   <- ts / (ts + td)                                                 
    kst     <- suppressWarnings(ks.test(rep(score, ts), rep(score, td)))
    SPECKS[i] <- kst$statistic
    ## Wilcoxon statistic 
    Ut      <- suppressWarnings(wilcox.test(rep(score, ts), rep(score, td))) 
    U[i]    <- Ut$statistic
    ## Calculate PO50
    predsyn <- (ptabs > ptabd)                                        
    PO50[i] <- (sum(ts[predsyn]) + sum(td[!predsyn])) / (sum(ts) + sum(td)) * 100 - 50
    
    MabsDD[i]  <- sum(abs(diff))/sum(ts)
    WMabsDD[i] <- sum(abs(diff)/sqrt(expect))/sqrt(2/pi)
    S_WMabsDD[i] <- WMabsDD[i]/df[i]
    
    dBhatt[i] <- sqrt(1 - sum(sqrt(ptabd*ptabs)))
  }
  
  tab.obs <- tab.obs[[1]]  
  
  if (m == 1) {
    tab.syn <- tab.syn[[1]]
    tab.zdiff <- tab.zdiff[[1]]
  }
  
  res <- list(m = m,
              VW = VW, 
              FT = FT,
              JSD = JSD, 
              SPECKS = SPECKS,
              WMabsDD = WMabsDD,
              U = U,
              G = G,
              pMSE = pMSE, 
              PO50 = PO50,
              MabsDD = MabsDD,
              dBhatt = dBhatt, 
              S_VW = S_VW,
              S_FT = S_FT,
              S_JSD = S_JSD,
              S_WMabsDD = S_WMabsDD,
              S_G = S_G,
              S_pMSE = S_pMSE,
              df = df,
              dfG = dfG,
              nempty = unlist(nempty),
              tab.obs = tab.obs,
              tab.syn = tab.syn,
              tab.zdiff = tab.zdiff,
              digits = digits,
              print.stats = print.stats,
              print.zdiff = print.zdiff,
              print.tables = print.tables,
              n = sum(object$n),
              k.syn = k.syn)
  
  class(res) <- "utility.tab"
  return(res)
}


###-----group_num----------------------------------------------------------
# function to categorise continuous variables

group_num <- function(x1, x2, xsyn, n = 5, style = "quantile", cont.na = NA, ...) {
  
  # Categorise 2 continuous variables into factors of n groups
  # with same groupings determined by the first one
  # xsyn - all synthetic values (for m syntheses)
  
  if (!is.numeric(x1) | !is.numeric(x2) | !is.numeric(xsyn)) 
    stop("x1, x2, and xsyn must be numeric.\n", call. = FALSE)
  
  # Select non-missing(nm) values
  x1nm <- x1[!(x1 %in% cont.na) & !is.na(x1)]
  x2nm <- x2[!(x2 %in% cont.na) & !is.na(x2)]
  xsynnm <- xsyn[!(xsyn %in% cont.na) & !is.na(xsyn)]
  
  # Derive breaks
  my_breaks <- unique(suppressWarnings(classIntervals(c(x1nm, xsynnm),
                                                      n = n, style = style, ...))$brks)
  
  my_levels <- c(levels(cut(x1nm, breaks = my_breaks,
                            dig.lab = 8, right = FALSE, include.lowest = TRUE)),
                 cont.na[!is.na(cont.na)])
  
  # Apply groupings to non-missing data
  x1[!(x1 %in% cont.na) & !is.na(x1)] <- as.character(cut(x1nm,
                                                          breaks = my_breaks, dig.lab = 8, right = FALSE, include.lowest = TRUE))
  x2[!(x2 %in% cont.na) & !is.na(x2)] <- as.character(cut(x2nm,
                                                          breaks = my_breaks, dig.lab = 8, right = FALSE, include.lowest = TRUE))
  x1 <- factor(x1, levels = my_levels)
  x2 <- factor(x2, levels = my_levels)
  
  return(list(x1,x2))
}

###-----utility.tables-----------------------------------------------------
utility.tables <- function(object, data, ...) UseMethod("utility.tables")


###-----utility.tables.default---------------------------------------------
utility.tables.default <- function(object, ...)
  stop("No compare method associated with class ", class(object), call. = FALSE)


###-----utility.tables.data.frame---utility.tables.list--------------------
utility.tables.data.frame <- utility.tables.list <- 
  function(object, data, 
           cont.na = NULL, not.synthesised = NULL, 
           tables = "twoway", maxtables = 5e4, 
           vars = NULL, third.var = NULL,
           useNA = TRUE, ngroups = 5,
           tab.stats = c("pMSE", "S_pMSE", "df"), 
           plot.stat = "S_pMSE", plot = TRUE,
           print.tabs = FALSE, digits.tabs = 4,
           max.scale = NULL, min.scale = 0, plot.title = NULL,
           nworst = 5, ntabstoprint = 0, k.syn = FALSE,
           low = "grey92", high = "#E41A1C",
           n.breaks = NULL, breaks = NULL, ...){
    if (is.null(data)) stop("Requires parameter 'data' to give name of the real data.\n", call. = FALSE)
    if (is.null(object)) stop("Requires parameter 'object' to give name of the synthetic data.\n", call. = FALSE)   
    
    if (is.list(object) & !is.data.frame(object)) m <- length(object)
    else if (is.data.frame(object)) m <- 1
    else stop("object must be a data frame or a list of data frames.\n", call. = FALSE)
    
    # sort out cont.na to make it into a complete named list
    cna <- cont.na
    cont.na <- as.list(rep(NA, length(data)))
    names(cont.na) <- names(data)
    if (!is.null(cna)) {
      if (!is.list(cna) | any(names(cna) == "") | is.null(names(cna))) 
        stop("Argument 'cont.na' must be a named list with names of selected variables.", call. = FALSE)  
      if (any(!names(cna) %in% names(data))) stop("Names of the list cont.na must be variables in data.\n", call. = FALSE)
      for (i in 1:length(cna)) {
        j <- (1:length(data))[names(cna)[i] == names(data)]
        cont.na[[j]] <- unique(c(NA,cna[[i]]))
      }
    }
    
    syn.method = rep("ok", length(data))
    if  (!is.null(not.synthesised)) {
      if ( !is.null(not.synthesised) && !all(not.synthesised %in% names(data))) stop("not.synthesised must be names of variables in data\n\n", call. = FALSE)
      syn.method[names(data) %in% not.synthesised] <- ""
    }
    
    object <- list(syn = object, m = m, strata.syn = NULL, 
                   method = syn.method, cont.na = cont.na)
    class(object ) <- "synds"
    
    res <- utility.tables.synds(object = object, data = data, 
                                tables = tables, maxtables = maxtables, 
                                vars = vars, third.var = third.var, 
                                useNA = useNA, ngroups = ngroups, 
                                tab.stats = tab.stats, plot.stat = plot.stat, 
                                plot = plot, print.tabs = print.tabs, 
                                digits.tabs = digits.tabs, max.scale = max.scale, 
                                min.scale = min.scale, plot.title = plot.title, 
                                nworst = nworst, ntabstoprint = ntabstoprint, 
                                k.syn = k.syn, low = low, high = high,
                                n.breaks = n.breaks, breaks = breaks, ...)
    
    res$call <- match.call()
    return(res)
  }


###-----utility.tables-----------------------------------------------------
utility.tables.synds <- function(object, data, 
                                 tables = "twoway", maxtables = 5e4, 
                                 vars = NULL, third.var = NULL,
                                 useNA = TRUE, ngroups = 5,  
                                 tab.stats = c("pMSE", "S_pMSE", "df"), 
                                 plot.stat = "S_pMSE", plot = TRUE,
                                 print.tabs = FALSE, digits.tabs = 4,
                                 max.scale = NULL, min.scale = 0, plot.title = NULL,
                                 nworst = 5, ntabstoprint = 0, k.syn = FALSE, 
                                 low = "grey92", high = "#E41A1C",
                                 n.breaks = NULL, breaks = NULL, ...){
  
  if (is.null(object)) stop("Requires parameter 'object' to give name of the synthetic data object.\n", call. = FALSE)   
  if (is.null(data)) stop("Requires parameter 'data' to give name of the original data.\n", call. = FALSE)
  if (!inherits(object, "synds")) stop("'object' must be of class 'synds', a synthetic data object", call. = FALSE)
  if (!is.data.frame(data)) stop("'data' must be a data.frame.\n", call. = FALSE)
  if (!(tables %in% c("twoway", "oneway", "threeway"))) stop("Argument tables must be 'oneway', 'twoway' or 'threeway.'\n", call. = FALSE)
  
  if (is.null(vars) ) {
    if (object$m == 1) vars <- names(object$syn)
    else vars <- names(object$syn[[1]])
    vno <- 1:length(vars)
  } else if (is.character(vars)) {
    if (!all(vars %in% names(data))) stop("vars must be in names of original data.\n", call. = FALSE)
    if (object$m == 1){
      if (!all(vars %in% names(object$syn))) stop("vars must be in names of synthetic data.\n", call. = FALSE)
      vno <- match(vars, names(object$syn))
    } else if (object$m > 1){
      if (!all(vars %in% names(object$syn[[1]]))) stop("vars must be in names of synthetic data.\n", call. = FALSE)
      vno <- match(vars, names(object$syn[[1]]))
    }
  } else if (is.numeric(vars)) {
    vno <- vars
    if (object$m == 1)  {
      if(!all(vars %in% 1:length(object$syn))) stop("vars must be in 1:length(object$syn).\n", call. = FALSE)
      vars <- names(object$syn)[vno]
    } else if (object$m > 1) {  
      if(!all(vars %in% 1:length(object$syn[[1]]))) stop("vars must be in 1:length(object$syn[[1]]).\n", call. = FALSE)
      vars <- names(object$syn[[1]])[vno]
    }
  }
  ord <- order(vno)
  vars <- vars[ord]; vno <- vno[ord]  ## get in right order in syn
  # now vars is character, vno is numeric, both refer to position in synthetic data
  
  # check tab.stats (for tabs) and plot.stat (for plots)
  if (is.null(tab.stats)) tab.stats <- c("pMSE", "S_pMSE", "df")
  else if(!all(tab.stats %in%  c("VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", 
                                 "MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE", "df", "dfG", "all"))) stop('Parameter tab.stats can only include:
"VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", 
"MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE", "df", "dfG", "all".\n', call. = FALSE, sep = "")
  if (any(tab.stats == "all")) tab.stats <- c("VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", 
                                              "MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE", "df", "dfG")
  
  if (is.null(plot.stat)) plot.stat <- "S_pMSE"
  else if (!(length(plot.stat) == 1 && plot.stat %in%  c("VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", 
                                                         "MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE"))) 
    stop('Parameter plot.stat must be just one of:\n"VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", 
"MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE".\n', call. = FALSE)
  # End of check of input parameters 
  
  # Add labels to names.vars to get plots in right order
  nv <- length(vars)
  names.vars <- paste(ntoc(vno), vars, sep = ".")
  
  # Code to make list of table utility values table of results and values to plot
  if (tables %in% c("twoway")) {
    npairs    <- nv*(nv - 1)/2
    pairs.num <- combn(1:nv, 2)
    pairs <- combn(vars, 2)
  } else if (tables == "threeway") {
    npairs <- nv*(nv - 1)*(nv - 2)/6
    pairs  <- combn(vars, 3)
    pairs.num <- combn(1:nv, 3)
    pair.names <- apply(pairs, 2, function(x) paste(as.vector(x), collapse=":"))
  } else if (tables == "oneway") {
    utility.list <- as.list(1:nv)
    npairs <- nv
    pairs <- matrix(vars, 1, nv)
    pairs.num <- matrix(1:nv, 1, nv)
  }
  
  if (npairs > maxtables) {
    sel <- sample(1:npairs, maxtables)
    pairs <- pairs[, sel, drop = FALSE]
    pairs.num <- pairs.num[, sel, drop = FALSE]
    npairs <- maxtables
    cat("Total tables requested exceeds limit of maxtables = ", maxtables,
        ", so utility is measured for a sample of ", maxtables, ".\n", sep = "")
    if (tables == "threeway") {
      plot <- FALSE    
      cat("You cannot produce plots of threeway tables from sampled tables.\n", sep = "")
    }
  }
  
  X1 <- X2 <- X3 <- rep("", npairs)
  for (i in 1:npairs){
    X1[i] <- names.vars[pairs.num[1, i]]
    if (!tables == "oneway")  X2[i] <- names.vars[pairs.num[2, i]]
    if (tables == "threeway") X3[i] <- names.vars[pairs.num[3, i]]
  }
  
  utility.list <- as.list(1:npairs)
  # now make tabs of results
  tabs <- matrix(NA, npairs, length(tab.stats)) 
  colnames(tabs) <- tab.stats
  if (tables == "oneway") dimnames(tabs)[[1]] <- paste(X1)
  else if (tables == "twoway") dimnames(tabs)[[1]] <- paste(X1, X2, sep = ":")
  else if (tables == "threeway") dimnames(tabs)[[1]] <- paste(X1, X2, X3, sep = ":")
  
  val <- rep(NA, npairs)
  
  for (i in 1:npairs) {
    utility.list[[i]] <- utility.tab(object, data, vars = pairs[, i], 
                                     ngroups = ngroups, k.syn = k.syn, 
                                     useNA = useNA, print.flag = FALSE, ...)
    if (i == 1) {
      tab.ind <- match(tab.stats, names(utility.list[[i]]))
      val.ind <- match(plot.stat, names(utility.list[[i]]) )
    }
    tabs[i,] <- sapply(utility.list[[i]][tab.ind], mean)
    val[i]      <- sapply(utility.list[[i]][val.ind], mean) # value for plotting
  }  
  
  # find worst set of variables
  nworst <- min(npairs, nworst)
  worstn <- val[order(-val)][1:nworst]
  names(worstn) <- dimnames(tabs)[[1]][order(-val)][1:nworst]
  worstnind <- (1:npairs)[order(-val)][1:nworst]
  worsttabs <- as.list(1:nworst)
  for (i in 1:nworst) {
    worsttabs[[i]] <- list(tab.obs = utility.list[[tab.obs = worstnind[i]]]$tab.obs, 
                           tab.syn = utility.list[[tab.obs = worstnind[i]]]$tab.syn)
  }
  names(worsttabs) <- names(worstn)
  
  # calculate vars with highest scores 
  var.scores <- NULL
  nn <- dimnames(tabs)[[1]]
  num <- score <- rep(NA, nv)
  names(num) <- names(score) <- names.vars
  for (i in 1:nv) {  
    num[i]   <- sum(grepl(names.vars[i], nn))
    score[i] <- sum(val[grepl(names.vars[i], nn)])
  }
  var.scores <- sort(score/num, decreasing = TRUE)
  if (tables == "threeway"){
    if (is.null(third.var)) third.var <- names(var.scores)[1]  # select worst as third.var if not specified
    else third.var <- names.vars[vars == third.var]
  } 
  
  if (tables == "oneway") toplot <- data.frame(X1 = X1, X2 = X2, val = val)
  else if (tables == "twoway") toplot <- rbind(data.frame(X1 = X1, X2 = X2, val = val), 
                                               data.frame(X1 = X2, X2 = X1, val = val))
  else if (tables == "threeway") {
    toplot <- data.frame(X1 = X1, X2 = X2, X3 = X3, val = val)
    toplot <- toplot[(third.var == toplot$X1 | third.var == toplot$X2 | third.var == toplot$X3), ]
    toplot[third.var == toplot$X1, 1:2] <- toplot[third.var == toplot$X1, 2:3]
    toplot[third.var == toplot$X2,   2] <- toplot[third.var == toplot$X2,   3]
    v2 <- toplot[, -3]
    v2$X1 <- toplot$X2
    v2$X2 <- toplot$X1
    toplot <- rbind(toplot[, c(1, 2, 4)], v2)
  }
  
  if (is.null(plot.title)){  
    if (tables == "twoway") plot.title <- bquote("Two-way utility:"~bold(.(plot.stat))~"for pairs of variables")
    else if (tables == "oneway") plot.title <- bquote("One-way utility:"~bold(.(plot.stat))~"for each variable")
    else if (tables == "threeway") plot.title <- bquote("Three-way utility:"~bold(.(plot.stat))~"for pairs along with"~bold(.(third.var)))
  }
  
  if (!is.null(max.scale)) {
    if (max(toplot$val, na.rm = TRUE) >  max.scale) {
      cat("Maximum of plot scale set to ", max.scale, 
          " (lower than maximum in results ",
          max(toplot$val, na.rm = TRUE), ").\n", sep = "")
      toplot$val[toplot$val > max.scale] <- max.scale
    }
  } else max.scale <- max(toplot$val, na.rm = TRUE) 
  
  if (!is.null(min.scale) & min.scale != 0 & (min(toplot$val, na.rm = TRUE) < min.scale)) {
    cat("Minimum of plot scale set to ", min.scale, " (greater than minimum in results ",
        min(toplot$val, na.rm = TRUE), ").\n", sep = "")
    toplot$val[toplot$val < min.scale] <- min.scale
  }
  
  if (!is.null(n.breaks)){
    plot.scale <- scale_fill_steps(
      n.breaks = n.breaks,
      low = low, high = high,
      limits = c(0, max.scale), ...)
  } else if (!is.null(breaks)){
    plot.scale <- scale_fill_steps(
      breaks = breaks,
      low = low, high = high,
      limits = c(0, max.scale), ...)
  } else {
    plot.scale <- scale_fill_gradient(
      low = low, high = high, 
      limits = c(0, max.scale), ...) 
  }
  
  
  utility.plot <- ggplot(toplot, aes(x = X2, y = X1)) + 
    geom_raster(aes(fill = val)) + 
    plot.scale + 
    labs(x = "", y = "", title = plot.title) +
    theme_minimal() + 
    theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 0.9, vjust = 0.2), 
          axis.text.y = element_text(size = 10, margin = margin(r = 0)),
          title = element_text(size = 11),
          legend.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  res <- list(tabs = tabs,
              plot.stat = plot.stat, 
              tables = tables, 
              third.var = third.var, 
              utility.plot = utility.plot,  
              var.scores = var.scores, plot = plot,
              print.tabs = print.tabs, 
              digits.tabs = digits.tabs, 
              plot.title = plot.title, 
              max.scale = max.scale, min.scale = min.scale,
              ntabstoprint = ntabstoprint, nworst = nworst, 
              worstn = worstn, worsttabs = worsttabs)
  
  class(res) <- "utility.tables"
  return(res)
}


###-----ntoc---------------------------------------------------------------
# to make labels for variables of constant length from an integer
ntoc <- function(x)
{
  nch <- nchar(max(x))
  res <- as.character(x)
  res <- str_pad(res, nch, "left", "0")
  return(res)
}


###-----print.utility.tables-----------------------------------------------
print.utility.tables <- function(x, print.tabs = NULL, digits.tabs = NULL,  
                                 plot = NULL, plot.title = NULL, 
                                 max.scale = NULL, min.scale = NULL,
                                 nworst = NULL, ntabstoprint = NULL, ...) {
  
  if (is.null(print.tabs)) print.tabs <- x$print.tabs
  if (is.null(digits.tabs)) digits.tabs <- x$digits.tabs
  if (is.null(plot)) plot <- x$plot
  if (is.null(plot.title)) plot.title <- x$plot.title 
  if (is.null(ntabstoprint)) ntabstoprint <- x$ntabstoprint
  if (is.null(nworst)) nworst <- x$nworst 
  if (is.null(max.scale)) max.scale <- x$max.scale
  if (is.null(min.scale)) min.scale <- x$min.scale
  
  if (x$tables == "twoway") cat("\nTwo-way utility: ", x$plot.stat, " value plotted for ", 
                                dim(x$tabs)[[1]], " pairs of variables.\n", sep = "")
  if (x$tables == "oneway") cat("\nUnivariate utility: ", x$plot.stat, " value plotted for ", 
                                dim(x$tabs)[[1]], " variables.\n", sep = "")
  if (x$tables == "threeway") {
    cat("\nThree-way utility (total of ", dim(x$tabs)[[1]]," variable combinations):\n", sep = "")
    cat("\nAverage of 3-way scores ", x$plot.stat,
        " (ordered) for 3-way tables including each variable.\n", sep = "")
    print(x$var.scores)
    cat("\nVariable with highest average score, ", x$third.var,
        ", selected to make plots.\nTo see others, set parameter 'third.var'.\n", sep = "")
  }
  
  cat("\nVariable combinations with worst ", nworst ,
      " utility scores (", x$plot.stat,"):\n", sep = "")
  print(round(x$worstn, digits.tabs))
  
  if (ntabstoprint > nworst) {
    ntabstoprint <- nworst
    cat("Only ", nworst, 
        " tables can be printed. For more rerun with 'nworst' set to a larger number (up to ", 
        dim(x$tabs)[[1]], ").\n", sep = "")
  }  
  if (ntabstoprint > 0) {
    cat("\nPrinting table of observed and synthetic data for the ", 
        ntabstoprint, " table with the worst utility\n", sep = "")
    for (i in 1:ntabstoprint) {
      cat("Tables of ", names(x$worsttabs[i]), "\nOriginal data:\n", sep = "")
      print(x$worsttabs[[i]]$tab.obs)
      cat("Synthetic data:\n")
      print(x$worsttabs[[i]]$tab.syn)
    }
  }
  
  if (plot) {
    if (!is.null(max.scale)) x$utility.plot$scales$scales[[1]]$limits[2] <- max.scale
    if (!is.null(min.scale)) x$utility.plot$scales$scales[[1]]$limits[1] <- min.scale
    print(x$utility.plot)
  }
  
  if (!print.tabs) {
    cat("\nMedians and maxima of selected utility measures for all tables compared\n")
    Medians <- apply(x$tabs, 2, function(x) median(x, na.rm = TRUE))
    Maxima  <- apply(x$tabs, 2, function(x) max(x, na.rm = TRUE))
    print(round(data.frame(Medians, Maxima), digits.tabs))
    cat("\nFor more details of all scores use print.tabs = TRUE.\n")
  } else {
    cat("\nTable of selected utility measures\n")
    print(round(x$tabs, digits.tabs))
  }
  
  invisible(x)
}

.onAttach <- function(...) {
  info <- "Find out more at https://www.synthpop.org.uk/"
  packageStartupMessage(paste(strwrap(info), collapse = "\n"))
}



###-----compare.CI---------------------------------------------------------
compare.CI <- function(synthetic, observed, ci.level, intercept, ...)
{
  CI <- qnorm(1- (1 - ci.level)/2)
  ##Initiate
  if(nrow(observed) > nrow(synthetic)){
    numVar <- nrow(synthetic); rNames <- rownames(synthetic)
  } else{
    numVar <- nrow(observed); rNames <- rownames(observed)
  } 
  CIoverlap <- matrix(NA, nrow = numVar, ncol = 1, 
                      dimnames = list(rNames, "CI overlap"))
  
  ##Calculate CIoverlap
  for(i in 1:numVar) {
    
    ##Store CIs
    syn.upper <- synthetic[i, 1] + (CI * synthetic[i, 2])
    syn.lower <- synthetic[i, 1] - (CI * synthetic[i, 2])
    obs.upper <- observed[i, 1]  + (CI * observed[i, 2])
    obs.lower <- observed[i, 1]  - (CI * observed[i, 2])
    
    ## CI overlap
    overlap.lower <- max(obs.lower, syn.lower)       
    overlap.upper <- min(obs.upper, syn.upper)       
    
    CIoverlap[i, 1] <- 0.5 * 
      (((overlap.upper - overlap.lower) / (obs.upper - obs.lower)) + 
         ((overlap.upper - overlap.lower) / (syn.upper - syn.lower)))
  }
  
  if(intercept == FALSE) {
    CIoverlap = CIoverlap[-1, , drop = FALSE]
  }
  
  return(as.data.frame(CIoverlap))
}




codebook.syn <- function(data, maxlevs = 3)
{
  # function to decribe variables in a data frame 
  
  if (!(is.data.frame(data))) stop("codebook.syn() requires a data frame as a parameter.\n", call. = FALSE)
  n <- dim(data)[[1]]
  p <- dim(data)[[2]]
  
  # calculate number and % of missing and non-missing values
  nmiss <- sapply(data, function(x) length(x[is.na(x)]))
  perctmiss <- round(nmiss/n * 100, 2)
  nok <- sapply(data, function(x) length(x[!is.na(x)]))
  ndistinct <- sapply(data, function(x) length(table(x)))
  dfclass <- sapply(data, class)
  
  fortab2 <- details <- rep("", length(nmiss))
  
  for (i in 1:p) {
    if (any(class(data[,i]) == "character")) details[i] <- paste("Max length: ", 
                                                                 max(nchar(data[,i])), sep = "")
    if (any(class(data[,i]) == "numeric")) details[i] <- paste("Range: ", 
                                                               min(data[!is.na(data[,i]), i]), " - ", max(data[!is.na(data[,i]),i]), 
                                                               sep = "")
    if (any(class(data[,i]) == "factor") & ndistinct[i] > maxlevs ) { 
      details[i] <- "See table in labs"
      fortab2[i] <- paste("'", paste(names(table(data[,i])), collapse = "' '"), 
                          "'", sep = "")
    }
    if (any(class(data[,i]) == "factor") & ndistinct[i] <= maxlevs ) details[i] <-
        paste("'", paste(names(table(data[,i])), collapse = "' '"),"'", sep = "")
  }
  
  if (any(grepl("factor", sapply(data, class)) & ndistinct > maxlevs )) {
    vnum <- (1:p)[grepl("factor", sapply(data, class)) & ndistinct > maxlevs]
    tabs2 <- vector("list", length(vnum))
    names(tabs2) <- names(data)[vnum]
    for (i in 1:length(vnum)) {
      tabs2[[i]] <- data.frame(label = names(table(data[, vnum[i]])))
    }
  } else tabs2 <- NULL
  
  result <- data.frame(variable = names(data), class = sapply(dfclass, paste, collapse = "-"), 
                       nmiss = nmiss, perctmiss = perctmiss, 
                       ndistinct = ndistinct, details = details)
  rownames(result) <- 1:p
  list(tab = result, labs = tabs2)
}


extended_syn <- function(data, method = "cart",
                visit.sequence = (1:ncol(data)),
                predictor.matrix = NULL,  
                m = 1, k = nrow(data), proper = FALSE, 
                minnumlevels = 1, maxfaclevels = 60,
                rules = NULL, rvalues = NULL,                                  
                cont.na = NULL, semicont = NULL,
                smoothing = NULL, event = NULL, denom = NULL,
                drop.not.used = FALSE, drop.pred.only = FALSE,                                            
                default.method = c("normrank", "logreg", "polyreg", "polr"),
                numtocat = NULL, catgroups = rep(5, length(numtocat)),
                models = FALSE,
                print.flag = TRUE,
                seed = "sample",
                ...)
{
  
  #----------------------------------------------------------------------
  # the code for the following checking functions is included within syn
  # so as to allow them to access the local objects
  #----------------------------------------------------------------------
  obs.vars <- names(data)
  
  #if (k == 0) m <- 0 
  
  # set default method for everything to cart and to blank (which will get defaults) if method is "parametric"
  # if (all(method == "")) method = "cart"  # change to "ctree"?
  # else if (length(method)==1 && method=="parametric") method=rep("",dim(data)[2])
  
  if (!is.null(attr(data,"filetype")) && attr(data,"filetype") == "sav") {
    var.lab <- attr(data,"variable.labels")
    val.lab <- attr(data,"label.table")
  } else {
    var.lab <- val.lab <- NULL
  }
  
  # check problematic characters in varaiable names
  has_space  <- grepl(" ", obs.vars) + 
    sapply(data, function(x) {res <- any(is.na(x))}) +
    sapply(data, is.numeric)
  if (any(has_space == 3)) stop(paste("Your data have numeric variable(s) with missing values with names that include spaces:\n  ",
                                      paste0("`", paste(obs.vars[has_space == 3], collapse = "`, `"), "`"),
                                      "\nThese should be renamed for synthpop to work correctly."), call. = FALSE)
  
  bad_char <- "[^\\w_.]"
  has_bad_char  <- str_detect(obs.vars, bad_char)
  if (any(has_bad_char)) cat("WARNING: Some variable names include special characters",
                             unlist(str_extract_all(obs.vars, bad_char)),"\n  ",
                             paste0(str_subset(obs.vars, bad_char), collapse = ", "),
                             "\nYou must rename them for synthpop to work correctly.")
  
  # if visit sequence includes variable names change them into column indecies 
  if (is.character(visit.sequence)) {
    nametoind <- match(visit.sequence, colnames(data))
    if (any(is.na(nametoind))) stop("Unrecognized variable(s) in visit.sequence: ", 
                                    paste(visit.sequence[is.na(nametoind)], collapse = ", "), call. = FALSE)
    else visit.sequence <- nametoind
  } else if (!all(visit.sequence %in% 1:length(data))) stop("Column indices in visit.sequence must be between 1 and ", 
                                                            length(data), sep = "", call. = FALSE)  
  
  # expand user's smoothing method (single string) to all numeric variables
  if (length(smoothing) == 1 & is.character(smoothing)) {
    numeric.vars <- which(sapply(data, is.numeric))
    smoothing <- as.list(rep(smoothing,length(numeric.vars)))
    names(smoothing) <- names(numeric.vars)
  }
  
  size.warn <- 100 + 10*length(visit.sequence)
  if (dim(data)[1] <  size.warn & print.flag == TRUE) {
    cat("CAUTION: Your data set has fewer observations (",  dim(data)[1], 
        ") than we advise.\nWe suggest that there should be at least ", 
        size.warn, " observations\n(100 + 10 * no. of variables used in modelling the data).\n",
        "Please check your synthetic data carefully with functions\ncompare(), utility.tab(), and utility.gen().\n\n", sep = "")
  }
  
  ##-----------------------code for numtocat------------------------------------------
  if (!is.null(numtocat)) {
    
    # if numtocat is numeric change to column names and check names if not
    if (is.numeric(numtocat)) {
      if (!all(numtocat %in% 1:length(data))) stop("Column indices in visit.sequence must be between 1 and ", 
                                                   length(data), sep = "", call. = FALSE)  
      numtocat <- names(data)[numtocat]
    } else if (!all(numtocat %in% names(data))) stop("numtocat variable(s): ",
                                                     paste(numtocat[!numtocat %in% names(data)], collapse = ", ")," are not in data.", sep = "", call. = FALSE)
    
    # check if numtocat in visit sequence and remove from numtocat if not
    if (!all(numtocat %in% names(data)[visit.sequence])) {
      cat("\nVariable(s) in numtocat (", paste(numtocat[!numtocat %in% names(data)[visit.sequence]], collapse = ", "),
          ") not in visit.sequence and have been removed from numtocat.\n" , sep = "")
      if (length(catgroups) > 1) catgroups <- catgroups[numtocat %in% names(data)[visit.sequence]]
      numtocat <- numtocat[numtocat %in% names(data)[visit.sequence]]
      if (length(numtocat) == 0) stop("\nAll variables in numtocat removed, perhaps try without this parameter.\n", call. = FALSE)
    }
    
    # get cont.na for numtocat vars
    cna <- cont.na
    if (!any(names(cna) %in% numtocat)) cna <- NULL
    if (!is.null(cna) & any(names(cna) %in% numtocat)) {
      cna <- cna[names(cna) %in%  numtocat]
      if (length(cna) == 0) cna <- NULL
    }
    
    # make numtocat object incorporating some checks
    numtocat.obj <- numtocat.syn(data, numtocat = numtocat, cont.na = cna, 
                                 catgroups = catgroups, print.flag = FALSE)
    
    # adjust visit.sequence
    visit.sequence <- c(visit.sequence, length(data) + 1:length(numtocat))
    
    # replace data with categorised with orig vars at end
    data <- cbind(numtocat.obj$data, numtocat.obj$orig)
    
    if (length(method) == 1) {
      if (method == "parametric") stop("'parametric' not available, when numtocat is used, methods must be specified in detail.\n", call. = FALSE )
      else method <- rep(method, length(data) - length(numtocat))
    }
    
    # assign method for final columns as nested and give changed parameters correct names
    method <- c(method, paste("nested", numtocat, sep = "."))
    # checks on methods for numtocat variables
    if (any(method[numtocat.obj$ind] %in% 
            c("norm", "normrank", "lognorm", "sqrtnorm", "cubertnorm"))) {
      nummth <- method[numtocat.obj$ind] %in% 
        c("norm", "normrank", "lognorm", "sqrtnorm", "cubertnorm")
      stop("Method(s) (", paste(method[numtocat.obj$ind][nummth], collapse = ", "),
           ") assigned to variable(s) in numcat (", 
           paste(numtocat[nummth], collapse = ", "),
           ")\nunsuitable for categorical data.", sep = "", call. = FALSE)
    }  
    
    # modify visit sequence and predictor matrix to synthesis numtocat 
    # original variables after the others
    # and adjust other parameters to match if not null
    if (!is.null(predictor.matrix)) {
      predictor.matrix <- cbind(predictor.matrix, matrix(0, dim(data)[2], length(numtocat)))
      predictor.matrix <- rbind(predictor.matrix, matrix(0, length(numtocat), dim(data)[2] + length(numtocat)))
      for (i in 1:length(numtocat)) { 
        predictor.matrix[dim(data)[2] + i, numtocat.obj$ind[i]] <- 1
      }
    }
    
    # move these parameters to numeric versions
    newnames <- function(x) { 
      # mini function to change names of lists
      xn <- names(x)
      ind <- match(numtocat, xn)
      ind <- ind[!is.na(ind)]
      names(x)[ind] <- paste("orig", names(x)[ind], sep = ".")
      return(x)
    }
    if (!is.null(rules))     rules     <- newnames(rules)
    if (!is.null(rvalues))   rvalues   <- newnames(rvalues)
    if (!is.null(cont.na))   cont.na   <- newnames(cont.na)
    if (!is.null(smoothing)) smoothing <- newnames(smoothing) 
    
    if (print.flag == TRUE) cat(
      "**************************************************************
 The numeric variable(s): ", 
      paste(names(data)[numtocat.obj$ind], collapse = ", "), 
      "\n will been synthesised as grouped variables and their numeric
 values generated from boostrap samples within categories.
**************************************************************\n", sep = "") 
  }
  ##-----------------end of--code for numtocat----------------------------
  
  
  ##-----------------------check.visit.sequence.syn-----------------------
  
  check.visit.sequence.syn <- function(setup){
    
    vis      <- setup$visit.sequence
    nvar     <- setup$nvar
    varnames <- setup$varnames
    method   <- setup$method
    
    # visit.sequence can include column indices only
    # not all columns have to be given - variables
    # not in the visit sequence won't be synthesised
    
    # stop if variable in visit.sequence more than once
    if (any(duplicated(vis))) stop("Visit sequence includes repeated variable names/numbers.\n", call. = FALSE)
    
    # remove any visitSequnce members outside 1:nvar
    if (any(!(vis %in% 1:nvar))) {
      cat("Element(s): {",paste(vis[!(vis %in% 1:nvar)],
                                collapse = ", "),"} of the 'visit.sequence' removed as not valid. No such column.\n\n", sep = "")
      vis <- as.numeric(vis[vis %in% 1:nvar])
    }
    
    # remove event indicator(s) from visit.sequence, if present
    event.in.vis <- !is.na(match(vis,event))
    if (!is.null(event) & any(event.in.vis) && method[which(event == vis[event.in.vis])] == "survctree") {
      cat("Variable(s) {", paste0(names(data)[vis][event.in.vis], collapse = ", "),
          "} with method(s) {",paste0(method[vis[event.in.vis]], collapse = ", "),
          "} removed from 'visit.sequence'\nbecause used as event indicator(s).\nAny event indicators will be synthesised along with the corresponding survival time(s). \n\n",
          sep = "")
      vis <- vis[!event.in.vis]
      if (length(vis) < 2) stop("Visit sequence now of length ",
                                length(vis),". No synthesis done.", call. = FALSE)
    } 
    #GRdenom new code
    #!BN adjusted to allow visit sequence with selected vars only 
    #! have to add a condition when denominator is not in visit seq at all;
    #! sampler has to be changed still
    #!---                                                                   
    #  check that denominator comes before the count for a binomial with denom
    #if (any(denom>0)) {
    #   denom.in.vis<-(1:nvar)[denom>0]
    #       for (j in denom.in.vis){
    #          posj<-(1:length(vis))[match(j,vis)]
    #          kj <-denom[j]
    #          posk<-(1:length(vis))[match(kj,vis)]
    #      if (posj<=posk) 
    #         stop("\n Denominator ",varnames[j]," for ",varnames[kj]," must be synthesisied before it\n")
    #   }
    #
    # }                                                               
    
    # check that denominator comes before the count for a binomial with denom 
    for (j in which(denom[vis] > 0)) { 
      denom.pos <- match(denom[vis][j],vis)
      if (j < denom.pos) stop("Denominator ",varnames[denom[vis][j]]," for ",
                              varnames[vis[j]]," must be synthesisied before it\n",
                              call. = FALSE)
    }
    #!---
    
    # stop if visit.sequence is empty
    if (length(vis) == 0) stop(paste("Seems no variables being synthesised.\nCheck parameter 'visit.sequence'."), call. = FALSE)
    
    # All variables used in passive synthesis have to be synthesised BEFORE
    # the variables they apply to
    for (j in which(is.passive(method[vis]))) {  #  
      var.present <- match(all.vars(as.formula(method[vis][j])),varnames) 
      var.in.vis  <- match(var.present,vis)
      if (j < max(var.in.vis) | any(is.na(var.in.vis))) stop("Variable(s) used in passive synthesis for ",
                                                             varnames[vis][j]," has/have to be synthesised BEFORE the variables they apply to.", call. = FALSE)
    }
    
    setup$visit.sequence <- vis
    return(setup)
  }
  ##-----------------end of--check.visit.sequence.syn---------------------
  
  
  ##-----------------------check.predictor.matrix.syn---------------------
  
  check.predictor.matrix.syn <- function(setup){
    ## checks the predictor.matrix
    ## makes consistency edits of the predictor.matrix
    
    pred     <- setup$predictor.matrix
    nvar     <- setup$nvar
    varnames <- setup$varnames
    vis      <- setup$visit.sequence
    method   <- setup$method
    denom    <- setup$denom                     #GRdenom new
    
    # set up default predictor matrix (if not provided by a user)
    # to lower diagonal in order of visitSequnce but with
    # elements for variables not to be synthesised set to 0
    
    pred.dt           <- matrix(0, nvar, nvar)
    pred.dt[vis, vis] <- outer(1:length(vis), 1:length(vis), ">")
    if (is.null(pred)) pred <- pred.dt
    
    # basic corrections for a default matrix or the one provided by a user
    dimnames(pred)   <- list(varnames, varnames)
    diag(pred)       <- 0
    
    # select from visit.sequence variables that are synthesised
    # (=method different than "")
    vis.syn <- vis
    if (!all(method == "") & length(method) > 1) vis.syn <- intersect(vis, which(method != ""))
    # removing predictors for variables with "" method
    if (length(vis.syn) < length(vis)) {
      vis.blank        <- setdiff(vis,vis.syn)
      pred[vis.blank,] <- 0
    }
    # removing predictors for variables not in visit.sequence
    pred[setdiff(1:nvar, vis),] <- 0
    
    # removing predictors for variables with "sample" method
    for (j in which(method == "sample")) pred[j,] <- 0
    
    # removing survival time from predictors
    for (j in which(method == "survctree")) pred[,j] <- 0
    
    # removing event indicator from predictors
    for (j in which(method == "survctree" & event > 0)) pred[,event[j]] <- 0
    #GRdenom new lines
    #  remove denom from prediction of its numerator
    for  (j in which(method == "logreg")) {
      if (denom[j] > 0) pred[j, denom[j]] <- 0
    }
    # to here
    # checking consistency between visit.sequence and predictor matrix
    # provided by a user: dropping from predictors variables that are
    # not already synthesised
    preddel <- which((pred[, vis.syn, drop = FALSE] == 1 & 
                        pred.dt[, vis.syn, drop = FALSE] == 0), arr.ind = TRUE)
    if (length(vis) > 1) {
      pred[,vis.syn] <- ifelse((pred[,vis.syn] == 1 & pred.dt[, vis.syn] == 0),
                               0, pred[, vis.syn])
      if (nrow(preddel) > 0) cat(paste("Not synthesised predictor ",
                                       varnames[vis.syn][preddel[, 2]],
                                       " removed from predictor.matrix for variable ",
                                       varnames[preddel[, 1]], ".\n", sep = ""))
    }
    setup$predictor.matrix <- pred
    setup$visit.sequence   <- vis
    return(setup)
  }
  ##-----------------end of--check.predictor.matrix.syn-------------------
  
  
  ##-----------------------check.method.syn------------------------------
  
  check.method.syn <- function(setup, data, proper) {
    ## check method, set defaults if appropriate
    
    method         <- setup$method
    default.method <- setup$default.method
    vis            <- setup$visit.sequence
    nvar           <- setup$nvar
    varnames	      <- setup$varnames
    pred		        <- setup$predictor.matrix
    event          <- setup$event
    denom          <- setup$denom                    
    
    # check that all ipf and allcat are at start of visit sequence   
    mcatall <- (method %in% "catall")[vis]
    mipf    <- (method %in% "ipf")[vis]
    if (any(mipf) & any(mcatall)) stop("Methods 'ipf' and 'catall' cannot both be used.\nIf you want all margins fitted for a set of variables,\nthen you could use 'ipf' and specify othmargins appropriately.\n", call. = FALSE)
    
    if (any(mcatall)) {
      if (any(mcatall != mcatall[order(!mcatall)])) stop("All variables with method 'catall' must be together at start of visit sequence.\n", call. = FALSE)
      if (sum(mcatall) == 1) {
        method[1] <- "sample"
        cat("First method changed to 'sample' from 'catall' as set for a single variable only.\n", call. = FALSE)
      }
    }
    if (any(mipf)) {
      if (any(mipf != mipf[order(!mipf)])) stop("All variables with method 'ipf' must be together at start of visit sequence.\n", call. = FALSE)
      if (sum(mipf) == 1) {
        method[1] <- "sample"
        cat("First method changed to 'sample' from 'ipf' as set for a single variable only.\n", call. = FALSE)
      }
    }
    
    # change method for constant variables but leave passive variables untouched
    # factors and character variables with missing data won't count,
    # as NA is made into an additional level
    for (j in 1:nvar) {
      if (!is.passive(method[j]) & method[j] != "ipf" & method[j] != "catall") {
        if (is.numeric(data[,j])) {
          v <- var(data[,j], na.rm = TRUE)
          if (!is.na(v)) constant <- (v < 1000 * .Machine$double.eps) else
            constant <- is.na(v) | v < 1000 * .Machine$double.eps
        } else {
          constant <- all(duplicated(data[,j])[-1L])
        }
        
        if (constant) {
          if (any(vis == j)) {
            method[j] <- "constant" 
            if (print.flag == T) cat('Variable ', varnames[j], 
                                     ' has only one value so its method has been changed to "constant".\n', sep = "")
            pred[j, ] <- 0
          }
          if (any(pred[, j] != 0)) { 
            if (print.flag == T) cat("Variable ", varnames[j], 
                                     " removed as predictor because only one value.\n", sep = "")
            pred[, j] <- 0
          }
        }
      }
    } 
    
    # check that passive relationship holds in original data
    #---
    passive.idx <- grep("~", method)
    for (i in passive.idx) {
      data.val <- data[,i]
      passive.val <- syn.passive(data, method[i])$res[[1]]
      
      if (is.factor(data.val)) {
        levels(data.val)[levels(data.val) == "NAtemp"] <- NA
        if (!all(levels(data.val) == levels(passive.val))) stop("Levels of passively created factor ",
                                                                varnames[i], " differ from original.\n",
                                                                sep = "", call. = FALSE)
      }
      
      NAderived <- sum( is.na(passive.val) & !is.na(data.val))
      NAorig    <- sum(!is.na(passive.val) &  is.na(data.val))
      nonmiss   <- sum(abs(as.numeric(passive.val)[!is.na(data.val) & !is.na(passive.val)] -
                             as.numeric(data.val)[!is.na(data.val) & !is.na(passive.val)]) > 1e-8 )
      
      if (sum(NAderived + NAorig + nonmiss) > 0) {
        cat("\nVariable(s) ", varnames[i]," with passive synthesis: relationship does not hold in data.\n", sep = "")
        if (NAderived > 0) cat("Total of ", NAderived," case(s) where value in data but some predictors missing.\n", sep = "")
        if (NAorig > 0) cat("Total of ", NAorig," case(s) where no predictors missing but no value in data.\n", sep = "")
        if (nonmiss > 0) cat("Total of ", nonmiss," case(s) where predictors do not give value present in data.\n", sep = "")
        cat("You might want to recompute the variable(s) in the original data.\n")
      }
      
      if (is.numeric(data.val) & any(is.na(data.val))) cat("\nVariable ", varnames[i],
                                                           " with passive synthesis has missing values\nso it will not be used to predict other variables.\n", sep = "")
    }
    
    
    # # check that passive variables obey rule in original data  ##GR0621
    # #---
    # passive.idx <- grep("~", method)
    # for (i in passive.idx) {
    #   test <- syn.passive(data, method[i])$res
    #   NAderived <- sum(is.na(test[[1]]) & !is.na(data[,i]))
    #   NAorig <-  sum(!is.na(test[[1]]) & is.na(data[,i]))
    #   nonmiss <- sum(abs(as.numeric(test[[1]])[!is.na(data[,i]) & !is.na(test)] - as.numeric(data[,i])[!is.na(data[,i]) & !is.na(test)]) > 1e-8 )
    #   
    #   if (sum(NAderived + NAorig + nonmiss) >0) {
    #     cat("\n\nVariable(s) ", varnames[i]," with passive synthesis: relationship does not hold in data\n" )
    #     if (NAderived >0 ) cat("Total of ", NAderived," cases where value in data but some predictors missing\n")
    #     if (NAorig >0 ) cat("Total of ", NAorig," cases no predictors missing but no value in  data \n")
    #     if (nonmiss >0 ) cat("Total of ", nonmiss," cases where predictors do not give value in data\n")
    #     stop("You must recompute the variables in the original data in order for the synthesis to run.", call. = FALSE)
    #   }
    #   if (is.factor(data[,i])) {
    #     resor <- data[,i]
    #     resor <- addNA(resor, ifany = TRUE)
    #     levels(resor)[is.na(levels(resor))] <- "NAtemp" 
    #     if (  !all(levels(resor) == levels(test[[1]]))) stop("Levels of passively created factor ", varnames[i]," differ from original\n", call. = FALSE)
    #   }
    #   if (is.numeric(data[,i]) & any(is.na(data[,i]))) cat("\n\nVariable ",varnames[i], " with passive synthesis has missing values\n so it will not be used to predict later variables.")
    # }
    
    
    
    
    # check nested variables
    #---
    nestmth.idx <- grep("nested", method)
    gr.vars <- vector("character", length(method))
    gr.vars[nestmth.idx] <- substring(method[nestmth.idx], 8)
    
    if (length(nestmth.idx) > 0) { 
      for (i in nestmth.idx) {
        # check if provided grouped var exists
        if (gr.vars[i] == "") stop("No variable name provided for 'nested' method for ", 
                                   varnames[i] ,".\nSet method as 'nested.varname' instead of 'nested'.\n", call. = FALSE)
        if (!(gr.vars[i] %in% varnames)) stop("Unrecognized variable ", gr.vars[i], 
                                              " provided for 'nested' method for ", varnames[i] ,"\n", call. = FALSE)
        if (gr.vars[i] == varnames[i]) stop("Variable ", varnames[i], 
                                            " can not be predicted by itself.\n", call. = FALSE) 
        
        # check if var nested in gr.var
        #? tabvars   <- table(data[,i], data[,gr.vars[i]]) 
        tabvars <- table(data[,i], data[,gr.vars[i]], useNA = "ifany") 
        tabvars01 <- ifelse(tabvars > 0, 1, 0)
        ymulti <- rowSums(tabvars01) > 1
        if ("NAtemp" %in% names(ymulti)) ymulti["NAtemp"] <- FALSE  
        ymulti[names(ymulti) %in% cont.na[[i]]] <- FALSE  # missing values and cont.na are excluded 
        if (any(ymulti)) cat("\nNOTE: Variable ", varnames[i], 
                             " is not nested within its predictor ", gr.vars[i], ".\nCheck values of ", 
                             varnames[i], ": ", paste0(rownames(tabvars01)[ymulti], collapse = ", "), 
                             "\n\n", sep = "")
        
        # adjust predictor matrix
        pred[i, -match(gr.vars[i], varnames)] <- 0  # remove all predictors except the group var
        pred[-match(varnames[i], gr.vars), i] <- 0  # exclude detailed var from predictors except when used for nested method
      }
      if (m > 0) method[nestmth.idx] <- "nested"
    }
    #---
    
    # check if var has predictors
    if (sum(pred) > 0) has.pred <- apply(pred != 0, 1, any)   # GR condition added
    else has.pred <- rep(0, nvar) 
    
    if (any(method == "parametric")) {
      # set method for first in visit.sequence to "sample"
      # change to default methods for variables with predictors
      
      if (length(vis) > 1) {
        for (j in vis[-1]) {
          if (has.pred[j]) {
            y <- data[,j]
            if (is.numeric(y))        method[j] <- default.method[1]
            else if (nlevels(y) == 2) method[j] <- default.method[2]
            else if (is.ordered(y) & nlevels(y) > 2) method[j] <- default.method[4]
            else if (nlevels(y) > 2)  method[j] <- default.method[3]
            else if (is.logical(y))   method[j] <- default.method[2]
            else if (nlevels(y) != 1) stop("Variable ",j," ",varnames[j],
                                           " type not numeric or factor.", call. = FALSE) # to prevent a constant values failing
          } else if (method[j] != "constant") method[j] <- "sample" 
        }
      }
    }
    
    # check whether the elementary synthesising methods are available
    # on the search path
    active    <- !is.passive(method) & !(method == "") & !(method == "constant") 
    if (sum(active) > 0) {
      # fullNames <- paste("syn", method[active], sep=".") #!GR-29/8/16
      fullNames <- method[active]                                #!GR-29/8/16
      if (m == 0) fullNames[grep("nested",fullNames)] <- "nested"  #!GR-29/8/16
      fullNames <- paste("syn", fullNames, sep = ".")               #!GR-29/8/16
      notFound  <- !(fullNames %in% c('syn.bag', 'syn.cart', 'syn.cartbboot', 'syn.collinear', 
                                      'syn.ctree', 'syn.cubertnorm', 'syn.lognorm', 'syn.logreg', 'syn.nested', 'syn.norm', 
                                      'syn.normrank', 'syn.pmm', 'syn.polr', 'syn.polyreg', 'syn.ranknorm', 'syn.rf', 
                                      'syn.sample', 'syn.satcat', 'syn.smooth', 'syn.sqrtnorm', 'syn.survctree') | 
                       sapply(fullNames, exists, mode = "function", inherit = TRUE)) 
      if (any(notFound)) stop(paste("The following functions were not found:",
                                    paste(unique(fullNames[notFound]), collapse = ", ")), call. = FALSE)
    }
    
    # type checks on built-in  methods 
    
    for (j in vis) {
      y     <- data[,j]
      vname <- colnames(data)[j]
      mj    <- method[j]
      mlist <- list(m1 = c("logreg","polyreg","polr","ipf","catall"), #GRBN
                    m2 = c("norm","normrank","survctree"),
                    m3 = c("norm","normrank","survctree","logreg"))
      # In case of type mismatch stop execution
      
      #                                                       #GRdenom lines changed
      # check for logistic with denominator
      # 
      if (denom[j] > 0) {
        if (!(mj %in% c("logreg"))) {
          method[j] <- "logreg"
          cat("Variable ", vname," has denominator (", colnames(data[denom[j]]), 
              ") and method ", mj, " has been changed to logreg\n", sep = "")
        }
        #if (!(mj %in% c("logreg"))) stop("Variable ", vname," has denominator (",
        #  colnames(data[denom[j]]), ") and method should be set to logreg and not ",mj,"\n", 
        #  call. = FALSE)
        #
        #  check all integers
        
        if (!((is.integer(y) | all((y - round(y)) == 0, na.rm = TRUE)) &  #!!!!!! address missing data issue
              (is.integer(data[denom[j]]) | all((data[denom[j]] - round(data[denom[j]]) == 0), na.rm = TRUE)))) #!!!!!! address missing data issue   
          stop("Variable ", vname," and denominator ", colnames(data[denom[j]]),
               " must be integers\n", call. = FALSE)
        if (any((data[denom[j]] - y) < 0, na.rm = TRUE)) stop("Variable ", vname,    #!!!!!! address missing data issue
                                                              " must be less than or equal denominator ",
                                                              colnames(data[denom[j]]),"\n", call. = FALSE)
      } else {
        if (is.numeric(y) & (mj %in% mlist$m1) & !(j %in% numtocat)) {                                             #!GRipf    numtocat added
          stop('Type mismatch for variable ', vname,
               '.\nSynthesis method "', mj, 
               '" is for categorical data unless grouped with numtocat.',
               sep = "", call. = FALSE)
        } else if (is.factor(y) & nlevels(y) == 2 & (mj %in% mlist$m2)) {
          stop('Type mismatch for variable ', vname,
               '.\nSyhthesis method "', mj, '" is not for factors.',
               sep = "", call. = FALSE)
        } else if (is.factor(y) & nlevels(y) > 2 & (mj %in% mlist$m3)) {
          stop('Type mismatch for variable ', vname,
               '.\nSyhthesis method "', mj,
               '" is not for factors with three or more levels.',
               sep = "", call. = FALSE)
        }                                               
      }
    }
    
    # check method for variables without predictors
    # set to "sample" if method is not valid
    
    # check if var has predictors (re-compute it)
    if (sum(pred) > 0) has.pred <- apply(pred != 0, 1, any)  #  GR condition added
    else has.pred <- rep(0, sqrt(length(pred)))            # this needed in case pred now has dimension 1
    
    for (j in vis) {
      if (!has.pred[j] & substr(method[j], 1, 6) != "nested" & is.na(any(match(method[j],
                                                                               c("", "constant", "sample", "sample.proper", "catall", "ipf"))))) 
      {
        if (print.flag == TRUE) cat('\nMethod "', method[j],
                                    '" is not valid for a variable without predictors (',
                                    names(data)[j],')\nMethod has been changed to "sample"\n\n', sep = "")
        method[j] <- "sample"
      }
    }
    
    # check survival method and events are consistent
    error.message <- "Invalid event value, must be logical, factor (2-level), or numeric (1/0)." 
    if (any(method == "survctree")) {
      for (j in vis) {   # checks for survival variables
        vname <- colnames(data)[j]
        if (method[j] == "survctree") {
          if (!is.numeric(data[,j])) stop("Variable ", vname,
                                          " should be a numeric survival time.", call. = FALSE)
          if (any(!is.na(data[,j]) & data[,j] < 0)) stop("Variable ", vname,          
                                                         " should be a non-negative survival time.", call. = FALSE)
          
          if (is.na((match(event[j], 1:nvar)))) {
            cat("Variable ", vname, " is a survival time. Corresponding event not in data, assuming no censoring.\n\n", sep = "")
            event[j] <- -1 # used to indicate no censoring
          } else {
            if (any(is.na(data[, event[j]]))) stop("Missing values in event indicator '", colnames(data)[event[j]], 
                                                   "' for survival time '", vname, "'. No data synthesised.", call. = FALSE)
            
            if (is.character(data[, event[j]])) {
              stop(error.message, call. = FALSE)
            } else if (is.logical(data[, event[j]])) {
              tabEI <- table(as.numeric(data[, event[j]]))
            } else if (is.factor(data[, event[j]])) {
              tabEI <- table(as.numeric(data[, event[j]]) - 1)
              cat("Value", levels(data[, event[j]])[2], "of event indicator",
                  colnames(data)[event[j]], "assumed to indicate an event.\n\n")
            } else {
              tabEI <- table(data[, event[j]])
            }
            
            if (length(tabEI) != 2) {
              if (length(tabEI) == 1 & all(tabEI == 1)) cat("Variable ", vname,
                                                            " is a survival time with all cases having events.\n", sep = "")
              else if (length(tabEI) == 1 & all(tabEI == 0)) stop("Variable ",
                                                                  vname," is a survival time with no cases having events.\n",
                                                                  "Estimation not possible.", sep = "", call. = FALSE)
              else stop(error.message, call. = FALSE)
            }
            if (!all(as.character(names(tabEI)) == c("0","1"))) {
              stop(error.message, call. = FALSE)
            }
          }
        } else {
          # checks for non-survival variables with events
          if (event[j] != 0) {
            cat("Variable ", vname, " has event set to ", colnames(data)[event[j]],
                ' although method is "', method[j], '". Event indicator reset to none.\n', sep = "")
            event[j] <- 0
          }
        }
      }
    } else if (!all(event == 0)) {
      cat("No variables have a survival method, so all event indicators are ignored.\n")
      event <- rep(0, nvar)
    }
    
    ## change names for proper imputations and check
    #for(j in unique(vis)){
    #  if(proper==T & method[j]!="") method[j] <- paste(method[j],
    #                                                   ".proper",sep="")
    #}
    
    # check collinearity of variables
    if (sum(pred > 0)  & m > 0) {                                     
      inpred <- apply(pred != 0, 1, any) | apply(pred != 0, 2, any)
      if (any(inpred)) {
        collout <- collinear.out(data[, inpred, drop = FALSE])
        if (length(collout) > 0) {
          for (i in 1:length(collout)) {
            if (print.flag) cat("Variables ", paste(collout[[i]], collapse = ", "),
                                " are collinear.", sep = "")
            vars <- match(collout[[i]], varnames[vis])
            vfirst <- collout[[i]][vars == min(vars)]
            nfirst <- match(vfirst,varnames)
            nall   <- match(collout[[i]],varnames)
            if (print.flag) cat(" Variables later in 'visit.sequence'\nare derived from ",
                                vfirst, ".\n\n", sep = "")
            for (ii in nall) {
              if (ii != nfirst) {
                method[ii] <- "collinear"
                pred[ii,]  <- 0
                pred[,ii]  <- 0
                pred[ii, nfirst] <- 1
              }
            }
          }
        }                 
      } 
    } 
    
    setup$event            <- event
    setup$method           <- method
    setup$predictor.matrix <- pred
    setup$visit.sequence   <- vis
    setup$denom            <- denom                  
    
    return(setup)
  }
  ##--------------------end-of--check.method.syn-------------------------
  
  
  ##------------------check.rules.syn------------------------------------
  
  check.rules.syn <- function(setup, data) {
    
    rules      <- setup$rules
    rvalues    <- setup$rvalues
    pred       <- setup$predictor.matrix
    nvar       <- setup$nvar                                                               
    varnames   <- setup$varnames
    method     <- setup$method
    vis        <- setup$visit.sequence
    #browser()  
    # Check the syntax
    #------------------
    # check the length of the rules and corresponding values
    if (any(sapply(rules,length) != sapply(rvalues,length)))
      stop("The number of data rules for each variable should equal the number of corresponding values.\n  Check variable(s): ",
           paste(varnames[sapply(rules,length) != sapply(rvalues,length)], collapse = ", "), ".", call. = FALSE)
    
    # special characters 
    char.allowed <- c("","|","||","&","&&","==",">=","<=","<",">",
                      "!=","==-",">=-","<=-","<-",">-","!=-","=='",".",")","(",";","-",
                      "'","\"","\"(",")\"","'(",")'") #### . ( and ) added
    char.present <- paste(gsub("\\w"," ",unlist(rules)),collapse = " ") # remove word characters and concatenate
    char.present <- strsplit(char.present,"[[:space:]]+")[[1]]    # split into seperate characters
    char.wrong   <- !(char.present %in% char.allowed)             # identify unxepected characters
    #if (any(char.wrong)) stop("Unexpected character(s) in rules: ",paste(char.present[char.wrong],collapse=" "),".")
    
    # variables names (=a string before a special character) must be in varnames 
    rule.sep <- lapply(sapply(rules, strsplit, "[|&]"), unlist)       # split into seperate conditions
    get.vars <- lapply(rule.sep, function(x) gsub("[<>=!].*", "", x)) # remove everything after a special character
    #get.vars <- lapply(get.vars,function(x) gsub(" ","",x))          # remove spaces
    get.vars <- lapply(get.vars, trimws)                              # Remove leading and trailing spaces
    get.vars <- lapply(get.vars, function(x) gsub("[\\(\\)]", "", x)) # remove brackets  
    get.vars <- lapply(get.vars, function(x) gsub("is.na", "", x))    # remove function name
    get.vars <- lapply(get.vars, function(x) gsub("`", "", x))        # remove `
    get.vars <- lapply(get.vars, function(x) x[x != ""])              # remove empty strings  ?? why this
    
    vars.in.rules <- unique(unlist(get.vars))
    vars.wrong <- !(vars.in.rules %in% varnames)                   # identify unxepected variables
    
    if (any(vars.wrong)) stop("Unexpected variable(s) in rules: ",
                              paste(vars.in.rules[vars.wrong], collapse = ", "), ".", call. = FALSE)
    
    # remove rules with warning for ipf and catall
    vars.with.rules <- varnames[rules != ""]
    if (any(method[varnames %in% vars.with.rules] %in% c("catall","ipf"))){
      cat("\nRules cannot be used for variables synthesised by ipf or catall")
      cat("\nbut values can be restricted by defining structural zero cells\nwith ipf.structzero or catall.structzero parameter.\n")
      rules[method %in% c("catall","ipf")] <-  rvalues[method %in% c("catall","ipf")] <- ""
      cat("\nRules defined for variable(s) ",
          paste0(varnames[method %in% c("catall","ipf") & varnames %in% vars.with.rules], collapse = ", "),
          " have been deleted.\n\n", sep = "")
      setup$rules <- rules
      setup$rvalues <- rvalues
      if (all(rules == "")) {
        return(setup)
      }
    }
    
    if (any(char.wrong)) {
      cat("One of rules may not be correct. If this is the case compare your rules and Error below.\nOtherwise rules have been applied.\n") 
      rs <- unlist(rules); names(rs) <- varnames
      rs <- cbind(rs[rs != ""]); colnames(rs) <- ""
      cat("\nYour rules are:")
      print(rs); cat("\n")
    }
    
    # Check that missingness in the data obeys the rules in rules
    nonmissing <- vector("list", nvar)
    isfactor   <- sapply(data, is.factor)
    yes.rules <- sapply(rules, function(x) any(x != ""))
    lth.rules <- sapply(rules, length)
    for (i in 1:nvar) {
      if (yes.rules[i]) {
        for (r in 1:lth.rules[i]) {
          if (is.na(rvalues[[i]][r]) & !isfactor[i]) {
            nonmissing[[i]][r] <- with(data,sum(!is.na(data[eval(parse(text = rules[[i]][r])), i])))
          } else if (is.na(rvalues[[i]][r]) & isfactor[i]) {    # different for factors because <NA> is treated as a level
            #nonmissing[[i]][r] <- with(data,sum(!is.na(as.character(data[eval(parse(text=rules[[i]][r])),i]))))
            nonmissing[[i]][r] <- with(data,sum(as.character(data[eval(parse(text = rules[[i]][r])),i]) != "NAtemp" &
                                                  as.character(data[eval(parse(text = rules[[i]][r])),i]) != "NAlogical"))       
          } else {
            nonmissing[[i]][r] <- with(data,sum(data[eval(parse(text = rules[[i]][r])),i] != rvalues[[i]][r] |
                                                  is.na(data[eval(parse(text = rules[[i]][r])),i])))
          }
        }
      }
    }
    any.nonmissing <- sapply(nonmissing, function(x) any(x > 0))
    if (any(any.nonmissing) > 0) cat("\nUnexpected values (not obeying the rules) found for variable(s): ",
                                     paste(varnames[any.nonmissing > 0], collapse = ", "),
                                     ".\nRules have been applied but make sure they are correct.\n", sep = "")
    
    # Check visit sequence 
    # all variables used in missing data rules have to be synthesised BEFORE 
    # the variables they apply to
    var.position <- lapply(get.vars, function(x) match(unique(x),varnames))
    var.in.vis   <- lapply(var.position, function(x) if (length(x) == 0) {
      x <- 0
    } else if (any(is.na(match(x,vis)))) {
      x[!is.na(match(x, vis))] <- match(x, vis)
      x[is.na(match(x, vis))]  <- nvar
    } else {
      x <- match(x,vis)})
    max.seq      <- sapply(var.in.vis, max, na.rm = T)
    not.synth    <- match(1:nvar,vis)[!is.na( match(1:nvar,vis))] <= max.seq[!is.na( match(1:nvar,vis))]
    if (any(not.synth,na.rm = TRUE)) stop("Variable(s) used in missing data rules for ",
                                          paste(varnames[!is.na( match(1:nvar,vis))][not.synth & !is.na(not.synth)], collapse = " "),
                                          " have to be synthesised BEFORE the variables they apply to.", call. = FALSE)
    
    # Check if a variable with missing values predicts other variables only if its
    # missing values are a subset of the missing values of the predicted variables
    # and remove from a prediction matrix if not. 
    # It refers to missing values coded as NA, otherwise variable can be used as 
    # a predictor without restrictions.
    
    #for (i in 1:nvar){
    #  if (!is.na(rvalues[i])) data[with(data,eval(parse(text=rules[i]))),i] <- NA
    #}
    patternRules <- matrix(0, nrow = nrow(data), ncol = ncol(data))
    for (i in 1:nvar) {
      if (yes.rules[i]) {
        for (r in 1:lth.rules[i]) {
          if (is.na(rvalues[[i]][r])) patternRules[with(data,eval(parse(text = rules[[i]][r]))), i] <- 1
        }
      }
    }
    patternNA <- is.na(data) + 0
    patternNA <- ifelse(patternRules == patternNA, patternNA, 0)
    diffNAij  <- function(i, j, dataNA) sum(dataNA[, i] - dataNA[, j] < 0)
    diffNA    <- Vectorize(diffNAij, vectorize.args = list("i", "j"))
    predNA    <- outer(1:nvar, 1:nvar, diffNA, dataNA = patternNA)
    
    # predNAwrong <- which ((pred==1 & predNA>0),arr.ind=TRUE)
    # pred        <- ifelse((pred==1 & predNA>0),0,pred)
    # if(nrow(predNAwrong)>0) cat(paste("Missing values of variable ",
    # varnames[predNAwrong[,2]]," are not a subset of missing values of variable ",
    # varnames[predNAwrong[,1]]," and cannot be used as its predictor (removed).\n",sep=""),
    # "\n",sep="")
    
    setup$predictor.matrix <- pred
    
    return(setup)
  }
  ##-----------------end of--check.rules.syn----------------------------
  
  
  ##------------------namedlist------------------------------------
  # check args that should be provided as a named list 
  # and create list with elements for each variable 
  namedlist <- function(x, varnames = colnames(data), 
                        nvars = length(varnames), 
                        missval = NA, argname, argdescription = "", 
                        asvector = FALSE){
    if (is.null(x)) {
      x <- as.list(rep(missval, nvars))
    } else if (!is.list(x) | any(names(x) == "") | is.null(names(x))) {
      stop("Argument '", argname,"' must be a named list with names of selected ", 
           argdescription, " variables.", call. = FALSE)  
    } else {
      x.missval <- as.list(rep(missval,nvars))
      x.ind <- match(names(x), varnames)
      if (any(is.na(x.ind))) stop("Unrecognized variable names in '", 
                                  argname,"': ",paste(names(x)[is.na(x.ind)], collapse = ", "), call. = FALSE)
      # For 'event' and 'denom' check if denominators' name exist and 
      # change them to column indecies   
      if (argname %in% c("denom", "event") & is.character(argname)) {
        denom.ind <- lapply(x,match,varnames)
        if (any(is.na(denom.ind))) stop("Unrecognized variable(s) provided as ", argname, "(s): ", 
                                        paste(unlist(x)[is.na(denom.ind)], collapse = ", "), call. = FALSE) 
        x <- denom.ind
      }
      x.missval[x.ind] <- x
      x <- x.missval 
    }
    names(x) <- varnames
    if (asvector) x <- unlist(x)
    return(x)  
  }
  ##-----------------end of--namedlist-----------------------------
  
  
  #----------------------- now syn continues here ----------------------
  # Basic checks of provided parameters:
  # dimensions, consistency, replication, ...
  
  call <- match.call()
  nvar <- ncol(data)
  if (!is.na(seed) & seed == "sample") {
    seed <- sample.int(1e9, 1)
    # cat("No seed has been provided and it has been set to ", seed,".\n\n", sep="")
  }
  if (!is.na(seed)) set.seed(seed)
  
  if (!(is.matrix(data) | is.data.frame(data)))
    stop("Data should be a matrix or data frame.")
  if (nvar < 2) stop("Data should contain at least two columns.", call. = FALSE)
  
  # S U B S A M P L E   S I Z E
  if (k != nrow(data) & print.flag == TRUE & m > 0) {
    # if (k > nrow(data)) {
    #   cat("Warning: Subpopulation size (k=",k,") cannot be greater than the population size (",
    #       nrow(data),").\n","Synthetic data sets of same size as data will be produced.\n\n",sep="")
    #       k <- nrow(data)
    # } else
    cat("Sample(s) of size ", k, " will be generated from original data of size ",
        nrow(data),".\n\n", sep = "")
  }
  
  # M E T H O D S
  method <- gsub(" ", "", method) # remove any spaces in or around method
  # # must be the same length as visit.sequence
  # if (length(method) > 1 & length(method) != length(visit.sequence)) 
  #   stop(paste("The length of method (", length(method),
  #              ") must be the same length as the visit.sequence (",length(visit.sequence),").", sep = ""), 
  #        call. = FALSE) 
  
  # expand user's syhthesising method (single string) to all variables
  if (length(method) == 1) {
    if (is.passive(method)) stop("Cannot have a passive syhthesising method for every column.", call. = FALSE)
    method <- rep(method, nvar)
    if (!(method[1] %in% c("catall", "ipf"))) method[visit.sequence[1]] <- "sample"
    # set method to "" for vars not in visit.sequence
    method[setdiff(1:nvar, visit.sequence)] <- ""
  }
  # if user specifies multiple methods, check the length of the argument
  # methods must be given for all columns in the data
  if (length(method) != nvar) stop(paste("The length of method (", length(method),
                                         ") does not match the number of columns in the data (", nvar, ").", sep = ""), 
                                   call. = FALSE)
  
  
  # P R E D I C T O R   M A T R I X
  if (!is.null(predictor.matrix)) {
    if (!is.matrix(predictor.matrix)) {
      stop("Argument 'predictor.matrix' is not a matrix.", call. = FALSE)
    } else if (nvar != nrow(predictor.matrix) | nvar != ncol(predictor.matrix))
      stop(paste("The 'predictor.matrix' has ",nrow(predictor.matrix),
                 " row(s) and ", ncol(predictor.matrix),
                 " column(s). \nBoth should match the number of columns in the data (",
                 nvar, ").", sep = ""), call. = FALSE)
  }
  
  data     <- as.data.frame(data)
  varnames <- dimnames(data)[[2]]
  
  # Named lists: check args and create list with elements for each variables 
  # C O N T I N O U S  V A R S  W I T H  M I S S I N G  D A T A  C O D E S
  # S E M I - C O N T I N O U S  V A R S 
  semicont <- namedlist(semicont, missval = NA, argname = "semicont", 
                        argdescription = "semi-continuous")
  cont.na  <- namedlist(cont.na, missval = NA, argname = "cont.na", 
                        argdescription = "")
  # combine cont.na and semicont lists  
  cont.na.ini <- cont.na
  cont.na <- mapply(c, cont.na, semicont, SIMPLIFY = FALSE)
  cont.na <- lapply(cont.na, unique)
  # R U L E S  and  R V A L U E S 
  rules   <- namedlist(rules, missval = "", argname = "rules", 
                       argdescription = "")
  rvalues <- namedlist(rvalues, missval = NA, argname = "rvalues", 
                       argdescription = "")
  # S M O O T H I N G
  smoothing <- namedlist(smoothing, missval = "", argname = "smoothing", 
                         argdescription = "", asvector = TRUE)
  if (any(smoothing != "")) {
    varsmoothind <- which(smoothing != "")
    varnumind    <- which(sapply(data, is.numeric))
    smoothnumind <- match(varsmoothind, varnumind)
    if (any(is.na(smoothnumind)) & print.flag == TRUE)
      cat("\nSmoothing can only be applied to numeric variables.\nNo smoothing will be applied to variable(s): ",
          paste(varnames[varsmoothind[is.na(smoothnumind)]], collapse = ", "), "\n", sep = "")   
    smoothing[varsmoothind[is.na(smoothnumind)]] <- "" 
  }
  
  # D E N O M   
  denom <- namedlist(denom, missval = 0, argname = "denom", asvector = TRUE)
  # E V E N T
  event <- namedlist(event, missval = 0, argname = "event", asvector = TRUE)
  
  # Perform various validity checks on the specified arguments
  setup <- list(visit.sequence = visit.sequence,
                method = method,
                default.method = default.method,
                predictor.matrix = predictor.matrix,
                nvar = nvar,
                varnames = varnames, 
                rules = rules,
                rvalues = rvalues,
                cont.na = cont.na,                                         
                event = event,
                denom = denom)                   #GRdenom new
  
  setup <- check.visit.sequence.syn(setup)
  setup <- check.predictor.matrix.syn(setup)
  
  
  # C H A N G E  D A T A  T Y P E  &  M O D I F Y  F A C T O R  L E V E L S
  #---
  # apply only if in predictor matrix 
  # GR added condition and else
  if (!is.null(setup$predictor.matrix) & sum(setup$predictor.matrix > 0)) {
    inpred <- apply(setup$predictor.matrix != 0, 1, any)*(!(method %in% c("","sample"))) |                # GR added to allow null methods not affected
      apply(setup$predictor.matrix != 0, 2, any)  # if anywhere in predictor.matrix
  } else {
    inpred <- rep(FALSE, sqrt(length(setup$predictor.matrix)))
  }
  notevent <- is.na(match(1:nvar,setup$event))       # if not in event list
  
  # Convert any character variables into factors for variables in pred
  ischar    <- sapply(data,is.character)
  chartofac <- (ischar * inpred) > 0
  if (sum(chartofac) > 0) {
    cat("\nVariable(s):",paste0(varnames[chartofac], collapse = ", "),
        "have been changed for synthesis from character to factor.\n", sep = " ")
    for (j in (1:nvar)[chartofac]) data[,j] <- as.factor(data[,j]) 
  }
  
  # Changing numeric variables with fewer than 'minnumlevels' into factors
  #  Default for this now set to 1 (only those with a single level are changed)                       
  # this allows correct synthesis of any with only one value as well as missing values
  #  Warning if numeric vars with < 5 levels (20 too many as picks up months)
  #  Also only need to do this if variable in predictionMatrix
  #  and any inappropriate methods are changed to the default for factors
  nlevel      <- sapply(data, function(x) length(table(x)))
  ifnum       <- sapply(data, is.numeric)
  innumtocat  <- rep(FALSE,length(data)) 
  
  if (minnumlevels < 5 & any(nlevel > minnumlevels & nlevel <= 5 & ifnum & inpred & notevent)) {
    cat("Warning: In your synthesis there are numeric variables with 5 or fewer levels: ",                     
        paste0(varnames[nlevel > minnumlevels & nlevel <= 5 & ifnum & inpred & notevent], collapse = ", "), ".",
        "\nConsider changing them to factors. You can do it using parameter 'minnumlevels'.\n", sep = "")  
  }
  
  vartofactor <- which(nlevel <= minnumlevels & ifnum & inpred & notevent)
  for (j in vartofactor) data[,j] <- as.factor(data[,j])
  if (length(vartofactor) > 0) {
    cat("\nVariable(s): ", paste0(varnames[vartofactor], collapse = ", "),
        " numeric but with only ", minnumlevels, 
        " or fewer distinct values turned into factor(s) for synthesis.\n\n", sep = "")
    for (j in vartofactor) {
      if (setup$method[j] %in% c("norm","norm.proper",
                                 "normrank","normrank.proper")) {
        if (nlevel[j] == 2) setup$method[j] <- default.method[2]
        else setup$method[j] <- default.method[3]
        cat("Method for ",varnames[j]," changed to ",setup$method[j],"\n\n")
      }
    }
  }
  
  # Modifies a factor by turning NA into an extra level
  isfactor  <- sapply(data,is.factor)
  for (j in (1:nvar)[isfactor & inpred & notevent]) {
    data[,j] <- addNA(data[,j], ifany = TRUE)
    levels(data[,j])[is.na(levels(data[,j]))] <- "NAtemp"          
  } 
  
  islogicalNA <- sapply(data, function(x) (is.logical(x) & any(is.na(x))))   
  for (j in (1:nvar)[islogicalNA & inpred & notevent]) {
    data[,j] <- addNA(data[,j], ifany = TRUE)
    levels(data[,j])[is.na(levels(data[,j]))] <- "NAlogical"          
  }
  
  
  #---
  setup  <- check.method.syn(setup, data, proper)
  if (any(rules != "")) setup <- check.rules.syn(setup, data)
  
  method           <- setup$method
  predictor.matrix <- setup$predictor.matrix
  visit.sequence   <- setup$visit.sequence
  event            <- setup$event
  rules            <- setup$rules
  rvalues          <- setup$rvalues
  cont.na          <- setup$cont.na
  default.method   <- setup$default.method
  denom            <- setup$denom  
  
  ############################################################
  
  method[!(1:length(method) %in% visit.sequence)] <- ""  
  
  # Identify any factors with > maxfaclevels levels that are in 'visit.sequence'
  no.fac.levels   <- sapply(data, function(x) length(levels(x)))
  too.many.levels <- no.fac.levels > maxfaclevels
  notsampling <- !(grepl("nested", method) | grepl("sample", method) | grepl("satcat", method) | grepl("constant", method))
  if (any(inpred & too.many.levels & notsampling)) {
    stop("We have stopped your synthesis because you have factor(s) with more than\n",
         maxfaclevels," levels: ", paste0(varnames[inpred & too.many.levels]," (",
                                          no.fac.levels[inpred & too.many.levels],"). ", collapse = ", "), 
         "This may cause computational problems that lead to failures\nand/or long running times. ",
         "You can try continuing by increasing 'maxfaclevels'\nand perhaps by trying one or more of the following:\n",
         " - omitting this variable as a predictor in the 'predictor.matrix' for some\n   or all variables,\n",
         " - leaving it/them until the end of the 'visit.sequence',\n",
         " - combining categories for these variables to make fewer categories,\n",
         " - using or creating a grouping of each variable (as above) and setting the\n",
         "   method for the variable with many levels to 'nested' within the groups.\n\n", 
         call. = FALSE)
  }
  
  # Not used variables are identified and dropped if drop.not.used==T
  # reclculate inpred & notevent in case they have changed after
  # check.method and check.data
  if (sum(predictor.matrix) > 0) {                            # GR condition added
    inpred      <- apply(predictor.matrix != 0, 1, any) | apply(predictor.matrix != 0, 2, any) # if anywhere in predictor.matrix
    ispredictor <- apply(predictor.matrix != 0, 2, any)    # if used as predictor
  }
  else inpred <- ispredictor <- rep(0, sqrt(length(predictor.matrix))) 
  
  notinvs     <- is.na(match(1:nvar,visit.sequence))  # if not in visit.sequence
  notsynth    <- notinvs | (!notinvs & method == "")  # if not synthesised
  notevent    <- is.na(match(1:nvar,event))           # if not in event list
  
  # identify columns not used as events or predictors or in visitSequnce
  out <- !inpred & notevent & notsynth
  
  if (any(out) & print.flag == TRUE) {
    cat("\nVariable(s):", paste0(varnames[out], collapse = ", "),
        "not synthesised or used in prediction.\n", sep = " ")
    if (drop.not.used == T) cat("The variable(s) will be removed from data and not saved in synthesised data.\n\n")
    else cat("CAUTION: The synthesised data will contain the variable(s) unchanged.\n\n")
  }
  
  
  # remove columns not used from data and replace predictor matrix, visit sequence, nvar and others
  if (any(out) & drop.not.used == T) {
    if (sum(!out) == 0) stop("No variables left to be synthesised", call. = FALSE) ######to stop if all data excluded 
    newnumbers        <- rep(0,nvar)
    newnumbers[!out]  <- 1:sum(!out)
    visit.sequence    <- newnumbers[visit.sequence]
    visit.sequence    <- visit.sequence[!visit.sequence == 0]
    predictor.matrix  <- predictor.matrix[!out,!out]
    
    event[event != 0] <- newnumbers[event[event != 0]]
    event             <- event[!out]
    denom[denom != 0] <- newnumbers[denom[denom != 0]]
    denom             <- denom[!out]
    
    data              <- data[,!out]
    nvar              <- sum(!out)
    method            <- method[!out]
    varnames          <- varnames[!out]
    
    if (nvar == 1) {                             #  GR added  note having to reassign character vector
      cl <- class(data)
      data <- data.frame(data)
      if ( any(cl == "character")) data[,1] <- as.character(data[,1])
      names(data) <- varnames
    }
    cont.na          <- cont.na[!out]
    cont.na.ini      <- cont.na.ini[!out]      #BN13/11  
    semicont         <- semicont[!out]         #BN13/11  
    smoothing        <- smoothing[!out]
    rules            <- rules[!out]
    rvalues          <- rvalues[!out]
    var.lab          <- var.lab[!out]
    val.lab          <- val.lab[!out]
    
    # recalculate these
    if (sum(predictor.matrix > 0)) {                                 # GR condition added
      inpred <- apply(predictor.matrix != 0, 1, any) |
        apply(predictor.matrix != 0, 2, any)      # if anywhere in predictor.matrix
      ispredictor <- apply(predictor.matrix != 0, 2, any) # if used as predictor
    }
    else inpred <- ispredictor <- rep(0,sqrt(length(predictor.matrix))) 
    
    notinvs  <- is.na(match(1:nvar,visit.sequence)) # if not in visit.sequence
    notsynth <- notinvs | (!notinvs & method == "") # if not synthesised
    notevent <- is.na(match(1:nvar,event))          # if not in event list
  }
  
  # Print out info on variables not synthesised but used in prediction
  pred.not.syn <- (ispredictor & notsynth)
  if (sum(pred.not.syn ) > 0 & drop.pred.only == FALSE) pred.not.syn[pred.not.syn == TRUE] <- FALSE
  
  if (sum(pred.not.syn ) > 0 & print.flag == TRUE) {
    cat("\nVariable(s):", paste0(varnames[ispredictor & notsynth], collapse = ", "),
        "used as predictors but not synthesised.\n", sep = " ")
    if (drop.pred.only == T) {
      cat("The variable(s) will not be saved with the synthesised data.\n")
    } else {
      cat("CAUTION: The synthesised data will contain the variable(s) unchanged.\n")
    }
  } 
  
  if (sum(predictor.matrix) > 0){
    
    pm <- padMis.syn(data, method, predictor.matrix, visit.sequence,
                     nvar, rules, rvalues, default.method, cont.na, smoothing, event, denom)
    
    # Pad the Syhthesis model with dummy variables for the factors
    # p <- padModel.syn(data, method, predictor.matrix, visit.sequence,
    #                   nvar, rules, rvalues)
    p  <- padModel.syn(pm$data, pm$method, pm$predictor.matrix, pm$visit.sequence,
                       pm$nvar, pm$rules, pm$rvalues, pm$factorNA, pm$smoothing, pm$event, pm$denom)
    if (k != dim(data)[1]) {
      # create a non-empty data frame in case some variables are kept unsynthesised
      p$syn <- p$syn[sample(1:nrow(data), k, replace = TRUE),]
      dimnames(p$syn)[[1]] <- 1:k
    }
    if (sum(duplicated(names(p$data))) > 0)
      stop("Column names of padded data should be unique.", call. = FALSE)
    
    p$cont.na <- pm$cont.na  
    
  } else {
    
    p <- list(data = data,                        
              syn = data,
              predictor.matrix = predictor.matrix, 
              method = method, 
              visit.sequence = visit.sequence, 
              rules = rules,
              rvalues = rvalues,
              cont.na = cont.na, 
              event = event,                
              denom = denom,
              categories = NULL,
              smoothing = smoothing)         
    
    if (k != dim(data)[1]) {
      p$syn <- p$syn[sample(1:nrow(data), k, replace = TRUE),]   
      dimnames(p$syn)[[1]] <- 1:k
    }
  }
  
  if (m > 0) {
    #  syn <- list(m)
    #	 for(i in 1:m){
    #     syn[[i]] <- data
    #     if (k != dim(data)[1]) syn[[i]] <- syn[[i]][sample(1:dim(data)[1], k, replace = TRUE), ]
    #   }    
    syn <- vector("list",m)                        
    for (i in 1:m) {                               
      syn[[i]] <- data[0, ]                         
      syn[[i]][1:k, ] <- NA
      #if (k > 0){                                 #!BN-12/08/2016 - for syn.strata
      #syn[[i]] <- syn[[i]][1:k,]                 
      #dimnames(syn[[i]])[[1]] <- 1:k             
      #}                                           
    } 
  }
  else syn <- NULL
  
  synall <- sampler.syn(p, data, m, syn, visit.sequence, rules, rvalues, 
                        event, proper, print.flag, k, pred.not.syn, models, numtocat, ...)
  
  syn <- synall$syn
  fits <- synall$fits
  
  if (m == 1) {
    syn <- syn[[1]]
    fits <- fits[[1]]
  } 
  
  
  #-----------------------------------------------------------------------
  # restore the original NA's in the data
  # for(j in p$visit.sequence) p$data[(!r[,j]),j] <- NA
  
  names(method)         <- varnames
  names(visit.sequence) <- varnames[visit.sequence]
  
  # Put grouped numeric variables and their synthesising parameters 
  # back into their correct positions
  #---
  if (!is.null(numtocat)) { 
    out <- (length(method) - length(numtocat) + 1:length(numtocat))
    if (m == 1) {
      tocat <- match(numtocat, names(syn))
      syn[, tocat] <- syn[, out]
      syn <- syn[, -out]
    } else {
      for (i in 1:m) {
        tocat <- match(numtocat, names(syn[[1]]))
        syn[[i]][, tocat] <- syn[[i]][, out]
        syn[[i]] <- syn[[i]][,-out]
      }
    }
    
    # move their parameters to numeric versions
    method             <- method[-out]
    visit.sequence     <- visit.sequence[-out]
    smoothing[tocat]   <- smoothing[out]
    smoothing          <- smoothing[-out]
    cont.na.ini[tocat] <- cont.na.ini[out]
    cont.na.ini        <- cont.na.ini[-out]
    rules[tocat]       <- rules[out]
    rules              <- rules[-out]
    rvalues[tocat]     <- rvalues[out]
    rvalues            <- rvalues[-out]
    predictor.matrix   <- predictor.matrix[-out,-out]
  }
  
  # ---------------------------------------------------------------------- 
  # #!GR added 14/01/21
  # put numtofac variables back to numeric and chartofac back to character
  
  if (length(chartofac) > 0 ) {
    if (m == 1) {
      tochange <- (1:length(syn)) [chartofac]
      for (i in tochange) syn[,i] <- as.character(syn[,i])}
    else {for (j in 1:m) {
      tochange <- (1:length(syn[[1]])) [chartofac]
      for (i in tochange) syn[[j]][,i] <- as.character(syn[[j]][,i])
    }
    }
  }
  
  if (any(vartofactor)) {
    if (m == 1)  {
      tochange <- (1:length(syn))[vartofactor]
      for (i in tochange) syn[,i] <- as.numeric(as.character(syn[,i]))
    }
    else { for(j in 1:m) {
      tochange <- (1:length(syn[[1]]))[vartofactor]
      for (i in tochange) syn[[j]][,i] <- as.numeric(as.character(syn[[j]][,i]))
    }
    }
  }
  
  #---------------------------------------------------------------------
  #Tidy models
  
  #---
  if (models) {
    if (m == 1) { 
      # drop fits from dummies and nulls
      fitout <- sapply(fits, function(x) is.null(x) || (is.character(x) && x =="dummy")) | 
        grepl("orig.", names(fits))
      if (any(fitout)) fits <- fits[!fitout]
      # move the models for non-missing values to original position
      n.0 <- grepl("\\.0", names(fits))      
      if (any(n.0)) {
        fits[gsub("\\.0", "", names(fits[n.0]))] <- fits[n.0]
        fits <- fits[!n.0]
      }
    }
    if (m > 1) { 
      for (j in 1:m) {
        fitout <- sapply(fits[[j]], function(x) is.null(x) || (is.character(x) && x =="dummy")) |
          grepl("orig.", names(fits))
        if (any(fitout)) fits[[j]] <- fits[[j]][!fitout]
        n.0 <- grepl("\\.0",names(fits[[j]]))
        if ( any(n.0)) {
          fits[[j]][gsub("\\.0","",names(fits[[j]][n.0]))] <- fits[[j]][n.0]
          fits[[j]] <- fits[[j]][!n.0]
        }
      }
    }
  }
  
  # Save, and return, but don't return data
  #---
  syndsobj <- list(call = call,
                   m = m,
                   syn = syn,
                   method = method,
                   visit.sequence = visit.sequence,
                   predictor.matrix = predictor.matrix,
                   smoothing = smoothing,
                   event = event,
                   denom = denom,                  
                   #                 minbucket = minbucket,
                   proper = proper,
                   n = nrow(data),
                   k = k,
                   rules = rules,
                   rvalues = rvalues,
                   cont.na = cont.na.ini,          
                   semicont = semicont,            
                   drop.not.used = drop.not.used,
                   drop.pred.only = drop.pred.only,
                   models = fits,
                   seed = seed,
                   var.lab = var.lab,
                   val.lab = val.lab,
                   obs.vars = obs.vars,
                   numtocat = numtocat,                      #!GRipf 2 new parameters
                   catgroups = catgroups)
  # if (diagnostics) syndsobj <- c(syndsobj, list(pad = p))
  class(syndsobj) <- "synds"
  return(syndsobj)
}
#-----------------------------sampler.syn-------------------------------
# The sampler controls the generation of conditional distributions
# This function is called by syn()

sampler.syn <- function(p, data, m, syn, visit.sequence,
                        rules, rvalues, event, proper,
                        print.flag, k, pred.not.syn, 
                        models, numtocat,  ...)
{
  #--- Assign optional parameters (...) to appropriate synthesising function   
  dots  <- as.list(substitute(list(...)))[-1L]         
  meth.with.opt <- paste(c("cart", "cartbboot", "ctree", "survctree", "polyreg", 
                           "norm", "lognorm", "sqrtnorm", "cubertnorm", "normrank", "pmm",
                           "polr", "rf", "ranger", "bag", "ipf", "catall"), collapse = "\\.|")
  meth.check <- grep(meth.with.opt, names(dots), value = TRUE)
  args.err <- !(names(dots) %in% meth.check)
  if (any(args.err)) stop("Unknown optional parameter(s): ", 
                          paste(names(dots)[args.err], collapse = ", "),
                          "\nNote that they have to be method specific, e.g. 'ctree.minbucket' and NOT 'minbucket'\n", 
                          call. = FALSE)
  if (length(dots) == 0) {
    mth.args <- NULL
  } else {  
    #mth.args.dots <- strsplit(names(dots), "\\.")
    mth.args.dots <- regmatches(names(dots), regexpr("\\.", names(dots)), invert = TRUE)
    mth.dots  <- unique(lapply(mth.args.dots, "[[", 1))
    args.dots <- lapply(mth.args.dots, "[[", -1)
    mth.args  <- setNames(vector("list", length(mth.dots)), unlist(mth.dots))
    
    for (i in 1:length(mth.dots)) { 
      ind <- grep(mth.dots[[i]], names(dots))
      mth.args[[i]] <- setNames(dots[ind], args.dots[ind])
    } 
  } 
  #---
  
  fits <- NULL 
  
  if (m > 0) {
    if (models) fits <- rep(list(setNames(vector("list", length(p$method)),
                                          names(p$method))), m) 
    for (i in 1:m) {  # Synthesising loop
      if (print.flag & m > 1) cat("\nSynthesis number ", i, 
                                  "\n--------------------\n", sep = "")  
      if (print.flag & m == 1) cat("\nSynthesis\n-----------\n", sep = "")  
      
      # Code for methods that take more than one variable together: ipf & catall      
      #--------------------------------------------------------------------------
      rest.visit.sequence <- p$visit.sequence  # when no grouped methods used
      
      if (any(p$method %in% c("catall", "ipf"))) {
        ordmethod <- p$method[p$visit.sequence]
        grind <- (1:length(p$visit.sequence))[ordmethod %in% ordmethod[1]]
        
        ## to reorder any dummies for grouped variables
        if (any(names(p$visit.sequence) %in% 
                paste(names(p$visit.sequence[grind]), "1", sep = "."))) {  
          dumind <- (1:length(p$visit.sequence))[names(p$visit.sequence) %in% 
                                                   paste(names(p$visit.sequence[grind]), "1", sep = ".")]
          othind <- (1:length(p$visit.sequence))[-c(grind, dumind)]
          p$visit.sequence <- p$visit.sequence[c(grind, dumind, othind)]
          ordmethod <- p$method[p$visit.sequence]
        }
        
        grouped <- p$visit.sequence[ordmethod %in% ordmethod[1]]
        
        if (print.flag == TRUE) {
          if (length(rest.visit.sequence) > 0  && 
              ncol(data) - length(numtocat) > length(grouped)) {
            
            cat("First ", length(grouped), " variables (", 
                paste(names(grouped), collapse = ", "),
                ") synthesised together by method '", ordmethod[1], "'\n", sep = "")
            if (ordmethod[1] == "catall" && !is.null(mth.args) && 
                "epsilon" %in% names(mth.args$catall) && mth.args$catall$epsilon > 0)
              cat("Synthesis made differentially private with parameter epsilon of ",
                  mth.args$catall$epsilon,"\n",
                  "Note that only these first variables will be made differentially private.\n")
            if (ordmethod[1] == "ipf" && !is.null(mth.args) && 
                "epsilon" %in% names(mth.args$ipf) && mth.args$ipf$epsilon > 0) 
              cat("Synthesis made differentially private with parameter epsilon of ",
                  mth.args$ipf$epsilon,"\n",
                  "Note that only these first variables will be made differentially private.\n")
          } else {
            cat("All ", length(grouped), 
                " variables in the data synthesised together by method '", 
                ordmethod[1], "'\n", sep = "")
            
            if (ordmethod[1] == "catall" && !is.null(mth.args) && 
                "epsilon" %in% names(mth.args$catall) && 
                mth.args$catall$epsilon > 0) 
              cat("Synthesis made differentially private with parameter epsilon of ",
                  mth.args$catall$epsilon,"\n")
            if (ordmethod[1] == "ipf" && !is.null(mth.args) && 
                "epsilon" %in% names(mth.args$ipf) && 
                mth.args$ipf$epsilon > 0) 
              cat("Synthesis made differentially private with parameter epsilon of ",
                  mth.args$ipf$epsilon,"\n")
          }   
        }   
        x <- p$data[, grouped]
        if (!(ordmethod[1] %in% names(mth.args))) fun.args <- NULL else
          fun.args  <- mth.args[[ordmethod[1]]]
        f <- paste("syn", ordmethod[1], sep = ".")
        synfun <- do.call(f, args = c(list(x = x, k = k, 
                                           proper = proper), fun.args))
        p$syn[, grouped] <- synfun$res
        if (models) {
          fits[[i]][[grouped[1]]] <- synfun$fit
          for (j in 2:length(grouped)) fits[[i]][[grouped[j]]] <- 
              paste("See first in group:", names(grouped)[1])
        }
        
        rest.visit.sequence <- p$visit.sequence[-(1:length(grouped))]
        if (length(rest.visit.sequence) > 0 & print.flag & 
            ncol(data) - length(numtocat) > length(grouped)) cat("\nRemaining variables:\n")
      }
      
      # Other variables 
      #--------------------------------------------------------------------------
      if (length(rest.visit.sequence) > 0) {           
        prcount <- 0 # to get new lines come out right
        for (j in rest.visit.sequence) {
          
          theMethod <- p$method[j]
          # get optional parameters for theMethod if they are provided
          if (!(theMethod %in% names(mth.args))) fun.args <- NULL else            
            fun.args  <- mth.args[[theMethod]]                                     
          
          vname <- dimnames(p$data)[[2]][j]
          if (print.flag & theMethod != "dummy"  
              & j <= (ncol(data) - length(numtocat))) { 
            cat(" ", vname, sep = "")
            prcount <- prcount + 1
          }  
          if (print.flag & prcount %% 10 == 0 & 
              j <= (ncol(data) - length(numtocat))) cat("\n")                                                  
          
          ya  <-  1:nrow(p$data) 
          ypa <- 1:k    
          
          # ya = yavailable, ym = ymissing                                            
          if (any(p$rules[[j]] != "")) {
            com.rules  <- paste(p$rules[[j]], collapse = " | ")
            evalrul.y  <- with(p$data,eval(parse(text = com.rules)))
            ym         <- which(evalrul.y == TRUE & !is.na(evalrul.y))
            ya         <- setdiff(1:nrow(p$data), ym)                                  
            evalrul.yp <- with(p$syn,eval(parse(text = com.rules)))         
            ypm        <- which(evalrul.yp == TRUE & !is.na(evalrul.yp))        
            ypa        <- setdiff(1:nrow(p$syn), ypm)       
          }                                                                       
          
          # != "", != "dummy", != "passive"
          if (theMethod != "" & (!is.passive(theMethod)) & theMethod != "dummy" ) {
            
            if (theMethod %in% c("sample", "sample.proper", "constant")) {
              
              y   <- p$data[ya, j]
              if (is.factor(y)) y <- y[, drop = TRUE]
              xp  <- length(ypa)
              x   <- length(ya)
              nam <- vname
              f   <- paste("syn", theMethod, sep = ".")
              if (theMethod == "constant") {
                synfun <- do.call(f, args = list(y = y, xp = xp, ...))    
              } else if (is.numeric(y)) {
                synfun <- do.call(f, args = list(y = y, xp = xp,       
                                                 smoothing = p$smoothing[j], cont.na = p$cont.na[[j]], 
                                                 proper = proper, ...)) 
              } else {
                synfun <- do.call(f, args = list(y = y, xp = xp, 
                                                 proper = proper, ...)) 
              }
              p$syn[ypa, j]  <- synfun$res
              if (models) fits[[i]][[j]] <- synfun$fit
              
            } else {
              
              x    <- p$data[ya, p$predictor.matrix[j, ] == 1, drop = FALSE]
              xp   <- p$syn[ypa, p$predictor.matrix[j, ] == 1, drop = FALSE]
              y    <- p$data[ya, j]
              if (is.factor(y)) y <- y[, drop = TRUE]
              nam  <- vname
              f    <- paste("syn", theMethod, sep = ".") 
              if (!theMethod %in% c("collinear", "nested")) {   # nested needs added to allow missing values 
                #if(theMethod!="collinear"){                  
                keep <- remove.lindep.syn(x, y, ...)
                x    <- x[, keep, drop = FALSE]
                xp   <- xp[, keep, drop = FALSE]
              }                                            
              if (theMethod == "survctree") {
                if (p$event[j] == -1) yevent <- rep(1,length(y))                   
                else yevent  <- p$data[ya,p$event[j]]
                survres      <- do.call(f, args = c(list(y = y, yevent = yevent,
                                                         x = x, xp = xp, proper = proper), 
                                                    fun.args))
                p$syn[ypa, j] <- survres[[1]]                                # synthetic data survival goes to p$syn
                if (p$event[j] != -1) p$syn[ypa,p$event[j]] <- survres[[2]]  # synthetic data event goes to p$syn
                if (models) fits[[i]][[j]] <- survres$fit 
              } else if (theMethod == "logreg" & p$denom[j] != 0) {                 
                synfun <- do.call(f, args = list(y = y, x = x, xp = xp,
                                                 denom = p$data[ya,p$denom[j]], 
                                                 denomp = p$syn[ypa, p$denom[j]],       
                                                 proper = proper, ...))
                p$syn[ypa, j] <- synfun$res
                if (models) fits[[i]][[j]] <- synfun$fit           
                
              } else if (theMethod == "nested") {
                if (is.numeric(y)) {
                  synfun <- do.call(f, args = c(list(y = y, x = x, xp = xp,
                                                     smoothing = p$smoothing[j], cont.na = p$cont.na[[j]], 
                                                     proper = proper), fun.args))
                } else {
                  synfun <- do.call(f, args = c(list(y = y, x = x, xp = xp,
                                                     proper = proper), fun.args))
                }
                p$syn[ypa, j] <- synfun$res
                if (models) fits[[i]][[j]] <- synfun$fit
                
              } else {
                if (is.numeric(y)) {
                  synfun <- do.call(f, args = c(list(y = y, x = x, xp = xp,
                                                     smoothing = p$smoothing[j],
                                                     proper = proper), fun.args))
                } else {
                  synfun <- do.call(f, args = c(list(y = y, x = x, xp = xp,
                                                     proper = proper), fun.args))
                }
                p$syn[ypa, j] <- synfun$res
                if (models) fits[[i]][[j]] <- synfun$fit
              }
            }
            
            if (any(p$rules[[j]] != "")) {
              if (length(p$rules[[j]]) == 1 & length(ypm) > 0) {
                p$syn[ypm,j] <- p$rvalues[[j]] 
              } else {
                for (r in 1:length(p$rules[[j]])) {
                  revalrul.yp  <- with(p$syn,eval(parse(text = p$rules[[j]][r])))  
                  rypm <- which(revalrul.yp == TRUE & !is.na(revalrul.yp))
                  if (length(rypm) > 0) p$syn[rypm,j] <- p$rvalues[[j]][r]
                }
              }                 
            }  
          } # end of !="", !="dummy", !="passive"
          
          else if (is.passive(theMethod)) {
            class0 <- class(p$syn[,j])
            synfun <- syn.passive(data = p$syn, func = theMethod) 
            
            if (is.factor(synfun$res[[1]]) & any(is.na(synfun$res[[1]]))) {
              synfun$res[[1]] <- addNA(synfun$res[[1]], ifany = TRUE)
              levels(synfun$res[[1]])[is.na(levels(synfun$res[[1]]))] <- "NAtemp"
            }
            
            p$syn[, j] <- synfun$res
            class(p$syn[,j]) <- class0
            if (models) fits[[i]][[j]] <- synfun$fit 
          }
          
          else if (theMethod == "dummy") {    # replace dummy variables in p$syn
            # getting dummy values from a synthesised categorical variable
            cat.columns <- p$syn[, p$categories[j, 4]]  # this is the single column with the data for which this is the dummy
            model.frame(~cat.columns - 1, data = p$syn) 
            p$syn[, (j:(j + p$categories[p$categories[j, 4], 2] - 1))] <- # replaces all the dummies for this variable with
              matrix((model.matrix(~cat.columns - 1)[, -1]),                # dummies calculated from the synthesised data
                     ncol = p$categories[p$categories[j, 4], 2],
                     nrow = nrow(p$syn))
            p$syn[,j] <- as.numeric(p$syn[, j])
            remove("cat.columns")
            if (models) fits[[i]][[j]] <- "dummy"    
          }
        } # end j loop 
      } # end other variables
      if (print.flag) cat("\n")  
      
      #if (k==dim(data)[1]) syn[[i]] <- p$syn[,1:dim(data)[2]]
      #else syn[[i]] <- p$syn[sample(1:dim(data)[1],k),1:dim(data)[2]]
      
      syn[[i]] <- p$syn[, 1:dim(data)[2], drop = FALSE]
      nms <- names(data)
      # exclude unsynthesised if drop.pred.only set to true
      if (sum(pred.not.syn ) > 0) {
        syn[[i]] <- syn[[i]][, !pred.not.syn]
        nms <- nms[!pred.not.syn]  # GR save names to use below if data just one column
      }
      # Prevent a single character column being changed to a factor
      chgetochar <- (sum(!pred.not.syn) == 1 & any(class(syn[[i]][, 1]) == "character"))       
      syn[[i]] <- as.data.frame(syn[[i]])
      if (chgetochar) {
        syn[[i]][, 1] <- as.character(syn[[i]][, 1])
        names(syn[[i]]) <- nms
      }
      
      #turn NA level in factors / logical to missing NA's
      # and remove contrasts 
      for (j in (1:ncol(syn[[i]]))) {
        if (is.factor(syn[[i]][,j])) {                                    #!BN-20/04/16
          if ("NAlogical" %in% levels(syn[[i]][,j])) {
            levels(syn[[i]][,j])[levels(syn[[i]][,j]) == "NAlogical"] <- NA
            syn[[i]][,j] <- as.logical(syn[[i]][,j])
          } else {                                                            
            # syn[[i]][,j] <- factor(syn[[i]][,j],exclude=NA,levels=levels(syn[[i]][,j]))
            levels(syn[[i]][,j])[levels(syn[[i]][,j]) == "NAtemp"] <- NA  #!BN 10/08/15 
          }
          #! attributes(syn[[i]][,j])$contrasts <- NULL                    #!BN-28/04/16   UNCOMMENT???? 
        }                                                                       
      }
    } # end i loop (m)
  } # end synthesising (m > 0)
  
  return(list(syn = syn, fits = fits))
}


###-----remove.lindep.syn--------------------------------------------------

remove.lindep.syn <- function(x, y, eps = 0.00001, maxcor = 0.99999, 
                              allow.na = FALSE, ...) 
{
  if (ncol(x) == 0) return(NULL) 
  if (eps <= 0) stop("\n Argument 'eps' must be positive.", call. = FALSE)
  xobs <- sapply(x, as.numeric)                                       
  yobs <- as.numeric(y)
  keep <- unlist(apply(xobs, 2, var) > eps)
  keep[is.na(keep)] <- FALSE
  keep <- keep & suppressWarnings((unlist(apply(xobs, 2, cor, yobs)) < maxcor)) # if y includes NA -> NAs error
  if (all(!keep)) warning("\nAll predictors are constant or have too high correlation.\n")
  ksum <- sum(keep)
  cx   <- cor(xobs[, keep, drop = FALSE], use = "all.obs")
  eig  <- eigen(cx, symmetric = TRUE)
  ncx  <- cx
  while (eig$values[ksum]/eig$values[1] < eps) {
    j   <- (1:ksum)[order(abs(eig$vectors[, ksum]), decreasing = TRUE)[1]]
    keep[keep][j] <- FALSE
    ncx  <- cx[keep[keep], keep[keep], drop = FALSE]
    ksum <- ksum - 1
    eig  <- eigen(ncx)
  }
  # if (!all(keep)) cat("\tVariable(s): ", paste(dimnames(x)[[2]][!keep], collapse = ", "),
  #   " removed due to linear dependency",sep="")
  return(keep)
}
