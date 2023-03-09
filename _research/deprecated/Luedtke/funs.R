## NOTE TO SELF: Need to update for observational setting -- or even clinical trials where g0 != 1/2 (I wasn't careful about this)

require(SuperLearner)
require(hitandrun) # get uniform simplex sample from this package
require(nloptr)
# source('~/optimal_tx/code/boundedlogistic.R')
# source('~/optimal_tx/code/SL_fits.R')

## Qbar0.fun
## INPUTS
## A : treatment
## W1,W2,W3 : covariates
## OUTPUT
## An mean outcome Y in this covariate-treatment strata

Qbar0.fun = function(A,W){
	W1 = W[,1]
	W2 = W[,2]
	W3 = W[,3]
	W4 = W[,4]
	return(0.5*plogis(1-W1^2 + 3*W2 + 5*W3^2*A - 4.45*A)+0.5*plogis(-0.5- W3 + 2*W1*W2 + 3*abs(W2)*A - 1.5*A))}
#	return(0.3*plogis(3.5-W1^2 + 3*W2 + 2*exp(W3)*A - 5*A)+0.7*plogis(-1.5 - W3 + 2*W1*W2 + 2.5*abs(W2)*A - 1.07*A))}
# have lots of data for the one below
#	return(0.3*plogis(2-W1^2 + 3*W2 + exp(W3)*A - 5*A)+0.7*plogis(-1 - W3 + 2*W1*W2 + 3*abs(W2)*A - 1.03*A))}

# return(0.3*plogis(-W1^2 + 3*W2 + exp(0.5*W3)*A - 6*A)+0.7*plogis(-1 - W3 + 2*W1*W2 + 3*abs(W2)*A - 1.03*A))}

Qbar0.fun.easy2 = function(A,W){
	W1 = W[,1]
	W2 = W[,2]
	W3 = W[,3]
	W4 = W[,4]
	return(0.45 + 0.05*W1 -0.05*W2 + 0.2*A*W2*W3 + 0.1*A*W2^2 + 0.15*A*W3-0.033*A)}
	
## sim.data
## INPUTS
## n : sample size
## Qbar0.sim : a function to generate the mean outcome in a covariate-treatment strata
## g : treatment probability in an RCT
## binom : generates binomial Y if true, Gaussian Y if false
## W.fun : a function to generate W's which takes n (sample size) as input
## OUTPUT
## A data frame containing:
##	ObsData : the observed data
##	FullData : the full data (ObsData + counterfactual outcomes)

sim.data = function(n,Qbar0=Qbar0.fun,g=0.5,binom=TRUE,W.fun=rnorm){
	W1=W.fun(n)
	W2=W.fun(n)
	W3=W.fun(n)
	W4=W.fun(n)
	W=data.frame(W1=W1,W2=W2,W3=W3,W4=W4)
	A = rbinom(n,1,g)
	if(binom){
		gen.Y = function(W1,W2,W3,A){
			u = runif(n)
			return(cbind(as.numeric(u<Qbar0(A,W)),as.numeric(u<Qbar0(0,W)),as.numeric(u<Qbar0(1,W))))}
	} else {
		gen.Y = function(W1,W2,W3,A){
			z = rnorm(n)
			return(cbind(Qbar0(A,W1,W2,W3)+z,Qbar0(0,W1,W2,W3)+z,Qbar0(1,W1,W2,W3))+z)}
	}
	Y.mat = gen.Y(W1,W2,W3,A)
	Y = Y.mat[,1]
	Y0 = Y.mat[,2]
	Y1 = Y.mat[,3]
	return(list(ObsData = data.frame(W1=W1,W2=W2,W3=W3,W4=W4,A=A,Y=Y),FullData = data.frame(W1=W1,W2=W2,W3=W3,W4=W4,A=A,Y=Y,Y0=Y0,Y1=Y1)))
}


## qbarV.mse
## Fit MSE estimate of qbarV, where we take V to be W1 and W2
## INPUTS
## W : data frame of covariates
## A : vector of treatments
## Y : vector of outcomes
## SL.library: SuperLearner library for estimating Qbar(V)
## Qbar0 : a function to generate the mean outcome in a covariate-treatment strata
## g : treatment probability in an RCT
## V : a vector of column names in ObsData upon which treatment decisions will be made
## newX : Include a different newX in the SL fit?
## vR : Use this to pass pre-specified folds to SuperLearner. See SL documentation
## wgts : Weights to be passed to SuperLearner for the multiple time point case
## OUTPUT
## 

qbarV.mse = function(W,A,Y,SL.library,Qbar0=Qbar0.fun,g=0.5,V=c('W1','W2'),newX=NULL,folds=NULL,wgts=NULL){
	D1 = (2*A-1)/g * (Y-Qbar0(A,W)) + Qbar0(1,W) - Qbar0(0,W)
	if(is.null(folds)){ SL.fit = SuperLearner(Y=D1,X=subset(W,select=V),SL.library=SL.library,family=gaussian(),newX=newX,obsWeights=wgts)
	} else {
		SL.fit = SuperLearner(Y=D1,X=subset(W,select=V),SL.library=SL.library,family=gaussian(),newX=newX,id=folds,obsWeights=wgts,cvControl=list(V=max(folds)))
	}
	return(SL.fit)
}


## qbarV.logit
## Fit MSE estimate of qbarV, where we take V to be W1 and W2
## INPUTS
## W : data frame of covariates
## A : vector of treatments
## Y : vector of outcomes
## SL.library: SuperLearner library for estimating Qbar(V)
## Qbar0 : a function to generate the mean outcome in a covariate-treatment strata
## g : treatment probability in an RCT
## bd : bounds on QbarV (called (a,b) in the write-up)
## V : a vector of column names in ObsData upon which treatment decisions will be made
## newX : Include a different newX in the SL fit?
## vR : Use this to pass pre-specified folds to SuperLearner. See SL documentation
## wgts : Weights to be passed to SuperLearner for the multiple time point case
## OUTPUT
## 

qbarV.logit = function(W,A,Y,SL.library,Qbar0=Qbar0.fun,g=0.5,bd=c(-1,1),V=c('W1','W2'),newX=NULL,folds=NULL,wgts=NULL){
	D1 = (2*A-1)/g * (Y-Qbar0(A,W)) + Qbar0(1,W) - Qbar0(0,W)
	D1.ab = (D1-bd[1])/(bd[2]-bd[1])
	if(sum(D1.ab<0)>0 | sum(D1.ab>1)>0) warning('Some D1 values out of (0,1). Only the linear logistic regression fits will work.')
	if(is.null(folds)){
		SL.fit = SuperLearner(Y=D1.ab,X=subset(W,select=V),SL.library=SL.library,family=binomial(),newX=newX,obsWeights=wgts)
	} else {
		SL.fit = SuperLearner(Y=D1.ab,X=subset(W,select=V),SL.library=SL.library,family=binomial(),newX=newX,id=folds,obsWeights=wgts,cvControl=list(V=max(folds)))
	}
	return(SL.fit)
}

## cl.setup
## Creates response and weights for the cl problem.
## INPUTS
## W : data frame of covariates
## A : vector of treatments
## Y : vector of outcomes
## Qbar0 : a function to generate the mean outcome in a covariate-treatment strata
## g : treatment probability in an RCT
## OUTPUT
## zA list containing:
## data : the outcome in the classification algorithm
## wgts : the weight in the classification algorithm.

cl.setup = function(W,A,Y,Qbar,g){
	n <- nrow(W)
	g.fun = function(AA,WW){ # **** change this line once generalized to non-RT context
		if(length(AA)==1) AA=rep(AA,nrow(WW))
		g.out = rep(NA,length(AA))
		g.out[AA==0] = 1-g
		g.out[AA==1] = g
		return(g.out)}
	S = (Y>=g.fun(0,W)*Qbar(rep(1,n),W) + g.fun(1,W)*Qbar(rep(0,n),W))
	Z = A*S + (1-A)*(1-S)
	gAW = g.fun(A,W)
	K = (1-2*S)*(-Y/gAW+Qbar(1-A,W) + (1/gAW-1)*Qbar(A,W))
	return(list(Z=Z,wgts=K))
}


## cl.SL
## Creates a SL object on the data as created by cl.setup.
## INPUTS
## W : data frame of covariates
## A : vector of treatments
## Y : vector of outcomes
## SL.library : SuperLearner library consisting of cl algorithms
##	that accept observation weights
## V : a vector of column names in ObsData upon which treatment decisions will be made
## Qbar0 : a function to generate the mean outcome in a covariate-treatment strata
## g : treatment probability in an RCT
## folds : If specified, a vector specifying the folds of the observations in W,A,Y
## vR : Use this to pass pre-specified folds to SuperLearner. See SL documentation
## mtp.wgts : Weights to be passed to SuperLearner for the multiple time point case
## OUTPUT
## A SuperLearner object

cl.SL = function(W,A,Y,SL.library,V,Qbar=Qbar0.fun,g=0.5,folds=NULL,mtp.wgts=NULL){
	cs.out = cl.setup(W,A,Y,Qbar,g)
	if(!is.null(mtp.wgts)) cs.out$wgts = cs.out$wgts * mtp.wgts
	if(is.null(folds)){
		out = list(SL.fit=SuperLearner(Y=cs.out$Z,X=subset(W,select=V),family=binomial(),SL.library=SL.library,obsWeights=cs.out$wgts))
	} else {
		out = list(SL.fit=SuperLearner(Y=cs.out$Z,X=subset(W,select=V),family=binomial(),SL.library=SL.library,obsWeights=cs.out$wgts,id=folds,cvControl=list(V=max(folds))))
	}
	class(out) = 'clSL'
	return(out)
}


## d.tmle
## Fits a TMLE given a decision vector.
## INPUT
## A : treatment
## Y : vector of outcomes
## Qb0.dopt : Qbar0(d,W)
## g : g(A|W), i.e. probability of an individual receiving the treatment they received
## d : Optimal treatment for each individual
## bound.Qb0 : should we use the logistic fluctuation model?
## OUTPUT
## 

d.tmle = function(A,Y,Qb0.dopt,g,d,bound.Qb0=TRUE,wgts=NULL){
	if(is.null(wgts)) wgts = rep(1,length(A))
	g.inv = (A==d)/g
	if(bound.Qb0){
		logit.Qb0 = qlogis(Qb0.dopt)
		eps = glm(Y~-1+offset(logit.Qb0) + g.inv,weights=wgts,family=binomial())$coef
		Qb0.star = plogis(logit.Qb0+eps*g.inv)
	} else {
		eps = glm(Y~-1+offset(Qb0.dopt)+g.inv,weights=wgts)$coef
		Qb0.star = Qb0.dopt + eps*g.inv
	}
	tmle.var = var(wgts*g.inv * (Y-Qb0.star) + wgts*Qb0.star - mean(wgts*Qb0.star))/length(A)
	return(list(psi=mean(wgts*Qb0.star),var.est=tmle.var))
}


## qbarV.alpha.cl
## Fits alpha according to a cl setup with a convex approximation to the 0-1 loss.
## INPUTS
## W : data frame of covariates
## A : vector of treatments
## Y : vector of outcomes
## jointSL.fit : a jointSL object
## Qbar : a function to generate the mean outcome in a covariate-treatment strata. Default gives an IPTW estimator
## g : treatment probability in an RCT
## bound.Qb0 : Does Qbar0(A,W) fall between 0 and 1?
## convex.loss : Loss to use in cl fit -- either 'hinge' or 'log'
## risk.est.method : How should the risk be estimated? A list containing a string for the method:
##	empirical risk ('emp.risk') or CV empirical risk ('cv.emp.risk')
## wgts : Weights to be passed to SuperLearner for the multiple time point case
## truncate.latent : either a cutoff at which the absolute value of latent functions should be truncated,
##	or NULL for no truncation
## OUTPUT
## 

qbarV.alpha.cl = function(W,A,Y,jointSL.fit,Qbar,g,bound.Qb0=TRUE,convex.loss='log',risk.est.method=list(method='emp.risk'),wgts=NULL,truncate.latent=NULL){
	if(is.null(wgts)) wgts=1
	if(is.null(jointSL.fit$SL.QbV)&is.null(jointSL.fit$SL.cl)) stop('At least one of SL.fit and SL.cl.fit must be specified.')
	if(!(convex.loss%in%c('hinge','log'))) stop('qbarV.alpha.cl currently only supports hinge or log loss.')

	if(risk.est.method$method=='emp.risk'){
		CV.preds = NULL
		if(!is.null(jointSL.fit$SL.QbV)) {
			if(is.null(jointSL.fit$SL.QbV$bd)) { QbV.CV.preds = jointSL.fit$SL.QbV$SL.fit$Z
			} else QbV.CV.preds = jointSL.fit$SL.QbV$SL.fit$Z*(jointSL.fit$SL.QbV$bd[2]-jointSL.fit$SL.QbV$bd[1]) + jointSL.fit$SL.QbV$bd[1]
			colnames(QbV.CV.preds) = paste('QbV_',names(jointSL.fit$SL.QbV$SL.fit$coef),sep='')
			CV.preds = cbind(CV.preds,QbV.CV.preds)
			num.QbV = ncol(QbV.CV.preds)
		} else num.QbV = 0
		if(!is.null(jointSL.fit$SL.cl)){
			cl.CV.preds = jointSL.fit$SL.cl$SL.fit$Z-1/2
			colnames(cl.CV.preds) = paste('cl_',names(jointSL.fit$SL.cl$SL.fit$coef),sep='')
			CV.preds = cbind(CV.preds,cl.CV.preds)
			num.cl = ncol(cl.CV.preds)
		} else num.cl = 0

		if(length(truncate.latent)>0) CV.preds = pmin(pmax(CV.preds,-truncate.latent),truncate.latent)

		cs.out = cl.setup(W,A,Y,Qbar,g)
		Z = cs.out$Z
		if(convex.loss=='hinge') { risk.fun = function(b){mean(wgts*cs.out$wgts*pmax(0,1-(2*Z-1)*(c(CV.preds%*%cbind(b)))))}
		} else if(convex.loss=='log') risk.fun = function(b){mean(wgts*cs.out$wgts*(-plogis((2*Z-1)*(c(CV.preds%*%cbind(b))),log.p=TRUE)))}
			#function(b){mean(wgts*cs.out$wgts*log(1+exp(-(2*Z-1)*(c(CV.preds%*%cbind(b))))))}
		num.alg = num.QbV + num.cl
		init = rep(1/num.alg,num.alg)
		names(init) = colnames(CV.preds)
		## **** rep(2,num.QbV) is inappropraiate when g is not 1/2 -- trying to make sure loss is bounded here
		alpha.out = optim(init,risk.fun,method='L-BFGS-B',lower=rep(0,num.alg),upper=rep(1,num.alg))$par
		# used to be: alpha.out = optim(init,risk.fun,method='L-BFGS-B',lower=rep(0,num.alg),upper=qlogis(0.975)*c(rep(2,num.QbV),rep(1/2,num.cl)))$par
		if(sum(alpha.out)!=0){
			return(alpha.out/sum(alpha.out))
		} else {
			return(alpha.out)}
	} else if(risk.est.method$method=='cv.emp.risk') {
		stop('Currently not supported. Not estimating Q or g in simulations, so irrelevant anyway.')
	} else stop('Invalid risk estimating method.')
}


## make.jointSL
## Makes a jointSL object
## INPUTS
## W : data frame of covariates
## A : vector of treatments
## Y : vector of outcomes
## QbV.SL.library : SuperLearner library for estimating Qbar(V). Can be NULL
## cl.SL.library : SuperLearner library for classification. Can be NULL
## Qbar0 : a function to generate the mean outcome in a covariate-treatment strata. Default gives an IPTW estimator
## g : treatment probability in an RCT
## V : a vector of column names in ObsData upon which treatment decisions will be made
## alpha.risk.est : Options are: empirical risk ('emp.risk') or TMLE ('tmle') or DR-IPCW
## bd : bounds on QbarV (called (a,b) in the write-up); to use MSE to fit the SuperLearner, let bd=NULL
## num.folds : number of folds for CV **** (SL or CV-TMLE?)
## convex.loss : convex surrogate loss function to use for classification
## wgts : Weights to be passed to SuperLearner for the multiple time point case
## OUTPUT
## A joint SL object

make.jointSL = function(W,A,Y,QbV.SL.library,cl.SL.library,Qbar0,g,V,alpha.risk.est,bd,num.folds=10,convex.loss='log',wgts=NULL){
	require('SuperLearner')
	num.alg = length(QbV.SL.library) + length(cl.SL.library)
	if(num.alg==0) stop('Must specify either QbV or classification SuperLearner library.')
	folds = (sample(1:nrow(W)) - 1) %% num.folds + 1

	if(is.null(QbV.SL.library)){
		QbVSL.obj = NULL
	} else {
		if(is.null(bd)){
			QbVSL.obj = list(SL.fit=qbarV.mse(W,A,Y,QbV.SL.library,V=V,Qbar0=Qbar0,newX=NULL,folds=folds,wgts=wgts),bd=bd)
		} else if (bd[1]<bd[2]) {
			QbVSL.obj = list(SL.fit=qbarV.logit(W,A,Y,QbV.SL.library,V=V,bd=bd,Qbar0=Qbar0,newX=NULL,folds=folds,wgts=wgts),bd=bd)
		} else {stop('Invalid entry for bd')}
		class(QbVSL.obj) = 'QbVSL'
	}

	if(is.null(cl.SL.library)){
		clSL.obj = NULL
	} else {
		clSL.obj = cl.SL(W,A,Y,cl.SL.library,V,Qbar=Qbar0,g=g,folds=folds,mtp.wgts=wgts)
	}

	jointSL.obj = list(SL.QbV=QbVSL.obj,SL.cl=clSL.obj,alpha=rep(1/num.alg,num.alg),V=V)
	class(jointSL.obj) = 'jointSL'

	if(num.alg>1){
		jointSL.obj$alpha = qbarV.alpha.cl(W,A,Y,jointSL.obj,Qbar0,g,bound.Qb0=TRUE,convex.loss=convex.loss,risk.est.method=list(method='emp.risk'),wgts=wgts)

		if(alpha.risk.est=='tmle' | alpha.risk.est=='dripcw'){
			jointSL.obj$alpha = qbarV.alpha2(W,A,Y,jointSL.obj,V,find.init=0,grid.size=30,Qbar0=Qbar0,unbd=bd,risk.est.method=list(method=alpha.risk.est),wgts=wgts)
		}
	} else {
		jointSL.obj$alpha = 1
		if(length(QbV.SL.library)==1){ names(jointSL.obj$alpha) = paste('QbV_',names(jointSL.obj$SL.QbV$SL.fit$coef),sep='')
		} else names(jointSL.obj$alpha) = paste('cl_',names(jointSL.obj$SL.cl$SL.fit$coef),sep='')
	}

	return(jointSL.obj)
}


## predict.jointSL
## Prediction function for a jointSL object
## INPUTS
## jointSL.obj : a jointSL object (a list containing SL.QbV, SL.cl, alpha)
	## SL.QbV : SuperLearner fit that takes A, W1, and W2 as predictors
	## SL.cl : SuperLearner fit that takes A, W1, and W2 as predictors and returns a probability estimate
	##	**** NOTE: in principle, doesn't need to be a probability, but any function f such that f(w)>1/2 returns
	##		a prediction where d=1 and d=0 otherwise
	## alpha : vector to combine prediction algoirthms from SL.QbV and SL.cl
	## V : the subset of covariates on which treatment decisions are to be made
## newdata : a W data frame on which to make predictions
## OUTPUT
## a list containing:
## preds : predictions for individual algorithms
## library.predict : matrix containing predictions for individual algorithms

predict.jointSL = function(jointSL.obj,newdata){
	V = jointSL.obj$V
	newdata = subset(newdata,select=V)
	library.predict = matrix(nrow=nrow(newdata),ncol=0)
	if(!is.null(jointSL.obj$SL.QbV)){
		QbV.preds = predict.QbVSL(jointSL.obj$SL.QbV,newdata=newdata)$library.predict
		colnames(QbV.preds) = paste('QbV_',colnames(QbV.preds),sep='')
		library.predict = cbind(library.predict,QbV.preds)
	}
	if(!is.null(jointSL.obj$SL.cl)){
		cl.preds = predict.clSL(jointSL.obj$SL.cl,newdata=newdata)$library.predict
		colnames(cl.preds) = paste('cl_',colnames(cl.preds),sep='')
		library.predict = cbind(library.predict,cl.preds)
	}
	pred.out = c(library.predict%*%cbind(jointSL.obj$alpha))
	return(list(pred=pred.out,library.predict=library.predict,d.pred=as.numeric(pred.out>0)))}


## predict.QbVSL
## Prediction function for a QbVSL object
## INPUTS
## QbVSL.obj : a QbVSL object (a list containing SL.fit and bd)
## newdata : a W data frame on which to make predictions
## OUTPUT
## a list containing:
## preds : predictions for individual algorithms
## library.predict : matrix containing predictions for individual algorithms

predict.QbVSL = function(QbVSL.obj,newdata){
	if(is.null(QbVSL.obj$bd)) { return(predict.SuperLearner(QbVSL.obj$SL.fit,newdata))
	} else {
		library.predict = predict.SuperLearner(QbVSL.obj$SL.fit,newdata=newdata)$library.predict*(QbVSL.obj$bd[2]-QbVSL.obj$bd[1]) + QbVSL.obj$bd[1]
		return(list(pred=c(library.predict%*%cbind(QbVSL.obj$SL.fit$coef)),library.predict=library.predict)) }
}


## predict.clSL
## Prediction function for a clSL object
## INPUTS
## clSL.obj : a list containing SL.fit, a SuperLearner fit, but will be rescaled by subtracting 1/2 from predictions
## OUTPUT
## a list containing:
## preds : predictions for individual algorithms
## library.predict : matrix containing predictions for individual algorithms

predict.clSL = function(clSL.obj,newdata){
	out = predict.SuperLearner(clSL.obj$SL.fit,newdata)
	out$pred = c(out$pred-1/2)
	out$library.predict = out$library.predict-1/2
	return(out)
}


## qbarV.alpha
## NOTE : THE NON-CV ALGORITHMS ARE ALREADY DOING WHAT THEIR CV COUNTERPARTS SHOULD BE DOING.
##		THE DISTINCTION IN THIS CODE ACTUALLY MAKES NO DIFFERENCE -- 4/2/2014
##		LEAVING CODE BECAUSE CAN MODIFY FOR WHEN Q/G ESTIMATED IF NEEDED
## Directly targets the choice of alpha on the simplex to target the resulting QbarV at EYd.
## INPUTS
## W : data frame of covariates
## A : vector of treatments
## Y : vector of outcomes
## SL.fit : SuperLearner fit that takes A, W1, and W2 as predictors
## V : a vector of column names in ObsData upon which treatment decisions will be made
## bound.Qb0 : Does Qbar0(A,W) fall between 0 and 1?
## nloptr.method : Method to use to optimize over the simplex with nloptr. If want
##	to optimize using Monte Carlo draws, make NULL. Probably want either NLOPT_LN_SBPLX or NLOPT_LN_COBYLA
## find.init : if nloptr.method is not NULL and find.init>0, then uses the best find.init
##	points drawn at random from the simplex (out of grid.size points)
## grid.size : size of the grid on which to do the grid search
## Qbar0 : a function to generate the mean outcome in a covariate-treatment strata
## g : treatment probability in an RCT
## unbd : Were the variables transformed to respect a bound? NULL if not, vector
##	containing bound limits otherwise
## risk.est.method : How should the risk be estimated? A list containing a string for the method:
##	empirical risk ('emp.risk'), CV-TMLE ('cv.tmle'), CV empirical risk ('cv.emp.risk'), or TMLE ('tmle').
##	CV-TMLE/CV empirical risk also require: a function of (W,A,Y,newX) for
##	estimating Qbar via SuperLearner, to be put in the list as CV.SL.fun (should have same SL library as SL.fit!);
##	a number of folds, to be put into the list as num.folds
## OUTPUT
## 

qbarV.alpha = function(W,A,Y,SL.fit,V,bound.Qb0=TRUE,nloptr.method="NLOPT_LN_SBPLX",find.init=20,grid.size=2e3,Qbar0=Qbar0.fun,g=0.5,unbd=NULL,risk.est.method=list(method='emp.risk')){
	lib.errors = (SL.fit$errorsInLibrary | SL.fit$errorsInCVLibrary)
	cumsum.lib.errors = cumsum(lib.errors)
	num.alg = ncol(SL.fit$Z)-sum(lib.errors)
	if(num.alg==0){stop('All algorithms returned errors.')}

	# Function which adds zeros to coefficients that had errors
	add.errors = function(alph){
		sapply(1:ncol(SL.fit$Z),function(i){
		if(lib.errors[i]){
			return(0)
		} else {
			return(alph[i-cumsum.lib.errors[i]])}})}

	# create SL fits on training samples if using CV-TMLE or CV empirical risk
	if(risk.est.method$method=='cv.tmle' | risk.est.method$method=='cv.emp.risk'){
		num.samples = nrow(W)
		if(risk.est.method$num.folds>num.samples){
			warning('More folds than data points. Leave-one-out cross-validation will be used.')
			risk.est.method$num.folds = num.samples}
		fold.remainder = num.samples%%risk.est.method$num.folds
		if(fold.remainder==0){
			folds = rep(1:risk.est.method$num.folds,each=num.samples/risk.est.method$num.folds)
		} else {
			folds = c(rep(1:fold.remainder,each=ceiling(num.samples/risk.est.method$num.folds)),
				rep((fold.remainder+1):risk.est.method$num.folds,each=floor(num.samples/risk.est.method$num.folds)))}

		Q.Z = lapply(1:risk.est.method$num.folds,function(i){
			tmp.SL = risk.est.method$CV.SL.fun(W[folds!=i,],A[folds!=i],Y[folds!=i],newX=subset(W[folds==i,],select=V))
			if(sum(tmp.SL$SL.library$library$predAlgorithm!=SL.fit$SL.library$library$predAlgorithm)>0) stop('SuperLearner fit in CV-TMLE should use the same library as the SuperLearner fit in SL.fit.')
			return(tmp.SL$library.predict[,which(!SL.fit$errorsInLibrary)])})
		Q.alpha = function(alph){
			do.call(rbind,lapply(Q.Z,function(Qf){
					return(Qf%*%alph)}))}
	}

	# evaluate the loss function at all points in a grid
	eval.grid = function(a.grid){
		if(risk.est.method$method=='emp.risk' | risk.est.method$method=='cv.emp.risk'){
			nrow.grid = nrow(a.grid)
			if(risk.est.method$method=='cv.emp.risk'){
				Q.est = c(Q.alpha(t(a.grid)))
			} else {
				Q.est = c(SL.fit$Z[,which(!SL.fit$errorsInLibrary)]%*%t(a.grid))
			}
			if(!is.null(unbd)){
				Q.est = Q.est*(unbd[2]-unbd[1]) + unbd[1]
			}

			W.rep = W[rep(1:nrow(W),nrow.grid),]
			A.rep = rep(A,nrow.grid)
			Y.rep = rep(Y,nrow.grid)

			d.opt = (Q.est>0)
			Qbopt = Qbar0(d.opt,W.rep)
			QbA = rep(Qbar0(A,W),nrow.grid)
			Qb0 = rep(Qbar0(0,W),nrow.grid)
			Qb1 = rep(Qbar0(1,W),nrow.grid)

			risk.est = colMeans(matrix(-((A.rep==d.opt)/g * (Y.rep-QbA) + Qbopt),ncol=nrow.grid))
		} else if(risk.est.method$method=='cv.tmle' | risk.est.method$method=='tmle') {
			risk.est = -apply(a.grid,1,function(alph){
				if(risk.est.method$method=='cv.tmle'){
					Q.est.curr = c(Q.alpha(cbind(alph)))
				} else {
					Q.est.curr = SL.fit$Z[,which(!SL.fit$errorsInLibrary)]%*%cbind(alph)
				}
				if(!is.null(unbd)){
					Q.est.curr = Q.est.curr*(unbd[2]-unbd[1]) + unbd[1] }
				d.opt.curr = (Q.est.curr>0)
				g.vec = rep(g,length(Q.est.curr)) ## **** note : if g becomes a function, this needs to change
				g.vec[A==0] = 1-g.vec[A==0]
				Qb0.dopt = Qbar0(d.opt.curr,W)
				return(d.tmle(A,Y,Qb0.dopt,g.vec,d.opt.curr,bound.Qb0=TRUE)$psi)})
		} else{
			stop('Invalid choice of risk.est.method. Currently only empirical risk and CV-TMLE risk estimates supported.')
		}
		return(list(risk.est=risk.est,sols=a.grid))}

	if(!is.null(nloptr.method) & (find.init>0)){
		if(find.init>grid.size){
			warning('find.init should be at least as large as grid.size. Setting find.init=grid.size.')
			find.init = grid.size}
		init.out = eval.grid(rbind(simplex.sample(num.alg,grid.size)$samples,SL.fit$coef))
		simplex.grid = init.out$sols[order(init.out$risk.est)[1:find.init],]
	} else {
		simplex.grid = rbind(simplex.sample(num.alg,grid.size)$samples,SL.fit$coef)
	}

	## optimize if desired
	if(!is.null(nloptr.method)){
		f = function(x){
			eval.grid(rbind(x))$risk.est}	
		nloptr.out=apply(simplex.grid,1,function(i){
			nloptr(x0=i,
			eval_f=f,
			lb=rep(0,num.alg),
			opts=list("algorithm"=nloptr.method,
				"ftol_rel"=1.0e-7,
				"maxeval"=5e4))})
		# sols = t(sapply(nloptr.out,function(xx){xx$solution}))
		sols = t(sapply(nloptr.out,function(xx){xx$solution/sum(xx$solution)}))
		risk.est = sapply(nloptr.out,function(xx){xx$objective})
		nonopt.grid = diag(num.alg)
	} else {
		risk.est = NULL
		sols = NULL
		nonopt.grid = rbind(simplex.grid,diag(num.alg))
	}

	nog = eval.grid(nonopt.grid)
	risk.est = c(risk.est,nog$risk.est)
	sols = rbind(sols,nog$sols)
	alpha = sols[which.min(risk.est),]
	names(alpha) = names(SL.fit$coef)
	# give algorithms with errors alpha of zero
	alpha.w.errors = add.errors(alpha)

	SL.fit$coef = alpha.w.errors
	SL.fit$SL.predict = SL.fit$Z%*%cbind(SL.fit$coef)
	return(SL.fit)
}


## qbarV.alpha2
## Similar to qbarV.alpha but takes jointSL fit and weights (weights for multiple time point case)
## replaced emp.risk with dripcw

qbarV.alpha2 = function(W,A,Y,jointSL.fit,V,bound.Qb0=TRUE,nloptr.method="NLOPT_LN_SBPLX",find.init=20,grid.size=2e3,Qbar0=Qbar0.fun,g=0.5,unbd=NULL,risk.est.method=list(method='dripcw'),wgts=NULL){
	Z = NULL
	if(!is.null(jointSL.fit$SL.QbV)) {
		if(is.null(jointSL.fit$SL.QbV$bd)) { QbV.Z = jointSL.fit$SL.QbV$SL.fit$Z
		} else QbV.Z = jointSL.fit$SL.QbV$SL.fit$Z*(jointSL.fit$SL.QbV$bd[2]-jointSL.fit$SL.QbV$bd[1]) + jointSL.fit$SL.QbV$bd[1]
		colnames(QbV.Z) = paste('QbV_',names(jointSL.fit$SL.QbV$SL.fit$coef),sep='')
		Z = cbind(Z,QbV.Z) }
	if(!is.null(jointSL.fit$SL.cl)){
		cl.Z = jointSL.fit$SL.cl$SL.fit$Z-1/2
		colnames(cl.Z) = paste('cl_',names(jointSL.fit$SL.cl$SL.fit$coef),sep='')
		Z = cbind(Z,cl.Z) }
	# predict(jointSL.fit,subset(W,select=V))$library.predict
	num.alg = ncol(Z)
	if(num.alg==0){stop('All algorithms returned errors.')}
	if(is.null(wgts)) wgts = rep(1,nrow(Z))

	# evaluate the loss function at all points in a grid
	eval.grid = function(a.grid){
		if(risk.est.method$method=='dripcw'){
			nrow.grid = nrow(a.grid)
			Q.est = c(Z%*%t(a.grid))

			if(!is.null(unbd)){
				Q.est = Q.est*(unbd[2]-unbd[1]) + unbd[1]
			}

			W.rep = W[rep(1:nrow(W),nrow.grid),]
			A.rep = rep(A,nrow.grid)
			Y.rep = rep(Y,nrow.grid)
			wgts.rep = rep(wgts,nrow.grid)

			d.opt = (Q.est>0)
			Qbopt = Qbar0(d.opt,W.rep)
			QbA = rep(Qbar0(A,W),nrow.grid)
			Qb0 = rep(Qbar0(0,W),nrow.grid)
			Qb1 = rep(Qbar0(1,W),nrow.grid)

			risk.est = colMeans(matrix(-(wgts.rep*((A.rep==d.opt)/g * (Y.rep-QbA) + Qbopt)),ncol=nrow.grid))
		} else if(risk.est.method$method=='tmle') {
			risk.est = -apply(a.grid,1,function(alph){
				Q.est.curr = Z%*%cbind(alph)
				if(!is.null(unbd)){
					Q.est.curr = Q.est.curr*(unbd[2]-unbd[1]) + unbd[1] }
				d.opt.curr = as.numeric(Q.est.curr>0)
				g.vec = rep(g,length(Q.est.curr)) ## **** note : if g becomes a function, this needs to change
				g.vec[A==0] = 1-g.vec[A==0]
				Qb0.dopt = Qbar0(d.opt.curr,W)
				return(d.tmle(A,Y,Qb0.dopt,g.vec,d.opt.curr,bound.Qb0=TRUE,wgts=wgts)$psi)})
		} else{
			stop('Invalid choice of risk.est.method. Currently only empirical risk and CV-TMLE risk estimates supported.')
		}
		return(list(risk.est=risk.est,sols=a.grid))}

	if(!is.null(nloptr.method) & (find.init>0)){
		if(find.init>grid.size){
			warning('find.init should be at least as large as grid.size. Setting find.init=grid.size.')
			find.init = grid.size}
		init.out = eval.grid(rbind(simplex.sample(num.alg,grid.size)$samples,jointSL.fit$alpha))
		simplex.grid = init.out$sols[order(init.out$risk.est)[1:find.init],]
	} else {
		simplex.grid = rbind(simplex.sample(num.alg,grid.size)$samples,jointSL.fit$alpha)
	}

	## optimize if desired
	if(!is.null(nloptr.method)){
		f = function(x){
			eval.grid(rbind(x))$risk.est}	
		nloptr.out=apply(simplex.grid,1,function(i){
			nloptr(x0=pmax(i,0),
			eval_f=f,
			lb=rep(0,num.alg),
			opts=list("algorithm"=nloptr.method,
				"ftol_rel"=1.0e-7,
				"maxeval"=5e4))})
		# sols = t(sapply(nloptr.out,function(xx){xx$solution}))
		sols = t(sapply(nloptr.out,function(xx){xx$solution/sum(xx$solution)}))
		risk.est = sapply(nloptr.out,function(xx){xx$objective})
		nonopt.grid = diag(num.alg)
	} else {
		risk.est = NULL
		sols = NULL
		nonopt.grid = rbind(simplex.grid,diag(num.alg))
	}

	nog = eval.grid(nonopt.grid)
	risk.est = c(risk.est,nog$risk.est)
	sols = rbind(sols,nog$sols)
	alpha = sols[which.min(risk.est),]
	names(alpha) = names(jointSL.fit$alpha)

	# jointSL.fit$alpha = alpha
	# jointSL.fit$SL.predict = Z%*%cbind(jointSL.fit$alpha)
	return(alpha)
}



## EYd.tmle
## Given an algorithm for estimating Qbar(V), estimate EYd using the TMLE,
## CV-TMLE, IPTW, and G-comp estimators.
## INPUTS
## W : data frame of covariates
## A : vector of treatments
## Y : vector of outcomes
## SL.library: SuperLearner library for estimating Qbar(V)
## Qbar0 : a function to generate the mean outcome in a covariate-treatment strata. Default gives an IPTW estimator
## g : treatment probability in an RCT
## V : a vector of column names in ObsData upon which treatment decisions will be made
## alpha.risk.est : Options are: (i) nonnegative least-squares ('nnls'); (ii) empirical risk ('emp.risk'); (iii) CV-TMLE ('cv.tmle'); (iv) TMLE ('tmle')
## bd : bounds on QbarV (called (a,b) in the write-up); to use MSE to fit the SuperLearner, let bd=NULL
## num.folds : number of folds to use when running CV-TMLE (both for the alpha risk estimate and the CV-TMLE EYd estimate)
## OUTPUT
## ****

EYd.tmle = function(W,A,Y,SL.library,Qbar0=function(AA,WW){rep(1/2,length(AA))},g=0.5,V=c('W1','W2'),alpha.risk.est='nnls',bd=c(-1,1),num.folds=10,cv.ests=FALSE){
	num.samples = length(A)
	# Get initial SL fits
	if(is.null(bd)){
		SL.fun.init = function(WW,AA,YY,newX){qbarV.mse(WW,AA,YY,SL.library,V=V,Qbar0=Qbar0,newX=newX)}
	} else if (bd[1]<bd[2]) {
		SL.fun.init = function(WW,AA,YY,newX){qbarV.logit(WW,AA,YY,SL.library,V=V,bd=bd,Qbar0=Qbar0,newX=newX)}
	} else {stop('Invalid entry for bd')}

	# Update alpha in initial SL fits if necessary
	if(alpha.risk.est=='emp.risk' | alpha.risk.est=='tmle'){
		SL.fun = function(WW,AA,YY){qbarV.alpha(WW,AA,YY,SL.fun.init(WW,AA,YY,NULL),V,Qbar0=Qbar0,risk.est.method=list(method=alpha.risk.est),unbd=bd,find.init=0,grid.size=30)}
	} else if(alpha.risk.est=='cv.tmle' | alpha.risk.est=='cv.emp.risk'){
		SL.fun = function(WW,AA,YY,newX){qbarV.alpha(WW,AA,YY,SL.fun.init(WW,AA,YY,NULL),V,Qbar0=Qbar0,risk.est.method=list(method=alpha.risk.est,num.folds=num.folds,CV.SL.fun=SL.fun.init),unbd=bd,find.init=0,grid.size=30)}
	} else if (alpha.risk.est=='nnls'){
		SL.fun = function(WW,AA,YY){SL.fun.init(WW,AA,YY,NULL)}
	} else {stop('Invalid entry for alpha.risk.est.')}
	QbV.fit = SL.fun(W,A,Y)

	# Get pieces of likelihood needed to estimate EYd
	if(is.null(bd)){d.est = (QbV.fit$SL.predict>0)
	} else {d.est = (QbV.fit$SL.predict*(bd[2]-bd[1])+bd[1]>0)}
	Qbd = Qbar0(d.est,W)
	g.vec = rep(g,num.samples)

	# Estimate EYd
	ic = (A==d.est)/g * (Y-Qbd) + Qbd
	dripcw.est = mean(ic)
	tmle.out = d.tmle(A,Y,Qbd,g,d.est)
	tmle.est = tmle.out$psi
	EYd.est = c(dripcw.est,tmle.est)
	names(EYd.est) = c('dripcw','tmle')
	ci.widths = c(qnorm(0.975)*sd(ic)/sqrt(num.samples),qnorm(0.975)*sqrt(tmle.out$var.est))
	names(ci.widths) = c('dr.ci.width','tmle.ci.width')

	if(cv.ests){
		# Get SL fits for CV-TMLE
		if(num.folds>num.samples){
			warning('More folds than data points. Leave-one-out cross-validation will be used.')
			num.folds = num.samples}
		fold.remainder = num.samples%%num.folds
		if(fold.remainder==0){	
			folds = rep(1:num.folds,each=num.samples/num.folds)
		} else {
			folds = c(rep(1:fold.remainder,each=ceiling(num.samples/num.folds)),
				rep((fold.remainder+1):num.folds,each=floor(num.samples/num.folds)))}
		QbV.CV.vec = unlist(lapply(1:num.folds,function(i){
			tmp = SL.fun(W[folds!=i,],A[folds!=i],Y[folds!=i])
			return(c(predict.SuperLearner(tmp,newdata=W[folds==i,])$library.predict%*%cbind(tmp$coef)))}))
		if(is.null(bd)){d.CV = (QbV.CV.vec>0)
		} else {d.CV = (QbV.CV.vec*(bd[2]-bd[1])+bd[1]>0)}

		# Estimate CV EYd
		Qbd.CV = Qbar0(d.CV,W)
		cv.ic = (A==d.CV)/g * (Y-Qbd.CV) + Qbd.CV
		cvdripcw.est = mean(cv.ic)
		cvtmle.out = d.tmle(A,Y,Qbd.CV,g,d.CV)
		cvtmle.est = cvtmle.out$psi
		EYd.est = c(EYd.est,cvdripcw.est,cvtmle.est)
		names(EYd.est)[3:4] = c('cvdripcw','cvtmle')
 		ci.widths = c(ci.widths,qnorm(0.975)*sd(cv.ic)/sqrt(num.samples),qnorm(0.975)*sqrt(cvtmle.out$var.est))
 		names(ci.widths)[3:4] = c('cvdr.ci.width','cvtmle.ci.width')
	}

	return(list(EYd.est=EYd.est,ci.widths=ci.widths,QbV.fit=QbV.fit))
}


## EYd.tmle2
## Given an algorithm for estimating Qbar(V), estimate EYd using the TMLE,
## CV-TMLE, IPTW, and G-comp estimators when classification is included in SL libraries.
## INPUTS
## W : data frame of covariates
## A : vector of treatments
## Y : vector of outcomes
## QbV.SL.library : SuperLearner library for estimating Qbar(V). Can be NULL
## cl.SL.library : SuperLearner library for classification. Can be NULL
## jointSL.list : If specified, overrides the two SuperLearner library arguments -- mostly used for simulation purposes
##	A list containing named elements titled:
##		all : jointSL object fit on entire data set
##		cv : a list of length num.folds, specifying the cross-validated jointSL fits
##			elements need to correspond to folds as specified in function
## Qbar0 : a function to generate the mean outcome in a covariate-treatment strata. Default gives an IPTW estimator
## g : treatment probability in an RCT
## V : a vector of column names in ObsData upon which treatment decisions will be made
## alpha.risk.est : Options are: (i) nonnegative least-squares ('nnls'); (ii) empirical risk ('emp.risk'); (iii) CV-TMLE ('cv.tmle'); (iv) TMLE ('tmle')
## bd : bounds on QbarV (called (a,b) in the write-up); to use MSE to fit the SuperLearner, let bd=NULL
## num.folds : number of folds to use when running CV-TMLE (both for the alpha risk estimate and the CV-TMLE EYd estimate)
## OUTPUT
## ****

EYd.tmle2 = function(W,A,Y,QbV.SL.library,cl.SL.library,jointSL.list=NULL,Qbar0=function(AA,WW){rep(1/2,length(AA))},g=0.5,V=c('W1','W2'),alpha.risk.est='emp.risk',bd=c(-1,1),num.folds=10,cv.ests=TRUE){
	num.samples = length(A)
	# Get initial SL fits
	if(is.null(jointSL.list)){
		SL.fun = function(WW,AA,YY){make.jointSL(WW,AA,YY,QbV.SL.library,cl.SL.library,Qbar0=Qbar0,g=g,V=V,alpha.risk.est=alpha.risk.est,bd=bd,num.folds=10)}
		QbV.fit = SL.fun(W,A,Y)
		d.est = (predict(QbV.fit,newdata=subset(W,select=V))$pred>0)
	} else {
		QbV.fit = jointSL.list[[1]]
		d.est = (predict(jointSL.list[[1]],newdata=subset(W,select=V))$pred>0)
	}

	# Get pieces of likelihood needed to estimate EYd
	Qbd = Qbar0(d.est,W)
	g.vec = rep(g,num.samples)

	# Estimate EYd
	ic = (A==d.est)/g * (Y-Qbd) + Qbd
	dripcw.est = mean(ic)

	tmle.out = d.tmle(A,Y,Qbd,g,d.est)
	tmle.est = tmle.out$psi
	EYd.est = c(dripcw.est,tmle.est)
	names(EYd.est) = c('dripcw','tmle')
	ci.widths = c(qnorm(0.975)*sd(ic)/sqrt(num.samples),qnorm(0.975)*sqrt(tmle.out$var.est))
	names(ci.widths) = c('dr.ci.width','tmle.ci.width')

	# Get CV estimates if desired
	if(cv.ests){
		# Get SL fits for CV-TMLE
		if(num.folds>num.samples){
			warning('More folds than data points. Leave-one-out cross-validation will be used.')
			num.folds = num.samples}
		fold.remainder = num.samples%%num.folds
		if(fold.remainder==0){	
			folds = rep(1:num.folds,each=num.samples/num.folds)
		} else {
			folds = c(rep(1:fold.remainder,each=ceiling(num.samples/num.folds)),
				rep((fold.remainder+1):num.folds,each=floor(num.samples/num.folds)))}
		QbV.CV.vec = unlist(lapply(1:num.folds,function(i){
			if(is.null(jointSL.list)){
				tmp = SL.fun(W[folds!=i,],A[folds!=i],Y[folds!=i])
			} else {
				tmp = jointSL.list[[2]][[i]]
			}
			return(c(predict(tmp,newdata=subset(W[folds==i,],select=V))$library.predict%*%cbind(tmp$alpha)))}))
		d.CV = (QbV.CV.vec>0)

		Qbd.CV = Qbar0(d.CV,W)
		cv.ic = (A==d.CV)/g * (Y-Qbd.CV) + Qbd.CV
		cvdripcw.est = mean(cv.ic)
		cvtmle.out = d.tmle(A,Y,Qbd.CV,g,d.CV)
		cvtmle.est = cvtmle.out$psi
		EYd.est = c(EYd.est,cvdripcw.est,cvtmle.est)
		names(EYd.est)[3:4] = c('cvdripcw','cvtmle')
 		ci.widths = c(ci.widths,qnorm(0.975)*sd(cv.ic)/sqrt(num.samples),qnorm(0.975)*sqrt(cvtmle.out$var.est))
 		names(ci.widths)[3:4] = c('cvdr.ci.width','cvtmle.ci.width')
	}

	return(list(EYd.est=EYd.est,ci.widths=ci.widths,QbV.fit=QbV.fit))
}


## eval.performance
## Evaluates the performance of each SL fit according to the logistic, squared error, and EYa loss functions using Monte Carlo simulations.
## **** add the other true failures as we go
## INPUTS
## tx.alg : SuperLearner fit that takes A, W1, and W2 as predictors, jointSL object, param.blip object, 
##	or 0 indicating no treatment for anyone or 1 indicating treatment for everyone
## n.test : Number of values to simulate from the true observed data distribution
## Qbar0 : a function to generate the mean outcome in a covariate-treatment strata
## g : treatment probability in an RCT
## V : variables upon which treatment decisions will be based
## unbd : Were the variables transformed to respect a bound? NULL if not, vector
##	containing bound limits otherwise
## W.fun : a function to generate W's which takes n (sample size) as input
## OUTPUT
## A vector containing the estimated risk from each of the three outputs.

eval.performance = function(tx.alg,n.test=1e6,Qbar0=Qbar0.fun,g=0.5,V=c('W1','W2'),unbd=NULL,W.fun=rnorm){
	all.data = sim.data(n.test,W.fun=W.fun)
	test.set = all.data$ObsData
	W1 = test.set$W1
	W2 = test.set$W2
	W3 = test.set$W3
	W4 = test.set$W4
	W=data.frame(W1=W1,W2=W2,W3=W3,W4=W4)
	A = test.set$A
	Y = test.set$Y
	Qb0 = Qbar0(0,W)
	Qb1 = Qbar0(1,W)
	QbA = Qbar0(A,W)
	D1 = (2*A-1)/g * (Y-QbA) + Qb1 - Qb0
	D1.ipcw = (2*A-1)/g * Y
	if(is.numeric(tx.alg)){
		if(tx.alg==0 | tx.alg==1){
			mse.true = NA
			mse.ipcw = NA
			d.est = rep(tx.alg,n.test)
		} else {
			stop('Numeric input for tx.alg must be 0 or 1.')}
	} else if(is(tx.alg,'SuperLearner')) {
		SL.preds = predict.SuperLearner(tx.alg,subset(test.set,select=V),onlySL=TRUE)$pred
		if(!is.null(unbd)) SL.preds = SL.preds*(unbd[2]-unbd[1]) + unbd[1]
		## Mean-squared error
		mse.true = mean((SL.preds - D1)^2)
		mse.ipcw = mean((SL.preds - D1.ipcw)^2)
		d.est = as.numeric(SL.preds>0)
	} else if(is(tx.alg,'param.blip')) {
		pb.preds = predict(tx.alg,newdata=W)
		mse.true = mean((pb.preds - D1)^2)
		mse.ipcw = mean((pb.preds - D1.ipcw)^2)
		d.est = as.numeric(pb.preds>0)
	} else if(is(tx.alg,'jointSL')) {
		SL.preds = predict(tx.alg,subset(test.set,select=V))$pred
		## Mean-squared error
		mse.true = mean((SL.preds - D1)^2)
		mse.ipcw = mean((SL.preds - D1.ipcw)^2)
		d.est = as.numeric(SL.preds>0)
	} else {
		stop('tx.alg must either be a SuperLearner fit, param.blip object, jointSL, or a numeric 0 or 1.')}

	## EYd
	Q.d = Qbar0(d.est,W)
	EYd = mean(Q.d)

	out.vec = c(mse.true=mse.true,mse.ipcw=mse.ipcw,EYd=EYd,mean.d=mean(d.est))
	names(out.vec) = c('mse.true','mse.ipcw','EYd','mean.d')
	return(out.vec)
}


## EYd.est.eval
## Wrapper function for EYd.tmle and eval.performance.
## INPUTS
## See EYd.tmle, with the distinction that:
## Qbar0 : an estimated function to generate the mean outcome in a covariate-treatment strata
## Qbar0.true : the true Qbar0 function
## OUTPUT
## 

EYd.est.eval = function(W,A,Y,SL.library,Qbar0=function(AA,WW){rep(1/2,length(AA))},Qbar0.true=Qbar0.fun,g=0.5,V=c('W1','W2'),alpha.risk.est='nnls',bd=c(-1,1),num.folds=10,W.fun=rnorm){
	EYd.out = EYd.tmle(W,A,Y,SL.library,Qbar0=Qbar0,g=g,V=V,alpha.risk.est=alpha.risk.est,bd=bd,num.folds=num.folds)
	perf = tryCatch(colMeans(do.call(rbind,lapply(1:10,function(i){eval.performance(EYd.out$QbV.fit,n.test=1e5,Qbar0=Qbar0.true,g=g,unbd=bd,V=V,W.fun=W.fun)}))),
		error=function(e){NULL})
	out.list = list(EYd.est=EYd.out$EYd.est,ci.width=EYd.out$ci.width,cv.ci.width=EYd.out$cv.ci.width,perf=perf,alpha=EYd.out$QbV.fit$coef)
	return(out.list)
}


## EYd.est.eval2
## Wrapper function for EYd.tmle and eval.performance when classification functions are included.
## INPUTS
## See EYd.tmle, with the distinction that:
## Qbar0 : an estimated function to generate the mean outcome in a covariate-treatment strata
## Qbar0.true : the true Qbar0 function
## just.perf : Just do a performance test, or also check inference?
## alpha.methods : what methods should be tried for fitting alpha? Valid options are:
##	'QbV_nnls', 'QbV_log', 'cl_log', and 'Qbv.cl_log'
## OUTPUT
## 
	
EYd.est.eval2 = function(W,A,Y,QbV.SL.library,cl.SL.library,Qbar0=function(AA,WW){rep(1/2,length(AA))},Qbar0.true=Qbar0.fun,g=0.5,V=c('W1','W2'),bd=c(-1,1),num.folds=10,W.fun=rnorm,convex.loss='log',just.perf=FALSE,eval.methods=c('QbV.cl_log')){
	
	num.samples = nrow(W)
	fold.remainder = num.samples%%num.folds
	if(fold.remainder==0){	
		folds = rep(1:num.folds,each=num.samples/num.folds)
	} else {
		folds = c(rep(1:fold.remainder,each=ceiling(num.samples/num.folds)),
			rep((fold.remainder+1):num.folds,each=floor(num.samples/num.folds)))}

	if(!just.perf){
		jointSL.CV = lapply(1:num.folds,function(i){
			make.jointSL(W[folds!=i,],A[folds!=i],Y[folds!=i],QbV.SL.library,cl.SL.library,Qbar0=Qbar0,g=g,V=V,alpha.risk.est='emp.risk',bd=bd,num.folds=num.folds,convex.loss=convex.loss)})
	} else {
		jointSL.CV = NULL
	}
	jointSL.all = make.jointSL(W,A,Y,QbV.SL.library,cl.SL.library,Qbar0=Qbar0,g=g,V=V,alpha.risk.est='emp.risk',bd=bd,num.folds=num.folds,convex.loss=convex.loss)

	jointSL.list = list(all=jointSL.all,cv=jointSL.CV)
	jsl.lol = list() # list of lists containing the jointSL list for each coefficient strategy
	ind = 1

	# QbV_nnls
	if('QbV_nnls'%in%eval.methods){
		curr.jsl = jointSL.list
		curr.jsl[[1]]$SL.cl = NULL
		curr.jsl[[1]]$alpha = curr.jsl[[1]]$SL.QbV$SL.fit$coef
		if(!just.perf){
			curr.jsl[[2]] = lapply(1:num.folds,function(i){
				x = curr.jsl[[2]][[i]]
				x$SL.cl = NULL
				x$alpha = x$SL.QbV$SL.fit$coef
				return(x)})
		}
		jsl.lol[[ind]] = curr.jsl
		names(jsl.lol)[ind] = 'QbV_nnls'
		ind = ind+1
	}

	# QbV_log
	if('QbV_log'%in%eval.methods){
		curr.jsl = jointSL.list
		curr.jsl[[1]]$SL.cl = NULL
		curr.jsl[[1]]$alpha = qbarV.alpha.cl(W,A,Y,curr.jsl[[1]],Qbar0,g,bound.Qb0=TRUE,convex.loss=convex.loss,risk.est.method=list(method='emp.risk'))
		if(!just.perf){
			curr.jsl[[2]] = lapply(1:num.folds,function(i){
				x = curr.jsl[[2]][[i]]
				x$SL.cl = NULL
				x$alpha = qbarV.alpha.cl(W[folds!=i,],A[folds!=i],Y[folds!=i],x,Qbar0,g,bound.Qb0=TRUE,convex.loss=convex.loss,risk.est.method=list(method='emp.risk'))
				return(x)})
		}
		jsl.lol[[ind]] = curr.jsl
		names(jsl.lol)[ind] = paste('Qbv_',convex.loss,sep='')
		ind = ind+1
	}

	# # QbV_tmle
	if('Qbv_tmle'%in%eval.methods){
		curr.jsl = jointSL.list
		curr.jsl[[1]]$SL.cl = NULL
		curr.jsl[[1]]$alpha = qbarV.alpha(W,A,Y,curr.jsl[[1]]$SL.QbV$SL.fit,V,Qbar0=Qbar0,risk.est.method=list(method='tmle'),unbd=bd,find.init=0,grid.size=30)$coef
		curr.jsl[[2]] = lapply(1:num.folds,function(i){
			x = curr.jsl[[2]][[i]]
			x$SL.cl = NULL
			x$alpha = qbarV.alpha(W[folds!=i,],A[folds!=i],Y[folds!=i],x$SL.QbV$SL.fit,V,Qbar0=Qbar0,risk.est.method=list(method='tmle'),unbd=bd,find.init=0,grid.size=30)$coef
			return(x)})
		jsl.lol[[ind]] = curr.jsl
		names(jsl.lol)[ind] = 'QbV_tmle'
		ind = ind+1
	}

	# cl_log
	if('cl_log'%in%eval.methods){
		curr.jsl = jointSL.list
		curr.jsl[[1]]$SL.QbV = NULL
		curr.jsl[[1]]$alpha = qbarV.alpha.cl(W,A,Y,curr.jsl[[1]],Qbar0,g,bound.Qb0=TRUE,convex.loss=convex.loss,risk.est.method=list(method='emp.risk'))
		if(!just.perf){
			curr.jsl[[2]] = lapply(1:num.folds,function(i){
				x = curr.jsl[[2]][[i]]
				x$SL.QbV = NULL
				x$alpha = qbarV.alpha.cl(W[folds!=i,],A[folds!=i],Y[folds!=i],x,Qbar0,g,bound.Qb0=TRUE,convex.loss=convex.loss,risk.est.method=list(method='emp.risk'))
				return(x)})
		}
		jsl.lol[[ind]] = curr.jsl
		names(jsl.lol)[ind] = paste('cl_',convex.loss,sep='')
		ind = ind+1
	}

	## QbV.cl_log
	if('QbV.cl_log'%in%eval.methods){
		jsl.lol[[ind]] = jointSL.list
		names(jsl.lol)[ind] = paste('Qbv.cl_',convex.loss,sep='')
		ind = ind + 1
	}

	if(ind==1){
		stop('Invalid eval.methods specification.')
	}

	out.list = lapply(jsl.lol,function(jointSL.curr){
		if(!just.perf){
			EYd.out = EYd.tmle2(W,A,Y,QbV.SL.library=NULL,cl.SL.library=NULL,jointSL.list=jointSL.curr,Qbar0=Qbar0,g=g,V=V,alpha.risk.est=NULL,bd=bd,num.folds=num.folds)
			cv.perf = lapply(jointSL.curr[[2]],function(cv.QbV.fit){
				tryCatch(colMeans(do.call(rbind,lapply(1:10,function(i){eval.performance(cv.QbV.fit,n.test=1e5,Qbar0=Qbar0.true,g=g,unbd=bd,V=V,W.fun=W.fun)}))),error=function(e){NULL})})
		} else {
			EYd.out = list(QbV.fit = jointSL.curr[[1]],EYd.est=NULL,ci.width=NULL)
			cv.perf = NULL
		}
		perf = tryCatch(colMeans(do.call(rbind,lapply(1:10,function(i){eval.performance(EYd.out$QbV.fit,n.test=1e5,Qbar0=Qbar0.true,g=g,unbd=bd,V=V,W.fun=W.fun)}))),error=function(e){NULL})
		return(list(EYd.est=EYd.out$EYd.est,ci.width=EYd.out$ci.width,perf=perf,cv.perf=cv.perf,alpha=EYd.out$QbV.fit$alpha))
	})
	names(out.list) = names(jsl.lol)
	return(out.list)
}


## EYd.mc
## Returns an estimate of the true EYd for a given simulation, using W3 as V.

EYd.mc = function(n.mc=1e5,n.eval=1e5,Qbar0=Qbar0.fun,g=0.5,W.fun=rnorm){
	Obs = sim.data(n.mc,Qbar0=Qbar0,g=0.5,binom=TRUE,W.fun=W.fun)$ObsData
	W1.mc = Obs$W1
	W2.mc = Obs$W2
	W4.mc = Obs$W4
	eval.data = sim.data(n.eval,Qbar0=Qbar0,g=0.5,binom=TRUE,W.fun=W.fun)$FullData
	d.opt = sapply(1:n.eval,function(i){
		curr.W3 = eval.data$W3[i]
		return(mean(Qbar0(1,cbind(W1.mc,W2.mc,curr.W3,W4.mc))-Qbar0(0,cbind(W1.mc,W2.mc,curr.W3,W4.mc)))>0)})
	# Y.opt = rep(NA,n.eval)
	# Y.opt[d.opt==0] = eval.data$Y0[d.opt==0]
	# Y.opt[d.opt==1] = eval.data$Y1[d.opt==1]
	# return(mean(Y.opt))
	Q.d = Qbar0(d.opt,subset(eval.data,select=c(W1,W2,W3,W4)))
	return(mean(Q.d))
}


## EYd.mc2
## Returns an estimate of the true EYd for a given simulation, using W1, W2, W3, W4 as V.

EYd.mc2 = function(n.eval=1e5,Qbar0=Qbar0.fun,g=0.5,W.fun=rnorm){
	eval.data = sim.data(n.eval,Qbar0=Qbar0,g=0.5,binom=TRUE,W.fun=W.fun)$FullData
	d.opt = (Qbar0(1,cbind(eval.data$W1,eval.data$W2,eval.data$W3,eval.data$W4))-Qbar0(0,cbind(eval.data$W1,eval.data$W2,eval.data$W3,eval.data$W4))>0)
	Q.d = Qbar0(d.opt,subset(eval.data,select=c(W1,W2,W3,W4)))
	return(mean(Q.d))
}


## param.blip
## Fits a parametric blip function for the given formula, link function (implied by bd), and choice of risk function.
## INPUTS
## W : data frame of covariates
## A : vector of treatments
## Y : vector of outcomes
## form : a formula or an object that can be coerced to a formula object, with no response and
##	(a subset of the) variables in W as the predictors
## risk.est.method : either 'QbV' (simply try to fit Qbar(V) directly), 'emp.risk', or 'tmle'
##	**** should implement cv.tmle as well eventually
## nloptr.method : Method to use to optimize over the simplex with nloptr. If want
##	to optimize using Monte Carlo draws, make NULL. Probably want either NLOPT_LN_SBPLX or NLOPT_LN_COBYLA
## num.init : number of initial points to use when searching for the optimal beta for risk.est.method 'emp.risk' or 'tmle'
## bd : bounds on QbarV (called (a,b) in the write-up). If NULL then specifies Gaussian family with identity link,
##	otherwise binomial with logit link
## Qbar0 : a function to generate the mean outcome in a covariate-treatment strata
## g : treatment probability in an RCT
## OUTPUT
## Coefficients to be used in model.

fit.param.blip = function(W,A,Y,form,risk.est.method='QbV',nloptr.method="NLOPT_LN_SBPLX",num.init=40,bd=NULL,Qbar0=Qbar0.fun,g=0.5){
	if(!is(form,'formula')) form = tryCatch(as.formula(form),error=function(e){'Cannot coerce user-supplied form to a formula object.'})
	if(!(risk.est.method%in%c('QbV','emp.risk','tmle'))) stop('Invalid user-supplied risk.est.method. Please select from QbV, emp.risk, or tmle.')

	if(is.null(bd)){
		fam='gaussian'
		inv.link = identity
	} else {
		fam = 'binomial'
		inv.link = plogis}
	g.vec = rep(g,length(A)) ## **** note : if g becomes a function, this needs to change. Also would need to change for CV-TMLE.
	g.vec[A==0] = 1-g.vec[A==0]


	D1 = (2*A-1)/g.vec * (Y-Qbar0(A,W)) + Qbar0(1,W) - Qbar0(0,W)
	if(!is.null(bd)) D1 = (D1-bd[1])/(bd[2]-bd[1])
	
	if(max(abs(D1))<=1 | fam=='gaussian') { init.fit = glm(reformulate(labels(terms(form)),response='D1'),data=data.frame(W,D1),family=get(fam,mode='function',envir=parent.frame()))
	} else { init.fit = boundedlogistic(reformulate(labels(terms(form)),response='D1'),data=data.frame(W,D1))}

	if(risk.est.method=='QbV'){
		beta.out = coef(init.fit)
	} else {
		if(num.init%%1!=0 | num.init<1){
			stop('num.init must be a positive integer.')
		} else {
			num.coef = length(coef(init.fit))
			# create a set of initial betas to feed to the optimization algorithm
			b.grid = rbind(coef(init.fit))
			if(num.init>1){
				b.grid = rbind(b.grid,matrix(rnorm(num.coef*(num.init-1)),ncol=num.coef))}
			# model matrix to use when evaluating Qbar(V) for a given beta
			mm = model.matrix(form,data=W)
			# function to optimize over
			f = function(beta){
				if(is.null(bd)){ QbV.est = inv.link(c(mm%*%cbind(beta)))
				} else { QbV.est = inv.link(c(mm%*%cbind(beta)))*(bd[2]-bd[1]) + bd[1] }
				d.est = QbV.est>0
				Qb0.dopt = Qbar0(d.est,W)
				if(risk.est.method=='emp.risk'){
					return(-mean((A==d.est)/g.vec * (Y-Qb0.dopt) + Qb0.dopt))
				} else {
					return(-d.tmle(A,Y,Qb0.dopt,g.vec,d.est,bound.Qb0=FALSE)$psi)}}
			num.coef = ncol(b.grid)
			nloptr.out = apply(b.grid,1,function(i){
				i = i/sqrt(sum(i^2))
				nloptr(x0=i,
					eval_f=f,
					lb=rep(-1,num.coef),
					ub=rep(1,num.coef),
					opts=list("algorithm"=nloptr.method,
						"ftol_rel"=1.0e-7,
						"maxeval"=5e4))})
			sols = t(sapply(nloptr.out,function(xx){xx$solution}))
			risk.est = sapply(nloptr.out,function(xx){xx$objective})
			beta.out = sols[which.min(risk.est),]
		}
	}
	out = list(beta=beta.out,bd=bd,form=form,inv.link=inv.link,W=W)
	class(out) = 'param.blip'
	return(out)
}


## pred.param.blip
## Takes as input the output of the param.blip function and generates predictions for a newdata object.
## Note: In the future I should create an object class unifying the parametric and nonparametric algorithms.
## INPUT
## pb : an output from the parametric blip function
## newdata : a new W to use to estimate QbarV
## OUTPUT
## 

predict.param.blip = function(pb,newdata=NULL){
	beta = pb$beta
	bd = pb$bd
	form = pb$form
	inv.link = pb$inv.link

	if(is.null(newdata)) newdata = pb$W
	mm = model.matrix(form,data=newdata)
	if(is.null(bd)){ QbV.est = c(inv.link(mm%*%cbind(beta)))
	} else { QbV.est = inv.link(c(mm%*%cbind(beta)))*(bd[2]-bd[1]) + bd[1] }
	return(QbV.est)}



################################### Two time point ######################################################################
#########################################################################################################################

## Qbar0.mtp
## INPUT
## A : data frame of treatments containing columns A1, A2
## W : data frame of covariates containing columns W1, W2, W3, W4 (though W4 not used)
## OUTPUT
## E[Y|A,W]

Qbar0.mtp = function(A,W){
	0.4 + (-0.4*A$A1- A$A2*W$W3 - 1.5*A$A1*sign(W$W1) + 0.08*A$A1*A$A2 + A$A2*(W$W3)^2 -1.5*A$A1*W$W1 -0.5*A$A1*(W$W2)^2-0.1*A$A2)*0.069}

# c(-.4,0) +
# c(-1.5,1.5) +
# c(-1.5,1.5) +
# c(0,0.08) +
# c(0,1.5^2) +
# c(-1.5,1.5) +
# c(-0.5,0) +
# c(-0.1,0)

Qbar0.mtp.noW34 = function(A,W){
	0.4 + (-0.4*A$A1 - 1.5*A$A1*sign(W$W1) + 0.08*A$A1*A$A2 + A$A2*(A$A1*1.5^2 + (1-A$A1)*0.25^2)/3 -1.5*A$A1*W$W1 - 0.5*A$A1*(W$W2)^2-0.1*A$A2)*0.069}

# Qbar0.mtp.noW34 = function(A,W){0.5 + (0.25*A$A1 + 1.5*A$A1*cos(2*pi*W$W1) + 0.3*A$A1*A$A2 + 0.75*A$A2/3 -1.5*A$A1*W$W1 -A$A1*(W$W2)^2-0.25*A$A2-0.22*A$A1*A$A2)*0.5/5.75}

# Qbar0.mtp = function(A,W){
# 	0.5 + (0.25*A$A1 - A$A2*W$W3 + 1.5*A$A1*cos(2*pi*W$W1) + 0.3*A$A1*A$A2 + 0.75*A$A2*(W$W3)^2 -1.5*A$A1*W$W1 -A$A1*(W$W2)^2-0.25*A$A2-0.22*A$A1*A$A2)*0.5/5.75}

## Qbar0.mtp
## INPUT
## Conditioning on different things.

Qbar0.mtp.ipcw1 = function(A1,W){rep(1/2,length(A1))}
Qbar0.mtp.ipcw2 = function(A2,W){rep(1/2,length(A2))}

Qbar0.mtp.true1 = function(A1,W){
	Qbar0.mtp.noW34(data.frame(A1=A1,A2=W$A2),subset(W,select=c('W1','W2')))}
Qbar0.mtp.true2 = function(A2,W){
	Qbar0.mtp(data.frame(A1=W$A1,A2=A2),subset(W,select=c('W1','W2','W3','W4')))}


# Qbar0.mtp.true1 = function(A1,W){
	# 0.4 + (-0.4*A$A1 - 1.5*A$A1*sign(W$W1) + 0.3*A$A1/2 + (A$A1*1.5^2 + (1-A$A1)*0.25^2)/6 -1.5*A$A1*W$W1 - 0.5*A$A1*(W$W2)^2-0.1/2-0.22*A$A1/2)*0.069}

	

# Qbar0.mtp.true1 = function(A1,W){0.5 + (0.25*A1 + 1.5*A1*cos(2*pi*W$W1) + 0.3*A1*W$A2 + 0.75*W$A2/3 -1.5*A1*W$W1 -A1*(W$W2)^2-0.25*W$A2-0.22*A1*W$A2)*0.5/5.75}
# Qbar0.mtp.true2 = function(A2,W){0.5 + (0.25*W$A1 - A2*W$W3 + 1.5*W$A1*cos(2*pi*W$W1) + 0.3*W$A1*A2 + 0.75*A2*(W$W3)^2 -1.5*W$A1*W$W1 -W$A1*(W$W2)^2-0.25*A2-0.22*W$A1*A2)*0.5/5.75}


## sim.data.mtp
## INPUT
## n : Sample size
## OUTPUT
## A list containing ObsData, the observed data frame,
## and FullData, the data frame containing counterfactuals

sim.data.mtp = function(n){
	W1 = 2*runif(n)-1
	W2 = 2*runif(n)-1
	A1 = rbinom(n,1,1/2)
	u3 = runif(n)
	W3 = (2*u3-1)*(1.25*A1+0.25)
	W3.0 = (2*u3-1)*(0+0.25)
	W3.1 = (2*u3-1)*(1.25+0.25)
	u4 = runif(n)
	W4 = (2*u4-1)*(1.25*A1+0.25)
	W4.0 = (2*u4-1)*(0+0.25)
	W4.1 = (2*u4-1)*(1.25+0.25)
	A2 = rbinom(n,1,1/2)

	A = data.frame(A1=A1,A2=A2)
	W = data.frame(W1=W1,W2=W2,W3=W3,W4=W4)

	u = runif(n)
	Y = as.numeric(u<Qbar0.mtp(A,W))

	Y.00 = as.numeric(u<Qbar0.mtp(data.frame(A1=rep(0,n),A2=rep(0,n)),data.frame(W1=W1,W2=W2,W3=W3.0,W4=W4.0)))
	Y.01 = as.numeric(u<Qbar0.mtp(data.frame(A1=rep(0,n),A2=rep(1,n)),data.frame(W1=W1,W2=W2,W3=W3.0,W4=W4.0)))
	Y.10 = as.numeric(u<Qbar0.mtp(data.frame(A1=rep(1,n),A2=rep(0,n)),data.frame(W1=W1,W2=W2,W3=W3.1,W4=W4.1)))
	Y.11 = as.numeric(u<Qbar0.mtp(data.frame(A1=rep(1,n),A2=rep(1,n)),data.frame(W1=W1,W2=W2,W3=W3.1,W4=W4.1)))

	return(list(ObsData = data.frame(W1=W1,W2=W2,W3=W3,W4=W4,A1=A1,A2=A2,Y=Y),FullData=data.frame(W1=W1,W2=W2,W3.0=W3.0,W3.1=W3.1,W4.0=W4.0,W4.1=W4.1,A1=A1,A2=A2,Y.00=Y.00,Y.01=Y.01,Y.10=Y.10,Y.11=Y.11)))
}

# sim.data.mtp = function(n){
# 	W1 = 2*runif(n)-1
# 	W2 = 2*runif(n)-1
# 	W3 = 2*runif(n)-1
# 	W4 = 2*runif(n)-1
# 	W = data.frame(W1=W1,W2=W2,W3=W3,W4=W4)

# 	A1 = rbinom(n,1,1/2)
# 	A2 = rbinom(n,1,1/2)
# 	A = data.frame(A1=A1,A2=A2)

# 	u = runif(n)
# 	Y = as.numeric(u<Qbar0.mtp(A,W))

# 	return(list(ObsData = data.frame(W1=W1,W2=W2,W3=W3,W4=W4,A1=A1,A2=A2,Y=Y)))
# }


## make.mtp.jointSL
## **** see notes in code. This doesn't extend well to: non-RCT case, multiple time point case, censoring case, 
##	"quadruple" (or doubly doubly) robust case
## Multiple time point case. (**** Currently only implemented for two time points, but can extend to multiple time points)
## INPUT
## L : data frame of covariates
## A : data frame of treatments, where the first column is the first time point, second is second time point, etc.
## Y : vector of outcomes
## QbV.SL.library : SuperLearner library for estimating Qbar(V). Can be NULL
## cl.SL.library : SuperLearner library for classification. Can be NULL
## Qbar0.list : a list of functions that generate the mean outcome in a covariate-treatment strata for each time point
## g : treatment probability in an RCT
## V.list : a list of vectors of column names in W upon which treatment decisions will be made. First one for first time point, etc.
## L.list : a list of vectors of column names in W and earlier versions of A upon which treatment decisions will be made.
##	Must be an increasing list
## alpha.risk.est : Options are: (i) nonnegative least-squares ('nnls'); (ii) empirical risk ('emp.risk'); (iii) CV-TMLE ('cv.tmle'); (iv) TMLE ('tmle')
## bd : bounds on QbarV (called (a,b) in the write-up); to use MSE to fit the SuperLearner, let bd=NULL
## num.folds : number of folds for CV **** (SL or CV-TMLE?)
## convex.loss : convex surrogate loss function to use for classification
## OUTPUT
## A mtp.jointSL object, which is given by a list containing:
##	fit : a list of jointSL fits
##	V.list : same as input

make.mtp.jointSL = function(L,A,Y,QbV.SL.library,cl.SL.library,Qbar0.list,g,V.list,L.list,alpha.risk.est,bd,num.folds=10,convex.loss='log'){
	n = nrow(L)

	# change this if add more time points
	VA.list = V.list
	VA.list[[2]] = c(V.list[[2]],colnames(A)[1])
	# change this if switch to observational -- should be different for each time point
	g.vec = rep(g,n)

	# second time point
	# **** for observational setting need to weight observations according to g!
	# **** for censoring need to restrict to uncensored observations!
	tp2 = make.jointSL(W=data.frame(subset(L,select=L.list[[2]]),A1=A$A1),A=A$A2,Y=Y,QbV.SL.library=QbV.SL.library,cl.SL.library=cl.SL.library,Qbar0=Qbar0.list[[2]],g=g,V=VA.list[[2]],alpha.risk.est=alpha.risk.est,bd=bd,convex.loss=convex.loss,num.folds=num.folds,wgts=1/g.vec)
	d.a0.fixed = predict(tp2,newdata=data.frame(subset(L,select=L.list[[2]]),A1=A$A1))$d.pred

	# first time point
	# **** This won't work when g is different from 0.5, nor for "quadruple" robust case.
	# Qbar0 part probably needs to change, but otherwise good?
	tp1 = make.jointSL(W=data.frame(subset(L,select=L.list[[1]]),A2=A$A2),A=A$A1,Y=Y,QbV.SL.library=QbV.SL.library,cl.SL.library=cl.SL.library,Qbar0=Qbar0.list[[1]],g=g,V=VA.list[[1]],alpha.risk.est=alpha.risk.est,bd=bd,convex.loss=convex.loss,num.folds=num.folds,wgts=(d.a0.fixed==A$A2)/g.vec)
	# tp1 = make.jointSL(W=subset(L,select=L.list[[1]]),A=A$A1,Y=Y,QbV.SL.library=QbV.SL.library,cl.SL.library=cl.SL.library,Qbar0=Qbar0.list[[1]],g=g,V=VA.list[[1]],alpha.risk.est=alpha.risk.est,bd=bd,convex.loss=convex.loss,num.folds=num.folds,wgts=(d.a0.fixed==A$A2)/g.vec)

	out = list(fits=list(tp1,tp2),VA.list=VA.list)
	names(out$fits) = c('Time1','Time2')
	class(out) = 'mtp.jointSL'
	return(out)
}


## predict.mtp.jointSL
## Make predictions based on an mtp.jointSL object.
## INPUT
## obj : mtp.jointSL object, chich contains a list of jointSL objects
## newdata : a W data frame on which to make predictions
## OUTPUT
## A data frame containing treatment decisions in the rows with treatment
## at different time points in the columns.

predict.mtp.jointSL = function(obj,newdata){
	out = data.frame(tmp=rep(NA,nrow(newdata)))[,-1]
	length.obj = length(obj$fits)
	d.df = data.frame()
	for(i in 1:length.obj){
		d.curr = predict(obj$fits[[i]],newdata=subset(cbind(out,newdata),select=obj$VA.list[[i]]))$d.pred
		out = data.frame(out,d.curr)
		colnames(out)[i] = paste('A',i,sep='')
	}
	return(out)
}


## EYd.mtp.tmle
## Given an algorithm for estimating Qbar(V), estimate EYd using the TMLE,
## CV-TMLE, IPTW, and G-comp estimators when classification is included in SL libraries.
## INPUTS
## L : data frame of covariates
## A : vector of treatments
## Y : vector of outcomes
## QbV.SL.library : SuperLearner library for estimating Qbar(V). Can be NULL
## cl.SL.library : SuperLearner library for classification. Can be NULL
## Qbar0.list : ****
## V.list : ****
## L.list : ****
## mtp.jointSL.list : If specified, overrides the two SuperLearner library arguments -- mostly used for simulation purposes
##	A list containing named elements titled:
##		all : mtp.jointSL object fit on entire data set
##		cv : a list of length num.folds, specifying the cross-validated mtp.jointSL fits
##			elements need to correspond to folds as specified in function
## g : treatment probability in an RCT
## alpha.risk.est : Options are: (i) nonnegative least-squares ('nnls'); (ii) empirical risk ('emp.risk'); (iii) CV-TMLE ('cv.tmle'); (iv) TMLE ('tmle')
## bd : bounds on QbarV (called (a,b) in the write-up); to use MSE to fit the SuperLearner, let bd=NULL
## num.folds : number of folds to use when running CV-TMLE (both for the alpha risk estimate and the CV-TMLE EYd estimate)
## non.cv.ests : inference using TMLE and DR-IPCW?
## cv.ests : inference using CV-TMLE and CV-DR-IPCW?
## Q.trunc : truncate Q at Q.trunc and 1-Q.trunc for calculation of fluctuations
## OUTPUT
## ****

EYd.mtp.tmle = function(L,A,Y,QbV.SL.library,cl.SL.library,Qbar0.list,V.list,L.list,Q.SL.library=cl.SL.library,mtp.jointSL.list=NULL,g=0.5,alpha.risk.est='emp.risk',bd=c(-1,1),num.folds=10,non.cv.ests=TRUE,cv.ests=TRUE,Q.trunc=1e-2,use.Qbar0.list=FALSE){
	num.samples = nrow(A)
	# change for non-RCT setting
	g.vec = rep(g,num.samples)

	if(is.null(mtp.jointSL.list)){
		SL.fun = function(LL,AA,YY){make.mtp.jointSL(LL,AA,YY,QbV.SL.library,cl.SL.library,Qbar0.list,g,V.list,L.list,alpha.risk.est,bd,num.folds=num.folds)}
		QbV.fit = SL.fun(L,A,Y)
		d.est = predict(QbV.fit,newdata=L)
	} else {
		QbV.fit = mtp.jointSL.list[[1]]
		d.est = predict(mtp.jointSL.list[[1]],newdata=L)
	}

	if(non.cv.ests){
		# QbV.fit = c(0,1)
		# d.est = data.frame(A1=rep(0,num.samples),A2=rep(1,num.samples))

		if(!is.null(Q.SL.library)){
			Q2.fit = SuperLearner(Y=Y,X=data.frame(subset(L,select=L.list[[2]]),A),SL.library=Q.SL.library,family=binomial())
			Q2.init = c(predict(Q2.fit,newdata=data.frame(subset(L,select=L.list[[2]]),d.est))$pred)
		} else Q2.init = Qbar0.list[[2]](d.est[,2],data.frame(L,A1=d.est[,1]))
		# } else Q2.init = rep(1/2,num.samples)

		tx.match2 = as.numeric(rowSums(abs(as.matrix(d.est)-as.matrix(A)))==0)/(g.vec)^2 # change for non-RCT setting
		eps2 = coef(glm(Y~-1+offset(qlogis(pmin(pmax(Q2.init,Q.trunc),1-Q.trunc)))+tx.match2,family=binomial))[1]
		print(eps2)
		Q2.star = plogis(qlogis(pmin(pmax(Q2.init,Q.trunc),1-Q.trunc)) + eps2*tx.match2)

		if(!is.null(Q.SL.library)){
			Q1.fit = SuperLearner(Y=Q2.init,X=data.frame(subset(L,select=L.list[[1]]),A1=A[,1]),SL.library=Q.SL.library,family=binomial())
			Q1.init = c(predict(Q1.fit,newdata=data.frame(subset(L,select=L.list[[1]]),A1=d.est[,1]))$pred)
			# Q1.dripcw.fit = SuperLearner(Y=Q2.init,X=data.frame(subset(L,select=L.list[[1]]),A1=A[,1]),SL.library=Q.SL.library,family=binomial())
			# Q1.dripcw = c(predict(Q1.dripcw.fit,newdata=data.frame(subset(L,select=L.list[[1]]),A1=d.est[,1]))$pred)
			Q1.dripcw = Q1.init
		} else Q1.init <- Q1.dripcw <- Qbar0.list[[1]](d.est[,1],data.frame(L,A2=d.est[,2]))
		# } else Q1.init <- Q1.dripcw <- rep(1/2,num.samples)

		tx.match1 = as.numeric(abs(d.est[,1]-A[,1])==0)/g.vec
		eps1 = coef(glm(Q2.star~-1+offset(qlogis(pmin(pmax(Q1.init,Q.trunc),1-Q.trunc)))+tx.match1,family=binomial))[1]
		print(eps1)
		Q1.star = plogis(qlogis(pmin(pmax(Q1.init,Q.trunc),1-Q.trunc)) + eps1*tx.match1)

		ic.dripcw = tx.match2*(Y-Q2.init) + tx.match1*(Q2.init-Q1.dripcw) + Q1.dripcw
		dripcw.est = mean(ic.dripcw)
		ci.width.dripcw = qnorm(0.975)*sd(ic.dripcw)/sqrt(num.samples)
		
		tmle.est = mean(Q1.star)
		ci.width.tmle = qnorm(0.975)*sd(tx.match2*(Y-Q2.star) + tx.match1*(Q2.star-Q1.star) + Q1.star-tmle.est)/sqrt(num.samples)

		EYd.est = c(dripcw.est,tmle.est)
		names(EYd.est) = c('dripcw','tmle')
		ci.widths = c(ci.width.dripcw,ci.width.tmle)
		names(ci.widths) = c('dr.ci.width','tmle.ci.width')
	} else {
		EYd.est = NULL
		ci.widths = NULL
	}

	if(cv.ests){
		# Get SL fits for CV-TMLE
		if(num.folds>num.samples){
			warning('More folds than data points. Leave-one-out cross-validation will be used.')
			num.folds = num.samples}
		fold.remainder = num.samples%%num.folds
		if(fold.remainder==0){	
			folds = rep(1:num.folds,each=num.samples/num.folds)
		} else {
			folds = c(rep(1:fold.remainder,each=ceiling(num.samples/num.folds)),
				rep((fold.remainder+1):num.folds,each=floor(num.samples/num.folds)))}
		d.CV.fit = lapply(1:num.folds,function(i){
			if(is.null(mtp.jointSL.list)){
				return(SL.fun(L[folds!=i,],A[folds!=i,],Y[folds!=i]))
			} else {
				return(mtp.jointSL.list[[2]][[i]])}})
		d.CV.1 = do.call(rbind,lapply(1:num.folds,function(i){predict(d.CV.fit[[i]],newdata=L[folds==i,])}))
		d.CV.0 = lapply(1:num.folds,function(i){predict(d.CV.fit[[i]],newdata=L[folds!=i,])})

		if(!is.null(Q.SL.library)){
			Q2.fit.cv = lapply(1:num.folds,function(i){
				return(SuperLearner(Y=Y[folds!=i],X=data.frame(subset(L[folds!=i,],select=L.list[[2]]),A[folds!=i,]),SL.library=Q.SL.library,family=binomial()))})
			Q2.init.1 = c(unlist(lapply(1:num.folds,function(i){
				return(predict(Q2.fit.cv[[i]],newdata=data.frame(subset(L[folds==i,],select=L.list[[2]]),d.CV.1[folds==i,]))$pred)})))
			Q2.init.0 = lapply(1:num.folds,function(i){
				predict(Q2.fit.cv[[i]],newdata=data.frame(subset(L[folds!=i,],select=L.list[[2]]),d.CV.0[[i]]))$pred
				})
		} else {
			Q2.init.1 = c(unlist(lapply(1:num.folds,function(i){
				Qbar0.list[[2]](d.CV.1[folds==i,2],data.frame(L[folds==i,],A1=d.CV.1[folds==i,1]))})))
			Q2.init.0 = lapply(1:num.folds,function(i){
				Qbar0.list[[2]](d.CV.0[[i]][,2],data.frame(L[folds!=i,],A1=d.CV.0[[i]][,1]))})
		}
		# else {
		# 	Q2.init.1 = rep(1/2,num.samples)
		# 	Q2.init.0 = lapply(1:num.folds,function(i){rep(1/2,nrow(d.CV.0[[i]]))})
		# }

		tx.match2.1 = as.numeric(rowSums(abs(as.matrix(d.CV.1)-as.matrix(A)))==0)/(g.vec)^2 # change for non-RCT setting
		eps2.cv = coef(glm(Y~-1 + offset(qlogis(pmin(pmax(Q2.init.1,Q.trunc),1-Q.trunc))) + tx.match2.1,family=binomial,weights=as.numeric(rowSums(abs(as.matrix(d.CV.1)-as.matrix(A)))==0)))[1]
		# eps2.cv = coef(glm(Y~-1 + offset(qlogis(pmin(pmax(Q2.init.1,Q.trunc),1-Q.trunc))) + tx.match2.1,family=binomial))[1]
		print(eps2.cv)

		Q2.star.0 = lapply(1:num.folds,function(i){
			tx.match2.curr = as.numeric(rowSums(abs(as.matrix(d.CV.0[[i]])-as.matrix(A[folds!=i,])))==0)/g^2 # change for non-RCT
			return(plogis(qlogis(pmin(pmax(Q2.init.0[[i]],Q.trunc),1-Q.trunc)) + eps2.cv*tx.match2.curr))})
			# tx.match2.curr = (rowSums(abs(as.matrix(d.CV.0[[i]])-as.matrix(A[folds!=i,])))==0)/(g.vec[folds!=i])^2
			# return(plogis(qlogis(pmin(pmax(Q2.init.0[[i]],Q.trunc),1-Q.trunc)) + eps2.cv*tx.match2.curr))})

		Q2.star.1 = plogis(qlogis(pmin(pmax(Q2.init.1,Q.trunc),1-Q.trunc)) + eps2.cv*tx.match2.1)

		if(!is.null(Q.SL.library)){
			Q1.cvtmle.init = c(unlist(lapply(1:num.folds,function(i){
				tmp = SuperLearner(Y=Q2.star.0[[i]],X=data.frame(subset(L[folds!=i,],select=L.list[[1]]),A1=A[folds!=i,1]),SL.library=Q.SL.library,family=binomial())
				predict(tmp,newdata=data.frame(subset(L[folds==i,],select=L.list[[1]]),A1=d.CV.1[folds==i,1]))$pred})))
			# Q1.cvdripcw.1 = c(unlist(lapply(1:num.folds,function(i){
			# 	tmp = SuperLearner(Y=Q2.init.0[[i]],X=data.frame(subset(L[folds!=i,],select=L.list[[1]]),A1=A[folds!=i,1]),SL.library=Q.SL.library,family=binomial())
			# 	predict(tmp,newdata=data.frame(subset(L[folds==i,],select=L.list[[1]]),A1=d.CV.1[folds==i,1]))$pred})))
			Q1.cvdripcw.1 = Q1.cvtmle.init
		} else {
			Q1.cvtmle.init <- Q1.cvdripcw.1 <- c(unlist(lapply(1:num.folds,function(i){
				Qbar0.list[[1]](d.CV.1[folds==i,1],data.frame(L[folds==i,],A2=d.CV.1[folds==i,2]))
				})))
		}
		# else Q1.cvtmle.init <- Q1.cvdripcw.1 <- rep(1/2,num.samples)

		tx.match1.1 = (abs(as.matrix(d.CV.1[,1])-as.matrix(A)[,1])==0)/g.vec # change for non-RCT setting
		eps1.cv = coef(glm(Q2.star.1 ~ -1 + offset(qlogis(pmin(pmax(Q1.cvtmle.init,Q.trunc),1-Q.trunc))) + tx.match1.1,family=binomial))[1]
		# eps1.cv = coef(glm(Q2.star.1 ~ -1 + offset(qlogis(pmin(pmax(Q1.cvtmle.init,Q.trunc),1-Q.trunc))) + tx.match1.1,family=binomial))[1]
		print(eps1.cv)

		Q1.cvtmle.star = plogis(qlogis(pmin(pmax(Q1.cvtmle.init,Q.trunc),1-Q.trunc)) + eps1.cv*tx.match1.1)
		# Q1.cvtmle.star = plogis(qlogis(pmin(pmax(Q1.cvtmle.init,Q.trunc),1-Q.trunc)) + eps1.cv*tx.match1.1)

		cv.ic.dripcw = tx.match2.1*(Y-Q2.init.1) + tx.match1.1*(Q2.init.1-Q1.cvdripcw.1) + Q1.cvdripcw.1
		cv.dripcw.est = mean(cv.ic.dripcw)
		cv.ci.width.dripcw = qnorm(0.975)*sd(cv.ic.dripcw)/sqrt(num.samples)
		
		cv.tmle.est = mean(Q1.cvtmle.star)
		print(paste('First part:',mean(tx.match2.1*(Y-Q2.star.1)),sep=''))
		print(paste('Second part:',mean(tx.match1.1*(Q2.star.1-Q1.cvtmle.star)),sep=''))
		cv.ci.width.tmle = qnorm(0.975)*sd(tx.match2.1*(Y-Q2.star.1) + tx.match1.1*(Q2.star.1-Q1.cvtmle.star) + Q1.cvtmle.star-cv.tmle.est)/sqrt(num.samples)

		EYd.est = c(EYd.est,cv.dripcw.est,cv.tmle.est)
		names(EYd.est)[3:4] = c('cvdripcw','cvtmle')

		ci.widths = c(ci.widths,cv.ci.width.dripcw,cv.ci.width.tmle)
 		names(ci.widths)[3:4] = c('cvdr.ci.width','cvtmle.ci.width')
	} else d.CV.fit=NULL

	return(list(d = d.est, EYd.est=EYd.est,ci.widths=ci.widths,mtp.jointSL.list=list(all=QbV.fit,cv=d.CV.fit)))
}


## eval.performance.mtp
## Evaluates the performance of each SL fit according to the logistic, squared error, and EYa loss functions using Monte Carlo simulations.
## **** add the other true failures as we go
## INPUTS
## tx.alg : SuperLearner fit that takes A, W1, and W2 as predictors, jointSL object, param.blip object, 
##	or a vector of static treatments
## n.test : Number of values to simulate from the true observed data distribution
## Qbar0 : a function to generate the mean outcome in a covariate-treatment strata
## g : treatment probability in an RCT
## V : variables upon which treatment decisions will be based
## unbd : Were the variables transformed to respect a bound? NULL if not, vector
##	containing bound limits otherwise
## OUTPUT
## A vector containing the estimated risk from each of the three outputs.

eval.performance.mtp = function(tx.alg,n.test=1e6,Qbar0=Qbar0.mtp,g=0.5,unbd=NULL){
	all.data = sim.data.mtp(n.test)
	fd = all.data$FullData
	test.set = all.data$ObsData
	W1 = test.set$W1
	W2 = test.set$W2
	W3 = test.set$W3
	W4 = test.set$W4
	W=data.frame(W1=W1,W2=W2,W3=W3,W4=W4)
	A = subset(test.set,select=c('A1','A2'))
	Y = test.set$Y
	if(is.numeric(tx.alg)){
		if(sum(!(tx.alg%in%c(0,1)))==0){
			mse.true = NA
			mse.ipcw = NA
			d.est = as.data.frame(rbind(tx.alg))[rep(1,n.test),]
			colnames(d.est) = paste('A',1:ncol(d.est),sep='')
		} else {
			stop('Numeric inputs for tx.alg must be 0 or 1.')}
	} else if(is(tx.alg,'mtp.jointSL')) {
		d.est.start = predict(tx.alg,W)
		d1.1 = which(d.est.start$A1==1)
		W.counterfactual = data.frame(W1=W1,W2=W2,W3=fd$W3.0,W4=fd$W4.0)
		W.counterfactual$W3[d1.1] = fd$W3.1[d1.1]
		W.counterfactual$W4[d1.1] = fd$W4.1[d1.1]
		d.est = predict(tx.alg,W.counterfactual)
		## Mean-squared error
		mse.true = NA
		mse.ipcw = NA
	} else {
		stop('Bad tx.alg specification.')}

	Yd = fd$Y.00
	Yd[which((d.est[,1]==0)&(d.est[,2]==1))] = fd$Y.01[which((d.est[,1]==0)&(d.est[,2]==1))]
	Yd[which((d.est[,1]==1)&(d.est[,2]==0))] = fd$Y.10[which((d.est[,1]==1)&(d.est[,2]==0))]
	Yd[which((d.est[,1]==1)&(d.est[,2]==1))] = fd$Y.11[which((d.est[,1]==1)&(d.est[,2]==1))]

	EYd = mean(Yd)

	out.vec = c(mse.true=mse.true,mse.ipcw=mse.ipcw,EYd=EYd,mean.d=mean(unlist(d.est)))
	names(out.vec) = c('mse.true','mse.ipcw','EYd','mean.d')
	return(out.vec)
}

# eval.performance.mtp = function(tx.alg,n.test=1e6,Qbar0=Qbar0.mtp,g=0.5,unbd=NULL){
# 	test.set = sim.data.mtp(n.test)$ObsData
# 	W1 = test.set$W1
# 	W2 = test.set$W2
# 	W3 = test.set$W3
# 	W4 = test.set$W4
# 	W=data.frame(W1=W1,W2=W2,W3=W3,W4=W4)
# 	A = subset(test.set,select=c('A1','A2'))
# 	Y = test.set$Y
# 	if(is.numeric(tx.alg)){
# 		if(sum(!(tx.alg%in%c(0,1)))==0){
# 			mse.true = NA
# 			mse.ipcw = NA
# 			d.est = as.data.frame(rbind(tx.alg))[rep(1,n.test),]
# 			colnames(d.est) = paste('A',1:ncol(d.est),sep='')
# 		} else {
# 			stop('Numeric inputs for tx.alg must be 0 or 1.')}
# 	} else if(is(tx.alg,'mtp.jointSL')) {
# 		d.est = predict(tx.alg,W)
# 		## Mean-squared error
# 		mse.true = NA
# 		mse.ipcw = NA
# 	} else {
# 		stop('Bad tx.alg specification.')}

# 	## EYd
# 	Q.d = Qbar0(d.est,W)
# 	EYd = mean(Q.d)

# 	# EYd = mean(4*(rowSums(abs(as.matrix(d.est)-as.matrix(A)))==0)*Y)

# 	out.vec = c(mse.true=mse.true,mse.ipcw=mse.ipcw,EYd=EYd,mean.d=mean(unlist(d.est)))
# 	names(out.vec) = c('mse.true','mse.ipcw','EYd','mean.d')
# 	return(out.vec)
# }


# ## EYd.est.eval.mtp
# ## Wrapper function for EYd.mtp.tmle and eval.performance.mtp.
# ## INPUTS
# ## See EYd.mtp.tmle
# ## OUTPUT
# ## 

# EYd.est.eval.mtp = function(W,A,Y,SL.library,Qbar0=function(AA,WW){rep(1/2,length(AA))},Qbar0.true=Qbar0.fun,g=0.5,V=c('W1','W2'),alpha.risk.est='nnls',bd=c(-1,1),num.folds=10,W.fun=rnorm){
# 	tmle.out = EYd.mtp.tmle(L,A,Y,QbV.SL.library,cl.SL.library,Qbar0.list=list(Qbar0.mtp.ipcw1,Qbar0.mtp.ipcw2),V.list=V.list,L.list=L.list,Q.SL.library=NULL,mtp.jointSL.list=NULL,g=0.5,alpha.risk.est='emp.risk',bd=NULL,num.folds=10,cv.ests=TRUE)




# 	EYd.out = EYd.tmle(W,A,Y,SL.library,Qbar0=Qbar0,g=g,V=V,alpha.risk.est=alpha.risk.est,bd=bd,num.folds=num.folds)
# 	perf = tryCatch(colMeans(do.call(rbind,lapply(1:10,function(i){eval.performance(EYd.out$QbV.fit,n.test=1e5,Qbar0=Qbar0.true,g=g,unbd=bd,V=V,W.fun=W.fun)}))),
# 		error=function(e){NULL})
# 	out.list = list(EYd.est=EYd.out$EYd.est,ci.width=EYd.out$ci.width,cv.ci.width=EYd.out$cv.ci.width,perf=perf,alpha=EYd.out$QbV.fit$coef)
# 	return(out.list)
# }


## EYd.mc.mtp
## Returns an estimate of the true EYd for a given simulation, using list(c('W1','W2'),c('W1','W2','W3','W4')) as V.list.

EYd0.mtp = function(n.eval=1e5){

	dat = sim.data.mtp(n.eval)
	test.set = dat$ObsData
	W = subset(test.set,select=paste('W',1:4,sep=''))
	A = subset(test.set,select=c('A1','A2'))

	A.11 <- A.00 <- A.01 <- A.10 <- A
	A.11$A1 <- A.11$A2 <- A.10$A1 <- A.01$A2 <- 1
	A.00$A1 <- A.00$A2 <- A.10$A2 <- A.01$A1 <- 0

	Q.00 = Qbar0.mtp(A.00,W)
	Q.01 = Qbar0.mtp(A.01,W)
	Q.10 = Qbar0.mtp(A.10,W)
	Q.11 = Qbar0.mtp(A.11,W)

	Q.0 = Qbar0.mtp.noW34(data.frame(A1=0,A2=as.numeric(Q.01>Q.00)),subset(W,select=c('W1','W2')))
	Q.1 = Qbar0.mtp.noW34(data.frame(A1=1,A2=as.numeric(Q.11>Q.10)),subset(W,select=c('W1','W2')))
	
	first.treat = as.numeric(Q.1>Q.0)
	second.treat = rep(NA,length(first.treat))
	second.treat[which(first.treat==0)] = as.numeric(Q.01>Q.00)[which(first.treat==0)]
	second.treat[which(first.treat==1)] = as.numeric(Q.11>Q.10)[which(first.treat==1)]

	# Q.d = Qbar0.mtp(data.frame(A1=first.treat,A2=second.treat),W)
	# return(mean(Q.d))

	Y.d = rep(NA,n.eval)
	Y.d[which((first.treat==0) & (second.treat==0))] = dat$FullData$Y.00[which((first.treat==0) & (second.treat==0))]
	Y.d[which((first.treat==0) & (second.treat==1))] = dat$FullData$Y.01[which((first.treat==0) & (second.treat==1))]
	Y.d[which((first.treat==1) & (second.treat==0))] = dat$FullData$Y.10[which((first.treat==1) & (second.treat==0))]
	Y.d[which((first.treat==1) & (second.treat==1))] = dat$FullData$Y.11[which((first.treat==1) & (second.treat==1))]

	return(mean(Y.d))
}


## EYd.mc.mtp
## Returns an estimate of the true EYd for a given simulation, using list(c('W1'),c('W1','W3')) as V.list.

# EYd0.mtp2 = function(n.eval=1e5){
# 	Qbar0.mtp = function(A,W){0.5 + (0.25*A$A1 - A$A2*W$W3 + 1.5*A$A1*cos(2*pi*W$W1) + 0.3*A$A1*A$A2 + 0.75*A$A2*(W$W3)^2 -1.5*A$A1*W$W1 -A$A1*(W$W2)^2-0.25*A$A2-0.22*A$A1*A$A2)*0.5/5.75}
# 	Qbar0.mtp.noW24 = function(A,W){0.5 + (0.25*A$A1 - A$A2*W$W3 + 1.5*A$A1*cos(2*pi*W$W1) + 0.3*A$A1*A$A2 + 0.75*A$A2*(W$W3)^2 -1.5*A$A1*W$W1 -A$A1/3-0.25*A$A2-0.22*A$A1*A$A2)*0.5/5.75}
# 	Qbar0.mtp.noW234 = function(A,W){0.5 + (0.25*A$A1 + 1.5*A$A1*cos(2*pi*W$W1) + 0.3*A$A1*A$A2 + 0.75*A$A2/3 -1.5*A$A1*W$W1 -A$A1/3-0.25*A$A2-0.22*A$A1*A$A2)*0.5/5.75}

# 	test.set = sim.data.mtp(n.eval)$ObsData
# 	W = subset(test.set,select=paste('W',1:4,sep=''))
# 	A = subset(test.set,select=c('A1','A2'))

# 	A.11 <- A.00 <- A.01 <- A.10 <- A
# 	A.11$A1 <- A.11$A2 <- A.10$A1 <- A.01$A2 <- 1
# 	A.00$A1 <- A.00$A2 <- A.10$A2 <- A.01$A1 <- 0

# 	Q.00 = Qbar0.mtp.noW24(A.00,subset(W,select=c('W1','W3')))
# 	Q.01 = Qbar0.mtp.noW24(A.01,subset(W,select=c('W1','W3')))
# 	Q.10 = Qbar0.mtp.noW24(A.10,subset(W,select=c('W1','W3')))
# 	Q.11 = Qbar0.mtp.noW24(A.11,subset(W,select=c('W1','W3')))

# 	Q.0 = Qbar0.mtp.noW234(data.frame(A1=0,A2=as.numeric(Q.01>Q.00)),subset(W,select=c('W1')))
# 	Q.1 = Qbar0.mtp.noW234(data.frame(A1=1,A2=as.numeric(Q.11>Q.10)),subset(W,select=c('W1')))
	
# 	first.treat = as.numeric(Q.1>Q.0)
# 	second.treat = rep(NA,length(first.treat))
# 	second.treat[which(first.treat==0)] = as.numeric(Q.01>Q.00)[which(first.treat==0)]
# 	second.treat[which(first.treat==1)] = as.numeric(Q.11>Q.10)[which(first.treat==1)]

# 	Q.d = Qbar0.mtp(data.frame(A1=first.treat,A2=second.treat),W)

# 	return(mean(Q.d))
# }




# predict.mtp.jointSL = function(obj,newdata){
# 	length.obj = length(obj)
# 	all.possible.A = as.data.frame(t(sapply(0:(2^(length.obj-1)-1),function(j){as.numeric(as.character(intToBits(j)))})[1:(length.obj-1),]))
# 	colnames(all.possible.A) = paste('A',1:ncol(all.possible.A),sep='')
# 	out = as.data.frame(t(apply(newdata,function(curr.W){
# 		curr.W.df = as.data.frame(matrix(x,nrow=1))
# 		colnames(curr.W.df) = colnames(newdata)
# 		A.d=list()
# 		d.opt = rep(NA,length.obj)
# 		# forward step
# 		for(i in 1:length.obj){
# 			if(i>1){
# 				all.A = all.possible.A[1:2^(i-1),1:(i-1)]
# 				newdata.curr = cbind(all.A,curr.W.df[rep(1,2^(i-1)),])
# 				d.predict = predict(obj[[i]],newdata=newdata.curr)$d.pred
# 				A.d[[i]] = cbind(all.A,d.predict)
# 			} else d.opt[1] = predict(obj[[i]],newdata=newdata)$d.pred
			
# 		}
# 		# backward step
# 		for(i in seq(length.obj-1,1,by=-1)){
# 			d.opt[length.obj-i+1] = tail(A.d[[i]][which(apply(A.d[[i]],1,function(x){return(sum(d.opt[1:(length.obj-i)]!=head(x,n=-1))==0)})),],n=1)
# 		}
# 		return(d.opt)})))
# 	colnames(out) = paste('A',1:length.obj,sep='')
# 	return(out)
# }








# ############# alternative mtp simulation

# ## Qbar0.mtp
# ## INPUT
# ## A : data frame of treatments containing columns A1, A2
# ## W : data frame of covariates containing columns W1, W2, W3, W4 (though W4 not used)
# ## OUTPUT
# ## E[Y|A,W]

# ff1 = function(x){(x+1)^1.5-1.25*(x+1)}
# ff2 = function(x){(x+0.5)^2-1.25*(x+0.5)}

# Qbar0.mtp = function(A,W){
# 	0.4 + (- 2*A$A1*W$W1*(W$W2+1)/2 + 1.5*A$A2*(W$W3)^2*ff1(W$W1) + A$A1*ff2(W$W2)*plogis(3*W$W1))*0.14}

# # c(-2,2) +
# # 1.5*c(ff1(-11/36),ff1(1)) + 
# # c(ff2(0.125),ff2(-1)) 

# Qbar0.mtp.noW34 = function(A,W){
# 	0.4 + (- 2*A$A1*W$W1*(W$W2+1)/2 + 1.5*A$A2*(A$A1*1^2 + (1-A$A1)*0.25^2)*ff1(W$W1)/3 + A$A1*ff2(W$W2)*plogis(3*W$W1))*0.14}
