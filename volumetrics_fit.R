#! /usr/bin/Rscript --vanilla
rm(list=ls(pos=.GlobalEnv), pos=.GlobalEnv)
Sys.setenv(TZ="Europe/Dublin")
library(foreach)
library(splines)
library(ic.infer)
library(Matrix)
library(plot3D)
library(oce)
library(colorRamps)
dyn.load("./Cpack/libdlam.so")
dyn.load("./Cpack/libmaths.so")
options(CBoundsCheck = TRUE)
source("utils/impro.R")
source("utils/amide.R")
source("utils/tictoc.r")
source("funs/tubemod-funs.R")
source("funs/LM_optim.R")
source("funs/C_dlam.R")
source("funs/C_maths.R")
source("funs/exsarcroi.R")
source("funs/foos.R")
source("funs/plot_phases.R")
source("funs/regul_foos.R")
source("funs/S_rfuns.R")
source("funs/tmi2015_funs.R")

volumetrics <- function(xr,xt,eroi){
# Input:
# - xr is a 5-column data frame (uptake,weigth,x1,x2,x3) in the raw coordinates
# - xt is a 5-column data frame (uptake,weigth,x1,x2,x3) in the PCA coordinates
# - eroi is the output of extract.roi() from package mia 
#   (available at https://github.com/ericwol/mia/releases/tag/1.0.4) 
# An example of use of extract.roi():
# 	eroi = extract.roi(xr,pcut=.25,alpha=5,nb=26,nres=25)
#
	# main settings:
	if(nrow(xr)>3000){
		nb = Jh = 26 
		nres = Jphi = 25 
	} else {
		nb = Jh = 6
		nres = Jphi = dimz
	}
	gam0 = .75
	modcon = list(astep=TRUE,axs=2,bxs=2,dorep=0,
				gammet = 1, gamval = gam0, lm.method = 2,
				pe.min = 1, # minimum percent error to be achieved
				maxit = 10, nbloops.b = 1, nbloops.tau = 1)
	thres = eroi$thres
	#
	xt = xt[,c(1,3:5,2)] # (w,x1,x2,x3,z)
	ab = hab(xt)
	zb = eroi$z.pasbins
	for(k in 1:nb){
		ii = which((xt[,4]<=zb[k+1])&(xt[,4]>zb[k]))
		if(length(ii)){
			xt[ii,2] = xt[ii,2]-as.numeric(eroi$mu1[k])
			xt[ii,3] = xt[ii,3]-as.numeric(eroi$mu2[k])
		}
	}
	#
	Omega = spline.omega(ab)
	#
	# -----------------------------------------------------------------------------------
	# info from PAS domain...
	xr.init = xr
	xa = xr[eroi$icut,]
	
	nx = length(unique(xa[,4]))
	ny = length(unique(xa[,3]))
	nz = length(unique(xa[,5]))
	
	i0 = which(xa[,1]<=0)
	if(length(i0)){
		xa[i0,1] = abs(.001*rnorm(length(i0)))
	}
	
	eroi$s = 1+0*eroi$s
	z = xr[,1]
	hv = xr[,5]
	w = xr[,2]
	xx = xr[,3:5]
	xt = xt[,2:4]

	# ------------------------------------------------------------------------ VOI features
	xx.o = xx
	z.o = z
	xxh = xx

	gx = sort(unique(xx[,2]))
	gy = sort(unique(xx[,1]))
	gz = sort(unique(xx[,3]))
	gxh = sort(unique(xxh[,2]))
	gyh = sort(unique(xxh[,1]))
	gzh = sort(unique(xxh[,3]))
	
	nghbr.xx = list.neighborhood(xx,gx,gy,gz)
	xth = xt
	xinds = c(1:nrow(xx))
	nghbr.xxh = nghbr.xx

	# ------------------------------------------------------------------------ spline stuff
	phia = 0
	phib = 2*pi
	alpha = 3
	
	# -----------------------------------------------------------------------------------
	# Convert initial profile to Gamma profile
	# Reconstruct estimated profiles from fitg() 
	v = eroi$efit[,1] 
	zuy=eroi$zuy
	zz=zuy[,1]   # h-range values
	uu=zuy[,2]   # voxel phases
	yy=zuy[,3]   # uptakes
	nz=Jh
	nn=length(yy) 
	nv=round(max(1,nn/nz)) 
	zb=sort(zz)[c(1,c(1:(nz-1))*nv,nn)] 		# discretized z-grid
	zb[1]=zb[1]-max(.1e-9,.001*(zb[2]-zb[1]))   # so as to include boundaries
	ghat = zhat = uhat = uout = vtst = vout = numeric(nn)
	ghat2 = zhat2 = numeric(nn)
	zza = xt[,4]
	for(k in 1:nz) { 
		kinds = ((zz<=zb[k+1])&(zz>zb[k]))
		yu=yy[kinds]
		ghat[kinds] = yu/eroi$a[k]
		vk=v[kinds]
		ghat2[kinds] = exp(-vk^2/2)
	}
	# from Gaussian to Gamma profile
	g2g = gauss.to.gam(v,alpha,doplot=FALSE,gquant=.992)
	fac = g2g$fac
	
	# -----------------------------------------------------------------------------------
	all.inds = get.inds(Jphi,Jh)
	a.inds = all.inds$a
	tau.inds = all.inds$tau
	b.inds = all.inds$b
	c.inds = all.inds$c
	s.inds = all.inds$s
	xi.inds = all.inds$xi
	mu1.inds = all.inds$mu1
	mu2.inds = all.inds$mu2	
	# structural parameters
	io = init.gpars(eroi,g2g)
	a.init = io$a0
	b.init = io$b0
	tau.init = io$t0
	th0 = numeric(length(unlist(all.inds)))
	th0[c.inds]   = eroi$c
	th0[s.inds]   = eroi$s
	th0[xi.inds]  = eroi$xi
	th0[mu1.inds] = c(matrix(eroi$mu1,nc=1))
	th0[mu2.inds] = c(matrix(eroi$mu2,nc=1))
	th0[a.inds]   = rep(eroi$a,each=Jphi)*fac
	th0[b.inds]   = b.init
	th0[tau.inds] = tau.init
	nb = length(eroi$z.midpoints)
	n = length(z)
	nv = round(max(1,n/nb)) 
	zb = sort(xt[,3])[c(1,c(1:(nb-1))*nv,n)] 
	zb[1] = zb[1]-max(.1e-9,.001*(zb[2]-zb[1])) 	# pas bins
	zb.new = zb 
	zb = eroi$z.pasbins
	zz = xt[,3]	
	rs = us = phis = numeric(n)
	muk = matrix(0,nr=n,nc=2)
	for(k in 1:nb){
		kinds = which((zz<=zb[k+1])&(zz>zb[k]))
		p.out = polar(xt[kinds,1:2],c(0,0),eroi$msg.out$sg[k,],thcorr=2)
		us[kinds] = p.out[,2]
		rs[kinds] = p.out[,3]
		phis[kinds] = p.out[,1]
		muk[kinds,] = matrix(c(eroi$mu1[k],eroi$mu2[k]),nc=2,nr=length(kinds),byr=T)
	}
	rphs = cbind(rs,phis,xt[,3])
	ab = range(rphs[,3])
	x1s=matrix(0,nr=n,nc=Jh)
	x2s=matrix(0,nr=n,nc=(Jh*Jphi))
	x3s=matrix(0,nr=n,nc=Jphi)		
	x = unlist(xt)
	for(i in 1:n) {
		phi = rphs[i,2]
		h = x[i,3]
		x1 = ecbs(h,ab[1],ab[2],Jh)
		x1s[i,]=x1
		x2 = etpb12(c(phi,h),phia,phib,Jphi,ab[1],ab[2],Jh) 
		x2s[i,]=x2
		x3 = epcbs(phi,phia,phib,Jphi)
		x3s[i,]=x3
	}
	orph = list(rphs=rphs,x1s=x1s,x2s=x2s,x3s=x3s)	
	beta0 = eroi$b*g2g$sc.adj
	b0 = rep(beta0,each=Jphi)	
	a0vec = b0vec = tau0vec = numeric(n)
	for(k in 1:nb){
		kinds = which((zz<=zb[k+1])&(zz>zb[k]))
		tau0vec[kinds] = eroi$tau[k]
		a0vec[kinds] = c(eroi$a)[k]
		b0vec[kinds] = beta0[k]
	}
	ps.r = zh.r = zh.ar = numeric(n)
	for(i in 1:n){
		bterm = sum(x2s[i,]*th0[b.inds])
		ps.r[i] = sum(x1s[i,]*th0[tau.inds]) + rs[i]/bterm
		zh.r[i] = ps.r[i]^(alpha-1)*exp(-ps.r[i])*a0vec[i]
		zh.ar[i] = ps.r[i]^(alpha-1)*exp(-ps.r[i])*sum(x2s[i,]*th0[a.inds])
	}

	# -----------------------------------------------------------------------------------
	# Initialize new parametrization using priors from extractroi()
	gyt = sort(unique(xt[,1]))
	gxt = sort(unique(xt[,2]))
	gzt = sort(unique(xt[,3]))
	ntx = length(gxt)
	nty = length(gyt)
	ntz = length(gzt)
	ab = range(rphs[,3])
	
	# Initial params
	ztrue = z
	
	ii = sort(c(c.inds,s.inds,xi.inds,mu1.inds,mu2.inds)) # optim wrt a+b+tau
	exinds = c(1:length(unlist(all.inds)))[ii]
	pis = c(1:length(unlist(all.inds)))[-ii]
	fixstruct = prod(is.element(c(c.inds,s.inds,xi.inds,mu1.inds,mu2.inds),exinds))
	if((length(pis)==length(a.inds))&&(sum(abs(pis-a.inds))==0)){modcon$maxit=1}
	
	theta.true = th0
	# Estimation bounds on theta
	thUL = set.bounds(th0,all.inds)#,method=0,dor=dorep)
	thL = thUL$lower
	thU = thUL$upper
	# reset bounds for non-optimised params
	thL[exinds] = th0[exinds]-1e-12
	thU[exinds] = th0[exinds]+1e-12
	thL[-exinds] = pmin(thL[-exinds],c(th0[-exinds]-1.5))
	thU[-exinds] = pmax(thU[-exinds],c(th0[-exinds]+1.5))
	# mu's
	mu.inds = sort(c(mu1.inds,mu2.inds))
	# a
	if(!mmatch(exinds,a.inds)){
		thL[a.inds] = 1e-2 #pmax(1e-3,theta.true[a.inds]-10) 
		thU[a.inds] = theta.true[a.inds]+10
	}
	# b
	if(!mmatch(exinds,b.inds)){
		thU[b.inds] = theta.true[b.inds]+4*abs(theta.true[b.inds])
		thL[b.inds] = pmax(1e-3,theta.true[b.inds]/10)
	}
	# tau
	if(!mmatch(exinds,tau.inds)){
		thL[tau.inds] = rep(.2,length(tau.inds))
		thU[tau.inds] = rep(5,length(tau.inds))
	}	
	th0 = pmax(thL,pmin(th0,thU))
	theta.true = th0 
	
	# LX.R = eval.lam.R(orph,th0[tau.inds],th0[a.inds],th0[b.inds],alpha)
	LX.C = eval.lam.C(orph,th0[tau.inds],th0[a.inds],th0[b.inds],alpha)
	LX = LX.C
	r1 = orph$rphs[xinds,1]
	
	# ----------------------------------------------------------------------------
	NPA = length(unlist(all.inds))-length(exinds) # nb of params being optimised
	fixstruct = as.logical(prod(is.element(c(c.inds,s.inds,xi.inds,mu1.inds,mu2.inds),exinds)))

	# ----------------------------------------------------------------------------
	# run main optim code (this part is replicated from optim_routine_exsarcroi_parallel.R)
	tic = Sys.time()
	gam = modcon$gamval
	nocritchange = 0
	RSSs.nog  = RSSs = RSSs.loops = added = ks = ns = ginds = NULL
	theta.all = theta.vec = trhat = diffths = NULL
	gamma.grid = gamma.crit = rss.critg  = NULL
	lap.vec = normk = NULL
	grid.gammas = c(seq(5e-3,0.995,.05),1.25)
	theta.init = th0
	theta1 = theta.init
	
	rph.i = orph
	rphs = orph$rphs[xinds,]
	x1s = orph$x1s[xinds,]
	x2s = orph$x2s[xinds,]
	x3s = orph$x3s[xinds,]

	Omega.tau = spline.omega1(hab(xth),Jh)

	n = nrow(xt)
	M0s = M1s = M01s = RSS.all.steps = NULL
	gamma.grid = gamma.crit = gams = NULL
	thetak.a = thetak.b = thetar.a = thetar.b = thetar.t = thetak.t = NULL
	thetas.all = NULL
	thetas.all.b = NULL
	thetas.all.mu = NULL
	thetas.all.tau = NULL
	
	pe.curr = 100
	kk=0
	nochange = FALSE
	
	while((kk<modcon$maxit)&&(abs(pe.curr)>modcon$pe.min)&&(!nochange)){
		kk = kk+1

		theta00 = theta1 # keep this in case we decide to reject the update
		theta0 = theta1

		LXH.0 = eval.lam.C(rph.i,theta0[tau.inds],theta0[a.inds],theta0[b.inds],alpha)
		GLXH.0 = eval.dlam.C(nghbr.xxh,LXH.0)
		M0s = c(M0s, (sum((z-LXH.0[xinds])^2 + gam*GLXH.0^2)))
		if((kk==1)){
			if(modcon$gammet>=0){
				eg = eval.grads.rd(xx,orph,nghbr.xxh,xinds,exinds,theta0,Jh,Jphi,alpha)
			} else {
				eg = list(DG=NULL,DGL=NULL)
			}
		}

		gamg = sort(c( gam, gam*(seq(0.9,0.1,l=4)), gam*(seq(1.1,5,l=5))))
		# print(paste("Selecting gamma... current value is",round(gam,5)))
		sgam.o1 = select.gam(eg$DG,eg$DGL,z,gammas=gamg,method=modcon$gammet,gamval=gam)
		gam = sgam.o1$gamma
		# print(paste("done selecting gamma:",round(gam,5)))
		gamma.grid = rbind(gamma.grid, sgam.o1$gammas)
		gamma.crit = rbind(gamma.crit, sgam.o1$crit)
		if(gam){
			ccrit = (sum((z-LXH.0[xinds])^2 + gam*GLXH.0^2))
		} else {
			ccrit = (sum((z-LXH.0[xinds])^2))
		}
		RSS.all.steps = c(RSS.all.steps, ccrit)
		
		### b-step
		if(!mmatch(exinds,b.inds)){
			# print(paste("... b-step"))
			for(i in 1:modcon$nbloops.b){
				theta0 = theta1			
				exbinds = sort(as.numeric(unlist(all.inds)[-b.inds]))
				LXH.a = eval.lam.C(orph,theta0[tau.inds],theta0[a.inds],theta0[b.inds],alpha)
				GLXH.a = eval.dlam.C(nghbr.xxh,LXH.a)[xinds]
				if(gam==0){
					ok=0
					try({
						Xb = bgrad(orph,theta0[tau.inds],theta0[a.inds],theta0[b.inds],alpha)[xinds,];
						X2.0 = Xb; # dummy
						ok=1;
					})
				}
				if((gam!=0)||(sum(is.na(c(Xb))))){
					eg = eval.grads.rd(xx,orph,nghbr.xxh,xinds,exbinds,theta0,Jh,Jphi,alpha)
					Xb = eg$DG
					X2.0 = eg$DGL
				}
				# pseudo-values:
				zstar = z-LXH.a[xinds]+c(Xb%*%theta0[b.inds])
				# matrix to ML-pseudo-inverse:
				# Ma  = crossprod(Xb)+gam*crossprod(X2.0) 
				Ma = fxuprod(Xb)+gam*fxuprod(X2.0)
				ob = regul.beta.wrap(Xb,zstar,Omega,bnd=thL[b.inds],
					gam=gam,XD=X2.0,G0=GLXH.a,b0=theta0[b.inds])
				thetak.b = rbind(thetak.b, ob$beta.ols)
				thetar.b = rbind(thetar.b, ob$beta.orls)
				theta1[b.inds] = ob$beta.orls
				# theta1 = pmin(thU,theta1)
				thetas.all.b = rbind(thetas.all.b,theta1)
			}
		}
		LXH.1 = eval.lam.C(orph,theta1[tau.inds],theta1[a.inds],theta1[b.inds],alpha)
		GLXH.1 = eval.dlam.C(nghbr.xxh,LXH.1)[xinds]
		if(gam){
			ccrit = (sum((z-LXH.1[xinds])^2 + gam*GLXH.1^2))
		} else {
			ccrit = (sum((z-LXH.1[xinds])^2))
		}
		RSS.all.steps = c(RSS.all.steps, ccrit)
		theta0.b = theta0
		theta1.b = theta1
		
		### tau-step (finite-diffs, usual Marquardt step)
		if(0){
			if(!mmatch(exinds,tau.inds)){
				# print(paste("... tau-step"))
				extauinds = sort(as.numeric(unlist(all.inds)[-tau.inds]))
				for(i in 1:modcon$nbloops.tau){
					theta0 = theta1
					LXH.0 = eval.lam.C(orph,theta0[tau.inds],theta0[a.inds],theta0[b.inds],alpha)
					GLXH.0 = eval.dlam.C(nghbr.xxh,LXH.0)
					eg = eval.grads.rd(xx,orph,nghbr.xxh,xinds,extauinds,theta0,Jh,Jphi,alpha)
					X1.0 = eg$DG
					X2.0 = eg$DGL
					delta = 5e-5
					# lmtau0 = ifelse(kk==1,delta*max(diag(crossprod(X1.0)+gam*crossprod(X2.0))),ml.tau$tau)
					lmtau0 = delta*max(diag(crossprod(X1.0)+gam*crossprod(X2.0)))
					ml.tau = marquardt.beta.rd(z,theta0,LXH.0,GLXH.0,X1.0,X2.0,gam,lmtau0,modcon$lm.method)
					theta1 = ml.tau$theta
					theta1[tau.inds] = pmin(thU[tau.inds],pmax(thL[tau.inds],theta1[tau.inds]))
					thetas.all.tau = rbind(thetas.all.tau,theta1)
				}
			}
		} else {
			if(!mmatch(exinds,tau.inds)){
				# print(paste("... tau-step"))
				extauinds = sort(as.numeric(unlist(all.inds)[-tau.inds]))
				for(i in 1:modcon$nbloops.tau){
					theta0 = theta1
					LXH.0 = eval.lam.C(orph,theta0[tau.inds],theta0[a.inds],theta0[b.inds],alpha)
					GLXH.0 = eval.dlam.C(nghbr.xxh,LXH.a)[xinds]
					eg = eval.grads.rd(xx,orph,nghbr.xxh,xinds,extauinds,theta0,Jh,Jphi,alpha)
					X1.0 = eg$DG
					X2.0 = eg$DGL
					# pseudo-values:
					zstar = z-LXH.0[xinds]+c(X1.0%*%theta0[tau.inds])
					# matrix to ML-pseudo-inverse:
					# Ma  = crossprod(X1.0)+gam*crossprod(X2.0) 
					Ma  = fxuprod(X1.0)+gam*fxuprod(X2.0) 
					ob = regul.tau.wrap(X1.0,zstar,Omega.tau,bnd=thL[tau.inds],
						gam=gam,XD=X2.0,G0=GLXH.0,b0=theta0[tau.inds])
					thetak.t = ob$beta.orls
					thetak.t = rbind(thetak.t, ob$beta.ols)
					thetar.t = rbind(thetar.t, ob$beta.orls)				
					theta1[tau.inds] = ob$beta.orls
					theta1[tau.inds] = pmin(thU[tau.inds],pmax(thL[tau.inds],theta1[tau.inds]))
					thetas.all.tau = rbind(thetas.all.tau,theta1)
				}
			}
		}
		LXH.1 = eval.lam.C(orph,theta1[tau.inds],theta1[a.inds],theta1[b.inds],alpha)
		GLXH.1 = eval.dlam.C(nghbr.xxh,LXH.1)[xinds]
		if(gam){
			ccrit = (sum((z-LXH.1[xinds])^2 + gam*GLXH.1^2))
		} else {
			ccrit = (sum((z-LXH.1[xinds])^2))
		}
		RSS.all.steps = c(RSS.all.steps, ccrit)
		
		### a-step
		if(!mmatch(exinds,a.inds)){
			# print(paste("... a-step"))
			theta0 = theta1
			exainds = sort(as.numeric(unlist(all.inds)[-a.inds]))
			LXH.0 = eval.lam.C(orph,theta0[tau.inds],theta0[a.inds],theta0[b.inds],alpha)
			GLXH.0 = eval.dlam.C(nghbr.xxh,LXH.0)[xinds]
			if(gam==0){
				ok=0
				try({
				Xa = agrad(orph,theta0[tau.inds],theta0[a.inds],theta0[b.inds],alpha)$gvals[xinds,];
				X2.0 = Xa; # dummy setting
				ok=1;
				})
			} 
			if((gam!=0)||(sum(is.na((c(Xa)))))){
				eg = eval.grads.rd(xx,orph,nghbr.xxh,xinds,exainds,theta0,Jh,Jphi,alpha)
				Xa = eg$DG
				X2.0 = eg$DGL
			}
			# matrix to ML-pseudo-inverse:
			# Ma  = crossprod(Xa)+gam*crossprod(X2.0)
			Ma  = fxuprod(Xa)+gam*fxuprod(X2.0)
			oas = regul.beta.wrap(Xa,z,Omega,bnd=thL[a.inds],gam=gam,
						XD=X2.0,G0=GLXH.0,b0=theta0[a.inds])
			thetak.a = rbind(thetak.a, oas$beta.ols)
			thetar.a = rbind(thetar.a, oas$beta.orls)
			theta1[a.inds] = oas$beta.orls
		}
		LXH.1 = eval.lam.C(orph,theta1[tau.inds],theta1[a.inds],theta1[b.inds],alpha)
		GLXH.1 = eval.dlam.C(nghbr.xxh,LXH.1)[xinds]
		if(gam){
			ccrit = (sum((z-LXH.1[xinds])^2 + gam*GLXH.1^2))
		} else {
			ccrit = (sum((z-LXH.1[xinds])^2))
		}
		RSS.all.steps = c(RSS.all.steps, ccrit)
		theta1.a = theta1

		if((kk>1)&&(ccrit>(rev(M1s)[1]))){
			nochange=TRUE
			theta1 = theta00
		} else {
			gams = c(gams, gam)
			M1s = c(M1s, ccrit)
			pe.curr = ((M1s[kk]-M0s[kk])/M0s[kk]*100)
			M01s = c(M01s, pe.curr)
			### record current estimate 
			thetas.all = rbind(thetas.all,theta1)
		}
	} # end of kk-for-loop
	toc = Sys.time()
	tictoc = (toc-tic)
	((kk<modcon$maxit)&&(abs(pe.curr)>modcon$pe.min)&&(!nochange))
	flag = 3 # "some other exit cause"
	if(nochange){
		flag = 1 #print("--- Terminated due to increasing RSS")
	} else {
		if(kk==modcon$maxit){
			flag = 2 # print("--- Terminated after reaching maxiter")
		}
	}
	if(abs(pe.curr)<modcon$pe.min){
		flag = 0 #print("--- Terminated after achieving minimal %-change")
	}

	M1 = ccrit
	zhat = eval.lam.C(orph,theta1[tau.inds],theta1[a.inds],theta1[b.inds],alpha)[xinds]
		
	# ----------------------------------------------------------------------------
	# hack for foreach combining...
	gpT = get.phases(theta.true,z,orph,xinds,alpha,all.inds)
	gp0 = get.phases(th0,z,orph,xinds,alpha,all.inds)
	gp1 = get.phases(theta1,z,orph,xinds,alpha,all.inds)
	# --- hets and phases ---
	nb = length(eroi$z.midpoints)
	n = nrow(xt)
	nv = round(max(1,n/nb))
	zb = eroi$z.pasbins
	zz = xt[,3]
	res0 = res1 = resT = 0*zz
	hets0 = hets1 = hetsT = mp = core.mp = numeric(Jh)
	for(k in 1:Jh){
	    kinds = which((zz<=zb[k+1])&(zz>zb[k]))
	    resT[kinds] = z[kinds]-gpT$zhat[kinds]
	    res0[kinds] = z[kinds]-gp0$zhat[kinds]
	    res1[kinds] = z[kinds]-gp1$zhat[kinds]
	    hetsT[k] = max(0, 1-(mean(resT[kinds]^2)/var(z[kinds])))
	    hets0[k] = max(0, 1-(mean(res0[kinds]^2)/var(z[kinds])))
	    hets1[k] = max(0, 1-(mean(res1[kinds]^2)/var(z[kinds])))
	    mp[k] = median(gp1$phases[kinds])
	}
	HT = 100-100*hetsT
	H0 = 100-100*hets0
	H1 = 100-100*hets1

	# ----------------------------------------------------------------------------
	regmod.out = list(thT=theta.true,thU=thU,thL=thL,th0=th0,th1=theta1,
			orph=orph,thetas.all=thetas.all,
			thetak.a=thetak.a,thetak.b=thetak.b,
			Jh=Jh,
			z=z,zz=zz,zhat=zhat,xinds=xinds,
			gp0=gp0,gp1=gp1,
			hets0=hets0,hets1=hets1,
			H0=H0,H1=H1,HT=HT,
			tictoc=tictoc)
	return(regmod.out)
}

# ------------------------------------------------------------------------ funs...
minmax.scale <- function(yh){
	return(	(yh-min(yh)) / (max(yh)-min(yh)))
}
get.gradients <- function(u,yh,resign=TRUE,rescale=TRUE){
	# u = u/median(ifelse(u>.001,u,.001))
	if(rescale){yh = minmax.scale(yh)}
	grad = (approx(u,yh,xout=u*1.1,rule=2)$y-approx(u,yh,xout=u*.9,rule=2)$y)/(ifelse(u>.01,u,.01)*.2)
	if(resign){grad = -grad}
	return(grad)
}
add.alpha <- function(COLORS, ALPHA){
	if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
	RGB <- col2rgb(COLORS, alpha=TRUE)
	RGB[4,] <- round(RGB[4,]*ALPHA)
	NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
	return(NEW.COLORS)
}
interp.slice <- function(slice){
	gx = c(1:dim(slice)[2])
	gy = c(1:dim(slice)[1]) 
	gz = 1
	tix = set.grid.from.midpoint(gx,.5)
	tiy = set.grid.from.midpoint(gy,.5) 
	tiz = gz
	slicea = array(slice,dim=c(dim(slice),1))
	return(trilin.approx(gx,gy,gz,slicea,tix,tiy,tiz)[,,1])
}
na.rm <- function(x){
	return(x[which(!is.na(x))])
}
quantize <- function(x,bins){
	N = length(x)
	B = length(bins)
	X = matrix(rep(x,each=B), byrow=TRUE, ncol=B, nrow=N)
	G = matrix(rep(bins,each=N), byrow=FALSE, ncol=B, nrow=N)
	d = unlist(apply(abs(X-G),1,which.min))
	o = numeric(N)*NA
	o[which(!is.na(x))] = bins[d]
	return(list(values=c(o),indices=d))
}

