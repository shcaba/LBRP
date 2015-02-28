#Assn5.R
## ___________________________________________________________________________
## Management Procedure Evalualtion
## Programer: Steve Martell
## Date: 09/07/2009
## Requires: 
##			Estimators: Schaefer.exe, LRGS.exe
## Notes:
## 1) The operating model is an age structured model (ASAM from class)
## 2) This is modified from the MPE assignment for the OSU crowd for MPE 
##
##
## ___________________________________________________________________________
## ___________________________________________________________________________

## DATA SECTION
setwd("~/Documents/CURRENT PROJECTS/OSU MPE Complexity/Martell Model/")
hakedata <-read.table("NamibianHake.data",header=T)
yr		<-	hakedata[,1]
yt		<-	hakedata[,2] #CPUE (units are in tons per effort)
ct		<-	hakedata[,3]/1000 #Catch (1000 tons)
n		<-	length(yr)
yrs	<- c(yr,max(yr)+1)
A		<-	25
age		<-	1:A
winf		<-	5.0
aa		<- 1e-5
bb		<- 3
linf	<- exp(log(winf/aa)/bb)
k		<-	0.14
adot		<-	3.5
tau.adot	<-	3.75
Phi		<-	c(winf=winf,k=k,adot=adot,tau.adot=tau.adot)
wa	<-	winf*(1-exp(-k*age))^3
la 	<- 	linf*(1-exp(-k*age))
ma	<-	plogis(age,adot,1/tau.adot)
dx	<- seq(as.integer(0.5*min(la)),as.integer(max(la)*1.35),by=1)
fa	<-	wa*ma
ALK = sapply(dx+0.5,pnorm,mean=la,sd=0.1*la)-sapply(dx-0.5,pnorm,mean=la,sd=0.1*la)

## PARAMETER SECTION
theta=c(Bo=3.7018792,kap=33.2988473,m=0.2622931,ahat=4.0758066,tau.ahat=2.3093062,tau.y=8.4240941)



## PROCEDURE SECTION
Asam	<-	function(theta,yr,ct)
{
	with(as.list(theta),
	{
		va	<-	plogis(age,ahat,1/(tau.ahat))
		iota	<-	exp(-m)^(age-1)
		iota[A] <- iota[A]/(1-exp(-m))
		phiE	<- 	sum(iota*wa*ma)
		phiB	<-	sum(iota*wa)
		so	<-	kap/phiE
		b	<-	(kap-1)/(Bo/phiB*phiE)
		N	<-	matrix(0,nrow=length(yr)+1,ncol=A)
		N[1,]	<- 	Bo/phiB*iota
		for(i in 1:length(yr))
		{
			#if(ct[i]>sum(N[i,]*wa*va)) ct[i]=0.9*sum(N[i,]*wa*va)
			Ut	<-	ct[i]/sum(N[i,]*wa*va)
			
			#Could add bioeconomic routine here to dynamically modify fishing mortality rates
			#based on economic variables, depreciation and reinvestment.

			Et	<-	sum(N[i,]*wa*ma)
			N[i+1,1]	<-	so*Et/(1+b*Et)
			N[i+1,2:A]	<-	N[i,1:(A-1)]*exp(-m)*(1-Ut*va[1:(A-1)])
			N[i+1,A]	<-	N[i+1,A]+N[i,A]*exp(-m)*(1-Ut*va[A])
		}
		Bt	<-	as.vector((wa*va) %*%t(N))
		SBt	<-	as.vector((wa*ma) %*%t(N))

		#Observations for Froese's metrics Pmat, Popt, Pmega
		xx = iota*wa*va
		optage=which(xx==max(xx))
		lopt= linf*(1-exp(-k*optage))
		lmat = 0.65*lopt  #based on Cope and Punt, 2009.
		iopt = seq(min(which(dx>=0.9*lopt)),min(which(dx>=1.1*lopt)),by=1)
		imat = which(dx>=lmat)
		imega = which(dx>=1.1*lopt)
		Pmat = Popt = Pmega = vector()
		for(i in 1:length(yr))
		{
			px = as.vector(t(ALK)%*%(N[i,]*va))  #this is the probability of sampling a fish of length dx
			fx = as.vector(rmultinom(1,100,px))	#this is a random sample of 500 lengths.
			fx = fx/sum(fx)
			Pmat[i] = sum(fx[imat])
			Popt[i] = sum(fx[iopt])
			Pmega[i] = sum(fx[imega])
			
			
		}
		Pobj = rowSums(cbind(Pmat,Popt,Pmega))
#		par(mfcol=c(3,2))
#		plot(dx,fx,type="s")
#		lines(dx,fx,type="h")
#		matplot(yr,cbind(Pmat,Popt,Pmega,Pobj),type="l",col=1:4,lty=1:4,ylim=c(0,2))
#		lines(yr,Bt[1:length(yr)]/Bt[1],lwd=2)
#		
#		legend("top",c("Pmat","Popt","Pmega","Pobj"),lty=1:4,col=1:4,ncol=4,bty="n")		
#		plot(SBt[1:length(yr)],Pmat,type="o"); abline(h=0.9)
#		plot(SBt[1:length(yr)],Pobj,type="o"); abline(h=0.9)
#		plot(SBt[1:length(yr)],Popt,type="o")
#		plot(SBt[1:length(yr)],Pmega,type="o")
#		
#		print(lm(SBt[1:length(yr)]~Pobj))
		#Observation model
		#zt		<-log(yt)-log(Bt[1:n])
		#epsilon	<-zt-mean(zt)

		#Statistical criterion
		#nloglike	<- -1.*sum(dnorm(epsilon,0,1/tau.y,log=T))
		
		#We want to minimize the following.
		#f <- nloglike
		#print(c(length(Bt),length(ct),length(yr)))
		return(list(Bt=Bt, ut=ct/Bt[1:length(yr)], Pobj = Pobj[length(yr)]))
	})
}
write.data	<- function(d)
{
	fn <- "MPEdata.dat"
	write(file=fn,dim(d)[1])
	write.table(d,file=fn,append=T,row.names=F,col.names=F)
}




operating.model	<- function(iseed=12345,np=20,estimator="schaefer",scenario="S1",harvest.rule="H1")
{
	## This function is the operating model that runs for
	## np years using the true parameter values.
	## RM is the  reference model
	#1) Generate random variates for current iseed
	set.seed(iseed)
	pyrs		<- c(yr,max(yr)+1:np)
	n		<- length(pyrs)
	epsilon	<- rnorm(n,0,0.3)
	
	#Objects for storing information
	bt	<- matrix(NA,nrow=np+1,ncol=n+1)
	o.ct	<- vector(length=n)
	
	#2) Loop over projection years (pyrs)
	for(i in 0:np)
	{
		#2) Use reference model to generate new data.
		yr	<- seq(1965,1987+i,by=1)
		RM	<- Asam(theta,yr,ct)
		switch(scenario,
			S1={ yt <- RM$Bt[1:length(yr)]*exp(epsilon[1:length(yr)])},
			S2={yt <- (RM$Bt[1:length(yr)]^0.7)*exp(epsilon[1:length(yr)])}
			)
		
		
		#3) Run estimator to determine biomass for HR
		write.data(cbind(yr,round(yt,3),round(ct,3)))
		
		switch(estimator,
			schaefer={
				#system("Schaefer.exe -ind MPEdata.dat -nox -est",show.output=F);
				system("./Schaefer -ind MPEdata.dat -nox -est",show.output=F);
				fmsy	<-scan("Schaefer.rep",skip=1,nlines=1,quiet=T)
				btmp	<-scan("Schaefer.rep",skip=3,nlines=1,quiet=T)
				bo	<-scan("Schaefer.rep",skip=5,nlines=1,quiet=T)},
			lrgs={
				#system("LRGS.exe -ind MPEdata.dat -nox -est",show.output=F);
				system("./LRGS -ind MPEdata.dat -nox -est",show.output=F);
				fmsy	<-scan("LRGS.rep",skip=1,nlines=1,quiet=T)
				btmp	<-scan("LRGS.rep",skip=3,nlines=1,quiet=T)
				bo	<-scan("LRGS.rep",skip=5,nlines=1,quiet=T)},
			datarule={
				fmsy <- 0.2
				btmp	<-rep(0,length=length(yr)+1)
				btmp[length(yr)+1] <- 6.167*RM$Pobj - 6.573
				bo <- 6.167*1.55 - 6.573		#the 1.55 is the Pobj at unfished conditions
						}
			)
		
		#4) Use harvest control rule to get next years catch
		switch(harvest.rule,
			H1 = {tac	<- fmsy*btmp[length(yr)+1]},
			H2 = {tac	<- 0.9*fmsy*btmp[length(yr)+1]},
			H3 ={	ft	<- 0
				d	<- btmp[length(yr)+1]/bo
				if(d>=0.4)ft=fmsy
				if(d>=0.1 & d<0.4) ft=fmsy*(d-0.1)/0.3
				tac	<- ft*btmp[length(yr)+1]}
			)
		
		#The following line assume the fishery can only harvest 90% of the stock
		#if the tac actually exceeds the true stock size.
		if(tac>=RM$Bt[length(yr)+1] | is.na(tac)) tac=0.9*RM$Bt[length(yr)+1]
		ct	<- c(ct,tac)
		
		#5) Store retrospective statistics
		bt[i+1,1:length(btmp)]	<- btmp
		
	}
	#6) Calculate performance measures
	o.ct	<-	ct
	st	<-	(length(yr)-np+1):(length(yr)-np+6)
	lt	<-	(length(yr)-np+7):(length(yr))
	cbar	<-	c(mean(ct[st]),mean(ct[lt]))
	aav	<-	c(sum(abs(ct[st]-ct[st-1]))/sum(ct[st]),sum(abs(ct[lt]-ct[lt-1]))/sum(ct[lt]))
	dbar	<-	c(mean(RM$Bt[st]/RM$Bt[1]),mean(RM$Bt[lt]/RM$Bt[1]))
	
	##optional graphics to display while running
	yrs<-c(yr,max(yr)+1)
	matplot(yrs,t(bt),type="l",lty=1,col="grey",ylim=c(0,6))	
	lines(yrs,RM$Bt,lwd=2)
	plot(yrs,ct,type="h",lwd=2)
	
	return(list(e.bt=bt,bt=bt[np+1,],ct=o.ct,Bt=RM$Bt,cbar=cbar,aav=aav*100,dbar=dbar))
}

## MAIN SECTION
par(mfcol=c(2,2),las=1)
nrep	<-	50
np	<-	25
if(!exists("S1")){
	CMP1.1<-lapply(1:nrep,operating.model,np=np,estimator="schaefer",scenario="S1",harvest.rule="H1")
	CMP2.1<-lapply(1:nrep,operating.model,np=np,estimator="schaefer",scenario="S1",harvest.rule="H2")
	CMP3.1<-lapply(1:nrep,operating.model,np=np,estimator="schaefer",scenario="S1",harvest.rule="H3")
	CMP4.1<-lapply(1:nrep,operating.model,np=np,estimator="lrgs",scenario="S1",harvest.rule="H1")
	CMP5.1<-lapply(1:nrep,operating.model,np=np,estimator="lrgs",scenario="S1",harvest.rule="H2")
	CMP6.1<-lapply(1:nrep,operating.model,np=np,estimator="lrgs",scenario="S1",harvest.rule="H3")
	CMP7.1<-lapply(1:nrep,operating.model,np=np,estimator="datarule",scenario="S1",harvest.rule="H1")
	CMP8.1<-lapply(1:nrep,operating.model,np=np,estimator="datarule",scenario="S1",harvest.rule="H2")
	CMP9.1<-lapply(1:nrep,operating.model,np=np,estimator="datarule",scenario="S1",harvest.rule="H3")
	S1	<- list(MP1=CMP1.1,MP2=CMP2.1,MP3=CMP3.1,MP4=CMP4.1,MP5=CMP5.1,MP6=CMP6.1,MP7=CMP7.1,MP8=CMP8.1,MP9=CMP9.1)
}
if(!exists("S2")){
	CMP1.2<-lapply(1:nrep,operating.model,np=np,estimator="schaefer",scenario="S2",harvest.rule="H1")
	CMP2.2<-lapply(1:nrep,operating.model,np=np,estimator="schaefer",scenario="S2",harvest.rule="H2")
	CMP3.2<-lapply(1:nrep,operating.model,np=np,estimator="schaefer",scenario="S2",harvest.rule="H3")
	CMP4.2<-lapply(1:nrep,operating.model,np=np,estimator="lrgs",scenario="S2",harvest.rule="H1")
	CMP5.2<-lapply(1:nrep,operating.model,np=np,estimator="lrgs",scenario="S2",harvest.rule="H2")
	CMP6.2<-lapply(1:nrep,operating.model,np=np,estimator="lrgs",scenario="S2",harvest.rule="H3")
	CMP7.2<-lapply(1:nrep,operating.model,np=np,estimator="datarule",scenario="S2",harvest.rule="H1")
	CMP8.2<-lapply(1:nrep,operating.model,np=np,estimator="datarule",scenario="S2",harvest.rule="H2")
	CMP9.2<-lapply(1:nrep,operating.model,np=np,estimator="datarule",scenario="S2",harvest.rule="H3")
	S2	<- list(MP1=CMP1.2,MP2=CMP2.2,MP3=CMP3.2,MP4=CMP4.2,MP5=CMP5.2,MP6=CMP6.2,MP7=CMP8.2,MP8=CMP8.2,MP9=CMP9.2)
}
#SET UP GRAPHICS DEVICE FOR PLOTTING (page-up page-down to scroll through figs)
graphics.off()
#if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
#windows(record=T); par(las=1,ps=12,cex.main=14/12)
par(mfcol=c(2,2))

worm.plot	<- function(CMP)
{
	#This function draws the worm plots given a CMP object.
	get.Bt <- function(CMP) CMP$Bt
	Bt<-sapply(CMP,get.Bt)
	B90<-apply(Bt,1,quantile,probs=c(0.05,0.5,0.95))
	ni <- dim(Bt)[2]
	yrs<-seq(1967,1967+dim(Bt)[1]-1)
	matplot(yrs,t(B90),type="n",ylim=c(0,max(Bt)))
	polygon(c(yrs,rev(yrs)),c(B90[1,],rev(B90[3,])),col="grey",lty=0)
	matlines(yrs,Bt[,c(sample(1:ni,3))],col=1,lty=1)
	lines(yrs,B90[2,],type="l",lwd=2,pch=19)
	
	get.Ct <- function(CMP) CMP$ct
	Ct<-sapply(CMP,get.Ct)
	C90<-apply(Ct,1,quantile,probs=c(0.05,0.5,0.95))
	ni <- dim(Ct)[2]
	yrs<-seq(1967,1967+dim(Bt)[1]-1)
	matplot(yrs,t(C90),type="n",ylim=c(0,max(Ct)))
	polygon(c(yrs,rev(yrs)),c(C90[1,],rev(C90[3,])),col="grey",lty=0)
	matlines(yrs,Ct[,c(sample(1:ni,3))],col=1,lty=1)
	lines(yrs,C90[2,],type="l",lwd=2,pch=19)
}

plot.stats	<-	function(S)
{	#This function plot summary statistics for a given scenario
	par(mfcol=c(3,2),mar=c(4,4,1,1),oma=c(1,1,1,1))
	bxplts<-function(S)
	{
		fn1<-function(MP,t=1){
			get.cbar <-function(MP) MP$cbar[t]
			c.bar<-sapply(MP,get.cbar)
			return(c.bar)
		}
		Cbar1<-sapply(S,fn1,t=1)
		Cbar2<-sapply(S,fn1,t=2)
		boxplot(as.data.frame(cbind(Cbar1,Cbar2)),col=2,ylim=c(0,0.5),ylab="Average catch")
		abline(v=9.5,lty=3); text(5,0,"Short term"); text(14,0,"Long term")
		
		fn2<-function(MP,t=1){
			get.aav <-function(MP) MP$aav[t]
			aa.v<-sapply(MP,get.aav)
			return(aa.v)
		}
		AAV1<-sapply(S,fn2,t=1)
		AAV2<-sapply(S,fn2,t=2)
		boxplot(as.data.frame(cbind(AAV1,AAV2)),col=3,ylim=c(0,100),ylab="AAV (%)")
		abline(v=9.5,lty=3); text(5,0,"Short term"); text(14,0,"Long term")
		
		fn3<-function(MP,t=1){
			get.dbar <-function(MP) MP$dbar[t]
			d.bar<-sapply(MP,get.dbar)
			return(d.bar)
		}
		Dbar1<-sapply(S,fn3,t=1)
		Dbar2<-sapply(S,fn3,t=2)
		boxplot(as.data.frame(cbind(Dbar1,Dbar2)*100),col="cyan",ylim=c(0,100),ylab="Average depletion")
		abline(v=9.5,lty=3); text(5,0,"Short term"); text(14,0,"Long term")
		abline(h=25,lty=2)
	}
	sapply(S,bxplts)
}

performance.stats<-function(MP)
{
	get.cbar	<-function(MP) MP$cbar
	c.bar<-sapply(MP,get.cbar)
	cbar<-apply(c.bar,1,quantile,prob=0.5)
	
	get.aav	<-function(MP) MP$aav
	aa.v<-sapply(MP,get.aav)
	aav<-apply(aa.v,1,quantile,prob=0.5)
	
	get.dbar	<-function(MP) MP$dbar
	d.bar<-sapply(MP,get.dbar)
	dbar<-apply(d.bar,1,quantile,prob=0.5)

	stats<-cbind(cbar=cbar,aav=aav,dbar=dbar)
	row.names(stats)=c("Short term","Long term")
	return(round(stats,3))
}

if(exists("S1")){
	par(mfrow=c(9,2),mar=c(0.2,2,0,0),oma=c(4,3,2,4),xaxt="n")
	sapply(S1,worm.plot); mtext(c("Year","Biomass","Scenario 1","Catch"),side=1:4,las=0,line=c(1,1,0,0),outer=T,font=c(0,0,2,0))
	sapply(S2,worm.plot); mtext(c("Year","Biomass","Scenario 2","Catch"),side=1:4,las=0,line=c(1,1,0,0),outer=T,font=c(0,0,2,0))
	
	#Performance statistics in the form of list of tables
	P1<-lapply(S1,performance.stats)
	P2<-lapply(S2,performance.stats)
	plot.stats(list(S1,S2))
	dev.copy2eps(file="MPEstats.eps")
}