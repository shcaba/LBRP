### needed functions ###
HCR.VBGF<-function(Linf,k,t0,ages)
{
Lengths_exp<-Linf*(1-exp(-k*(ages-t0)))
return(Lengths_exp)
}


HCR.LtWt.fit<-function(p,Lts)
{
  exp.wts<-p[1]*(Lts)^p[2]
  return(exp.wts)
}

Calc_Pobjs<-function(Comps_in,L_bins,Lmat,Lopt0.9,Lopt1.1,year.pick,gender.pick)
{
 Lt.comp.fl.gen.yr<-subset(Comps_in,Year == year.pick)
 if(gender.pick<2){Lt.comp.prep<-Lt.comp.fl.gen.yr[7:(length(L_bins)+6)]}
 if(gender.pick==2){Lt.comp.prep<-Lt.comp.fl.gen.yr[(length(L_bins)+7):ncol(Lt.comp.fl.gen.yr)]}
 Lmat_ind<-L_bins>=Lmat
 Lopt_ind<-L_bins>Lopt0.9 &  L_bins<Lopt1.1
 Lmega_ind<-L_bins>=Lopt1.1

 if(nrow(Lt.comp.prep)>0)
 {
 prop.comp<-as.numeric(colSums(Lt.comp.prep))/as.numeric(sum(Lt.comp.prep))
 Pmat<-sum(prop.comp[Lmat_ind])
 Popt<-sum(prop.comp[Lopt_ind])
 Pmega<-sum(prop.comp[Lmega_ind])
 Pobj<-Pmat+Popt+Pmega
 Pmat.opt.mega<-matrix(c(Pmat,Popt,Pmega,Pobj),nrow=4,ncol=1)
 rownames(Pmat.opt.mega)<-c("Pmat","Popt","Pmega","Pobj")
 colnames(Pmat.opt.mega)<-"Value"
 }
 else {Pmat.opt.mega<-matrix(NA,nrow=4,ncol=1)}
 return(Pmat.opt.mega)
}

Froese.plus<-function(LtComp.dat,L_bins,Linf,k,t0,M,a.wl,b.wl,Lmat,ages,fleet.pick,gender.pick=0,adj.neff="T") # add Catch.series to weight by catch series.
{
  #Calculate Lopt
  Lts<-HCR.VBGF(Linf,k,t0,ages)
  Wts<-HCR.LtWt.fit(c(a.wl,b.wl),Lts)
  stable.biomass<-(1*exp(-M*ages))*Wts
  Lopt<-Lts[stable.biomass==max(stable.biomass)]
  Lopt0.9<-Lopt*0.9
  Lopt1.1<-round(1.1*Lopt,1)
  #Re-weight according to effective sample sizes
  if(adj.neff=="T")
    {
      Tot.samps<-rowSums(LtComp.dat[,7:ncol(LtComp.dat)])
      Samp.perc<-LtComp.dat[,7:ncol(LtComp.dat)]/Tot.samps
#      if(exists("Catch.series")==T)
#        {
#          for(i in 1:length(fleet.pick))
#            {
#
#            }
#          LtComp.dat[,7:ncol(LtComp.dat)]<-Samp.perc*LtComp.dat$Neff*rowSums(Catch.series[,fleet.pick])
#        }
      LtComp.dat[,7:ncol(LtComp.dat)]<-Samp.perc*LtComp.dat$Neff
    }
  Lt.comp.fleet<-subset(LtComp.dat,Fleet == fleet.pick[1])
  if(exists("Catch.series")==T)
    {
     Lt.comp.fleet[,7:ncol(LtComp.dat)]<-Lt.comp.fleet[,7:ncol(LtComp.dat)]*Catch.series[Catch.series[,ncol(Catch.series)] %in% Lt.comp.fleet[,1],fleet.pick[1]]
    }
  if(length(fleet.pick)>1)
  {
    for(f in 2:length(fleet.pick))
      {
       Lt.comp.fleet.temp<-subset(LtComp.dat,Fleet == fleet.pick[f])
       if(exists("Catch.series")==T)
        {
        Lt.comp.fleet.temp[,7:ncol(LtComp.dat)]<-Lt.comp.fleet.temp[,7:ncol(LtComp.dat)]*Catch.series[Catch.series[,ncol(Catch.series)] %in% Lt.comp.fleet.temp[,1],fleet.pick[f]]
        }
       Lt.comp.fleet<-rbind(Lt.comp.fleet,Lt.comp.fleet.temp)
      }
  }
  Lt.comp.fl.gen<-subset(Lt.comp.fleet,Gender == gender.pick[1])
  if(length(gender.pick)>1)
  {
    for(g in 2:length(gender.pick))
      {
       Lt.comp.fl.gen.temp<-subset(Lt.comp.fleet,Gender == gender.pick[g])
       Lt.comp.fl.gen<-rbind(Lt.comp.fl.gen,Lt.comp.fl.gen.temp)
      }
  }
  yrs<-sort(unique(Lt.comp.fl.gen$Year))
  Pobjs.out<-matrix(NA,nrow=4,ncol=length(yrs),dimnames=list(c("Pmat","Popt","Pmega","Pobj"),yrs))
  for(i in 1:length(yrs))
  {
    cab.Pobjs.yr<-Calc_Pobjs(Lt.comp.fl.gen,L_bins,Lmat,Lopt0.9,Lopt1.1,yrs[i],gender.pick)
    Pobjs.out[,i]<-as.numeric(cab.Pobjs.yr[,1])
  }
#  Pobjs<-Calc_Pobjs(LtComp.dat,Lbins,Lmat,Lopt0.9,Lopt1.1,year.pick,Fleet.pick,gender.pick)
  Pcalcs.out<-list()
  Pcalcs.out[[1]]<-Pobjs.out
  Pcalcs.out[[2]]<-c(Lmat,Lopt,Lopt1.1)
  names(Pcalcs.out[[2]])<-c("Lmat","Lopt","Lmega")
  Pcalcs.out[[3]]<-Lmat/Lopt
  names(Pcalcs.out)[[1]]<-"Pout"
  names(Pcalcs.out)[[2]]<-"Lx"
  names(Pcalcs.out)[[3]]<-"Lmat/Lopt"
  return(Pcalcs.out)
}

plot.Pobj<-function(Pobjs.in,y.in,y.int,RPs_h,h.in=T,h.spp.in,RP.comp="TP",plot.leg=T)
{
#Determine and collect decision tree RPs
  RPs.plot.in<-as.data.frame(matrix(NA,nrow=5,ncol=length(Pobjs.in$Pout[4,])))
  rownames(RPs.plot.in)<-c("Px","RPtarget","RPlimit","Pt.bg","Px_type")
  colnames(RPs.plot.in)<-colnames(Pobjs.in$Pout)
  #Define column index based on Lmat/Lopt
  if(Pobjs.in$`Lmat/Lopt`<=0.65+(0.75-0.65)/2){col.i<-c(2:4)}
  if(Pobjs.in$`Lmat/Lopt`>0.65+(0.75-0.65)/2 & Pobjs.in$`Lmat/Lopt`<=0.75+(0.9-0.75)/2){col.i<-c(5:7)}
  if(Pobjs.in$`Lmat/Lopt`>=0.75+(0.9-0.75)/2){col.i<-c(8:10)}
  if(RP.comp=="TP"){rp.i<-2}
  if(RP.comp=="LP"){rp.i<-3}
  #Extract Px and RP values
  for(i in 1:length(Pobjs.in$Pout[4,]))
  {
    if(Pobjs.in[[1]][4,i]<=1) #Poj <=1
    {
      RPs.plot.in[1,i]<-Pobjs.in$Pout[1,i]
      if(h.in==F)(RPs.plot.in[2,i]<-RPs.plot.in[3,i]<-RPs_h$triggers[i,col.i[1]]-1)
      if(h.in==T)
      {
        RPs.plot.in[2,i]<-RPs_h$target[RPs_h$target$h %in% round(h.spp.in,1),col.i[1]]
        RPs.plot.in[3,i]<-RPs_h$limit[RPs_h$limit$h %in% round(h.spp.in,1),col.i[1]]
      }
      if(RPs.plot.in[1,i]>=RPs.plot.in[rp.i,i]){RPs.plot.in[4,i]<-"black"}
      if(RPs.plot.in[1,i]<RPs.plot.in[rp.i,i]){RPs.plot.in[4,i]<-"white"}
      RPs.plot.in[5,i]<-"Pmat"
    }
    if(Pobjs.in$Pout[4,i]>1 & Pobjs.in$Pout[4,i]<2)   #1<Poj<1.99
    {
      RPs.plot.in[1,i]<-Pobjs.in$Pout[1,i]
      if(h.in==F)(RPs.plot.in[2,i]<-RPs.plot.in[3,i]<-RPs_h$triggers[i,col.i[2]]-1)
      if(h.in==T)
      {
        RPs.plot.in[2,i]<-RPs_h$target[RPs_h$target$h %in% round(h.spp.in,1),col.i[2]]
        RPs.plot.in[3,i]<-RPs_h$limit[RPs_h$limit$h %in% round(h.spp.in,1),col.i[2]]
      }
      if(RPs.plot.in[1,i]>=RPs.plot.in[rp.i,i]){RPs.plot.in[4,i]<-"black"}
      if(RPs.plot.in[1,i]<RPs.plot.in[rp.i,i]){RPs.plot.in[4,i]<-"white"}
      RPs.plot.in[5,i]<-"Pmat"
    }
    if(Pobjs.in$Pout[4,i]>=2) #Poj >=1.99
    {
      RPs.plot.in[1,i]<-Pobjs.in$Pout[1,i]
      if(h.in==F)(RPs.plot.in[2,i]<-RPs.plot.in[3,i]<-RPs_h$triggers[i,col.i[3]]-1)
      if(h.in==T)
      {
        RPs.plot.in[2,i]<-RPs_h$target[RPs_h$target$h %in% round(h.spp.in,1),col.i[3]]
        RPs.plot.in[3,i]<-RPs_h$limit[RPs_h$limit$h %in% round(h.spp.in,1),col.i[3]]
      }
      if(RPs.plot.in[1,i]<=RPs.plot.in[rp.i,i]){RPs.plot.in[4,i]<-"black"}
      if(RPs.plot.in[1,i]>RPs.plot.in[rp.i,i]){RPs.plot.in[4,i]<-"white"}
      RPs.plot.in[5,i]<-"Popt"
    }
  }
#Plot
  Years<-as.numeric(colnames(RPs.plot.in))
  plot(Years,RPs.plot.in[1,],xlim=c(min(Years),max(Years)),ylim=y.in,xlab="",ylab=expression(P[x]),type="p",lwd=2,pch=21,col="black",bg=as.character(RPs.plot.in[4,]),cex=1.25,axes=F)
  box()
  axis(1)
  axis(2,at=seq(y.in[1],y.in[2],by=y.int),label=c(0,seq(y.in[1]+y.int,y.in[2],by=y.int)))
  if(y.in[1]!=0){axis.break(2,y.in[1]+y.int/2,style="slash",brw=0.05)}
#  axis(3,at=colnames(RPs.plot.in)),label=Dep.in)
  points(Years,RPs.plot.in[2,],pch=3,col="darkgreen")
  points(Years,RPs.plot.in[3,],pch=4,col="red")
  if(plot.leg == T){smartlegend("left","bottom",c(expression(paste(P[x], "> RP",sep="")),expression(paste(P[x], "< RP",sep="")),"Target RP","Limit RP"),lty=c(NA,NA,NA,NA),lwd=c(NA,NA,NA,NA),pch=c(21,21,3,4),col=c("black","black","green","red"),pt.bg=c("black","white",NA,NA),pt.cex=1.25,bty="n")}
  return(RPs.plot.in)
}

Spp.bp.plot<-function(spp.dep.in,spp.px.in,Dep.targ=0.4,plot.leg=T)
{
  spp.bp.dat<-cbind(spp.dep.in,t(spp.px.in[4,]))
  spp.bp.plot<-matrix(NA,nrow=nrow(spp.bp.dat),ncol=ncol(spp.bp.dat))
  colnames(spp.bp.plot)<-c("Year","Dep","F","LtRP")
  spp.bp.plot[,1]<-spp.bp.dat[,1]
  spp.bp.plot[spp.bp.dat[,2]>=Dep.targ,2]<-1
  spp.bp.plot[spp.bp.dat[,2]<Dep.targ,2]<-(-1)
  spp.bp.plot[spp.bp.dat[,3]<=1,3]<-1
  spp.bp.plot[spp.bp.dat[,3]>1,3]<-(-1)
  spp.bp.plot[spp.bp.dat[,4]=="black",4]<-1
  spp.bp.plot[spp.bp.dat[,4]=="white",4]<-(-1)
  plot(spp.bp.plot[,1],spp.bp.plot[,2]*0.4,ylim=c(-1,1),xlab="Year",ylab="",type="p",pch=21,col="black",bg="gray",axes="F")
  points(spp.bp.plot[,1],spp.bp.plot[,3]*0.5,pch=22,col="black",bg="gray")
  points(spp.bp.plot[,1],spp.bp.plot[,4]*0.6,pch=24,col="black",bg="gray")
  abline(h=0)
  box()
  axis(1)
  if(plot.leg == T){smartlegend("right","bottom",c("Depletion","Fishing mortality","Lt RP"),pch=c(21,22,24),col="black",pt.bg="gray",pt.cex=1.25,bty="n")}
  return(spp.bp.plot)
}

#Get number of agreement and number of conservative designations, first for Depletion, then Fishing mortality
count.agree<-function(Comp.in,perc="T")
{
  count.out<-c(NA,NA,NA,NA,0)
  count.x<-rep(1,nrow(Comp.in))
  count.out[1]<-sum(count.x[Comp.in[,2]==Comp.in[,4]])
  count.out[2]<-sum(count.x[Comp.in[,2]>Comp.in[,4]])
  count.out[3]<-sum(count.x[Comp.in[,3]==Comp.in[,4]])
  count.out[4]<-sum(count.x[Comp.in[,3]>Comp.in[,4]])
  if(perc=="T"){count.out<-count.out/nrow(Comp.in)}
  count.out[5]<-nrow(Comp.in)
  return(count.out)
}

###############################################################
#Calculate some needed VBGF parameters (mostly, t0)
VBGF.dat.spp<-read.table("clipboard", header=T)
fitvbgf.black.g1<-nls(lts~Linf*(1-exp(-k*(ages-t0))),data=list(ages=VBGF.dat.spp$Ages,lts=VBGF.dat.spp$Black_lt_g1),start=list(Linf=max(VBGF.dat.spp$Black_lt_g1),k=0.1,t0=0),control = list(reltol=0.00000000001))
fitvbgf.black.g2<-nls(lts~Linf*(1-exp(-k*(ages-t0))),data=list(ages=VBGF.dat.spp$Ages,lts=VBGF.dat.spp$Black_lt_g2),start=list(Linf=max(VBGF.dat.spp$Black_lt_g2),k=0.1,t0=0),control = list(reltol=0.00000000001))
fitvbgf.blue<-nls(lts~Linf*(1-exp(-k*(ages-t0))),data=list(ages=c(0:(length(na.omit(VBGF.dat.spp$Blue_Lt))-1)),lts=na.omit(VBGF.dat.spp$Blue_Lt)),start=list(Linf=max(VBGF.dat.spp$Blue_Lt,na.rm=T),k=0.1,t0=0),control = list(reltol=0.00000000001))
fitvbgf.boc<-nls(lts~Linf*(1-exp(-k*(ages-t0))),data=list(ages=c(0:(length(na.omit(VBGF.dat.spp$Boc_lt))-1)),lts=na.omit(VBGF.dat.spp$Boc_lt)),start=list(Linf=max(VBGF.dat.spp$Boc_lt,na.rm=T),k=0.1,t0=0),control = list(reltol=0.00000000001))
fitvbgf.widow.g1<-nls(lts~Linf*(1-exp(-k*(ages-t0))),data=list(ages=c(0:(length(na.omit(VBGF.dat.spp$Widow_lt_g1))-1)),lts=na.omit(VBGF.dat.spp$Widow_lt_g1)),start=list(Linf=max(VBGF.dat.spp$Widow_lt_g1,na.rm=T),k=0.1,t0=0),control = list(reltol=0.00000000001))
fitvbgf.widow.g2<-nls(lts~Linf*(1-exp(-k*(ages-t0))),data=list(ages=c(0:(length(na.omit(VBGF.dat.spp$Widow_lt_g2))-1)),lts=na.omit(VBGF.dat.spp$Widow_lt_g2)),start=list(Linf=max(VBGF.dat.spp$Widow_lt_g2,na.rm=T),k=0.1,t0=0),control = list(reltol=0.00000000001))

####################################
### Calculate spcecies Px values ###
####################################
CabN.LtComps<-read.table("clipboard", header=T)
CabN.ages<-c(0:30)
CabN.Linf<-58.97
CabN.k<-0.21
CabN.t0<-(-1.28)
CabN.M<-0.25
CabN.Lmat<-34.6
CabN.a.wl<-0.0000092
CabN.b.wl<-3.187
CabN.Pobjs.All<-Froese.plus(CabN.LtComps,seq(6,92,by=2),CabN.Linf,CabN.k,CabN.t0,CabN.M,CabN.a.wl,CabN.b.wl,CabN.Lmat,CabN.ages,fleet.pick=c(1:6),gender.pick=0)

CabS.LtComps<-read.table("clipboard", header=T)
CabS.ages<-c(0:30)
CabS.Linf<-58.97
CabS.k<-0.21
CabS.t0<-(-1.28)
CabS.M<-0.25
CabS.Lmat<-34.6
CabS.a.wl<-0.0000092
CabS.b.wl<-3.187
CabS.Pobjs.All<-Froese.plus(CabS.LtComps,seq(6,92,by=2),CabS.Linf,CabS.k,CabS.t0,CabS.M,CabS.a.wl,CabS.b.wl,CabS.Lmat,CabS.ages,fleet.pick=c(2:6),gender.pick=0)

CabOR.LtComps<-read.table("clipboard", header=T)
CabOR.ages<-c(0:30)
CabOR.Linf<-68.83
CabOR.k<-0.2
CabOR.t0<-(-2.26)
CabOR.M<-0.25
CabOR.Lmat<-43.67
CabOR.a.wl<-0.00000918
CabOR.b.wl<-3.188
CabOR.Pobjs.All<-Froese.plus(CabOR.LtComps,seq(6,92,by=2),CabOR.Linf,CabOR.k,CabOR.t0,CabOR.M,CabOR.a.wl,CabOR.b.wl,CabOR.Lmat,CabOR.ages,fleet.pick=c(1:6),gender.pick=0)

Black.LtComps<-read.table("clipboard", header=T)
Black.ages<-c(0:60)
Black.Linf<-49.5345
Black.k<-0.17074
Black.t0<-(-2)
Black.M<-0.24
Black.Lmat<-39.53
Black.a.wl<-0.00001677
Black.b.wl<-3
BlackN.Pobjs.All<-Froese.plus(Black.LtComps,seq(20,56,by=2),Black.Linf,Black.k,Black.t0,Black.M,Black.a.wl,Black.b.wl,Black.Lmat,Black.ages,fleet.pick=c(1:3),gender.pick=0)
BlackS.Pobjs.All<-Froese.plus(Black.LtComps,seq(20,56,by=2),Black.Linf,Black.k,Black.t0,Black.M,Black.a.wl,Black.b.wl,Black.Lmat,Black.ages,fleet.pick=c(7:12,15,21),gender.pick=0)

Blue.LtComps<-read.table("clipboard", header=T)
Blue.ages<-c(0:60)
Blue.Linf<-38.4573
Blue.k<-0.111716
Blue.t0<-(-2.6)
Blue.M<-0.1
Blue.Lmat<-26
Blue.a.wl<-0.00003408
Blue.b.wl<-2.874
Blue.Pobjs.All<-Froese.plus(Blue.LtComps,seq(10,52,by=2),Blue.Linf,Blue.k,Blue.t0,Blue.M,Blue.a.wl,Blue.b.wl,Blue.Lmat,Blue.ages,fleet.pick=c(1:4),gender.pick=0)

Bocaccio.LtComps<-read.table("clipboard", header=T)
Bocaccio.ages<-c(0:70)
Bocaccio.Linf<-67.738
Bocaccio.k<-0.21958
Bocaccio.t0<-(-0.7023)
Bocaccio.M<-0.15
Bocaccio.Lmat<-39.9
Bocaccio.a.wl<-0.000007355
Bocaccio.b.wl<-3.11359
Bocaccio.Pobjs.All<-Froese.plus(Bocaccio.LtComps,seq(16,76,by=2),Bocaccio.Linf,Bocaccio.k,Bocaccio.t0,Bocaccio.M,Bocaccio.a.wl,Bocaccio.b.wl,Bocaccio.Lmat,Bocaccio.ages,fleet.pick=c(1,3,4,5,7,8,9,10),gender.pick=0)

Gstriped.LtComps<-read.table("clipboard", header=T)
Gstriped.ages<-c(0:80)
Gstriped.Linf<-33.4
Gstriped.k<-0.105
Gstriped.t0<-(-2.51)
Gstriped.M<-0.08
Gstriped.Lmat<-20.97
Gstriped.a.wl<-0.0000074
Gstriped.b.wl<-3.167
Gstriped.Pobjs.All<-Froese.plus(Gstriped.LtComps,seq(10,40,by=1),Gstriped.Linf,Gstriped.k,Gstriped.t0,Gstriped.M,Gstriped.a.wl,Gstriped.b.wl,Gstriped.Lmat,Gstriped.ages,fleet.pick=c(1,2,4,5),gender.pick=0)

Widow.LtComps<-read.table("clipboard", header=T)
Widow.ages<-c(0:80)
Widow.Linf<-46.29
Widow.k<-0.2
Widow.t0<-(-1)
Widow.M<-0.125
Widow.Lmat<-37
Widow.a.wl<-0.00545
Widow.b.wl<-3.28781
Widow.Pobjs.All<-Froese.plus(Widow.LtComps,seq(10,64,by=2),Widow.Linf,Widow.k,Widow.t0,Widow.M,Widow.a.wl,Widow.b.wl,Widow.Lmat,Widow.ages,fleet.pick=c(1,2,3,4,11,12),gender.pick=0)

#############################################

#Read in steepness values for more specific Px evaluations
h.spp<-c(0.6,0.6,0.58,0.7,0.7,0.7,0.573,0.692,0.406)
names(h.spp)<-c("BlackN","BlackS","Blue","CabOR","CabN","CabS","Bocaccio","Gstriped","Widow")
#Read in Px trigger tables (Table 4 of Cope and Punt 2009)
RP_targs_h<-read.table("clipboard", header=T)
RP_lims_h<-read.table("clipboard", header=T)
RP_targs_trig<-read.table("clipboard", header=T)
RPs_h<-list(RP_targs_h,RP_lims_h,RP_targs_trig)
names(RPs_h)<-c("target","limit","triggers")

#Determine species-specific h trigger points
#TP
par(mfrow=c(2,2))
BlackN.px.TP<-plot.Pobj(BlackN.Pobjs.All,c(0.3,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[1],plot.leg=F)
BlackS.px.TP<-plot.Pobj(BlackS.Pobjs.All,c(0,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[2],plot.leg=F)
par(mfrow=c(2,2))
Blue.px.TP<-plot.Pobj(Blue.Pobjs.All,c(0,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[3],plot.leg=F)
Boc.px.TP<-plot.Pobj(Bocaccio.Pobjs.All,c(0,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[7],plot.leg=F)
par(mfrow=c(2,2))
CabOR.px.TP<-plot.Pobj(CabOR.Pobjs.All,c(0,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[6],plot.leg=F)
CabN.px.TP<-plot.Pobj(CabN.Pobjs.All,c(0.3,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[4],plot.leg=F)
CabS.px.TP<-plot.Pobj(CabS.Pobjs.All,c(0,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[5],plot.leg=F)
par(mfrow=c(2,2))
Gstr.px.TP<-plot.Pobj(Gstriped.Pobjs.All,c(0.6,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[8],plot.leg=F)
Widow.px.TP<-plot.Pobj(Widow.Pobjs.All,c(0,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[9],plot.leg=F)

#LP
par(mfrow=c(2,2))
BlackN.px.LP<-plot.Pobj(BlackN.Pobjs.All,c(0.3,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[1],RP.comp="LP",plot.leg=F)
BlackS.px.LP<-plot.Pobj(BlackS.Pobjs.All,c(0,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[2],RP.comp="LP",plot.leg=F)
par(mfrow=c(2,2))
Blue.px.LP<-plot.Pobj(Blue.Pobjs.All,c(0,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[3],RP.comp="LP",plot.leg=F)
Boc.px.LP<-plot.Pobj(Bocaccio.Pobjs.All,c(0,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[7],RP.comp="LP",plot.leg=F)
par(mfrow=c(2,2))
CabOR.px.LP<-plot.Pobj(CabOR.Pobjs.All,c(0,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[6],RP.comp="LP",plot.leg=F)
CabN.px.LP<-plot.Pobj(CabN.Pobjs.All,c(0.3,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[4],RP.comp="LP",plot.leg=F)
CabS.px.LP<-plot.Pobj(CabS.Pobjs.All,c(0,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[5],RP.comp="LP",plot.leg=F)
par(mfrow=c(2,2))
Gstr.px.LP<-plot.Pobj(Gstriped.Pobjs.All,c(0.6,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[8],RP.comp="LP",plot.leg=F)
Widow.px.LP<-plot.Pobj(Widow.Pobjs.All,c(0,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[9],RP.comp="LP",plot.leg=F)


#Read in depletion and 1-SPR/1-SPR50% values
BlackN.dep<-read.table("clipboard", header=T)
BlackS.dep<-read.table("clipboard", header=T)
Blue.dep<-read.table("clipboard", header=T)
Bocaccio.dep<-read.table("clipboard", header=T)
CabOR.dep<-read.table("clipboard", header=T)
CabN.dep<-read.table("clipboard", header=T)
CabS.dep<-read.table("clipboard", header=T)
Gstriped.dep<-read.table("clipboard", header=T)
Widow.dep<-read.table("clipboard", header=T)
Dep.in<-list(BlackN.dep,BlackS.dep,Blue.dep,Bocaccio.dep,CabOR.dep,CabN.dep,CabS.dep,Gstriped.dep,Widow.dep)
names(Dep.in)<-c("BlackN","BlackS","Blue","Bocaccio","CabOR","CabN","CabS","Gstriped","Widow")

#TP
par(mfrow=c(2,2))
BlackN.TGcomp.TP<-Spp.bp.plot(Dep.in$BlackN,BlackN.px.TP,Dep.targ=0.4,plot.leg=F)
BlackS.TGcomp.TP<-Spp.bp.plot(Dep.in$BlackS,BlackS.px.TP,Dep.targ=0.4,plot.leg=F)
par(mfrow=c(2,2))
Blue.TGcomp.TP<-Spp.bp.plot(Dep.in$Blue,Blue.px.TP,Dep.targ=0.4,plot.leg=F)
Bocaccio.TGcomp.TP<-Spp.bp.plot(Dep.in$Bocaccio,Boc.px.TP,Dep.targ=0.4,plot.leg=F)
par(mfrow=c(2,2))
Gstriped.TGcomp.TP<-Spp.bp.plot(Dep.in$Gstriped,Gstr.px.TP,Dep.targ=0.4,plot.leg=F)
Widow.TGcomp.TP<-Spp.bp.plot(Dep.in$Widow,Widow.px.TP,Dep.targ=0.4,plot.leg=F)
par(mfrow=c(2,2))
CabOR.TGcomp.TP<-Spp.bp.plot(Dep.in$CabOR,CabOR.px.TP,Dep.targ=0.4,plot.leg=F)
CabN.TGcomp.TP<-Spp.bp.plot(Dep.in$CabN,CabN.px.TP,Dep.targ=0.4,plot.leg=F)
CabS.TGcomp.TP<-Spp.bp.plot(Dep.in$CabS,CabS.px.TP,Dep.targ=0.4,plot.leg=F)

#LP
par(mfrow=c(2,2))
BlackN.TGcomp.LP<-Spp.bp.plot(Dep.in$BlackN,BlackN.px.LP,Dep.targ=0.25,plot.leg=F)
BlackS.TGcomp.LP<-Spp.bp.plot(Dep.in$BlackS,BlackS.px.LP,Dep.targ=0.25,plot.leg=F)
par(mfrow=c(2,2))
Blue.TGcomp.LP<-Spp.bp.plot(Dep.in$Blue,Blue.px.LP,Dep.targ=0.25,plot.leg=F)
Bocaccio.TGcomp.LP<-Spp.bp.plot(Dep.in$Bocaccio,Boc.px.LP,Dep.targ=0.25,plot.leg=F)
par(mfrow=c(2,2))
Gstriped.TGcomp.LP<-Spp.bp.plot(Dep.in$Gstriped,Gstr.px.LP,Dep.targ=0.25,plot.leg=F)
Widow.TGcomp.LP<-Spp.bp.plot(Dep.in$Widow,Widow.px.LP,Dep.targ=0.25,plot.leg=F)
par(mfrow=c(2,2))
CabOR.TGcomp.LP<-Spp.bp.plot(Dep.in$CabOR,CabOR.px.LP,Dep.targ=0.25,plot.leg=F)
CabN.TGcomp.LP<-Spp.bp.plot(Dep.in$CabN,CabN.px.LP,Dep.targ=0.25,plot.leg=F)
CabS.TGcomp.LP<-Spp.bp.plot(Dep.in$CabS,CabS.px.LP,Dep.targ=0.25,plot.leg=F)
#Get legend
CabN.TGcomp.LP<-Spp.bp.plot(Dep.in$CabN,CabN.px.LP,Dep.targ=0.25,plot.leg=T)

#Get agreement counts
BlackN.agree<-count.agree(BlackN.TGcomp.TP,perc="T")
BlackS.agree<-count.agree(BlackS.TGcomp.TP,perc="T")
Blue.agree<-count.agree(Blue.TGcomp.TP,perc="T")
Bocaccio.agree<-count.agree(Bocaccio.TGcomp.TP,perc="T")
Gstriped.TGcomp.agree<-count.agree(Gstriped.TGcomp.TP,perc="T")
Widow.agree<-count.agree(Widow.TGcomp.TP,perc="T")
CabOR.agree<-count.agree(CabOR.TGcomp.TP,perc="T")
CabN.agree<-count.agree(CabN.TGcomp.TP,perc="T")
CabS.agree<-count.agree(CabS.TGcomp.TP,perc="T")
TP.agreements<-rbind(BlackN.agree,BlackS.agree,Blue.agree,Bocaccio.agree,Gstriped.TGcomp.agree,Widow.agree,CabOR.agree,CabN.agree,CabS.agree)
colnames(TP.agreements)<-c("Dep.agree","Dep.conser","F.agree","F.conserv","N")

BlackN.agree<-count.agree(BlackN.TGcomp.LP,perc="T")
BlackS.agree<-count.agree(BlackS.TGcomp.LP,perc="T")
Blue.agree<-count.agree(Blue.TGcomp.LP,perc="T")
Bocaccio.agree<-count.agree(Bocaccio.TGcomp.LP,perc="T")
Gstriped.TGcomp.agree<-count.agree(Gstriped.TGcomp.LP,perc="T")
Widow.agree<-count.agree(Widow.TGcomp.LP,perc="T")
CabOR.agree<-count.agree(CabOR.TGcomp.LP,perc="T")
CabN.agree<-count.agree(CabN.TGcomp.LP,perc="T")
CabS.agree<-count.agree(CabS.TGcomp.LP,perc="T")
LP.agreements<-rbind(BlackN.agree,BlackS.agree,Blue.agree,Bocaccio.agree,Gstriped.TGcomp.agree,Widow.agree,CabOR.agree,CabN.agree,CabS.agree)
colnames(LP.agreements)<-c("Dep.agree","Dep.conser","F.agree","F.conserv","N")

save(VBGF.dat.spp,fitvbgf.black.g1,fitvbgf.black.g2,fitvbgf.blue,fitvbgf.boc,fitvbgf.widow.g1,fitvbgf.widow.g2,CabN.Pobjs.All,CabS.Pobjs.All,CabOR.Pobjs.All,BlackN.Pobjs.All,BlackS.Pobjs.All,Blue.Pobjs.All,
Bocaccio.Pobjs.All,Gstriped.Pobjs.All,Widow.Pobjs.All,h.spp,RP_targs_h,RP_lims_h,RP_targs_trig,RPs_h,BlackN.px.TP,BlackS.px.TP,Blue.px.TP,Boc.px.TP,CabOR.px.TP,CabN.px.TP,CabS.px.TP,Gstr.px.TP,Widow.px.TP,
BlackN.px.LP,BlackS.px.LP,Blue.px.LP,Boc.px.LP,CabOR.px.LP,CabN.px.LP,CabS.px.LP,Gstr.px.LP,Widow.px.LP,BlackN.dep,BlackS.dep,Blue.dep,Bocaccio.dep,CabOR.dep,CabN.dep,CabS.dep,Gstriped.dep,Widow.dep,Dep.in,
BlackN.TGcomp.TP,BlackS.TGcomp.TP,Blue.TGcomp.TP,Bocaccio.TGcomp.TP,Gstriped.TGcomp.TP,Widow.TGcomp.TP,CabOR.TGcomp.TP,CabN.TGcomp.TP,CabS.TGcomp.TP,BlackN.TGcomp.LP,BlackS.TGcomp.LP,Blue.TGcomp.LP,
Bocaccio.TGcomp.LP,Gstriped.TGcomp.LP,Widow.TGcomp.LP,CabOR.TGcomp.LP,CabN.TGcomp.LP,CabS.TGcomp.LP,TP.agreements,LP.agreements,file="C:/Users/copeja/Documents/Publications/In review/Rags2Fishes/DatNOutputs_2011.DMP")

load("C:/Users/copeja/Documents/Publications/In review/Rags2Fishes/DatNOutputs.DMP")
load("C:/Users/copeja/Documents/Publications/In review/Rags2Fishes/DatNOutputs_2011.DMP")

#For Rags 2 Fishes agreement, cautious, and risky categories
Spp.bp.plot.R2F()
#TP
par(mfrow=c(2,2))
BlackN.TGcomp.TP.R2F<-Spp.bp.plot.R2F(Dep.in$BlackN,BlackN.px.TP,Dep.targ=0.4,plot.leg=F)
BlackS.TGcomp.TP.R2F<-Spp.bp.plot.R2F(Dep.in$BlackS,BlackS.px.TP,Dep.targ=0.4,plot.leg=F)
par(mfrow=c(2,2))
Blue.TGcomp.TP.R2F<-Spp.bp.plot.R2F(Dep.in$Blue,Blue.px.TP,Dep.targ=0.4,plot.leg=F)
Bocaccio.TGcomp.TP.R2F<-Spp.bp.plot.R2F(Dep.in$Bocaccio,Boc.px.TP,Dep.targ=0.4,plot.leg=F)
par(mfrow=c(2,2))
Gstriped.TGcomp.TP.R2F<-Spp.bp.plot.R2F(Dep.in$Gstriped,Gstr.px.TP,Dep.targ=0.4,plot.leg=F)
Widow.TGcomp.TP.R2F<-Spp.bp.plot.R2F(Dep.in$Widow,Widow.px.TP,Dep.targ=0.4,plot.leg=F)
par(mfrow=c(2,2))
CabOR.TGcomp.TP.R2F<-Spp.bp.plot.R2F(Dep.in$CabOR,CabOR.px.TP,Dep.targ=0.4,plot.leg=F)
CabN.TGcomp.TP.R2F<-Spp.bp.plot.R2F(Dep.in$CabN,CabN.px.TP,Dep.targ=0.4,plot.leg=F)
CabS.TGcomp.TP.R2F<-Spp.bp.plot.R2F(Dep.in$CabS,CabS.px.TP,Dep.targ=0.4,plot.leg=F)

#LP
par(mfrow=c(2,2))
BlackN.TGcomp.LP.R2F<-Spp.bp.plot.R2F(Dep.in$BlackN,BlackN.px.LP,Dep.targ=0.25,plot.leg=F)
BlackS.TGcomp.LP.R2F<-Spp.bp.plot.R2F(Dep.in$BlackS,BlackS.px.LP,Dep.targ=0.25,plot.leg=F)
par(mfrow=c(2,2))
Blue.TGcomp.LP.R2F<-Spp.bp.plot.R2F(Dep.in$Blue,Blue.px.LP,Dep.targ=0.25,plot.leg=F)
Bocaccio.TGcomp.LP.R2F<-Spp.bp.plot.R2F(Dep.in$Bocaccio,Boc.px.LP,Dep.targ=0.25,plot.leg=F)
par(mfrow=c(2,2))
Gstriped.TGcomp.LP.R2F<-Spp.bp.plot.R2F(Dep.in$Gstriped,Gstr.px.LP,Dep.targ=0.25,plot.leg=F)
Widow.TGcomp.LP.R2F<-Spp.bp.plot.R2F(Dep.in$Widow,Widow.px.LP,Dep.targ=0.25,plot.leg=F)
par(mfrow=c(2,2))
CabOR.TGcomp.LP.R2F<-Spp.bp.plot.R2F(Dep.in$CabOR,CabOR.px.LP,Dep.targ=0.25,plot.leg=F)
CabN.TGcomp.LP.R2F<-Spp.bp.plot.R2F(Dep.in$CabN,CabN.px.LP,Dep.targ=0.25,plot.leg=F)
CabS.TGcomp.LP.R2F<-Spp.bp.plot.R2F(Dep.in$CabS,CabS.px.LP,Dep.targ=0.25,plot.leg=F)

Spp.bp.plot.R2F<-function(spp.dep.in,spp.px.in,Dep.targ=0.4,plot.leg=T)
{
  spp.bp.dat<-cbind(spp.dep.in,t(spp.px.in[4,]))
  spp.bp.plot<-matrix(NA,nrow=nrow(spp.bp.dat),ncol=ncol(spp.bp.dat))
  colnames(spp.bp.plot)<-c("Year","Dep","F","LtRP")
  spp.bp.plot[,1]<-spp.bp.dat[,1]
  spp.bp.plot[spp.bp.dat[,2]>=(0.8*Dep.targ)&spp.bp.dat[,2]<=Dep.targ|spp.bp.dat[,2]>=Dep.targ&spp.bp.dat[,2]<=(1.2*Dep.targ),2]<-0
  #spp.bp.plot[spp.bp.dat[,2]>=Dep.targ&spp.bp.dat[,2]<=(1.2*Dep.targ),2]<-1
  spp.bp.plot[spp.bp.dat[,2]<0.8*Dep.targ,2]<-(-1)
  spp.bp.plot[spp.bp.dat[,2]>1.2*Dep.targ,2]<-1
  spp.bp.plot[spp.bp.dat[,3]<=1,3]<-1
  spp.bp.plot[spp.bp.dat[,3]>1,3]<-(-1)
  spp.bp.plot[spp.bp.dat[,4]=="black",4]<-1
  spp.bp.plot[spp.bp.dat[,4]=="white",4]<-(-1)
  plot(spp.bp.plot[,1],spp.bp.plot[,2]*0.4,ylim=c(-1,1),xlab="Year",ylab="",type="p",pch=21,col="black",bg="gray",axes="F")
  points(spp.bp.plot[,1],spp.bp.plot[,3]*0.5,pch=22,col="black",bg="gray")
  points(spp.bp.plot[,1],spp.bp.plot[,4]*0.6,pch=24,col="black",bg="gray")
  abline(h=0)
  box()
  axis(1)
  if(plot.leg == T){smartlegend("right","bottom",c("Depletion","Fishing mortality","Lt RP"),pch=c(21,22,24),col="black",pt.bg="gray",pt.cex=1.25,bty="n")}
  return(spp.bp.plot)
}

count.ag.caut.risk<-function(Comp.in,perc="T")
{
  count.out<-c(NA,NA,NA,NA,NA,NA,0)
  count.x<-rep(1,nrow(Comp.in))
  count.out[1]<-sum((count.x[Comp.in[,2]+1)==Comp.in[,4])
  count.out[2]<-sum(count.x[Comp.in[,2]>Comp.in[,4]])
  count.out[3]<-sum(count.x[Comp.in[,2]>Comp.in[,4]])
  count.out[4]<-sum(count.x[Comp.in[,3]==Comp.in[,4]])
  count.out[5]<-sum(count.x[Comp.in[,3]>Comp.in[,4]])
  count.out[6]<-sum(count.x[Comp.in[,3]>Comp.in[,4]])
  if(perc=="T"){count.out<-count.out/nrow(Comp.in)}
  count.out[7]<-nrow(Comp.in)
  return(count.out)
}

