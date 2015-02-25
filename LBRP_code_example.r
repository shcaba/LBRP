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


#Example
Dir<-"C:/Users/copeja/Desktop"
load(file=pastes(Dir"/RPs_in.DMP"),sep="")
load(file=pastes(Dir"/Blue_LtComps.DMP"),sep="")
Blue.ages<-c(0:60)
Blue.Linf<-38.4573
Blue.k<-0.111716
Blue.t0<-(-2.6)
Blue.M<-0.1
Blue.Lmat<-26
Blue.a.wl<-0.00003408
Blue.b.wl<-2.874
Blue.Pobjs.All<-Froese.plus(Blue.LtComps,seq(10,52,by=2),Blue.Linf,Blue.k,Blue.t0,Blue.M,Blue.a.wl,Blue.b.wl,Blue.Lmat,Blue.ages,fleet.pick=c(1:4),gender.pick=0)
#Plot relative to reference points
Blue.px.TP<-plot.Pobj(Blue.Pobjs.All,c(0,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[3],plot.leg=F)
Blue.px.LP<-plot.Pobj(Blue.Pobjs.All,c(0,1),0.1,RPs_h,h.in=T,h.spp.in=h.spp[3],RP.comp="LP",plot.leg=F)
