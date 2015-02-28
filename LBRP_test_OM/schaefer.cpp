#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <schaefer.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  n.allocate("n");
  d.allocate(1,n,1,3,"d");
  yr.allocate(1,n);
  yt.allocate(1,n);
  ct.allocate(1,n);
		yr=ivector(column(d,1));
		yt=column(d,2);
		ct=column(d,3);
		cout<<"Finished reading data\n"<<endl;
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  msy.allocate(0,2.,"msy");
  fmsy.allocate(0,1,"fmsy");
  tau.allocate(0,100,"tau");
  f.allocate("f");
  k.allocate("k");
  #ifndef NO_AD_INITIALIZE
  k.initialize();
  #endif
  r.allocate("r");
  #ifndef NO_AD_INITIALIZE
  r.initialize();
  #endif
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  pt.allocate(1,n+1,"pt");
  #ifndef NO_AD_INITIALIZE
    pt.initialize();
  #endif
}

void model_parameters::userfunction(void)
{
	schaefer_model();
	
}

void model_parameters::schaefer_model(void)
{
	r=2.0*fmsy;
	k=4.0*msy/r;
	pt(1)=1.0; fpen=0;
	for(int i=1;i<=n;i++)
	{
		pt(i+1)=posfun(pt(i)+r*pt(i)*(1.-pt(i))-ct(i)/k,1.e-3,fpen);
	}
	dvar_vector zt=log(yt)-log(pt(1,n)*k);
	zt=zt-mean(zt);
	f=0.5*size_count(zt)*log(1/tau)+0.5*tau*norm2(zt);
	f+=0.5*log(0.25)+log(msy)+1./(2.*0.25)*square(log(msy)-log(mean(ct)));
	f+= -(1.1-1.)*log(fmsy) - (1.9-1.)*log(1.-fmsy);
	f+=1000.*fpen;
}

void model_parameters::mpe_output()
{
	ofstream ofs("Schaefer.rep");
	ofs<<"#fmsy\n"<<fmsy<<endl;
	ofs<<"#Bt\n"<<pt*k<<endl;
	ofs<<"#Bo\n"<<k<<endl;
	
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
		mpe_output();
		
		
}

void model_parameters::preliminary_calculations(void){
  admaster_slave_variable_interface(*this);
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
  #if defined(__GNUDOS__) || defined(DOS386) || defined(__DPMI32__)  || \
     defined(__MSVC32__)
      if (!arrmblsize) arrmblsize=150000;
  #else
      if (!arrmblsize) arrmblsize=25000;
  #endif
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    // so we can stop here
    exit(i);
  }
}
