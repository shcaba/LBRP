//******************************************************
//	Programmer: Steve Martell
//	Project Name:Schaefer Model
//	Date:
//	Version:
//	Comments:
//	
//******************************************************/
DATA_SECTION
	init_int n;
	init_matrix d(1,n,1,3);
	ivector yr(1,n);
	vector yt(1,n);
	vector ct(1,n);
	LOC_CALCS
		yr=ivector(column(d,1));
		yt=column(d,2);
		ct=column(d,3);
		cout<<"Finished reading data\n"<<endl;
	END_CALCS
	
	
PARAMETER_SECTION
	init_bounded_number msy(0,2.);
	init_bounded_number fmsy(0,1);
	init_bounded_number tau(0,100);
	
	objective_function_value f;

	number k;
	number r;
	number fpen;
	vector pt(1,n+1);

PROCEDURE_SECTION
	schaefer_model();

	

FUNCTION schaefer_model
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

FUNCTION void mpe_output()
	ofstream ofs("Schaefer.rep");
	ofs<<"#fmsy\n"<<fmsy<<endl;
	ofs<<"#Bt\n"<<pt*k<<endl;
	ofs<<"#Bo\n"<<k<<endl;
	

REPORT_SECTION
		mpe_output();
		
		

