//******************************************************
//	Programmer: Steve Martell
//	Project Name: LRGS Model
//	Date:
//	Version:
//	Comments: Uses Fmsy and MSY as leading parameters
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
	init_bounded_number msy(0.01,2.0);
	init_bounded_number fmsy(0.01,0.95);
	init_bounded_number s(0.4,0.99);
	init_bounded_number tau(0,100);
	!! fmsy=0.2;
	!! s = 0.8;
	objective_function_value f;

	number bo;
	number a;
	number b;
	number fpen;
	
	vector bt(1,n+1);
	
PROCEDURE_SECTION
	lrgs_model();
	
FUNCTION lrgs_model
	fpen=0;
	dvariable et;
	a = square((1.0-s+fmsy))/(1.0-s);
	b = posfun(fmsy*(a/(1.0-s+fmsy)-1.0)/msy,1.e-5,fpen);
	bo = (1-a-s)/(b*(s-1)); 

	bt(1)=bo;
	for(int i=1;i<=n;i++)
	{
		if(i>4) et=bt(i-3); else et=bo;
		bt(i+1)=posfun(s*bt(i)+a*et/(1+b*et) - ct(i),1.e-3,fpen);
	}
	dvar_vector zt=log(yt)-log(bt(1,n));
	
	zt=zt-mean(zt);
	f=0.5*size_count(zt)*log(1/tau)+0.5*tau*norm2(zt);
	f+=0.5*log(0.5)+log(msy)+1./(2.*0.5)*square(log(msy)-log(mean(ct)));
	f+= -(1.1-1.)*log(fmsy) - (1.9-1.)*log(1.-fmsy); //beta prior for fmsy
	f+=1000.*fpen;

	
FUNCTION void mpe_output()
	ofstream ofs("LRGS.rep");
	ofs<<"#fmsy\n"<<fmsy<<endl;
	ofs<<"#Bt\n"<<bt<<endl;
	ofs<<"#Bo\n"<<bo<<endl;

	//Algebra check, cmsy should correspond to msy
	//Checks out OK;
	//double bmsy =value( (1-s-sqrt(-a*(s-1)))/(b*(s-1)));
	//double cmsy = value(bmsy*(s-1)+a*bmsy/(1+b*bmsy));
	//cout<<cmsy<<endl;

REPORT_SECTION
	mpe_output();
	


