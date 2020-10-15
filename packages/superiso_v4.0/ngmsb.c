#include "src/include.h"
#include "src/nmssmtools.h"


#define USE_NMSSMTOOLS /* to be commented if NMSSMTOOLS is unavailable */

/*--------------------------------------------------------------------*/
/* Calculation of the observables, corresponding to an NGMSB point generated by NMSSMTools */
/*--------------------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[200];
	double Lambda,Mmess,tanb,N5,lambda,AL,sgnmu,del_h,mtop,mbot,alphas_mz,delta0;
	double obs[Nobs_BKsll+1],omega;

  	if(argc<6) 
  	{ 
    		printf(" This program needs 4 parameters:\n"
           	"   Lambda  scale of the SUSY breaking (10000-100000 GeV)\n"
           	"   Mmess   messenger mass scale > Lambda\n"
	   	"   N5      equivalent number of 5+5bar messenger fields\n"
           	"   tanb    tan(beta) \n"
          	"   lambda    Yukawa coupling\n");
    		printf(" Auxiliary parameters are:\n"
          	"   sgnmu    +/-1,  sign of Higgsino mass term (default 1)\n"   
           	"   AL       (if non-zero)\n"
           	"   Del_h    (if non-zero)\n"
             	"   mtop     top quark pole mass\n"
           	"   mbot     Mb(Mb) scale independent b-quark mass\n"
           	"   alphas_mz  strong coupling at MZ\n");                        
      		exit(1); 
  	} 
	else  
  	{  
		sscanf(argv[1],"%lf",&Lambda);
     		sscanf(argv[2],"%lf",&Mmess);
     		sscanf(argv[3],"%lf",&N5);
     		sscanf(argv[4],"%lf",&tanb);
     		sscanf(argv[5],"%lf",&lambda);
     		if(argc>6) sscanf(argv[6],"%lf",&sgnmu); else sgnmu=1;
     		if(argc>7) sscanf(argv[7],"%lf",&AL); else AL=0.;
     		if(argc>8) sscanf(argv[8],"%lf",&del_h); else del_h=0.;
     		if(argc>8) sscanf(argv[8],"%lf",&mtop); else mtop=173.34;   
     		if(argc>9) sscanf(argv[9],"%lf",&mbot); else mbot=4.18;
     		if(argc>10) sscanf(argv[10],"%lf",&alphas_mz); else alphas_mz=0.1184;
  	}	

	int filesOK=1;
#ifdef USE_NMSSMTOOLS
	if(!test_file(NMSSMTools)) 
	{
		printf("\"%s\" absent. Please check the NMSSMTOOLS path or comment \"#define USE_NMSSMTOOLS\" in ngmsb.c\n",NMSSMTools);
		filesOK=0;
	}
#endif
	if(!filesOK) return 1;
	
	if(!test_file("tmp")) system("mkdir tmp");
	chdir("tmp");
	
	printf("\n");
	
	printf("SuperIso v4.0 - F. Mahmoudi\n\n");

#ifdef USE_NMSSMTOOLS	
	sprintf(name,"ngmsb_nmssmtools%d.tmplha",getpid());
	nmssmtools_ngmsb(Lambda,Mmess,tanb,(int)N5,lambda,AL,del_h,sgnmu,mtop,mbot,alphas_mz,name);
	printf("NGMSB - SLHA file generated by NMSSMTOOLS\n\n");

	delta0=delta0_calculator(name);
	if(delta0 !=0.)
	{

		printf("Observable\t\t\tValue\n\n");


		printf("BR(b->s gamma)\t\t\t%.3e\n",bsgamma_calculator(name));
		printf("delta0(B->K* gamma)\t\t%.3e\n\n",delta0);

		printf("BR(Bs->mu mu)\t\t\t%.3e\n",Bsmumu_calculator(name));
		printf("BR(Bs->mu mu)_untag\t\t%.3e\n",Bsmumu_untag_calculator(name));
		printf("BR(Bd->mu mu)\t\t\t%.3e\n\n",Bdmumu_calculator(name));
	
		printf("BR(B->K* mu mu)_low\t\t%.3e\n",BRobs_BKstarmumu_lowq2_calculator(name,obs));
		printf("AFB(B->K* mu mu)_low\t\t%.3e\n",obs[1]);
		printf("FL(B->K* mu mu)_low\t\t%.3e\n",obs[2]);
		printf("P1(B->K* mu mu)_low\t\t%.3e\n",obs[5]);
		printf("P2(B->K* mu mu)_low\t\t%.3e\n",obs[14]);
		printf("P4'(B->K* mu mu)_low\t\t%.3e\n",obs[17]);
		printf("P5'(B->K* mu mu)_low\t\t%.3e\n",obs[18]);
		printf("P6'(B->K* mu mu)_low\t\t%.3e\n",obs[19]);
		printf("P8'(B->K* mu mu)_low\t\t%.3e\n",obs[21]);
		printf("AI(B->K* mu mu)_low\t\t%.3e\n\n",AI_BKstarmumu_lowq2_calculator(name));
	
		printf("BR(B->K* mu mu)_high\t\t%.3e\n",BRobs_BKstarmumu_highq2_calculator(name,obs));
		printf("AFB(B->K* mu mu)_high\t\t%.3e\n",obs[1]);
		printf("FL(B->K* mu mu)_high\t\t%.3e\n",obs[2]);
		printf("P1(B->K* mu mu)_high\t\t%.3e\n",obs[5]);
		printf("P2(B->K* mu mu)_high\t\t%.3e\n",obs[14]);
		printf("P4'(B->K* mu mu)_high\t\t%.3e\n",obs[17]);
		printf("P5'(B->K* mu mu)_high\t\t%.3e\n",obs[18]);
		printf("P6'(B->K* mu mu)_high\t\t%.3e\n",obs[19]);
		printf("P8'(B->K* mu mu)_high\t\t%.3e\n",obs[21]);
		printf("AI(B->K* mu mu)_high\t\t%.3e\n\n",AI_BKstarmumu_highq2_calculator(name));


		printf("BR(B->Xs mu mu)_low\t\t%.3e\n",BRBXsmumu_lowq2_calculator(name));
		printf("BR(B->Xs mu mu)_high\t\t%.3e\n",BRBXsmumu_highq2_calculator(name));
		printf("q0^2(AFB(B->Xs mu mu)\t\t%.3e\n",A_BXsmumu_zero_calculator(name));
		printf("BR(B->Xs tau tau)_high\t\t%.3e\n\n",BRBXstautau_highq2_calculator(name));
	
		printf("BR(B->tau nu)\t\t\t%.3e\n",Btaunu_calculator(name));
		printf("R(B->tau nu)\t\t\t%.3e\n",RBtaunu_calculator(name));
		printf("BR(B->D tau nu)\t\t\t%.3e\n",BDtaunu_calculator(name));
		printf("BR(B->D tau nu)/BR(B->D e nu)\t%.3e\n",BDtaunu_BDenu_calculator(name));
		printf("BR(Ds->tau nu)\t\t\t%.3e\n",Dstaunu_calculator(name));
		printf("BR(Ds->mu nu)\t\t\t%.3e\n",Dsmunu_calculator(name));
		printf("BR(D->mu nu)\t\t\t%.3e\n",Dmunu_calculator(name));
		printf("BR(K->mu nu)/BR(pi->mu nu)\t%.3e\n",Kmunu_pimunu_calculator(name));
		printf("Rmu23(K->mu nu)\t\t\t%.3e\n\n",Rmu23_calculator(name));

		printf("a_muon\t\t\t\t%.3e\n\n",muon_gm2_calculator(name));

		printf("theory_excluded\t\t\t%d\n\n",NMSSM_theory_excluded(name));

				
		flha_generator(name,"../output.flha");
		printf("output.flha generated\n\n");	
	}
	else printf("Invalid point\n\n");
	sprintf(name,"rm ngmsb_nmssmtools%d.tmplha",getpid());
	system(name);		
#endif
	return 1;
}
