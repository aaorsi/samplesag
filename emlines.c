/**
 * @file emlines.c
 * @brief Here the calculation of the different emission lines, powered by 
 * star formation activity, is performed 	
 *
 * The emission-line luminosity is read from a binary file where 
 *  line-ratios are defined as a function of Z and q, the ionization 
 * parameter. 
 * @author Alvaro Orsi
 */

#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"
#include <stdio.h>
#include <strings.h>

/*
float * LinesArray;
float * ZArr;
float * QArr;
int NZ, NQ;
int NLines;
char *LineName;
*/
	
void get_emlinesinfo(char *info_file)
{
	FILE * ff;
	char buf[SED_MAXC],info[SED_MAXC];
	char *p;
	int id,j; 
	float cw;

	if (!(ff = fopen(info_file,"r")))
	{
		printf("DUMPHDF5_SED.C: Emission lines info file %s not found. Check!\n",info_file);
		exit(0);
	} 		

	j = 0;

	// The first 2 lines are header+blank space, so they are irrelevant here.
	// If by any chance this changes, then this should be modified accordingly.

	fgets(buf,SED_MAXC,ff);
	fgets(buf,SED_MAXC,ff);

	while(!feof(ff))
	{			
		fscanf(ff,"%d %s \t \t \t %f",&id,Emlines_name[j],&cw);
		j++;
	}	
	fclose(ff);

	printf("Emission Lines to be computed:\n");
	for (j=0;j<NLines;j++)
	{
		printf("%s, Central wavelength= %f\n",Emlines_name[j],CentralWavelength[j]);
	}
}		



void read_emlines(char *emlines_file)
{
	
	int i,j,k;
	int NTot;
	int NStr;
	FILE *fel;

//	Hdisk = 1.686;
//	hbulge= 1.350;
	
	if (! (fel = fopen(emlines_file,"r")))
	{
		printf(" Cannot open file %s. Check!\n",emlines_file);
		exit(0);
	}

	printf("Reading emission lines from file %s\n",emlines_file);	
	
	
	fread(&NLines,sizeof(int),1,fel);
	
	for (i = 0; i<NLines; i++)
		fread(&CentralWavelength[i],sizeof(float),1,fel);

	fread(&SED_NZEL,sizeof(int),1,fel);	
	fread(&SED_NQEL,sizeof(int),1,fel);
	fread(&SED_NAEL,sizeof(int),1,fel);

	printf("NLines %d SED_NZEL %d SED_NQEL %d SED_NAEL %d\n",NLines,SED_NZEL,SED_NQEL,SED_NAEL);

	for (i = 0;i<SED_NZEL;i++)
		fread(&ZArr[i],sizeof(float),1,fel);

	for (i = 0;i<SED_NQEL;i++)
		fread(&QArr[i],sizeof(float),1,fel);

	for (i = 0;i<SED_NAEL;i++)
		fread(&Age_emlines_Arr[i],sizeof(float),1,fel);

	NTot = SED_NZEL * SED_NQEL * SED_NAEL;

	for (i = 0;i<SED_NZEL; i++)
		printf("ZArr[%d] = %f\n",i,ZArr[i]);

	for (i = 0;i<SED_NQEL; i++)
		printf("QArr[%d] = %f\n",i,QArr[i]);

	for (i = 0;i<SED_NAEL;i++)		
		printf("Age_emlines[%d] = %f\n",i,Age_emlines_Arr[i]);

	k = 0;
	for (i = 0; i< NTot; i++)
	{
		for (j = 0; j< NLines; j++)		
		{
			fread(&LinesArray[k],sizeof(double),1,fel);
			k++;
		}
	}
	
	fclose(fel);

}


	

double calc_emlines(float QGas, float ZGas, float age, int idl)
{
	
	// QGas corresponds to the Q parameter of the galaxy
	// ZGas is the cold gas metallicity
    // line is the pointer to the array where the line intensities are stored.

	int i,j,iz,iq;
	float Z0,Z1, Q0,Q1,A0,A1;
	int iz0,iz1,iq0,iq1,ia0,ia1;
	int ia;
	float f00,f10,f01,f11;
	float t1,t2,t3,t4;
	float den;
	float line;

	age = 0;	
	ia = 0;

// Bilinear interpolation
// (age is approximated to its nearest value)
/*
	if (age < Age_emlines_Arr[0])
	{
		printf("calc_emlines(): Age of SSP %f lower than minimum age of photoionization model %f. Check!\n",age,Age_emlines_Arr[0]);
		exit(0);
	}


	i = 0;
	while (age >= Age_emlines_Arr[i] && i < SED_NAEL-1)
	{
		A0 = Age_emlines_Arr[i];
		ia0 = i;
		i++;
	}
	A1 = Age_emlines_Arr[i];
	ia1 = i;
	ia = (fabs(age-A0) < fabs(age-A1)) ? ia0 : ia1;	

*/

	if (QGas < QArr[0]) // || ZGas < ZArr[0])	
	{
		line = 0.0;
		return;
	}

	if (ZGas < ZArr[0])
	{
		iz0	= 0;
		iz1 = 1;
		Z0 	= ZArr[iz0];
		Z1	= ZArr[iz1];
	}
	else
	{
		i = 0;
		while(ZGas >= ZArr[i] && i < SED_NZEL-1)
		{
			Z0 = ZArr[i];
			iz0 = i;
			i++;
		}	
		iz1 = i;
		Z1 = ZArr[iz1];
	}

	if (QGas < QArr[0])
	{
		iq0	= 0;
		iq1	= 1;
		Q0	= QArr[0];		
		Q1	= QArr[1];
	}
	else
	{
		i = 0;
		while(QGas >= QArr[i] && i < SED_NQEL-1)
		{
			Q0 = QArr[i];	
			iq0 = i;
			i++;
		}	
	
		Q1 = QArr[i];
		iq1 = i;
	}

//	for (i = 0;i<NLines;i++)
//	{
/*	f00 = LinesArray[idl + NLines*iz0 + SED_NZEL*NLines*iq0 + NLines*SED_NZEL*SED_NQEL*ia];
	f10 = LinesArray[idl + NLines*iz1 + SED_NZEL*NLines*iq0 + NLines*SED_NZEL*SED_NQEL*ia];
	f01 = LinesArray[idl + NLines*iz0 + SED_NZEL*NLines*iq1 + NLines*SED_NZEL*SED_NQEL*ia];
	f11 = LinesArray[idl + NLines*iz1 + SED_NZEL*NLines*iq1 + NLines*SED_NZEL*SED_NQEL*ia];
*/
	f00 = LinesArray[idl + NLines*iz0 + SED_NZEL*NLines*iq0];
	f10 = LinesArray[idl + NLines*iz1 + SED_NZEL*NLines*iq0];
	f01 = LinesArray[idl + NLines*iz0 + SED_NZEL*NLines*iq1];
	f11 = LinesArray[idl + NLines*iz1 + SED_NZEL*NLines*iq1];

	den = (Z1-Z0)*(Q1-Q0);
				
	t1 = f00*(Z1-ZGas)*(Q1-QGas)/den;
	t2 = f10*(ZGas-Z0)*(Q1-QGas)/den;
	t3 = f01*(Z1-ZGas)*(QGas-Q0)/den;
	t4 = f11*(ZGas-Z0)*(QGas-Q0)/den;

	line = t1 + t2 + t3 + t4;
	line = (line <= 0.0) ? 1e-30 : line;

	return(line);

//	return (line);

//	}	

}	


double integ_line(float ZGal, double nlyc, int iline, int ised)
{

/* This routine should provide different ways of computing the ionization parameter.
   From which then it's trivial to compute the resulting, integrated emission line intensity. 

   To date (27/06/12) the simplest one is implemented: Spreading the ionizing flux
   over an exponential surface density profile, and, by assuming constant number density, 
   integrating for different radii, assuming that all the ionizing flux 
   is concentrated within the half-mass radius of the disk. 

   A similar approach should be implemented for the bulge.
*/


/*	LineMod = 			
								0 q = q0, scaling relation
		   					    1 q(r) and an exponential disk profile 
							    2 q(M), with a molecular cloud power-law distribution 
*/

	float rmin,rmax,r_i,rbin,Mbin;
	int i,Nrsteps,iHalpha,Nmsteps;
	double line,line_i,q_r,m_i,q_m,N_m,L_m;
	double fline,FracHyda;
	float mq;
	float hconst;
	float nconst;
	float Zscale,rscaleZ;
	float mgas_r,mmet_r; 
	double alpha,beta,betaQ,eta,kappa,gamm,Mmax,Mmin,Mtot,Ntot;
	double log_qconst,qconst,logMmin,logMmax,nHconst,nH2_i,nlyc_i;	
	double eps_ff;
	double K,M_,Gamm,QMax,QMin;
	double nlyc_40,halphaK,pow40_third;
	double cel;

	float eps_k,eps_g;
	float h2_h1;

	int el_model;

	
	line = -30;
	mq = 0.0;

	if (isnan(nlyc))
	{
		Mean_q = 0.0;	
		return(line);
	}


	if (ised == 1)
		el_model = LineMod;
	if (ised ==2)
		el_model = BulgeMod;

		QMax	= QArr[SED_NQEL-1];
		QMin    = QArr[0];	

	if (el_model ==0)
	{
		iHalpha = 2;
		alpha = log10(1.37)-12;
//		Mean_q = 5.e7;
		Mean_q = Qval;
		FracHyda   = log10(calc_emlines(Mean_q,ZGal,0,2));
		cel = calc_emlines(Mean_q,ZGal,0,iline);
		line = (cel < 1e-20) ? -30 : nlyc + alpha + log10(cel) - FracHyda -40.0;
	}
	
	if (el_model == 1)
	{
//	Q-M model
		iHalpha = 2;
		alpha = log10(1.37)-12;


		if (Qvar == 1)
		{

//			K 		= 7.5e6;	
//			M_ 		= 1.7;
//			Gamm	= -0.5; 

			K		= Kqz;
			M_		= Zqz;
			Gamm	= Gqz;
		
			if ( mstellar < 1e-3)
				Mean_q = QMax;
			else if (mstellar > 1.0)
				Mean_q = QMin;
			else
				Mean_q = K * pow(mstellar/M_,Gamm);
		}

// Q-Z Model
		if (Qvar ==2) 
		{
//			K 		= 9.42e6;
//			M_ 		= 0.0066;
//			Gamm	= -1.06;

			K		= Kqz;
			M_		= Zqz;
			Gamm	= Gqz;
		
//			if ( log10(ZGal) < -3.7)
//				Mean_q = QMax;
//			else if (log10(ZGal) > -2.2)
//				Mean_q = QMin;
//			else

			Mean_q = K * pow(ZGal/M_,Gamm);
			if (Mean_q < QMin)
				Mean_q = QMin;
			if (Mean_q > QMax)
				Mean_q = QMax;
		}

// Fit to Nagao+06 data
	
		if (Qvar ==3)
		{
			K 		= 2.0e7;
			M_ 		= 0.0054;
			Gamm	= -0.4;
		
/*			if ( ZGal < 0.0013)
				Mean_q = 7.5e7;
			else if (ZGal > 0.04)
				Mean_q = 1.9e7;
			else
*/
				Mean_q = K * pow(ZGal/M_,Gamm);
				if (Mean_q > QMax)
					Mean_q = QMax;
				if (Mean_q < QMin)
					Mean_q = QMin;
		}
	
			ZGal *= pow(10,Zoff);
			FracHyda   = calc_emlines(Mean_q,ZGal,0,2);

			cel = calc_emlines(Mean_q,ZGal,0,iline);

			line = (cel <= 1e-29 || FracHyda <= 1e-29) ? -30 : nlyc + alpha + log10(cel) - log10(FracHyda) -40.0;
			


	}
			
/*	
	OBSOLETE MODEL, UNUSED BUT LEFT HERE ONLY IN CASE IS USEFUL IN THE FUTURE.
	if (LineMod == 1) 
	{

		hconst =  Hdisk;
 		rscale_phy = log10(rscale * UnitLength_in_cm);
//		Sqr_radius = pow(hconst*rscale_phy,2);
		Sqr_radius = 2*(rscale_phy + log10(hconst));

		rmin = -2.0;
		rmax = 1.0;

		rmin *= 1./(log10(exp(1))); //It's more intuitive to express this in base 10, and then convert to base e.
		rmax *= 1./(log10(exp(1)));

		Nrsteps = 10;	
	
		rbin = (rmax - rmin)/(Nrsteps-1);
		Sig0 = -(log10(2*PI)+2*rscale_phy);

		nconst = log10(10.0);

// Radial profile	
//		printf("Line: %s\n",Emlines_name[iline]);
//		printf("ZGal = %f\n******\n",ZGal);
		for (i= 0;i<Nrsteps;i++)
		{

			r_i  = exp(i*rbin + rmin);	
			hc = hconst*r_i;
			q_r  = pow(10,nlyc-nconst + Sig0) * exp(-hc);
			line += r_i*r_i*calc_emlines(q_r,ZGal,age,iline);
//			printf("Line_%d = %g r_i %g q_r  %g\n",i,line,r_i,q_r);
		}

		if (line <= 0)
			line = -30.0;
		else
			line  = log10(line)+log10(2*PI*rbin) + Sqr_radius - 40.0;

		Mean_q = pow(10,(nlyc-nconst)+Sig0) * exp(-Hdisk);
		
//		printf("line tot = %g\n",line);

	}
*/

	if (el_model == 2) 
	{	
		if (mgas <= 0.0)
		{	
			line = -30;
			Mean_q = 0.0;
			return(line);
		}

		eps_ff		= 0.1;  		// Cloud's filling factor (-1 means use scaling relation).
		beta		= BetaGMC;		// Relation between H2/H1
		betaQ		= 1.0;			// Exponent of M-Q relation

		if (BetaGMC == -1.0)		// This activates Lagos+11 scaling relations
		{
			h2_h1 		= 0.09 * pow(mgas,0.9)*pow((1+ZZ[Snapshot]),2.4);// Following Eq. 7 from Lagos et al. (2011)	
			beta		= h2_h1/(1 + h2_h1);
		}

		if (beta > 1.0)
			beta = 1.0;

		alpha		= -2.0;			// Exponent of dn/dM of GMCs
		Mmax 		= 3.0e7;		// both in Solar mass 
		Mmin 		= 500.0; 
		eta	 		= Eta_GMC/(PC2CM*PC2CM);		// M_{sun}/pc^2
		Mtot  		= FracHyd * beta* mgas* 1.0e10 /Hubble_h;	//mgas is in units of M_sun
		kappa 		= Mtot / (log(Mmax/Mmin));
		Ntot 		= kappa/(alpha+1) * (pow(Mmax,alpha+1) - pow(Mmin,alpha+1));

		Nmsteps = 10;
		logMmin = log(Mmin);
		logMmax = log(Mmax);
		Mbin    = (logMmax - logMmin)/(Nmsteps-1);

		nHconst = (3 * pow(eta,1.5) / (4*PI*MassHyd)) * SunMass;
		qconst  = (3./4) * pow(AlphaB,2./3) * pow(3/(4*PI),1./3);
		// Check out the units of eta!
	
		nlyc_40 = pow(10,nlyc-40.0);
		pow40_third = pow(1e40,1./3);

		halphaK = 1.37e-12;
		line_i  = 0.0;

		for (i = 0;i<Nmsteps; i++)
		{
			m_i   		= exp(i*Mbin + logMmin);
			nH2_i   	= nHconst * pow(m_i,-0.5);	
			nlyc_i 		= nlyc_40 * ( m_i / Mtot);
			q_m			= pow(eps_ff,2./3) * qconst * pow( nlyc_i * nH2_i,1./3) * pow40_third;
			if (q_m > QMax)
				q_m = QMax;
			if (q_m < QMin)
				continue;

			fline 		= calc_emlines(q_m,ZGal,0,iline);
//			N_m   = kappa * pow(m_i,alpha);

			FracHyda   = calc_emlines(q_m,ZGal,0,2);
			
			L_m	  = ( fline <= 1e-20 || FracHyda <= 1e-20 ) ? 0.0 : \
					1e40 * nlyc_i * halphaK  * fline / FracHyda;

			line_i += pow(m_i,alpha+1)*L_m;
			mq   += pow(m_i,alpha+1)*q_m;

			if (i ==0) 
				Qlow = q_m;
		}
		
		Qhigh = q_m;

		if (mq > 0.0) 
		{
			Mean_q = kappa * mq * Mbin/Ntot;
			line = log10(kappa*line_i*Mbin * (Hubble_h*Hubble_h)) - 40.0; // [correct mass units]
		} 
		else 
		{	
			line 	= -30.;
			Mean_q	= 0.;
		}

	}
	return(line);
}
