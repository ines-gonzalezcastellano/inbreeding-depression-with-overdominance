/* SIMOVERDOM_NEW.c (19/07/2024) */

/* ***************************************************** */

#include "libhdr"
#define NN 10001  /* max number of NIND */
#define MM 8001  /* max number of NCRO */
#define SS 100000  /* max number of SNPs */
#define GG 4000  /* max number of genes */

int NIND, NCRO, NLOCI, TOTLOCI, numSNPs, numSNPsneu, numSNPsdel, numSNPslet, numSNPsod, numSNPNP;
int i, j, k, l, rep, classes, replicates, g, b;
int RM[NN], ran_i;
int gm[NN][MM][2];
int ss, cc, crom[SS], loc[SS], genes, bb;

int NINDNP, gmNP[NN][MM][2];
int chromNP[MM][31], chrom[MM][31];
unsigned long long int x, posNP[MM][31], pos[MM][31], posplink[SS];
double sNP[MM][31], atNP[MM][31], hsNP[MM][31], hatNP[MM][31], qNP[MM][31];
double s[MM][31], s1[MM][31], s2[MM][31], hs[MM][31], at[MM][31], hat[MM][31], freqNP[MM][31], LEQ_NP;
double w, genvalue[NN], FhatI[NN], FhatII[NN], FhatIII[NN], Fhom[NN], Fexh[NN], Froh100[NN], Froh1000[NN], Froh5000[NN];

double pm_s[NN], MAF;

double d_a_del, alfa_a_del, va_del, vd_del, id_del;
double d_a_let, alfa_a_let, va_let, vd_let, id_let;
double d_a_od, alfa_a_od, va_od, vd_od, id_od;

double AA, Aa, aa, q[MM][31];
double FITN[NN], FITN_F[NN];
double mean_FITN, mean_FITN_F;

double numLOCI, EHom, EHet, Hom[NN], H_AA[NN], H_Aa[NN], H_aa[NN], sum_VR1[NN], sum_VR2[NN], sum_YANG1[NN], sum_YANG2[NN], sum_yang2_delta[NN], sum_LH2[NN];
double fvr1[NN], fvr2[NN], fyang1[NN], fyang2[NN], fyang2_delta[NN], fLH1[NN], fLH2[NN], fhom[NN];
double IDFvr1, IDFyang1, SqE_Fvr1, SqE_Fyang1;

double IDFhatI, IDFhatII, IDFhatIII, IDFhom, IDFexh, IDFroh100, IDFroh1000, IDFroh5000;
double SqE_FhatI, SqE_FhatII, SqE_FhatIII, SqE_Fhom, SqE_Fexh, SqE_Froh100, SqE_Froh1000, SqE_Froh5000;

struct acc gmean_s, gvar_s, AVE_VA, AVE_VD, AVE_ID, AVE_VA_del, AVE_VD_del, AVE_ID_del, AVE_VA_let, AVE_VD_let, AVE_ID_let, AVE_VA_od, AVE_VD_od, AVE_ID_od;
struct acc AVE_numSNPsneu[SS], AVE_numSNPsdel[SS], AVE_numSNPslet[SS], AVE_numSNPsod[SS];

struct acc q_0[NN];
struct acc AVE_q_0[NN];

struct acc qdel[SS], qlet[SS], qod[SS], qneu[SS];

struct acc AVE_mean_FITN, AVE_mean_FITN_F, AVE_ID_FITN;

struct acc AFhatI, AFhatII, AFhatIII, AFhom, AFexh, AFroh100, AFroh1000, AFroh5000, Phe;
struct acc AVE_IDFhatI, AVE_IDFhatII, AVE_IDFhatIII, AVE_IDFhom, AVE_IDFexh, AVE_IDFroh100, AVE_IDFroh1000, AVE_IDFroh5000;
struct acc AVE_SqE_FhatI, AVE_SqE_FhatII, AVE_SqE_FhatIII, AVE_SqE_Fhom, AVE_SqE_Fexh, AVE_SqE_Froh100, AVE_SqE_Froh1000, AVE_SqE_Froh5000;

struct acc AFvr1, AFyang1;
struct acc AVE_IDFvr1, AVE_IDFyang1;
struct acc AVE_SqE_Fvr1, AVE_SqE_Fyang1;
struct covacc PheFvr1, PheFyang1;

struct covacc PheFhatI, PheFhatII, PheFhatIII, PheFhom, PheFexh, PheFroh100, PheFroh1000, PheFroh5000;
struct covacc FhatIFhatII, FhatIFhatIII, FhatIFhom, FhatIFexh, FhatIFroh100, FhatIFroh1000, FhatIFroh5000;
struct covacc FhatIIFhatIII, FhatIIFhom, FhatIIFexh, FhatIIFroh100, FhatIIFroh1000, FhatIIFroh5000;
struct covacc FhatIIIFhom, FhatIIIFexh, FhatIIIFroh100, FhatIIIFroh1000, FhatIIIFroh5000;
struct covacc FhomFexh, FhomFroh100, FhomFroh1000, FhomFroh5000;
struct covacc FexhFroh100, FexhFroh1000, FexhFroh5000;
struct covacc Froh100Froh1000, Froh100Froh5000;
struct covacc Froh1000Froh5000;

FILE *fptr, *fgen, *frep, *fdat, *fpop, *ffmap, *ffreq, *ffped, *ffPphen, *fFfile, *fsumID, *fsumM, *fsumV, *fsumR, *fsumSqE, *fsumB, *fvaluesF;

/* ***************************************************** */

main()
{
	fptr = fopen ("dfilename.dat","w");
	fgen = fopen ("genfile.dat","w");
	ffreq = fopen ("frequencies.dat","w");
	fsumID = fopen ("summaryID","w");
	fsumM = fopen ("summaryM","w");
	fsumV = fopen ("summaryV","w");
	fsumSqE = fopen ("summarySqE","w");
	fsumR = fopen ("summaryR","w");
	fsumB = fopen ("summaryB","w");
	fvaluesF = fopen ("valuesF","w");

	getinputs();
	recombination_masks();
	natural_population();

	fprintf(ffreq, "0.0 ");
	for (j=1; j<classes; j++)	fprintf(ffreq, "%4.2f ", (double)j/100);

	for (rep=1; rep<=replicates; rep++)
	{
		numSNPs = 0.0;

		frep = fopen ("repfile.dat","a");
		fprintf (frep,"\nreplicate %d\n", rep);
		fclose(frep);

		fprintf(ffreq, "\n");

	//	if (tracelevel!=0) fprintf (fptr,"\n\nreplicate %d\n\n", rep);

		sample();
		frequency_genes();
		genotypic_values();
		phenotypeB();
//		if (tracelevel!=0) dumpphenotypes();
		plink_files();
		int status1 = system("bash shell_F_values");
		estimates_of_F();
		regression_pm_F();
		ID_HOMOZYGOTES();
		settozero();
	 	printout();
		writeseed();
  	}
}

/* ***************************************************** */

getinputs()
{
	tracestart();
	getseed();

	getintandskip("NINDNP (max 10000):",&NINDNP,2,10000);
	getintandskip("NIND (max 1000):",&NIND,2,10000);
	NLOCI=30;
	getintandskip("Number of frequeny classes :",&classes,1,1000);
	getrealandskip("MAF :",&MAF,0.0,1.0);
	getintandskip("Number of replicates :",&replicates,1,infinity);
}

/* **************************************************** */

recombination_masks ()
{
	for (l=0; l<NLOCI; l++)   RM[l]=pow(2.0,(double)l);
}

/* ***************************************************** */

natural_population ()
{
	int x, ds, g0, g1;
	unsigned long long int dpos;

	double dps, da, dh, dq;

	/* ***** take effects of genes ***** */

	fdat=fopen("list_allsnps","r");

	fscanf(fdat,"%d", &x);
	numSNPNP = x;

	NCRO = numSNPNP/NLOCI;
	TOTLOCI = NCRO * NLOCI;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		fscanf(fdat,"%d%llu%lf%lf%lf%lf", &ds, &dpos, &dps, &da, &dh, &dq);
		chromNP[k][l] = ds;
		//fprintf(fptr,"chromNP[%d][%d]=%d\n", k, l, chromNP[k][l]);
		posNP[k][l] = dpos;
		if (dps < -1.0) dps=(-1.0);
		if (da == -99.0) da=0.0;
		sNP[k][l] = dps;
		atNP[k][l] = da;
		hsNP[k][l] = dh;
		hatNP[k][l] = dh;
		freqNP[k][l] = dq;
//		if((tracelevel!=0)&&(k==0)&&(l==0)) fprintf(fptr,"k=%d l=%d posNP=%llu chromNP=%d sNP=%f hNP=%f freqNP=%f\n", k, l, posNP[k][l], chromNP[k][l], sNP[k][l], hsNP[k][l], freqNP[k][l]);
	}

	/* ***** take genotypic values of natural population ***** */

	fpop=fopen("dataBP.ped","r");

	for (i=0; i<NINDNP; i++)
	{
		lookfortext("IND");

		for (j=1; j<=5; j++)	fscanf(fpop,"%d", &x);

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		{
			fscanf(fpop,"%d%d", &g0, &g1);

			if (g0 == 2)	gmNP[i][k][0]=(gmNP[i][k][0] | RM[l]);
			if (g1 == 2)	gmNP[i][k][1]=(gmNP[i][k][1] | RM[l]);
		}
	}

	fclose(fpop);

//	if (tracelevel!=0)
	{
		fprintf(fptr, "\nNatural population genotypes\n");
		for (i=0; i<NIND; i++)
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((i==0)&&(k==0))
		{
			if ((gmNP[i][k][0] & RM[l])==RM[l])	fprintf(fptr, "1 ");
			else								fprintf(fptr, "0 ");
			if ((gmNP[i][k][1] & RM[l])==RM[l])	fprintf(fptr, "1 ");
			else								fprintf(fptr, "0 ");
		}
		fprintf(fptr, "\n");
	}

	/* ***** estimate LEQ in the natural population ***** */

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NINDNP; i++)
		{
			if (((gmNP[i][k][0] & RM[l])==RM[l])&&((gmNP[i][k][1] & RM[l])==RM[l]))		aa+=1.0;
	    		else if (((gmNP[i][k][0] & RM[l])!=RM[l])&&((gmNP[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     else	Aa+=1.0;
		}

		qNP[k][l] = (aa/(double)NINDNP)+(Aa/(2.0*(double)NINDNP));

//		if ((tracelevel!=0)&&(k==0))	fprintf(fptr, "k=%d l=%d AA=%f Aa=%f aa=%f q=%f\n", k, l, AA, Aa, aa, qNP[k][l]);

		LEQ_NP += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);
	}

//	if (tracelevel!=0)
	{
		fprintf(fptr, "\n LEQ_NP = %f\n", LEQ_NP);
//		for (k=0; k<NCRO; k++)
//		for (l=0; l<NLOCI; l++)
//		if (k==0)
//		fprintf(fptr, "\n k=%d l=%d   sNP=%f  hsNP=%f  qNP=%f", k, l, sNP[k][l], hsNP[k][l], qNP[k][l]);
	}

	fclose(fdat);
}

/* ***************************************************** */

sample ()
{
	int g;

	/* ***** sample the first NIND individuals from the Base Population ***** */

	/* ** Randomise NINDNP individuals ** */

	for (i=0; i<NINDNP; i++)
	{
		ran_i = (int)(uniform() * NINDNP);

		for (k=0; k<NCRO; k++)
		{
			g=gmNP[i][k][0]; gmNP[i][k][0]=gmNP[ran_i][k][0]; gmNP[ran_i][k][0]=g;
			g=gmNP[i][k][1]; gmNP[i][k][1]=gmNP[ran_i][k][1]; gmNP[ran_i][k][1]=g;

//			if ((tracelevel!=0)&&(i==0)&&(k==0))	fprintf (fptr," i = %d  k = %d  gmNP0 = %d   gmNP1 = %d\n", i, k, gmNP[i][k][0], gmNP[i][k][1]);
		}
	}

//	if (tracelevel!=0)
	{
		fprintf(fptr, "\nFirst individual natural population\n");
		for (l=0; l<NLOCI; l++)
		{
			if ((gmNP[0][0][0] & RM[l])==RM[l])	fprintf(fptr, "1 ");
			else								fprintf(fptr, "0 ");
			if ((gmNP[0][0][1] & RM[l])==RM[l])	fprintf(fptr, "1 ");
			else								fprintf(fptr, "0 ");
		}
		fprintf(fptr, "\n");
	}

	for (i=0; i<NIND; i++)
	for (k=0; k<NCRO; k++)
	{
		gm[i][k][0]=gmNP[i][k][0];
		gm[i][k][1]=gmNP[i][k][1];

//		if ((tracelevel!=0)&&(i==0)&&(k==0))	fprintf (fptr," i = %d  k = %d  gm0 = %d   gm1 = %d\n", i, k, gm[i][k][0], gm[i][k][1]);
	}

//	if (tracelevel!=0)
	{
		fprintf(fptr, "\nFirst individual sampled population\n");
		for (l=0; l<NLOCI; l++)
		{
			if ((gm[0][0][0] & RM[l])==RM[l])	fprintf(fptr, "1 ");
			else							fprintf(fptr, "0 ");
			if ((gm[0][0][1] & RM[l])==RM[l])	fprintf(fptr, "1 ");
			else							fprintf(fptr, "0 ");
		}
		fprintf(fptr, "\n");
	}

	/* ***** take effects of genes from Base Population ***** */

	if (rep == 1)
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		chrom[k][l] = chromNP[k][l];
		//fprintf(fptr,"chrom[%d][%d]=%d\n", k, l, chrom[k][l]);
		pos[k][l] = posNP[k][l];
		s[k][l] = sNP[k][l];
		at[k][l] = atNP[k][l];
		hs[k][l] = hsNP[k][l];
		hat[k][l] = hatNP[k][l];
//		if((k==4000)&&(l==0)) fprintf(fptr,"\n\nk=%d l=%d %llu chrom=%d\n\n", k, l, pos[k][l], chrom[k][l]);

//		if (tracelevel!=0)    if ((k<4)&&(l<5)) fprintf(fptr,"\n k=%d l=%d s=%f hs=%f", k, l, s[k][l], hs[k][l]);
	}
}

/* ***************************************************** */

frequency_genes ()
{
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NIND; i++)
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		aa+=1.0;
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
			else	Aa+=1.0;
		}

		q[k][l] = (aa/(double)NIND)+(Aa/(2.0*(double)NIND));

//		if ((tracelevel!=0)&&(k==0))	fprintf(fptr, "k=%d l=%d s=%f  hs=%f  AA=%f Aa=%f aa=%f q=%f\n", k, l, s[k][l], hs[k][l], AA, Aa, aa, q[k][l]);

		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(q[k][l] >= MAF)&&(s[k][l]==0.0)) numSNPs ++;
//		if (tracelevel!=0)	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(q[k][l] >= MAF)) fprintf(fptr,"\nk=%d l=%d numSNPs=%d q=%f\n", k, l, numSNPs, q[k][l]);

		if (q[k][l] == 0.0) accum (&q_0[0], 1.0);
		for (j=1; j<classes; j++)	if ((q[k][l] > ((double)(j-1)/classes))&&(q[k][l] <= ((double)j/classes))) accum (&q_0[j], 1.0);
	}

	fprintf(ffreq, "%f ", accsum(&q_0[0]));
	accum (&AVE_q_0[0], accsum(&q_0[0]));
	for (j=1; j<classes; j++)
	{
		fprintf(ffreq, "%f ", accsum(&q_0[j]));
		accum (&AVE_q_0[j], accsum(&q_0[j]));
	}

	/* ******************* ADDITIVE AND DOMINANCE VARIANCE AND INBREEDING DEPRESSION RATE ***************** */

	va_del = 0.0;
	vd_del = 0.0;
	id_del = 0.0;
    
	va_let = 0.0;
	vd_let = 0.0;
	id_let = 0.0;
    
	va_od = 0.0;
	vd_od = 0.0;
	id_od = 0.0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s[k][l] != 0.0))
	{
		if (s[k][l]<= (-0.9))
		{
			d_a_let = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
			alfa_a_let = (-s[k][l]/2.0) + ( d_a_let * (2.0*q[k][l] - 1.0) );
			va_let += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_let * alfa_a_let;
			vd_let += (2.0 * d_a_let * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_let * q[k][l] * (1.0-q[k][l]));
			id_let += (2.0 * d_a_let * q[k][l] * (1.0-q[k][l]));
		}
		else if (s[k][l]>0.0)
		{
			s1[k][l]=(s[k][l]*hs[k][l])/(1+s[k][l]*hs[k][l]);
			s2[k][l]=(s[k][l]*(1+hs[k][l]))/(1+s[k][l]*hs[k][l]);
            
			d_a_od = (s[k][l]*(hs[k][l]-0.5))/(1+s[k][l]*hs[k][l]);
			alfa_a_od = (q[k][l] * s2[k][l]) - ((1-q[k][l]) * s1[k][l]);
            
			va_od += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_od * alfa_a_od;
			vd_od += (2.0 * d_a_od * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_od * q[k][l] * (1.0-q[k][l]));
			id_od += (2.0 * d_a_od * q[k][l] * (1.0-q[k][l]));
		}
		else
		{
			d_a_del = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
			alfa_a_del = (-s[k][l]/2.0) + ( d_a_del * (2.0*q[k][l] - 1.0) );
 			va_del += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_del * alfa_a_del;
 			vd_del += (2.0 * d_a_del * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_del * q[k][l] * (1.0-q[k][l]));
			id_del += (2.0 * d_a_del * q[k][l] * (1.0-q[k][l]));
		}
	}

	accum (&AVE_VA_let, va_let);
	accum (&AVE_VD_let, vd_let);
	accum (&AVE_ID_let, id_let);

	accum (&AVE_VA_od, va_od);
	accum (&AVE_VD_od, vd_od);
	accum (&AVE_ID_od, id_od);

	accum (&AVE_VA_del, va_del);
	accum (&AVE_VD_del, vd_del);
	accum (&AVE_ID_del, id_del);

	accum (&AVE_VA, va_let + va_od + va_del);
	accum (&AVE_VD, vd_let + vd_od + vd_del);
	accum (&AVE_ID, id_let + id_od + id_del);
    
    
	/* ******************* del, let, od  ***************** */

	numSNPslet = 0.0;
	numSNPsod = 0.0;
	numSNPsdel = 0.0;
	numSNPsneu = 0.0;
    
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s[k][l] != 0.0))
	{
		if (s[k][l]<= (-0.9))
		{
			numSNPslet++;
			accum (&qlet, q[k][l]);
		}
		else if (s[k][l]>0.0)
		{
			numSNPsod ++;
			accum (&qod, q[k][l]);
		}
		else
		{
			numSNPsdel ++;
 			accum (&qdel, q[k][l]);
		} 
	}
    
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(q[k][l] >= MAF)&&(s[k][l]==0.0))
	{
		numSNPsneu ++;
		accum (&qneu, q[k][l]);
	}
    
	accum (&AVE_numSNPslet, numSNPslet);
	accum (&AVE_numSNPsod, numSNPsod);
	accum (&AVE_numSNPsdel, numSNPsdel);
	accum (&AVE_numSNPsneu, numSNPsneu);
}

/* ***************************************************** */

genotypic_values ()
{
//	if (tracelevel!=0)	fprintf(fptr,"\n\n Genotypic values \n");

	for (i=0; i<NIND; i++)
	{
		pm_s[i] = 1.0;

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	pm_s[i] *= (1.0 + s[k][l]) /* aa */;
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
			else	pm_s[i] *= (1.0 + (s[k][l]*hs[k][l])) /* Aa */;
		}
		//if (tracelevel!=0)
		//{
			//if (i <= 10)
			//fprintf(fptr,"genotypic_values %d    pm_s = %f\n", i, pm_s[i]);
			//for (k=0; k<NCRO; k++) if ((i==0)&&(k<=2)&&(gm[i][k][0]!=0)) fprintf (fptr,"%d   gm0=%d   gm1=%d\n", i, gm[i][k][0], gm[i][k][1]);
		//}
	}
}

/* ***************************************************** */

phenotypeB ()
{
	int ii, it;
	double gsum_s=0.0, gsum2_s=0.0, gsum_f=0.0;

	for (i=0; i<NIND; i++)
	{
		gsum_s += pm_s[i];
		gsum2_s += (pm_s[i]*pm_s[i]);
	}

	accum (&gmean_s, gsum_s/(double)NIND);
	accum (&gvar_s, (gsum2_s - (gsum_s*gsum_s / (double)NIND)) / ((double)NIND - 1.0));

//	if (tracelevel!=0)   fprintf(fptr,"\ngmean_s = %f  gvar_s = %f C2 = %f\n",
//	gsum_s/(double)NIND, (gsum2_s - (gsum_s*gsum_s / (double)NIND)) / (double)NIND, ( (gsum2_s*(double)NIND) / (gsum_s*gsum_s) ) - 1.0);

}

/* ***************************************************** */

dumpphenotypes()
{
//	if (tracelevel==0)   return (0);

//	fprintf(fptr,"\n Fitness values\n");
//	for (i=0; i<NIND; i++)   fprintf(fptr,"i=%d pm_s=%f\n", i, pm_s[i]);
}

/* ***************************************************** */

plink_files()
{
	int last, lastpos;

	ffped = fopen ("data.ped","w");
	ffmap = fopen ("data.map","w");
	ffPphen = fopen ("qt.phe","w");
    
	ffped = fopen ("data.ped","a");
	ffmap = fopen ("data.map","a");
	ffPphen = fopen ("qt.phe","a");

	// data.map

	last = 1;
	lastpos = 1;
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(q[k][l] >= MAF)&&(s[k][l]==0.0))
	{
		if (chrom[k][l] != last)
		{
			last = chrom[k][l];
			lastpos = pos[k][l];
		}
		fprintf(ffmap,"%d SNP%llu 0 %llu\n", chrom[k][l], pos[k][l], pos[k][l]-lastpos);
	}

	for (i=0; i<NIND; i++)
	{
		// data.ped

		fprintf(ffped,"1 IND%d 0 0 1 -9 ", i);
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(q[k][l] >= MAF)&&(s[k][l]==0.0))
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		fprintf(ffped,"T T ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffped,"A A ");
	    		else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffped,"T A ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	fprintf(ffped,"A T ");
		}
		fprintf(ffped,"\n");
	}

	for (i=0; i<NIND; i++)
	{
		// qt.phe
        
		fprintf(ffPphen,"%f\n", pm_s[i]);
	}
    
	fclose (ffped);
	fclose (ffmap);
	fclose (ffPphen);
}

/* ***************************************************** */
//int status1 = system("bash shell_F_values");
/* ***************************************************** */

estimates_of_F ()
{
	// FREQUENCIES FROM CURRENT GENERATION

	numLOCI = 0.0;
	EHom = 0.0;
	EHet = 0.0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(q[k][l] >= MAF)&&(s[k][l]==0.0))
	{
		numLOCI ++;
		EHom += (1.0 - 2.0*q[k][l]*(1.0-q[k][l]));
		EHet += 2.0*q[k][l]*(1.0-q[k][l]);
	}

	for (i=0; i<NIND; i++)
	{
		// FREQUENCIES FROM CURRENT GENERATION

		Hom[i] = 0.0;
		H_AA[i] = 0.0;
		H_Aa[i] = 0.0;
		H_aa[i] = 0.0;
		sum_VR1[i] = 0.0;
		sum_VR2[i] = 0.0; 
		sum_YANG1[i] = 0.0; 
		sum_YANG2[i] = 0.0; 
		sum_LH2[i] = 0.0;
		sum_yang2_delta[i] = 0.0;

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(q[k][l] >= MAF)&&(s[k][l]==0.0))
		{
	    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
			{
				Hom[i] ++;
				H_aa[i] ++;
				sum_VR1[i] += (2.0-(2.0*q[k][l]))*(2.0-(2.0*q[k][l]));
				sum_VR2[i] += ( ( (2.0-2.0*q[k][l])*(2.0-2.0*q[k][l]) ) / (2.0*q[k][l]*(1.0-q[k][l])) ) - 1.0; 
				sum_YANG1[i] += ( (4.0)-((1.0+2.0*q[k][l])*2.0)+(2.0*q[k][l]*q[k][l]) ); 
				sum_YANG2[i] += ( (4.0)-((1.0+2.0*q[k][l])*2.0)+(2.0*q[k][l]*q[k][l]) ) / (2.0*q[k][l]*(1.0-q[k][l]));
				sum_yang2_delta[i] += (1.0/q[k][l]);
			}
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))
			{
				Hom[i] ++;
				H_AA[i] ++;
				sum_VR1[i] += (0.0-(2.0*q[k][l]))*(0.0-(2.0*q[k][l]));
				sum_VR2[i] += ( ( (0.0-(2.0*q[k][l]))*(0.0-(2.0*q[k][l])) ) / (2.0*q[k][l]*(1.0-q[k][l])) ) - 1.0; 
				sum_YANG1[i] += (2.0*q[k][l]*q[k][l]); 
				sum_YANG2[i] += (2.0*q[k][l]*q[k][l]) / (2.0*q[k][l]*(1.0-q[k][l])); 
				sum_yang2_delta[i] += (1.0/(1.0-q[k][l]));
			}
			else
			{
				H_Aa[i] ++;
				sum_VR1[i] += (1.0-(2.0*q[k][l]))*(1.0-(2.0*q[k][l]));
				sum_VR2[i] += ( ( (1.0-(2.0*q[k][l]))*(1.0-(2.0*q[k][l])) ) / (2.0*q[k][l]*(1.0-q[k][l])) ) - 1.0; 
				sum_YANG1[i] += ( (1.0)-((1.0+2.0*q[k][l])*1.0)+(2.0*q[k][l]*q[k][l]) ); 
				sum_YANG2[i] += ( (1.0)-((1.0+2.0*q[k][l])*1.0)+(2.0*q[k][l]*q[k][l]) ) / (2.0*q[k][l]*(1.0-q[k][l])); 
				sum_LH2[i] += ( 1.0 / (2.0*q[k][l]*(1.0-q[k][l])) ); 
			}
		}

		fvr1[i] = (sum_VR1[i]/EHet)-1.0;
		fvr2[i] = sum_VR2[i] / numLOCI;
		fyang1[i] = sum_YANG1[i] / EHet;
		fyang2[i] = sum_YANG2[i] / numLOCI;
		fLH1[i] = (Hom[i]-EHom) / (numLOCI-EHom);
		fLH2[i] = 1.0 - ( sum_LH2[i] / numLOCI);
		fhom[i] = Hom[i] / numLOCI;
		fyang2_delta[i] = (sum_yang2_delta[i] / numLOCI) - 1.0;
	}
}

/* ***************************************************** */

regression_pm_F ()
{
	int n;

	fFfile = fopen ("data.F","r");

	for (i=0; i<NIND; i++)
	{
		fscanf(fFfile,"%lf", &w);
		genvalue[i] = log(w);

		fscanf(fFfile,"%lf", &w);
		FhatI[i] = w;

		fscanf(fFfile,"%lf", &w);
		FhatII[i] = w;

		fscanf(fFfile,"%lf", &w);
		FhatIII[i] = w;

		fscanf(fFfile,"%lf", &w);
		Fhom[i] = w;

		fscanf(fFfile,"%lf", &w);
		Fexh[i] = w;

		fscanf(fFfile,"%lf", &w);
		Froh100[i] = w;

		fscanf(fFfile,"%lf", &w);
		Froh1000[i] = w;

		fscanf(fFfile,"%lf", &w);
		Froh5000[i] = w;
        
		//if (tracelevel!=0) fprintf(fptr,"%d    genvalue[i] = %f\n", i, genvalue[i]);
	}

	close(fFfile);

	/* ********************* New estimates ********************************* */

	for (i=0; i<NIND; i++)
	{
 // 		fprintf(fvaluesF, "phe= %f     FhatI[i]= %f  fvr2[i]= %f     FhatII[i]= %f  fLH2[i]= %f      FhatIII[i]= %f  fyang2[i]= %f  fyang2_d[i]= %f     Fexh[i]= %f  fLH1[i]= %f\n",
//		genvalue[i], FhatI[i], fvr2[i], FhatII[i], fLH2[i], FhatIII[i], fyang2[i], fyang2_delta[i], Fexh[i], fLH1[i]); 

 		fprintf(fvaluesF, "phe= %f     fyang2[i]= %f  fyang2_d[i]= %f  H_AA[i]= %f  H_Aa[i]= %f  H_aa[i]= %f\n",
		genvalue[i], fyang2[i], fyang2_delta[i], H_AA[i], H_Aa[i], H_aa[i]); 

		FhatI[i] = fvr2[i];
		FhatII[i] = fLH2[i];
		FhatIII[i] = fyang2[i];
		Fexh[i] = fLH1[i];
	}

	/* ******************************************************************* */

	for (i=0; i<NIND; i++)
	{
		if (genvalue[i] > -99999.0)
		{
			accum (&Phe, genvalue[i]);
			accum (&AFhatI, FhatI[i]);
			accum (&AFhatII, FhatII[i]);
			accum (&AFhatIII, FhatIII[i]);
			accum (&AFhom, Fhom[i]);
			accum (&AFexh, Fexh[i]);
			accum (&AFroh100, Froh100[i]);
			accum (&AFroh1000, Froh1000[i]);
			accum (&AFroh5000, Froh5000[i]);
   			accum (&AFvr1, fvr1[i]);
   			accum (&AFyang1, fyang1[i]);
 		}
		if (genvalue[i] > -99999.0)
		{
   			covaccum (&PheFhatI, genvalue[i], FhatI[i]);
   			covaccum (&PheFhatII, genvalue[i], FhatII[i]);
   			covaccum (&PheFhatIII, genvalue[i], FhatIII[i]);
   			covaccum (&PheFhom, genvalue[i], Fhom[i]);
   			covaccum (&PheFexh, genvalue[i], Fexh[i]);
   			covaccum (&PheFroh100, genvalue[i], Froh100[i]);
   			covaccum (&PheFroh1000, genvalue[i], Froh1000[i]);
   			covaccum (&PheFroh5000, genvalue[i], Froh5000[i]);
   			covaccum (&PheFvr1, genvalue[i], fvr1[i]);
   			covaccum (&PheFyang1, genvalue[i], fyang1[i]);
		}

		covaccum (&FhatIFhatII, FhatI[i], FhatII[i]);
		covaccum (&FhatIFhatIII, FhatI[i], FhatIII[i]);
		covaccum (&FhatIFhom, FhatI[i], Fhom[i]);
		covaccum (&FhatIFexh, FhatI[i], Fexh[i]);
		covaccum (&FhatIFroh100, FhatI[i], Froh100[i]);
		covaccum (&FhatIFroh1000, FhatI[i], Froh1000[i]);
		covaccum (&FhatIFroh5000, FhatI[i], Froh5000[i]);

		covaccum (&FhatIIFhatIII, FhatII[i], FhatIII[i]);
		covaccum (&FhatIIFhom, FhatII[i], Fhom[i]);
		covaccum (&FhatIIFexh, FhatII[i], Fexh[i]);
		covaccum (&FhatIIFroh100, FhatII[i], Froh100[i]);
		covaccum (&FhatIIFroh1000, FhatII[i], Froh1000[i]);
		covaccum (&FhatIIFroh5000, FhatII[i], Froh5000[i]);

		covaccum (&FhatIIIFhom, FhatIII[i], Fhom[i]);
		covaccum (&FhatIIIFexh, FhatIII[i], Fexh[i]);
		covaccum (&FhatIIIFroh100, FhatIII[i], Froh100[i]);
		covaccum (&FhatIIIFroh1000, FhatIII[i], Froh1000[i]);
		covaccum (&FhatIIIFroh5000, FhatIII[i], Froh5000[i]);

		covaccum (&FhomFexh, Fhom[i], Fexh[i]);
		covaccum (&FhomFroh100, Fhom[i], Froh100[i]);
		covaccum (&FhomFroh1000, Fhom[i], Froh1000[i]);
		covaccum (&FhomFroh5000, Fhom[i], Froh5000[i]);

		covaccum (&FexhFroh100, Fexh[i], Froh100[i]);
		covaccum (&FexhFroh1000, Fexh[i], Froh1000[i]);
		covaccum (&FexhFroh5000, Fexh[i], Froh5000[i]);

		covaccum (&Froh100Froh1000, Froh100[i], Froh1000[i]);
		covaccum (&Froh100Froh5000, Froh100[i], Froh5000[i]);

		covaccum (&Froh1000Froh5000, Froh1000[i], Froh5000[i]);
        
        //if (tracelevel!=0) fprintf(fptr,"%d    genvalue[i] = %f  FhatI[i]=%f   PheFhatI=%f   AFhatI=%f\n", i, genvalue[i], FhatI[i], covariance(&PheFhatI), variance(&AFhatI));
	}

	IDFhatI = covariance(&PheFhatI) / variance(&AFhatI);
	IDFhatII = covariance(&PheFhatII) / variance(&AFhatII);
	IDFhatIII = covariance(&PheFhatIII) / variance(&AFhatIII);
	IDFhom = covariance(&PheFhom) / variance(&AFhom);
	IDFexh = covariance(&PheFexh) / variance(&AFexh);
	IDFroh100 = covariance(&PheFroh100) / variance(&AFroh100);
	IDFroh1000 = covariance(&PheFroh1000) / variance(&AFroh1000);
	if (accmean(&AFroh5000) > 0.0)  IDFroh5000 = covariance(&PheFroh5000) / variance(&AFroh5000);

	IDFvr1 = covariance(&PheFvr1) / variance(&AFvr1);
	IDFyang1 = covariance(&PheFyang1) / variance(&AFyang1);
    
	if (tracelevel!=0) fprintf(fptr,"rep=%d    ID(2dpq)=%f    IDFhatI=%f    IDFhatII=%f    IDFhatIII=%f    IDFhom=%f    IDFexh=%f    IDFroh100=%f    IDFroh1000=%f    IDFroh5000=%f\n", rep, (id_let + id_od + id_del), -IDFhatI, -IDFhatII, -IDFhatIII, -IDFhom, -IDFexh, -IDFroh100, -IDFroh1000, -IDFroh5000);
    
    //if (tracelevel!=0) fprintf(fptr,"rep=%d    AFhatI=%f    AFhatII=%f    AFhatIII=%f    AFhom=%f    AFexh=%f    AFroh100=%f    AFroh1000=%f    AFroh5000=%f\n", rep, accmean(&AFhatI), accmean(&AFhatII), accmean(&AFhatIII),accmean(&AFhom), accmean(&AFexh), accmean(&AFroh100), accmean(&AFroh1000), accmean(&AFroh5000));
    
	accum (&AVE_IDFhatI, -IDFhatI);
	accum (&AVE_IDFhatII, -IDFhatII);
	accum (&AVE_IDFhatIII, -IDFhatIII);
	accum (&AVE_IDFhom, -IDFhom);
	accum (&AVE_IDFexh, -IDFexh);
	accum (&AVE_IDFroh100, -IDFroh100);
	accum (&AVE_IDFroh1000, -IDFroh1000);
	accum (&AVE_IDFroh5000, -IDFroh5000);
    	accum (&AVE_IDFvr1, -IDFvr1);
    	accum (&AVE_IDFyang1, -IDFyang1);

	SqE_FhatI = pow((id_let + id_od + id_del) + IDFhatI, 2.0);
	SqE_FhatII = pow((id_let + id_od + id_del) + IDFhatII, 2.0);
	SqE_FhatIII = pow((id_let + id_od + id_del) + IDFhatIII, 2.0);
	SqE_Fhom = pow((id_let + id_od + id_del) + IDFhom, 2.0);
	SqE_Fexh = pow((id_let + id_od + id_del) + IDFexh, 2.0);
	SqE_Froh100 = pow((id_let + id_od + id_del) + IDFroh100, 2.0);
	SqE_Froh1000 = pow((id_let + id_od + id_del) + IDFroh1000, 2.0);
	SqE_Froh5000 = pow((id_let + id_od + id_del) + IDFroh5000, 2.0);
	SqE_Fvr1 = pow((id_let + id_od + id_del) + IDFvr1, 2.0);
	SqE_Fyang1 = pow((id_let + id_od + id_del) + IDFyang1, 2.0);
    
	if (tracelevel!=0) fprintf(fptr,"rep=%d    SqE_FhatI =%f    SqE_FhatII =%f    SqE_FhatIII =%f    SqE_Fhom =%f    SqE_Fexh =%f    SqE_Froh100 =%f    SqE_Froh1000 =%f    SqE_Froh5000 =%f\n", rep, SqE_FhatI, SqE_FhatII, SqE_FhatIII, SqE_Fhom, SqE_Fexh, SqE_Froh100, SqE_Froh1000, SqE_Froh5000);
    
	accum (&AVE_SqE_FhatI, SqE_FhatI);
	accum (&AVE_SqE_FhatII, SqE_FhatII);
	accum (&AVE_SqE_FhatIII, SqE_FhatIII);
	accum (&AVE_SqE_Fhom, SqE_Fhom);
	accum (&AVE_SqE_Fexh, SqE_Fexh);
	accum (&AVE_SqE_Froh100, SqE_Froh100);
	accum (&AVE_SqE_Froh1000, SqE_Froh1000);
	accum (&AVE_SqE_Froh5000, SqE_Froh5000);
	accum (&AVE_SqE_Fvr1, SqE_Fvr1);
	accum (&AVE_SqE_Fyang1, SqE_Fyang1);
}

/* ***************************************************** */

ID_HOMOZYGOTES()
{
	int rnd;

	for (i=0; i<NIND; i++)
	{
		FITN[i] = 1.0;
		FITN_F[i] = 1.0;

		if (uniform() < 0.5) rnd = 1;
		else				 rnd = 0;

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		{
	    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))
			{
				FITN[i] *= (1.0 + s[k][l]*hs[k][l]);

				if (rnd == 1)
				{
					FITN_F[i] *= (1.0 + s[k][l]);
				}
				else		/*11*/;
			}
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
			{
				FITN[i] *= (1.0 + s[k][l]*hs[k][l]);

				if (rnd == 1)	/*11*/;
				else
				{
					FITN_F[i] *= (1.0 + s[k][l]);
				}
			}
			else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
			{
				FITN[i] *= (1.0 + s[k][l]);
				FITN_F[i] *= (1.0 + s[k][l]);
			}
		}
	}

	mean_FITN = 0.0;
	mean_FITN_F = 0.0;

	for (i=0; i<NIND; i++)
	{
		mean_FITN += FITN[i]/(double)NIND;
		mean_FITN_F += FITN_F[i]/(double)NIND;
	}

	accum(&AVE_mean_FITN, mean_FITN);
	accum(&AVE_mean_FITN_F, mean_FITN_F);

	if (mean_FITN_F > 0.0) accum(&AVE_ID_FITN, -log(mean_FITN_F/mean_FITN));

	return(0);
}

/* ***************************************************** */

settozero()
{
	for (j=0; j<classes; j++)	initacc (&q_0[j]);
}

/* ***************************************************** */

printout()
{
        
	fgen = fopen ("genfile.dat","w");

	fprintf(fgen, "NINDNP=%d  NIND=%d  NCRO=%d  TOTLOCI=%d  numSNPNP=%d  numSNPsneu=%f  numSNPsdel=%f  numSNPslet=%f  numSNPod=%f  LEQ_NP=%f\n", NINDNP, NIND, NCRO, TOTLOCI, numSNPNP, accmean(&AVE_numSNPsneu), accmean(&AVE_numSNPsdel), accmean(&AVE_numSNPslet), accmean(&AVE_numSNPsod), LEQ_NP);

	fprintf(fgen, "\n*********** Mean fitness and B (2dpq) ***********\n\n");

	fprintf(fgen, "W          B          Bdel        Blet        Bod\n");
		fprintf(fgen, "%10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n", accmean(&gmean_s), accmean(&AVE_ID), accmean(&AVE_ID_del), accmean(&AVE_ID_let), accmean(&AVE_ID_od));

	fprintf(fgen, "\n*********** ID from homozygotes ***********\n\n");

	fprintf(fgen, "mean_FITN=%6.4f    mean_FITN_F=%6.4f    ID_F1=%6.4f\n", accmean(&AVE_mean_FITN), accmean(&AVE_mean_FITN_F), accmean(&AVE_ID_FITN));
           
	fprintf(fgen, "\n*********** F summary ***********\n\n");
    
	fprintf(fgen, "AFhatI=%6.4f    AFhatII=%6.4f    AFhatIII=%6.4f\n", accmean(&AFhatI), accmean(&AFhatII), accmean(&AFhatIII)); 
	fprintf(fgen, "AFhom=%6.4f    AFexh=%6.4f\n", accmean(&AFhom), accmean(&AFexh)); 
	fprintf(fgen, "AFvr1=%6.4f    AFyang1=%6.4f\n", accmean(&AFvr1), accmean(&AFyang1)); 
	fprintf(fgen, "AFroh100=%6.4f    AFroh1000=%6.4f    AFroh5000=%6.4f\n", accmean(&AFroh100), accmean(&AFroh1000), accmean(&AFroh5000)); 
    
	fprintf(fgen, "\n*********** correlations ***********\n\n");
    
    	fprintf(fgen, "             FhatI FhatII    FhatIII Fhom     Fexh      Froh100  Froh1000 Froh5000\n");
	fprintf(fgen, "Phe        %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f   %6.4f\n",
		correlation(&PheFhatI), correlation(&PheFhatII), correlation(&PheFhatIII), correlation(&PheFhom), correlation(&PheFexh), correlation(&PheFroh100), correlation(&PheFroh1000), correlation(&PheFroh5000)); 
	fprintf(fgen, "FhatI              %6.4f   %6.4f  %6.4f  %6.4f   %6.4f   %6.4f   %6.4f\n",
		correlation(&FhatIFhatII), correlation(&FhatIFhatIII), correlation(&FhatIFhom), correlation(&FhatIFexh), correlation(&FhatIFroh100), correlation(&FhatIFroh1000), correlation(&FhatIFroh5000)); 
	fprintf(fgen, "FhatII                       %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n",
		correlation(&FhatIIFhatIII), correlation(&FhatIIFhom), correlation(&FhatIIFexh), correlation(&FhatIIFroh100), correlation(&FhatIIFroh1000), correlation(&FhatIIFroh5000)); 
	fprintf(fgen, "FhatIII                               %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n",
		correlation(&FhatIIIFhom), correlation(&FhatIIIFexh), correlation(&FhatIIIFroh100), correlation(&FhatIIIFroh1000), correlation(&FhatIIIFroh5000)); 
	fprintf(fgen, "Fhom                                           %6.4f   %6.4f   %6.4f   %6.4f\n",
		correlation(&FhomFexh), correlation(&FhomFroh100), correlation(&FhomFroh1000), correlation(&FhomFroh5000)); 
	fprintf(fgen, "Fexh                                                    %6.4f   %6.4f   %6.4f\n",
		correlation(&FexhFroh100), correlation(&FexhFroh1000), correlation(&FexhFroh5000)); 
	fprintf(fgen, "Froh100                                                          %6.4f   %6.4f\n", correlation(&Froh100Froh1000), correlation(&Froh100Froh5000)); 
	fprintf(fgen, "Froh1000                                                                  %6.4f\n", correlation(&Froh1000Froh5000));
    
    
	fprintf(fgen, "\n*********** B summary ***********\n\n");

	fprintf(fgen, "ID_2dpq=%6.4f\n", accmean(&AVE_ID)); 
	fprintf(fgen, "ID_F1=%6.4f\n", accmean(&AVE_ID_FITN));
	fprintf(fgen, "IDFhatI=%6.4f    IDFhatII=%6.4f    IDFhatIII=%6.4f\n", accmean(&AVE_IDFhatI), accmean(&AVE_IDFhatII), accmean(&AVE_IDFhatIII)); 
	fprintf(fgen, "IDFhom=%6.4f    IDFexh=%6.4f\n", accmean(&AVE_IDFhom), accmean(&AVE_IDFexh)); 
	fprintf(fgen, "IDFvr1=%6.4f    IDFyang1=%6.4f\n", accmean(&AVE_IDFvr1), accmean(&AVE_IDFyang1)); 
	fprintf(fgen, "IDFroh100=%6.4f    IDFroh1000=%6.4f    IDFroh5000=%6.4f\n", accmean(&AVE_IDFroh100), accmean(&AVE_IDFroh1000), accmean(&AVE_IDFroh5000)); 
    
	fprintf(fgen, "\n*********** SqE summary ***********\n\n");
    
	fprintf(fgen, "%10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n", accmean(&AVE_SqE_Fvr1), accmean(&AVE_SqE_FhatI), accmean(&AVE_SqE_FhatII), accmean(&AVE_SqE_Fyang1), accmean(&AVE_SqE_FhatIII), accmean(&AVE_SqE_Fhom), accmean(&AVE_SqE_Fexh), accmean(&AVE_SqE_Froh100), accmean(&AVE_SqE_Froh1000), accmean(&AVE_SqE_Froh5000)); 

	fclose(fgen);
    
	//	Summary ID file
    
	fsumID = fopen ("summaryID","w");
    
	fprintf(fsumID, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", accmean(&AVE_ID), accmean(&AVE_ID_FITN), accmean(&AVE_IDFvr1), accmean(&AVE_IDFhatI), accmean(&AVE_IDFhatII), accmean(&AVE_IDFyang1), accmean(&AVE_IDFhatIII), accmean(&AVE_IDFhom), accmean(&AVE_IDFexh), accmean(&AVE_IDFroh100), accmean(&AVE_IDFroh1000), accmean(&AVE_IDFroh5000)); 
    
	fclose(fsumID);
    
	// Summary Mean file
    
	fsumM = fopen ("summaryM","w");
    
	fprintf(fsumM, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", accmean(&AFvr1), accmean(&AFhatI), accmean(&AFhatII), accmean(&AFyang1), accmean(&AFhatIII), accmean(&AFhom), accmean(&AFexh), accmean(&AFroh100), accmean(&AFroh1000), accmean(&AFroh5000)); 
    
	fclose(fsumM);
    
	// Summary Var file
        
	fsumV = fopen ("summaryV","w");
    
	fprintf(fsumV, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", variance(&AFvr1), variance(&AFhatI), variance(&AFhatII), variance(&AFyang1), variance(&AFhatIII), variance(&AFhom), variance(&AFexh), variance(&AFroh100), variance(&AFroh1000), variance(&AFroh5000)); 
    
	fclose(fsumV);
    
	//	Summary SqE file
    
	fsumSqE = fopen ("summarySqE","w");

	fprintf(fsumSqE, "%6.4f  %6.4f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n", accmean(&AVE_ID), accmean(&AVE_ID_FITN), accmean(&AVE_SqE_Fvr1), accmean(&AVE_SqE_FhatI), accmean(&AVE_SqE_FhatII), accmean(&AVE_SqE_Fyang1), accmean(&AVE_SqE_FhatIII), accmean(&AVE_SqE_Fhom), accmean(&AVE_SqE_Fexh), accmean(&AVE_SqE_Froh100), accmean(&AVE_SqE_Froh1000), accmean(&AVE_SqE_Froh5000)); 

	fclose(fsumSqE);
    
	//	SummaryR (correlations) file
    
	fsumR = fopen ("summaryR","w");
    
	fprintf(fsumR, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
		correlation(&PheFvr1), correlation(&PheFhatI), correlation(&PheFhatII), correlation(&PheFyang1), correlation(&PheFhatIII), correlation(&PheFhom), correlation(&PheFexh), correlation(&PheFroh100), correlation(&PheFroh1000), correlation(&PheFroh5000)); 
	fprintf(fsumR, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
		correlation(&FhatIFhatII), correlation(&FhatIFhatIII), correlation(&FhatIFhom), correlation(&FhatIFexh), correlation(&FhatIFroh100), correlation(&FhatIFroh1000), correlation(&FhatIFroh5000)); 
	fprintf(fsumR, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
		correlation(&FhatIIFhatIII), correlation(&FhatIIFhom), correlation(&FhatIIFexh), correlation(&FhatIIFroh100), correlation(&FhatIIFroh1000), correlation(&FhatIIFroh5000)); 
	fprintf(fsumR, "%6.4f %6.4f %6.4f %6.4f %6.4f ",
		correlation(&FhatIIIFhom), correlation(&FhatIIIFexh), correlation(&FhatIIIFroh100), correlation(&FhatIIIFroh1000), correlation(&FhatIIIFroh5000)); 
	fprintf(fsumR, "%6.4f %6.4f %6.4f %6.4f ",
		correlation(&FhomFexh), correlation(&FhomFroh100), correlation(&FhomFroh1000), correlation(&FhomFroh5000)); 
	fprintf(fsumR, "%6.4f %6.4f %6.4f ",
		correlation(&FexhFroh100), correlation(&FexhFroh1000), correlation(&FexhFroh5000)); 
	fprintf(fsumR, "%6.4f %6.4f ", correlation(&Froh100Froh1000), correlation(&Froh100Froh5000)); 
	fprintf(fsumR, "%6.4f\n", correlation(&Froh1000Froh5000));
    
	fclose(fsumR);
    
	//	Summary B file
    
	fsumB = fopen ("summaryB","w");

	fprintf(fsumB, "%10.8f  %10.8f  %10.8f  %10.8f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", accmean(&AVE_ID), accmean(&AVE_ID_del), accmean(&AVE_ID_let), accmean(&AVE_ID_od), accmean(&AVE_numSNPsneu), accmean(&AVE_numSNPsdel), accmean(&AVE_numSNPslet), accmean(&AVE_numSNPsod), accmean(&qneu), accmean(&qdel), accmean(&qlet), accmean(&qod)); 

	fclose(fsumB); 
}

/* ***************************************************** */

lookfortext(s)
char *s;
{
   int len, i, curchar;
   char c;

   curchar = 0;
   len = 0;

   for (i=0; i<=100; i++)
   {
      if (s[i] == '\0') break;
      len++;
   }
   do
   {
      c = getc(fpop);

      if (c==s[curchar])
      {
         curchar++;
         if (curchar==len) return(0);
      }
      else curchar = 0;
   }
   while (c != EOF);
}

/* ********************************************************************************************* */
