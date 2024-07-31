// SNP_BP_SLIM3_2.c

#include "libhdr"

#define NN 100000 // Maximum 100000 SNPs segregating
#define CC 20001   // Maximum N=10000 individuals
#define MC 25000 // Maximum 25000 SNPs per chromosome

int i, j, s, ss, m, x, b, nind, HOM1, HOM2, crossover, mut_n;
unsigned long long int pos[NN], xx, a, MUT_pos;
int mut[NN], crom[CC][MC], num_mut_crom[CC];
double w, ps[NN], ef[NN], h[NN], q[NN], FITN[CC];
double FITNparents[CC], MUT_h, MUT_s;
int numSNPchip, SNPchip[NN], numSNP;
char ch;

FILE *fdat, *fgen, *fPpos, *fPphen, *fallsnps, *fLISTQTL, *fcoll;

main()
{
	getseed();
	getintandskip("NIND :",&nind,10,10000);

	readfiles();
	PLINK_files();
	writeseed();

	return(0);
}

/* **************************************************************************** */

readfiles()
{
	fgen = fopen ("dataBP.ped","w");
	fPpos = fopen ("dataBP.map","w");
	fPphen = fopen ("qtBP.phe","w");

	fallsnps = fopen ("list_allsnps","w");
	fLISTQTL = fopen ("list_qtls","w");

	// ********** Read slimout to get SNP positions, effects and frequencies ********** 

	fdat = fopen ("slimout","r");

	// ********** get the position, effects and frequencies of all mutations simulated (numSNP) **********

	lookfortext("Mutations:");

	while (!feof(fdat))
	{
		s ++;
		fscanf(fdat,"%d", &x);
		mut[s] = x;
		fscanf(fdat,"%d", &x);
		for (j=1; j<=4; j++)	fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%llu", &xx);
		pos[s] = xx;
		if (pos[s] == pos[s-1])
		{
			pos[s] = pos[s-1] + 1;
		}
		fscanf(fdat,"%lf", &w);
		ps[s] = w;
		if (ps[s] < -1.0)   ps[s] = (-1.0);
		ef[s] = fabs(w);
		fscanf(fdat,"%lf", &w);
		h[s] = w;
		if (ps[s] == 0.0)	h[s] = 0.0;
		for (j=1; j<=4; j++)	fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%d", &x);
		fscanf(fdat,"%d", &x);
		q[s] = x / (2.0*nind);
		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%c", &ch);
		if (ch != 'I')	ungetc(ch, fdat);
		else		break;
	}
	numSNP = s;

	fclose(fdat);

	// ********** reorder by genome position **********

	for (s=1; s<numSNP; s++)
	for (ss=s; ss<=numSNP; ss++)
	{
		if (pos[ss] < pos[s])
		{
			a=pos[s]; pos[s]=pos[ss]; pos[ss]=a;
			b=mut[s]; mut[s]=mut[ss]; mut[ss]=b;
			w=ps[s]; ps[s]=ps[ss]; ps[ss]=w;
			w=ef[s]; ef[s]=ef[ss]; ef[ss]=w;
			w=h[s]; h[s]=h[ss]; h[ss]=w;
			w=q[s]; q[s]=q[ss]; q[ss]=w;
		}
	}

	// ********** read collect_par_mutations.txt and correct h value **********

	fcoll = fopen ("collect_par_mutations.txt","r");

	while (!feof(fcoll))
	{
		mut_n ++;

		fscanf(fcoll,"%lf", &w);
		MUT_h = w;
		fscanf(fcoll,"%lf", &w);
		MUT_s = w;
		fscanf(fcoll,"%llu", &xx);
		MUT_pos = xx;
	
//		if (mut_n < 10)	printf("\n MUT_h = %f    MUT_s = %f     MUT_pos = %llu", MUT_h, MUT_s, MUT_pos);

		for (s=1; s<=numSNP; s++)
       	if (ps[s] < 0.0)
		{
			if (pos[s] == MUT_pos)
			{
				h[s] = MUT_h;
//				printf("\n pos_s = %llu    MUT_pos = %llu    h_s = %f", pos[s], MUT_pos, h[s]);
//				goto nextmut;
			}
		}
//		nextmut: /**/;
	}

	fclose(fcoll);

	// ********** get genotypes **********

	fdat = fopen ("slimout","r");
	lookfortext("Genomes:");

	fscanf(fdat,"%c", &ch);

	for (i=1; i<=(2*nind);i++)
	{
		m = 0;

		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%d", &x);
		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%c", &ch);
		if (ch == '\n')	goto next;
		else	ungetc(ch, fdat);
		while (!feof(fdat))
		{
			fscanf(fdat,"%c", &ch);

			if (ch == '\n')	break;
			else
			{
				ungetc(ch, fdat);
				m ++;
				fscanf(fdat,"%d", &x);
				crom[i][m] = x;
			}			
		}
		num_mut_crom[i] = m;

		next: /***/;
	}

	fclose(fdat);

	return(0);
}

/* **************************************************************************** */

PLINK_files()
{
	double sum, sum2;

	for (i=1; i<(2*nind);i+=2)
	{
		if (i%2 != 0)
		{
			fprintf(fgen,"1 IND%d 0 0 1 0   ", (i+1)/2);
			FITN[(i+1)/2] = 1.0;
		}
		for (s=1; s<=numSNP; s++)
		{
			if (i == 1)
			{
				// PLINK POSFILE
				fprintf(fPpos,"1 SNP%d 0 %llu\n", s, pos[s]);
			}
			HOM1 = 0; HOM2 = 0;

			for (m=1; m<=num_mut_crom[i]; m++)
			{
				if (crom[i][m] == mut[s])
				{
					HOM1 = 1;
					break;
				}
			}
			for (m=1; m<=num_mut_crom[i+1]; m++)
			{
				if (crom[i+1][m] == mut[s])
				{
					HOM2 = 1;
					break;
				}
			}

			if ((HOM1==0) && (HOM2==0))		fprintf(fgen,"1 1  ");
			else if ((HOM1==1) && (HOM2==0))
			{
				fprintf(fgen,"2 1  ");
				FITN[(i+1)/2] *= (1.0 + ps[s]*h[s]);
			}
			else if ((HOM1==0) && (HOM2==1))
			{
				fprintf(fgen,"1 2  ");
				FITN[(i+1)/2] *= (1.0 + ps[s]*h[s]);	
			}
			else
			{
				fprintf(fgen,"2 2  ");
				FITN[(i+1)/2] *= (1.0 + ps[s]);
			}
		}
		fprintf(fgen,"\n");

		if (i%2 != 0)
		{
			// PLINK PHENOFILE
			fprintf(fPphen,"1 IND%d %f\n", (i+1)/2, FITN[(i+1)/2]);
		}
	}

	// List all SNPs
//	fprintf(fallsnps,"chr  pos  s  a  h  q\n");
	fprintf(fallsnps,"%d\n", numSNP);
	for (s=1; s<=numSNP; s++)	fprintf(fallsnps,"1  %llu  %f  %f  %f  %f\n", pos[s], ps[s], ef[s], h[s], q[s]);

	//FILE WITH QTLs
	fprintf(fLISTQTL,"SNP  pos  s  a  h  q\n");
	for (s=1; s<=numSNP; s++)
	if (ps[s] != 0)
	fprintf(fLISTQTL,"%d  %llu  %f  %f  %f  %f\n", s, pos[s], ps[s], ef[s], h[s], q[s]);

	return(0);
}

/* **************************************************************************** */

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
      c = getc(fdat);

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

