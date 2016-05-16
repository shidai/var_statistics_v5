#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
////#include <gsl/gsl_rng.h>
////#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include "T2toolkit.h"
//#include "cpgplot.h"
//#include <omp.h>
#define WINSIZE 2.5
#define PRECISION 0.005

typedef struct acfStruct {
	int n; // number of dynamic spectrum

	double cFlux;  // pulsar flux density
	double whiteLevel;  // white noise level

	double phaseGradient;
	double cFreq; // observing central frequency

	double bw; // observing bandwidth
	double tint; // integration time
	int nchn;
	int nsubint;
	int nsubint_win;

	double f0; // scintillation bandwidth
	double t0; // scintillation time-scale

	int ns; // sampling number of spatial scale 
	int nf; // sampling number of frequency scale
	double size[2]; // sampling boundary
	double steps;
	double stepf;
	double *s; // spatial scale
	double *f; // bw scale
	double *acf2d;  // ACF 
	double *psrt;  // power spactrum
	fftw_complex *eField; // complex electric field
	fftw_complex *intensity;  // intensity 
	double **dynSpec; // dynamic spectrum 
	double **dynSpecWindow; // dynamic spectrum window, nchn*nsubint dimension

	float *dynPlot; // dynamic spectrum for pgplot

	double winsize;
} acfStruct;

typedef struct noiseStruct {
	int n; // number of dynamic spectrum
	int nchn;
	int nsubint;

	float *noisePlot; // dynamic spectrum for pgplot
	double whiteLevel;  // white noise level
	double detection; // detection threshold
} noiseStruct;

typedef struct controlStruct {
	int n; // number of dynamic spectrum
	//char oname[1024]; // output file name
	char dname[1024]; // graphics device

	double cFreq;  // observing central frequency
	double bw;
	double tint;
	double tsub;  // subintegration time
	double chanBW;  // subchannel bandwidth

	int nsub;  // number of subintegrations
	int nsub_win;  // number of subintegrations in the window
	int nchan; // number of subchannels

	double scint_freqbw;   // scintillation bandwidth
	double scint_ts;       // scintillation timescale

	double whiteLevel;   // white noise level, mJy
	double cFlux;        // flux density of pulsars, mJy

	int noplot;

	double winsize;
}controlStruct;

int idft2d (acfStruct *acfStructure);
int dft2d (acfStruct *acfStructure, fftw_complex *out);
int calACF (acfStruct *acfStructure);
int power (acfStruct *acfStructure);
void deallocateMemory (acfStruct *acfStructure);
void allocateMemory (acfStruct *acfStructure);
int simDynSpec (acfStruct *acfStructure, long seed);
int winDynSpec (acfStruct *acfStructure, long seed);
int calculateScintScale (acfStruct *acfStructure, controlStruct *control);
float calculateNDynSpec (acfStruct *acfStructure);
void preAllocateMemory (acfStruct *acfStructure);
//void preAllocateMemory (acfStruct *acfStructure, controlStruct *control);
float find_peak_value (int n, float *s);
int calSize (acfStruct *acfStructure, double *size, double *ratio);
int windowSize (acfStruct *acfStructure, double *size);
int readParams(char *fname, char *dname, int n, controlStruct *control);
//int readParams(char *fname, char *oname, char *dname, int n, controlStruct *control);
void initialiseControl(controlStruct *control);

//void heatMap (acfStruct *acfStructure, char *dname);
//void palett(int TYPE, float CONTRA, float BRIGHT);
//int plotDynSpec (char *pname, char *dname);

//int qualifyVar (acfStruct *acfStructure, noiseStruct *noiseStructure, controlStruct *control);
float chiSquare (float *data, int n, float noise);
float moduIndex (float *data, int n);
float variance (float *data, int n);
int histogram (float *data, int n, float *x, float *val, float low, float up, int step);
float find_max_value (int n, float *s);
float find_min_value (int n, float *s);

float calNoise (noiseStruct *noiseStructure, controlStruct *control);
int simNoise (noiseStruct *noiseStructure, long seed);
int calThreshold (noiseStruct *noiseStructure, float *var_n);

void readDiss (char *Tname, char *Fname, double *tdiss, double *fdiss);
int readDissNum (char *Tname);

float calNoise (noiseStruct *noiseStructure, controlStruct *control)
{
	int nchn, nsubint;
	long seed;
	int i;
	int n = control->n;

	float std_n;
	float *var_n;
	var_n = (float *)malloc(sizeof(float)*n);

	noiseStructure->n = control->n; 

	noiseStructure->nchn = control->nchan; 
	noiseStructure->nsubint = control->nsub_win; 
	nchn = noiseStructure->nchn;
	nsubint = noiseStructure->nsubint;

	noiseStructure->whiteLevel = sqrt(nchn*nsubint)*control->whiteLevel; // mJy

	// allocate memory
	noiseStructure->noisePlot = (float *)malloc(sizeof(float)*nsubint*nchn);

	// simulate noise
	for (i=0; i<noiseStructure->n; i++)
	{
		seed = TKsetSeed();
		simNoise (noiseStructure, seed);

		var_n[i] = variance (noiseStructure->noisePlot, nsubint*nchn);
	}

	std_n = sqrt(variance(var_n, n));
	//calThreshold (noiseStructure, var_n);

	// deallocate memory
	free(noiseStructure->noisePlot);
	free(var_n);

	return std_n;
}

int calThreshold (noiseStruct *noiseStructure, float *var_n)
{
	int i;
	int n = noiseStructure->n;

	float max;
	float percent;
	float threshold;
	int num;

	max = find_max_value (n, var_n);
	//printf ("Results: %f %f %f %f\n", meanVar, varVar, meanVarN, varVar_n);
	
	percent = 1.0;
	threshold = max;
	while (percent >= 0.95)
	{
		threshold -= 0.01;
		num = 0;
		for (i=0; i<n; i++)
		{
			if (var_n[i] <= threshold)
			{
				num++;
			}
		}
		percent = (float)(num)/n;
	}

	noiseStructure->detection = threshold;

	return 0;
}

float calculateNDynSpec (acfStruct *acfStructure)
{
	long seed;
	int i;

	int n = acfStructure->n;
	int nsub = acfStructure->nsubint_win;
	int nchan = acfStructure->nchn;

	float meanVar;
	//float *var;
	//var = (float*)malloc(sizeof(float)*n);

	acfStructure->dynPlot = (float *)malloc(sizeof(float)*nsub*nchan);

	meanVar = 0.0;
	for (i=0; i<acfStructure->n; i++)
	{
		seed = TKsetSeed();
		//acfStructure->phaseGradient = TKgaussDev(&seed);
		//printf ("Phase gradient: %lf\n", acfStructure->phaseGradient);

		winDynSpec (acfStructure, seed);
		//printf ("Make DynSpec %d\n", i);
		
		// calculate the variance
		meanVar += variance (acfStructure->dynPlot, nsub*nchan);

	}
	meanVar = meanVar/n;
	//printf ("%d %f\n", nchan, meanVar/n);
	
	free(acfStructure->dynPlot);
				
	return meanVar;
}

int calculateScintScale (acfStruct *acfStructure, controlStruct *control)
//int calculateScintScale (acfStruct *acfStructure, controlStruct *control, long seed)
{
	//FILE *fin;
	long seed;

	//printf ("Starting simulating dynamic spectrum\n");
	// moved to preAllocateMemory
	acfStructure->n = control->n; 
	acfStructure->whiteLevel = control->whiteLevel; // mJy
	acfStructure->cFreq = control->cFreq; // MHz

	acfStructure->nchn = control->nchan;
	acfStructure->nsubint = control->nsub;
	acfStructure->bw = control->bw; // MHz
	acfStructure->tint = control->tint;  // s
	//acfStructure->bw = control->chanBW*control->nchan; // MHz
	//acfStructure->tint = control->nsub*control->tsub;  // s

	acfStructure->t0 = control->scint_ts; // s
	acfStructure->f0 = control->scint_freqbw;  // MHz

	acfStructure->winsize = control->winsize;

	acfStructure->cFlux = control->cFlux; // mJy

	//printf ("Scintillation bandwidth: %lf (MHz)\n", acfStructure->f0);
	//printf ("Scintillation time-scale: %lf (s)\n", acfStructure->t0);

	// moved to main
	//preAllocateMemory (acfStructure, control);
	preAllocateMemory (acfStructure);
	allocateMemory (acfStructure);

	calACF (acfStructure);
	power (acfStructure);
		
	seed = TKsetSeed();
		
	acfStructure->phaseGradient = 0.0;
	simDynSpec (acfStructure, seed);

	// deallocate memory
	free(acfStructure->s); 
	free(acfStructure->f); 
	free(acfStructure->acf2d);
	free(acfStructure->psrt);

	fftw_free(acfStructure->eField); 
	fftw_free(acfStructure->intensity); 

	// allocate memory
	//acfStructure->dynSpecWindow = (double **)malloc(sizeof(double *)*nchn);
	//for (i = 0; i < nchn; i++)
	//{
	//	acfStructure->dynSpecWindow[i] = (double *)malloc(sizeof(double)*nsubint);
	//}

	//for (i=0; i<acfStructure->n; i++)
	//{
	//	seed = TKsetSeed();
	//	//acfStructure->phaseGradient = TKgaussDev(&seed);
	//	acfStructure->phaseGradient = 0.0;
	//	//printf ("Phase gradient: %lf\n", acfStructure->phaseGradient);

	//	winDynSpec (acfStructure, seed, i);
	//}

	/*
	if (acfStructure->n == 1)
	{
		printf ("Dynamic spectrum is output into %s\n", control->oname);

		if ((fin=fopen(control->oname,"w"))==NULL)
		{
			printf ("Can't open output file!\n");
			exit(1);
		}

		fprintf(fin,"INFO nsub nchn bandwidth tint cFreq\n");
		fprintf(fin,"START %d %d %lf %lf %lf\n",acfStructure->nsubint,acfStructure->nchn,acfStructure->bw,acfStructure->tint,acfStructure->cFreq);

		for (i=0;i<acfStructure->nchn;i++)
		{
			for (j=0;j<acfStructure->nsubint;j++)
			{
				fprintf(fin,"%d %d %lf\n", i, j, acfStructure->dynPlot[0][i*acfStructure->nsubint+j]);
			}
		}

		if (fclose(fin))
		{
			printf ("Can't close output file!\n");
			exit(1);
		}
	}
	else
	{
		printf ("%d dynamic spectra are simulated\n", acfStructure->n);
	}
	*/

	return 0;
}

int calACF (acfStruct *acfStructure)
{
	//printf ("Calculating ACF\n");
	int i,j;
	int ns = acfStructure->ns;
	int nf = acfStructure->nf;
	double *acf;
	acf = (double *)malloc(sizeof(double)*ns*nf);

	// moved to allocateMemory
	//double steps = acfStructure->steps;
	//double stepf = acfStructure->stepf;
	//for (i = 0; i < ns; i++)
	//{
	//	acfStructure->s[i] = -acfStructure->size[1]+i*steps;
	//}

      	//for (i = 0; i < nf; i++)
	//{
	//	acfStructure->f[i] = -acfStructure->size[0]+i*stepf;
	//}

	double rand;
	rand = acfStructure->phaseGradient;
	//printf ("%lf\n",rand);
	
	int n = 0;
	for (i = 0; i < nf; i++)
	{
		for (j = 0; j < ns; j++)
		{
			//acf[n] = exp(-pow((pow(acfStructure->s[j],2.5)+pow(acfStructure->f[i],1.5)),2.0/3.0));
			//acf[n] = exp(-pow((pow(fabs(acfStructure->s[j]),2.5)+pow(fabs(acfStructure->f[i]),1.5)),2.0/3.0));
			acf[n] = exp(-pow((pow(fabs(acfStructure->s[j]+2.0*rand*0.4*acfStructure->f[i]),2.5)+pow(fabs(acfStructure->f[i]),1.5)),2.0/3.0));
			n++;
		}
	}

	//////////////////////////////////////////////////////////////////

	for (i = 0; i < nf; i++)
	{
		for (j = 0; j < ns; j++)
		{
			if (i >= (int)ceil(nf/2.0) && j >= (int)ceil(ns/2.0))
			{
				acfStructure->acf2d[i*ns+j] = acf[ns*(i-(int)ceil(nf/2.0))+(j-(int)ceil(ns/2.0))];
			}
			else if  (i >= (int)ceil(nf/2.0) && j < (int)ceil(ns/2.0))
			{
				acfStructure->acf2d[i*ns+j] = acf[ns*(i-(int)ceil(nf/2.0))+(j+(int)floor(ns/2.0))];
			}
			else if  (i < (int)ceil(nf/2.0) && j < (int)ceil(ns/2.0))
			{
				acfStructure->acf2d[i*ns+j] = acf[ns*(i+(int)floor(nf/2.0))+(j+(int)floor(ns/2.0))];
			}
			else if  (i < (int)ceil(nf/2.0) && j >= (int)ceil(ns/2.0))
			{
				acfStructure->acf2d[i*ns+j] = acf[ns*(i+(int)floor(nf/2.0))+(j-(int)ceil(ns/2.0))];
			}
		}
	}

	free(acf);

	return 0;
}


int dft2d (acfStruct *acfStructure, fftw_complex *out)
{
	int i;
	int n0 = acfStructure->nf;
	int n1 = acfStructure->ns;
	double *in;
	in = (double *)malloc(sizeof(double)*n0*n1);

	fftw_plan p;
	
	p = fftw_plan_dft_r2c_2d (n0, n1, in, out, FFTW_ESTIMATE);
	//p = fftw_plan_dft_r2c_2d (n0, n1, in, out, FFTW_MEASURE);

	for (i = 0; i < n0*n1; i++)
	{
		in[i] = acfStructure->acf2d[i];
	}

	fftw_execute(p);

	fftw_destroy_plan(p);
	free (in);
  
	return 0;
}

int idft2d (acfStruct *acfStructure)
{
	int i;
	int n0 = acfStructure->nf;
	int n1 = acfStructure->ns;

	fftw_complex *in;
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n0*n1);

	fftw_plan p;
	
	p = fftw_plan_dft_2d (n0, n1, in, acfStructure->intensity, FFTW_BACKWARD, FFTW_ESTIMATE);
	//p = fftw_plan_dft_2d (n0, n1, in, acfStructure->intensity, FFTW_BACKWARD, FFTW_MEASURE);

	for (i = 0; i < n0*n1; i++)
	{
		in[i][0] = acfStructure->eField[i][0];
		in[i][1] = acfStructure->eField[i][1];
	}

	fftw_execute(p);

	fftw_destroy_plan(p);
	fftw_free (in);
  
	return 0;
}

int power (acfStruct *acfStructure)
{
	int i;
	int nf = acfStructure->nf;
	int ns = acfStructure->ns;
	/////////////////////////////////////////////////////////////////////////////////

	fftw_complex *out;
	
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*((int)(ns/2)+1));
	
	dft2d (acfStructure, out);

	int n, j;
	n = 0;
	for (i = 0; i < nf; i++)
	{
		for (j = 0; j < ns; j++)
		{
			if (j < (int)(ns/2)+1)
			{
				acfStructure->psrt[n] = sqrt(sqrt(pow(out[i*((int)(ns/2)+1)+j][0],2.0)+pow(out[i*((int)(ns/2)+1)+j][1],2.0)));
			}
			else
			{
				if (i == 0)
				{
					acfStructure->psrt[n] = sqrt(sqrt(pow(out[i*((int)(ns/2)+1)+ns-j][0],2.0)+pow(out[i*((int)(ns/2)+1)+ns-j][1],2.0)));
				}
				else
				{
					acfStructure->psrt[n] = sqrt(sqrt(pow(out[(nf-i)*((int)(ns/2)+1)+ns-j][0],2.0)+pow(out[(nf-i)*((int)(ns/2)+1)+ns-j][1],2.0)));
				}
			}
			//printf ("%lf ", acfStructure->psrt[n]);
			n++;
		}
		//printf ("\n");
	}

	fftw_free(out); 

	return 0;
}

void allocateMemory (acfStruct *acfStructure)
{
	int i;
	int ns, nf;

	double steps = acfStructure->steps;
	double stepf = acfStructure->stepf;

	ns = acfStructure->ns;
	nf = acfStructure->nf;
	
	acfStructure->s = (double *)malloc(sizeof(double)*ns);
	acfStructure->f = (double *)malloc(sizeof(double)*nf);
	acfStructure->acf2d = (double *)malloc(sizeof(double)*ns*nf);
	acfStructure->psrt = (double *)malloc(sizeof(double)*ns*nf);

	acfStructure->eField = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*ns);
	acfStructure->intensity = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*ns);
	
	acfStructure->dynSpec = (double **)malloc(sizeof(double *)*nf);

	for (i = 0; i < nf; i++)
	{
		acfStructure->dynSpec[i] = (double *)malloc(sizeof(double)*ns);
	}

	for (i = 0; i < ns; i++)
	{
		acfStructure->s[i] = -acfStructure->size[1]+i*steps;
	}

      	for (i = 0; i < nf; i++)
	{
		acfStructure->f[i] = -acfStructure->size[0]+i*stepf;
	}

}

void deallocateMemory (acfStruct *acfStructure)
{
	//int n = acfStructure->n; // number of dynamic spectrum
	int nf = acfStructure->nf;
	//int nchn = acfStructure->nchn;
	
	int i;
	for (i = 0; i < nf; i++)
	{
		free(acfStructure->dynSpec[i]);
	}

	//for (i = 0; i < nchn; i++)
	//{
	//	free(acfStructure->dynSpecWindow[i]);
	//}

	//free(acfStructure->dynPlot);
}

int simNoise (noiseStruct *noiseStructure, long seed)
{
	int nchn = noiseStructure->nchn;
	int nsubint = noiseStructure->nsubint;

	int i, j;

	for (i = 0; i < nchn; i++)
	{
		for (j = 0; j < nsubint; j++)
		{
			noiseStructure->noisePlot[i*nsubint+j] = (float)(noiseStructure->whiteLevel*TKgaussDev(&seed));   // create noise image pixels
		}
	}

	return 0;
}

int simDynSpec (acfStruct *acfStructure, long seed)
{
	int nf = acfStructure->nf;
	int ns = acfStructure->ns;

	int i;
	int j;

	int n = 0;
	double sum = 0.0;

	//printf ("Simulating dynamic spectrum\n");
	//seed = TKsetSeed();
	//printf ("seed %ld\n",seed);

	for (i = 0; i < nf*ns; i++)
	{
		//acfStructure->eField[i][0] = acfStructure->psrt[i];
		//acfStructure->eField[i][1] = acfStructure->psrt[i];
		acfStructure->eField[i][0] = acfStructure->psrt[i]*TKgaussDev(&seed);
		acfStructure->eField[i][1] = acfStructure->psrt[i]*TKgaussDev(&seed);
		//printf ("%lf\n",TKgaussDev(&seed));
		//printf ("######### \n");
	}

	/////////////////////////////////////////////////////////////////////////////////

	// ifft
	idft2d (acfStructure);

	// form the matrix and normalize
	for (i = 0; i < nf; i++)
	{
		for (j = 0; j < ns; j++)
		{
			acfStructure->dynSpec[i][j] = pow(acfStructure->intensity[n][1]/(nf*ns),2.0)+pow(acfStructure->intensity[n][0]/(nf*ns),2.0);
			//fprintf (fp, "%lf  ", acfStructure->dynSpec[i][j]);
			sum += pow(acfStructure->intensity[n][1]/(nf*ns),2.0)+pow(acfStructure->intensity[n][0]/(nf*ns),2.0);
			n++;
		}
		//fprintf (fp, "\n");
	}

	sum = sum/n;

	//printf ("Normalization %.10lf\n",sum);
	for (i = 0; i < nf; i++)
	{
		for (j = 0; j < ns; j++)
		{
			acfStructure->dynSpec[i][j] = acfStructure->dynSpec[i][j]/sum;
			//printf ("%d %d\n", i, j);
		}
	}

	return 0;
}

int winDynSpec (acfStruct *acfStructure, long seed)
{
	int nf = acfStructure->nf;
	int ns = acfStructure->ns;

	int nchn = acfStructure->nchn;

	int nsubint = acfStructure->nsubint;
	int nsubint_win = acfStructure->nsubint_win;

	double bw = acfStructure->bw;
	double tint = acfStructure->tint;
	double f0 = acfStructure->f0;
	double t0 = acfStructure->t0;

	int i, ii;
	int j;
	float temp;
	int tempt = (int)(nsubint/nsubint_win);  // number of pixels to average

	double rand1, rand2;
	double rand11, rand22;
	int nf0;
	int ns0;

	double dynSpecWindow;

	//printf ("Simulating dynamic spectrum\n");
	//seed = TKsetSeed();
	//printf ("seed %ld\n",seed);

	rand1 = TKgaussDev(&seed);
	rand2 = TKgaussDev(&seed);
	rand11 = rand1 - floor(rand1);
	rand22 = rand2 - floor(rand2);
	//printf ("rand %lf\n",rand);

	if (f0 == (bw/nchn) && t0 >= (tint/nsubint))
	{
		//printf ("f0 < (bw/nchn) && t0 >= (tint/nsubint)\n");
		nf0 = (int)(rand11*(nf-nchn));
		ns0 = (int)(rand22*(ns-nsubint));

		for (i = 0; i < nchn; i++)
		{
			for (j = 0; j < nsubint_win; j++)
			{
				temp = 0.0;
				for (ii=0; ii<tempt; ii++)
				{
						temp += acfStructure->dynSpec[i+nf0][j*tempt+ns0+ii];
				}
				dynSpecWindow = temp/tempt;
				acfStructure->dynPlot[i*nsubint+j] = (float)(dynSpecWindow*acfStructure->cFlux);   // add in noise here
			}
			//printf ("%d\n", tempf);
		}
	}
	else
	{
		printf ("Number of subintegration too small!\n");
	}

	return 0;
}

int windowSize (acfStruct *acfStructure, double *size)
{
	double bw, f0, tint, t0;
	bw = acfStructure->bw;
	f0 = acfStructure->f0;
	tint = acfStructure->tint;
	t0 = acfStructure->t0;

	if ( (bw/f0) > 6 )
	{
		size[0] = bw/f0;
	}
	else
	{
		size[0] = 6.0;
	}
			  
	if ( (tint/t0) > 6 )
	{
		size[1] = tint/t0;
	}
	else
	{
		size[1] = 6.0;
	}

	double ratio[2];
	calSize (acfStructure, size, ratio);
	//printf ("f0 ratio: %lf\n", ratio[0]);
	//printf ("s0 ratio: %lf\n", ratio[1]);

	while (ratio[0] >= 1e-6 || ratio[1] >= 1e-6)
	{
		size[0] = 1.05*size[0];
		size[1] = 1.05*size[1];
		calSize (acfStructure, size, ratio);
		//printf ("f0 ratio: %lf\n", ratio[0]);
		//printf ("s0 ratio: %lf\n", ratio[1]);
	}

	acfStructure->size[0] = acfStructure->winsize*size[0];
	acfStructure->size[1] = acfStructure->winsize*size[1];
	//acfStructure->size[0] = size[0];
	//acfStructure->size[1] = size[1];

	return 0;
}

int calSize (acfStruct *acfStructure, double *size, double *ratio)
{
	int i;
	double rand = acfStructure->phaseGradient;
	double steps = acfStructure->steps;
	double stepf = acfStructure->stepf;

	// for calculate box size
	int nf = (int)(size[0]*2/stepf)+1;
	int ns = (int)(size[1]*2/steps)+1;

	float  smax;
	float  fmax;
	float  *s, *acfs;
	float  *f, *acff;
	//float  s[ns], acfs[ns];
	//float  f[nf], acff[nf];
	//double c; // value at the center

	s = (float *)malloc(sizeof(float)*ns);
	f = (float *)malloc(sizeof(float)*nf);
	acfs = (float *)malloc(sizeof(float)*ns);
	acff = (float *)malloc(sizeof(float)*nf);

	for (i = 0; i < ns; i++)
	{
		s[i] = -size[1]+i*steps;
	}

	for (i = 0; i < nf; i++)
	{
		f[i] = -size[0]+i*stepf;
	}

	//c = exp(-pow((pow(fabs(s[(int)(ns/2)]+2.0*rand*0.4*f[(int)(nf/2)]),2.5)+pow(fabs(f[(int)(nf/2)]),1.5)),2.0/3.0));

	for (i = 0; i < nf; i++)
	{
		acff[i] = exp(-pow((pow(fabs(s[0]+2.0*rand*0.4*f[i]),2.5)+pow(fabs(f[i]),1.5)),2.0/3.0));
	}

	for (i = 0; i < ns; i++)
	{
		acfs[i] = exp(-pow((pow(fabs(s[i]+2.0*rand*0.4*f[0]),2.5)+pow(fabs(f[0]),1.5)),2.0/3.0));
	}

	smax = find_peak_value (ns, acfs);
	fmax = find_peak_value (nf, acff);

	ratio[0] = fmax;
	ratio[1] = smax;

	free(s);
	free(f);
	free(acfs);
	free(acff);

	return 0;
}

float find_peak_value (int n, float *s)
{
	int i;
	float temp[n];

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	float a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
															
		c = (a >= b ? a : b);
		
		temp[i+1] = c;
	
	}
	
	return temp[n-1];
}
						
void preAllocateMemory (acfStruct *acfStructure)
//void preAllocateMemory (acfStruct *acfStructure, controlStruct *control)
{
	double bw, f0, tint, t0;
	int nchn, nsubint;

	double steps;
	double stepf;

	double size[2]; // sampling boundary

	int nf, ns;

	//acfStructure->n = control->n; 
	//acfStructure->cFlux = control->cFlux; // mJy
	//acfStructure->whiteLevel = control->whiteLevel; // mJy
	//acfStructure->cFreq = control->cFreq; // MHz
	//acfStructure->bw = fabs(control->chanBW*control->nchan); // MHz
	//acfStructure->f0 = control->scint_freqbw;  // MHz
	//acfStructure->tint = control->nsub*control->tsub;  // s
	//acfStructure->t0 = control->scint_ts; // s
	//acfStructure->nchn = control->nchan;
	//acfStructure->nsubint = control->nsub;
	//printf ("Scintillation bandwidth: %lf (MHz)\n", acfStructure->f0);
	//printf ("Scintillation time-scale: %lf (s)\n", acfStructure->t0);

	bw = acfStructure->bw;
	f0 = acfStructure->f0;
	tint = acfStructure->tint;
	t0 = acfStructure->t0;

	nchn = acfStructure->nchn;
	nsubint = acfStructure->nsubint;

	// take f0 < bw_chn into account
	if (f0 >= (bw/nchn))
	{
		stepf = (bw/f0)/nchn;
	}
	else
	{
		stepf = 1.0;
	}
	
	if (t0 >= (tint/nsubint))
	{
		steps = (tint/t0)/nsubint;
	}
	else
	{
		steps = 1.0;
	}
	
	acfStructure->steps = steps;
	acfStructure->stepf = stepf;
	//printf ("%lf\n", steps);
	//printf ("%lf\n", stepf);
	
	windowSize (acfStructure, size);
	//printf ("f0 size: %lf\n", size[0]);
	//printf ("s0 size: %lf\n", size[1]);

	nf = (int)(acfStructure->size[0]*2/stepf)+1;
	ns = (int)(acfStructure->size[1]*2/steps)+1;
	//nf = (int)(size[0]*2/stepf)+1;
	//ns = (int)(size[1]*2/steps)+1;
	//printf ("ns: %d\n", ns);
	//printf ("nf: %d\n", nf);

	acfStructure->ns = ns;
	acfStructure->nf = nf;
}

int readParams(char *fname, char *dname, int n, controlStruct *control)
//int readParams(char *fname, char *oname, char *dname, int n, controlStruct *control)
{
	FILE *fin;
	char param[1024];
	int endit=-1;
	int finished=0;

	// define the output file name
	//strcpy(control->oname,oname);
	strcpy(control->dname,dname);

	control->n = n;

	///////////////////////////////////////
	if ((fin=fopen(fname,"r"))==NULL)
	{
		printf ("Can't open file!\n");
		exit(1);
	}

	//printf("Reading parameters...\n");

      	// Find the start observation
	while (!feof(fin))
	{
		if (fscanf(fin,"%s",param)==1)
		{
			if (strcasecmp(param,"START_OBS")==0)
			{
				endit=0;
				break;
			}
		}
		else 
			return 1;
	}

	if (endit==-1)
		return 1;

	do
	{
		fscanf(fin,"%s",param);
		if (strcasecmp(param,"END_OBS")==0)
			endit=1;
		else
		{
			if (strcasecmp(param,"SCINT_TS")==0)
				fscanf(fin,"%lf",&(control->scint_ts));
			else if (strcasecmp(param,"SCINT_BW")==0)
				fscanf(fin,"%lf",&(control->scint_freqbw));	  
			else if (strcasecmp(param,"T")==0)
      				fscanf(fin,"%lf",&(control->tint));
			else if (strcasecmp(param,"BW")==0)
      				fscanf(fin,"%lf",&(control->bw));
			else if (strcasecmp(param,"CFREQ")==0)
      				fscanf(fin,"%lf",&(control->cFreq));
			//else if (strcasecmp(param,"TSUB")==0)
      			//	fscanf(fin,"%lf",&(control->tsub));
			//else if (strcasecmp(param,"CHAN_BW")==0)
      			//	fscanf(fin,"%lf",&(control->chanBW));
			else if (strcasecmp(param,"NCHAN")==0)
				fscanf(fin,"%d",&(control->nchan));
			else if (strcasecmp(param,"NSUB")==0)
      				fscanf(fin,"%d",&(control->nsub));
			else if (strcasecmp(param,"WHITE_LEVEL")==0)
      				fscanf(fin,"%lf",&(control->whiteLevel));
			else if (strcasecmp(param,"CFLUX")==0)
			      	fscanf(fin,"%lf",&(control->cFlux));
			else if (strcasecmp(param,"WINSIZE")==0)
			      	fscanf(fin,"%lf",&(control->winsize));
		}
	} while (endit==0);

	//////////////////////////////////////////////////////////////////////
	//if (control->tsys != 0.0 && control->tsky != 0.0 && control->gain != 0.0 && control->whiteLevel == 0)
	//{
	//	control->radioNoise = (control->tsys+control->tsky)/(control->gain)/sqrt(2.0*(control->tsubRequested/control->nbin)*(fabs(control->obsBW)/control->nchan));
	//}
	//else if (control->tsys == 0.0 && control->tsky == 0.0 && control->gain == 0.0 && control->whiteLevel != 0)
	//{
	//	control->radioNoise = control->whiteLevel;
	//}
	//else 
	//{
	//	printf ("Double definiation of radio-meter noise!\n");
	//	exit (1);
	//}
	//printf ("Nchan: %d; Tsys: %lf; Tsky: %lf; Gain: %lf; Radio-meter noise: %lf mJy\n", control->nchan, control->tsys, control->tsky, control->gain, control->radioNoise);

	if (fclose(fin))
	{
		printf ("Can't close file.\n");
		exit(1);
	}

	return finished;
}

void initialiseControl(controlStruct *control)
{
	//strcpy(control->primaryHeaderParams,"UNKNOWN");
	//strcpy(control->exact_ephemeris,"UNKNOWN");
	//strcpy(control->src,"UNKNOWN");
	//strcpy(control->oname,"UNKNOWN");
	strcpy(control->dname,"1/xs");
	control->tsub = 0;
	
	control->n = 1; // simulate 1 dynamic spectrum by default

	// Standard defaults
	control->tint = 10000;
	control->bw = 10000;
	control->nchan = 100;
	control->nsub = 10000;
	control->nsub_win = 100;
	control->chanBW = control->bw/control->nchan;   // MHz
	control->tsub = control->tint/control->nsub;  // second

	control->cFreq = 1400.0;
	control->whiteLevel = 0;   // mJy
	control->scint_ts  = 0.0;  // second
	//control->scint_ts0  = 0.0;  // second
	//control->scint_ts1  = 0.0;  // second
	//control->scint_ts_step  = 0.0;  // second
	control->scint_freqbw = 0.0;   // MHz
	//control->scint_freqbw0 = 0.0;   // MHz
	//control->scint_freqbw1 = 0.0;   // MHz
	//control->scint_f_step = 0.0;   // MHz
	control->cFlux = 0.0;   // mJy

	control->winsize = 2.5;   
}

//void heatMap (acfStruct *acfStructure, char *dname)
//{
//	//int i,j;                     
//	//int dimx = acfStructure.ns;
//	//int dimy = acfStructure.nf; // dimensions 
//	//float tab[dimx*dimy];       // value
//	char caption[1024];
//	sprintf (caption, "%s %.2f %s %s %.2f %s %s %.2f %s", "Freq:", acfStructure->cFreq, "MHz", "BW:", acfStructure->bw, "MHz", "Length:", acfStructure->tint, "s");
//  
//	float zmin,zmax;            /* min et max des valeurs de la fonction */
//	float tr[6];                /* matrice utilisee par pgimag */
//
//	int dimx = acfStructure->nsubint;
//	int dimy = acfStructure->nchn;
//	double bw = acfStructure->bw;
//  
//
//	float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
//	float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
//	float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
//	float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
//
//	double f1 = acfStructure->cFreq-bw/2.0-2.0*bw/dimy; // MHz
//	double f2 = acfStructure->cFreq+bw/2.0-2.0*bw/dimy; // MHz
//	//printf ("f1 f2: %lf %lf\n", f1, f2);
//
//	zmin=0; 
//	zmax=find_peak_value (dimx*dimy, acfStructure->dynPlot[0]);
//	//double f1 = 1241; // MHz
//	//double f2 = 1497; // MHz
//	/*The transformation matrix TR is used to calculate the world
//	coordinates of the center of the "cell" that represents each
//	array element. The world coordinates of the center of the cell
//	corresponding to array element A(I,J) are given by:
//	X = TR(1) + TR(2)*I + TR(3)*J
//	Y = TR(4) + TR(5)*I + TR(6)*J
//	Usually TR(3) and TR(5) are zero -- unless the coordinate
//	transformation involves a rotation or shear.  The corners of the
//	quadrilateral region that is shaded by PGIMAG are given by
//	applying this transformation to (I1-0.5,J1-0.5), (I2+0.5, J2+0.5).*/
//  
//	//tr[0]=0;
//	//tr[1]=(float)(dimy)/dimx;
//	//tr[2]=0;
//	//tr[3]=0;
//	//tr[4]=0;
//	//tr[5]=1;
//  
//	tr[0]=-0.5;
//       	tr[1]=1;
//       	tr[2]=0;
//      	tr[3]=f2+0.5;
//      	tr[4]=0;
//      	tr[5]=-bw/dimy;
//
//	// plot 
//	//cpgbeg(0,"?",1,1);
//	cpgbeg(0,dname,1,1);
//	//cpgbeg(0,"2/xs",1,1);
//      	cpgsch(1.2); // set character height
//      	cpgscf(2); // set character font
//	//cpgswin(0,dimx,164,132); // set window
//	cpgswin(0,dimx,f2,f1); // set window
//	//cpgsvp(0.1,0.9,0.1,0.9); // set viewport
//      	//cpgenv(1,dimx,f1,f2,0,0); // set window and viewport and draw labeled frame
//	cpgbox("BCTSIN",4,4,"BCTSIN",16,8);
//	//cpgbox("BCTSIN",10,5,"BCTSIN",50,5);
//      	cpglab("Subintegration","Frequency (MHz)",caption);
//      	//cpglab("Subintegration","Frequency (MHz)","Freq: 150.0 MHz BW: -32.000 MHz Length: 960.0 s");
//	//cpgtick(0,0,1024,0,1.0/64,0.1,0.2,0,0,"1");
//	//palett(3, -0.4, 0.3);
//	//cpgimag(tab,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//	cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
//	cpgimag(acfStructure->dynPlot[0],dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//	//cpggray(acfStructure->dynPlot,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//
//	cpgend();
//} 
//
//void palett(int TYPE, float CONTRA, float BRIGHT)
//{
////-----------------------------------------------------------------------
//// Set a "palette" of colors in the range of color indices used by
//// PGIMAG.
////-----------------------------------------------------------------------
//	float GL[] = {0.0, 1.0};
//	float GR[] = {0.0, 1.0};
//	float GG[] = {0.0, 1.0};
//	float GB[] = {0.0, 1.0};
//	float RL[] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
//	float RR[] = { 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
//	float RG[] = { 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
//	float RB[] = { 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};
//	float HL[] = {0.0, 0.2, 0.4, 0.6, 1.0};
//	float HR[] = {0.0, 0.5, 1.0, 1.0, 1.0};
//	float HG[] = {0.0, 0.0, 0.5, 1.0, 1.0};
//	float HB[] = {0.0, 0.0, 0.0, 0.3, 1.0};
//	float WL[] = {0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0};
//	float WR[] = {0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0};
//	float WG[] = {0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0};
//	float WB[] = {0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0};
//	float AL[] = {0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0};
//	float AR[] = {0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
//	float AG[] = {0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0};
//	float AB[] = {0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//      
//	if (TYPE == 1)
//	{   
//		//-- gray scale
//		cpgctab(GL, GR, GG, GB, 2, CONTRA, BRIGHT);
//	}
//	else if (TYPE == 2) 
//	{
//		//-- rainbow
//		cpgctab(RL, RR, RG, RB, 9, CONTRA, BRIGHT);
//	}
//	else if (TYPE == 3) 
//	{
//		//-- heat
//		cpgctab(HL, HR, HG, HB, 5, CONTRA, BRIGHT);
//	}
//	else if (TYPE == 4) 
//	{
//		//-- weird IRAF
//		cpgctab(WL, WR, WG, WB, 10, CONTRA, BRIGHT);
//	}
//	else if (TYPE == 5) 
//	{
//		//-- AIPS
//		cpgctab(AL, AR, AG, AB, 20, CONTRA, BRIGHT);
//	}
//}

//int plotDynSpec (char *pname, char *dname)
//{
//	FILE *fin;
//	char start[128];
//	double tint,bw,cFreq;
//	int nsub,nchn;
//	float val;
//	int n1, n2;
//	float *dynSpec;
//
//	int i;
//
//	int dimx;
//	int dimy;
//	float zmin,zmax;            /* min et max des valeurs de la fonction */
//	float tr[6];                /* matrice utilisee par pgimag */
//
//	float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
//	float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
//	float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
//	float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
//	char caption[1024];
//
//	if ((fin=fopen(pname,"r"))==NULL)
//	{
//		printf ("Can't open dynamic spectrum!\n");
//		exit(1);
//	}
//
//      	// Find the start of dynamic spectrum, which contains basic info
//	while (!feof(fin))
//	{
//		if (fscanf(fin,"%s %d %d %lf %lf %lf",start,&nsub,&nchn,&bw,&tint,&cFreq)==6)
//		{
//			if (strcasecmp(start,"START")==0)
//			{
//				break;
//			}
//		}
//		//else 
//		//	return 1;
//	}
//
//	sprintf (caption, "%s %.2f %s %s %.2f %s %s %.2f %s", "Freq:", cFreq, "MHz", "BW:", bw, "MHz", "Length:", tint, "s");
//
//	dynSpec = (float*)malloc(sizeof(float)*nsub*nchn);
//
//	i = 0;
//	while (fscanf(fin,"%d %d %f", &n1, &n2, &val)==3)
//	{
//		dynSpec[i] = val;
//		i++;
//	}
//
//	if (fclose(fin))
//	{
//		printf ("Can't close dynamic spectrum!\n");
//		exit(1);
//	}
//
//	dimx = nsub;
//	dimy = nchn;
//  
//
//	zmin=0; 
//	zmax=find_peak_value (dimx*dimy, dynSpec);
//
//	double f1 = cFreq-bw/2.0-2.0*bw/dimy; // MHz
//	double f2 = cFreq+bw/2.0-2.0*bw/dimy; // MHz
//	//double f1 = 1241; // MHz
//	//double f2 = 1497; // MHz
//	/*The transformation matrix TR is used to calculate the world
//	coordinates of the center of the "cell" that represents each
//	array element. The world coordinates of the center of the cell
//	corresponding to array element A(I,J) are given by:
//	X = TR(1) + TR(2)*I + TR(3)*J
//	Y = TR(4) + TR(5)*I + TR(6)*J
//	Usually TR(3) and TR(5) are zero -- unless the coordinate
//	transformation involves a rotation or shear.  The corners of the
//	quadrilateral region that is shaded by PGIMAG are given by
//	applying this transformation to (I1-0.5,J1-0.5), (I2+0.5, J2+0.5).*/
//  
//	//tr[0]=0;
//	//tr[1]=(float)(dimy)/dimx;
//	//tr[2]=0;
//	//tr[3]=0;
//	//tr[4]=0;
//	//tr[5]=1;
//  
//	tr[0]=-0.5;
//       	tr[1]=1;
//       	tr[2]=0;
//      	tr[3]=f2+0.5;
//      	tr[4]=0;
//      	tr[5]=-bw/dimy;
// 
//	// plot 
//	//cpgbeg(0,"?",1,1);
//	cpgbeg(0,dname,1,1);
//	//cpgbeg(0,"2/xs",1,1);
//      	cpgsch(1.2); // set character height
//      	cpgscf(2); // set character font
//	cpgswin(0,dimx,f2,f1); // set window
//	//cpgsvp(0.1,0.9,0.1,0.9); // set viewport
//      	//cpgenv(1,dimx,f1,f2,0,0); // set window and viewport and draw labeled frame
//	cpgbox("BCTSIN",4,4,"BCTSIN",16,8);
//      	cpglab("Subintegration","Frequency (MHz)",caption);
//	//cpgtick(0,0,1024,0,1.0/64,0.1,0.2,0,0,"1");
//	//palett(3, -0.4, 0.3);
//	//cpgimag(tab,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//	//cpggray(dynSpec,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//	cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
//	cpgimag(dynSpec,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//
//	cpgend();
//
//	free(dynSpec);
//
//	return 0;
//} 

//int qualifyVar (acfStruct *acfStructure, noiseStruct *noiseStructure, controlStruct *control)
//{
//	int i;
//	int n = acfStructure->n;
//	int n_n = noiseStructure->n;
//
//	float m0;
//	//float *m;
//
//	float *var; 
//	//float *var_n;
//	//float varVar, meanVar;
//	//float varVar_n, meanVar_n;
//
//	int nsub = control->nsub;
//	int nchan = control->nchan;
//
//	int num;
//
//	//m = (float*)malloc(sizeof(float)*n);
//	var = (float*)malloc(sizeof(float)*n);
//	//var_n = (float*)malloc(sizeof(float)*n_n);
//
//	for (i=0; i<n; i++)
//	{
//		//m[i] = moduIndex (acfStructure->dynPlot[i], nsub*nchan);
//		var[i] = variance (acfStructure->dynPlot[i], nsub*nchan);
//		//printf ("%f\n", var[i]);
//		//printf ("%f \n", m[i]);
//	}
//
//	num = 0;
//	for (i=0; i<n; i++)
//	{
//		if (var[i] >= noiseStructure->detection)
//		{
//			num++;
//		}
//	}
//	acfStructure->probability = (float)(num)/n;
//	//printf ("%d %d %f %f\n", num, n, (float)(num)/n, acfStructure->probability);
//
//	//for (i=0; i<n_n; i++)
//	//{
//	//	var_n[i] = variance (noiseStructure->noisePlot[i], nsub*nchan);
//	//	//fprintf (fin2, "%f\n", var_n[i]);
//	//}
//
//	/*
//	if (fclose(fin1))
//	{
//		printf ("Can't close file...\n");
//		exit(1);
//	}
//
//	if (fclose(fin2))
//	{
//		printf ("Can't close file...\n");
//		exit(1);
//	}
//	*/
//
//	//m0 = 0.0;
//	//meanVar = 0.0;
//	//for (i=0; i<n; i++)
//	//{
//	//	m0 += m[i];
//	//	meanVar += var[i];
//	//}
//	//m0 = m0/n;
//	//meanVar = meanVar/n;
//
//	//meanVar_n = 0.0;
//	//for (i=0; i<n_n; i++)
//	//{
//	//	meanVar_n += var_n[i];
//	//}
//	//meanVar_n = meanVar_n/n_n;
//
//	//varVar = variance (var, n);
//	//varVar_n = variance (var_n, n_n);
//
//	//printf ("Results: %f %f %f %f %f\n", m0, meanVar, varVar, meanVar_n, varVar_n);
//
//	/*
//	maxV = find_max_value(n,var);
//	minV = find_min_value(n,var);
//	maxVN = find_max_value(n_n,var_n);
//	minVN = find_min_value(n_n,var_n);
//	////////////////////////////
//	// make histogram
//	if (control->noplot != 1)
//	{
//		cpgbeg(0,control->dname,1,1);
//
//		xHis = (float*)malloc(sizeof(float)*step);
//		val = (float*)malloc(sizeof(float)*step);
//		xHisN = (float*)malloc(sizeof(float)*step);
//		valN = (float*)malloc(sizeof(float)*step);
//		histogram (var, n, xHis, val, minV-2*(maxV-minV)/step, maxV+2*(maxV-minV)/step, step);
//		histogram (var_n, n_n, xHisN, valN, minVN-2*(maxVN-minVN)/step, maxVN+2*(maxVN-minVN)/step, step);
//
//      		cpgsch(1); // set character height
//      		cpgscf(1); // set character font
//
//		// find the max and min
//		max1 = find_max_value(step,val);
//		max2 = find_max_value(step,valN);
//		max = (max1 >= max2 ? max1 : max2);
//      		//cpgenv(-5,5,0,4500,0,1); // set window and viewport and draw labeled frame
//      		cpgenv(minVN-1, maxV+1, 0, max+0.1*max, 0, 1); // set window and viewport and draw labeled frame
//
//		sprintf(caption, "%s", "Variance histogram");
//      		cpglab("Variance","Number",caption);
//		cpgbin(step,xHis,val,0);
//		cpgsci(2);
//		cpgbin(step,xHisN,valN,0);
//		cpgsci(3);
//		cpgline(100, xThreshold, yThreshold);
//		///////////////////////////////////////////////////////
//		cpgend();
//
//		free(xHis);
//		free(val);
//		free(xHisN);
//		free(valN);
//	}
//	*/
//
//	//free(m);
//	free(var);
//	//free(var_n);
//	//free(flux0);
//	//free(flux);
//
//	return 0;
//}

float chiSquare (float *data, int n, float noise)
{
	int i;

	float ave;
	float chiS;
	
	ave = 0.0;
	for (i=0; i<n; i++)
	{
		ave += data[i];
	}
	ave = ave/n;

	chiS = 0.0;
	for (i=0; i<n; i++)
	{
		chiS += pow(data[i]-ave,2)/pow(noise,2);
	}

	return chiS/(n-1);
}

float moduIndex (float *data, int n)
{
	int i;

	float ave, devi;
	float m;

	ave = 0.0;
	for (i=0; i<n; i++)
	{
		ave += data[i];
	}
	ave = ave/n;

	devi = 0.0;
	for (i=0; i<n; i++)
	{
		devi += pow(data[i]-ave,2);
	}
	devi = sqrt(devi/n);

	m = devi/ave;

	return m;
}

float variance (float *data, int n)
{
	int i;

	float ave, devi;

	ave = 0.0;
	for (i=0; i<n; i++)
	{
		ave += data[i];
	}
	ave = ave/n;

	devi = 0.0;
	for (i=0; i<n; i++)
	{
		devi += pow(data[i]-ave,2);
	}
	devi = devi/n;

	return devi;
}

int histogram (float *data, int n, float *x, float *val, float low, float up, int step)
{
	int i,j,count;
	float width;
	float *temp;

	temp = (float*)malloc(sizeof(float)*(step+1));

	width = (up-low)/step;
	for (i=0; i<step; i++)
	{
		x[i] = low + i*width + width/2.0;
	}

	for (i=0; i<=step; i++)
	{
		temp[i] = low + i*width;
	}

	for (i=0; i<step; i++)
	{
		count = 0;
		for (j=0; j<n; j++)
		{
			if (data[j]>=temp[i] && data[j]<temp[i+1])
			{
				count += 1;
			}
		}
		//val [i] = count;
		val[i] = (float)(count)/n;
	}

	free(temp);
	return 0;
}

float find_max_value (int n, float *s)
{
	int i;
	float *temp;

	temp = (float *)malloc(sizeof(float)*n);

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	float a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
															
		c = (a >= b ? a : b);
		
		temp[i+1] = c;
	
	}
	
	c = temp[n-1];
	free(temp);

	return c;
}

float find_min_value (int n, float *s)
{
	int i;
	float *temp;

	temp = (float *)malloc(sizeof(float)*n);

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	float a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
															
		c = (a <= b ? a : b);
		
		temp[i+1] = c;
	
	}
	
	c = temp[n-1];
	free(temp);

	return c;
}

int readDissNum (char *Tname)
{
	int nt;
	FILE *fint;
	double temp;

	if ((fint=fopen(Tname, "r"))==NULL)
	{
		printf ("Can't open scintillation time-scale...\n");
		exit(1);
	}

	nt = 0;
	while (fscanf(fint, "%lf", &temp) == 1)
	{
		nt++;
	}

	if (fclose(fint))
	{
		printf ("Can't close scintillation time-scale...\n");
		exit(1);
	}

	return nt;
}



void readDiss (char *Tname, char *Fname, double *tdiss, double *fdiss)
{
	int i;
	FILE *fint, *finf;
	double temp;

	///////////////////////////////////////////////////////////////////////
	if ((fint=fopen(Tname, "r"))==NULL)
	{
		printf ("Can't open scintillation time-scale...\n");
		exit(1);
	}

	if ((finf=fopen(Fname, "r"))==NULL)
	{
		printf ("Can't open scintillation bandwidth...\n");
		exit(1);
	}

	i = 0;
	while (fscanf(fint, "%lf", &temp) == 1)
	{
		tdiss[i] = temp;
		i++;
	}

	i = 0;
	while (fscanf(finf, "%lf", &temp) == 1)
	{
		fdiss[i] = temp;
		i++;
	}

	if (fclose(fint))
	{
		printf ("Can't close scintillation time-scale...\n");
		exit(1);
	}

	if (fclose(finf))
	{
		printf ("Can't close scintillation time-scale...\n");
		exit(1);
	}
}
