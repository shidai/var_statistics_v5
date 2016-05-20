// make a color map of sensitivity as a function scintillation time-scale and bandwidth
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "detection_matched.h"
#include <mpi.h>

int main (int argc, char* argv[])
{
	int id;  //  process rank
	int p;   //  number of processes

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &p);

	//FILE *fin;
	controlStruct control;
	acfStruct acfStructure;
	noiseStruct noiseStructure;

	char fname[1024];   // read in parameter file
	char Tname[1024];   // read in noise

	int nNoise;
	double *noise;

	//char oname[1024];   // output file name
	//char pname[1024];   // file to plot
	char dname[1024];   // graphics device

	int i;
	int n = 1;  // number of dynamic spectrum to simulate

	int plotMode = 0;  // plot existing dynamic spectrum, off by default
	control.noplot = 1;  // don't show dynamic spectrum while simulationg, on by default

	// read options
	for (i=0;i<argc;i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
			strcpy(fname,argv[++i]);
			strcpy(Tname,argv[++i]);
			printf ("Parameters are in %s\n", fname);
		}
		//else if (strcmp(argv[i],"-o")==0)
		//{
		//	strcpy(oname,argv[++i]);
		//	//printf ("Dynamic spectrum is output into %s\n", oname);
		//}
		else if (strcmp(argv[i],"-n")==0)
		{
			n = atoi (argv[++i]);
			printf ("Number of dynamic spectrum to simulate %d\n", n);
		}
		//else if (strcmp(argv[i],"-p")==0) // just plot 
		//{
		//	strcpy(pname,argv[++i]);
		//	plotMode = 1;
		//	printf ("Plotting exiting dynamic spectrum.\n");
		//}
		//else if (strcmp(argv[i],"-noplot")==0)
		//{
		//	control.noplot = 1; // Don't show dynamic spetrum while simulationg
		//}
		//else if (strcmp(argv[i],"-dev")==0)
		//{
		//	strcpy(dname,argv[++i]);
		//	printf ("Graphic device: %s\n", dname);
		//}
	}

	nNoise = readNum (Tname);
	noise = (double *)malloc(sizeof(double)*nNoise);
	readNoise (Tname, noise);

	sprintf(dname, "%s", "1/xs");
	/*
	if ((fin=fopen(oname, "w"))==NULL)
	{
		printf ("Can't open file...\n");
		exit(1);
	}
	*/

	if (plotMode==0)
	{
		// Simulate dynamic spectrum
		// read parameters
		initialiseControl(&control);
		readParams (fname,dname,n,&control);
		printf ("Finished reading parameters.\n");

		allocateNoise (&noiseStructure, &control);
		calculateScintScale (&acfStructure, &control);

		for (i=id; i<nNoise; i+=p)
		//for (i=0; i<nNoise; i++)
		{
			control.whiteLevel = pow(10.0, noise[i]);
			calNoise (&noiseStructure, &control);
					
			calculateNDynSpec (&acfStructure, &control, &noiseStructure);
				
			//printf ("%lf %lf %lf\n", control.whiteLevel, acfStructure.meanM/sqrt(acfStructure.varM), acfStructure.meanV/sqrt(acfStructure.varV));
			//fflush (stdout);
		}
		deallocateNoise (&noiseStructure);
		deallocateMemory (&acfStructure);
	}
	//else
	//{
	//	plotDynSpec(pname, dname);
	//}

	/*
	if (fclose(fin))
	{
		printf ("Can't close file...\n");
		exit(1);
	}
	*/
		
	MPI_Finalize ();

	free(noise);
	return 0;
}

