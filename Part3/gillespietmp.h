// this header file contains a rather more general framework for gillespie simulation than in the previous parts

#define MAXR 100
#define MAXN 100

#define RND drand48()

// structure containing descriptions of reactions
typedef struct
{
  int nreacs, nprods;
  int reacs[5], prods[5];
  int rstoc[5], pstoc[5];
} Reaction;

// structure containing parameterisation for a given gillespie simulation (reactions, number of iterations, time window, etc)
typedef struct tagGillespie
{
  int N, M;
  Reaction R[MAXR];
  int init[MAXN];
  double c[MAXR];
  int nrun;
  int numt;
} Gillespie;

// quick factorial function
int Factorial(int n)
{
  if(n <= 1) return 1;
  else return n*Factorial(n-1);
}

// quick binomial function
int Binomial(int n, double p)
{
  int x;
  int i;

  x = 0;
  for(i = 0; i < n; i++)
    if(RND < p) x++;
  return x;
}

// work out combinatoric terms in concentration and stochiometry for gillespie simulation
double Combinatoric(double conc, int stoc)
{
  double tmp;
  int i;

  if(conc <= 0) 
    { 
      if(conc < 0) 
	{ 
	  return 0;
	  printf("odd error early1\n"); 
	  exit(0); 
	} 
      return 0; 
    }
  tmp = 1;
  for(i = 0; i < stoc; i++)
    tmp *= (conc-i);
  return tmp/Factorial(stoc);
}

// run a gillespie simulation without cell divisions
double GillespieSimulate(Gillespie G, int verbose, FILE *fp)
{
  int x[MAXN], oldx[MAXN];
  double h[MAXR];
  double a[MAXR];
  double r1, r2;

  double tau;
  double dt = 10, nextt;
  int tmpt, lastrec;
  double t;
  int i, j, k;
  double tmp;
  double az0;
  double sum;
  int n;
  int reaction;

  int run;
  double mean, sd;
  double mean0, mean1;

  char string[40];
  double *rec;

  double stat, total;

  // allocate memory for array recording copy number statistics
  rec = (double*)malloc(sizeof(double)*(G.numt+10)*G.nrun*G.N);

  // loop through iterations
  for(run = 0; run < G.nrun; run++)
    {
      printf("Traj %i\n", run);

      // initialise the system
      n = 0; 
      lastrec = 0;
      t = 0; 
      for(i = 0; i < G.N; i++)
	{
	  x[i] = G.init[i];
	  rec[((G.numt+10)*G.nrun)*i+(G.numt+10)*run+0] = x[i];
	}

      // this loop encloses time progression
      for(;t < G.numt;)
	{
	  for(i = 0; i < G.N; i++)
	    oldx[i] = x[i];
	  for(j = 0; j < G.M; j++)
	    {
	      h[j] = 1.0;
	      for(k = 0; k < G.R[j].nreacs; k++)
		{
		  h[j] *= Combinatoric(x[G.R[j].reacs[k]], G.R[j].rstoc[k]);
		}
	    }
	  // compute statistics for each reaction propensity
	  az0 = 0;
	  for(j = 0; j < G.M; j++)
	    {
	      a[j] = h[j]*G.c[j];
	      az0 += a[j];
	    }
	  // make our random choice
	  r1 = RND; r2 = RND;
	  tau = (1.0/az0)*log(1.0/r1);
	  tmp = r2*az0; sum = 0;
	  for(j = 0; sum < tmp && j < G.M; j++)
	    sum += a[j];

	  // we're not finished: so update populations accordingly
	  reaction = j-1;

	  if(reaction == -1)
	    {
	      tau = G.numt;
	    }
	  else
	    {
	      for(j = 0; j < G.R[reaction].nreacs; j++)
		x[G.R[reaction].reacs[j]] -= G.R[reaction].rstoc[j];

	      for(j = 0; j < G.R[reaction].nprods; j++)
		x[G.R[reaction].prods[j]] += G.R[reaction].pstoc[j];	
	    }

	  // increment time according to reaction delay
	  //	 	 	  	  	  printf("Was at %i %i  (%.4f), did reaction %i, now at %i %i  (%.4f)\n", oldx[0], oldx[1], t, reaction, x[0], x[1], t+tau);
	  //	  printf("Was at %i %i [%i %i %i %i] (%.4f), did reaction %i, now at %i %i [%i %i %i %i] (%.4f)\n", oldx[0], oldx[1], oldx[4], oldx[5], oldx[24], oldx[25], t, reaction, x[0], x[1], x[4], x[5], x[24], x[25], t+tau);


	  t += tau;
	  if(t > G.numt) t = G.numt;

	  // record copy number statistics
	  for(tmpt = floor(lastrec+1); tmpt < t; tmpt ++)
	    {
	      for(i = 0; i < G.N; i++)
		rec[((G.numt+10)*G.nrun)*i+(G.numt+10)*run+tmpt] = oldx[i];
	    }
	  lastrec = t;
	}
    }


  // the simulation code ends here -- we're left with arrays of copy number records
  // output to file if desired (set with verbose)
  for(t = 0; t < G.numt; t++)
    {
      if(verbose)   fprintf(fp, "%.4f ", t);
      for(i = 0; i < G.N; i++)
	{
	  mean = sd = 0;
	  for(run = 0; run < G.nrun; run++)
	    mean += rec[(G.numt+10)*G.nrun*i + (G.numt+10)*run + (int)t];
	  mean /= G.nrun;
	  for(run = 0; run < G.nrun; run++)
	    sd += pow(rec[(G.numt+10)*G.nrun*i + (G.numt+10)*run + (int)t]-mean, 2);
	  sd /= G.nrun-1; sd = sqrt(sd);
	  if(verbose)	  fprintf(fp, "%.4f %.4f ", mean, sd);
	}
      if(verbose)     fprintf(fp, "\n");
    }

  free(rec);

  return mean;
}

// run a gillespie simulation including binomial cell divisions
double GillespieSimulateDivisions(Gillespie G, int verbose, FILE *fp, double celltau, int crazyverbose)
{
  int x[MAXN], oldx[MAXN];
  double h[MAXR];
  double a[MAXR];
  double r1, r2;

  double tau;
  double dt = 10, nextt;
  int tmpt, lastrec;
  double t;
  int i, j, k;
  double tmp;
  double az0;
  double sum;
  int n;
  int reaction;

  int run;
  double mean, sd;
  double mean0, mean1;
  double nextdivision;

  char string[40];
  double *rec;
  FILE *fp1;

  double stat, total;
  int nhist = 2000;
  double *hist;
  int thisrec;

  // allocate memory for array recording copy number statistics
  rec = (double*)malloc(sizeof(double)*(G.numt+10)*G.nrun*G.N);
  hist = (double*)malloc(sizeof(double)*nhist);

  for(run = 0; run < G.nrun; run++)
    {
      if(run % 100 == 0)
      printf("Traj %i\n", run);

      // initialise the system
      n = 0; 
      lastrec = 0;
      t = 0; 
      nextdivision = celltau;
      for(i = 0; i < G.N; i++)
	{
	  x[i] = G.init[i];
	  rec[((G.numt+10)*G.nrun)*i+(G.numt+10)*run+0] = x[i];
	}

      // this loop encloses time progression
      for(;t < G.numt;)
	{
	  for(i = 0; i < G.N; i++)
	    oldx[i] = x[i];
	  for(j = 0; j < G.M; j++)
	    {
	      h[j] = 1.0;
	      for(k = 0; k < G.R[j].nreacs; k++)
		{
		  h[j] *= Combinatoric(x[G.R[j].reacs[k]], G.R[j].rstoc[k]);
		}
	    }
	  // compute statistics for each reaction propensity
	  az0 = 0;
	  for(j = 0; j < G.M; j++)
	    {
	      a[j] = h[j]*G.c[j];
	      az0 += a[j];
	    }
	  // make our random choice
	  r1 = RND; r2 = RND;
	  tau = (1.0/az0)*log(1.0/r1);
	  tmp = r2*az0; sum = 0;
	  for(j = 0; sum < tmp && j < G.M; j++)
	    sum += a[j];

	  reaction = j-1;

	  // if the next reaction won't occur before the next scheduled cell division
	  if(t + tau > nextdivision)
	    {
	      reaction = -2;
	      // increment time minutely past the division
	      tau = nextdivision-t+1e-8;
	      nextdivision += celltau;
	      // binomially partition cellular elements
	      for(i = 0; i < G.N; i++)
		{
		  x[i] = Binomial(x[i], 0.5);
		}
	    }
	  else
	    {
	      // perform the reaction if appropriate and update copy numbers
	  if(reaction == -1)
	    {
	      tau = G.numt;
	    }
	  else
	    {
	      for(j = 0; j < G.R[reaction].nreacs; j++)
		x[G.R[reaction].reacs[j]] -= G.R[reaction].rstoc[j];

	      for(j = 0; j < G.R[reaction].nprods; j++)
		x[G.R[reaction].prods[j]] += G.R[reaction].pstoc[j];	
	    }
	    }
	  // increment time according to reaction delay
	  //	  printf("Was at %i %i  (%.4f), did reaction %i, now at %i %i  (%.4f)\n", oldx[0], oldx[1], t, reaction, x[0], x[1], t+tau);
	  //	  printf("Was at %i %i [%i %i %i %i] (%.4f), did reaction %i, now at %i %i [%i %i %i %i] (%.4f)\n", oldx[0], oldx[1], oldx[4], oldx[5], oldx[24], oldx[25], t, reaction, x[0], x[1], x[4], x[5], x[24], x[25], t+tau);


	  t += tau;
	  if(t > G.numt) t = G.numt;

	  // record copy number statistics
	  for(tmpt = floor(lastrec+1); tmpt < t; tmpt ++)
	    {
	      for(i = 0; i < G.N; i++)
		rec[((G.numt+10)*G.nrun)*i+(G.numt+10)*run+tmpt] = oldx[i];
	    }
	  lastrec = t;
	}
    }


  // the simulation code ends here -- we're left with arrays of copy number records

  // output to file if desired (with verbose)
  // output distributional detail if desired (with crazyverbose)
  if(crazyverbose)
    fp1 = fopen("crazyverbosebase.txt", "w");

  for(t = 0; t < G.numt; t++)
    {
      if(verbose)   fprintf(fp, "%.4f ", t);
      if(crazyverbose)
	{
	  // build up and output a histogram for distributional detail
	  fprintf(fp1, "%.4f ", t);
	  for(i = 0; i < nhist; i++)
	    hist[i] = 0;
	  for(run = 0; run < G.nrun; run++)
	    {
	      thisrec = (int)rec[(G.numt+10)*G.nrun*0 + (G.numt+10)*run + (int)t];
	      //if(thisrec < 0)
		//
	      if(thisrec>0 && thisrec<nhist)
	      hist[thisrec]++;
	    }
	  for(i = 0; i < nhist; i++)
	    fprintf(fp1, "%.5f ", hist[i]/G.nrun);
	  fprintf(fp1, "\n");
	}

      for(i = 0; i < G.N; i++)
	{
	  // output summary statistics
	  mean = sd = 0;
	  for(run = 0; run < G.nrun; run++)
	    mean += rec[(G.numt+10)*G.nrun*i + (G.numt+10)*run + (int)t];
	  mean /= G.nrun;
	  for(run = 0; run < G.nrun; run++)
	    sd += pow(rec[(G.numt+10)*G.nrun*i + (G.numt+10)*run + (int)t]-mean, 2);
	  sd /= G.nrun-1; sd = sqrt(sd);
	  if(verbose)	  fprintf(fp, "%.4f %.4f ", mean, sd);
	}
      if(verbose)     fprintf(fp, "\n");
    }
      if(crazyverbose)
	fclose(fp1);

  free(rec);
  free(hist);

  return mean;
}

// run a gillespie simulation including subtractive cell divisions
// basically the same as above except for the partitioning step
double GillespieSimulateSubtractions(Gillespie G, int verbose, FILE *fp, double celltau, int eta, double tstep, int crazyverbose)
{
  int x[MAXN], oldx[MAXN];
  double h[MAXR];
  double a[MAXR];
  double r1, r2;

  double tau;
  double dt = 10, nextt;
  double tmpt, lastrec;
  double t;
  int i, j, k;
  double tmp;
  double az0;
  double sum;
  int n;
  int reaction;

  int run;
  double mean, sd;
  double mean0, mean1;
  double nextdivision;

  char string[40];
  double *rec;

  double stat, total;
  int postdiv[G.nrun];
  int ndiv;
  int nsteps = G.numt/tstep;

  FILE *fp1;
  double hist[30];
  int thisrec;

  // allocate memory for array recording copy number statistics
  rec = (double*)malloc(sizeof(double)*(nsteps+10)*G.nrun*G.N);

  srand48(34);

  for(run = 0; run < G.nrun; run++)
    {
      if(run % 100 == 0)
      printf("Traj %i\n", run);

      // initialise the system
      n = 0; 
      ndiv = 0;
      lastrec = 0;
      t = 0; 
      nextdivision = celltau;
      for(i = 0; i < G.N; i++)
	{
	  x[i] = G.init[i];
	  rec[((nsteps+10)*G.nrun)*i+(nsteps+10)*run+0] = x[i];
	}

      // this loop encloses time progression
      for(;t < G.numt;)
	{
	  for(i = 0; i < G.N; i++)
	    oldx[i] = x[i];
	  for(j = 0; j < G.M; j++)
	    {
	      h[j] = 1.0;
	      for(k = 0; k < G.R[j].nreacs; k++)
		{
		  h[j] *= Combinatoric(x[G.R[j].reacs[k]], G.R[j].rstoc[k]);
		}
	    }
	  // compute statistics for each reaction propensity
	  az0 = 0;
	  for(j = 0; j < G.M; j++)
	    {
	      a[j] = h[j]*G.c[j];
	      az0 += a[j];
	    }
	  // make our random choice
	  r1 = RND; r2 = RND;
	  tau = (1.0/az0)*log(1.0/r1);
	  tmp = r2*az0; sum = 0;
	  for(j = 0; sum < tmp && j < G.M; j++)
	    sum += a[j];

	  // we're not finished: so update populations accordingly
	  reaction = j-1;

	  if(t + tau > nextdivision)
	    {
	      reaction = -2;
	      tau = nextdivision-t+1e-8;
	      nextdivision += celltau;

	      // here is the difference from the binomial partitioning case
	      for(i = 0; i < G.N; i++)
		{
		  x[i] -= Binomial(2.*eta, 0.5);
		}
	      ndiv++;
	      if(ndiv == 3) 
		postdiv[run] = x[0];
	    }
	  else
	    {
	  if(reaction == -1)
	    {
	      tau = G.numt;
	    }
	  else
	    {
	      for(j = 0; j < G.R[reaction].nreacs; j++)
		x[G.R[reaction].reacs[j]] -= G.R[reaction].rstoc[j];

	      for(j = 0; j < G.R[reaction].nprods; j++)
		x[G.R[reaction].prods[j]] += G.R[reaction].pstoc[j];	
	    }
	    }
	  // increment time according to reaction delay
	  //	  printf("Was at %i %i  (%.4f), did reaction %i, now at %i %i  (%.4f)\n", oldx[0], oldx[1], t, reaction, x[0], x[1], t+tau);
	  //	  printf("Was at %i %i [%i %i %i %i] (%.4f), did reaction %i, now at %i %i [%i %i %i %i] (%.4f)\n", oldx[0], oldx[1], oldx[4], oldx[5], oldx[24], oldx[25], t, reaction, x[0], x[1], x[4], x[5], x[24], x[25], t+tau);


	  t += tau;
	  if(t > G.numt) t = G.numt;

	  // record copy number statistics
	  //	  for(tmpt = floor(lastrec+1); tmpt < t; tmpt ++)
	  for(tmpt = lastrec; tmpt < t; tmpt += tstep)
	    {
	      for(i = 0; i < G.N; i++)
		rec[((nsteps+10)*G.nrun)*i+(nsteps+10)*run+(int)(tmpt/tstep)] = oldx[i];
	    }
	  lastrec = t;
	}
    }


  // the simulation code ends here -- we're left with arrays of copy number records
  if(crazyverbose)
    fp1 = fopen("crazyverbose.txt", "w");

  for(t = 0; t < nsteps; t++)
    {
      if(verbose)   fprintf(fp, "%.4f ", (t+1)*tstep);
      if(crazyverbose)
	{
	  fprintf(fp1, "%.4f ", (t+1)*tstep);
	  for(i = 0; i < 30; i++)
	    hist[i] = 0;
	  for(run = 0; run < G.nrun; run++)
	    {
	      thisrec = (int)rec[(nsteps+10)*G.nrun*0 + (nsteps+10)*run + (int)t];
	      //if(thisrec < 0)
		//		printf("here");
	      hist[thisrec+10]++;
	    }
	  for(i = 0; i < 30; i++)
	    fprintf(fp1, "%.5f ", hist[i]/G.nrun);
	  fprintf(fp1, "\n");
	}

      for(i = 0; i < G.N; i++)
	{
	  mean = sd = 0;
	  for(run = 0; run < G.nrun; run++)
	    mean += rec[(nsteps+10)*G.nrun*i + (nsteps+10)*run + (int)t];
	  mean /= G.nrun;
	  for(run = 0; run < G.nrun; run++)
	    sd += pow(rec[(nsteps+10)*G.nrun*i + (nsteps+10)*run + (int)t]-mean, 2);
	  sd /= G.nrun-1; sd = sqrt(sd);
	  if(verbose)	  fprintf(fp, "%.4f %.4f ", mean, sd);
	}
      if(verbose)     fprintf(fp, "\n");
    }

  if(crazyverbose)
    fclose(fp1);

  free(rec);

  mean = sd = 0;
  for(run = 0; run < G.nrun; run++)
    {
      mean += postdiv[run];
    }
  mean /= G.nrun;
  for(run = 0; run < G.nrun; run++)
    {
      sd += pow(mean-postdiv[run],2);
    }
  sd /= G.nrun-1;

  return sd;
}

// add a reaction to the gillespie framework within G
// call this with labels of reactants and products
// see e.g. relaxed.h for examples of use
void SetupGillespie(int r, Gillespie *G, int nreacs, int nprods, ...)
{
  va_list valist;
  int i;

  G->R[r].nreacs = nreacs; G->R[r].nprods = nprods;
  va_start(valist, nprods);
  for(i = 0; i < nreacs; i++)
    G->R[r].reacs[i] = va_arg(valist, int);
  for(; i < nreacs+nprods; i++)
    G->R[r].prods[i-nreacs] = va_arg(valist, int);
  va_end(valist);
}
