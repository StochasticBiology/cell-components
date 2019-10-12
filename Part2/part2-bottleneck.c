// this code performs numeric stochastic simulation for comparison with the analytic theory concerning the bottleneck model in Fig. 3

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// flags for different partitioning regimes
#define BINOMIAL 0
#define CLUSTER 1
#define DETERMINISTIC 2
#define TELOMERE 3
#define TELOMERE_RANDOM 4

#define RND drand48()

// structure containing reaction parameterisation
typedef struct
{
  int nreacs, nprods;
  int reacs[5], prods[5];
  int rstoc[5], pstoc[5];
} Reaction;

// quick factorial function
int Factorial(int n)
{
  if(n <= 1) return 1;
  else return n*Factorial(n-1);
}

// used to compute reaction propensities in gillespie algorithm
double Combinatoric(double conc, int stoc)
{
  double tmp;
  int i;

  tmp = 1;
  for(i = 0; i < stoc; i++)
    tmp *= (conc-i);
  return tmp/Factorial(stoc);
}

// gaussian random number generator
double gsl_ran_gaussian(const double sigma)
{
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      x = -1 + 2 * drand48();
      y = -1 + 2 * drand48();

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
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

// number of timesteps and iterations
#define NUMT 240
#define NUMR 10000

// structure containing parameterisation for the bottleneck model
typedef struct tagP
{
  double lambda1, lambda2;
  int m0;
  double tau;
  int ndiv;
} Params;

// use gillespie algorithm to simulate bottleneck model
void Simulate(char *fname, Params P, int flag)
{
  double c[100];
  int x[100];
  int lastx;
  int M, N;
  double h[100];
  double a[100];
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
  Reaction r[100];
  int reaction;
  double hetrec[NUMR];
  double hethist[100];
  
  int run;

  char string[40];
  FILE *fp;
  double *rec0, *rec1;
  double mean0[NUMT], sd0[NUMT];
  double nclus;
  int watch = 0;

  srand48(495);

  // allocate memory
  rec0 = (double*)malloc(sizeof(double)*NUMT*(NUMR+10));

  // set up BID reactions

  // WTa -> WTa + WTa
  r[0].nreacs = 1; r[0].nprods = 2;
  r[0].reacs[0] = 0; r[0].prods[0] = 0; r[0].prods[1] = 0;
  r[0].rstoc[0] = 1; r[0].pstoc[0] = 1; r[0].pstoc[1] = 1;

  // WTa -> 0
  r[1].nreacs = 1; r[1].nprods = 1;
  r[1].reacs[0] = 0; r[1].prods[0] = 2;
  r[1].rstoc[0] = 1; r[1].pstoc[0] = 1;

  // 0 -> WTa 
  r[2].nreacs = 0; r[2].nprods = 1;
  r[2].reacs[0] = 2; r[2].prods[0] = 0;
  r[2].rstoc[0] = 0; r[2].pstoc[0] = 1;

  double nu1 = 0.01, nu2 = 0.01;
  double mu = 0;//0.001;
  double lambda1 = 0.1, lambda2 = 0.1;
  int taucc = P.tau;

  double nextdivision;
  double cellcyclelength;
  double p2onset, p3onset;
  int phase2, phase3;
  int n1 = 5;
  int n2 = 3;

  double residual;
  double res;
  double endsim;
  int ndiv;
  int eta;

  // loop through iterations
  for(run = 0; run < NUMR; run++)
    {
      printf("%i\n", run);

      // initialise system and timings
      n = 0; ndiv = 0;
      lastrec = 0;
      N = 2; M = 3;
      x[0] = P.m0;
      c[0] = P.lambda1; c[1] = 0; c[2] = 0;

      t = 0; 
      cellcyclelength = taucc;
      nextdivision = cellcyclelength;

      lastx = x[0];
      rec0[NUMT*run+0] = lastx;

      // start simulation
      for(i = 0; t < NUMT; i++)
	{
	  // catch pathological behaviour
	  if(x[0] > 2000000) {t = endsim; printf("went nuts\n"); break;}
	  az0 = 0;

	  // compute propensities for each reaction
	  for(j = 0; j < M; j++)
	    {
	      h[j] = 1.0;
	      for(k = 0; k < r[j].nreacs; k++)
		h[j] *= Combinatoric(x[r[j].reacs[k]], r[j].rstoc[k]);
	    }

	  for(j = 0; j < M; j++)
	    {
	      a[j] = h[j]*c[j];
	      az0 += a[j];
	    }

	  // draw random numbers to determine next reaction
	  r1 = drand48(); r2 = drand48();
	  tau = (1.0/az0)*log(1.0/r1);

	  // if next reaction won't occur before next division
	  if(t + tau > nextdivision)
	    {
	      t = nextdivision;

	      // record copy number
	      for(tmpt = lastrec+1; tmpt < t; tmpt ++)
		{
		  rec0[NUMT*run + tmpt] = lastx;
		}
	      nextdivision += cellcyclelength;

	      // apply binomial partitioning
	      x[0] = Binomial(x[0], 0.5);

	      x[0] = (x[0] < 0 ? 0 : x[0]);
	      lastx = x[0];
	      lastrec = t-1e-7;
	     
	      ndiv++;
	      // if we've now passed the minimum copy number time
	      if(ndiv >= P.ndiv)
		{
		  c[0] = P.lambda2;
		}

	    }
	  else      //perform the reaction
	    {
	      tmp = r2*az0; sum = 0;
	      for(j = 0; sum < tmp && j < M; j++)
		sum += a[j];
	      t += tau;

	      // record copy number
	      for(tmpt = lastrec+1; tmpt < t; tmpt ++)
		{
		  rec0[NUMT*run + tmpt] = lastx;
		}

	      reaction = j-1;

	      // make changes according to the reaction chosen
	      for(j = 0; j < r[reaction].nreacs; j++)
		x[r[reaction].reacs[j]] -= r[reaction].rstoc[j];

	      for(j = 0; j < r[reaction].nprods; j++)
		x[r[reaction].prods[j]] += r[reaction].pstoc[j];		

	      lastrec = t;
	      lastx = x[0];
	    }

	}


    }

  // output statistics to file
  fp = fopen(fname, "w");
  for(i = 0; i < NUMT; i++)
    {
      mean0[i] = sd0[i] = 0; 
      // compute mean
      for(j = 0; j < NUMR; j++)
	{
	  mean0[i] += rec0[NUMT*j+i];
	}
      mean0[i] /= NUMR; 
      // compute s.d.
      for(j = 0; j < NUMR; j++)
	{
	  sd0[i] += pow(mean0[i]-rec0[NUMT*j+i], 2);
	}
      sd0[i] /= NUMR; 
      sd0[i] = sqrt(sd0[i]);
      // output
      fprintf(fp, "%.1f %.5f %.5f\n", (double)i, mean0[i], sd0[i]*sd0[i]);
    }
  fclose(fp);


  free(rec0);

}

 
int main(int argc, char *argv[])
{
  double t[100], v[100];
  int n = 0;
  char string[40];
  Params P;
  double mlambda;
  double nuplus;
  int flag = 0;
  double rho;
  int expt;

  // initialise some parameters
  rho = 0.0693147;

  P.m0 = 100000;
  P.tau = 10;
  P.ndiv = 12;

  // parse command-line arguments to determine which parameterisation to run
  if(argc != 2)
    {
      printf("Which expt?\n");
      return 0;
    }
  expt = atoi(argv[1]);
  // apply appropriate output file and parameters for chosen experiment, and run
  switch(expt)
    {
    case 0:
      sprintf(string, "bneckbignuma.dat");
      P.lambda1 = rho; P.lambda2 = (2.*rho - P.lambda1);
      Simulate(string, P, 0);
      break;
    case 1:
      sprintf(string, "bneckbignumb.dat");
      P.lambda1 = 0.02; P.lambda2 = (2.*rho - P.lambda1);
      Simulate(string, P, 0); 
      break;
    case 2:
      sprintf(string, "bneckbignumc.dat");
      P.lambda1 = 0.03; P.lambda2 = (2.*rho - P.lambda1);
      Simulate(string, P, 0);
      break;
    case 3:
      sprintf(string, "bneckbignumd.dat");
      P.lambda1 = 0.04; P.lambda2 = (2.*rho - P.lambda1);
      Simulate(string, P, 0);
      break;
    case 4:
      sprintf(string, "bneckbignume.dat");
      P.lambda1 = 0.05; P.lambda2 = (2.*rho - P.lambda1);
      Simulate(string, P, 0);
      break;
    }

  return 0;
}

