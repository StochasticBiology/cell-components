/*** this code runs Gillespie simulations of birth-death-partitioning models with various partitioning regimes, in order to compare analytic results ****/
/*** takes one command-line parameter: partitioning regime to simulate (according to #defines below) -- facilitates running 4 in parallel ****/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// define references for each of the partitioning regimes
#define BINOMIAL 0
#define CLUSTER 1
#define DETERMINISTIC 2
#define TELOMERE 3
#define TELOMERE_ALT 4

// define the timespan of a simulation and the number of repeat simulations to compute statistics
#define NUMT 50
#define NUMR 100000

// structure containing the information required for a reaction in stochastic simulation
// the length-5 arrays are lazy: we have no reactions that need near this number of reactants/products
typedef struct
{
  int nreacs, nprods;
  int reacs[5], prods[5];
  int rstoc[5], pstoc[5];
} Reaction;

// compute a factorial
int Factorial(int n)
{
  if(n <= 1) return 1;
  else return n*Factorial(n-1);
}

// function to perform the combinatoric computation involved in determining reaction propensities in the gillespie algorithm
double Combinatoric(double conc, int stoc)
{
  double tmp;
  int i;

  tmp = 1;
  for(i = 0; i < stoc; i++)
    tmp *= (conc-i);
  return tmp/Factorial(stoc);
}

// function to produce a gaussian-distributed random number
// from the gsl library
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

// use gaussian random number and normal approximation to binomial to produce a binomial-like random variate 
int Binomial(double x)
{
  int max, min;
  x = nearbyint(x);
  double ans = (0.5*x + gsl_ran_gaussian(0.5*sqrt(x)));

  return nearbyint(ans);
  max = ((int)ans) + 1;
  min = ((int)ans);

  if(max - ans < ans - min) return max;
  else return min;
}

// structure to contain the parameterisation for a particular model: birth and death rates, size of a partitioned cluster, partitioning regime, initial copy number, cell cycle length, subtractive copy number loss
typedef struct tagP
{
  double lambda, nu;
  int nclus;
  int partition;
  int m0;
  double tau;
  int eta;
} Params;

// perform gillespie simulation of the model system with a given parameterisation
void Simulate(char *fname, Params P, int flag)
{
  double c[100];
  int x[100];
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

  srand48(1495);

  // allocate memory to record copy numbers
  rec0 = (double*)malloc(sizeof(double)*NUMT*(NUMR+10));

  // set up reactions in the system:
  // birth: m -> m + m
  r[0].nreacs = 1; r[0].nprods = 2;
  r[0].reacs[0] = 0; r[0].prods[0] = 0; r[0].prods[1] = 0;
  r[0].rstoc[0] = 1; r[0].pstoc[0] = 1; r[0].pstoc[1] = 1;

  // death: m -> 0
  r[1].nreacs = 1; r[1].nprods = 1;
  r[1].reacs[0] = 0; r[1].prods[0] = 2;
  r[1].rstoc[0] = 1; r[1].pstoc[0] = 1;

  double nu1 = 0.01, nu2 = 0.01;
  double mu = 0;
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

  // loop through number of repeat experiments
  for(run = 0; run < NUMR; run++)
    {
      printf("%i\n", run);

      // initialise parameters and initial conditions
      n = 0; ndiv = 0;
      lastrec = 0;
      N = 2; M = 2;
  x[0] = P.m0;
 c[0] = P.lambda; c[1] = P.nu;
 endsim = NUMT;

 // initialise time and cell cycle timer
      t = 0; 
      cellcyclelength = taucc;
      nextdivision = cellcyclelength;

      // this loop encloses time in the stochastic simulation
      for(i = 0; t < endsim; i++)
	{
	  // if copy number has diverged, output a message and jump to end of simulation
	  if(x[0] > 50000 || x[0] == 0)
	    {
	  if(x[0] > 50000) {printf("went nuts\n");}
	  if(x[0] == 0) {printf("bottomed out\n");}
	  t = endsim; 
	      for(tmpt = lastrec; tmpt < t; tmpt ++)
		{
		  rec0[NUMT*run + tmpt] = (x[0]);
		}
	      lastrec = t;
	      break;
	    }

	  // choose a reaction according to gillespie rules 
	  az0 = 0;
	  // loop through reactions calculating propensities
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
	  // pick random numbers to determine time until next reaction
	  r1 = drand48(); r2 = drand48();
	  tau = (1.0/az0)*log(1.0/r1);

	  if(t + tau > nextdivision)
	    {
	      // if a cell divisions happens sooner than the next reaction is projected to occur
	      // shift time to next division time and update records and cell cycle timer
	      t = nextdivision;
	      for(tmpt = lastrec; tmpt < t; tmpt ++)
		{
		  rec0[NUMT*run + tmpt] = (x[0]);
		}
	      lastrec = t;
	      nextdivision += cellcyclelength;

	      // partition existing population according to the partitioning regime specified in the parameter set
	      switch(P.partition)
		{
		case BINOMIAL: x[0] = Binomial(x[0]); break;
		case DETERMINISTIC: x[0] = floor(x[0]/2); break;
		case CLUSTER: x[0] = P.nclus*Binomial((double)x[0]/P.nclus); break;
		case TELOMERE: x[0] -= P.eta; break;
		case TELOMERE_ALT: x[0] -= P.eta; break;
		}

	      // this shouldn't happen: but regularise copy number to zero if we've fallen below
	      x[0] = (x[0] < 0 ? 0 : x[0]);
	     
	      ndiv++;
	    }
	  else 
	    {
     //perform the gillespie reaction
	      // pick the reaction according to propensity list and random number choice
	      tmp = r2*az0; sum = 0;
	      for(j = 0; sum < tmp && j < M; j++)
		sum += a[j];
	      reaction = j-1;

	      // update time and copy number record
	      t += tau;
	      for(tmpt = lastrec; tmpt < t; tmpt ++)
		{
		  rec0[NUMT*run + tmpt] = (x[0]);
		}
	      lastrec = t;

	      // update copy number according to this reaction
	      for(j = 0; j < r[reaction].nreacs; j++)
		x[r[reaction].reacs[j]] -= r[reaction].rstoc[j];

	      for(j = 0; j < r[reaction].nprods; j++)
		x[r[reaction].prods[j]] += r[reaction].pstoc[j];		

	    }

	}


    }

  /***** we've now completed the set of stochastic simulations *********/

  // open a file for output
  fp = fopen(fname, "w");
  // loop through the timespan of the simulation
  for(i = 0; i < NUMT; i++)
    {
      // compute mean and s.d. of copy number at this time
      mean0[i] = sd0[i] = 0; 
      for(j = 0; j < NUMR; j++)
	{
	  mean0[i] += rec0[NUMT*j+i];
	}
      mean0[i] /= NUMR; 
      for(j = 0; j < NUMR; j++)
	{
	  sd0[i] += pow(mean0[i]-rec0[NUMT*j+i], 2);
	}
      sd0[i] /= NUMR; 
      sd0[i] = sqrt(sd0[i]);
      // output these stats to file
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

  // use command-line argument to specify partitioning regime
  if(argc != 2)
    {
      printf("wrong number of arguments: need partitioning regime (0-4)\n");
      return 0;
    }
  P.partition = atoi(argv[1]);
  if(!(P.partition >= BINOMIAL && P.partition <= TELOMERE_ALT))
    P.partition = BINOMIAL;

  // initialise other parameters
  P.tau = 10;
  P.m0 = 100;
  P.nclus = 10; 
  P.eta = 100;
  P.nu = 0.01;
  P.lambda = P.nu + log(2.)/P.tau;

  if(P.partition == TELOMERE_ALT)
    {
      P.lambda = 0.02;
      P.m0 = 10;
      P.eta = 2;
    }

  // set up output filename and simulate
  sprintf(string, "partitioning-regime-%i.dat", P.partition);
  Simulate(string, P, 0);

  return 0;
}
