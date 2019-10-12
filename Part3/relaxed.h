// this requires gillespietmp.h

#define _MTDNA 0

// structure for relaxed replication parameterisation
typedef struct tag
{
  double alpha;
  double Noptimal;
  double tau;
} ChinneryParams;

// initialise simulation framework for relaxed replication simulation
void InitialiseRelaxed(Gillespie *G, ChinneryParams C, int n0)
{
  int n = 0;
  int i, j;
  double alpha, lambda, nu;

  // deal with subtleties of alpha parameter
  if(C.alpha <= 1)
    {
      alpha = C.alpha*C.Noptimal/C.tau;
      lambda = (1.-C.alpha)/C.tau;
      nu = 1./C.tau;
    }
  else
    {
      alpha = C.alpha*C.Noptimal/C.tau;
      lambda = 0;
      nu = 1./C.tau - (1.-C.alpha)/C.tau;
    }

  // introduce reactions
  SetupGillespie(n++, G, 0, 1, /**/            /* -> */ _MTDNA);
  SetupGillespie(n++, G, 1, 2, /**/ _MTDNA,    /* -> */ _MTDNA, _MTDNA);
  SetupGillespie(n++, G, 1, 0, /**/ _MTDNA     /* -> */       );

  G->N = 1; G->M = n;

  // set up stoichiometries and rates
  for(i = 0; i < G->M; i++)
    {
      for(j = 0; j < G->R[i].nreacs; j++)
	G->R[i].rstoc[j] = 1;
      for(j = 0; j < G->R[i].nprods; j++)
	G->R[i].pstoc[j] = 1;
      G->c[i] = 1;
    }

  G->c[0] = alpha; G->c[1] = lambda; G->c[2] = nu;

  // initialise
  for(i = 0; i < G->N; i++)
    G->init[i] = 0;

  G->init[_MTDNA] = n0;

  G->nrun = 10000;
  G->numt = 10;
}
