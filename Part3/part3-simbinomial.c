#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>

#include "gillespietmp.h"
#include "relaxed.h"

// simple wrapper for gillespie simulation of relaxed replication model with binomial partitioning
int main(void)
{
  Gillespie G;
  ChinneryParams C;
  double celltau;
  FILE *fp;
  double nend;

  celltau = 5;
  C.tau = 10; C.Noptimal = 1000; C.alpha = 5;

  InitialiseRelaxed(&G, C, 1000);
  G.numt = 40;
  G.nrun = 100000;

  fp = fopen("trythisnormal.txt", "w");
  nend = GillespieSimulateDivisions(G, 1, fp, celltau, 1);
  fclose(fp);

  return 0;
}
