#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>

#include "gillespietmp.h"
#include "relaxed.h"

// simple wrapper for gillespie simulation of relaxed replication model with subtractive partitioning
int main(void)
{
  Gillespie G;
  ChinneryParams C;
  double celltau;
  FILE *fp;
  double postdiv;
  celltau = 5;
  C.tau = 10; C.Noptimal = 100; C.alpha = 5;

  InitialiseRelaxed(&G, C, 100);
  G.numt = 40;
  G.nrun = 500000;

  fp = fopen("trythisyeast.txt", "w");
  postdiv = GillespieSimulateSubtractions(G, 1, fp, celltau, 15, 0.1, 0);
  fclose(fp);

  printf("%f\n", postdiv);

  return 0;
}
