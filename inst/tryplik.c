/*

 tryplik.c - Calculates likelihoods of tryptase copy numbers (alpha,
 active beta, beta frameshift) using tryptase read count information. Usage:

 TrypLik WT FS Areads Breads Dreads (-pop POP)

 WT = number of reads supporting active (non-frameshift) beta alleles
 FS = number of reads supporting beta3 frameshift alleles
 Areads = number of reads uniquely mapping to tryptase alpha
 Breads = number of reads uniquely mapping to tryptase beta
 Dreads = number of reads uniquely mapping to tryptase delta
 POP = Optional population identifier for genotype likelihood priors and
 frameshift probability FSfreq (default is EUR for European ancestry
 individuals)

 Note we round FS<=2 down to 0 (i.e., FS<=2 is called as 0 beta frameshift
 alleles).  Similarly, FS>2 implies at least one beta frameshift allele.
 Finally, we assume all individuals have at least 1 active beta allele.

 To compile, type

 cc -O -o TrypLik tryplik.c -lm

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Multiply likelihoods by arbitrary factor to reduce rounding issues */

#define FACTOR 1E250

/* Mean expected A/D and B/D values due to unknown biases */

#define expA 1.059
#define expB 1.027

void do_tryptase_calc (int WT, int FS, int Areads, int Breads, int Dreads, int POP, long double ****FinalLiks, long double * TotalLik)
{

  int Aeffective, Beffective, Deffective;
  int a, b, c;
  long double **GenoPriors, **GenoLiks, **FSLiks;
  double fact(int), FSfreq;

  if (POP==1) /* Set prior FS allele frequency among beta alleles */
    FSfreq = 0.0436;
  else if (POP==2)
    FSfreq = 0.0014;
  else if (POP==3)
    FSfreq = 0.0575;
  else
    FSfreq = 0.1349;

  Aeffective = (int) (Areads/2.75 + 0.5);
  Beffective = (int) (Breads/2.75 + 0.5);
  Deffective = (int) (Dreads/2.75 + 0.5);

  GenoPriors = (long double **) malloc (7 * sizeof (long double *));
  for (a=0; a<7; ++a) {
    GenoPriors[a] = (long double *) malloc (9 * sizeof (long double));
    for (b=1; b<9; ++b)
      GenoPriors[a][b] = 0.;
  }
  GenoLiks = (long double **) malloc (7 * sizeof (long double *));
  for (a=0; a<7; ++a) {
    GenoLiks[a] = (long double *) malloc (9 * sizeof (long double));
    for (b=1; b<9; ++b)
      GenoLiks[a][b] = 0.;
  }
  FSLiks = (long double **) malloc (9 * sizeof (long double *));
  for (a=0; a<9; ++a) {
    FSLiks[a] = (long double *) malloc (9 * sizeof (long double));
    for (b=0; b<9; ++b)
      FSLiks[a][b] = 0.;
  }

  if (POP==0) { /* Sets genotype priors for EUR ancestry populations */
    GenoPriors[0][2] = 0.0001;
    GenoPriors[0][3] = 0.013;
    GenoPriors[0][4] = 0.219;
    GenoPriors[0][5] = 0.020028;
    GenoPriors[0][6] = 0.0004;
    GenoPriors[0][7] = 0.00005;
    GenoPriors[0][8] = 0.00001;

    GenoPriors[1][1] = 0.0001;
    GenoPriors[1][2] = 0.013;
    GenoPriors[1][3] = 0.434;
    GenoPriors[1][4] = 0.03;
    GenoPriors[1][5] = 0.0004;
    GenoPriors[1][6] = 0.00001;
    GenoPriors[1][7] = 0.000001;

    GenoPriors[2][1] = 0.001;
    GenoPriors[2][2] = 0.208;
    GenoPriors[2][3] = 0.03;
    GenoPriors[2][4] = 0.004;
    GenoPriors[2][5] = 0.0001;
    GenoPriors[2][6] = 0.00004;

    GenoPriors[3][1] = 0.0001;
    GenoPriors[3][2] = 0.018;
    GenoPriors[3][3] = 0.005;
    GenoPriors[3][4] = 0.0004;
    GenoPriors[3][5] = 0.00001;

    GenoPriors[4][2] = 0.002;
    GenoPriors[4][3] = 0.0008;
    GenoPriors[4][4] = 0.00005;
    GenoPriors[4][5] = 0.00001;

    GenoPriors[5][2] = 0.0002;
    GenoPriors[5][3] = 0.0001;
    GenoPriors[5][4] = 0.00004;
    GenoPriors[5][5] = 0.00001;

    GenoPriors[6][2] = 0.00003;
    GenoPriors[6][3] = 0.00001;
    GenoPriors[6][4] = 0.000001;
  }
  else if (POP==1) { /* Sets genotype priors for AFR ancestry populations */
    GenoPriors[0][2] = 0.001;
    GenoPriors[0][3] = 0.015;
    GenoPriors[0][4] = 0.274;
    GenoPriors[0][5] = 0.1;
    GenoPriors[0][6] = 0.008;
    GenoPriors[0][7] = 0.0003;
    GenoPriors[0][8] = 0.0001;

    GenoPriors[1][1] = 0.0001;
    GenoPriors[1][2] = 0.012;
    GenoPriors[1][3] = 0.39;
    GenoPriors[1][4] = 0.06;
    GenoPriors[1][5] = 0.005;
    GenoPriors[1][6] = 0.0001;
    GenoPriors[1][7] = 0.00001;

    GenoPriors[2][1] = 0.0007;
    GenoPriors[2][2] = 0.126536;
    GenoPriors[2][3] = 0.003;
    GenoPriors[2][4] = 0.0001;
    GenoPriors[2][5] = 0.00001;
    GenoPriors[2][6] = 0.000001;

    GenoPriors[3][1] = 0.0001;
    GenoPriors[3][2] = 0.003;
    GenoPriors[3][3] = 0.0003;
    GenoPriors[3][4] = 0.0001;
    GenoPriors[3][5] = 0.0001;

    GenoPriors[4][2] = 0.0001;
    GenoPriors[4][3] = 0.0001;
    GenoPriors[4][4] = 0.0001;
    GenoPriors[4][5] = 0.0001;

    GenoPriors[5][2] = 0.00001;
    GenoPriors[5][3] = 0.00001;
    GenoPriors[5][4] = 0.00001;
    GenoPriors[5][5] = 0.00001;

    GenoPriors[6][2] = 0.000001;
    GenoPriors[6][3] = 0.000001;
    GenoPriors[6][4] = 0.000001;
  }
  else if (POP==2) { /* Sets genotype priors for EAS ancestry populations */
    GenoPriors[0][2] = 0.025;
    GenoPriors[0][3] = 0.037;
    GenoPriors[0][4] = 0.032;
    GenoPriors[0][5] = 0.004;
    GenoPriors[0][6] = 0.0002;
    GenoPriors[0][7] = 0.00002;
    GenoPriors[0][8] = 0.000002;

    GenoPriors[1][1] = 0.0021;
    GenoPriors[1][2] = 0.22;
    GenoPriors[1][3] = 0.23;
    GenoPriors[1][4] = 0.012;
    GenoPriors[1][5] = 0.0005;
    GenoPriors[1][6] = 0.00001;
    GenoPriors[1][7] = 0.000003;

    GenoPriors[2][1] = 0.004;
    GenoPriors[2][2] = 0.42;
    GenoPriors[2][3] = 0.005;
    GenoPriors[2][4] = 0.0005;
    GenoPriors[2][5] = 0.00001;
    GenoPriors[2][6] = 0.000002;

    GenoPriors[3][1] = 0.0001;
    GenoPriors[3][2] = 0.006;
    GenoPriors[3][3] = 0.001;
    GenoPriors[3][4] = 0.0001;
    GenoPriors[3][5] = 0.00001;

    GenoPriors[4][2] = 0.0001;
    GenoPriors[4][3] = 0.0001;
    GenoPriors[4][4] = 0.0001;
    GenoPriors[4][5] = 0.0001;

    GenoPriors[5][2] = 0.00001;
    GenoPriors[5][3] = 0.00001;
    GenoPriors[5][4] = 0.00001;
    GenoPriors[5][5] = 0.00001;

    GenoPriors[6][2] = 0.000001;
    GenoPriors[6][3] = 0.000001;
    GenoPriors[6][4] = 0.000001;
  }
  else if (POP==3) { /* Sets genotype priors for SAS ancestry populations */
    GenoPriors[0][2] = 0.004;
    GenoPriors[0][3] = 0.053;
    GenoPriors[0][4] = 0.165;
    GenoPriors[0][5] = 0.0325;
    GenoPriors[0][6] = 0.0006;
    GenoPriors[0][7] = 0.00005;
    GenoPriors[0][8] = 0.000005;

    GenoPriors[1][1] = 0.0005;
    GenoPriors[1][2] = 0.017;
    GenoPriors[1][3] = 0.444;
    GenoPriors[1][4] = 0.01;
    GenoPriors[1][5] = 0.0005;
    GenoPriors[1][6] = 0.00005;
    GenoPriors[1][7] = 0.000005;

    GenoPriors[2][1] = 0.0003;
    GenoPriors[2][2] = 0.256129;
    GenoPriors[2][3] = 0.004;
    GenoPriors[2][4] = 0.001;
    GenoPriors[2][5] = 0.0001;
    GenoPriors[2][6] = 0.00001;

    GenoPriors[3][1] = 0.0001;
    GenoPriors[3][2] = 0.01;
    GenoPriors[3][3] = 0.0003;
    GenoPriors[3][4] = 0.0001;
    GenoPriors[3][5] = 0.00001;

    GenoPriors[4][2] = 0.0003;
    GenoPriors[4][3] = 0.0001;
    GenoPriors[4][4] = 0.0001;
    GenoPriors[4][5] = 0.0001;

    GenoPriors[5][2] = 0.0001;
    GenoPriors[5][3] = 0.00001;
    GenoPriors[5][4] = 0.00001;
    GenoPriors[5][5] = 0.00001;

    GenoPriors[6][2] = 0.00001;
    GenoPriors[6][3] = 0.000001;
    GenoPriors[6][4] = 0.000001;
  }

  if ((40*Areads) < Dreads) /* Genotype likelihoods when #alpha = 0 */
    for (b=2; b<9; ++b)
      GenoLiks[0][b] = FACTOR * GenoPriors[0][b] *
        pow ((b*expB)/(b*expB+2.), (double) Beffective) *
        pow (2./(b*expB+2.), (double) Deffective);
  else if ((4*Areads) >= Dreads) /* Genotype likelihoods when A/D>=0.5 */
      for (a=1; a<7; ++a)
        for (b=1; b<9; ++b)
          GenoLiks[a][b] = FACTOR * GenoPriors[a][b] *
            pow ((a*expA)/(a*expA+b*expB+2.), (double) Aeffective) *
            pow ((b*expB)/(a*expA+b*expB+2.), (double) Beffective) *
            pow (2./(a*expA+b*expB+2.), (double) Deffective);

  FSLiks[1][0] = 1.;
  if (FS<=2) /* Frameshift likelihoods when #BetaFS=0 */
for (a=2; a<9; ++a)
  FSLiks[a][0] = 1.;
  else /* Frameshift likelihoods when #BetaFS>0 */
for (a=2; a<9; ++a)
  for (b=1; b<a; ++b)
    FSLiks[a][b] = (fact (a) / (fact (b) * fact (a-b))) *
      pow (FSfreq, (double) b) * pow ((1. - FSfreq), (double) a-b) *
      pow ((double) b/a, (double) FS) *
      pow ((double) (a-b)/a, (double) WT);

  for (a=0; a<7; ++a)
    for (b=1; b<9; ++b)
      for (c=0; c<b; ++c) {
        (*FinalLiks)[a][b][c] = GenoLiks[a][b] * FSLiks[b][c];
        (*TotalLik) += (*FinalLiks)[a][b][c];
      }
}

int main (int argc, char *argv[])
{
  int WT, FS, Areads, Breads, Dreads;
  int POP;
  int a, b, c;
  long double ***FinalLiks, TotalLik = 0.;

  if (argc == 1) {
    printf("Usage: TrypLik WT FS Areads Breads Dreads (-pop POP)\n");
    exit(1);
  }

  WT = atoi (argv[1]);
  FS = atoi (argv[2]);
  Areads = atoi (argv[3]);
  Breads = atoi (argv[4]);
  Dreads = atoi (argv[5]);
  POP = 0;
  if (argc>6 && strcmp(argv[6], "-pop")==0) {
    if (strcmp (argv[7], "AFR") == 0)
      POP = 1;
    else if (strcmp (argv[7], "EAS") == 0)
      POP = 2;
  }

  // note: memory leak since these aren't free()d
  FinalLiks = (long double ***) malloc (7 * sizeof (long double **));
  for (a=0; a<7; ++a) {
    FinalLiks[a] = (long double **) malloc (9 * sizeof (long double *));
    for (b=1; b<9; ++b) {
      FinalLiks[a][b] = (long double *) malloc (b * sizeof (long double));
      for (c=0; c<b; ++c)
        FinalLiks[a][b][c] = 0.;
    }
  }

  do_tryptase_calc(WT, FS, Areads, Breads, Dreads, POP, &FinalLiks, &TotalLik);

  printf("Alpha count\tBeta count\tBeta FS count\tPosterior likelihood\n");
  for (a=0; a<7; ++a)
    for (b=1; b<9; ++b)
      for (c=0; c<b; ++c)
	if (FinalLiks[a][b][c] / TotalLik > 0.000005)
	  printf("%d\t\t%d\t\t%d\t\t%.5Lf\n", a, b, c,
		 FinalLiks[a][b][c] / TotalLik);
}

double fact (int a)
{
  int c;
  double b=1.;

  if (a<2)
    return (1.);
  else
    for (c=2; c<=a; ++c)
      b *= (double) c;
  return (b);
}

