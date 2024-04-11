#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include "Rtypes.h"
#include "TROOT.h"
#include "TFile.h"
#include "TRandom.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TKey.h"
#include "TSystem.h"
#include "TCutG.h"
#include "TTree.h"
#include "gdecomp.h"

#include "GEBSort.h"
#include "GTMerge.h"
#define NGSGE NGE
#include "gsang.h"


#define TAPE_MOVED 10
#define EBIS_CLOCK 11
#define BETA_FIRED 12

#define ALL2DS 0
#define TRACE 0
#define SZ_EXTRA 0

/* Gain Calibration */

float ehigain[NGE + 1];
float ehioffset[NGE + 1];
float ehibase[NGE + 1];
void SetBeta ();

/* Other variables */

unsigned long long int EvTimeStam0 = 0;

/* pointers to ROOT spectra */

TH1D *hEventCounter;
TH2F *hGErate, *hBGOrate;
TH2F *hEhiRaw, *hEhiRawRaw, *hEhiCln, *hEhiDrty, *hEhiCln_nodop;
TH2F *hGeBGO_DT;
TH2F *pzraw;
TH2F *base1_diff;
TH2F *base2_diff;
TH1D *hrr;
TH1D *dgs_ds_singles;
TH2F *dgs_gg_1keV;
TH2F *dgs_gg;
TH2F *dgs_ggdt;
TH2F *dgs_gg_decay1;
TH2F *dgs_gg_decay1_beta;
TH2F *dgs_gg_decay2;
TH2F *dgs_gg_decay2_beta;
TH2F *dgs_gg_decay3;
TH2F *dgs_gg_decay3_beta;
TH2F *baseXid;
TH2F *dtevsp;
#if(SZ_EXTRA)
TH2F *dtevXbaseline[NGE + 1];
TH2F *spvdcXsum1[NGE + 1];
TH2F *spvdcXsum1_a[NGE + 1];
TH2F *spvdcXsum1_b[NGE + 1];
TH2F *spbaseXbase[NGE + 1];
TH2F *sampled_baselineXgid;
TH2F *sampled_baselinesum1Xgid;
TH2F *cor_factorXgid;
TH2F *eXdtev_det[NGE + 1];
TH2F *eXsum1[NGE + 1];
TH2F *sum1Xdtev[NGE + 1];
#endif
#if(TRACE)
TH2F *traces[NGE + 1];
#define NTRACES 200
int traceno[NTRACES + 1];
#endif
TH1D *rate_dgs_min;
TH2F *rate_dgs_sec;
TH2F *spfactor;
TH2F *eXdtev;
//TH1D *spcfdfrac1;
//TH1D *spcfdfrac2;
TH2F *sptt1;
TH2F *sptt2;
TH2F *dgs_sumtraceXge1;
TH2F *dgs_sumtraceXge2;

extern TH1D *hBetaCounter;
extern TH2F *hGeBeta_DT;
extern TH2F *hEhiBeta;
extern TH2F *hGeTape_DT;

#if(ALL2DS)
TH2F *e2_e1vse1[NGE + 1];
TH2F *SZe2_e1vse1[NGE + 1];
TH2F *hbase1, *hbase2;
#endif

long long int tpehits[100];

/* parameters */

extern DGSEVENT DGSEvent[MAXCOINEV];
extern int ng;
extern PARS Pars;
int tlkup[NCHANNELS];
int tid[NCHANNELS];

extern double angle[NGSGE], anglephi[NGSGE];

extern double DopCorFac[NGSGE], ACFac[NGSGE];

float factor[NGE], Base;

int en = 0;
float bl_array[5000];

float base1, base2, base1_av = 0;
float ave_base[NGE + 1];

int firstbl = 1;
int firstbl2[NGE + 1];
long long int dgsHeaderID[20];

float PZ[NGE + 1];

long long int ngood_e[NGE + 1];
long long int nbad_e[NGE + 1];

float hibaselim[NGE + 1], lobaselim[NGE + 1];
float base[NGE + 1] = { 0.0 };

#define NAVETRACE 1000

long long int nn_all[NGE + 1];
long long int nn_badsum1[NGE + 1];
long long int nn_badsum12[NGE + 1];

/*-----------------------------------------------------*/

int
exit_dgs ()
{
  /* declarations */

  int i, j, i1, i2, nn;
  float r1;
  long long int ill, nnt;
  double d1;
  FILE *fp;
  float meanrate, norm;
  TH2F *h2;

  printf ("\nbegin exit_dgs\n");

  /* header ID statistics */

  ill = 0;
  for (i = 0; i < 20; i++)
    ill += dgsHeaderID[i];
  for (i = 0; i < 20; i++)
    if (dgsHeaderID[i] > 0)
      printf ("DGS header ID %2i seen %12lli times, %6.2f%%\n", i, dgsHeaderID[i], 100.0 * dgsHeaderID[i] / ill);

  /* good/bad energies */

  printf ("\n");
  printf ("good/bad energy statistics\n");
  printf ("\n");
  for (i = 1; i <= NGE + 1; i++)
    if (ngood_e[i] > 0)
      {
        d1 = ngood_e[i] + nbad_e[i];
        d1 = ngood_e[i] / d1;
        d1 *= 100.0;
        printf ("ge %3i good/total fraction: %7.1f %%\n", i, d1);
      };
  printf ("\n");

  /* list avarage base lines we had at the end */

#if(0)
  fp = fopen ((char *) "last_baseline.txt", "w");
  printf ("last running baseline values:\n");
  for (i = 0; i < NGE; i++)
    if (ave_base[i] > 0)
      {
        printf ("ave_base[%3i]=%8.1f\n", i, ave_base[i]);
        fprintf (fp, "%i %f\n", i, ave_base[i]);
      };
  printf ("\n");
  fclose (fp);
#endif

  fp = fopen ((char *) "last_baseline.txt", "r");
  if (fp != NULL)
    {
      printf ("last running baseline values:\n");
      i2 = fscanf (fp, "%i %f", &i1, &r1);
      while (i2 == 2)
        {
          ave_base[i1] = r1;
          printf ("start ave_base[%3i]=%8.1f\n", i1, ave_base[i1]);
          i2 = fscanf (fp, "%i %f", &i1, &r1);
        };
      printf ("\n");
      fclose (fp);
    }
  else
    {
      for (i = 0; i < NGE; i++)
        ave_base[i] = 0;
      firstbl2[i] = 1;
    };

  /* tpe hits */

  printf ("\n");
  for (i = 0; i < 100; i++)
    if (tpehits[i] > 0)
      printf ("tpe= %i had %li hits\n", i, tpehits[i]);

  /* last baseline */

  printf ("\n");
  printf ("last baselines:\n");
  printf ("\n");
  for (i = 1; i <= NGE; i++)
    printf ("base[%3.3i]=%8.1f\n", i, base[i]);
  printf ("\n");

  /* rates in ge's */

  fp = fopen ("d2.cmd", "w");

  h2 = (TH2F *) gROOT->FindObject ("dtev");
  if (h2 == NULL)
    {
      printf ("could not find dtev\n");
    }
  else
    {
      printf ("found dtev\n");

      for (i = 1; i <= NGE; i++)
        {

          /* find mean from dtev spectrum */

          meanrate = 0;
          norm = 0;;
          for (j = 0; j < 1500; j++)
            {
              meanrate += j * h2->GetBinContent (i, j);
              norm += h2->GetBinContent (i, j);
            }
          if (norm > 0)
            meanrate /= norm;
          else
            meanrate = 0;

          if (meanrate > 0)
            {
              printf ("ge %3i: ", i);
              printf ("mean dtev %6.1f us, ", meanrate);
              printf ("mean rate: %5.1f kHz", 1.0 / meanrate * 1000.0);
              printf ("\n");
              fprintf (fp, "d2(\"vdcXsum1_%3.3i\")\n", i);
            };
        };
      fclose (fp);

    };

  /* bad sum statistics */

  printf ("\n");
  printf ("dgs statistics:\n");
  nnt = 0;
  nn = 0;
  for (i = 0; i <= NGE; i++)
    if (nn_all[i] > 0)
      {
        nnt += nn_all[i];
        nn++;
      }
  for (i = 0; i <= NGE; i++)
    if (nn_all[i] > 0)
      if (nn_badsum1[i] > 0 || nn_badsum12[i] > 0)
        {
          printf ("ge%3.3i ", i);
          printf (" fired %5.1f %% of all [rel %5.1f %%], ", 100.0 * nn_all[i] / nnt, nn * 100.0 * nn_all[i] / nnt);
          printf ("[sum1 <1] = %5.2f %% ", 100.0 * nn_badsum1[i] / nn_all[i]);
          printf ("[sum2<sum1] = %5.2f %%\n", 100.0 * nn_badsum12[i] / nn_all[i]);
        };
  printf ("we saw %3i dgs germanum detectors\n", nn);

  /* done */

  printf ("done exit_dgs\n");
  return (0);

};

/*-----------------------------------------------------*/

int
sup_dgs ()
{
  /* declarations */

  char str[STRLEN];
  int i, i1, i2, i7, i8;
  int imod, ichan;
  FILE *fp;
  float r1, r2;
  int maperrors=0;

  void getcal (char *);

// functions for making root histograms 

  TH2F *make2D (const char *, int, int, int, int, int, int);
  TH1D *make1D (const char *, int, int, int);

/* base diff spectra */

  /* -------------------- */
  /* read in the map file */
  /* -------------------- */

  for (i = 0; i < NCHANNELS; i++)
    {
      tlkup[i] = NOTHING;
      tid[i] = NOTHING;
    };

  fp = fopen ("map.dat", "r");
  if (fp == NULL)
    {
      printf ("need a \"map.dat\" file to run\n");
      exit (1);
    };

  maperrors = 0;
  printf ("\nmapping - started\n");

  /* map file format */
  /*  */
  /* nnnn  mm  kk  str */
  /*   |    |   |   +--- human readable; not used in program  */
  /*   |    |   +------- id number for that type */
  /*   |    +----------- type id, 1:GE 2:BGO 3:side 4:AUX etc */
  /*   +---------------- channel label (digId*10+ch) */

  i2 = fscanf (fp, "\n%i %i %i %s", &i1, &i7, &i8, str);
  while (i2 == 4)
    {
      /* trap duplicate entries */

      if (tlkup[i1] != NOTHING)
        {
          printf ("problem in mapfile: ");
          printf ("channel label %4i has already been assigned\n", i1);
          maperrors++;
          fflush (stdout);
        };

      tlkup[i1] = i7;
      tid[i1] = i8;
      i2 = fscanf (fp, "\n%i %i %i %s", &i1, &i7, &i8, str);
      ichan = i1 % 10;
      imod = (i1 - ichan) / 10;
      printf ("map: %4i ( %3i %1i ) %3i %3i %7s -- (official) ", i1, imod, ichan, i7, i8, str);
      switch (i7)
        {
        case GE:
          printf ("GE\n");
          break;
        case BGO:
          printf ("BGO\n");
          break;
        case SIDE:
          printf ("SIDE\n");
          break;
        case AUX:
          printf ("AUX\n");
          break;
        case DSSD:
          printf ("DSSD\n");
          break;
        case FP:
          printf ("FP\n");
          break;
        case XARRAY:
          printf ("XARRAY\n");
          break;
        case CHICO2:
          printf ("CHICO2\n");
          break;
        case SSD:
          printf ("SSD\n");
          break;
        case CLOVER:
          printf ("CLOVER\n");
          break;
        case SPARE:
          printf ("SPARE\n");
          break;
        case SIBOX:
          printf ("SIBOX\n");
          break;
        default:
          printf ("dont know what this is\n");
          break;

        };
    };
  fclose (fp);

  sprintf (str, "awk '{print $2, $3}' map.dat | uniq -d");
  system (str);

  if (maperrors > 0)
    {
      printf ("we had %i map.dat label errors\n", maperrors);
      printf ("Please fix map.dat and try again\n");
      fflush (stdout);
      exit (1);
    };

  printf ("\nmapping - complete!!\n");

#if(ALL2DS)
  base1_diff = make2D ("base1_diff", NGE, 1, NGE + 1, 2048, -1024, 1024);
  base2_diff = make2D ("base2_diff", NGE, 1, NGE + 1, 2048, -1024, 1024);
  hbase1 = make2D ("hbase1", NGE, 1, NGE + 1, 4096, 0, 16384);
  hbase2 = make2D ("hbase2", NGE, 1, NGE + 1, 4096, 0, 16384);
#endif


// 2-D's for Rate

  hEventCounter = make1D ("EvntCounter", 14400, 0, 14400);      // Good for 4 hours if Counts/sec
  hrr = make1D ("rr", 1012, 0, 2);
  hGErate = make2D ("GErate", 14400, 0, 14400, NGE, 1, NGE + 1);
  hBGOrate = make2D ("BGOrate", 14400, 0, 14400, NGE, 1, NGE + 1);
#if(ALL2DS)
  for (i = 1; i <= NGE; i++)
    {
      sprintf (str, "e2_e1vse1_%3.3i", i);
      e2_e1vse1[i] = make2D (str, 2048, 1, 10000, 1024, 1, 10000);
      sprintf (str, "SZe2_e1vse1_%3.3i", i);
      SZe2_e1vse1[i] = make2D (str, 2048, 1, 10000, 1024, 1, 10000);
    };
#endif

// 2-D's for Energy

  hEhiRaw = make2D ("EhiRaw", LENSP, 0, LENSP + 1, NGE, 1, NGE + 1);
  hEhiRawRaw = make2D ("EhiRawRaw", LENSP, 0, LENSP + 1, NGE, 1, NGE + 1);
  hEhiCln = make2D ("EhiCln", LENSP, 0, LENSP + 1, NGE, 1, NGE + 1);
  hEhiCln_nodop = make2D ("EhiCln_nodop", LENSP, 0, LENSP + 1, NGE, 1, NGE + 1);
  hEhiDrty = make2D ("EhiDrty", LENSP, 0, LENSP + 1, NGE, 1, NGE + 1);

// 2-D's for Tacs

  hGeBGO_DT = make2D ("GeBGO_DT", 400, -200, 200, NGE, 1, NGE + 1);

// 2-D's for PZ and Baseline

  pzraw = new TH2F ("pzraw", "pzraw spectra", NGE, 1, NGE + 1, 2000, 0, 2.0);

  dgs_ds_singles = make1D("DGS_DS_Singles",10000,0,10000);

  /* simple gg matix */

//  dgs_gg_1keV = make2D ("gg_1keV", 5000, 0, 5000, 5000, 0, 5000);
  dgs_gg_1keV = make2D ("gg_1keV", 10000, 0, 10000, 10000, 0, 10000);
  dgs_gg_1keV->SetXTitle ("g1");
  dgs_gg_1keV->SetYTitle ("g2");
//  dgs_gg = make2D ("gg", 2500, 0, 5000, 2500, 0, 5000);
  dgs_gg = make2D ("gg", 5000, 0, 10000, 5000, 0, 10000); 
  dgs_gg->SetXTitle ("g1");
  dgs_gg->SetYTitle ("g2");

  /* baseline */

  baseXid = new TH2F ("baseXid", "baseline vs geid", 4000, 1, 16000, NGE, 1, NGE);
  baseXid->SetXTitle ("base values from SZ algo 1 or 2");
  baseXid->SetYTitle ("ge id");

  if (Pars.have_decay_data)
    {
      /* tape station spectra */

      hBetaCounter = make1D ("BetaCounter", 14400, 0, 14400);
      hBetaCounter->SetXTitle ("sec since run start");
      hBetaCounter->SetYTitle ("counts");

      hGeBeta_DT = make2D ("GeBeta_DT", Pars.GGMAX, 1, Pars.GGMAX, 1000, -500, 500);
      hGeBeta_DT->SetXTitle ("ge energy");
      hGeBeta_DT->SetYTitle ("ge and beta time diff");

      hEhiBeta = make2D ("GeBeta_ID", Pars.GGMAX, 1, Pars.GGMAX, NGE, 1, NGE + 1);
      hEhiBeta->SetXTitle ("ge energy");
      hEhiBeta->SetYTitle ("ge id");

      hGeTape_DT = make2D ("GeTape_DT", Pars.GGMAX, 1, Pars.GGMAX, 1000, 0, 4000);
      hGeTape_DT->SetXTitle ("ge energy");
      hGeTape_DT->SetYTitle ("time since last tape move");

      /* simple gg decay matix, gates in chat file */

      dgs_gg_decay1 = new TH2F ("gg_decay1", "gg, in decay window 1", Pars.GGMAX, 1, Pars.GGMAX, Pars.GGMAX, 1, Pars.GGMAX);
      dgs_gg_decay1->SetXTitle ("g1");
      dgs_gg_decay1->SetYTitle ("g2");

      dgs_gg_decay2 = new TH2F ("gg_decay2", "gg, in decay window 2", Pars.GGMAX, 1, Pars.GGMAX, Pars.GGMAX, 1, Pars.GGMAX);
      dgs_gg_decay2->SetXTitle ("g1");
      dgs_gg_decay2->SetYTitle ("g2");

      dgs_gg_decay3 = new TH2F ("gg_decay3", "gg, in decay window 3", Pars.GGMAX, 1, Pars.GGMAX, Pars.GGMAX, 1, Pars.GGMAX);
      dgs_gg_decay3->SetXTitle ("g1");
      dgs_gg_decay3->SetYTitle ("g2");

      /* simple gg decay matix, gates in chat file, also require beta */

      dgs_gg_decay1_beta = new TH2F ("gg_decay1_beta", "gg with beta, in decay window 1", Pars.GGMAX, 1, Pars.GGMAX, Pars.GGMAX, 1, Pars.GGMAX);
      dgs_gg_decay1_beta->SetXTitle ("g1");
      dgs_gg_decay1_beta->SetYTitle ("g2");

      dgs_gg_decay2_beta = new TH2F ("gg_decay2_beta", "gg with beta, in decay window 2", Pars.GGMAX, 1, Pars.GGMAX, Pars.GGMAX, 1, Pars.GGMAX);
      dgs_gg_decay2_beta->SetXTitle ("g1");
      dgs_gg_decay2_beta->SetYTitle ("g2");

      dgs_gg_decay3_beta = new TH2F ("gg_decay3_beta", "gg with beta, in decay window 3", Pars.GGMAX, 1, Pars.GGMAX, Pars.GGMAX, 1, Pars.GGMAX);
      dgs_gg_decay3_beta->SetXTitle ("g1");
      dgs_gg_decay3_beta->SetYTitle ("g2");

    };

  dtevsp = new TH2F ("dtev", "time since last event", NGE, 1, NGE + 1, 10000, -10, 10000);
  dtevsp->SetXTitle ("gsid");
  dtevsp->SetYTitle ("time [us]");

#if(SZ_EXTRA)
  for (i = 1; i < NGE + 1; i++)
    {
      sprintf (str, "vdcXsum1_%3.3i", i);
      spvdcXsum1[i] = new TH2F ((char *) str, (char *) "SZ/fig6 factor plot", 1000, 0, 8000, 1000, 0, 6000);
      spvdcXsum1[i]->SetYTitle ("base[gsid]-baselast[gsid]");
      spvdcXsum1[i]->SetXTitle ("sum1-sum1last[gsid]");

      sprintf (str, "vdcXsum1_%3.3i_a", i);
      spvdcXsum1_a[i] = new TH2F ((char *) str, (char *) " SZ_base plot", 1000, 0, 8000, 1000, 0, 6000);
      spvdcXsum1_a[i]->SetYTitle ("SZ calculated base[gsid]");
      spvdcXsum1_a[i]->SetXTitle ("sum1");

      sprintf (str, "vdcXsum1_%3.3i_b", i);
      spvdcXsum1_b[i] = new TH2F ((char *) str, (char *) "sampled_baseline vs sum1 plot", 1000, 0, 8000, 1000, 0, 6000);
      spvdcXsum1_b[i]->SetYTitle ("sampled_baseline");
      spvdcXsum1_b[i]->SetXTitle ("sum1");

      sprintf (str, "eXsum1_%3.3i", i);
      eXsum1[i] = new TH2F ((char *) str, (char *) "energy vs sum1 plot", 5000, 0, 10000, 1500, 0, 1500);
      eXsum1[i]->SetYTitle ("SZ calculated energy");
      eXsum1[i]->SetXTitle ("sum1");

      sprintf (str, "baseXbase_%3.3i", i);
      spbaseXbase[i] = new TH2F ((char *) str, (char *) "SZ baseline vs sampled_baseline", 1000, 0, 5000, 1000, 0, 5000);
      spbaseXbase[i]->SetYTitle ("SZ baseline");
      spbaseXbase[i]->SetXTitle ("sampled_baseline");

      sprintf (str, "eXdtev_det_%3.3i", i);
      eXdtev_det[i] = new TH2F ((char *) str, (char *) "dtev vs corrected energy", 7000, 1, 7000, 500, 0, 500);
      eXdtev_det[i]->SetYTitle ("dtev");
      eXdtev_det[i]->SetXTitle ("Energy");

      sprintf (str, "sum1Xdtev_%3.3i", i);
      sum1Xdtev[i] = new TH2F ((char *) str, (char *) str, 1000, 0, 14000, 200, 1, 200);
      sum1Xdtev[i]->SetXTitle ("sum1");
      sum1Xdtev[i]->SetYTitle ("dtev");

    };

  sprintf (str, "sampled_baselineXgid");
  sampled_baselineXgid = new TH2F ((char *) str, (char *) "sampled_baseline vs gid", NGE, 1, NGE + 1, 1000, 0, 14000);
  sampled_baselineXgid->SetYTitle ("sampled_baseline");
  sampled_baselineXgid->SetXTitle ("gsid");

  sprintf (str, "sampled_baselinesum1Xgid");
  sampled_baselinesum1Xgid = new TH2F ((char *) str, (char *) "sampled_baseline - sum1 vs gid", NGE, 1, NGE + 1, 1000, -2000, 2000);
  sampled_baselinesum1Xgid->SetYTitle ("sampled_baseline-sum1");
  sampled_baselinesum1Xgid->SetXTitle ("gsid");

  sprintf (str, "cor_factorXgid");
  cor_factorXgid = new TH2F ((char *) str, (char *) "cor_factor vs gid", NGE, 1, NGE + 1, 1000, 0, 20);
  cor_factorXgid->SetYTitle ("energy corr factor");
  cor_factorXgid->SetXTitle ("gsid");

 sprintf (str, "dtevXbaseline_%3.3i", i);
 dtevXbaseline[i] = new TH2F ((char *) str, (char *) str, 5000, 0, 10000, 2000, 0, 500);
 dtevXbaseline[i]->SetYTitle ("dtev");
 dtevXbaseline[i]->SetXTitle ("calculated baseline");



#endif

  

#if(TRACE)
  for (i = 1; i < NGE + 1; i++)
    {
      sprintf (str, "traces_%3.3i", i);
      traces[i] = new TH2F ((char *) str, (char *) "traces", NTRACES, 1, NTRACES, 1001, 0, 1000);
      traceno[i] = 0;
    };
#endif

  rate_dgs_sec = new TH2F ((char *) "rate_dgs_sec", (char *) "rate of dgs", NGE, 1, NGE + 1, RATELEN, 0, RATELEN);
  rate_dgs_sec->SetXTitle ("sec");
  rate_dgs_sec->SetYTitle ("Hz");

  spfactor = new TH2F ((char *) "factor", (char *) "SZ2 extrapolation factor: DGS", NGE, 1, NGE + 1, 1000, -5, 5);
  spfactor->SetXTitle ("geid");
  spfactor->SetYTitle ("factor");

  eXdtev = new TH2F ((char *) "eXdtev", (char *) "Energies vs dtev", 1500, 1, 1500, 500, 0, 500);
  eXdtev->SetXTitle ("e gamma [channels]");
  eXdtev->SetYTitle ("dtev [us]");

//  spcfdfrac1 = new TH1D ((char *) "cfdfrac1", (char *) "additional CFD fraction [10 ns units] relatrive to nearest", 2000, -30, 30);
//  spcfdfrac2 = new TH1D ((char *) "cfdfrac2", (char *) "additional CFD fraction [10 ns units] relative to event TS", 2000, -30, 30);

  sptt1 = new TH2F ((char *) "tt1", (char *) "ge time vs ge time", 109 * 110 + 110, 0, 109 * 110 + 110, 100, -5, 5);
  sptt1->SetXTitle ("id1*110+id2, id1<id2");
  sptt1->SetYTitle ("10nsec units");

  sptt2 = new TH2F ((char *) "tt2", (char *) "ge time vs ref ge", 111, 0, 110, 50, -5, 5);
  sptt2->SetXTitle ("id2");
  sptt2->SetYTitle ("10nsec units");

/* list what we have */

  //printf (" we have define the following spectra:\n");

  Pars.wlist = gDirectory->GetList ();
//  Pars.wlist->Print ();


  /* Set Default Calibration and PZ parameters */

  for (i = 0; i <= NGE + 1; i++)
    {
      ehigain[i] = 1.0;
      ehioffset[i] = 0.0;
      PZ[i] = 1.0;
      factor[i] = 0.0;
      ehibase[i] = 0.0;
    };

  for (i = 0; i < 5000; i++)
    {
      bl_array[i] = 0;
    }

  for (i = 0; i < NGE + 1; i++)
    {
      ave_base[i] = 0;
    }


  /* get the DGS calibration file */

  getcal (Pars.dgs_ecalfn);

  /* dgs header ids */

  for (i = 0; i < 20; i++)
    dgsHeaderID[i] = 0;

#if(0)
  /* list enabled detectors */

  for (j = 0; j < MAXDETNO; j++)
    if (!Pars.enabled[j])
      printf ("bin_dgs: detector %3i is DISABLED\n", j);
  printf ("all others are enabled\n");
#endif

  /* good/bad energies */

  for (i = 0; i < NGE + 1; i++)
    {
      ngood_e[NGE + 1] = 0;
      nbad_e[NGE + 1] = 0;
    };

  fp = fopen ((char *) "last_baseline.txt", "r");
  if (fp != NULL)
    {
      printf ("last running baseline values:\n");
      i2 = fscanf (fp, "%i %f", &i1, &r1);
      while (i2 == 2)
        {
          ave_base[i1] = r1;
          printf ("start ave_base[%3i]=%8.1f\n", i1, ave_base[i1]);
          i2 = fscanf (fp, "%i %f", &i1, &r1);
        };
      printf ("\n");
      fclose (fp);
    }
  else
    {
      for (i = 0; i < NGE; i++)
        ave_base[i] = 0;
    };

  /* set default baseline limits */

  for (i = 0; i < NGE; i++)
    {
      hibaselim[i] = 4000;
      lobaselim[i] = 0;
    };

  /* get the base line limits to use */

  printf ("attempting to open \"baseline_limits.txt\"\n");
  fp = fopen ((char *) "baseline_limits.txt", "r");
  if (fp != NULL)
    {

      /* first set all values to default */

      for (i = 0; i < NGE; i++)
        {
          hibaselim[i] = 4000;
          lobaselim[i] = 0;
        };

      /* read any values in file */

      i2 = fscanf (fp, "%i %f %f", &i1, &r1, &r2);
      while (i2 == 3)
        {
          lobaselim[i1] = r1;
          hibaselim[i1] = r2;
          i2 = fscanf (fp, "%i %f %f", &i1, &r1, &r2);
        };
      printf ("base limits read\n");
      fclose (fp);

    }
  else
    {
      printf ("\"baseline_limits.txt\" not found, using defaults\n");

    };

  /* list base line limits */

  for (i = 0; i < NGE; i++)
    printf ("ge%3.3i: baseline limits from %5.0f to %5.0f\n", i, lobaselim[i], hibaselim[i]);

  /* init */

  for (i = 0; i < 100; i++)
    tpehits[i] = 0;

  printf ("Pars.dgs_MM= %f\n", Pars.dgs_MM);

  if (Pars.dgs_algo == 2)
    {
      assert (Pars.dgs_SZ_t2 > 0);
      assert (Pars.dgs_SZ_t1 > 0);
      assert (Pars.dgs_SZ_t1 > Pars.dgs_SZ_t2);
    };

  /* trace matrix for Pat and Co */

  dgs_sumtraceXge1 = new TH2F ("dgs_sumtraceXge1", "summed-traces vs geid", NGE, 1, NGE + 1, 4000, 0, 16000);
  dgs_sumtraceXge2 = new TH2F ("dgs_sumtraceXge2", "summed-traces-base vs geid", NGE, 1, NGE + 1, 4000, 0, 16000);

  for (i = 0; i <= NGE; i++)
    {
      nn_all[i] = 0;
      nn_badsum1[i] = 0;
      nn_badsum12[i] = 0;
    };

  /* done */

  printf ("\nsup_dgs done!!\n");

  return (0);

};


/* ----------------------------------------------------------------- */

int
bin_dgs (GEB_EVENT *GEB_event)
{

  /* declarations */

  char str[128];
  int i, j, gsid, e, st;
  double d1, d2, erawraw[NGE + 1];
  float pz1, pz2, pz3, pz4, sum1, sum2;
  signed short int ssi;
  double x1, x2, y1, y2, bb, aa, cc;
  float etmp[20];
  double ttmp[20];
  int idtmp[20];
  int jj, indx;
  FILE *fp;
  static int nwritetraces = 0;
  static int nsumtraces[NGE + 1] = { 0 };
  float sp[16000];

  double MicroSECOND = 100;
  float msample = 1024;         /* 10.24us = 1024*10ns */
  static float baselast[NGE] = { 0.0 }, sum1last[NGE] = {0.0};
  static long long int DTlast[NGE] = { 0 };

  int RelEvT = 0, tdiff = 0, tdiff1 = 0, tdiff2 = 0;
  float Energy, Energy_nodop, top, bot, r1, rr, ggdt;
  long long int ggdtll;
  unsigned long long int ulli, now, then;

  int is_beta = 0;
  static long long int tape_ts, ebis_ts, beta_ts;

  /* prototypes */

  int GebTypeStr (int type, char str[]);
  int DGSEvDecompose_v3 (unsigned int *ev, int len, DGSEVENT * DGSEvent, int tlkup[], int tid[]);
  int wr_spe (char *, int *, float *sp);

  /* Print debug */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("\nentered bin_dgs:\n");

  /* loop through the coincidence event and fish out DGS data */
  /* (gamma rays) count in ng */
  /* It will also handle type 15 data, data from the trigger */
  /* passed on the the digitizer */

  ng = 0;
  for (i = 0; i < GEB_event->mult; i++)
    {
      if (GEB_event->ptgd[i]->type == GEB_TYPE_DGS)
        {
          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              GebTypeStr (GEB_event->ptgd[i]->type, str);
              printf ("bin_dgs (header): %2i> %2i, %s, TS=%lli\n", i, GEB_event->ptgd[i]->type, str, GEB_event->ptgd[i]->timestamp);
            }

          st =
            DGSEvDecompose_v3 ((unsigned int *) GEB_event->ptinp[i], GEB_event->ptgd[i]->length / sizeof (unsigned int), &DGSEvent[ng], tlkup, tid);
          if (st != 0)
            return (0);

          ng++;
        }
    }

  /* process the three CFD values (if there) */
  /* They are 14 bit floats */

  for (i = 0; i < ng; i++)
    {
      DGSEvent[i].event_cfd_timestamp = DGSEvent[i].event_timestamp;
      if (DGSEvent->cfd_valid_flag)
        {

          /* extract the cfd values */

#if(1)
          /* TL's way */

          ssi = (signed short int) (DGSEvent[i].cfd_sample_0 << 2);
          DGSEvent[i].cfd_0 = (float) ssi / 4.0;
          ssi = (signed short int) (DGSEvent[i].cfd_sample_1 << 2);
          DGSEvent[i].cfd_1 = (float) ssi / 4.0;
          ssi = (signed short int) (DGSEvent[i].cfd_sample_2 << 2);
          DGSEvent[i].cfd_2 = (float) ssi / 4.0;
          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              printf ("DGSEvent[%i].cfd_sample_0|1|2 = %10i ", i, DGSEvent[i].cfd_sample_0);
              printf ("%10i ", DGSEvent[i].cfd_sample_1);
              printf ("%10i (unsigned)\n", DGSEvent[i].cfd_sample_2);
              printf ("[ tl] find zero crossing from: (in 10ns units)\n");
              printf ("  %lli %10.2f\n", DGSEvent[i].event_timestamp, DGSEvent[i].cfd_0);
              printf ("  %lli %10.2f\n", DGSEvent[i].event_timestamp - 1, DGSEvent[i].cfd_1);
              printf ("  %lli %10.2f\n", DGSEvent[i].event_timestamp - 2, DGSEvent[i].cfd_2);
            };
#endif
#if(0)
          /* JTA's way */

          ssi = (signed short int) DGSEvent[i].cfd_sample_0;
          if (DGSEvent[i].cfd_sample_0 > 8191)
            ssi = ssi - 16384;
          DGSEvent[i].cfd_0 = (float) ssi;
          ssi = (signed short int) DGSEvent[i].cfd_sample_1;
          if (DGSEvent[i].cfd_sample_1 > 8191)
            ssi = ssi - 16384;
          DGSEvent[i].cfd_1 = (float) ssi;
          ssi = (signed short int) DGSEvent[i].cfd_sample_2;
          if (DGSEvent[i].cfd_sample_2 > 8191)
            ssi = ssi - 16384;
          DGSEvent[i].cfd_2 = (float) ssi;
          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              printf ("DGSEvent[%i].cfd_sample_0|1|2 = %10i ", i, DGSEvent[i].cfd_sample_0);
              printf ("%10i ", DGSEvent[i].cfd_sample_1);
              printf ("%10i (unsigned)\n", DGSEvent[i].cfd_sample_2);
              printf ("[jta] find zero crossing from: (in 10ns units)\n");
              printf ("  %lli %10.2f\n", DGSEvent[i].event_timestamp, DGSEvent[i].cfd_0);
              printf ("  %lli %10.2f\n", DGSEvent[i].event_timestamp - 1, DGSEvent[i].cfd_1);
              printf ("  %lli %10.2f\n", DGSEvent[i].event_timestamp - 2, DGSEvent[i].cfd_2);
            };
#endif

          /* find zero cross */

          if (DGSEvent[i].cfd_0 < 0 && DGSEvent[i].cfd_1 > 0)
            {
              y2 = DGSEvent[i].cfd_1;
              y1 = DGSEvent[i].cfd_0;
              x1 = DGSEvent[i].event_timestamp;
              x2 = DGSEvent[i].event_timestamp - 1;
            }
          else if (DGSEvent[i].cfd_1 < 1 && DGSEvent[i].cfd_2 > 0)
            {
              y2 = DGSEvent[i].cfd_1;
              y1 = DGSEvent[i].cfd_2;
              x1 = DGSEvent[i].event_timestamp - 1.0;
              x2 = DGSEvent[i].event_timestamp - 2.0;
            }
          else
            {
              /* use the first two anyway, good move? */
              y2 = DGSEvent[i].cfd_1;
              y1 = DGSEvent[i].cfd_0;
              x1 = DGSEvent[i].event_timestamp;
              x2 = DGSEvent[i].event_timestamp - 1;
            }
          bb = (y2 - y1) / (x2 - x1);
          aa = y2 - x2 * bb;
          if (bb != 0)
            DGSEvent[i].event_cfd_timestamp = -aa / bb;

          /* bin the correction, note: x1 not .event_timestamp */

//          spcfdfrac1->Fill (DGSEvent[i].event_cfd_timestamp - x1);
//          spcfdfrac2->Fill (DGSEvent[i].event_cfd_timestamp - DGSEvent[i].event_timestamp);

          /* detailed debug list */

          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
//              printf ("DGSEvent[%i].cfd_sample_0|1|2 = %10i ", i, DGSEvent[i].cfd_sample_0);
//              printf ("%10i ", DGSEvent[i].cfd_sample_1);
//              printf ("%10i (unsigned)\n", DGSEvent[i].cfd_sample_2);
//              printf ("DGSEvent[%i].cfd_0|2|3 =        %10.2f", i, DGSEvent[i].cfd_0);
//              printf (" %10.2f ", DGSEvent[i].cfd_1);
//              printf (" %10.2f (floats)\n", DGSEvent[i].cfd_2);
//              printf ("find zero crossing from: (in 10ns units)\n");
//              printf ("  %lli %10.2f\n", DGSEvent[i].event_timestamp, DGSEvent[i].cfd_0);
//              printf ("  %lli %10.2f\n", DGSEvent[i].event_timestamp - 1, DGSEvent[i].cfd_1);
//              printf ("  %lli %10.2f\n", DGSEvent[i].event_timestamp - 2, DGSEvent[i].cfd_2);
              printf ("  %10.2f (cross zero) <-- is CFD time in 10 ns units\n", DGSEvent[i].event_cfd_timestamp);
            };

        };
    };


  /* debug list */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("\n\nCurEvNo: %i has %i DGS events (ng)\n", Pars.CurEvNo, ng);
      for (i = 0; i < ng; i++)
        {
          printf ("DGSEvent[%i].event_timestamp: %llu ", i, DGSEvent[i].event_timestamp);
          printf ("tpe= %i tid=%i ", DGSEvent[i].tpe, DGSEvent[i].tid);
          printf ("DGS header type: %i ", DGSEvent[i].header_type);
          if (DGSEvent[i].header_type < 15)
            printf (" (regular data)");
          else
            printf (" (trigger data)");
          printf ("\n");
        };
    };

  /* rate spectrum */

  if (EvTimeStam0 == 0)
    EvTimeStam0 = DGSEvent[0].event_timestamp;
  RelEvT = (int) ((DGSEvent[0].event_timestamp - EvTimeStam0) / 100000000);     // overflow?
  hEventCounter->Fill (RelEvT);

  /* Loop */

  for (i = 0; i < ng; i++)
    {
      Energy = 0;
      gsid = DGSEvent[i].tid;
      tpehits[DGSEvent[i].tpe]++;
      if (DGSEvent[i].tpe == GE)
        if (Pars.enabled[gsid])
          {

            /* this should be unnecessary to do; */
            /* but we sometimes crash if we don't. */
            /* Needs to be looked at */

            if (gsid < 1 || gsid > NGE)
              {
                printf ("bad gsid= %i\n", gsid);
                fflush (stdout);
                gsid = 0;
              };

            /* fill the summed trace matix */

            if (Pars.CurEvNo <= Pars.NumToPrint && 0)
              {
                printf ("bin_dgs: DGSEvent[%i].traceLen=%i\n", i, DGSEvent[i].traceLen);
                for (j = 0; j < DGSEvent[i].traceLen; j++)
                  printf ("bin_dgs: DGSEvent[%i].trace[%4i]=%i\n", i, j, DGSEvent[i].trace[j]);
              };
            if (nsumtraces[gsid] < NAVETRACE)
              {
                nsumtraces[gsid]++;
                for (j = 0; j < DGSEvent[i].traceLen; j++)
                  {
                    dgs_sumtraceXge1->Fill (gsid, j, (float) DGSEvent[i].trace[j] / NAVETRACE);
//              dgs_sumtraceXge2->Fill(gsid,(float)DGSEvent[i].trace[j]-(float)DGSEvent[i].sampled_baseline);
                    dgs_sumtraceXge2->Fill (gsid, j, (float) DGSEvent[i].trace[j] - base[gsid]);
                  };
              };

#if(TRACE)
            if (Pars.CurEvNo > 1000)
              if (traceno[gsid] < NTRACES)
                {
                  for (j = 0; j < DGSEvent[i].traceLen; j++)
                    traces[gsid]->SetBinContent (traceno[gsid], j, DGSEvent[i].trace[j]);
                  traceno[gsid]++;
                };
#endif


            hGErate->Fill ((int) ((DGSEvent[0].event_timestamp - EvTimeStam0) / 100000000), gsid);

            /* time since last event */
            /* easy for LED timing, little tricky for CFD timing */

            if (DGSEvent->header_type % 2 == 1)
              {
                if (Pars.CurEvNo <= Pars.NumToPrint)
                  {
                    printf ("__parameters from jta.c:\n");
                    printf ("__*LED timing;; has 47 bits 0x%16.16jx\n", (unsigned long long int) 1 << 47);
                    printf ("__DGSEvent[%i].pileup_count = %i\n", i, DGSEvent[i].pileup_count);
                    printf ("__DDGSEvent[%i].event_timestamp    = %lli\n", i, DGSEvent[i].event_timestamp);
                    printf ("__DGSEvent[%i].last_disc_timestamp = %lli\n", i, DGSEvent[i].last_disc_timestamp);
                    printf ("__dtev = %lli ", DGSEvent[i].event_timestamp - DGSEvent[i].last_disc_timestamp);
                    printf ("or %f usec\n", (DGSEvent[i].event_timestamp - DGSEvent[i].last_disc_timestamp) / 100.);
                    printf ("__DGSEvent[%i].p2_mode = %lli; ", i, DGSEvent[i].p2_mode);
                    printf ("__GSEvent[%i].capture_parst_ts = %lli\n", i, DGSEvent[i].capture_parst_ts);
                    printf ("__DGSEvent[%i].event_timestamp = %lli\n", i, DGSEvent[i].event_timestamp);
                    printf ("__DGSEvent[%i].sampled_baseline =%i\n", i, DGSEvent[i].sampled_baseline);
                    printf ("__DGSEvent[%i].sum1 | sum2 = %i | %i; ", i, DGSEvent[i].sum1, DGSEvent[i].sum2);
                    printf ("sum2-sum1 = %i\n", DGSEvent[i].sum2 - DGSEvent[i].sum1);
                    printf ("__DGSEvent[%i].peak_timestamp = %lli\n", i, DGSEvent[i].peak_timestamp);
                    printf ("__DGSEvent[%i].p_sum = %i\n", i, DGSEvent[i].p_sum);
                    printf ("__DGSEvent[%i].multiplex_data_field= %i\n", i, DGSEvent[i].multiplex_data_field);
                    printf ("__DGSEvent[%i].early_pre_rise= %i\n", i, DGSEvent[i].early_pre_rise);
                    printf ("__DGSEvent[%i].course_timestamp= %li\n", i, DGSEvent[i].course_timestamp);
                  };
                DGSEvent[i].dtev = ((double) (DGSEvent[i].event_timestamp) - (double) DGSEvent[i].last_disc_timestamp);
              }
            else
              {
                if (Pars.CurEvNo <= Pars.NumToPrint)
                  printf ("*CFD timing;; has 29 bits 0x%16.16jx\n", (unsigned long long int) 1 << 29);
                now = DGSEvent[i].event_timestamp & 0x3fffffff;
                then = DGSEvent[i].last_disc_timestamp & 0x3fffffff;
                DGSEvent[i].dtev = (double) now - (double) then;
//                DGSEvent[i].dtev += 78; /* <---- special TEST */
                if (DGSEvent[i].dtev <= 0)
                  DGSEvent[i].dtev += 0x40000000;
                if (Pars.CurEvNo <= Pars.NumToPrint)
                  {
                    printf ("tnow = %25lli 0x%16.16x %f sec\n", now, now, (double) now / 100000000.0);
                    printf ("then = %25lli 0x%16.16x %f sec\n", then, then, (double) then / 100000000.0);
                  };
              };
            DGSEvent[i].dtev /= MicroSECOND;
            if (Pars.CurEvNo <= Pars.NumToPrint)
              {
                printf ("DGSEvent[i].dtev= %f us ", DGSEvent[i].dtev);
                printf (" %f ms ", DGSEvent[i].dtev / 1000.0);
                printf (" %f s\n", DGSEvent[i].dtev / 1000000.0);
              };

            /* bin the DGSEvent[i].dtev. For the first few events */
            /* per det, you will see large DGSEvent[i].dtev because */
            /* previous event is not set yet */

            dtevsp->Fill (gsid, DGSEvent[i].dtev);

            /*-------------------------------------------------*/
            /* select the code to extract the gamma ray energy */
            /*-------------------------------------------------*/

            /* SZ_1.c and SZ_2.c are the new algos from SZ in 2021 */
            /* based on his paper */

            /* SZ_0_3456.c is the method implemented with m1begin...m2end */
            /* used for type 3/4 and 5/6 header types */

            /* make sure sum1 and sum2 have reasonable values */
            /* or we can get spikes and poluted last baselines */

            nn_all[gsid]++;
            if ((DGSEvent[i].sum1 > 0) && (DGSEvent[i].sum1 < DGSEvent[i].sum2))
              {

                switch (Pars.dgs_algo)
                  {
                  case 0:
#include "SZ_0_3456.stub"
                    break;
                  case 1:
#include "SZ_1.stub"
                    break;
                  case 2:
#include "SZ_2.stub"
                    break;
                  default:
                    printf ("dgs_algo %i not available\n", Pars.dgs_algo);
                    exit (1);
                    break;
                  }

                /* gain match */

                erawraw[i] = Energy;
                Energy = Energy * ehigain[gsid] + ehioffset[gsid];
                if (Pars.CurEvNo <= Pars.NumToPrint)
                  printf ("Energy = %f (calibrated), raw = %f\n", Energy, erawraw[i]);
              }
            else
              {
                if (DGSEvent[i].sum1 < 1)
                  nn_badsum1[gsid]++;
                if (DGSEvent[i].sum1 > DGSEvent[i].sum2)
                  nn_badsum12[gsid]++;
                Energy = 0;
                erawraw[i] = Energy;
              };

            /* count good and bad energies */

            if (Energy > 10.0 && Energy < 4000.0)
              {
                ngood_e[gsid]++;
                if (Pars.CurEvNo <= Pars.NumToPrint)
                  printf ("GOOD energy\n");
              }
            else
              {
                nbad_e[gsid]++;
                if (Pars.CurEvNo <= Pars.NumToPrint)
                  printf ("BAD energy\n");
              }

            /* Doppler correct the energy */

            Energy_nodop = Energy;
            if (Pars.beta != 0)
              {
                d1 = angtheta[gsid - 1] / 57.29577951;
                Energy = Energy * (1 - Pars.beta * cos (d1)) / sqrt (1 - Pars.beta * Pars.beta);
              }
            if (Pars.enabled[gsid])
              {
                DGSEvent[i].ehi = Energy;
                DGSEvent[i].ehi_nodop = Energy_nodop;
              }
            else
              {
                /* mark bad */
                DGSEvent[i].ehi = -1;
                DGSEvent[i].ehi_nodop = -1;
              }
            DGSEvent[i].id = gsid;

            /* fill the eXsum1 for PZ determination */

#if(SZ_EXTRA)
            eXsum1[DGSEvent[i].id]->Fill (sum1, DGSEvent[i].ehi);
            eXdtev_det[DGSEvent[i].id]->Fill (DGSEvent[i].ehi, DGSEvent[i].dtev);
#endif

#if(ALL2DS)
            if (Pars.enabled[gsid])
              {
                d1 = (DGSEvent[i].sum2 - DGSEvent[i].sum1) / Pars.dgs_MM;       /* uncorrected energy */
                d2 = DGSEvent[i].sum1 / Pars.dgs_MM;    /* ~ baseline */
                if (gsid <= NGE)
                  if (d1 < (double) 8192)
                    if (d2 < (double) 8192)
                      {
                        e2_e1vse1[gsid]->Fill (d1, d2, 1);
                        d1 = DGSEvent[i].ehi;   /* SZ corrected */
                        SZe2_e1vse1[gsid]->Fill (d1, d2, 1);
                      };
              };
#endif




            /* do the Compton suppression by setting a flag */
            /* loop through event for BGO. Note loop uses 'j' */

            for (j = 0; j < ng; j++)
              {
                if (DGSEvent[j].tpe == BGO && DGSEvent[j].tid == gsid)
                  {             // BGO & GE in coincidence
                    tdiff = (int) (DGSEvent[i].event_timestamp - DGSEvent[j].event_timestamp);
                    hGeBGO_DT->Fill (tdiff, gsid);
                    if (abs (tdiff) <= 50)
                      {
                        DGSEvent[i].flag = 1;   // Mark as Dirty Ge
//                    printf("BGO %i vetoed\n",DGSEvent[j].tid);
                      };
                  }
              }
          }


      if (DGSEvent[i].tpe == BGO)
        {
          if (Pars.enabled[gsid])
            {
              hBGOrate->Fill ((int) ((DGSEvent[0].event_timestamp - EvTimeStam0) / 100000000), DGSEvent[i].tid);
              DGSEvent[i].ehi = (float) (DGSEvent[i].sum2) - (float) (DGSEvent[i].sum1);
            }
          else
            {
              DGSEvent[i].ehi = -1;
            }
        };

      if (Pars.have_decay_data)
        {

          /* look for possible tape station signals */

          /* Note: normally the tape station data is piped */
          /* into the data stream through the DGS digitizers. */
          /* Thus, the only way to can detect tape station */
          /* data is through the map file */
          /* hunting for GEB_TYPE_XA headers does not work */

          /* here is an example of tape station entries in map file */
          /*   1010   10  111 TapeT0 */
          /*   1011   11  111 EBIS */
          /*   1018   12  111 BETA */

          /* Tape Moved */

          if (DGSEvent[i].tpe == TAPE_MOVED)
            {
              printf ("Tape Moved %llu %llu \n", DGSEvent[i].event_timestamp, tape_ts);
              tape_ts = DGSEvent[i].event_timestamp;
            };

          /* EBIS clock */

          if (DGSEvent[i].tpe == EBIS_CLOCK)
            {
              ebis_ts = DGSEvent[i].event_timestamp;
              //printf("EBIS Clock %i \n",DGSEvent[i].tpe);
            };


          /* Beta fired */

          if (DGSEvent[i].tpe == BETA_FIRED)
            {
              is_beta = 1;
              beta_ts = DGSEvent[i].event_timestamp;
              d1 = (double) (beta_ts - EvTimeStam0) / 100000000;
              if (d1 < 14400)
                hBetaCounter->Fill (d1);
            };
        };

    }                           /* for (i = 0; i < ng; i++) */






  /*-------------------------------*/
  /* Finally done with the primary */
  /* for (i = 0; i < ng; i++) loop */
  /* now start binning things */
  /*-------------------------------*/

  /* Energy Histogram loop */

  for (i = 0; i < ng; i++)
    {
      if (DGSEvent[i].tpe == GE)
        {
          if (erawraw[i] > 0 && erawraw[i] < LENSP)
            hEhiRawRaw->Fill (erawraw[i], gsid);
          e = (int) DGSEvent[i].ehi;
          if (e > 0 && e < LENSP)
            {
              gsid = DGSEvent[i].tid;
              hEhiRaw->Fill (e, gsid);
              if (DGSEvent[i].flag == 0)
                {
                  hEhiCln->Fill (e, gsid);
                  hEhiCln_nodop->Fill (DGSEvent[i].ehi_nodop, gsid);

                  /* clean energy vs time since last event */

                  eXdtev->Fill (DGSEvent[i].ehi, DGSEvent[i].dtev);

                  /* beta gated matrix */

                  if (Pars.have_decay_data)
                    if (is_beta)
                      {
                        /* beta-gated e_gamma vs time difference between Ge & beta */

                        tdiff = (int) (DGSEvent[i].event_timestamp - beta_ts);
                        hGeBeta_DT->Fill (e, tdiff);

                        /* beta-gated e_gamma vs id */

                        if (tdiff > Pars.gb_dt_lo && tdiff < Pars.gb_dt_hi)
                          hEhiBeta->Fill (e, gsid);

                        /* Generate tape histogram for beta-gated e_gamma */

                        if (tape_ts != 0)
                          {
                            tdiff = ((DGSEvent[i].event_timestamp - tape_ts) / 10000000);
                            hGeTape_DT->Fill (e, tdiff);
                          };

                      };

                };
              if (DGSEvent[i].flag == 1)
                hEhiDrty->Fill (e, gsid);
            };
        };
    }


  // Downscaled singles
  if (ng == 1){
    if (Pars.enabled[DGSEvent[0].tid])
      if (DGSEvent[0].tpe == GE)
        if (DGSEvent[0].flag == 0)
          if (DGSEvent[0].ehi > 10)
            dgs_ds_singles->Fill((double)DGSEvent[i].ehi); 
  // 1.7841523 is hardcoded correction for change in DGS gain range
  }
  /* gg matrices */

  if (ng >= 2)
    for (i = 0; i < ng; i++)
      if (Pars.enabled[DGSEvent[i].tid])
        if (DGSEvent[i].tpe == GE)
          if (DGSEvent[i].flag == 0)
//            if (DGSEvent[i].ehi > 0 && DGSEvent[i].ehi < 5000)
            if (DGSEvent[i].ehi > 0 && DGSEvent[i].ehi < 10000)
              for (j = i + 1; j < ng; j++)
                if (Pars.enabled[DGSEvent[j].tid])
                  if (DGSEvent[j].tpe == GE)
                    if (DGSEvent[j].flag == 0)
//                      if (DGSEvent[j].ehi > 0 && DGSEvent[j].ehi < 5000)
                      if (DGSEvent[j].ehi > 0 && DGSEvent[j].ehi < 10000)
                        {
                          /* regular gg matrix */

			  // 1.7841523 is hardcoded correction for change in DGS gain range
                          dgs_gg_1keV->Fill ((double) DGSEvent[i].ehi, (double) DGSEvent[j].ehi, 1);
                          dgs_gg_1keV->Fill ((double) DGSEvent[j].ehi, (double) DGSEvent[i].ehi, 1);
 
			  // 1.7841523 is hardcoded correction for change in DGS gain range                        
    			           dgs_gg->Fill ((double) DGSEvent[i].ehi, (double) DGSEvent[j].ehi, 1);
                          dgs_gg->Fill ((double) DGSEvent[j].ehi, (double) DGSEvent[i].ehi, 1);

                          if (Pars.CurEvNo <= Pars.NumToPrint)
                            printf ("filled dgs_gg with %f %f\n", DGSEvent[i].ehi, DGSEvent[j].ehi);

                          if (Pars.have_decay_data)
                            {

                              /* ge energy vs ge-ge-dt */

                              ggdtll = DGSEvent[i].event_timestamp - DGSEvent[j].event_timestamp;
                              ggdt = fabs ((double) ggdtll);
                              dgs_ggdt->Fill ((double) DGSEvent[i].ehi, (double) ggdt);

                              if (ggdt < Pars.ggdt)
                                {

                                  /* decay gamma-gamma matrix */

                                  tdiff1 = ((DGSEvent[i].event_timestamp - tape_ts) / 10000000);
                                  tdiff2 = ((DGSEvent[j].event_timestamp - tape_ts) / 10000000);

                                  if (tdiff1 > Pars.decay_lo1 && tdiff1 < Pars.decay_hi1)
                                    if (tdiff2 > Pars.decay_lo1 && tdiff2 < Pars.decay_hi1)
                                      {
                                        dgs_gg_decay1->Fill ((double) DGSEvent[i].ehi, (double) DGSEvent[j].ehi, 1);
                                        dgs_gg_decay1->Fill ((double) DGSEvent[j].ehi, (double) DGSEvent[i].ehi, 1);

                                        /* also require we saw beta in time window */

                                        if (is_beta)
                                          {
                                            tdiff = (int) (DGSEvent[i].event_timestamp - beta_ts);
                                            if (tdiff > Pars.gb_dt_lo && tdiff < Pars.gb_dt_hi)
                                              {
                                                dgs_gg_decay1_beta->Fill ((double) DGSEvent[i].ehi, (double) DGSEvent[j].ehi, 1);
                                                dgs_gg_decay1_beta->Fill ((double) DGSEvent[j].ehi, (double) DGSEvent[i].ehi, 1);
                                              };
                                          };
                                      };

                                  if (tdiff1 > Pars.decay_lo2 && tdiff1 < Pars.decay_hi2)
                                    if (tdiff2 > Pars.decay_lo2 && tdiff2 < Pars.decay_hi2)
                                      {
                                        dgs_gg_decay2->Fill ((double) DGSEvent[i].ehi, (double) DGSEvent[j].ehi, 1);
                                        dgs_gg_decay2->Fill ((double) DGSEvent[j].ehi, (double) DGSEvent[i].ehi, 1);

                                        /* also require we saw beta in time window */

                                        if (is_beta)
                                          {
                                            tdiff = (int) (DGSEvent[i].event_timestamp - beta_ts);
                                            if (tdiff > Pars.gb_dt_lo && tdiff < Pars.gb_dt_hi)
                                              {
                                                dgs_gg_decay2_beta->Fill ((double) DGSEvent[i].ehi, (double) DGSEvent[j].ehi, 1);
                                                dgs_gg_decay2_beta->Fill ((double) DGSEvent[j].ehi, (double) DGSEvent[i].ehi, 1);
                                              };
                                          };
                                      };

                                  if (tdiff1 > Pars.decay_lo3 && tdiff1 < Pars.decay_hi3)
                                    if (tdiff2 > Pars.decay_lo3 && tdiff2 < Pars.decay_hi3)
                                      {
                                        dgs_gg_decay3->Fill ((double) DGSEvent[i].ehi, (double) DGSEvent[j].ehi, 1);
                                        dgs_gg_decay3->Fill ((double) DGSEvent[j].ehi, (double) DGSEvent[i].ehi, 1);

                                        /* also require we saw beta in time window */

                                        if (is_beta)
                                          {
                                            tdiff = (int) (DGSEvent[i].event_timestamp - beta_ts);
                                            if (tdiff > Pars.gb_dt_lo && tdiff < Pars.gb_dt_hi)
                                              {
                                                dgs_gg_decay3_beta->Fill ((double) DGSEvent[i].ehi, (double) DGSEvent[j].ehi, 1);
                                                dgs_gg_decay3_beta->Fill ((double) DGSEvent[j].ehi, (double) DGSEvent[i].ehi, 1);
                                              };
                                          };
                                      };
                                };
                            };
                        };

  /* check the cfd timing (bi207 case) */

  jj = 0;
  if (ng >= 2)
    for (i = 0; i < ng; i++)
      if (Pars.enabled[DGSEvent[i].tid])
        if (DGSEvent[i].tpe == GE)
          if (DGSEvent[i].flag == 0)
            {
              etmp[jj] = DGSEvent[i].ehi;
              ttmp[jj] = DGSEvent[i].event_cfd_timestamp;
              idtmp[jj] = DGSEvent[i].tid;
              jj++;
            };

  if (jj == 2)
    if (fabs (etmp[0] - 569) < 5.0 || fabs (etmp[0] - 1063) < 5.0)
//   if (fabs (etmp[0] - 1173) < 5.0 || fabs (etmp[0] - 1333) < 5.0)
      {
        /* good candidate */

        if (idtmp[0] < idtmp[1])
          indx = idtmp[0] * 110 + idtmp[1];
        else
          indx = idtmp[1] * 110 + idtmp[0];
//      printf("xxx %i %f %f e %f %f %i %i \n", indx, ttmp[0], ttmp[1], etmp[0], etmp[1], idtmp[0], idtmp[1]);
        sptt1->Fill (indx, ttmp[0] - ttmp[1]);
        if (idtmp[0] == 9)
          sptt2->Fill (idtmp[1], ttmp[1] - ttmp[0]);

      };

/* debug list the gamma rays we found */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("final list: we have %i gamma rays\n", ng);
      for (i = 0; i < ng; i++)
        {
          printf ("%2i> ", i);
          printf ("chan_id=%i; ", DGSEvent[i].chan_id);
          printf ("board_id=%i; ", DGSEvent[i].board_id);
          printf ("id =%i; ", DGSEvent[i].id);
          printf ("tpe=%i; ", DGSEvent[i].tpe);
          printf ("tid=%i; ", DGSEvent[i].tid);
          printf ("EventTS=%llu; ", DGSEvent[i].event_timestamp);
          printf ("ehi=%8f ", DGSEvent[i].ehi);
          if (DGSEvent[i].flag == 1)
            printf ("(dirty)");
          else
            printf ("(clean)");
          if (Pars.enabled[gsid])
            printf ("(enabled)");
          else
            printf ("(disabled)");
          printf ("\n");
          fflush (stdout);
        };
    };

  /* done */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("exit bin_dgs\n");

  return (0);

}

/* ---------------------------------------------------------------------*/

TH2F *
make2D (const char *txt, int xln, int xlo, int xhi, int yln, int ylo, int yhi)
{
  char str[STRLEN];

  TH2F *mkTH2F (char *, char *, int, double, double, int, double, double);


  TH2F *h2D;

  sprintf (str, txt);
  h2D = mkTH2F (str, str, xln, xlo, xhi, yln, ylo, yhi);

  return h2D;
}

/*---------------------------------------------------------------------*/

TH1D *
make1D (const char *txt, int xln, int xlo, int xhi)
{
  char str[STRLEN];
  double xlod, xhid;
  TH1D *mkTH1D (char *, char *, int, double, double);
  TH1D *h1D;

  xlod = xlo;
  xhid = xhi;

  sprintf (str, txt);
  h1D = mkTH1D (str, str, xln, xlod, xhid);
  return h1D;
}

/*---------------------------------------------------------------------*/

void
getcal (char *file)
{
  int i, j, ret = 0;
  float b, c, d;
  char mystring[1000];
  FILE *fp;


  /* get pol zero */

  fp = fopen (Pars.dgs_PZfn, "r");      // read mode
  if (fp == NULL)
    {
      printf ("could not open the DGS PZ file: %s; use default == 1\n", Pars.dgs_PZfn);
    }
  else
    {

      // read file and parse

      while (fgets (mystring, 100, fp) != NULL)
        {
          ret = sscanf (mystring, "%i %f ", &i, &b);
          PZ[i] = b;
//      printf ("ge %3i has pz of %8.4f\n", i, PZ[i]);
        }
      fclose (fp);
    };


  /* get energy cal file */

  fp = fopen (Pars.dgs_ecalfn, "r");    // read mode
  if (fp == NULL)
    {
      printf ("could not open the DGS cal file: %s; set 0 offset and 1.0 gain\n", Pars.dgs_ecalfn);
      for (i = 1; i <= 110; i++)
        {
          ehigain[i] = 1;
          ehioffset[i] = 0;
        };
    }
  else
    {

      // read file and parse

      while (fgets (mystring, 100, fp) != NULL)
        {
          ret = sscanf (mystring, "%i %f %f ", &i, &c, &d);
          ehigain[i] = d;
          ehioffset[i] = c;
          printf ("ge %3i has offset and gain of: %8.2f %8.4f and a PZ of %8.4f ", i, ehioffset[i], ehigain[i], PZ[i]);
          if (Pars.enabled[i])
            printf ("enabled\n");
          else
            printf ("DISABLED!\n");
        }
      fclose (fp);
    };

  /* get extrapolation factors */

  fp = fopen (Pars.dgs_factorfn, "r");  // read mode
  if (fp == NULL)
    {
      printf ("could not open the DGS factor file: %s; use default == 0\n", Pars.dgs_factorfn);
      for (j = 1; j <= NGE; j++)
        factor[j] = 0;
    }
  else
    {

      // read file and parse

      while (fgets (mystring, 100, fp) != NULL)
        {
          ret = sscanf (mystring, "%i %f ", &i, &b);
          factor[i] = b;
//      printf ("ge %3i has factor of %8.4f\n", i, PZ[i]);
        }
      fclose (fp);
    };



}
